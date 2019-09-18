#include "Equations_conf.h"

// ============================================
#ifdef RANS_THERMAL_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "MGSolverRANS_thermal.h"

// configuration files -----------
#include "Printinfo_conf.h"

#include "MGGeomEl.h"
#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"
#include "MGFE.h"
#include "MGUtils.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"
#include "parallelM.h"

// ======================================================
/// This function constructs the 3d-2D MGSolTBK class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
 * This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
 * equation name "T", basic variable name "T"
 */

MGSolRANS_thermal::MGSolRANS_thermal(
    MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
    const int nvars_in[],                    ///< KLQ number of variables
    std::string eqname_in,                   ///< equation name
    std::string varname_in                   ///< basic variable name
    )
    : MGSolDA(mg_equations_map_in, nvars_in, eqname_in, varname_in),
      _offset(_mgmesh._NoNodes[_NoLevels - 1]),  // mesh nodes
      _dt(stod(_mgutils._sim_config["dt"])),     // parameter  dt
      _uref(_mgutils._mat_prop["Uref"]),         // parameter  u reference
      _lref(_mgutils._mat_prop["Lref"]),         // parameter  l reference
      _rhof(_mgutils._mat_prop["rho0"]),         // parameter density
      _muf(_mgutils._mat_prop["mu0"]),
      _Tref(_mgutils._mat_prop["Tref"]),  // parameter  temperature reference
      _cp0(_mgutils._mat_prop["cp0"]),    // parameter  Cp reference
      _kappa0(_mgutils._mat_prop["kappa0"]) {
  /// A) reading parameters  for field coupling (in _FF_idx[])
  _RANS_t_parameter.read_param(_mgutils);
  _nTKdim = DIMENSION;
  _IRe = _muf / (_rhof * _uref * _lref);
  _alpha = _kappa0 / (_rhof * _cp0);
  _IPrdl = _alpha / _IRe;
  _qs = _mgutils._mat_prop["qs"];

  for (int k_index = 0; k_index < 30; k_index++) { _FF_idx[k_index] = -1; }

  /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
  int _var_index = 0.;

  std::map<std::string, int> YesNo;
  YesNo["yes"] = 1;
  YesNo["no"] = 0;
  std::map<std::string, int> StabMap;
  StabMap["gls"] = 1;
  StabMap["supgc"] = 0;
  StabMap["sgs"] = -1;
  _SolveTTBK = YesNo[_mgutils._sim_config["SolveThermalTurbulence"]];
  _SUPG = _RANS_t_parameter._Supg;
  _UPWIND = _RANS_t_parameter._Upwind;
  _SCPG = _RANS_t_parameter._Scpg;
  _RNSTAB = _RANS_t_parameter._ReactionNumberBased;
  _NLITER = _RANS_t_parameter._MaxNonLinearIterations;
  _FlatProfile = _RANS_t_parameter._FlatProfile;
  _SolveSteady = _RANS_t_parameter._SolveSteady;
  _AxiSym = (int)(_mgutils._geometry["Axysim"]);
  _WallDist = _mgutils._geometry["Wall_dist"];
  _Restart = 0;

  if (_RANS_t_parameter._ModifiedSupg == "no") {
    _ModifiedSupg = false;
    _theta = 0;
  } else {
    _ModifiedSupg = true;
    _theta = StabMap[_RANS_t_parameter._ModifiedSupg];
  }

  return;
}

//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void MGSolRANS_thermal::GenMatRhs(
    const double /**< time (in) */, const int Level /**< discretization Level (in) */,
    const int mode /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
    )              // ===============================================
{
  /// a) Set up
  const int unsteady_flag = _RANS_t_parameter._SolveSteady;
  // geometry ---------------------------------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];    // element sides
  int el_conn[NDOF_FEM];                               // element connectivity
  int el_neigh[NDOF_FEM];                              // bd element connectivity
  int sur_toply[NDOF_FEMB];                            // boundary topology

  // gauss integration  -----------------------------------------------------------------------------
  double u_old[NDOF_FEM];
  const int el_ngauss = _fe[2]->_NoGauss1[_nTKdim - 1];   // elem gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[_nTKdim - 2];  // bd elem gauss points

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int el_ndof[3];
  el_ndof[0] = 1;
  int elb_ndof[3];
  elb_ndof[0] = 1;       // number of el dofs
  int el_mat_nrows = 0;  // number of mat rows (dofs)

  for (int ideg = 1; ideg < 3; ideg++) {
    el_ndof[ideg] = _fe[ideg]->_NoShape[_nTKdim - 1];
    elb_ndof[ideg] = _fe[ideg]->_NoShape[_nTKdim - 2];
    el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
  };

  const int el_ndof2 = _fe[2]->_NoShape[_nTKdim - 1];

  int el_mat_ncols = el_mat_nrows;  // square matrix

  std::vector<int> el_dof_indices(el_mat_ncols);  // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  for (int k = 0; k < 30; k++) {  // coupling  basic system fields
    const int idx = _data_eq[2].tab_eqs[k];
    _FF_idx[k] = (idx >= 0) ? _data_eq[2].indx_ub[idx] : -1;
  }

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
  A[Level]->zero();

  if (mode == 1) {
    b[Level]->zero();  // global matrix+rhs
  }

  _KeM.resize(el_mat_nrows, el_mat_ncols);
  _FeM.resize(el_mat_nrows);  // resize  local  matrix+rhs

  int ndof_lev = 0;

  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
    ndof_lev += delta;
  }

  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];      // stop element

  for (int iel = 0; iel < (nel_e - nel_b); iel++) {  // LOOP OVER MESH ELEMENTS

    // set to zero matrix and rhs and center
    _KeM.zero();
    _FeM.zero();
    _y_dist = _mgmesh._dist[iel + nel_b];
    // ----------------------------------------------------------------------------------
    /// 1. Geometry and element  fields
    // ----------------------------------------------------------------------------------
    _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, _xx_qnds);
    _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);

    for (int deg = 0; deg < 3; deg++) {  // OLD SOLUTION
      for (int eq = 0; eq < _data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(
            0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq], el_ndof[deg], el_conn, offset,
            _data_eq[deg].indx_ub[eq], _data_eq[deg].ub);
      }
    }

    // ----------------------------------------------------------------------------------
    /// 2. Boundary integration  (bc)
    // ----------------------------------------------------------------------------------

    _EquationNodes.clear();
    for (int k = 0; k < el_ndof[2]; k++) {
      _bc_el[k] = (_bc_vol[k] / 10 == 0) ? 0 : 1;
      if (_bc_el[k] != 0) { _EquationNodes.push_back(k); }
    }

    _WallElement = 0;

    for (int dim = 0; dim < _nTKdim; dim++) { _NormMidCell[dim] = 0.; }

    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int idof = 0; idof < elb_ndof[2]; idof++) {
          sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside];  // local nodes
          int idofb = sur_toply[idof];                                            // connectivity vector

          for (int idim = 0; idim < _nTKdim; idim++) {
            _xxb_qnds[idim * NDOF_FEMB + idof] = _xx_qnds[idim * NDOF_FEM + idofb];  // coordinates
          }
        }
        bc_set(sur_toply, el_ndof[2], elb_ndof[2], elb_ngauss);
      }
    }

    // ----------------------------------------------------------------------------------
    //   3. Volume integration
    // ----------------------------------------------------------------------------------

    vol_integral(el_ndof2, el_ngauss, mode, el_conn);

    // ----------------------------------------------------------------------------------
    //  4. add local to global
    // ----------------------------------------------------------------------------------
    A[Level]->add_matrix(_KeM, el_dof_indices);  // global matrix

    if (mode == 1) {
      b[Level]->add_vector(_FeM, el_dof_indices);  // global rhs
    }
  }  // END LOOP OVER MESH ELEMENTS

  /// 5. clean
  el_dof_indices.clear();
  A[Level]->close();

  if (mode == 1) { b[Level]->close(); }

    //   A[Level]->print(); b[Level]->print();
#ifdef PRINT_INFO
  std::cout << " Matrix Assembled(T)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

  return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolRANS_thermal::MGTimeStep_no_up(
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
  // =========================================================================================

  /// A) Set up the time step
  if (_SolveTTBK == 1) {
    std::cout << std::endl
              << "\033[038;5;" << KTT_F + 50 << ";1m "
              << "--------------------------------------------------- \n\t" << _eqname.c_str() << "   "
              << _var_names[0].c_str() << " solution of problem " << _mgutils.get_name() << " with dir "
              << _dir << "\n --------------------------------------------------- \n \033[0m";

    // SET UP XOOLD AND XNONL VECTORS AFTER RESTART
    if (_Restart == 0) {
      x_old[_NoLevels - 1]->localize(*x_oold[_NoLevels - 1]);
      x_old[_NoLevels - 1]->localize(*x_nonl[_NoLevels - 1]);
    }

    _Restart = 1;

    // EQUATION ASSEMBLY + SOLUTION
    StandardTimeStep(time);
  }

  return;
}  // =======================================================================================

void MGSolRANS_thermal::MGUpdateStep() {
  x_old[_NoLevels - 1]->localize(*x_oold[_NoLevels - 1]);

  int size = x[_NoLevels - 1]->size();

  for (int i = 0; i < size; i++) {
    double k[1] = {(*(x[_NoLevels - 1]))(i)};
    if (k[0] < _LowerLimit) {
      k[0] = _LowerLimit;
      x[_NoLevels - 1]->set(i, k[0]);
    }
  }

  x[_NoLevels - 1]->localize(*x_old[_NoLevels - 1]);

  return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolRANS_thermal::MGTimeStep(
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
  // =========================================================================================

  /// A) Set up the time step
  if (_SolveTTBK == 1) {
    std::cout << std::endl
              << "\033[038;5;" << KTT_F + 50 << ";1m "
              << "--------------------------------------------------- \n\t" << _eqname.c_str() << "   "
              << _var_names[0].c_str() << " solution of problem " << _mgutils.get_name() << " with dir "
              << _dir << "\n --------------------------------------------------- \n \033[0m";

    // SET UP XOOLD AND XNONL VECTORS AFTER RESTART
    if (_Restart == 0) {
      x_old[_NoLevels - 1]->localize(*x_oold[_NoLevels - 1]);
      x_old[_NoLevels - 1]->localize(*x_nonl[_NoLevels - 1]);
    }

    _Restart = 1;
    StandardTimeStep(time);

    // UPDATE XOLD AND XOOLD VECTORS
    MGUpdateStep();
  }

  return;
}  // =======================================================================================

void MGSolRANS_thermal::StandardTimeStep(const double time) {
  std::clock_t start_time = std::clock();
  GenMatRhs(time, _NoLevels - 1, 1);  // matrix and rhs
  for (int Level = 0; Level < _NoLevels - 1; Level++) {
    GenMatRhs(time, Level, 0);  // matrix
  }
  std::clock_t end_time = std::clock();
  MGSolve(1.e-6, 40);
  std::clock_t end_time2 = std::clock();

#if PRINT_TIME == 1
  std::cout << " Ass. time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC << "s "
            << " Ass. and sol. time: =" << double(end_time2 - start_time) / CLOCKS_PER_SEC << "s "
            << std::endl;
#endif

  return;
}

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;
