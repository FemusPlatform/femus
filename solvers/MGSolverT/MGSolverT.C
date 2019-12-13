#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS  // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverT.h"

// configuration files -----------
#include "Printinfo_conf.h"

// standard  library
#include <sstream>
#include "MGFE.h"
#include "MGGeomEl.h"
#include "MGUtils.h"
#include "MeshExtended.h"
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "linear_solverM.h"
#include "numeric_vectorM.h"
#include "parallelM.h"
#include "sparse_matrixM.h"

// ======================================================
/// This function constructs the 3d-2D MGSolT class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
 * This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
 * equation name "T", basic variable name "T"
 */
MGSolT::MGSolT(
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
      _muf(_mgutils._mat_prop["mu0"]),           // parameter viscosity
      _Tref(_mgutils._mat_prop["Tref"]),         // parameter  temperature reference
      _cp0(_mgutils._mat_prop["cp0"]),           // parameter  Cp reference
      _kappa0(_mgutils._mat_prop["kappa0"]) {    // parameter  conductivity reference
  //  =========================================================================

  // READ PARAMETERS FROM CLASS T_parameter
  _T_parameter.read_param(_mgutils);

  /// A) reading parameters  for field coupling (in _FF_idx[])
  _nTdim = DIMENSION;
  for (int k_index = 0; k_index < 30; k_index++) { _FF_idx[k_index] = -1; }
  /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
  _var_names[0] = varname_in;
  _refvalue[0] = _Tref;

  /// C ) Setting the  solver type (with _solver[l]->set_solver_type(SOLVERT))
  for (int l = 0; l < _NoLevels; l++) { _solver[l]->set_solver_type(_T_parameter._SolverType); }

  /// D) Setting nondimensional parameters _alpha _IPrdl _IRe .....
  _alpha = _kappa0 / (_rhof * _cp0);
  _IPrdl = _rhof * _alpha / _muf;
  _IRe = _muf / (_rhof * _uref * _lref);

  _IPrdl_turb = 1. / _T_parameter._Prt;
  _alpha_turb = 0.;

  _qheat = _mgutils._mat_prop["qheat"] * _lref / (_rhof * _cp0 * _Tref * _uref);
  _qs = _mgutils._mat_prop["qs"] / (_rhof * _cp0 * _Tref * _uref);
  _Wall_dist = _mgutils._geometry["Wall_dist"];

  _SolveT = (_mgutils._sim_config["SolveTemperature"].compare("yes") == 0) ? true : false;

  _Axisym = _mgutils._geometry["Axisym"];

  _NumRestartSol = _T_parameter._NumRestartSol;

  return;
}

//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void MGSolT::GenMatRhs(
    const double /**< time (in) */, const int Level /**< discretization Level (in) */,
    const int mode /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
) {                // ===============================================
                   //   double Crank_Nicolson =1.;
  /// a) Set up

  // geometry ---------------------------------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];    // element sides
  int el_conn[NDOF_FEM];                               // element connectivity
  int el_neigh[NDOF_FEM];                              // bd element connectivity
  int sur_toply[NDOF_FEMB];                            // boundary topology

  // gauss integration  -----------------------------------------------------------------------------
  double x_m[DIMENSION];
  double normal[DIMENSION];
  const int el_ngauss = _fe[2]->_NoGauss1[_nTdim - 1];   // elem gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[_nTdim - 2];  // bd elem gauss points

  double WallDist[NDOF_FEM], AlphaTurb[NDOF_FEM];

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int el_ndof[3];
  el_ndof[0] = 1;
  int elb_ndof[3];
  elb_ndof[0] = 1;       // number of el dofs
  int el_mat_nrows = 0;  // number of mat rows (dofs)
  for (int ideg = 1; ideg < 3; ideg++) {
    el_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 1];
    elb_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 2];
    el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
  };
  const int el_ndof2 = _fe[2]->_NoShape[_nTdim - 1];

  int el_mat_ncols = el_mat_nrows;                // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);  // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  for (int k = 0; k < 30; k++) {  // coupling  basic system fields
    const int idx = _data_eq[2].tab_eqs[k];
    _FF_idx[k] = (idx >= 0) ? _data_eq[2].indx_ub[idx] : -1;
  }
  double vel_g[DIMENSION];
  for (int idim = 0; idim < _nTdim; idim++) {
    vel_g[idim] = 0.;  // velocity not coupled
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
  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    // set to zero matrix and rhs and center
    _KeM.zero();
    _FeM.zero();

    /// 1. geometry and element  fields ------------------------------------
    // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
    _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, _xx_qnds);
    _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);
    // fill the node data vectors
    for (int deg = 0; deg < 3; deg++) {
      for (int eq = 0; eq < _data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(
            0, 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq], el_ndof[deg], el_conn, offset,
            _data_eq[deg].indx_ub[eq], _data_eq[deg].ub);
      }
    }

    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[T_F]]->get_el_sol(
        0, 0, 1, el_ndof[2], el_conn, offset, 0, _T_1ts);  // time step -1
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[T_F]]->get_el_sol(
        1, 0, 1, el_ndof[2], el_conn, offset, 0, _T_2ts);  // time step -2

    // ----------------------------------------------------------------------------------
    /// 2. Boundary integration  (bc)
    // ----------------------------------------------------------------------------------

    for (int k = 0; k < el_ndof[2]; k++) {
      _bc_el[k] = (_bc_vol[k] / 10 == 0) ? 0 : 1;  // boundary condition
      _bc_bd[k] = 1;
    }
    for (int idim = 0; idim < _nTdim; idim++) {
      x_m[idim] = 0.;
      for (int idofk = 0; idofk < el_ndof[2]; idofk++)
        x_m[idim] += _xx_qnds[idim * NDOF_FEM + idofk] / NDOF_FEM;
    }

    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int idof = 0; idof < elb_ndof[2]; idof++) {
          sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside];  // local nodes
          int idofb = sur_toply[idof];                                            // connectivity vector
          _bc_bd[idofb] = 0;
          for (int idim = 0; idim < _nTdim; idim++) {
            _xxb_qnds[idim * NDOF_FEMB + idof] = _xx_qnds[idim * NDOF_FEM + idofb];  // coordinates
          }
        }
        int sign_normal = 1;
        _fe[2]->normal_g(_xxb_qnds, x_m, normal, sign_normal);
        bc_set(sur_toply, el_ndof[2], elb_ndof[2], elb_ngauss, sign_normal);
      }
    }
    // ----------------------------------------------------------------------------------
    //   3. Volume integration
    // ----------------------------------------------------------------------------------

    if (_FF_idx[K_F] >= 0) _y_dist = _mgmesh._dist[iel + nel_b];
    if (_FF_idx[IB_F] >= 0)
      _wall_frac = _mgmesh._VolFrac[iel + nel_b];
    else
      _wall_frac = 0.;

    vol_integral(el_ndof2, el_ngauss, mode);

    A[Level]->add_matrix(_KeM, el_dof_indices);  // global matrix
    if (mode == 1) {
      b[Level]->add_vector(_FeM, el_dof_indices);  // global rhs
    }
  }  // end of element loop

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
void MGSolT::MGTimeStep(
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
  // =========================================================================================

  if (_SolveT) {
    std::cout << std::endl
              << "\033[038;5;" << 196 << ";1m "
              << "--------------------------------------------------- \n\t" << _eqname.c_str()
              << " solution of problem " << _mgutils.get_name()
              << "\n ---------------------------------------------------\n  \033[0m";

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

    x_old[0][_NoLevels - 1]->localize(*x_old[1][_NoLevels - 1]);
    x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);
  }
  return;
}  // =======================================================================================

#endif
// #endif // personal application

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
