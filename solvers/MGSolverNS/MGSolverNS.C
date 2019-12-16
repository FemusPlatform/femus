// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#include "MGSolverNS.h"  // Navier-Stokes class header file

#include "MGFE_conf.h"       // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

#include "EquationSystemsExtendedM.h"  // Equation map class
#include "MGFE.h"                      // Mesh class
#include "MeshExtended.h"

// local alg lib -----------------------------------------------
#include "linear_solverM.h"   // algebra solvers
#include "numeric_vectorM.h"  // algebra numerical vectors
#include "sparse_matrixM.h"   // algebra sparse matrices

// ==================================================================
/// This routine constructs the FSI class:
MGSolNS::MGSolNS(
    MGEquationsSystem& mg_equations_map_in, int nvars_in[], std::string eqname_in, std::string varname_in)
    : MGSolDA(mg_equations_map_in, nvars_in, eqname_in, varname_in),
      //
      //===================================================================================================//
      //                                    A) Reading parameters                                          //
      //===================================================================================================//
      //
      _offset(_mgmesh._NoNodes[_NoLevels - 1]),  // mesh nodes (top level)
      _dt(stod(_mgutils._sim_config["dt"])),     // parameter  dt
      _uref(_mgutils._mat_prop["Uref"]),         // parameter  u reference
      _lref(_mgutils._mat_prop["Lref"]),         // parameter  l reference
      _rhof(_mgutils._mat_prop["rho0"]),         // parameter density
      _muf(_mgutils._mat_prop["mu0"]) {          // parameter viscosity
  //===================================================================================================//
  //                                    B) Setting class variables                                     //
  //===================================================================================================//

  // READ PARAMETERS FROM CLASS NS_parameter
  _NS_parameter.read_param(_mgutils);
  _NumRestartSol = _NS_parameter._NumRestartSol;
  _nNSdim = DIMENSION;

  // class equation ---------------------------------------------------------------------
  for (int k_index = 0; k_index < 30; k_index++) { _FF_idx[k_index] = -1; }

  _pres_order = (_nvars[0] > 0) ? 0 : 1;
  // class variable names ------------------------------------------------------------
  //   _dir = 0;

  //===================================================================================================//
  //                                    C) Setting solver type                                         //
  //===================================================================================================//
  for (int l = 0; l < _NoLevels; l++) { _solver[l]->set_solver_type(_NS_parameter._SolverType); }

  _IRe = _muf / (_rhof * _lref * _uref);   // Reynolds number
  _IFr = 9.81 * _lref / (_uref * _uref);   // Froud number
  _dirg[0] = _mgutils._mat_prop["dirgx"];  // x-gravity
  _dirg[1] = _mgutils._mat_prop["dirgy"];  // y-gravity
  _dirg[2] = _mgutils._mat_prop["dirgz"];  // z-gravity

  for (int i = 0; i < _nNSdim; i++) { _xyzg[i] = 0.; }

  _SolveNS = (_mgutils._sim_config["SolveNavierStokes"].compare("yes") == 0) ? true : false;
  _Wall_dist = _mgutils._geometry["Wall_dist"];
  _AxiSym = (int)(_mgutils._geometry["Axisym"]);

  if (_AxiSym == 1 && _nNSdim == 3) {
    std::cout << "\033[1;31m MGSolNS: CANNOT USE AXISYM OPTION WITH 3D GEOMETRY \n \033[0m";
    std::abort();
  }

  _ImmersedBoundary = 0;

  return;
}  //****************************************************************************************************//

// ====================================================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================================================
void MGSolNS::GenMatRhs(const double /* time*/, const int Level, const int mode) {
  // ===================================================================================================
  //                                    A) Set up
  // ===================================================================================================

  // -------------------------- Geometry -------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];    // element nodes
  int el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];          // element connectivity
  int el_neigh[NDOF_FEM];                              // element connectivity
  int sur_toply[NDOF_FEMB];                            // boundary topology
  double normal[DIMENSION];                            // normal to the boundary

  // -------------------------- Gauss integration -----------------------------------------------------

  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION - 2];  // elem gauss points

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] ---------------------------------------
  int el_ndof[3];
  el_ndof[0] = NDOF_K;
  int elb_ndof[3];
  elb_ndof[0] = 0;       // number of el dofs
  int el_mat_nrows = 0;  // number of mat rows (dofs)

  for (int ideg = 1; ideg < 3; ideg++) {                                           //     ...
    el_ndof[ideg] = ((_nvars[ideg] > 0) ? _fe[ideg]->_NoShape[_nNSdim - 1] : 0);   //   computing
    elb_ndof[ideg] = ((_nvars[ideg] > 0) ? _fe[ideg]->_NoShape[_nNSdim - 2] : 0);  //     ...
    el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
  }

  el_ndof[1] = _fe[1]->_NoShape[_nNSdim - 1];

  el_mat_nrows += el_ndof[0] * _nvars[0];
  int el_mat_ncols = el_mat_nrows;                // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);  // element dof vector

  // fields -> Navier-Stokes ----------------------------------------------------------------------
  double x_m[DIMENSION];

  for (int k = 0; k < 30; k++) {  // coupling  basic system fields
    const int idx = _data_eq[2].tab_eqs[k];
    _FF_idx[k] = (idx >= 0) ? _data_eq[2].indx_ub[idx] : -1;
  }

  if (_FF_idx[IB_F] >= 0) _ImmersedBoundary = 1;

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
  A[Level]->zero();

  if (mode == 1) {
    b[Level]->zero();  // global matrix+rhs
  }

  _KeM.resize(el_mat_nrows, el_mat_ncols);
  _FeM.resize(el_mat_nrows);  // resize  local  matrix+rhs

  int ndof_lev = 0;  // for multilevel multi proc point ordering ------------------------------------------

  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    ndof_lev += _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
  }

  // ===================================================================================================
  //                                    B) Element  Loop over the volume (n_elem)
  // ===================================================================================================

  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];      // stop element

  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    // set to zero matrix and rhs and center01
    _KeM.zero();
    _FeM.zero();

    // ---------------------------------------------------------------------------------
    // Volume  geometry and element  fields
    // ---------------------------------------------------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, _xx_qnds);
    _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);
    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);
    // field data  ------------------------------------------------------
    get_el_field_data(iel, Level, el_conn, offset, el_ndof, ndof_lev);

    CalcVolume();

    // initializing  volume quantities
    for (int idim = 0; idim < _nNSdim; idim++) {  // quad loop entities (vector)
      x_m[idim] = 0.;

      for (int d = 0; d < el_ndof[2]; d++) {
        const int dnode = idim * el_ndof[2] + d;  // index local points
        x_m[idim] += _xx_qnds[dnode] / el_ndof[2];
        _normal_pt[dnode] = 0.;  // normal boundary
        _bc_el[dnode] = 0;
        _AlreadyWrittenDirBC[dnode] = false;
      }
    }

    for (int d = 0; d < el_ndof[1]; d++) {  // pressure  loop entities (scalar)
      if (_bc_bd[_nNSdim * el_ndof[2] + d] % 10 > 3) {
        _bc_el[_nNSdim * el_ndof[2] + d] = 0;
      } else {
        _bc_el[_nNSdim * el_ndof[2] + d] = 1;
      }
    }

    // SETTING _bc_el = -8 FOR EVERY ELEMENT BOUNDARY NODE
    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int lbnode = 0; lbnode < elb_ndof[2]; lbnode++) {
          int lnode = _mgmesh._GeomEl._surf_top[lbnode + elb_ndof[2] * iside];  // local nodes

          for (int idim = 0; idim < _nNSdim; idim++) { _bc_el[lnode + idim * el_ndof[2]] = _BdFlagId; }
        }
      }
    }

    _WallElement = 0;

    if (_FF_idx[IB_F] >= 0)
      _wall_frac = _mgmesh._VolFrac[iel + nel_b];
    else
      _wall_frac = 0.;

    // --------------------------------------------------------------------------
    //  Boundary and boundary bc
    // --------------------------------------------------------------------------
    //  New bc flags for Navier Stokes equation:
    //  _bc_el = -1 -2 -3 -> Neumann bc along tangential 1, tangential 2 (only 3d), normal directions
    //  _bc_el =  1  2  3 -> Dirichlet bc along tangential 1, tangential 2 (only 3d), normal directions
    //  _bc_el =  0 -> normal interior point
    //  Integration is performed for _bc_el <=0
    //  Test functions are set to zero when _bc_el >0
    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        // setup boundary element  ----------------------------------------------------------------
        for (int lbnode = 0; lbnode < elb_ndof[2]; lbnode++) {                  // quad quantities
          int lnode = _mgmesh._GeomEl._surf_top[lbnode + elb_ndof[2] * iside];  // local nodes
          sur_toply[lbnode] = lnode;                                            // lbnode -> lnode
          elb_conn[lbnode] = el_conn[lnode];  // connctivity el_conn->elb_conn

          for (int idim = 0; idim < _nNSdim; idim++) {  // coordinates
            _xxb_qnds[idim * elb_ndof[2] + lbnode] = _xx_qnds[idim * el_ndof[2] + lnode];
          }
        }

        _fe[2]->normal_g(_xxb_qnds, x_m, normal);
        int dir_maxnormal = (fabs(normal[0]) > fabs(normal[1])) ? 0 : 1;
        dir_maxnormal =
            (fabs(normal[dir_maxnormal]) > fabs(normal[DIMENSION - 1])) ? dir_maxnormal : DIMENSION - 1;
        //         std::cout<<" .........................................................\n Element "<<iel<<"
        //         side "<< iside <<std::endl;
        set_bc_matrix(dir_maxnormal, sur_toply, el_ndof, elb_ndof, elb_ngauss, normal, el_conn);
      }  // iside -1
    }    // -----------------------------  End Boundary -------------------------------------

    matrixrhsvol(el_ndof, mode, el_conn);

    // ----------------------------------------------------------------------------------
    //   E) Add them to the global matrix
    A[Level]->add_matrix(_KeM, el_dof_indices);

    if (mode == 1) { b[Level]->add_vector(_FeM, el_dof_indices); }
  }  //  =============== End of element loop =============================================

  // ===============================================================================
  el_dof_indices.clear();
  A[Level]->close();

  if (mode == 1) { b[Level]->close(); }

// ----------------------------------------------------------
#ifdef PRINT_INFO
  //     A[Level]->print();  b[Level]->print();
  std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) " << std::endl;
#endif

  return;
}  //****************************************************************************************************//

// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolNS::MGTimeStep_no_up(
    const double time,  // time
    const int max_iter  // Number of max inter
) {
  // =========================================================================================

  // ========================================================================================= //
  //              A) Set up the time step                                                      //
  // ========================================================================================= //
  if (_SolveNS) {
    std::cout << std::endl
              << "\033[038;5;" << NS_F + 50 << ";1m "
              << "--------------------------------------------------- \n\t" << _eqname.c_str()
              << " solution of problem " << _mgutils.get_name() << " with dir " << _dir
              << "\n ---------------------------------------------------\n\033[0m";

    // ========================================================================================= //
    //              B) Assemblying of the Matrix-Rhs                                             //
    // ========================================================================================= //

#if PRINT_TIME == 1
    std::clock_t start_time = std::clock();
#endif

    GenMatRhs(time, _NoLevels - 1, 1);  // matrix and rhs

    for (int Level = 0; Level < _NoLevels - 1; Level++) {
      GenMatRhs(time, Level, 0);  // matrix
    }

#if PRINT_TIME == 1
    std::clock_t end_time = std::clock();
    std::cout << "  Assembly time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC << " s "
              << std::endl;
#endif
    // ========================================================================================= //
    //              C) Solution of the linear MGsystem (MGSolFSI::MGSolve)                       //
    // ========================================================================================= //
    MGSolve(1.e-6, 15);
#if PRINT_TIME == 1
    end_time = std::clock();
    std::cout << " Assembly+solution time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC << "s "
              << std::endl;
#endif
    // ========================================================================================= //
    //              D) Non Linear Iterations                                                     //
    // ========================================================================================= //

    int iter_nl = _NS_parameter._MaxNonLinearIterations;

    for (int iter = 0; iter < iter_nl; iter++) {
      x[_NoLevels - 1]->localize(*x_nonl[_NoLevels - 1]);
      std::cout << std::endl
                << " === NON LINEAR ITERATOR: " << iter + 1 << "- " << _eqname.c_str() << " solution "
                << std::endl;
      GenMatRhs(time, _NoLevels - 1, 1);  // matrix and rhs

      for (int Level = 0; Level < _NoLevels - 1; Level++) {
        GenMatRhs(time, Level, 0);  // matrix
      }

      MGSolve(1.e-6, 15);  // solve
    }

    _TimeStep++;
  }

  return;
}  //****************************************************************************************************//

void MGSolNS::MGUpdateStep() {
  if (_SolveNS) {
    x_old[1][_NoLevels - 1]->localize(*x_old[2][_NoLevels - 1]);  // time step -2
    x_old[0][_NoLevels - 1]->localize(*x_old[1][_NoLevels - 1]);  // time step -1
    x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);         // time step
    x[_NoLevels - 1]->localize(*x_nonl[_NoLevels - 1]);
  }

  return;
}  //****************************************************************************************************//

// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolNS::MGTimeStep(
    const double time,  // time
    const int max_iter  // Number of max inter
) {
  // =========================================================================================
  MGTimeStep_no_up(time, max_iter);
  MGUpdateStep();
  return;
}  //****************************************************************************************************//

// ================================================================
// ====================== END velocity ============================

// ================================================================

// ========================================================================
// ========================================================================
//   NS EQUATION PARAMETERS
// ========================================================================
// ========================================================================

void MGSolNS::CalcTangDir(double Tang1[], double Tang2[], double normal[], TANG_CALC_TYPE Type) {
  CalcTangDir(0, NDOF_FEM, 10, 0, Tang1, Tang2, normal, Type);
  return;
}

void MGSolNS::CalcTangDir(
    int PhiNumber, int el_dof, int bc_el, int qp, double Tang1[], double Tang2[], TANG_CALC_TYPE Type) {
  double normal[DIMENSION];

  for (int kk = 0; kk < DIMENSION; kk++) { normal[kk] = _normal_pt[kk + PhiNumber * _nNSdim]; }

  CalcTangDir(PhiNumber, el_dof, bc_el, qp, Tang1, Tang2, normal, Type);
  return;
}

void MGSolNS::CalcTangDir(
    int PhiNumber, int el_dof, int bc_el, int qp, double Tang1[], double Tang2[], double normal[],
    TANG_CALC_TYPE Type) {
  // ---------------------------------------------------
  if (Type == VEL_BASED) {
    double scal_unorm = 0., mod_norm = 0., mod_vel = 0.;
    double Normal[DIMENSION], Vel[DIMENSION];

    for (int kk = 0; kk < DIMENSION; kk++) {
      Normal[kk] = normal[kk];
      mod_norm += Normal[kk] * Normal[kk];
      Vel[kk] = _data_eq[2].ub[(_FF_idx[NS_F] + kk) * el_dof + (NDOF_FEM - 1)];
      mod_vel += Vel[kk] * Vel[kk];
    }

    mod_norm = sqrt(mod_norm + 1.e-20);
    mod_vel = sqrt(mod_vel + 1.e-20);

    for (int kk = 0; kk < DIMENSION; kk++) {
      Normal[kk] = Normal[kk] / mod_norm;
      scal_unorm += Vel[kk] * Normal[kk];
    }

    double prod_u_t1 = 0.;

    for (int kk = 0; kk < DIMENSION; kk++) {
      Tang1[kk] = (Vel[kk] - scal_unorm * Normal[kk]) / mod_vel;
      prod_u_t1 += Vel[kk] * Tang1[kk];
    }

    // check if tang1 has same direction of velocity field
    if (prod_u_t1 < 0.)
      for (int kk = 0; kk < DIMENSION; kk++) { Tang1[kk] = -Tang1[kk]; }

    if (_nNSdim == 3) {  //
      for (int kk = 0; kk < DIMENSION; kk++) {
        int a = (kk + 1) % DIMENSION;
        int b = (kk + 2) % DIMENSION;
        Tang2[kk] = (Normal[a] * Tang1[b] - Normal[b] * Tang1[a]);
      }
    }
  }

  // -------------------------------------------------------------
  if (Type == GEOM_BASED) {
    Tang1[0] = -(normal[1] + normal[DIMENSION - 1] * (DIMENSION % 2));
    Tang1[1] = normal[0];
    Tang1[DIMENSION - 1] = normal[0];

    for (int k = 0; k < DIMENSION; k++) {
      int a = (k + 1) % DIMENSION;
      int b = (k + 2) % DIMENSION;
      Tang2[k] = (normal[a] * Tang1[b] - normal[b] * Tang1[a]);
    }
  }

  return;
}

// =======================================================================================================
void MGSolNS::SetBCFlags(
    double normal[],  // normal unit
    int sur_toply[],  // surface connectivity (respect to volume)
    int el_ndof,      // number of dofs
    int el_conn[]     // volume connectivity
) {  // -----------------------------------------------------------------------------------------------------
  // bc conditions
  // Dirichlet >5  (surface point setting ( no volume integration))
  // 6-> set zero component 7-> set component       ------------------------> _bc_el[] = 3,0,-3; _bc_el[]
  // =1,2,-1,-2; 8-> zero along dir (normal or tg) 9-> value along dir  ----------------> _bc_el[] = 3,0,-3;
  // _bc_el[] =1,2,-1,-2;
  //  Neumann (volume integration (+surface integration))   ---------------->
  // 1-> no surface integration                             ---------------->
  // 2-> pressure integration (normal only)  3-> pressure setting    ------->
  // 4-> tau=a.velocity                      4-> stau=a.(velocity-velocity_0)------>
  const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;  // element normal-tg bc flag
  const int norm_face =
      (abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 10000) / 1000 - 1;  // element max-normal-dir bc flag

  for (int lbnode = 0; lbnode < NDOF_FEMB; lbnode++)
    if (el_conn[lbnode] == 3) int a = 1;

  for (int lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {                // loop on surface nodes
    const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;         // point normal-tg  bc flag
    const int bc_var_check = abs(_bc_bd[sur_toply[lbnode]]) % 10000;  // bc_var
    int bc_var_normal = bc_var_check / 1000 - 1;                      // point max-normal-dir  bc flag
    const int bc_n_flag = bd_node / 10;                               // normal bc flag
    const int bc_tg_flag = bd_node % 10;                              // tg bc flag

    bool Norm_Dirichlet_Bound_Cond = (bc_n_flag > 5) ? true : false;
    bool Tang_Dirichlet_Bound_Cond = (bc_tg_flag > 5) ? true : false;
    int NormDir = bc_var_normal;

    // BC are imposed only if node bc equal to face bc
    if (bd_face == bd_node || (bd_node == 88 || bd_node == 66)) {
      int NormRow = sur_toply[lbnode] + NormDir * el_ndof;

      // Setting bc flags for normal direction --------------------------
      if (_bc_el[NormRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
        if (Norm_Dirichlet_Bound_Cond) {
          _bc_el[NormRow] = 3;
        } else {  // Neumann bc (surface integration)
          if (bc_n_flag == 1 || bd_node == 31) {
            _bc_el[NormRow] = 0;
          } else {
            _bc_el[NormRow] = -3;
          }
        }
      }  //  norm  ----------------------------------------------------------

      // Setting bc flags for tangential directions -------------------
      int dir_tg0 = 0, NeuTg = 0;

      for (int kdir = NormDir + 1; kdir < NormDir + _nNSdim; kdir++) {
        int idir = kdir % _nNSdim;
        int tg_dirRow = sur_toply[lbnode] + idir * el_ndof;

        if (_bc_el[tg_dirRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
          if (Tang_Dirichlet_Bound_Cond) {
            dir_tg0++;  // Dirichlet bc
            _bc_el[tg_dirRow] = dir_tg0;
          } else {  // Neumann bc (surface integration)
            if (bd_node == 11 || bd_node == 31) {
              _bc_el[tg_dirRow] =
                  0;  // for outflow bc we don't project the equation along normal-tangential directions
            } else {
              NeuTg++;
              _bc_el[tg_dirRow] = -NeuTg;
            }
          }
        }
      }  //  Tan ----------------------------------------------------------
    }    // if (bd_face==bd_node)
  }

  return;
}

// ===========================================================================
void MGSolNS::CorrectBCFlags(
    double Tang1[],   // tan in dir 1 (flow align)
    double Tang2[],   // tan in dir 2 (transverse flow)
    int MaxNormal,    // dir max normal
    int ElementNode,  // node id
    int ElDof,        // elemet dofs
    int BcType        // 0 dirichlet 1 neumann
) {
  int tg1 = (MaxNormal + 1) % DIMENSION;
  int tg2 = (MaxNormal + 2) % DIMENSION;

  int tang1_maxdir = (fabs(Tang1[tg1]) > fabs(Tang1[tg2])) ? tg1 : tg2;
  int tang2_maxdir = (fabs(Tang2[tg1]) > fabs(Tang2[tg2])) ? tg1 : tg2;

  if (tang1_maxdir == tang2_maxdir) {
    int last_dir = (tang1_maxdir == tg1) ? tg2 : tg1;

    if (fabs(Tang1[tang1_maxdir]) > fabs(Tang2[tang2_maxdir])) {
      tang2_maxdir = last_dir;
    } else {
      tang1_maxdir = last_dir;
    }
  }

  if (BcType == 1) {  // neumann
    if (_bc_el[ElementNode + (tg1)*ElDof] * _bc_el[ElementNode + tg2 * ElDof] == 2 &&
        _bc_el[ElementNode + tg1 * ElDof] + _bc_el[ElementNode + tg2 * ElDof] == -3) {
      _bc_el[ElementNode + (tang1_maxdir)*ElDof] = -1;
      _bc_el[ElementNode + (tang2_maxdir)*ElDof] = -2;
    }
  } else if (BcType == 0) {  // dirichlet
    if (_bc_el[ElementNode + (tg1)*ElDof] * _bc_el[ElementNode + tg2 * ElDof] == 2 &&
        _bc_el[ElementNode + tg1 * ElDof] + _bc_el[ElementNode + tg2 * ElDof] == 3) {
      _bc_el[ElementNode + (tang1_maxdir)*ElDof] = 1;
      _bc_el[ElementNode + (tang2_maxdir)*ElDof] = 2;
    }
  }

  return;
}

#endif  // ENDIF NS_EQUATIONS
// #endif  // NS_equation is personal
