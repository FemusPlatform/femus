// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include "MGSolverP.h"       // Navier-Stokes class header file
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "EquationSystemsExtendedM.h"  // Equation map class
#include "Solvertype_enum.h"

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers


// ================================================================
#if (NS_EQUATIONS%2)==0

// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================

MGSolP::MGSolP (
  MGEquationsSystem & mg_equations_map_in,
  const int             nvars_in[],    ///
  std::string     eqname_in,    ///< name equation (in)
  std::string     varname_in    ///< name variable (in)
) :  MGSolDA ( mg_equations_map_in, nvars_in, eqname_in, varname_in ),
  _offset ( _mgmesh._NoNodes[-1] ),          // mesh nodes (top level)
  _dt ( stod ( _mgutils._sim_config["dt"] ) ),                  // parameter  dt
  _uref ( _mgutils._mat_prop["Uref"] ),      // parameter  u reference
  _lref ( _mgutils._mat_prop["Lref"] ),      // parameter  l reference
  _rhof ( _mgutils._mat_prop["rho0"] ),      // parameter density
  _muf ( _mgutils._mat_prop["mu0"] )        // parameter viscosity
  {
  _nPdim = DIMENSION;
  _factor = 0.01;

  for ( int k_index = 0; k_index < 30; k_index++ ) {
      _FF_idx[k_index] = -1;
      }

  _var_names[0] = "p";
  _refvalue[0] = _rhof * _uref * _uref; // class variable names
  _IRe = _muf / _rhof;

  for ( int  l = 0; l < _NoLevels; l++ ) {// BICGSTABM  BICGM
      _solver[l]->set_solver_type ( GMRESM );
      }

  _P_parameter.read_param ( _mgutils, _mgmesh._iproc );
  _AssembleOnce = _P_parameter._AssembleOnce;
  _AlreadyAssembled = 0;
  return;
  }
//
//  ====================================================
//  This function assembles the matrix and the rhs:
//  ====================================================
//
void  MGSolP::GenMatRhs (
  const double /* time*/, // time  <-
  const int  Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
)    // ===============================================
  {
  int axis = ( int ) ( _mgutils._geometry["Axisym"] );
  const int   offset = _mgmesh._NoNodes[_NoLevels - 1];                               // mesh nodes
  // element connectivity
  const int  el_sides = _mgmesh._GeomEl._n_sides[0];                         // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                          // element connectivity
  int        el_neigh[NDOF_FEM];                                              // element connectivity

  /*----------------------------------- Gauss integration  --------------------------------------------*/

  double det, JxW_g, InvJac[3][DIMENSION * DIMENSION];                          // Jac, Jac*w Jacobean
  const int  el_ngauss = _fe[2]->_NoGauss1[ _nPdim - 1];                              //elem gauss points

  double dphijdx_g[3][DIMENSION];
  double dphiidx_g[3][DIMENSION];             // global derivatives at g point

  double u_1ts[DIMENSION * NDOF_FEM], u_2ts[DIMENSION * NDOF_FEM], u_3ts[DIMENSION * NDOF_FEM], u_div[DIMENSION * NDOF_FEM];
  double p_1ts[NDOF_P], p_2ts[NDOF_P], p_rhs[NDOF_P];

  double vel_gdx[DIMENSION * DIMENSION], p_gdx[DIMENSION];

  // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
  int el_ndof[3];
  el_ndof[0] = 1;                // number of el dofs
  int el_mat_nrows = 0;                           // number of mat rows (dofs)

  for ( int ideg = 1; ideg < 3; ideg++ ) {
      el_ndof[ideg] = _fe[ideg]->_NoShape[ _nPdim - 1];
      el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
      };

  int el_mat_ncols = el_mat_nrows;                     //square matrix

  std::vector<int > el_dof_indices ( el_mat_ncols );   // element dof vector

  /*---------------------------- Fields -> Navier-Stokes  [NS_F]  -------------------------------------*/
  for ( int k = 0; k < 30; k++ ) { // coupling  basic system fields
      const int idx = _data_eq[2].tab_eqs[k];
      _FF_idx[k] = ( idx >= 0 ) ? _data_eq[2].indx_ub[idx] : -1;
      }

  /* --------------------------- Element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --*/

  A[Level]->zero();

  if ( mode == 1 ) {
      b[Level]->zero();    // global matrix+rhs
      }

  DenseMatrixM KeM;
  DenseVectorM FeM;                // local  matrix+rhs
  KeM.resize ( el_mat_nrows, el_mat_ncols );
  FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

  int ndof_lev = 0;

  for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
      int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
      ndof_lev += delta;
      }

  const int  nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
  const int  nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

  double div_1ts, div_2ts, div_3ts;
  double pres_1ts, pres_2ts;

  if ( _P_parameter._TimeDisc == 2 ) {
      div_1ts = 3. / 2.;
      div_2ts = -2.;
      div_3ts = 1. / 2.;
      pres_1ts = 2.;
      pres_2ts = -1.;
      }

  if ( _P_parameter._TimeDisc == 1 ) {
      div_1ts = 1.;
      div_2ts = -1.;
      div_3ts = 0.;
      pres_1ts = 1.;
      pres_2ts = 0.;
      }

  for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {

      // ===================================================================================================
      //                          B) Element  Loop over the volume (n_elem)
      // ===================================================================================================

      // ------------------------ Set to zero matrix and rhs and center -----------------------------------
      KeM.zero();
      FeM.zero();
      // ------------------------ Geometry and element fields ---------------------------------------------
      // ------------------------ Element Connectivity (el_conn) and coordinates (xx_qnds) ----------------
      _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, _xx_qnds );
      _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh );
      get_el_dof_bc ( Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

      for ( int idim = 0; idim < _nPdim; idim++ ) {
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_1ts );
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_oldsol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_2ts );
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_oooldsol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_3ts );
          }

      _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol ( 0, 1, el_ndof[1], el_conn, offset, 0, p_1ts ); // dp pressure
      _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol ( 0, 1, el_ndof[1], el_conn, offset, 0, p_2ts ); // dp pressure

      for ( int dim = 0; dim < _nPdim; dim++ )
        for ( int node = 0; node < el_ndof[2]; node++ )
          u_div[node + dim * el_ndof[2]] = div_1ts * u_1ts[dim * NDOF_FEM + node] + div_2ts * u_2ts[dim * NDOF_FEM + node] + div_3ts * u_3ts[dim * NDOF_FEM + node];

      for ( int node = 0; node < el_ndof[1]; node++ )
        p_rhs[node] = pres_1ts * p_1ts[node] + pres_2ts * p_2ts[node];

      for ( int  j = 0; j < el_ndof[1]; j++ ) {
          _bc_el[j] = _bc_vol[j] / 10;
          }

            double linear_nodes[8*3];
  for(int dim=0; dim<3; dim++)
      for(int node=0; node<8; node++)
          linear_nodes[node + 8*dim] = _xx_qnds[node + dim*NDOF_FEM];

      // ===================================================================================================
      //                          C) Gaussian integration loop (n_gauss)
      // ===================================================================================================
      for ( int  qp = 0; qp < el_ngauss; qp++ ) {
          // ------------------------- Shape functions at gaussian points --------------------------------------
      const double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // quadratic Jacobian
      JxW_g = det2 * _fe[2]->_weight1[ _nPdim - 1][qp]; // quadratic weight
      _fe[2]->get_phi_gl_g ( _nPdim, qp, _phi_g[2] );                // quadratic shape function
      _fe[2]->get_dphi_gl_g ( _nPdim, qp, _InvJac2, _dphi_g[2] );    // global coord deriv
      
      const double det1      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac1 ); // quadratic Jacobian
      _fe[1]->get_dphi_gl_g ( _nPdim, qp, _InvJac1, _dphi_g[1] );    // global coord deriv
      _fe[1]->get_phi_gl_g ( _nPdim, qp, _phi_g[1] );    // global coord deriv

          double dtxJxW_g = JxW_g;

          interp_el_gdx ( u_div, 0, _nPdim, _dphi_g[2], el_ndof[2], vel_gdx ); // derivatives  vel_gdx[DIM][DIM]
          interp_el_gdx ( p_rhs, 0, 1, _dphi_g[1], el_ndof[1], p_gdx ); // derivatives  vel_gdx[DIM][DIM]

          double Div_g = 0.;

          for ( int dim = 0; dim < _nPdim; dim++ )
            Div_g += vel_gdx[dim + dim * DIMENSION];


          // ===================================================================================================
          //                          D) Local (element) assemblying pressure equation
          // ===================================================================================================
          for ( int  i = 0; i < el_ndof[1]; i++ )     {
              const double phii_g = _phi_g[1][i];

              if ( _bc_el[i] == 1 ) {

                  for ( int dim = 0; dim < _nPdim; dim++ )
                    for ( int ll = 0; ll < NDOF_FEM; ll++ )
                      FeM ( i ) -= dtxJxW_g * u_div[dim * NDOF_FEM + ll] * _dphi_g[2][ll + dim * el_ndof[2]] * phii_g / _dt;

                  for ( int dim = 0; dim < _nPdim; dim++ )
                    FeM ( i ) += dtxJxW_g * p_gdx[dim] * _dphi_g[1][i + dim * el_ndof[1]] ;

                  // ------------------------- Matrix Assemblying  ----------------------------
                  for ( int  j = 0; j < el_ndof[1]; j++ ) {
                      const double phij_g = _phi_g[1][j];
                      double Lap = 0.;

                      for ( int  idim = 0; idim < _nPdim; idim++ ) {
                          if ( axis == 1 ) {
                              Lap += ( 1 - ( idim ) ) * phii_g * phij_g / ( _xyzg[0] * _xyzg[0] ); // axysimmetry
                              }

                          Lap += _dphi_g[1][j + idim * el_ndof[1]] * _dphi_g[1][i + idim * el_ndof[1]]; // Laplacian
                          }

                      KeM ( i, j ) += dtxJxW_g * Lap;
                      } // ---------------------------------------------
                  }
              else {
                  if ( _bc_vol[i] == 4 || _bc_vol[i] == 0 )  {
                      FeM ( i ) += dtxJxW_g * p_1ts[i];
                      }

                  KeM ( i, i ) += dtxJxW_g;
                  }
              }
          } // end of the quadrature point qp-loop +++++++++++++++++++++++++

      A[Level]->add_matrix ( KeM, el_dof_indices );              // global matrix

      if ( mode == 1 ) {
          b[Level]->add_vector ( FeM, el_dof_indices ); // global rhs
          }
      } // end of element loop

  // clean and close
  el_dof_indices.clear();
  A[Level]->close();

  if ( mode == 1 ) {
      b[Level]->close();
      }

#ifdef PRINT_INFO
  std::cout << " Matrix Assembled(P)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif
  return;
  } /******************************************************************************************************/
//
//

void  MGSolP::GenRhs (
  const double /* time*/, // time  <-
  const int  Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
)    // ===============================================
  {
  int axis = ( int ) ( _mgutils._geometry["Axisym"] );
  const int   offset = _mgmesh._NoNodes[_NoLevels - 1];                               // mesh nodes
  // element connectivity
  const int  el_sides = _mgmesh._GeomEl._n_sides[0];                         // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                          // element connectivity
  int        el_neigh[NDOF_FEM];                                              // element connectivity

  /*----------------------------------- Gauss integration  --------------------------------------------*/

  double det, JxW_g;                          // Jac, Jac*w Jacobean
  const int  el_ngauss = _fe[2]->_NoGauss1[ _nPdim - 1];                              //elem gauss points

  double dphijdx_g[3][DIMENSION];
  double dphiidx_g[3][DIMENSION];             // global derivatives at g point

  double u_1ts[DIMENSION * NDOF_FEM], u_2ts[DIMENSION * NDOF_FEM], u_3ts[DIMENSION * NDOF_FEM], u_div[DIMENSION * NDOF_FEM];
  double p_1ts[NDOF_P], p_2ts[NDOF_P], p_rhs[NDOF_P];

  double vel_gdx[DIMENSION * DIMENSION], p_gdx[DIMENSION];

  // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
  int el_ndof[3];
  el_ndof[0] = 1;                // number of el dofs
  int el_mat_nrows = 0;                           // number of mat rows (dofs)

  for ( int ideg = 1; ideg < 3; ideg++ ) {
      el_ndof[ideg] = _fe[ideg]->_NoShape[ _nPdim - 1];
      el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
      };

  int el_mat_ncols = el_mat_nrows;                     //square matrix

  std::vector<int > el_dof_indices ( el_mat_ncols );   // element dof vector

  /*---------------------------- Fields -> Navier-Stokes  [NS_F]  -------------------------------------*/
  for ( int k = 0; k < 30; k++ ) { // coupling  basic system fields
      const int idx = _data_eq[2].tab_eqs[k];
      _FF_idx[k] = ( idx >= 0 ) ? _data_eq[2].indx_ub[idx] : -1;
      }

  /* --------------------------- Element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --*/
  if ( mode == 1 ) {
      b[Level]->zero();    // global matrix+rhs
      }

  DenseVectorM FeM;                // local  matrix+rhs
  FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

  int ndof_lev = 0;

  for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
      int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
      ndof_lev += delta;
      }

  const int  nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
  const int  nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

  double div_1ts, div_2ts, div_3ts;
  double pres_1ts, pres_2ts;

  if ( _P_parameter._TimeDisc == 2 ) {
      div_1ts = 3. / 2.;
      div_2ts = -2.;
      div_3ts = 1. / 2.;
      pres_1ts = 2.;
      pres_2ts = -1.;
      }

  if ( _P_parameter._TimeDisc == 1 ) {
      div_1ts = 1.;
      div_2ts = 0.;
      div_3ts = 0.;
      pres_1ts = 1.;
      pres_2ts = 0.;
      }
      
  for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {

      // ===================================================================================================
      //                          B) Element  Loop over the volume (n_elem)
      // ===================================================================================================

      // ------------------------ Set to zero matrix and rhs and center -----------------------------------
      FeM.zero();
      // ------------------------ Geometry and element fields ---------------------------------------------
      // ------------------------ Element Connectivity (el_conn) and coordinates (xx_qnds) ----------------
      _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, _xx_qnds );
      _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh );
      get_el_dof_bc ( Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

      for ( int idim = 0; idim < _nPdim; idim++ ) {
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_1ts );
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_oldsol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_2ts );
          _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_oooldsol ( 0, 1, el_ndof[2], el_conn, offset, idim, u_3ts );
          }

      _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol ( 0, 1, el_ndof[1], el_conn, offset, 0, p_1ts ); // dp pressure
      _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol ( 0, 1, el_ndof[1], el_conn, offset, 0, p_2ts ); // dp pressure

      for ( int dim = 0; dim < _nPdim; dim++ )
        for ( int node = 0; node < el_ndof[2]; node++ )
          u_div[node + dim * el_ndof[2]] =   div_1ts * u_1ts[dim * NDOF_FEM + node]
                                           + div_2ts * u_2ts[dim * NDOF_FEM + node] 
                                           + div_3ts * u_3ts[dim * NDOF_FEM + node];

      for ( int node = 0; node < el_ndof[1]; node++ )
        p_rhs[node] = pres_1ts * p_1ts[node]
                    + pres_2ts * p_2ts[node];

      for ( int  j = 0; j < el_ndof[1]; j++ ) {
          _bc_el[j] = _bc_vol[j] / 10;
          }
          
      // ===================================================================================================
      //                          C) Gaussian integration loop (n_gauss)
      // ===================================================================================================
      for ( int  qp = 0; qp < el_ngauss; qp++ ) {
          // ------------------------- Shape functions at gaussian points --------------------------------------
      const double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // quadratic Jacobian
      JxW_g = det2 * _fe[2]->_weight1[ _nPdim - 1][qp]; // quadratic weight
      _fe[2]->get_phi_gl_g ( _nPdim, qp, _phi_g[2] );                // quadratic shape function
      _fe[2]->get_dphi_gl_g ( _nPdim, qp, _InvJac2, _dphi_g[2] );    // global coord deriv
      
      const double det1      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac1 ); // quadratic Jacobian
      _fe[1]->get_dphi_gl_g ( _nPdim, qp, _InvJac1, _dphi_g[1] );    // global coord deriv
      _fe[1]->get_phi_gl_g ( _nPdim, qp, _phi_g[1] );    // global coord deriv

          double dtxJxW_g = JxW_g;

          interp_el_gdx ( u_div, 0, _nPdim, _dphi_g[2], el_ndof[2], vel_gdx ); // derivatives  vel_gdx[DIM][DIM]
          interp_el_gdx ( p_rhs, 0, 1, _dphi_g[1], el_ndof[1], p_gdx ); // derivatives  vel_gdx[DIM][DIM]


          // ===================================================================================================
          //                          D) Local (element) assemblying pressure equation
          // ===================================================================================================
          for ( int  i = 0; i < el_ndof[1]; i++ )     {
              const double phii_g = _phi_g[1][i];

              if ( _bc_el[i] == 1 ) {

                  for ( int dim = 0; dim < _nPdim; dim++ )
                    for ( int ll = 0; ll < NDOF_FEM; ll++ )
                      FeM ( i ) -= dtxJxW_g * u_div[dim * NDOF_FEM + ll] * _dphi_g[2][ll + dim * el_ndof[2]] * phii_g / _dt;

                  for ( int dim = 0; dim < _nPdim; dim++ )
                    FeM ( i ) += dtxJxW_g * p_gdx[dim] * _dphi_g[1][i + dim * el_ndof[1]] ;
                  }
              else {
                  if ( _bc_vol[i] == 4 || _bc_vol[i] == 0 )  {
                      FeM ( i ) += dtxJxW_g * p_1ts[i];
                      }
                  }
              }
          } // end of the quadrature point qp-loop +++++++++++++++++++++++++

      if ( mode == 1 ) {
          b[Level]->add_vector ( FeM, el_dof_indices ); // global rhs
          }
      } // end of element loop

  // clean and close
  el_dof_indices.clear();

  if ( mode == 1 ) {
      b[Level]->close();
      }

#ifdef PRINT_INFO
  std::cout << " RHS Assembled(P)  for  Level " << Level << "\n";
#endif
  return;
  } /******************************************************************************************************/
//
//


// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolP::MGTimeStep (
  const double time,  // time
  const int iter  // Number of max inter
)
  {
  MGTimeStep_no_up ( time, iter );
  MGUpdateStep();

  return;
  } /******************************************************************************************************/

// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolP::MGTimeStep_no_up (
  const double time,  // time
  const int iter  // Number of max inter
)
  {
  /* ========================================================================================= */
  /*              A) Set up the time step                                                      */
  /* ========================================================================================= */
  std::cout  << std::endl << "\033[038;5;" << 155 << ";1m "
             << "--------------------------------------------------- \n\t"
             <<  _eqname.c_str()
             << " solution of problem " << _mgutils.get_name()
             << "\n ---------------------------------------------------\n \033[0m";
  /* ========================================================================================= */
  /*              B) Assemblying of the Matrix-Rhs                                             */
  /* ========================================================================================= */
#if PRINT_TIME==1
  std::clock_t start_time = std::clock();
#endif


  if ( _AssembleOnce == 0 || _AlreadyAssembled == 0 ) {
      GenMatRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs

      for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
          GenMatRhs ( time, Level, 0 ); // matrix
          }

      if ( _AssembleOnce == 1 )
        _AlreadyAssembled = 1;
      }

  if ( _AssembleOnce == 1 && _AlreadyAssembled == 1 ) {
      GenRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs
      }

#if PRINT_TIME==1
  std::clock_t end_time = std::clock();
  std::cout << "  Assembly time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
  /* ========================================================================================= */
  /*               C) Solution of the linear MGsystem (MGSolP::MGSolve)                        */
  /* ========================================================================================= */
  MGSolve ( 1.e-6, 15 );
#if PRINT_TIME==1
  end_time = std::clock();
  std::cout << " Assembly+solution time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC
            << "s " << std::endl;
#endif

  return;
  } /******************************************************************************************************/

// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolP::MGUpdateStep()
  {

  x_old[_NoLevels - 1]->localize ( *x_oold[_NoLevels - 1] ); // p(n-1)
  x[_NoLevels - 1]->localize ( *x_old[_NoLevels - 1] );      // p(n)

  return;
  } /******************************************************************************************************/

#endif  //ENDIF NS_EQUATIONS%2==0
#endif  //ENDIF NS_EQUATIONS
// #endif  // NS_equation is personal



