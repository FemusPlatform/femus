// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include "MGSolverNS_coup.h"       // Navier-Stokes class header file
#include "MGFE_conf.h"        // FEM approximation
#include "MGFE.h"          // Mesh class


MGSolNS_coup::MGSolNS_coup(
  MGEquationsSystem & mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
) :  
MGSolNS ( mg_equations_map_in, nvars_in, eqname_in, varname_in ){
    _Coupled = 1;    
    SetVariableNames ( varname_in );
}


void MGSolNS_coup::SetVariableNames ( std::string varname_in ) {

  _dir=0;  
    
  _var_names[_nNSdim - 1] = "w";
  _refvalue[_nNSdim - 1] = _uref; // velocity 3D
  _var_names[0] = "u";
  _refvalue[0] = _uref; // velocity 2D
  _var_names[1] = "v";
  _refvalue[1] = _uref; // velocity 2D

  _var_names[_nNSdim] = "p";
  _refvalue[_nNSdim] = _rhof * _uref * _uref; // pressure

  return;
  }

void MGSolNS_coup::get_el_field_data (
  int iel, int Level,
  int el_conn [], int offset, int el_ndof[], int ndof_lev ) {
  int el_ndofp = el_ndof[1];

  for ( int deg = 0; deg < 3; deg++ )
    for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                               el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
        }

  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_nonl_sol ( 0, _nNSdim, el_ndof[2], el_conn, offset, 0, _u_nl );
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol ( 0, _nNSdim, el_ndof[2], el_conn, offset, 0, _u_1ts );
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_oldsol ( 0, _nNSdim, el_ndof[2], el_conn, offset, 0, _u_2ts );

  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol ( _nNSdim, 1, el_ndof[1], el_conn, offset, 0, _p_1ts ); // pressure
  _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_oldsol ( _nNSdim, 1, el_ndof[1], el_conn, offset, 0, _p_2ts ); // pressure

  for ( int i = 0; i < NDOF_P; i++ ) {
      if ( _NS_parameter._TimeDisc == 2 )
        _Pressure[i] = ( 2.*_p_1ts[i] - _p_2ts[i] );

      if ( _NS_parameter._TimeDisc == 1 )
        _Pressure[i] = ( _p_1ts[i] );
      }

  return;
  }



void MGSolNS_coup::ContinuityEquation ( double JxW_g2, int el_ndof[] ) {

  for ( int  ikl = 1; ikl < 2; ikl++ ) { // ikl=0 discontinuous ikl=1 continuous pressure
      for ( int   i = 0; i < el_ndof[ikl]; i++ ) {
          const int    indp      = i + el_ndof[2] * _nNSdim;
          const double psii_g    = _phi_g[ikl][i];
          const double dtxJxWp_g = ( _bc_el[indp] != 0 ) ? JxW_g2 : 0.;

          for ( int j = 0; j < el_ndof[2]; j++ )     {
              for ( int  jvar = 0; jvar < _nNSdim; jvar++ ) { // linear -quad
                  double Div = 0.;
                  CalcVelDivergence ( Div, ikl, i, j, jvar, el_ndof );

                  _KeM ( indp, j + jvar * el_ndof[2] ) +=  dtxJxWp_g * Div;
                  }// jvar
              }// j end linear-quad --------------------------------------------

              _KeM ( indp, i + _nNSdim * el_ndof[2] ) +=  dtxJxWp_g * 1.e-40;
              
          }// i
      }// ikl=0 discontinuous ikl=1 continuous pressure

  return;
  }



void MGSolNS_coup::MomentumEquation(double JxW_g2, int el_ndof[], int qp)
{
    double dphiidx_g2[DIMENSION], dphijdx_g2[DIMENSION];
    
    for ( int i = 0; i < el_ndof[2]; i++ ) { // LOOP OVER TEST FUNCTIONS --------------------------------------------------
          // set up test function phii_g, derivatives dphiidx_g2[dim], and supg _Phi_supg

          double phii_g = _phi_g[2][i];

          for ( int idim = 0; idim < _nNSdim; idim++ ) {
              dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof[2]];
              }

          for ( int  row_shift = 0; row_shift < _nvars[2]; row_shift++ )    { // LOOP OVER ROWS: i + row_shift*el_ndof2
              int indx     = i + row_shift * el_ndof[2];
              RowSetUp ( i, indx, qp, el_ndof )  ;


              // zero line for Dirichlet  bc -----------------------------------------------------------
              if ( _bc_el[indx] < 0 || ( _bc_el[indx] == 0 && ( _bc_vol[i] == 11 || _bc_vol[i] == 31 ) ) ) {
                  for ( int  ivarN = row_shift; ivarN < row_shift + _NComp; ivarN++ )     { // Loop over velocity components for equation written in row indx
                      const int ivar = ivarN % _nNSdim;

                      // computation of alpha beta (factor for adding equations) ----------------------------------------------
                      double alpha_beta = 1.;
                      if ( _BoundEquation ) alpha_beta = _ProjDir[abs ( _bc_el[indx] ) - 1][ivar];

                      double  dtxJxW_g = JxW_g2 * alpha_beta; // (factor for adding equations)

                      // -------------------------------------------------------------------------------------------------------------
                      // -------------------------------------- Assemblying rhs ------------------------------------------------------

                          if ( _NS_parameter._SolveSteady == 0 ) {
                              _FeM ( indx )  +=  dtxJxW_g *  _u_1ts_g[ivar] * ( phii_g + _Phi_supg ) / _dt;
                              }

                          _FeM ( indx )  +=  dtxJxW_g *  _IFr * _dirg[ivar] * ( phii_g + _Phi_supg );


                      // -------------------------------------------------------------------------------------------------------------
                      //-------------------------------------- Assemblying matrix ----------------------------------------------------
                      for ( int j = 0; j < el_ndof[2]; j++ ) {
                          const double phij_g = _phi_g[2][j];
                          double Lap_g = 0., Adv_g = 0., Supg_lap = 0., Turb_grad = 0., up = 0.;

                          if ( _AxiSym == 1 ) { // axysimmetry only --------------------------------------------------
                              Lap_g = 2.* ( 1 - ( ivar + _dir ) ) *
                                      ( _IReEff + _NS_parameter._Upwind * _f_upwind * _u_1ts_g[0] * _u_1ts_g[0] )
                                      * phij_g * ( phii_g + _Phi_supg ) / ( _xyzg[0] * _xyzg[0] );
                              } // --------------------------------------------------------------------------------

                          for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
                              dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof[2]];
                              up = _NS_parameter._Upwind * _f_upwind * _u_1ts_g[kdim] * _u_1ts_g[kdim];
                              Adv_g     += _u_1ts_g[kdim] * dphijdx_g2[kdim];                                                                // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
                              Lap_g     += ( _IReEff + up ) * dphijdx_g2[kdim] * dphiidx_g2[kdim];
                              Supg_lap  += ( _IReEff ) * _ddphi_g[2][j * _nNSdim * _nNSdim + kdim * _nNSdim + kdim];
                              }

                          //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
                          if ( _NS_parameter._SolveSteady == 0 ) {
                              _KeM ( indx, j + ivar * el_ndof[2] ) += dtxJxW_g *  ( phii_g + _Phi_supg ) * phij_g / _dt;
                              }

                          _KeM ( indx, j + ivar * el_ndof[2] ) += dtxJxW_g  * (
                                                                 Adv_g * ( phii_g + _Phi_supg ) // time + advection
                                                                 + Lap_g                       // viscous Laplacian
                                                                 - Supg_lap * _Phi_supg          // SUPG viscous Laplacian
                                                                 - Turb_grad * _Phi_supg
                                                               );

                          // non Diagonal blocks
                          for ( int  ivars = 0; ivars < _nvars[2]; ivars++ ) {
                              double part = dtxJxW_g *  _dphi_g[2][j + ivar * el_ndof[2]];
                              _KeM ( indx, j + ivars * el_ndof[2] )   += _IReEff * part * dphiidx_g2[ivars]; //
                              }
                          } // end A element matrix quad -quad (end loop on j)---------------------------------------------

                      // ------------------------------------------------------------------
                      // B^T element matrix ( p*div(v) )--------------------
                      for ( int  jp = 0; jp < el_ndof[1]; jp++ ) {
                          const double psij_g = _phi_g[1][jp];

                          //  B^T  axisimmetric term
                          if ( _AxiSym == 1 ) {
                              _KeM ( indx, jp + _nNSdim * el_ndof[2] ) -= dtxJxW_g * ( 1 - ( ivar) ) * psij_g * phii_g / _xyzg[0];
                              }

                          //  B^T diagonal term
                          _KeM ( indx, jp + _nNSdim * el_ndof[2] ) += dtxJxW_g * ( // MPascal
                                -1.*_phi_g[1][jp] * _dphi_g[2][i + ivar * el_ndof[2]]
                                +  _dphi_g[1][jp + ivar * el_ndof[1]] * ( 0.*phii_g + _Phi_supg )
                              );
                          //  B^T non-diagonal term
                          } // jp pressure linear
                      }
                  } // end loop ivar
              }
          } // end loop i
    
    return;
}




void MGSolNS_coup::set_bc_matrix (
  int dir_maxnormal,                          ///<  normal dir
  int sur_toply[],                            ///< boundary topology map
  int el_ndof[],                              ///< number of volume dofs
  int elb_ndof[],                             ///< number of boundary dofs
  int elb_ngauss,                             ///<  number of surface gaussian points
  double normal[],                            ///< normal
  int el_conn[]
) {
  const double delta = _Wall_dist;
  double pressure;
  double vel_bound[1];
  double xyz_g[DIMENSION];
  const int elb_dof0 = ( NDOF_K > 1 ) ? elb_ndof[_pres_order] * NDOF_K : elb_ndof[_pres_order];
  double det2 = _fe[2]->JacSur ( elb_ngauss - 1, _xxb_qnds, _InvJac2 ); // jacobian
  double Ipenalty = det2;                          // Dirichlet bc flagdouble xyz_bg[DIMENSION];
  Ipenalty = 1.e0;

  const int bd_face = abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 100;
  const int norm_face = ( abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 10000 ) / 1000 - 1 ;
  
  int bd_node, bc_n_flag, bc_tg_flag, bc_var_normal;
  for (int dir = 0; dir < _nNSdim; dir++) {
        _ProjDir[2][dir] = normal[dir];
    }
  
  // TANGENT VECTORS BASED ON BOUNDARY GEOMETRY
  CalcTangDir ( _ProjDir[0], _ProjDir[1], _ProjDir[2], GEOM_BASED );
  SetBCFlags ( _ProjDir[2], sur_toply, el_ndof[2], el_conn );

  // ...........................................................................................
  //     DIRICHLET BOUNDARY CONDITION
  // ...........................................................................................

  for ( int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++ ) {

      NodeBCFlags ( sur_toply[lbnode], bd_node, bc_tg_flag, bc_n_flag, bc_var_normal );
      
      if ( _nNSdim == 3 ) {
          CorrectBCFlags ( _ProjDir[0], _ProjDir[1], bc_var_normal, sur_toply[lbnode], el_ndof[2], 0 );
          }

      if ( bd_face == bd_node ) {
          for ( int row_shift = 0; row_shift < _nvars[2]; row_shift++ ) { // LOOP OVER SPACE DIRECTIONS =======================================
              const int  indx_row = sur_toply[lbnode] + row_shift * el_ndof[2]; // solution indx

              bool NormBC = ( abs ( _bc_el[indx_row] ) == 3 ) ? true : false;
              bool DirichletBC = ( _bc_el[indx_row] > 0 ) ? true : false;

              if ( DirichletBC /*!_AlreadyWrittenDirBC[indx_row]*/ ) {
//           _AlreadyWrittenDirBC[indx_row] = true;
                  int  bc_bc  = NormBC ? bc_n_flag : bc_tg_flag;
                  const int  bc_rhs = ( ( bc_bc % 2 ) == 1 ) ? 1 : 0; // bc_rhs -> rhs bc
                  bool OrDir = ( bc_bc > 7 ) ? true : false;
                  const double vart = bc_rhs * _u_1ts[indx_row];

                  double alpha_beta = 1.;

                  if ( bd_node == 88 || bd_node == 66 ) {
                      _KeM ( indx_row, indx_row ) = 1.;
                      _FeM ( indx_row )          = 0.;
                      }
                  else if ( OrDir ) {
                      for ( int   jvar = 0; jvar < _nNSdim; jvar++ )    { //  -beta_n*n*( u.n)
                          const int indj = sur_toply[lbnode] + jvar * el_ndof[2];
                          alpha_beta = -normal[jvar];

                          if ( !NormBC ) {
                              alpha_beta = ( _bc_el[indx_row] == 1 ) ? _ProjDir[0][jvar] : _ProjDir[1][jvar];
                              }

                          _KeM ( indx_row, indj ) += alpha_beta;
                          _FeM ( indx_row ) += _u_1ts[indj] * bc_rhs * alpha_beta;
                          }
                      }
                  else {
                      _KeM ( indx_row, indx_row ) = alpha_beta;
                      _FeM ( indx_row ) = vart * alpha_beta;
                      }// beta_0 u -> diagonal
                  }
              }// END DIRICHLET BLOCK ====================================================================================
          }
      }


  // ...........................................................................................
  //     PRESSURE BOUNDARY CONDITION
  // ...........................................................................................

  const int bc_side_check2 = _bc_bd[sur_toply[elb_ndof[2] - 1]]; // total  bc_var
  const int bc_side_normal = ( bc_side_check2 % 10000 ) / 1000 - 1; //
  const int bc_side_n_flag = ( ( bc_side_check2 % 10000 ) % 1000 ) / 10;

  if ( bc_side_n_flag < 4 && bc_side_n_flag > 1 ) {
      const int  rhs_p = bc_side_n_flag % 2;

      for ( int i = 0; i < elb_ndof[1]; i++ ) {
          const int   indx = _pres_order * sur_toply[i] + _nNSdim * NDOF_FEM; // volume dof index
          _bc_el[indx] = 0;                                                                                    // SET TO 1 -> CONTROL IF RIGHT
          _KeM ( indx, indx )  += Ipenalty;
          _FeM ( indx )  +=  rhs_p * 0.*Ipenalty * _data_eq[_pres_order].ub[_pres_order * sur_toply[i]];
//       std::cout<<_data_eq[_pres_order].ub[_pres_order*sur_toply[i]]<<std::endl;
          }
      }

  // ...........................................................................................
  //     INTEGRATION OVER BOUNDARY ELEMENT - STRESS CALCULATION
  // ...........................................................................................

  double VelOnGauss[DIMENSION], u_tau, y_plus, eff_stress, VelOnBound, veltg;
  CalcTangDir ( _ProjDir[0], _ProjDir[1], _ProjDir[2], VEL_BASED );

  for ( int  qp = 0; qp < elb_ngauss; qp++ ) { // START LOOP OVER GAUSS POINTS ===============================================
      double det   = _fe[2]->JacSur ( qp, _xxb_qnds, _InvJac2 ); // local coord _phi_g and jac
      double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp]; // weight
      _fe[2]->get_phi_gl_g ( _nNSdim - 1, qp, _phi_g[2] ); // global coord _phi_g
      _fe[1]->get_phi_gl_g ( _nNSdim - 1, qp, _phi_g[1] ); // global coord _phi_g

      if ( _AxiSym == 1 ) {
          interp_el_bd_sol ( _xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg );
          JxW_g  *= _xyzg[0];   // axisymmetric  (index ->0)
          }

      interp_el_bd_sol ( _u_1ts, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], VelOnGauss );


      for ( int i = 0; i < elb_ndof[2]; i++ ) { // LOOP OVER TEST FUNCTIONS ===================================================

          NodeBCFlags ( sur_toply[i], bd_node, bc_tg_flag, bc_n_flag, bc_var_normal );
          const int bc_varn_flag = bd_node / 10;
          const int bc_vartg_flag = bd_node % 10;

          double phii_g  = _phi_g[2][i];

          if ( bd_face == bd_node ) {


              const int pen_cond = ( bd_node == 44 ) ? 0 : 1;
              double lambda1 = ( 1 - pen_cond ) * _NS_parameter._Penalty_tg;
              double beta = lambda1 * _NS_parameter._Penalty_n;

              // Direction of maximum vector component of tangent 1 and 2 -> alignment of equations
              if ( _nNSdim == 3 ) {
                  CorrectBCFlags ( _ProjDir[0], _ProjDir[1], bc_var_normal, sur_toply[i], el_ndof[2], 1 );
                  }

              for ( int  row_shift = 0; row_shift < _nvars[2]; row_shift++ ) { // LOOP OVER ROWS: i + row_shift*el_ndof2
                  int indx     = sur_toply[i] + row_shift * el_ndof[2];
                  int NComp = ( _bc_el[indx] < 0 ) ? _nNSdim : 1;
                  int bc_el = ( _bc_el[indx] <= 0 ) ? 1 : 0;
                  const double dtJxW_g = bc_el * JxW_g;

                  bool NormBC = ( _bc_el[indx] == -3 ) ? true : false;
                  const int bc_flag = ( NormBC ) ? bd_node / 10 : bd_node % 10;

                  if ( bc_flag < 6 && bd_node != 11 ) { // only for neumann bc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      // NORMAL VECTOR FOR THE PROJECTION OF VOLUME PART EQUATION
                      // _normal_pt is used to project the equation (volume part) using the same normal used in the boundary
                      if ( qp == 0 ) for ( int dir = 0; dir < _nNSdim; dir++ ) {
                            _normal_pt[ sur_toply[i]*_nNSdim + dir] += normal[dir];
                            }

                      for ( int  ivarN = row_shift; ivarN < row_shift + NComp; ivarN++ ) { // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
                          const int ivar = ivarN % DIMENSION;

                          // Stress component along tangential directions (-1,-2) ---------------------------------------------------
                          if ( _bc_el[indx] == -1 || _bc_el[indx] == -2 ) {
                              if ( bc_flag != 1 ) {
                                  const int  bc_rhs = ( ( bc_vartg_flag % 2 ) == 1 ) ? 1 : 0;
                                  double alpha_beta = ( _bc_el[indx] == -1 ) ? _ProjDir[0][ivar]: _ProjDir[1][ivar];

                                  if ( bc_flag == 5 ) {
                                      alpha_beta = _ProjDir[0][ivar];
                                      }

                                  double tau = ( _bc_el[indx] == -1 ) ? ( _NS_parameter._Tg1_stress ) : _NS_parameter._Tg2_stress ; //utau*utau/vel_bound;
                                  double stress = alpha_beta * ( pen_cond * tau - lambda1 ) * dtJxW_g * _phi_g[2][i];

//                 FeM (indx) += /*bc_rhs**/ stress*VelOnGauss[ivar];
                                  for ( int j = 0; j < elb_ndof[2]; j++ ) {
                                      _KeM ( indx, sur_toply[j] + ivar * el_ndof[2] ) += -1.*stress * _phi_g[2][j] ; /// ( modulus[j] );
                                      }// END LOOP OVER MATRIX COMPONENTS - IMPLICIT STRESS
                                  }
                              }// End loop if equation is for tangential direction-----------------------------------------------------
                          else {
                              // Stress component along normal direction(-3,0) -----------------------------------------------------
                              double alpha_beta = 1.;

                              if ( _bc_el[indx] == -3 ) {
                                  alpha_beta = normal[ivar];

                                  // normal with penalty (u dot n) n dot phi
                                  if ( bc_varn_flag == 4 )
                                    for ( int j = 0; j < elb_ndof[2]; j++ )
                                      for ( int jvar = 0; jvar < _nNSdim; jvar++ ) {
                                          _KeM ( indx, sur_toply[j] + jvar * el_ndof[2] ) -= alpha_beta * normal[ivar] * normal[jvar] * beta * dtJxW_g * _phi_g[2][j] * _phi_g[2][i];
                                          }
                                  } //  ---------------------------------------------------------------------------------------------------

                              // surface integral p phi dot n   --> pressure_inlet - pressure_outlet - outflow_p ------------------------------------------------
                              if ( bc_varn_flag > 1 && bc_varn_flag < 4 ) {
                                  for ( int j = 0; j < elb_dof0; j++ ) {
                                      double det = dtJxW_g * alpha_beta * normal[ivar];
                                      _KeM ( indx, sur_toply[j] + _nNSdim * el_ndof[2] ) += det * _phi_g[2][i] * _phi_g[1][j];
//                       FeM(indx)  -= det*_phi_g[2][i]*_phi_g[1][j]*_data_eq[_pres_order].ub[(_FF_idx[NS_F]) *NDOF_FEM+_pres_order*sur_toply[j]];
                                      }
                                  } //  ---------------------------------------------------------------------------------------------------
                              }// contribution along normal direction

                          }// END LOOP IF BC FLAG IS FOR INTEGRATION
                      }// END LOOP OVER VARIABLES FOR LINEAR COMBINATION OF EQUATIONS ++++++++++++++++++++++++++++++++++++++++++++++
                  }// END LOOP OVER MATRIX ROWS RELATIVE TO TEST FUNCTION I
              }// END IF BC_NODE == BC_FACE
          }// END LOOP OVER BOUNDARY TEST FUNCTIONS =========================================================================
      }// END GAUSSIAN INTEGRATION ========================================================================================

  return;
  }

#endif









