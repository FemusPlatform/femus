// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
#if NS_EQUATIONS==2
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class


void MGSolNS::get_el_field_data (
  int iel, int Level,
  int el_conn [], int offset, int el_ndof[], int ndof_lev,
  double u_old[],  double u_oold[],   double u_nl[],
  double p_1ts[], double p_2ts[]
) {
  int el_ndofp = el_ndof[1];

  for ( int deg = 0; deg < 3; deg++ )
    for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                               el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
        }

  for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_nonl_sol ( 0, 1, el_ndof[2], el_conn, offset, kdim, u_nl );
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, kdim, u_old );
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_oldsol ( 0, 1, el_ndof[2], el_conn, offset, kdim, u_oold );
      }

  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol ( 0, 1, el_ndofp, el_conn, offset, 0, p_1ts );       // time step -1
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol ( 0, 1, el_ndofp, el_conn, offset, 0, p_2ts );   // time step -2


  return;
  }

// ==============================================================================================
void MGSolNS::matrixrhsvol (
  DenseMatrixM & KeM, DenseVectorM & FeM,
  int el_ndof[],
  double u_1ts[],
  double u_2ts[],
  double u_nl[],
  const int unsteady_flag,
  const int mode,
  int el_conn[],
  double p_1ts[],
  double p_2ts[]
) {
  double dphijdx_g2[DIMENSION], dphiidx_g2[DIMENSION];

  double vel_gddx[DIMENSION * DIMENSION * DIMENSION];
  double vel_gdx[DIMENSION * DIMENSION];

  double u_1ts_g[DIMENSION], u_2ts_g[DIMENSION];

  double rho = 1.;
  const int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim - 1];              //elem gauss points
  const int el_ndof2 = el_ndof[2];

  double ProjDir[3][DIMENSION];
  double Pressure[NDOF_P],   P_gdx[DIMENSION];

  for ( int i = 0; i < NDOF_P; i++ ) {
      if ( _NS_parameter._TimeDisc == 2 )
        Pressure[i] = ( 2.*p_1ts[i] - p_2ts[i] );

      if ( _NS_parameter._TimeDisc == 1 )
        Pressure[i] = ( p_1ts[i] );

      }

  double linear_nodes[8 * 3];

  for ( int dim = 0; dim < 3; dim++ )
    for ( int node = 0; node < 8; node++ )
      linear_nodes[node + 8 * dim] = _xx_qnds[node + dim * NDOF_FEM];


  for ( int qp = 0; qp <  el_ngauss; qp++ ) {
      // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
      // quadratic continuous (2)  (velocity)
      const double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // quadratic Jacobian
      double JxW_g2 = det2 * _fe[2]->_weight1[ _nNSdim - 1][qp]; // quadratic weight
      _fe[2]->get_phi_gl_g ( _nNSdim, qp, _phi_g[2] );                // quadratic shape function
      _fe[2]->get_dphi_gl_g ( _nNSdim, qp, _InvJac2, _dphi_g[2] );    // global coord deriv
      _fe[2]->get_ddphi_gl_g ( _nNSdim, qp, _InvJac2, _ddphi_g[2] );  // global second derivatives

      const double det1      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac1 ); // quadratic Jacobian
      _fe[1]->get_dphi_gl_g ( _nNSdim, qp, _InvJac1, _dphi_g[1] );    // global coord deriv

      // discontinuous (0) (disc pressure)
      if ( _nvars[0] > 0 ) _fe[0]->get_phi_gl_g ( _nNSdim, qp, _phi_g[0] ); // piecewise shape function

      // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------
      interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2] ); // field _ub_g[2][DIM]
      interp_el_sol ( u_2ts, 0, _nNSdim, _phi_g[2], el_ndof2, u_2ts_g ); // field _ub_g[2][DIM]
      interp_el_gdx ( u_1ts, 0, _nNSdim, _dphi_g[2], el_ndof2, vel_gdx ); // derivatives  vel_gdx[DIM][DIM]
      interp_el_gdx ( Pressure, 0, 1, _dphi_g[1], el_ndof[1], P_gdx );    // derivatives  vel_gdx[DIM][DIM]
      interp_el_gddx ( u_1ts, 0, _nNSdim, _ddphi_g[2], el_ndof2, vel_gddx );



      if ( _AxiSym == 1 ) {
          interp_el_sol ( _xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof2, _xyzg );
          JxW_g2  *= _xyzg[0];
          }

      for ( int idim = 0; idim <  _nNSdim; idim++ ) {
          u_1ts_g[idim] = _ub_g[2][_FF_idx[NS_F] + idim];
          }

      double IRe_eff = _IRe;

      if ( _FF_idx[MU_T] >= 0 ) {
          double muturb = _ub_g[2][_FF_idx[MU_T]];
          IRe_eff   += _IRe * muturb;
          }

      double f_upwind = CalcFUpwind ( u_1ts_g, _dphi_g[2], IRe_eff, _nNSdim, el_ndof2 );

      // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------

      // ===========================================================================================================
      //                                       D) Assembling NS equation
      // ===========================================================================================================

      for ( int i = 0; i < el_ndof2; i++ ) { // LOOP OVER TEST FUNCTIONS -------------------------------------------------
          // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
          int indx_eq     = i + _dir * el_ndof2;
          double phii_g = _phi_g[2][i];

          for ( int idim = 0; idim < _nNSdim; idim++ ) dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof2];

          if ( _bc_el[indx_eq] == -8 ) {
              _bc_el[indx_eq] = 0;

              if ( qp == 0 )
                std::cout << "ATTENTION!!! el node " << i << " glob node " << el_conn[i] << " row " << indx_eq << " bc " << _bc_el[indx_eq] << "   " << _dir << std::endl;
              }

          double Phi_supg = 0.;

          if ( _bc_el[indx_eq] == 0 && _NS_parameter._Supg == 1 ) {
              for ( int idim = 0; idim < _nNSdim; idim++ ) { // Phi_supg = 0. on boundaries with Dirichlet bc
                  Phi_supg += f_upwind * u_1ts_g[idim] * dphiidx_g2[idim]; //f_upwind
                  }
              }

          bool BoundEquation = ( _bc_el[indx_eq] <= -1 ) ? true : false;
          bool NormalEquation = false;

          if ( BoundEquation && _bc_el[indx_eq] == -3 ) {
              NormalEquation = true;

              for ( int dim = 0; dim < _nNSdim; dim++ ) {
                  ProjDir[2][dim] = _normal_pt[dim + i * _nNSdim];
                  }
              }

          if ( BoundEquation && !NormalEquation ) {
              CalcTangDir ( i, el_ndof[2], _bc_el[indx_eq], qp, ProjDir[0], ProjDir[1], VEL_BASED );
              }

          int NComp = ( BoundEquation ) ? _nNSdim : 1;                            // Default: 1 vel comp for equation in row indx -> ivar = row_shift

          // zero line for Dirichlet  bc -----------------------------------------------------------
          if ( _bc_el[indx_eq] < 0 || ( _bc_el[indx_eq] == 0 && ( _bc_vol[i] == 11 || _bc_vol[i] == 31 ) ) ) {
              for ( int  ivarN = _dir; ivarN < _dir + NComp; ivarN++ )     { // Loop over velocity components for equation written in row indx
                  const int ivar = ivarN % _nNSdim;
                  // computation of alpha beta (factor for adding equations) ----------------------------------------------
                  double alpha_beta = 1.;

                  if ( BoundEquation ) alpha_beta = ProjDir[abs ( _bc_el[indx_eq] ) - 1][ivar];

                  double  dtxJxW_g = JxW_g2 * alpha_beta; // (factor for adding equations)

                  // Regularization supg term ------------------------------------------------------------------------------------
                  if ( _NS_parameter._SolveSteady == 0 ) {
                      if ( _NS_parameter._TimeDisc == 1 )
                        FeM ( i )  +=  dtxJxW_g * rho * ( u_1ts_g[ivar] ) * ( phii_g + Phi_supg ) / ( _dt );

                      if ( _NS_parameter._TimeDisc == 2 )
                        FeM ( i )  +=  dtxJxW_g * rho * ( 2.*u_1ts_g[ivar] - 0.5 * u_2ts_g[ivar] ) * ( phii_g + Phi_supg ) / ( _dt );

                      }

                  FeM ( i )  +=  dtxJxW_g * rho * (
                                   _IFr * _dirg[ivar]   // force
                                   - P_gdx[ivar]        // pressure
                                 ) * ( phii_g + Phi_supg );

                  for ( int  ot = ivar; ot < ivar + _nNSdim; ot++ ) { // out of diagonal tensor components -> better explicit
                      const int ivars = ot % _nNSdim;
                      FeM ( i )   -= IRe_eff * dtxJxW_g * rho * vel_gdx[ivars * _nNSdim + ivar] * dphiidx_g2[ivars]; //
                      }

                  for ( int j = 0; j < el_ndof2; j++ ) {
                      const double phij_g = _phi_g[2][j];
                      double Lap_g = 0., Adv_g = 0., Supg_lap = 0., up = 0.;

                      if ( _AxiSym == 1 ) { // axysimmetry only --------------------------------------------------
                          Lap_g = 2.* ( 1 - ( _dir ) ) *
                                  ( IRe_eff + _NS_parameter._Upwind * f_upwind * u_1ts_g[0] * u_1ts_g[0] )
                                  * phij_g * ( phii_g + Phi_supg ) / ( _xyzg[0] * _xyzg[0] );
                          } // --------------------------------------------------------------------------------

                      for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
                          dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof2];

                          if ( _WallElement != 1 ) {
                              up = _NS_parameter._Upwind * f_upwind * u_1ts_g[kdim] * u_1ts_g[kdim];
                              }

                          if ( _NS_parameter._TimeDisc == 1 )
                            Adv_g     += ( u_1ts_g[kdim] ) * dphijdx_g2[kdim];

                          if ( _NS_parameter._TimeDisc == 2 )
                            Adv_g     += ( 2.*u_1ts_g[kdim] - u_2ts_g[kdim] ) * dphijdx_g2[kdim];

                          Lap_g     += ( IRe_eff + up ) * dphijdx_g2[kdim] * dphiidx_g2[kdim];
                          Supg_lap  += ( IRe_eff ) * _ddphi_g[2][j * _nNSdim * _nNSdim + kdim * _nNSdim + kdim];
                          }

                      //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
                      double TimeDerivative;

                      if ( _NS_parameter._TimeDisc == 2 )
                        TimeDerivative = dtxJxW_g * ( 3. / 2. ) * rho * ( phii_g + Phi_supg ) * phij_g / ( _dt );

                      if ( _NS_parameter._TimeDisc == 1 )
                        TimeDerivative = dtxJxW_g *  rho * ( phii_g + Phi_supg ) * phij_g / ( _dt );

                      const double KeM_ij = dtxJxW_g * rho * (
                                              Adv_g * ( phii_g + Phi_supg )  // advection
                                              + Lap_g                        // viscous Laplacian
                                              - Supg_lap * Phi_supg          // SUPG viscous Laplacian
                                            );
                      // non Diagonal blocks

                      if ( ivarN == _dir ) {
                          if ( _NS_parameter._SolveSteady == 0 ) KeM ( i, j ) += TimeDerivative;

                          KeM ( i, j ) += KeM_ij;
                          }
                      else {
                          if ( _NS_parameter._SolveSteady == 0 ) FeM ( i )  -=  TimeDerivative * u_1ts[j + ivar * NDOF_FEM];

                          FeM ( i )   -= KeM_ij * u_1ts[j + ivar * NDOF_FEM];
                          }
                      } // end A element matrix quad -quad (end loop on j)---------------------------------------------
                  }
              } // end loop ivar
          } // end loop i
      }

  return;
  }



// ==============================================================================================
void MGSolNS::matrixrhsvol_corr_step (
  DenseMatrixM & KeM, DenseVectorM & FeM,
  int el_ndof[],
  double u_old[],
  double u_nl[],
  const int unsteady_flag,
  const int mode,
  int el_conn[],
  double p_proj[],
  double p_old[]
) {
  double rho = 1.;
  int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim - 1];              //elem gauss points
  const int el_ndof2 = el_ndof[2];

  double Pressure[NDOF_P], Pressure_old[NDOF_P], P_gdx[DIMENSION], u_nlg[DIMENSION], P_ong[1];

  for ( int i = 0; i < NDOF_P; i++ ) {
      Pressure[i] = p_proj[i];
      Pressure_old[i] = p_old[i];
      }

  double ** ProjDir;
  ProjDir = new double*[3];

  for ( int l = 0; l < 3; l++ ) {
      ProjDir[l] = new double [DIMENSION];
      }

  for ( int qp = 0; qp <  el_ngauss; qp++ ) {
      // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
      // quadratic continuous (2)  (velocity)
      const double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // quadratic Jacobian
      double JxW_g2 = det2 * _fe[2]->_weight1[ _nNSdim - 1][qp];      // quadratic weight
      _fe[2]->get_phi_gl_g ( _nNSdim, qp, _phi_g[2] );                // quadratic shape function
      _fe[1]->get_phi_gl_g ( _nNSdim, qp, _phi_g[1] );                // quadratic shape function
      _fe[1]->get_dphi_gl_g ( _nNSdim, qp, _InvJac2, _dphi_g[1] );    // global coord deriv
      _fe[2]->get_dphi_gl_g ( _nNSdim, qp, _InvJac2, _dphi_g[2] );    // global coord deriv

      // discontinuous (0) (disc pressure)
      if ( _nvars[0] > 0 ) {
          _fe[0]->get_phi_gl_g ( _nNSdim, qp, _phi_g[0] ); // piecewise shape function
          }

      interp_el_sol ( _xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof2, _xyzg );
      interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2] ); // field _ub_g[2][DIM]
      interp_el_sol ( u_nl, 0, _nNSdim, _phi_g[2], el_ndof2, u_nlg ); // field _ub_g[2][DIM]

      const double TimeDisc = 1.;
      double alpha_beta = 1.;

      for ( int i = 0; i < el_ndof2; i++ ) { // LOOP OVER TEST FUNCTIONS --------------------------------------------------
          double phii_g = _phi_g[2][i];
          int indx_eq   = i + _dir * el_ndof2;
          int BoundEquation = ( _bc_el[indx_eq] <= -1 ) ? 1 : 0;
          int NormalEquation = ( BoundEquation == 1 && _bc_el[indx_eq] == -3 ) ? 1 : 0;

          if ( BoundEquation == 1 && NormalEquation == 0 ) {
              CalcTangDir ( i, el_ndof[2], _bc_el[indx_eq], qp, ProjDir[0], ProjDir[1], VEL_BASED );
              }

          for ( int dd = 0; dd < _nNSdim; dd++ ) {
              ProjDir[2][dd] = _normal_pt[dd + i * _nNSdim];
              }

          int NComp = ( BoundEquation == 1 ) ? _nNSdim : 1;                            // Default: 1 vel comp for equation in row indx -> ivar = row_shift

          if ( _bc_el[indx_eq] <= 0 ) {
              for ( int  ivarN = _dir; ivarN < _dir + NComp; ivarN++ )     { // Loop over velocity components for equation written in row indx
                  const int ivar = ivarN % _nNSdim;

                  if ( BoundEquation == 1 ) {
                      alpha_beta = ProjDir[abs ( _bc_el[indx_eq] ) - 1][ivar];
                      }

                  double  J_Wg_Alpha_phiig = JxW_g2 * alpha_beta * rho * phii_g;

                  for ( int f = 0; f < NDOF_P; f++ ) {
                      FeM ( i )  -=  J_Wg_Alpha_phiig * ( Pressure[f] - Pressure_old[f] ) * _dphi_g[1][f + ivar * el_ndof[1]];
//                         FeM(i)  +=  JxW_g2 * alpha_beta * rho * (Pressure[f] /*- p_old[f]*/) * _phi_g[1][f] * _dphi_g[2][i + ivar * el_ndof[2]]; // integration by parts (Pressure[f] /*- p_old[f]*/)
                      }

                  FeM ( i )  += ( 3. / 2. ) * TimeDisc * J_Wg_Alpha_phiig * _ub_g[2][_FF_idx[NS_F] + ivar] / _dt;

                  if ( ivarN == _dir ) {
                      for ( int j = 0; j < el_ndof2; j++ ) {
                          KeM ( i, j ) += ( 3. / 2. ) * TimeDisc * J_Wg_Alpha_phiig * _phi_g[2][j] / _dt;
                          } // end A element matrix quad -quad (end loop on j)---------------------------------------------
                      }
                  else {
                      FeM ( i )  -= ( 3. / 2. ) * TimeDisc * J_Wg_Alpha_phiig * ( u_nlg[ivar] ) / _dt;
                      }
                  }
              }
          }// end loop i
      }// end of the quadrature point qp-loop

  for ( int dd = 0; dd < 3; dd++ ) {
      delete [] ProjDir[dd];
      }

  delete [] ProjDir;

// ====================== end volume (element) =======================================
  return;
  }







#endif
#endif









