// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include "MGSolverNS_proj.h"       // Navier-Stokes class header file
#include "MGFE_conf.h"        // FEM approximation
#include "MGFE.h"          // Mesh class

// CONSTRUCTOR ===========================================================================
MGSolNS_proj::MGSolNS_proj (
  MGEquationsSystem & mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
) :
  MGSolNS ( mg_equations_map_in, nvars_in, eqname_in, varname_in ) {
  _Coupled = 0;
  SetVariableNames ( varname_in );
  }


// SET VARIABLE NAMES
void MGSolNS_proj::SetVariableNames ( std::string varname_in ) {

  if ( !varname_in.compare ( "u" ) ) {
      _dir = 0;   // u-equation
      }

  if ( !varname_in.compare ( "v" ) ) {
      _dir = 1;   // v-equation
      }

  if ( !varname_in.compare ( "w" ) ) {
      _dir = 2;   // w-equation
      }

  _var_names[0] = varname_in;
  _refvalue[0] = _uref;

  return;
  }

// GET OTHER SYSTEM DATA =================================================================
void MGSolNS_proj::get_el_field_data (
  int iel, int Level,
  int el_conn [], int offset, int el_ndof[], int ndof_lev ) {
  int el_ndofp = el_ndof[1];

  for ( int deg = 0; deg < 3; deg++ )
    for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                               el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
        }

  for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_nonl_sol ( 0, 1, el_ndof[2], el_conn, offset, kdim, _u_nl );
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, kdim, _u_1ts );
      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_oldsol ( 0, 1, el_ndof[2], el_conn, offset, kdim, _u_2ts );
      }

  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol ( 0, 1, el_ndofp, el_conn, offset, 0, _p_1ts );       // time step -1
  _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol ( 0, 1, el_ndofp, el_conn, offset, 0, _p_2ts );   // time step -2

  for ( int i = 0; i < NDOF_P; i++ ) {
      if ( _NS_parameter._TimeDisc == 2 )
        _Pressure[i] = ( 2.*_p_1ts[i] - _p_2ts[i] );

      if ( _NS_parameter._TimeDisc == 1 )
        _Pressure[i] = ( _p_1ts[i] );
      }

  return;
  }

// MOMENTUM EQUATION DISCRETIZATION ======================================================
void MGSolNS_proj::MomentumEquation ( double JxW_g2, int el_ndof[], int qp )
  {
  
  double  dphiidx_g2[DIMENSION], dphijdx_g2[DIMENSION];
  for ( int i = 0; i < el_ndof[2]; i++ ) { // LOOP OVER TEST FUNCTIONS -------------------------------------------------
      // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
      int indx_eq     = i + _dir * el_ndof[2];
      double phii_g = _phi_g[2][i];

      for ( int idim = 0; idim < _nNSdim; idim++ )
        dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof[2]];


      RowSetUp ( i, indx_eq, qp, el_ndof )  ;

      // zero line for Dirichlet  bc -----------------------------------------------------------
      if ( _bc_el[indx_eq] < 0 || ( _bc_el[indx_eq] == 0 && ( _bc_vol[i] == 11 || _bc_vol[i] == 31 ) ) ) {
          for ( int  ivarN = _dir; ivarN < _dir + _NComp; ivarN++ )     { // Loop over velocity components for equation written in row indx
              const int ivar = ivarN % _nNSdim;
              // computation of alpha beta (factor for adding equations) ----------------------------------------------
              double alpha_beta = 1.;

              if ( _BoundEquation ) alpha_beta = _ProjDir[abs ( _bc_el[indx_eq] ) - 1][ivar];

              double  dtxJxW_g = JxW_g2 * alpha_beta; // (factor for adding equations)

              // Regularization supg term ------------------------------------------------------------------------------------
              double TimeDer=0.;

              if ( _NS_parameter._SolveSteady == 0 ) {
                  if ( _NS_parameter._TimeDisc == 1 )
                    TimeDer =  ( _u_1ts_g[ivar] ) / ( _dt );

                  if ( _NS_parameter._TimeDisc == 2 )
                    TimeDer =  ( 2.*_u_1ts_g[ivar] - 0.5 * _u_2ts_g[ivar] ) / ( _dt );
                  }

              _FeM ( i )  +=  dtxJxW_g *  (
                               TimeDer
                               + _IFr * _dirg[ivar]   // force
                               - _P_gdx[ivar]        // pressure
                             ) * ( phii_g + _Phi_supg );

              for ( int  ot = ivar; ot < ivar + _nNSdim; ot++ ) { // out of diagonal tensor components -> better explicit
                  const int ivars = ot % _nNSdim;
                  _FeM ( i )   -= _IReEff * dtxJxW_g *  _vel_gdx[ivars * _nNSdim + ivar] * dphiidx_g2[ivars]; //
                  }

              for ( int j = 0; j < el_ndof[2]; j++ ) {
                  const double phij_g = _phi_g[2][j];
                  double Lap_g = 0., Adv_g = 0., Supg_lap = 0., up = 0.;

                  if ( _AxiSym == 1 ) { // axysimmetry only --------------------------------------------------
                      Lap_g = 2.* ( 1 - ( _dir ) ) *
                              ( _IReEff + _NS_parameter._Upwind * _f_upwind * _u_1ts_g[0] * _u_1ts_g[0] )
                              * phij_g * ( phii_g + _Phi_supg ) / ( _xyzg[0] * _xyzg[0] );
                      } // --------------------------------------------------------------------------------

                  for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
                      dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof[2]];

                      if ( _WallElement != 1 ) {
                          up = _NS_parameter._Upwind * _f_upwind * _u_1ts_g[kdim] * _u_1ts_g[kdim];
                          }

                      if ( _NS_parameter._TimeDisc == 1 )
                        Adv_g     += ( _u_1ts_g[kdim] ) * dphijdx_g2[kdim];

                      if ( _NS_parameter._TimeDisc == 2 )
                        Adv_g     += ( 2.*_u_1ts_g[kdim] - _u_2ts_g[kdim] ) * dphijdx_g2[kdim];

                      Lap_g     += ( _IReEff + up ) * dphijdx_g2[kdim] * dphiidx_g2[kdim];
                      Supg_lap  += ( _IReEff ) * _ddphi_g[2][j * _nNSdim * _nNSdim + kdim * _nNSdim + kdim];
                      }

                  //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
                  double TimeDerivative;

                  if ( _NS_parameter._TimeDisc == 2 )
                    TimeDerivative = dtxJxW_g * ( 3. / 2. ) *  ( phii_g + _Phi_supg ) * phij_g / ( _dt );

                  if ( _NS_parameter._TimeDisc == 1 )
                    TimeDerivative = dtxJxW_g *   ( phii_g + _Phi_supg ) * phij_g / ( _dt );

                  const double KeM_ij = dtxJxW_g  * (
                                          Adv_g * ( phii_g + _Phi_supg )  // advection
                                          + Lap_g                        // viscous Laplacian
                                          - Supg_lap * _Phi_supg          // SUPG viscous Laplacian
                                        );
                  
                  // non Diagonal blocks
                  if ( ivarN == _dir ) {
                      if ( _NS_parameter._SolveSteady == 0 ) _KeM ( i, j ) += TimeDerivative;
                      _KeM ( i, j ) += KeM_ij;
                      }
                  else {
                      if ( _NS_parameter._SolveSteady == 0 ) _FeM ( i )  -=  TimeDerivative * _u_1ts[j + ivar * NDOF_FEM];
                      _FeM ( i )   -= KeM_ij * _u_1ts[j + ivar * NDOF_FEM];
                      }
                  } // end A element matrix quad -quad (end loop on j)---------------------------------------------
              }
          } // end loop ivar
          
        
      } // end loop i

  return;
  }

// BOUNDARY CONDITIONS =================================================================== 
void MGSolNS_proj::set_bc_matrix(
    int dir_maxnormal,                          ///<  normal dir
    int sur_toply[],                            ///< boundary topology map
    int el_ndof[],                              ///< number of volume dofs
    int elb_ndof[],                             ///< number of boundary dofs
    int elb_ngauss,                             ///<  number of surface gaussian points
    double normal[],                            ///< normal
    int el_conn[]
)
{
    int bd_node, bc_n_flag, bc_tg_flag, bc_var_normal;
    for (int dir = 0; dir < _nNSdim; dir++) {
        _ProjDir[2][dir] = normal[dir];
    }

    const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;
    const int norm_face = (abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 10000) / 1000 - 1 ;

    // TANGENT VECTORS BASED ON BOUNDARY GEOMETRY
    CalcTangDir(_ProjDir[0], _ProjDir[1], _ProjDir[2], GEOM_BASED);  // for dirichlet bc

    
    SetBCFlags(normal, sur_toply, el_ndof[2], el_conn);
    double u_oold[DIMENSION * NDOF_FEM];
    for (int kdim = 0; kdim < _nNSdim; kdim++) {
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_oldsol(0, 1, el_ndof[2], el_conn, _offset, kdim, u_oold);
    }

    if (bd_face == 88 || bd_face == 84) {
        if (_FF_idx[K_F] >= 0) {
            if (_NS_parameter._WallFunctionApproach == 1) {
                _WallElement = 1;
            }
            for (int dir = 0; dir < _nNSdim; dir++) {
                _normal_pt[(NDOF_FEM - 1) *_nNSdim + dir] += normal[dir];
            }
        }
    }
        
    // ...........................................................................................
    //     DIRICHLET BOUNDARY CONDITION
    // ...........................................................................................

    for (int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {

        
        NodeBCFlags ( sur_toply[lbnode], bd_node, bc_tg_flag, bc_n_flag, bc_var_normal );
        const int norm_node = (abs(_bc_bd[sur_toply[lbnode]]) % 10000) / 1000 - 1 ;
        
        if(el_conn[lbnode]==3)
            int a=1;
        
        if (_nNSdim == 3) {
            CorrectBCFlags(_ProjDir[0], _ProjDir[1], bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
        }

        if (bd_face == bd_node) {
            const int  indx_row = sur_toply[lbnode] + _dir * el_ndof[2]; // solution indx
            const int  node = sur_toply[lbnode];

            bool NormBC = (abs(_bc_el[indx_row]) == 3) ? true : false;
            bool DirichletBC = (_bc_el[indx_row] > 0) ? true : false;


            if (DirichletBC && norm_face == norm_node/*!_AlreadyWrittenDirBC[indx_row]*/) {
                int  bc_bc  = NormBC ? bc_n_flag : bc_tg_flag;
                const int  bc_rhs = ((bc_bc % 2) == 1) ? 1 : 0; // bc_rhs -> rhs bc
                bool OrDir = (bc_bc > 7) ? true : false;
                const double vart = bc_rhs * _u_1ts[node + _dir * el_ndof[2]];

                double prod = 0., u_norm[DIMENSION];
                for (int i = 0; i < _nNSdim; i++) {
                    prod += _u_1ts[node + i * el_ndof[2]] * normal[i];
                }
                for (int i = 0; i < _nNSdim; i++) {
                    u_norm[i] = prod * normal[i];
                }

                double alpha_beta = 1.;
                if (bd_node == 88 || bd_node == 66) {
                    _KeM(node, node) = 1.;
                    _FeM(node)      = 0.;
                } else if (OrDir) {
                    //  IMPOSITION OF VELOCITY VECTOR ALONG GENERAL DIRECTION f
                    //  \vec(U)_f = (\vec(U)_{old} \dot \hat(f)) \hat(f)
                    //  f can be normal or tangential direction
                    if (bd_node == 98) {
                        _KeM(node, node) = 1.;
                        _FeM(node)       = u_norm[_dir];
                    } else {             
                        int direction = -2*(abs(_bc_el[indx_row])/3) + 2;                     
                        double scal_prod = 0.; double vel_proj[DIMENSION];
                        for (int i = 0; i < _nNSdim; i++) {
                                scal_prod += _u_1ts[node + i * el_ndof[2]] * _ProjDir[direction][i];
                            }
                            for (int i = 0; i < _nNSdim; i++) {
                                vel_proj[i] = scal_prod * _ProjDir[direction][i];
                            }
                        _KeM(node, node)  = 1.;
                        _FeM(node)        = vel_proj[_dir];
                    }
                } else {
                    _KeM(node, node) = alpha_beta;
                    _FeM(node) = vart * alpha_beta;
                }// beta_0 u -> diagonal
            }// END IF DIRICHLET BC
        }// END IF BD_NODE == BD_FACE
    }// END LOOP OVER BOUNDARY NODES - LBNODE

    // ...........................................................................................
    //     INTEGRATION OVER BOUNDARY ELEMENT - STRESS CALCULATION
    // ...........................................................................................
    for (int  qp = 0; qp < elb_ngauss; qp++) { // START LOOP OVER GAUSS POINTS ===============================================
        double det   = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);   // local coord _phi_g and jac
        double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp]; // weight
        
        _fe[2]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[2]);   // global coord _phi_g
        _fe[1]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[1]);   // global coord _phi_g

        if (_AxiSym == 1) {
            interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg);
            JxW_g  *= _xyzg[0];   // axisymmetric  (index ->0)
        }
        const double dtJxW_g = JxW_g;        
        
        double VelOnGauss[DIMENSION];
        interp_el_bd_sol(_u_1ts, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], VelOnGauss);

        for (int i = 0; i < elb_ndof[2]; i++) { // LOOP OVER TEST FUNCTIONS ===================================================

            const int bc_var_check  = abs(_bc_bd[sur_toply[i]]) % 10000;                         // bc_var
            const int bc_var_normal = bc_var_check / 1000 - 1;                            //
            const int bd_node       = abs(_bc_bd[sur_toply[i]]) % 100;
            const int indx          = sur_toply[i] + _dir * el_ndof[2];
            const int NComp         = (_bc_el[indx] < 0) ? _nNSdim : 1;
            const int bc_flag       = (_bc_el[indx] == -3) ? bd_node / 10 : bd_node % 10;  // normal or tangential bc flags
            const int pen_cond      = (bd_node == 44) ? 0 : 1;
            const double lambda1    = (1 - pen_cond) * _NS_parameter._Penalty_tg;
            const double beta       = lambda1 * _NS_parameter._Penalty_n;

            if (bd_face == bd_node) {
                // Direction of maximum vector component of tangent 1 and 2 -> alignment of equations
                if (_nNSdim == 3) {
                    CorrectBCFlags(_ProjDir[0], _ProjDir[1], bc_var_normal, sur_toply[i], el_ndof[2], 1);
                }
                if (bc_flag < 6 && bd_node != 11) { // only for neumann bc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // NORMAL VECTOR FOR THE PROJECTION OF VOLUME PART EQUATION
                    // _normal_pt is used to project the equation (volume part) using the same normal used in the boundary
                    if (qp == 0) for (int dir = 0; dir < _nNSdim; dir++) {
                            _normal_pt[ sur_toply[i]*_nNSdim + dir] += normal[dir];
                        }
                    for (int  ivarN = _dir; ivarN < _dir + NComp; ivarN++) { // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
                        const int ivar = ivarN % DIMENSION;
                        // Stress component along tangential directions (-1,-2) ---------------------------------------------------
                        if (_bc_el[indx] == -1 || _bc_el[indx] == -2) {
                            if (bc_flag != 1) {
//                                    _                              _
//                                   /   =>   ->    -->             /    ->   -->
//                                  /   (T *  n) *  phi_t ds  =    /  a (u *  phi_t) ds
//                                -s                             -s
                                double alpha_beta = _ProjDir[abs(_bc_el[indx]) - 1][ivar];
                                if (bc_flag == 5) {
                                    alpha_beta = _ProjDir[0][ivar];
                                }
                                double tau = (_bc_el[indx] == -1) ? (_NS_parameter._Tg1_stress) : _NS_parameter._Tg2_stress ;    //utau*utau/vel_bound;
                                double stress = alpha_beta * (pen_cond * tau - lambda1) * dtJxW_g * _phi_g[2][i];

                                if (ivarN == _dir)
                                    for (int j = 0; j < elb_ndof[2]; j++) {
                                        _KeM(sur_toply[i], sur_toply[j]) -= stress * _phi_g[2][j];
                                    }
                                else {
                                    _FeM(sur_toply[i]) += stress * VelOnGauss[ivar];
                                }
                            }
                        }// End loop if equation is for tangential direction-----------------------------------------------------
                        else if (_bc_el[indx] == -3) {
                            // normal with penalty (u dot n) n dot phi
//                                    _
//                                   /   ->   ->  ->   -->
//                                  /   (u *  n)  n *  phi_n ds = 0
//                                -s
                            if (bc_flag == 4) {
                                double alpha_beta = 1.;
                                alpha_beta = normal[_dir];
                                double norm_stress = dtJxW_g * alpha_beta * beta * _phi_g[2][i] * normal[_dir];
                                if (ivarN == _dir)
                                    for (int j = 0; j < elb_ndof[2]; j++) {
                                        _KeM(sur_toply[i], sur_toply[j]) -= norm_stress * normal[ivar] * _phi_g[2][j];
                                    }
                                else {
                                    _FeM(sur_toply[i]) += norm_stress * normal[ivar] * VelOnGauss[ivar];
                                }
                            } //  ---------------------------------------------------------------------------------------------------
                        }// contribution along normal direction
                    }// END LOOP IF BC FLAG IS FOR INTEGRATION
                }// END LOOP OVER VARIABLES FOR LINEAR COMBINATION OF EQUATIONS ++++++++++++++++++++++++++++++++++++++++++++++
            }// END IF BC_NODE == BC_FACE
        }// END LOOP OVER BOUNDARY TEST FUNCTIONS =========================================================================
    }// END GAUSSIAN INTEGRATION ========================================================================================
    return;
}


#endif









