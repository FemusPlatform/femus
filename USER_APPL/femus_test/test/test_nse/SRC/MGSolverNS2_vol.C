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

// Thermodinamical Properties ==========================================================================================
// constant
// double  _NSdensity(double T,double T_r) {return 1.;} // water (t K) T_ref
// double  _NSkviscosity(double T,double T_r) {return 1.;} // lead T K
// water
//     double  _NSdensity(double T,double T_r){return (1. - (T-3.9863)*(T-3.9863)*(T+288.9414)/(508929.2*(T+68.12963)))/
//         (1. - (T_r-3.9863)*(T_r-3.9863)*(T_r+288.9414)/(508929.2*(T_r+68.12963)));} // water (t C) T_ref
//      double  _NSkviscosity(double T){return 1.;} //  T K
//  // Lead
//     double  _NSdensity(double T,double T_r){return (11441-1.2795*T)/(11441-1.2795*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return return 4.55e-4*exp(1069/T);} // lead T K
//       // LBE
//     double  _NSdensity(double T,double T_r){return (11065-1.293*T)/(11065-1.293*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return 4.94e-4*exp(754.1/T);} // lead T K
//
//  ======================================================================================================================

void MGSolNS::get_el_field_data(
    int iel, int Level,
    int el_conn [], int offset, int el_ndof[], int ndof_lev,
    double u_old[],  double u_oold[],   double u_nl[],
    double p_proj[], double dp_proj[]
)
{
    // external fields (from constant 0 to quadratic 2) -------------------------------------------------------------------
//      for ( int deg=0; deg<3; deg++ )
//           for ( int eq=0; eq<_data_eq[deg].n_eqs; eq++ ) {
//                _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub );
//           }
//      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol ( _nNSdim,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub ); // pressure
//      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol ( 0,_nNSdim,el_ndof[2],el_conn, offset,0,u_old );          // old vel
//      _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_nonl_sol ( 0,_nNSdim,el_ndof[2],el_conn, offset,0,u_nl );      //  non linear vel

    int el_ndofp = el_ndof[1];

    // element nodes coordinates ----------------------------------------------------------------------------------------
//   for(int idim=0; idim<DIMENSION; idim++) for(int d=0; d< NDOF_FEM; d++) {
//       _data_eq[2].ub[idim*NDOF_FEM+d]=_xx_qnds[idim*NDOF_FEM+d];
//     }
    // external fields (from constant 0 to quadratic 2) -------------------------------------------------------------------
    for (int deg = 0; deg < 3; deg++)
        for (int eq = 0; eq < _data_eq[deg].n_eqs; eq++) {
            _data_eq[deg].mg_eqs[eq]->get_el_sol(0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                                 el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub);
        }

    for (int kdim = 0; kdim < _nNSdim; kdim++) {
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_nonl_sol(0, 1, el_ndof[2], el_conn, offset, kdim, u_nl);
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_sol(0, 1, el_ndof[2], el_conn, offset, kdim, u_old);
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_oldsol(0, 1, el_ndof[2], el_conn, offset, kdim, u_oold);
    }
    _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0, 1, el_ndofp, el_conn, offset, 0, p_proj);  // old pressure
    _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0, 1, el_ndofp, el_conn, offset, 0, dp_proj); // dp pressure


    return;
}

// ==============================================================================================
void MGSolNS::matrixrhsvol(
    DenseMatrixM &KeM, DenseVectorM &FeM,
    int el_ndof[],
    double u_old[],
    double u_nl[],
    const int unsteady_flag,
    const int mode,
    double NodeMuTurb[],
    double WallDist[],
    int el_conn[],
    double p_proj[],
    double dp_proj[]
)
{

    double dphijdx_g2[DIMENSION];
    double dphiidx_g2[DIMENSION];
    double vel_gddx[DIMENSION * DIMENSION * DIMENSION];
    double vel_gdx[DIMENSION * DIMENSION];
    double vel_g[DIMENSION], u_nlg[DIMENSION];
    double MuTurb[1];
    double x_m[DIMENSION];
    double MuTurbDxg[DIMENSION];
    double tang1[DIMENSION], tang2[DIMENSION];
    double rho = 1.;
    int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim - 1];              //elem gauss points
    const int el_ndof2 = el_ndof[2];
    double dt = _dt;
    if (dt > 1 || _NS_parameter._SolveSteady == 1) {
        dt = 1. ;
    }
    if (dt < 1.e-10) {
        dt = 1.e-10;
    }
    double Pressure[NDOF_P], P_gdx[DIMENSION];
    for (int i = 0; i < NDOF_P; i++) {
        Pressure[i] = p_proj[i] /*+ dp_proj[i]/_dt*/;
    }

    double utau = 0.;
    double WallOnG[1];
    if (_WallElement > 0) {
        double max_wd = -1.;
        int indx_max = 1;
        for (int k = 0; k < el_ndof2; k++) {
            if (max_wd < WallDist[k] + 1.e-10) {
                max_wd = WallDist[k];
                indx_max = k;
            }
        }
        const int indx_wf = indx_max; // */NDOF_FEM -1  ;
        double Vel[DIMENSION], VelNorm[DIMENSION], UMag = 1.e-20;
        double UNormMag = 1.e-20;
        for (int dir = 0; dir < DIMENSION; dir++) {
            Vel[dir] = _data_eq[2].ub[(_FF_idx[NS_F] + dir) * NDOF_FEM + indx_wf];
        }
        for (int dir = 0; dir < DIMENSION; dir++) {
            VelNorm[dir] = Vel[dir] * _normal_pt[(NDOF_FEM - 1) * _nNSdim + dir];
            UNormMag    += VelNorm[dir] * VelNorm[dir];
            UMag        += Vel[dir] * Vel[dir];
        }
        double UTangMag = sqrt(UMag - UNormMag);
        utau = _mgutils._TurbParameters->CalcUtau(UTangMag, WallDist[indx_wf]);
        _mgutils._TurbParameters->_utau = utau;
    }

    
//     double vel_dp[DIMENSION*NDOF_FEM];
//     for(int node = 0; node<NDOF_FEM; node++){
//       double det = _fe[1]->Jac3D_nodes(node, _xx_qnds, _InvJac1);
//       std::cout<<det<<std::endl;
//       _fe[1]->get_dphi_node(_nNSdim,node,_InvJac1,_dphi_g[1]);
//       interp_el_gdx(Pressure, 0, 1, _dphi_g[1], el_ndof[1], P_gdx); // derivatives  vel_gdx[DIM][DIM]
//       for (int idim = 0; idim <  _nNSdim; idim++) vel_dp[idim*NDOF_FEM + node] = -1.*_dt*P_gdx[idim];
//     }
    
    for (int qp = 0; qp <  el_ngauss; qp++) {
        // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
        // quadratic continuous (2)  (velocity)
        const double det2      = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);   // quadratic Jacobian
        const double det1      = _fe[1]->Jac(qp, _xx_qnds, _InvJac1);   // quadratic Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[ _nNSdim - 1][qp]/**dt*/;      // quadratic weight
        _fe[2]->get_phi_gl_g(_nNSdim, qp, _phi_g[2]);                   // quadratic shape function
        _fe[2]->get_dphi_gl_g(_nNSdim, qp, _InvJac2, _dphi_g[2]);       // global coord deriv
        _fe[2]->get_ddphi_gl_g(_nNSdim, qp, _InvJac2, _ddphi_g[2]);     // global second derivatives
        _fe[1]->get_phi_gl_g(_nNSdim, qp, _phi_g[1]);                   // linear shape funct
        _fe[1]->get_dphi_gl_g(_nNSdim, qp, _InvJac1, _dphi_g[1]);       // global coord deriv
        
        
        // discontinuous (0) (disc pressure)
        if (_nvars[0] > 0) {
            _fe[0]->get_phi_gl_g(_nNSdim, qp, _phi_g[0]);  // piecewise shape function
        }
        interp_el_sol(_xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof2, _xyzg);
        // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------
        // quadratic fields (velocity)
        interp_el_sol(_data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2]); // field _ub_g[2][DIM]
        interp_el_gdx(u_old, 0, _nNSdim, _dphi_g[2], el_ndof2, vel_gdx); // derivatives  vel_gdx[DIM][DIM]
        interp_el_gdx(Pressure, 0, 1, _dphi_g[1], el_ndof[1], P_gdx); // derivatives  vel_gdx[DIM][DIM]
        interp_el_gddx(u_old, 0, _nNSdim, _ddphi_g[2], el_ndof2, vel_gddx);

//         double vel_dp_g[DIMENSION];
//         interp_el_sol(vel_dp, 0, _nNSdim, _phi_g[2], el_ndof2, vel_dp_g);
//         
//         std::cout<<" ----------------------\n";
//         for(int f=0; f<_nNSdim; f++) std::cout<<" dir "<<f<<" from nodal der "<<vel_dp_g[f]<<" from P_gdx "<< -1.*_dt*P_gdx[f]<<std::endl;
        
        // linear fields (pressure) -> projection and segregated

        if (_AxiSym == 1) {
            JxW_g2  *= _xyzg[0];
        }

        // Velocity, Reynolds and upwind  --------------------------------------------------------------------------
        double mod2_vel = 1.e-20;
        _sP = 1.e-20;
        double mu = 1.; // Velocity
        for (int idim = 0; idim <  _nNSdim; idim++) {
            vel_g[idim] = 0.;
            x_m[idim] = 0.;
            u_nlg[idim] = 0.; // old and  non linear Velocity at gaussina point qp
            for (int k = 0; k < NDOF_FEM; k++) {
                const int dnode = k + idim * NDOF_FEM;
                x_m[idim] += _xx_qnds[dnode] / el_ndof[2];
                u_nlg[idim] += u_nl[dnode] * _phi_g[2][k];
                vel_g[idim] += u_old[dnode] * _phi_g[2][k];
            }
            mod2_vel += u_nlg[idim] * u_nlg[idim]; // module velocity
            for (int jdim = 0; jdim <  _nNSdim; jdim++) {
                _sP += 0.25 * (vel_gdx[ jdim * _nNSdim + idim] + vel_gdx[ idim * _nNSdim + jdim]) *
                       (vel_gdx[ jdim * _nNSdim + idim] + vel_gdx[ idim * _nNSdim + jdim]);
            }
        }
        
//         for (int idim = 0; idim <  _nNSdim; idim++) vel_g[idim] -= _dt* P_gdx[idim]   ; 
        
        
        mod2_vel = sqrt(mod2_vel);

        _mu_turb = 0.; // eff visc turb at g pt
        MuTurb[0] = 0.;
        for (int dir = 0; dir < _nNSdim; dir++) {
            MuTurbDxg[dir] = 0.;
        }
        // -------------------- Turbulence [K_F] -> (quad,_indx_eqs[K_F]) -------------------------------------
        if (_FF_idx[K_F] >= 0) { // turbulence _mu_turb evaluation
            _kappa_g[0] = _ub_g[2][_FF_idx[K_F]];
            _kappa_g[1] = _ub_g[2][_FF_idx[K_F] + 1];
            const int klog = _mgutils._TurbParameters->_klog;
            const double kappa = klog * exp(klog * _kappa_g[0]) + (1 - klog) * _kappa_g[0];
            _mu_turb = _mgutils._TurbParameters->CalcMuTurb(_kappa_g, _y_bcout, _sP);    // piece wise mu turb
            const double MuDurbin = kappa / (sqrt(_sP + 1.e-10) * _IRe);    // Durbin limit value
            if (_NS_parameter._InterpolatedMuTurb == 1 || _WallElement == 1) {
                if (_NS_parameter._InterpolatedMuTurb == 1) {
                    interp_el_gdx(NodeMuTurb, 0, 1, _dphi_g[2], el_ndof2, MuTurbDxg);
                }
                interp_el_sol(NodeMuTurb, 0, 1, _phi_g[2], el_ndof2, MuTurb); // muturb interpolated from the node values
                _mu_turb = MuTurb[0];
                if (_WallElement == 1) {
                    interp_el_sol(WallDist, 0, 1, _phi_g[2], el_ndof2, WallOnG); // muturb interpolated from the node values
                    double yplus2 = WallOnG[0] * utau / _IRe ;
//                          _mu_turb = 1./ ( 1./ ( 0.001093*yplus2*yplus2*yplus2 ) + 1./ ( 0.41*yplus2 ) );
                }
            }
            _mu_turb = min(MuDurbin, _mu_turb);
        }
        double IRe_eff = _IRe * (rho + _mu_turb);
        double f_upwind = CalcFUpwind(vel_g, _dphi_g[2], IRe_eff, _nNSdim, el_ndof2);

        double Div_g = 0.;
        for (int i = 0; i < _nNSdim; i++) {
            Div_g += vel_gdx[ i * _nNSdim + i];
        }

        // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
//     if (_FF_idx[T_F]>=0) {
//       double Temp_g=_ub_g[2][_FF_idx[T_F]];
//       rho = (11065-1.293* (Temp_g+273.15)) / (_rhof);
//     }

        // ===========================================================================================================
        //                                       D) Assembling NS equation
        // ===========================================================================================================



        for (int i = 0; i < el_ndof2; i++) { // LOOP OVER TEST FUNCTIONS --------------------------------------------------
            // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg

            double phii_g = _phi_g[2][i];
            for (int idim = 0; idim < _nNSdim; idim++) {
                dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof2];
            }
//                for ( int  row_shift=0; row_shift< _nvars[2]; row_shift++ )    { // LOOP OVER ROWS: i + row_shift*el_ndof2
            int indx_eq     = i + _dir * el_ndof2;

            // in the curved cylinder this is needed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (_bc_el[indx_eq] == -8) {
                _bc_el[indx_eq] = 0;
                std::cout << "ATTENTION!!! el node " << i << " glob node " << el_conn[i] << " row " << indx_eq << " bc " << _bc_el[indx_eq] << "   " << _dir << std::endl;
            }

            bool BoundEquation = (_bc_el[indx_eq] <= -1) ? true : false;
            bool NormalEquation = false;
            if (BoundEquation && _bc_el[indx_eq] == -3) {
                NormalEquation = true;
            }
            if (BoundEquation && !NormalEquation) {
                CalcTangDir(i, el_ndof[2], _bc_el[indx_eq], qp, tang1, tang2, 1);
            }
            int NComp = 1;                                 // Default: 1 vel comp for equation in row indx -> ivar = row_shift
            if (BoundEquation) {
                NComp = _nNSdim;     // For boundary sides (tangential or normal required) linear combination of _nnsdim equations
            }


            // zero line for Dirichlet  bc -----------------------------------------------------------
            if (_bc_el[indx_eq] <= 0) {
                for (int  ivarN = _dir; ivarN < _dir + NComp; ivarN++)     { // Loop over velocity components for equation written in row indx
                    const int ivar = ivarN % _nNSdim;

                    const int ImpEq = (ivarN > _dir) ? 0 : 1;

                    // computation of alpha beta (factor for adding equations) ----------------------------------------------
                    double alpha_beta = 1.;
                    if (BoundEquation) {
                        if (NormalEquation) {
                            alpha_beta = _normal_pt[ivar + i * _nNSdim];
                        } else {
                            alpha_beta = (_bc_el[indx_eq] == -1) ? tang1[ivar] : tang2[ivar];
                        }
                    }
                    double  dtxJxW_g = JxW_g2 * alpha_beta; // (factor for adding equations)

                    // Regularization supg term ------------------------------------------------------------------------------------
                    double Phi_supg = 0., Phi_supg2 = 0.;
                    if (_bc_el[indx_eq] == 0 && _NS_parameter._Supg == 1) {
                        for (int idim = 0; idim < _nNSdim; idim++) { // Phi_supg = 0. on boundaries with Dirichlet bc
                            Phi_supg += f_upwind * vel_g[idim] * dphiidx_g2[idim]; //f_upwind
                        }
                    }
                    // -------------------------------------------------------------------------------------------------------------
                    // -------------------------------------- Assemblying rhs ------------------------------------------------------

                    // -------------------------------------------------------------------------------------------------------------
                    //-------------------------------------- Assemblying matrix ----------------------------------------------------

                    if (mode == 1) {   // rsh  terms
                        if (_NS_parameter._SolveSteady == 0) {
                            FeM(i)  +=  dtxJxW_g * rho * vel_g[ivar] * (phii_g + Phi_supg) / _dt;
                        }
                       FeM(i)  +=  dtxJxW_g * rho * (
                                        _IFr * _dirg[ivar]   // force
//                                         - P_gdx[ivar]        // pressure 
                                    ) * (phii_g + Phi_supg);
                        for (int  ot = ivar; ot < ivar + _nNSdim; ot++) {// out of diagonal tensor components -> better explicit 
                            const int ivars = ot % _nNSdim;
                            FeM(i)   -= IRe_eff * dtxJxW_g * rho * vel_gdx[ivars * _nNSdim + ivar] * dphiidx_g2[ivars]; //
                        }           
                    }
                    if (ImpEq == 1) {
                        for (int j = 0; j < el_ndof2; j++) {
                            const double phij_g = _phi_g[2][j];
                            double Lap_g = 0., Adv_g = 0., Supg_lap = 0., Turb_grad = 0., up = 0.;
                            if (_AxiSym == 1) { // axysimmetry only --------------------------------------------------
                                Lap_g = 2.* (1 - (_dir)) *
                                        (IRe_eff + _NS_parameter._Upwind * f_upwind * u_nlg[0] * u_nlg[0])
                                        * phij_g * (phii_g + Phi_supg) / (_xyzg[0] * _xyzg[0]);
                            } // --------------------------------------------------------------------------------

                            for (int kdim = 0; kdim < _nNSdim; kdim++) {
                                dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof2];
                                if (_WallElement != 1) {
                                    up = _NS_parameter._Upwind * f_upwind * vel_g[kdim] * vel_g[kdim];
                                }
                                Adv_g     += (u_nlg[kdim] /*- _dt*P_gdx[kdim]*/)* dphijdx_g2[kdim];                                                                // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
                                Lap_g     += (IRe_eff + up) * dphijdx_g2[kdim] * dphiidx_g2[kdim];
                                Supg_lap  += (IRe_eff) * _ddphi_g[2][j * _nNSdim * _nNSdim + kdim * _nNSdim + kdim];
                            }

                            if (_FF_idx[K_F] >= 0 && _NS_parameter._InterpolatedMuTurb == 1) // grad(nut)*[grad(u) + grad(u)^T]*phi_supg, diagonal part
                                for (int kdim = 0; kdim < _nNSdim; kdim++) {
                                    Turb_grad += _IRe * MuTurbDxg[kdim] * _dphi_g[2][j + kdim * el_ndof2];
                                }

                            //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
                            if (_NS_parameter._SolveSteady == 0) {
                                KeM(i, j) += dtxJxW_g * rho * (phii_g + Phi_supg) * phij_g / _dt;
                            }

                            KeM(i, j) += dtxJxW_g * rho * (
                                             Adv_g * (phii_g + Phi_supg) // time + advection
                                             + Lap_g                       // viscous Laplacian
//                                              - Supg_lap * Phi_supg          // SUPG viscous Laplacian
//                                              - Turb_grad * Phi_supg
                                         );
//                             KeM(i, j) += IRe_eff * dtxJxW_g * rho * _dphi_g[2][j + ivar * el_ndof2] * dphiidx_g2[ivar]; //
                        } // end A element matrix quad -quad (end loop on j)---------------------------------------------
                    }
                    if (ImpEq == 0) {
                        double ExLap_g = 0., ExAdv_g = 0., ExSupg_lap = 0., ExTurb_grad = 0., Exup = 0.;
                        for (int kdim = 0; kdim < _nNSdim; kdim++) {
                            if (_WallElement != 1) {
                                Exup = _NS_parameter._Upwind * f_upwind * vel_g[kdim] * vel_g[kdim];
                            }
                            ExAdv_g     += (u_nlg[kdim] /*- _dt*P_gdx[kdim]*/)* vel_gdx[ ivar * _nNSdim + kdim];                                                            // Adv_g +=vel_g[kdim]*dphijdx_g2[kdim]*phii_g;
                            ExLap_g     += (IRe_eff + Exup) * vel_gdx[ ivar * _nNSdim + kdim] * dphiidx_g2[kdim];
//                                    ExSupg_lap  += ( IRe_eff ) *_ddphi_g[2][j* _nNSdim* _nNSdim+kdim* _nNSdim+kdim];
//                                    ExSupg_lap  += ( IRe_eff ) * vel_gddx[j*_nNSdim*_nNSdim   +kdim*_nNSdim  +kdim];
                        }
                        if (_NS_parameter._SolveSteady == 0) {
                            FeM(i)  -=  dtxJxW_g * rho * vel_g[ivar] * (phii_g + Phi_supg) / _dt;
                        }
                        FeM(i)   -= dtxJxW_g * rho * (
                                        ExAdv_g * (phii_g + Phi_supg) // advection
                                        + ExLap_g
                                        + IRe_eff * vel_gdx[ ivar * _nNSdim + ivar] * dphiidx_g2[ivar] //
                                    );
                    }
                }
            } // end loop ivar
//                }// END LOOP OVER MATRIX ROWS: i + row_shift*el_ndof2
        } // end loop i

        //------------------    QL    -----------------------------------------------------------

//           for ( int  ikl=1; ikl<2; ikl++ ) { // ikl=0 discontinuous ikl=1 continuous pressure
//                for ( int   i=0; i< el_ndof[ikl]; i++ ) {
//                     const int    indp      = i+el_ndof2*_nNSdim;
//                     const double psii_g    = _phi_g[ikl][i];
//                     const double dtxJxWp_g = ( _bc_el[indp] != 0 ) ? JxW_g2 : 0.;
//                     for ( int j=0; j<el_ndof2; j++ )     {
//                          for ( int  jvar=0; jvar< _nvars[2]; jvar++ ) { // linear -quad
//                               if ( _AxiSym ==1 && _nNSdim == 2 ) {
//                                    const double phij_g= _phi_g[2][j];
//                                    KeM ( indp,j+jvar*el_ndof2 ) += dtxJxWp_g* ( 1- jvar ) *psii_g*phij_g/_xyzg[0];
//                               }
//                               KeM ( indp,j+jvar*el_ndof2 ) +=  dtxJxWp_g*psii_g*_dphi_g[2][j+jvar*el_ndof2];
//                          }// jvar
//                     }// j end linear-quad --------------------------------------------
//                     // linear-linear Assemblying Matrix ------------------------------
//                     for ( int j=0; j<el_ndof[1]; j++ )  {
//                          const double psij_g=_phi_g[1][j];
//                          KeM ( indp,j+ _nNSdim*el_ndof2 )  += dtxJxWp_g*rho*psii_g*psij_g*1.e-40; //0.0* KOMP_NS;//
//                     } // end linear liner -----------------------------------------------------
//                }// i
//           }// ikl=0 discontinuous ikl=1 continuous pressure
    }// end of the quadrature point qp-loop

// ====================== end volume (element) =======================================
    return;
}



// ==============================================================================================
void MGSolNS::matrixrhsvol_corr_step(
    DenseMatrixM &KeM, DenseVectorM &FeM,
    int el_ndof[],
    double u_old[],
    double u_nl[],
    const int unsteady_flag,
    const int mode,
    double NodeMuTurb[],
    double WallDist[],
    int el_conn[],
    double p_proj[],
    double dp_proj[]
)
{
    double rho = 1.;
    int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim - 1];              //elem gauss points
    const int el_ndof2 = el_ndof[2];

    double Pressure[NDOF_P], P_gdx[DIMENSION];
    for (int i = 0; i < NDOF_P; i++) {
        Pressure[i] = p_proj[i];
    }

    for (int qp = 0; qp <  el_ngauss; qp++) {
        // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
        // quadratic continuous (2)  (velocity)
        const double det2      = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);   // quadratic Jacobian
        const double det1      = _fe[1]->Jac(qp, _xx_qnds, _InvJac1);   // quadratic Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[ _nNSdim - 1][qp];      // quadratic weight
        _fe[2]->get_phi_gl_g(_nNSdim, qp, _phi_g[2]);                   // quadratic shape function
        _fe[1]->get_dphi_gl_g(_nNSdim, qp, _InvJac2, _dphi_g[1]);       // global coord deriv

        // discontinuous (0) (disc pressure)
        if (_nvars[0] > 0) {
            _fe[0]->get_phi_gl_g(_nNSdim, qp, _phi_g[0]);  // piecewise shape function
        }
        interp_el_sol(_xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof2, _xyzg);
        interp_el_sol(_data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2]); // field _ub_g[2][DIM]
        interp_el_gdx(Pressure, 0, 1, _dphi_g[1], el_ndof[1], P_gdx); // derivatives  vel_gdx[DIM][DIM]

        double tang1[DIMENSION], tang2[DIMENSION];
        for (int i = 0; i < el_ndof2; i++) { // LOOP OVER TEST FUNCTIONS --------------------------------------------------

            double phii_g = _phi_g[2][i];
            int indx_eq   = i + _dir * el_ndof2;
            int BoundEquation = (_bc_el[indx_eq] <= -1) ? 1 : 0;
            int NormalEquation = (BoundEquation == 1 && _bc_el[indx_eq] == -3) ? 1 : 0;
            if (BoundEquation == 1 && NormalEquation == 0) {
                CalcTangDir(i, el_ndof[2], _bc_el[indx_eq], qp, tang1, tang2, 1);
            }
            int NComp = (BoundEquation == 1) ? _nNSdim : 1;                              // Default: 1 vel comp for equation in row indx -> ivar = row_shift
            double alpha_beta = 1.;
            if (_bc_el[indx_eq] <= 0) {
                for (int  ivarN = _dir; ivarN < _dir + NComp; ivarN++)     { // Loop over velocity components for equation written in row indx
                    const int ivar = ivarN % _nNSdim;
                    const int ImpEq = (ivarN > _dir) ? 0 : 1;
                    if (BoundEquation == 1) {
                        if (NormalEquation == 1) {
                            alpha_beta = _normal_pt[ivar + i * _nNSdim];
                        } else {
                            alpha_beta = (_bc_el[indx_eq] == -1) ? tang1[ivar] : tang2[ivar];
                        }
                    }
                    double  J_Wg_Alpha_phiig = JxW_g2 * alpha_beta * rho * phii_g;
                    for (int f = 0; f < NDOF_P; f++) {
                       FeM(i)  -=  J_Wg_Alpha_phiig * Pressure[f] * _dphi_g[1][f+ivar*el_ndof[1]];
                    }
                    if (ImpEq == 1) {
                        FeM(i)  +=  J_Wg_Alpha_phiig * _ub_g[2][_FF_idx[NS_F] + ivar] / _dt;
                        for (int j = 0; j < el_ndof2; j++) {
                            KeM(i, j) += J_Wg_Alpha_phiig * _phi_g[2][j] / _dt;
                        } // end A element matrix quad -quad (end loop on j)---------------------------------------------
                    }
                }
            } // end loop i
        }// end of the quadrature point qp-loop
    }
// ====================== end volume (element) =======================================
    return;
}




void MGSolNS::CalcTangDir(
    double Tang1[],
    double Tang2[],
    double normal[],
    int Type
)
{
    CalcTangDir(0, NDOF_FEM, 10, 0, Tang1, Tang2, normal, Type);
    return;
}


void MGSolNS::CalcTangDir(
    int PhiNumber,
    int el_dof,
    int bc_el,
    int qp,
    double Tang1[],
    double Tang2[],
    int Type
)
{
    double normal[DIMENSION];
    for (int kk = 0; kk < DIMENSION; kk++) {
        normal[kk]  = _normal_pt[kk + PhiNumber * _nNSdim];
    }
    CalcTangDir(PhiNumber, el_dof, bc_el, qp, Tang1, Tang2, normal, Type);
    return;
}


void MGSolNS::CalcTangDir(
    int PhiNumber,
    int el_dof,
    int bc_el,
    int qp,
    double Tang1[],
    double Tang2[],
    double normal[],
    int Type
)
{
    // ---------------------------------------------------
    if (Type == 1) {

        double scal_unorm = 0., mod_norm = 0., mod_vel = 0.;
        double Normal[DIMENSION], Vel[DIMENSION];

        for (int kk = 0; kk < DIMENSION; kk++) {
            Normal[kk]  = normal[kk];
            mod_norm   += Normal[kk] * Normal[kk];
            Vel[kk]     = _data_eq[2].ub[(_FF_idx[NS_F] + kk) * el_dof + (NDOF_FEM - 1)];
            mod_vel    += Vel[kk] * Vel[kk];
        }

        mod_norm = sqrt(mod_norm + 1.e-20);
        mod_vel  = sqrt(mod_vel + 1.e-20);

        for (int kk = 0; kk < DIMENSION; kk++) {
            Normal[kk] = Normal[kk] / mod_norm;
            scal_unorm += Vel[kk] * Normal[kk];
        }
        double prod_u_t1 = 0.;
        for (int kk = 0; kk < DIMENSION; kk++) {
            Tang1[kk]  = (Vel[kk] - scal_unorm * Normal[kk]) / mod_vel;
            prod_u_t1 += Vel[kk] * Tang1[kk];
        }
        // check if tang1 has same direction of velocity field
        if (prod_u_t1 < 0.)
            for (int kk = 0; kk < DIMENSION; kk++) {
                Tang1[kk] = -Tang1[kk];
            }
        if (_nNSdim == 3) { //
            for (int kk = 0; kk < DIMENSION; kk++) {
                int a = (kk + 1) % DIMENSION;
                int b = (kk + 2) % DIMENSION;
                Tang2[kk] = (Normal[a] * Tang1[b] - Normal[b] * Tang1[a]);
            }
        }
    }
    // -------------------------------------------------------------
    if (Type == 0) {
        Tang1[0] = - (normal[1] + normal[DIMENSION - 1] * (DIMENSION % 2));
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


#endif
#endif







