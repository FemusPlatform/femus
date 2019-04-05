#include "Equations_conf.h"

// ============================================
#ifdef RANS_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverRANS.h"       // Navier-Stokes class header file
#include "MGFE.h"          // Mesh class,double vel_g[]
#include <iomanip>

// ===============================================================================================
void MGSolRANS::vol_integral (
    const int el_ndof2,
    const int el_ngauss,
    double xx_qnds[],
    int ElId,
    int el_conn[]
)   // ==============================================================================================
{

    double xyz_g[DIMENSION], vel_dxg[DIMENSION * DIMENSION], x_2ts_g[1];

    const double T_0ts = ( _TimeDer == 2 ) ? 1.5 : 1.;
    const double T_1ts = ( _TimeDer == 2 ) ? 2.  : 1.;
    const double T_2ts = ( _TimeDer == 2 ) ? 0.5 : 0.;

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------

    for ( int qp = 0; qp < el_ngauss; qp++ ) { // INTERPOLATION OVER GAUSS NODE
        // shape functions at gaussian points -----------------------------------
        double det2      = _fe[2]->Jac ( qp, xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2    = det2 * _fe[2]->_weight1[_nTdim - 1][qp];   // weight

        _fe[2]->get_phi_gl_g ( _nTdim, qp, _phi_g[2] );                 // shape funct
        _fe[2]->get_dphi_gl_g ( _nTdim, qp, _InvJac2, _dphi_g[2] );     // global coord deriv
        _fe[2]->get_ddphi_gl_g ( _nTdim, qp, _InvJac2, _ddphi_g[2] );   // local second deriv

        interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2] ); // quadratic
        interp_el_sol ( _x_2ts, 0, 1, _phi_g[2], el_ndof2, x_2ts_g );

        // DERIVATIVES OF OLD SOLUTION
        interp_el_gdx ( _data_eq[2].ub, _FF_idx[K_F],      1, _dphi_g[2], el_ndof2, _T_dxg[0] );
        interp_el_gdx ( _data_eq[2].ub, _FF_idx[K_F] + 1,  1, _dphi_g[2], el_ndof2, _T_dxg[1] );
        interp_el_gddx ( _data_eq[2].ub, _FF_idx[K_F],     1, _ddphi_g[2], el_ndof2, _T_2dxg[0] );
        interp_el_gddx ( _data_eq[2].ub, _FF_idx[K_F] + 1, 1, _ddphi_g[2], el_ndof2, _T_2dxg[1] );

        interp_el_gdx ( _data_eq[2].ub, _FF_idx[MU_T], 1, _dphi_g[2], el_ndof2, _MuTurbDxg );
        for ( int l = 0; l < DIMENSION; l++ ) {
            _MuTurbDxg[l] *= _IRe * _InvSigma;
        }

        if ( _AxiSym == 1 ) {
            interp_el_sol ( _xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g );
            JxW_g2 *= xyz_g[0];
        }

        _kappa_g[0] = _ub_g[2][_FF_idx[K_F]];
        _kappa_g[1] = _ub_g[2][_FF_idx[K_F] + 1];
        
        double mod2_vel;
        // this function fills _Vel_g[] array and velocity tensor modulus _sP
        VelOnGaussAndTensorModulus ( mod2_vel, el_ndof2 );
        // this function calculates the source and diss terms for this equation
        CalcSourceAndDiss ( el_ndof2 );

        
        double tauc = 0., f_upwind = 0., h_eff = 1.e-21;
        double ParVel[DIMENSION] = {0., 0.};
        ParVel[DIMENSION-1]=0.;


        // STABILIZATION TERMS =====================================================
        f_upwind = CalcFUpwind ( _Vel_g, _dphi_g[2], _nueff, _nTdim, el_ndof2 );

        if ( ( _SUPG && _ModifiedSupg ) || _SCPG ) {
            CalcSUPG ( h_eff, f_upwind, mod2_vel, _Vel_g );
        }

        /// d) Local (element) assemblying energy equation
        double CylLap = 0, Phi_supg = 0., kem_ij;

        // =====================================================================================================================
        for ( int nEq=0; nEq<_EquationNodes.size(); nEq ++ ) {
            const int i = _EquationNodes[nEq];
            const double phii_g = _phi_g[2][i];
            Phi_supg = 0.;

            // NUMERICAL STABILIZATION FOR INTERIOR NODES
            Phi_supg = CalcPhiSupg ( i, _Vel_g, ParVel, tauc, f_upwind, _implicit_diss[_dir], el_ndof2 );

            // RIGHT HAND SIDE -----------------------------------------------------
            double SourceTerms = ( _explicit_source[_dir] - _explicit_diss[_dir] );
            double TimeDer = ( 1-_SolveSteady ) * ( T_1ts*_kappa_g[_dir] - T_2ts*x_2ts_g[0] ) / _dt ;
            _FeM ( i ) += JxW_g2 * ( phii_g + Phi_supg ) * ( SourceTerms + TimeDer);            

            for ( int j = 0; j < el_ndof2; j++ ) { // INTERPOLATION OVER ELEMENT NODES
                const double phij_g = _phi_g[2][j];
                // this function calculates the matrix coefficients -> advection, diffusion, cross diffusion
                CalcAdvectiveAndDiffusiveTerms ( i, j, el_ndof2, f_upwind );

                if ( _AxiSym == 1 ) {
                    CylLap = _nueff * ( /* phij_g/xyz_g[0] - */_dphi_g[2][j] ) / xyz_g[0];
                }

                double MatTimeDer = ( 1-_SolveSteady ) * T_0ts * ( phii_g + Phi_supg ) * phij_g / _dt;   // time term

                kem_ij = ( phii_g + Phi_supg ) * (
                             + _Adv                                 // advection term
                             - _Cross[_dir]                         // Cross diffusion, for wD
                             - _Log_Cross[_dir]                     // Cross term for logarithmic models
                             + _implicit_diss[_dir] * phij_g        // k dissipation
                         )
                         + (
                             _Lap                                   // diffusion
                             - Phi_supg * (
                                _LapSupg                   // lap_phi*Phi_supg
                               +_LapMuTurb                 // grad_phi*grad_nut*Phi_supg
                               + CylLap
                             )
                         );

                _KeM ( i, j ) += JxW_g2 * ( MatTimeDer + kem_ij );
            }// END LOOP OVER ELEMENT NODES - j
        }// END LOOP OVER TEST FUNCTIONS - i
    }// END QUADRADURE LOOP - qp


    return;
}


void MGSolRANS::wall_element (
    const int el_ndof2,
    const int el_ngauss,
    double xx_qnds[],
    int ElId,
    int el_conn[]
)   // ==============================================================================================
{

    double xyz_g[DIMENSION], vel_dxg[DIMENSION * DIMENSION], x_2ts_g[1];

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------

    // wall normal direction
    double grad[DIMENSION] = {0., 0.};
    grad[DIMENSION -1] = 0.;
    double area = 0.;
    for ( int qp = 0; qp < el_ngauss; qp++ ) {
        double det2      = _fe[2]->Jac ( qp, xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2    = det2 * _fe[2]->_weight1[_nTdim - 1][qp];   // weight
        area += JxW_g2;

        _fe[2]->get_dphi_gl_g ( _nTdim, qp, _InvJac2, _dphi_g[2] ); // global coord deriv

        double der_dist[DIMENSION] = {0., 0.};
        der_dist[DIMENSION -1] = 0.;


        for ( int n = 0; n < el_ndof2; n++ ) {
            for ( int dim = 0; dim<DIMENSION; dim++ ) {
                der_dist[dim] += JxW_g2 * _data_eq[2].ub[_FF_idx[DIST]*el_ndof2 + n] * _dphi_g[2][n + dim*el_ndof2];
            }
        }
        for ( int dim = 0; dim<DIMENSION; dim++ ) {
            grad[dim] += JxW_g2 * der_dist[dim];
        }
    }
    double norm[DIMENSION];
    double mod = 0;
    for ( int dim = 0; dim<DIMENSION; dim++ ) {
        mod += grad[dim]*grad[dim];
    }
    mod = sqrt ( mod ) + 1.e-10;
    for ( int dim = 0; dim<DIMENSION; dim++ ) {
        norm[dim] = -grad[dim]/mod;
    }


    // calculation of utau from mid cell point values
    double vel_mod =0., vel_norm = 0., vel_bound=0.;

    for ( int dim = 0; dim < _nTdim; dim ++ ) {
        double vel_d = _data_eq[2].ub[ ( _FF_idx[NS_F] + dim ) * NDOF_FEM + NDOF_FEM - 1];
        vel_mod  += vel_d * vel_d;
        vel_norm += vel_d * norm[dim];
    }

    vel_bound = sqrt ( vel_mod - vel_norm * vel_norm ) + 1.e-10;

    double WallDist = _data_eq[2].ub[_FF_idx[DIST]*el_ndof2 + el_ndof2 -1];

    double utau  = _mgutils._TurbParameters->CalcUtau ( vel_bound, WallDist );

    double kappa_omega[2];
    for ( int i=0; i<el_ndof2; i++ ) {
        double dist = _data_eq[2].ub[_FF_idx[DIST]*el_ndof2 + i];
        if ( dist < 1.e-8 ) {
            dist = WallDist;
        }
        _mgutils._TurbParameters->DynTurNearWallValues ( kappa_omega[0], kappa_omega[1], dist, utau );

        _KeM ( i, i ) = area;
        _FeM ( i )    = area*kappa_omega[_dir];
    }

    return;
}


void MGSolRANS::CalcSCPG (
    double ParVel[],
    double & tauc,
    double DefSource,
    double vel_g[],
    double mod2_vel,
    double h_eff
)   //Stabilized finite element method for heat transfer and turbulent flows inside industrial furnaces
{

    double GradMod = 0., ScalProd = 0., ParVelMod = 0.;

    for ( int i = 0; i < DIMENSION; i++ ) {
        GradMod += _dir * _T_dxg[1][i] * _T_dxg[1][i] + ( 1 - _dir ) * _T_dxg[0][i] * _T_dxg[0][i] ;
        ScalProd += vel_g[i] * ( _dir * _T_dxg[1][i] + ( 1 - _dir ) * _T_dxg[0][i] );
    }

    if ( GradMod < 1.e-10 ) {
        tauc = 0.;
    } else {
        GradMod = sqrt ( GradMod ) ;

        for ( int i = 0; i < DIMENSION; i++ ) {
            ParVel[i] = ScalProd * ( _dir * _T_dxg[1][i] + ( 1 - _dir ) * _T_dxg[0][i] ) / GradMod;
            ParVelMod += ParVel[i] * ParVel[i] + 1.e-20;
        }

        ParVelMod = sqrt ( ParVelMod );
        tauc = 0.5 * h_eff / mod2_vel;
    }

    return;
}

void MGSolRANS::CalcSUPG (
    double & h_eff,
    double & f_upwind,
    double mod2_vel,
    double vel_g[]
)   // NUMERICAL STABILIZATION - UPWIND AND SUPG ===============================
{
    const int  el_ndof2 = _fe[2]->_NoShape[_nTdim - 1];

    double VEL[DIMENSION];
    double vel_modulus = 1.e-10;

    VelocityForSUPG ( vel_modulus, vel_g, VEL );

    for ( int i = 0; i < el_ndof2; i++ ) {
        double hh = 1.e-20;

        for ( int idim = 0; idim < _nTdim; idim++ ) {
            hh += VEL[idim] * _dphi_g[2][i + idim * el_ndof2] / vel_modulus;
        }

        h_eff += fabs ( hh );
    }

    h_eff = 2. / h_eff;

    if ( h_eff < 1.e-10 ) {
        h_eff = 1. ;
        std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";
    }

    // STANDARD SUPG
    const double Pe_h     = 0.5 * vel_modulus * h_eff / ( _nueff );
    const double a_opt    = ( 1. / tanh ( Pe_h ) - 1. / Pe_h );

    if ( a_opt > 1. ) {
        std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n";
    }

    f_upwind = 0.5 * a_opt * h_eff / ( vel_modulus );

    if ( _ModifiedSupg ) { // MODIFIED SUPG - CROSS CONTRIBUTIONS AS ADDITIONAL ADVECTION ------------
        const double implicit_diss   = _implicit_diss[_dir];
        const double sigma = implicit_diss + 1. / _dt;
        {
            // TAU FORMULATION FROM ELIE HACHEM - STABILIZED FINITE ELEMENT METHOD FOR HEAT TRANSFER AND TURBULENT FLOWS INSIDE INDUSTRIAL FURNACES
            const double Pe1 = 6 * _nueff / ( sigma * h_eff * h_eff );
            const double Pe2 = vel_modulus * h_eff / ( 3.*_nueff );
            const double xi1 = ( Pe1 > 1. ) ? Pe1 : 1.;
            const double xi2 = ( Pe2 > 1. ) ? Pe2 : 1.;
            double f_up2;

            if ( _SolveSteady ) {
                f_up2 = h_eff * h_eff / ( ( implicit_diss ) * h_eff * h_eff * xi1 + 6.*_nueff * xi2 );
            } else {
                f_up2 = h_eff * h_eff / ( sigma * h_eff * h_eff * xi1 + 6.*_nueff * xi2 );
            }

            f_upwind = f_up2;
        }
        {
            // BOB   gamma = sigma/nueff
            double LapK = 0.;

            for ( int i = 0; i < _nTdim; i++ ) {
                LapK += _dir * ( 2.*_T_2dxg[0][i * _nTdim + i] + _T_2dxg[1][i * _nTdim + i] ) + ( 1. - _dir ) * _T_2dxg[0][i * _nTdim + i];
            }

            LapK *= _nueff;

            const double lambda = sigma - 0.5 * LapK;
            const double betak = lambda / ( h_eff * h_eff * sigma * sigma + 3.*_nueff * lambda );
            f_upwind = betak * h_eff * h_eff / 3.;
        }
    }// -----------------------------------------------------------------------------------------

//
//   f_upwind = 1/(sqrt(9.*(4.*_nueff/(h_eff*h_eff))*(4.*_nueff/(h_eff*h_eff)) + 4.*mod2_vel*mod2_vel/(h_eff*h_eff)));

    return;
}//==========================================================================

double MGSolRANS::CalcPhiSupg ( int i, double vel_g[], double ParVel[], double tauc, double f_upwind, double implicit_diss, int el_ndof2 )
{
    double SkewOperator, SymmOperator;
    double SAdv = 0., SDiff = 0., VelModulus = 0., SReact = 0., LapTurb = 0.;
    double ShockCapturing = 0.;
    SymmOperator = 0.;
    double Phi_supg = 0.;
    const double phii_g = _phi_g[2][i];

    double VelForSupg[DIMENSION];
    VelocityForSUPG ( VelModulus, vel_g, VelForSupg );

    for ( int idim = 0; idim < _nTdim; idim++ ) {
        const double  dphiidxg = _dphi_g[2][i + idim * el_ndof2];
        SAdv     += VelForSupg[idim] * dphiidxg;

        if ( _SCPG ) {
            ShockCapturing += tauc * ParVel[idim] * dphiidxg;
        }

        LapTurb += _T_2dxg[0][idim * _nTdim + idim];
        LapTurb += _dir * ( _T_2dxg[0][idim * _nTdim + idim] + _T_2dxg[1][idim * _nTdim + idim] );
    }// end loop over space dimension

    if ( _ModifiedSupg ) {
        if ( !_SolveSteady ) {
            SReact += ( 1. / _dt ) * phii_g ;
        }

        SReact  += implicit_diss * phii_g ;
        SDiff   *= _nueff;
        LapTurb *= _nueff;
    }

    SymmOperator = SReact - SDiff + LapTurb;
    SkewOperator = - SAdv;
    Phi_supg     = - f_upwind * ( SkewOperator + _theta * SymmOperator ) + ShockCapturing;

    return Phi_supg;
}




void MGSolRANS::CalcAdvectiveAndDiffusiveTerms ( int i, int j, int el_ndof2, double f_upwind )
{

    _Adv = _Lap = _LapSupg = _LapMuTurb = _Cross[0] = _Cross[1] = _Log_Cross[0] = _Log_Cross[1] = 0.;

    for ( int idim = 0; idim <  _nTdim; idim++ ) { // LOOP OVER SPACE COMPONENTS
        const double  dphiidxg = _dphi_g[2][i + idim * el_ndof2];
        const double  dphijdxg = _dphi_g[2][j + idim * el_ndof2];
        _Adv          += _Vel_g[idim] * dphijdxg;                                                      // advection
        _Lap          += _nueff * dphijdxg * dphiidxg;                                                // diffusion
    }

    return;
}


void MGSolRANS::CalcSourceAndDiss ( int el_ndof2 )
{

    double Div_g = 0.;

    double MuTurb[1], WallOnG[1];
    MuTurb[0] = 0.;

    for ( int k = 0; k < _nTdim; k++ ) {
        _MuTurbDxg[k] = 0.;
    }

    if ( _InterpolatedMuTurb == 1 ) {
        interp_el_gdx ( _NodeMuTurb, 0, 1, _dphi_g[2], el_ndof2, _MuTurbDxg );
        interp_el_sol ( _NodeMuTurb, 0, 1, _phi_g[2], el_ndof2, MuTurb ); // muturb interpolated from the node values
        interp_el_sol ( _NodeWallDist, 0, 1, _phi_g[2], el_ndof2, WallOnG );
        _y_dist = WallOnG[0];

        for ( int l = 0; l < DIMENSION; l++ ) {
            _MuTurbDxg[l] *= _IRe * _InvSigma;
        }
    }

    // function must be called using _kappa_g and not values from non linear iterations
    _mgutils._TurbParameters->CalcDynTurSourceAndDiss ( _kappa_g, _y_dist, _sP, _mu_turb, _source, _diss, Div_g );  // point wise source and diss
    _mu_turb = max ( 0., _mu_turb );

    if ( _InterpolatedMuTurb == 1 ) {
        _mu_turb = max ( 0., MuTurb[0] );
    }

    _explicit_source[0] = _source[0];
    _explicit_diss[0]   = _diss[0];
    _implicit_diss[0]   = 0.;
    _explicit_source[1] = _source[1];
    _explicit_diss[1]   = _diss[1];
    _implicit_diss[1]   = 0.;

    return;
}

void MGSolRANS::VelOnGaussAndTensorModulus ( double & mod2_vel, int el_ndof2 )
{

    double vel_dxg[DIMENSION * DIMENSION];
//             {
    // VELOCITY FIELD, VELOCITY MODULUS AND TURBULENCE SOURCE ==================
    _Vel_g[0] = _Vel_g[1] = _Vel_g[_nTdim - 1] = 0.;

    if ( _FF_idx[NS_F] > -1 )
        for ( int idim = 0; idim < _nTdim; idim++ ) {
            _Vel_g[idim] = _ub_g[2][_FF_idx[NS_F] + idim]; // velocity field
            mod2_vel += _Vel_g[idim] * _Vel_g[idim];
        }

    mod2_vel = sqrt ( mod2_vel );
    _sP = 1.e-20;
    interp_el_gdx ( _data_eq[2].ub, _FF_idx[NS_F], _nTdim, _dphi_g[2], el_ndof2, vel_dxg );

    for ( int idim = 0; idim < DIMENSION; idim++ )
        for ( int jdim = 0; jdim < DIMENSION; jdim++ ) {
            const double sss = ( vel_dxg[idim + jdim * DIMENSION] + vel_dxg[jdim + idim * DIMENSION] );
            _sP += sss * sss;
        }

    if ( _AxiSym == 1 ) {
        double xyz_g[DIMENSION];
        interp_el_sol ( _xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g );
        _sP += _Vel_g[0] * _Vel_g[0] / ( xyz_g[0] * xyz_g[0] );
    }

    return;
}

void MGSolRANS::VelocityForSUPG (
    double & mod2_vel,
    double vel_g[],
    double VEL[]
)   // NUMERICAL STABILIZATION - UPWIND AND SUPG ===============================
{

    for ( int i = 0; i < DIMENSION; i++ ) {
        VEL[i]  = vel_g[i];
        mod2_vel += VEL[i] * VEL[i];
    }

    mod2_vel = sqrt ( mod2_vel );

    return;
}//==========================================================================

void MGSolRANS::vol_stab (
    const int el_ndof2,
    const int el_ngauss,
    const int mode,
    int el_conn[]
)   // ==============================================================================================
{

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------

    double x_2ts_g[1];

    double t_0_c = ( _TimeDer == 1 ) ? 1.:1.5;
    double t_1_c = ( _TimeDer == 1 ) ? 1.:2.;
    double t_2_c = ( _TimeDer == 1 ) ? 0.:0.5;

    for ( int qp = 0; qp < el_ngauss; qp++ ) { // INTERPOLATION OVER GAUSS NODE
        // shape functions at gaussian points -----------------------------------
        double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];  // weight

        _fe[2]->get_phi_gl_g ( _nTdim, qp, _phi_g[2] );           // shape funct
        interp_el_sol ( _data_eq[2].ub, _FF_idx[NS_F], _nTdim, _phi_g[2], el_ndof2, _ub_g[2] ); // quadratic
        interp_el_sol ( _data_eq[2].ub, _FF_idx[K_F], 2, _phi_g[2], el_ndof2, _ub_g[2] ); // quadratic
        interp_el_sol ( _x_2ts, 0, 1, _phi_g[2], el_ndof2, x_2ts_g );

        _kappa_g[_dir] = _ub_g[2][_FF_idx[K_F]+_dir];

        /// d) Local (element) assemblying energy equation
        // =====================================================================================================================
        for ( int i = 0; i < el_ndof2; i++ ) { // LOOP OVER TEST FUNCTIONS PHI_I

            if ( _bc_el[i] != 0 ) {
                const double phii_g = _phi_g[2][i];
                double dtxJxW_g = JxW_g2; // area with bc and weight

                double TimeDer = phii_g * ( t_1_c * _kappa_g[_dir] - t_2_c * x_2ts_g[0] ) * 10. / _dt;
                _FeM ( i ) += dtxJxW_g * TimeDer;

                for ( int j = 0; j < el_ndof2; j++ ) { // INTERPOLATION OVER ELEMENT NODES
                    const double phij_g = _phi_g[2][j];
                    // LEFT HAND SIDE -----------------------------------------------------
                    double TimeDer = ( phii_g ) * phij_g * t_0_c * 10. / _dt;  // time term
                    _KeM ( i, j ) += dtxJxW_g * TimeDer ;
                }// END LOOP OVER ELEMENT NODES - j

            }// END LOOP OVER TEST FUNCTIONS - i
        }// END IF _WALL_ELEMENT < 2
    }// END QUADRADURE LOOP - qp

    return;
}


#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;


