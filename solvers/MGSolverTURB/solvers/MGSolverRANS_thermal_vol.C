#include "Equations_conf.h"

// ============================================
#ifdef RANS_THERMAL_EQUATIONS // 3D-2D Energy equation
// ============================================

#include <iomanip>      // std::setprecision
#include "MGSolverRANS_thermal.h"
#include "MGFE.h"          // Mesh class,double vel_g[]


// ===============================================================================================
void /**/MGSolRANS_thermal::vol_integral (
    const int el_ndof2,
    const int el_ngauss,
    const int mode,
    int el_conn[]
)   // ==============================================================================================
{

    double xyz_g[DIMENSION], T_dxg[DIMENSION];
    double rhocp = 1.;

    double ElemArea = 0.;

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------

    for ( int dim = 0; dim < _nTKdim; dim++ ) {
        _AlphaTurbDxg[dim] = 0.;
    }

    const int nEquations = _EquationNodes.size();

    for ( int qp = 0; qp < el_ngauss; qp++ ) { // INTERPOLATION OVER GAUSS NODE
        // shape functions at gaussian points -----------------------------------
        double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[_nTKdim - 1][qp];  // weight
        ElemArea += JxW_g2;

        // test functions
        _fe[2]->get_phi_gl_g ( _nTKdim, qp, _phi_g[2] );           // shape funct
        _fe[2]->get_dphi_gl_g ( _nTKdim, qp, _InvJac2, _dphi_g[2] ); // global coord deriv
        _fe[2]->get_ddphi_gl_g ( _nTKdim, qp, _InvJac2, _ddphi_g[2] ); // local second deriv

        // solution interpolation
        interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2] ); // quadratic

        // DERIVATIVES OF OLD SOLUTION
        interp_el_gdx ( _data_eq[2].ub, _FF_idx[KTT_F], 1, _dphi_g[2], el_ndof2, _KH_der );
        interp_el_gdx ( _data_eq[2].ub, _FF_idx[KTT_F] + 1, 1, _dphi_g[2], el_ndof2, _WH_der );
        interp_el_gddx ( _data_eq[2].ub, _FF_idx[KTT_F], 1, _ddphi_g[2], el_ndof2, _KH_2der );
        interp_el_gddx ( _data_eq[2].ub, _FF_idx[KTT_F] + 1, 1, _ddphi_g[2], el_ndof2, _WH_2der );

        if ( _AxiSym == 1 ) {
            interp_el_sol ( _xx_qnds, 0, _nTKdim, _phi_g[2], el_ndof2, xyz_g );
            JxW_g2 *= xyz_g[0];
        }

        // THERMAL TURBULENCE PRODUCTION
        _sT = 1.e-40;
        interp_el_gdx ( _data_eq[2].ub, _FF_idx[T_F], 1, _dphi_g[2], el_ndof2, T_dxg );

        for ( int idim = 0; idim < DIMENSION; idim++ ) {
            _sT += T_dxg[idim] * T_dxg[idim];
        }

        _alpha_turb =  max ( _ub_g[2][_FF_idx[ALPHA_T]], 0. );

        _Tkappa_g[0] = _ub_g[2][_FF_idx[KTT_F]];
        _Tkappa_g[1] = _ub_g[2][_FF_idx[KTT_F] + 1];
        _kappa_g[0]  = _ub_g[2][_FF_idx[K_F]];
        _kappa_g[1]  = _ub_g[2][_FF_idx[K_F] + 1];

        // ALPHA_TURB INTERPOLATION AND DERIVATIVES
        double mod2_vel = 1.e-20;
        VelOnGaussAndTensorModulus ( mod2_vel, el_ndof2 );

        if ( _RANS_t_parameter._InterpolatedAlphaTurb == 1 ) {
            interp_el_gdx ( _data_eq[2].ub, _FF_idx[ALPHA_T], 1, _dphi_g[2], el_ndof2, _AlphaTurbDxg );
            for ( int l = 0; l < DIMENSION; l++ ) {
                _AlphaTurbDxg[l] *= _IRe * _InvSigma;
            }
        }

        if ( _FF_idx[DIST] > -1 ) {
            _y_dist = _ub_g[2][_FF_idx[DIST]];
        }

        // function must be called using _kappa_g and not values from non linear iterations
        _mgutils._TurbParameters->CalcThermTurSourceAndDiss ( _kappa_g, _Tkappa_g, _y_dist, _sP, _sT, _alpha_turb, _source, _diss, _mecc_term );
        _alpha_eff = _alpha + max ( _ub_g[2][_FF_idx[ALPHA_T]], 0. ) * _InvSigma * _IRe;
        CalcSourceAndDiss ( );

//
//  Stabilized finite element method for heat transfer and turbulent flows inside industrial furnaces
//  On an improved unusual stabilized finite element method for the advective-reactive-diffusive equation
//
        double tauc = 0., f_upwind = 0., h_eff = 1.e-21;
        double ParVel[DIMENSION];

        for ( int i = 0; i < DIMENSION; i++ ) {
            ParVel[i] = 0.;
        }

        if ( _SUPG || _UPWIND || _SCPG ) {
            CalcSUPG ( h_eff, f_upwind, mod2_vel, _Vel_g );
        }

        const double SourceTerms = ( _explicit_source[_dir] - _explicit_diss[_dir] );
        const double TimeDer = ( 1-_SolveSteady ) * _Tkappa_g[_dir] / _dt ;

        /// d) Local (element) assemblying energy equation
        // =====================================================================================================================
        for ( int nEq = 0; nEq < nEquations; nEq++ ) { // LOOP OVER TEST FUNCTIONS PHI_I
            const int i = _EquationNodes[nEq];

            const double phii_g = _phi_g[2][i];
            double Phi_supg = 0.;

//             NUMERICAL STABILIZATION FOR INTERIOR ELEMENTS
            if ( /* _bc_vol[i] == 11 &&*/  _SUPG ) {
                Phi_supg = CalcPhiSupg ( i, _Vel_g, ParVel, tauc, f_upwind, _implicit_diss[_dir], el_ndof2 );
            }

            // RIGHT HAND SIDE -----------------------------------------------------
            _FeM ( i ) += JxW_g2 * ( phii_g + Phi_supg ) * ( TimeDer + SourceTerms );

            double CylLap=0.;

            for ( int j = 0; j < el_ndof2; j++ ) { // INTERPOLATION OVER ELEMENT NODES
                const double phij_g = _phi_g[2][j];
                CalcAdvectiveAndDiffusiveTerms ( i, j, el_ndof2, f_upwind );

                if ( _AxiSym == 1 ) {
                    CylLap = _alpha_eff * ( /* phij_g/xyz_g[0] -*/ _dphi_g[2][j] ) / xyz_g[0];
                }

                // LEFT HAND SIDE -----------------------------------------------------
                double TimeDer = ( 1-_SolveSteady ) * ( phii_g + Phi_supg ) * phij_g / _dt ; // time term

                double Eq_on_G = ( phii_g + Phi_supg ) * (
                                     + _Adv                               // advection term
                                     - _Cross[_dir]                       // Cross diffusion, for w
                                     - _Log_Cross[_dir]                   // Cross term for logarithmic models
                                     + _implicit_diss[_dir] * phij_g      // k dissipation
                                 )
                                 + (
                                     _Lap                                 // diffusion
                                     - Phi_supg * ( _LapSupg              // lap_phi*Phi_supg
                                                    + _LapMuTurb          // grad_phi*grad_nut*Phi_supg
                                                    + CylLap
                                                  )
                                 );
                _KeM ( i, j ) += JxW_g2 * ( TimeDer + Eq_on_G );
            }// END LOOP OVER ELEMENT NODES - j
        }// END LOOP OVER TEST FUNCTIONS - i
    }// END QUADRADURE LOOP - qp

    return;
}


void MGSolRANS_thermal::VelOnGaussAndTensorModulus ( double & mod2_vel, int el_ndof2 )
{

    double vel_dxg[DIMENSION * DIMENSION];
    _Vel_g[0] = _Vel_g[1] = _Vel_g[_nTKdim - 1] = 0.;

    if ( _FF_idx[NS_F] > -1 )
        for ( int idim = 0; idim < _nTKdim; idim++ ) {
            _Vel_g[idim] = _ub_g[2][_FF_idx[NS_F] + idim]; // velocity field
            mod2_vel += _Vel_g[idim] * _Vel_g[idim];
        }

    mod2_vel = sqrt ( mod2_vel );
    _sP = 1.e-20;
    interp_el_gdx ( _data_eq[2].ub, _FF_idx[NS_F], _nTKdim, _dphi_g[2], el_ndof2, vel_dxg );

    for ( int idim = 0; idim < DIMENSION; idim++ )
        for ( int jdim = 0; jdim < DIMENSION; jdim++ ) {
            const double sss = ( vel_dxg[idim + jdim * DIMENSION] + vel_dxg[jdim + idim * DIMENSION] );
            _sP += sss * sss;
        }

    if ( _AxiSym == 1 ) {
        double xyz_g[DIMENSION];
        interp_el_sol ( _xx_qnds, 0, _nTKdim, _phi_g[2], el_ndof2, xyz_g );
        _sP += _Vel_g[0] * _Vel_g[0] / ( xyz_g[0] * xyz_g[0] );
    }

    return;
}

void MGSolRANS_thermal::CalcAdvectiveAndDiffusiveTerms ( int i, int j, int el_ndof2, double f_upwind )
{

    _Adv = _Lap = _LapSupg = _LapMuTurb = _Cross[0] = _Cross[1] = _Log_Cross[0] = _Log_Cross[1] = 0.;

    for ( int idim = 0; idim <  _nTKdim; idim++ ) { // LOOP OVER SPACE COMPONENTS
        const double  dphiidxg = _dphi_g[2][i + idim * el_ndof2];
        const double  dphijdxg = _dphi_g[2][j + idim * el_ndof2];
        _Adv += _Vel_g[idim] * dphijdxg;                                                      // advection
        _Lap += _alpha_eff * dphijdxg * dphiidxg;                                             // diffusion
    }

    return;
}



double MGSolRANS_thermal::CalcPhiSupg ( int i, double vel_g[], double ParVel[], double tauc, double f_upwind, double implicit_diss, int el_ndof2 )
{
    double SkewOperator, SymmOperator;
    double SAdv = 0., SDiff = 0., VelModulus = 0., SReact = 0., LapTurb = 0.;
    double ShockCapturing = 0.;
    SymmOperator = 0.;
    double Phi_supg = 0.;
    const double phii_g = _phi_g[2][i];

    double VelForSupg[DIMENSION];
    VelocityForSUPG ( VelModulus, vel_g, VelForSupg );

    for ( int idim = 0; idim < _nTKdim; idim++ ) {
        const double  dphiidxg = _dphi_g[2][i + idim * el_ndof2];
        SAdv     += VelForSupg[idim] * dphiidxg;

        if ( _SCPG ) {
            ShockCapturing += tauc * ParVel[idim] * dphiidxg;
        }

        LapTurb += _KH_2der[idim * _nTKdim + idim];
        LapTurb += _dir * ( _KH_2der[idim * _nTKdim + idim] + _WH_2der[idim * _nTKdim + idim] );
    }// end loop over space dimension

    if ( _ModifiedSupg ) {
        if ( !_SolveSteady ) {
            SReact += ( 1. / _dt ) * phii_g ;
        }

        SReact  += implicit_diss * phii_g ;
        SDiff   *= _alpha_eff;
        LapTurb *= _alpha_eff;
    }

    SymmOperator = SReact - SDiff + LapTurb;
    SkewOperator = - SAdv;
    Phi_supg     = - f_upwind * ( SkewOperator + _theta * SymmOperator ) + ShockCapturing;

    return Phi_supg;
}

void MGSolRANS_thermal::CalcSUPG (
    double & h_eff,
    double & f_upwind,
    double mod2_vel,
    double vel_g[]
)   // NUMERICAL STABILIZATION - UPWIND AND SUPG ===============================
{
    const int  el_ndof2 = _fe[2]->_NoShape[_nTKdim - 1];

    double VEL[DIMENSION];
    double vel_modulus = 1.e-10;

    VelocityForSUPG ( vel_modulus, vel_g, VEL );

    for ( int i = 0; i < el_ndof2; i++ ) {
        double hh = 1.e-20;

        for ( int idim = 0; idim < _nTKdim; idim++ ) {
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
    const double Pe_h     = 0.5 * vel_modulus * h_eff / ( _alpha_eff );
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
            const double Pe1 = 6 * _alpha_eff / ( sigma * h_eff * h_eff );
            const double Pe2 = vel_modulus * h_eff / ( 3.*_alpha_eff );
            const double xi1 = ( Pe1 > 1. ) ? Pe1 : 1.;
            const double xi2 = ( Pe2 > 1. ) ? Pe2 : 1.;
            double f_up2;

            if ( _SolveSteady ) {
                f_up2 = h_eff * h_eff / ( ( implicit_diss ) * h_eff * h_eff * xi1 + 6.*_alpha_eff * xi2 );
            } else {
                f_up2 = h_eff * h_eff / ( sigma * h_eff * h_eff * xi1 + 6.*_alpha_eff * xi2 );
            }

            f_upwind = f_up2;
        }
        {
            // BOB   gamma = sigma/nueff
            double LapK = 0.;

            for ( int i = 0; i < _nTKdim; i++ ) {
                LapK += _dir * ( 2.*_KH_2der[i * _nTKdim + i] + _WH_2der[i * _nTKdim + i] ) + ( 1. - _dir ) * _KH_2der[i * _nTKdim + i];
            }

            LapK *= _alpha_eff;

            const double lambda = sigma - 0.5 * LapK;
            const double betak = lambda / ( h_eff * h_eff * sigma * sigma + 3.*_alpha_eff * lambda );
            f_upwind = betak * h_eff * h_eff / 3.;
        }
    }// -----------------------------------------------------------------------------------------

//
//   f_upwind = 1/(sqrt(9.*(4.*_alpha_eff/(h_eff*h_eff))*(4.*_alpha_eff/(h_eff*h_eff)) + 4.*mod2_vel*mod2_vel/(h_eff*h_eff)));

    return;
}//==========================================================================

void MGSolRANS_thermal::VelocityForSUPG (
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
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
