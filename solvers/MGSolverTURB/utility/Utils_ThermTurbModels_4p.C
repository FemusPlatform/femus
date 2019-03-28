#include "Utils_ThermTurbModels_4p.h"
#include <iostream>


double Therm4P::CalcAlphaTurb ( double * DynKappaAndOmega, double * ThermTurbVar, double dist )
{
    double alpha_t;

    _KappaAndOmega[0] = DynKappaAndOmega[0];
    _KappaAndOmega[1] = DynKappaAndOmega[1];
    CalcKappaHAndOmegaH ( ThermTurbVar );
    CalcDynCoefficients ( dist );

    const double omega = _KappaAndOmega[1], omegaT = _KappaHAndOmegaH[1];
    const double Rmax = 100.;
    const double R    = omega / omegaT;
    const double R_corr = ( R > Rmax ) ? Rmax / R : 1.;

    _rT      = R_corr * omega / omegaT;
    _F1t     = ( 1. - exp ( -_Rd / ( sqrt ( _IPr ) * 14. ) ) ) * ( 1. - exp ( -_Rd / 14. ) );
    _F2at    = exp ( -4.e-6 * _Rt * _Rt );
    _F2bt    = exp ( -2.e-5 * _Rt * _Rt );

    double  IPrdlT = ( 0.1 / __CMU ) * _F1t * (
                         1./__Prdl_inf                                                   /* Asymptotic contribution */
                         + 2.*_rT / ( _rT + 0.3 ) * _F2at                                /* Contribution far from wall */
                         + 1.3 * _IPr * sqrt ( 2.*_rT ) / ( pow ( _Rt, 0.75 ) ) * _F2bt  /* Near Wall contribution */
                     );

    alpha_t = fabs ( _Ret * IPrdlT ); //*_nu

    return alpha_t;
}

void Therm4P::ThermTurInitKappaAndOmegaValues(double *ThermTurbVariables, double vel, double diameter, double RefDeltaT)
{
    const double Re_h     = vel*diameter/_nu;
    const double In       = 0.16*pow ( Re_h,-0.125 );
    const double len      = 0.07*diameter;
    const double kappa    = 1.5* ( In*vel ) * ( In*vel );
    const double omega    = sqrt ( kappa ) /len;
    _KappaHAndOmegaH[1] = omega * _IPr;
    
    const double Trms = 0.005 * RefDeltaT;
    _KappaHAndOmegaH[0] = Trms * Trms;
    
    ConvertKHandWHtoLocal(ThermTurbVariables);
    
    return;
}

void Therm4P::ThermTurNearKappaHAndOmegaHValues(double* ThermTurInitKappaAndOmegaValues, double WallDist, double Utau)
{
    double turb[2];
    DynTurNearKappaAndOmegaValues(turb , WallDist, Utau);
    
    _KappaHAndOmegaH[0] = _IPr * _KappaAndOmega[0];
    _KappaHAndOmegaH[1] = _KappaAndOmega[1];
    
    ConvertKHandWHtoLocal(ThermTurInitKappaAndOmegaValues);
    return;
}


void Therm4P_logKHWH::CalcKappaHAndOmegaH ( double * TurbVariables )
{
    double kappah = exp ( TurbVariables[0] ), omegah = exp ( TurbVariables[1] );
    _KappaHAndOmegaH[0] = ( kappah > _LowerKappah ) ? kappah:_LowerKappah;
    _KappaHAndOmegaH[1] = ( omegah > _LowerOmegah ) ? omegah:_LowerOmegah;
    return;
}

void Therm4P_logKHWH::CalcThermSourceTerms (
    double * KappaAndOmega,
    double * ThermTurbVariables,
    double * SourceTerms,
    double * DissTerms,
    double * MeccTerms,
    double TempGradMod,
    double VelGradMod,
    double muT,
    double alphaT,
    double ydist )
{

    _KappaAndOmega[0] = KappaAndOmega[0];
    _KappaAndOmega[1] = KappaAndOmega[1];
    CalcKappaHAndOmegaH ( ThermTurbVariables );
    CalcDynCoefficients ( ydist );

    double prod_k_corr = 1.;    
    if(fabs(_KappaAndOmega[0])<_LowerKappa) prod_k_corr = 0.;
    
    const double f_exp  = ( __CD2 * ( 1. - 0.3 * exp ( -_Rt * _Rt / 42.25 ) ) - 1. ) * ( 1. - exp ( -_Rd / 5.7 ) ) * ( 1. - exp ( -_Rd / 5.7 ) );

    const double Rmax = 100.;
    const double R    = _KappaAndOmega[1] / _KappaHAndOmegaH[1];
    const double R_corr = ( R > Rmax ) ? Rmax / R : 1.;

    const double prod_k = 0.5 * VelGradMod * _nu * muT;
    const double prod_kt = _nu * alphaT * TempGradMod;

    // THERMAL SOURCE
    SourceTerms[0] = prod_kt / _KappaHAndOmegaH[0];
    SourceTerms[1] = ( __CP1 - 1. ) * prod_kt / _KappaHAndOmegaH[0];

    // THERMAL DISSIPATION
    DissTerms[0]   = __CMU * _KappaHAndOmegaH[1];
    DissTerms[1]   = ( __CD1 - 1. ) * __CMU * _KappaHAndOmegaH[1];

    // MECHANICAL DISSIPATION AND SOURCE
    MeccTerms[0]  = f_exp * __CMU * _KappaAndOmega[1] * R_corr;  // DISSIPATION
    MeccTerms[1]  = prod_k_corr*__CP2 * prod_k  / _KappaAndOmega[0];         // SOURCE
    
    return;
}


void Therm4P_logKHWH::ConvertKHandWHtoLocal(double* ThermTurbVariables)
{
    ThermTurbVariables[0] = log ( _KappaHAndOmegaH[0] );
    ThermTurbVariables[1] = log ( _KappaHAndOmegaH[1] );
    return;
}


void Therm4P_KHWH::CalcKappaHAndOmegaH ( double * TurbVariables )
{
    double kappah = ( TurbVariables[0] ), omegah = ( TurbVariables[1] );
    _KappaHAndOmegaH[0] = ( kappah > _LowerKappah ) ? kappah:_LowerKappah;
    _KappaHAndOmegaH[1] = ( omegah > _LowerOmegah ) ? omegah:_LowerOmegah;
    return;
}

void Therm4P_KHWH::CalcThermSourceTerms (
    double * KappaAndOmega,
    double * ThermTurbVariables,
    double * SourceTerms,
    double * DissTerms,
    double * MeccTerms,
    double TempGradMod,
    double VelGradMod,
    double muT,
    double alphaT,
    double ydist )
{

    _KappaAndOmega[0] = KappaAndOmega[0];
    _KappaAndOmega[1] = KappaAndOmega[1];
    CalcKappaHAndOmegaH ( ThermTurbVariables );
    CalcDynCoefficients ( ydist );

    double prod_k_corr = 1.;    
    if(fabs(_KappaAndOmega[0])<_LowerKappa) prod_k_corr = 0.;
    
    const double f_exp  = ( __CD2 * ( 1. - 0.3 * exp ( -_Rt * _Rt / 42.25 ) ) - 1. ) * ( 1. - exp ( -_Rd / 5.7 ) ) * ( 1. - exp ( -_Rd / 5.7 ) );

    const double Rmax = 100.;
    const double R    = _KappaAndOmega[1] / _KappaHAndOmegaH[1];
    const double R_corr = ( R > Rmax ) ? Rmax / R : 1.;

    const double prod_k = 0.5 * VelGradMod * _nu * muT;
    const double prod_kt = _nu * alphaT * TempGradMod;

    // THERMAL SOURCE
    SourceTerms[0] = prod_kt;
    SourceTerms[1] = ( __CP1 - 1. ) * prod_kt * _KappaHAndOmegaH[1] / _KappaHAndOmegaH[0];

    // THERMAL DISSIPATION
    DissTerms[0]   = __CMU * _KappaHAndOmegaH[1];
    DissTerms[1]   = ( __CD1 - 1. ) * __CMU * _KappaHAndOmegaH[1];

    // MECHANICAL DISSIPATION AND SOURCE
    MeccTerms[0]  = f_exp * __CMU * _KappaAndOmega[1] * R_corr;                       // DISSIPATION
    MeccTerms[1]  = prod_k_corr * __CP2 * prod_k * _KappaHAndOmegaH[1] / _KappaAndOmega[0];         // SOURCE

    return;
}


void Therm4P_KHWH::ConvertKHandWHtoLocal(double* ThermTurbVariables)
{
    ThermTurbVariables[0] = ( _KappaHAndOmegaH[0] );
    ThermTurbVariables[1] = ( _KappaHAndOmegaH[1] );
    return;
}

void Therm4P_KHEH::CalcKappaHAndOmegaH ( double * TurbVariables )
{
    double kappah = ( TurbVariables[0] ), omegah = ( TurbVariables[1] / (__CMU*TurbVariables[0]) );
    _KappaHAndOmegaH[0] = ( kappah > _LowerKappah ) ? kappah:_LowerKappah;
    _KappaHAndOmegaH[1] = ( omegah > _LowerOmegah ) ? omegah:_LowerOmegah;
    return;
}

void Therm4P_KHEH::ConvertKHandWHtoLocal(double* ThermTurbVariables)
{
    ThermTurbVariables[0] = ( _KappaHAndOmegaH[0] );
    ThermTurbVariables[1] = ( _KappaHAndOmegaH[1] * _KappaHAndOmegaH[0] * __CMU );
    return;
}

void Therm4P_KHEH::CalcThermSourceTerms (
    double * KappaAndOmega,
    double * ThermTurbVariables,
    double * SourceTerms,
    double * DissTerms,
    double * MeccTerms,
    double TempGradMod,
    double VelGradMod,
    double muT,
    double alphaT,
    double ydist )
{

    _KappaAndOmega[0] = KappaAndOmega[0];
    _KappaAndOmega[1] = KappaAndOmega[1];
    CalcKappaHAndOmegaH ( ThermTurbVariables );
    CalcDynCoefficients ( ydist );

    double prod_k_corr = 1.;    
    if(fabs(_KappaAndOmega[0])<_LowerKappa) prod_k_corr = 0.;
    
    const double f_exp  = ( __CD2 * ( 1. - 0.3 * exp ( -_Rt * _Rt / 42.25 ) ) - 1. ) * ( 1. - exp ( -_Rd / 5.7 ) ) * ( 1. - exp ( -_Rd / 5.7 ) );

    const double Rmax = 100.;
    const double R    = _KappaAndOmega[1] / _KappaHAndOmegaH[1];
    const double R_corr = ( R > Rmax ) ? Rmax / R : 1.;

    const double prod_k = 0.5 * VelGradMod * _nu * muT;
    const double prod_kt = _nu * alphaT * TempGradMod;

    // THERMAL SOURCE
    SourceTerms[0] = prod_kt;
    SourceTerms[1] = __CP1 * prod_kt * __CMU * _KappaHAndOmegaH[1];

    // THERMAL DISSIPATION
    DissTerms[0]   = __CMU * _KappaHAndOmegaH[1];
    DissTerms[1]   = __CD1 * __CMU * __CMU * _KappaHAndOmegaH[1];

    // MECHANICAL DISSIPATION AND SOURCE
    const double epsilonT = _KappaHAndOmegaH[0] * _KappaHAndOmegaH[1] * __CMU;
    MeccTerms[0]  = f_exp * __CMU * _KappaAndOmega[1] / _KappaHAndOmegaH[0];                       // DISSIPATION
    MeccTerms[1]  = prod_k_corr * __CP2 * prod_k * epsilonT / _KappaAndOmega[0];         // SOURCE

    return;
}



// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
