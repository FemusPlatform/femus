#include "TurbUtils.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <fstream>


//
using namespace std;

TurbUtils::TurbUtils(int MeshId, int nMeshes) :
 __Proc( 0 ), __MeshID(MeshId), __NMeshes(nMeshes) {
    FillModelMap();
    FillParameters(); // DEFAULT VALUES
}

TurbUtils::~TurbUtils() {}

void TurbUtils::FillParameters() {
    _BoundWallDist = 1.e-3;
    _nu = 1.e-4;
    _alpha = 1.e-4;
    _klog  = 0.;
    _wlog  = 0.;
    _khlog = 0.;
    _whlog = 0.;

    read_file();

    _klim = _wlim = _elim = _khlim = _whlim = _ehlim = 1.e-20;
    _SolveAlphaT = _SolveMuT = 0;
    
    std::cout<<"\n ========================================================= \n"
    << "\033[38;5;118m \t \t SETTING THE TURBULENCE MODEL \033[0m \n";

    _BoundWallDist = stod ( _FileMap["Wall_dist"] );

    const double rho = stod ( _FileMap["rho0"] );
    const double mu = stod ( _FileMap["mu0"] );
    const double kappa = stod ( _FileMap["kappa0"] );
    const double cp = stod ( _FileMap["cp0"] );

    _nu = mu/rho;
    _alpha = kappa/ ( rho*cp );
    __IPr  = _alpha/_nu;

//
    std::cerr<<"==================================================== \n";
    std::cerr<<" TURB_UTILS: available dynamic turbulence models: \n";
    for ( auto it = ::DynTurbModelMap.cbegin(); it != ::DynTurbModelMap.cend(); ++it )
        std::cerr <<"     "<< it->first << " " << it->second<< "\n";
    std::cerr<<" TURB_UTILS: available thermal turbulence models: \n";
    for ( auto it = ::ThermTurbModelMap.cbegin(); it != ::ThermTurbModelMap.cend(); ++it )
        std::cerr <<"     "<< it->first << " " << it->second<< "\n";
    std::cerr<<"    Chosen dynamic turb model: "<<_FileMap["RANS_dynamic"]<<std::endl;
    std::cerr<<"    Chosen themral turb model: "<<_FileMap["RANS_thermal"]<<std::endl;

    std::string DynModel = ( _FileMap["RANS_dynamic"] != "" ) ? _FileMap["RANS_dynamic"]:"default";
    int RANS_dynamic = ::DynTurbModelMap.at ( DynModel );
    switch ( RANS_dynamic ) {
    case default_dyn:
        _Nagano = 1;
        _Wilcox = _klog = _wlog = _emod = _numod = 0;
        break;
    case nagano_ke:
        _Nagano = _emod = 1;
        _Wilcox = _klog = _wlog = _numod =0;
        break;
    case nagano_kw:
        _Nagano = 1;
        _Wilcox = _klog = _wlog = _emod = _numod = 0;
        break;
    case nagano_log:
        _Nagano = _klog = _wlog = 1;
        _Wilcox = _numod = _emod =0;
        break;
    case wilcox:
        _Nagano = _klog = _wlog = _emod = _numod = 0;
        _Wilcox = 1;
        break;
    case wilcox_log:
        _Nagano = _emod = _numod = 0;
        _Wilcox = _klog = _wlog =1;
        break;
    default:
        std::cout<<"\033[1;31m\n=====================================================\n"
        <<"   TURB_UTILS: unknown dynamical turbulence model "<<_FileMap["RANS_dynamic"]
        <<"\n=====================================================\n\033[0m";
        abort();
        break;
    }

    std::string ThermModel = ( _FileMap["RANS_thermal"] != "" ) ? _FileMap["RANS_thermal"]:"default";
    int RANS_thermal = ::ThermTurbModelMap.at ( ThermModel );
    _Kays = 0;
    _4PwithKays = 0;
    
    if (_FileMap["UseKays"] != "" ) _4PwithKays = stoi(_FileMap["UseKays"]);
    
    switch ( RANS_thermal ) {
    case default_therm:
        _ehmod = _khlog = _whlog = 0;
        break;
    case nagano_keT:
        _ehmod = 1;
        _khlog = _whlog = 0;
        break;
    case nagano_kwT:
        _ehmod = _khlog = _whlog = 0;
        break;
    case nagano_logT:
        _ehmod = 0;
        _khlog = _whlog = 1;
        break;
    case kays:
        _Kays = 1;
        break;
    default:
        std::cout<<"\033[1;31m\n=====================================================\n"
        <<"   TURB_UTILS: unknown thermal turbulence model "<<_FileMap["RANS_thermal"]
        <<"\n=====================================================\n\033[0m";
        abort();
        break;
    }


    _IsFilled = true;

    if ( _FileMap ["YapCorrection"]!="" ) _YapCorr = stoi ( _FileMap["YapCorrection"] );
    else  _YapCorr = 0;
    if ( _FileMap ["DurbinConstrain"]!="" ) _Durbin = stoi ( _FileMap["DurbinConstrain"] );
    else _Durbin = 0;
    if ( _FileMap ["Park"]!="" ) _Park = stoi ( _FileMap["Park"] );
    else _Park = 1;
    if ( _FileMap ["Diameter"]!="" ) _diameter = stod ( _FileMap["Diameter"] );
    else _diameter = 0.1;
    if ( _FileMap ["AvVelocity"]!="" ) _vmid = stod ( _FileMap["AvVelocity"] );
    else _vmid = 0.1;
    if ( _FileMap ["utau"]!="" ) _InputUtau = stod ( _FileMap["utau"] );
    else _InputUtau = -1; // if negative then profiles will be calculated with default utau method

    
    if ( _FileMap ["SolveMuT"]!="" ) _SolveMuT = stoi ( _FileMap["SolveMuT"] );
    if ( _FileMap ["SolveAlphaT"]!="" ) _SolveAlphaT = stoi ( _FileMap["SolveAlphaT"] );
    
    if ( __Proc==0 ) {
        std::cerr <<" \n =====================================================\n";
        std::cerr <<"  TURB UTILS PARAMETERS   \n";
        std::cerr <<"   _YapCorr    "<<_YapCorr<<std::endl;
        std::cerr <<"   _Durbin     "<<_Durbin<<std::endl;
        std::cerr <<"   _Park       "<<_Park<<std::endl;
        std::cerr <<"   _diameter   "<<_diameter<<std::endl;
        std::cerr <<"   _vmid       "<<_vmid<<std::endl;
        std::cerr <<"   _InputUtau  "<<_InputUtau<<std::endl;
        std::cerr <<"   _Nagano      "<<_Nagano<<std::endl;
        std::cerr <<"   _Wilcox      "<<_Wilcox<<std::endl;
        std::cerr <<"   _klog        "<<_klog  <<std::endl;
        std::cerr <<"   _wlog        "<<_wlog   <<std::endl;
        std::cerr <<"   _emod        "<<_emod  <<std::endl;
        std::cerr <<"   _numod       "<<_numod  <<std::endl;
        std::cerr <<" =====================================================\n";
    }

    return;
}


void TurbUtils::FillModelMap() {

    _DynamicModel["nagano_k"]      = nagano_k    ;
    _DynamicModel["nagano_logk"]   = nagano_logk ;
    _DynamicModel["nagano_w"]      = nagano_w    ;
    _DynamicModel["nagano_logw"]   = nagano_logw ;
    _DynamicModel["nagano_e"]      = nagano_e    ;
    _DynamicModel["wilcox_k"]      = wilcox_k    ;
    _DynamicModel["wilcox_logk"]   = wilcox_logk ;
    _DynamicModel["wilcox_nut"]    = wilcox_nut  ;
    _DynamicModel["wilcox_w"]      = wilcox_w    ;
    _DynamicModel["wilcox_logw"]   = wilcox_logw ;

    _DynamicModel["nagano_ke"]  =     nagano_ke  ;
    _DynamicModel["nagano_kw"]  =     nagano_kw  ;
    _DynamicModel["nagano_log"] =     nagano_log ;
    _DynamicModel["wilcox"]     =     wilcox     ;
    _DynamicModel["wilcox_log"] =     wilcox_log ;

    _ThermalModel["natural_kh"] = natural_kh;
    _ThermalModel["logarithmic_kh"] = logarithmic_kh;
    _ThermalModel["natural_omegah"] = natural_omegah;
    _ThermalModel["logarithmic_omegah"] = logarithmic_omegah;
//     _ThermalModel["logarithmic_omegah"] = logarithmic_omegah;
    return;
}

double TurbUtils::CalcUtau ( double vel_bound, double dist ) {
    double umusk, ulog, ulin, utau, ulold, diff;
    double beta, vk;

    double yp;

    beta = 5.2;
    vk = 0.41;

    ulog = 0.;
    ulold = 11.6/vel_bound;

    // Calculation of utau through the linear relation
    ulin = sqrt ( vel_bound * _nu / dist );
    yp = ulin*dist/_nu;

    // Calculation of utau through the logarithmic relation
    diff = 100.;
    if ( yp > 5. ) {
        while ( diff > 1.e-6 ) {
            ulog = vel_bound*vk/ ( log ( exp ( beta*vk ) *dist*ulold/_nu ) );
            diff = fabs ( ulold - ulog );
            ulold = ulog;
        }
        yp = ulog*dist/_nu;
    } else {
        ulog = 0.;
    }

    // Calculation of utau through the musker relation

    if ( yp > 5. && yp < 40. ) {
        diff = 100.;
        int cont = 0;
        double umuskold = ulog;

        while ( diff > 1.e-6 ) {
            umusk = vel_bound/Musker ( dist, umuskold );
            diff = fabs ( umuskold - umusk );
            umuskold = umusk;
            cont ++;
            if ( cont > 3000 ) {
                umusk = 0.;
                break;
            }
        }
    }

    if ( yp > 5. )  {
        utau = max ( ulog,umusk );
        if ( umusk>2.*ulog ) {
            utau = ulog;
        }
    } else {
        utau = ulin;
    }

    return utau;
}

double TurbUtils::Musker ( double dist, double utau ) {
    double yplus = dist*utau/_nu;
    double vel = 5.424*atan ( ( 2.*yplus-8.15 ) /16.7 )
                 + 4.1693*log ( yplus + 10.6 )
                 - 0.8686*log ( yplus*yplus - 8.15*yplus+86 )
                 - 3.52;
    return vel;
}


double TurbUtils::CalcMuTurb ( double KappaAndOmega[], double dist, double vel_sp ) {

    if ( _numod==1 ) {
        _MuTurb = KappaAndOmega[0]/_nu;
    } else { // WILCOX AND NAGANO TURBULENCE MODELS - NOT FOR NUT EQUATION
        double kappa, omega, epsilon;

        // REAL K, OMEGA AND EPSILON CALCULATION
        kappa = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
        kappa = ( kappa>_klim ) ? kappa:_klim;
        if ( _emod==1 ) {
            epsilon = KappaAndOmega[1];
            omega   = epsilon/ ( kappa*__CMU );
        } else omega = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );
        omega = ( omega>_wlim ) ? omega:_wlim;

        // EDDY VISCOSITY CALCULATION
        __Ret   = kappa/ ( _nu*omega );            //  viscosity ratio
        __Rt    = __Ret/__CMU;                     // turbulent Reynolds number
        __Rd    = dist*sqrt ( kappa/sqrt ( __Rt ) ) /_nu;
        _MuTurb = __Ret;                           // turbulent viscosity ratio
        if ( _Nagano==1 ) {
            __fmu   = ( 1.-exp ( -1.*__Rd/14. ) ) * ( 1.-exp ( -1.*__Rd/14. ) );
            __fcorr = 1. + 5./pow ( __Rt,0.75 ) *exp ( -1.*__Rt*__Rt/40000. );
            _MuTurb *= __fcorr*__fmu;
        } else { // fmu and fcorr set to 1 -> in CalcAlphaTurb we calculate nut/__fmu*__fcorr
            __fmu = 1;
            __fcorr = 1;
        }
    }
//      std::cout<<_MuTurb*_nu<<"  "<<KappaAndOmega[0]*KappaAndOmega[0]/KappaAndOmega[1]<<"  "<<dist<<std::endl;
    return _MuTurb;
}

void TurbUtils::CalcDynTurSourceAndDiss ( double KappaAndOmega[], double dist, double vel_sp ,double &muturb, double source[2], double diss[2], double div_g ) {
    // REAL KAPPA, OMEGA AND EPSILON
    double kappa,omega,epsilon;
    kappa =  omega = epsilon =1.;
    kappa = KappaAndOmega[0];
    omega = KappaAndOmega[1];
    epsilon = KappaAndOmega[1];

    if ( _numod==0 ) kappa = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
    if ( _emod==0 )  omega = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );

    kappa = ( kappa>_klim ) ? kappa:_klim;
    omega = ( omega>_wlim ) ? omega:_wlim;
    if ( _emod==1 ) epsilon = ( epsilon>_elim ) ? epsilon:_elim;

    const double kCorr = ( _klog/kappa + ( 1.-_klog ) );
    const double wCorr = ( _wlog/omega + ( 1.-_wlog ) );

    // Values of __Rt and __Rd are assigned within TurbUtils::CalcMuTurb
    muturb =  CalcMuTurb ( KappaAndOmega,dist );
    double prod_k = 0.5 * _nu * muturb * vel_sp /*- 2.*kappa*div_g/3.*/;
    if ( _Park==1 )
        prod_k = min ( prod_k, sqrt ( 8./3. ) *kappa*sqrt ( vel_sp ) ); // k production limitation -> Park

    if ( _Nagano==1 ) { // NAGANO TURBULENCE MODEL-----------------------------------------------------------------------
        const double f_exp  = ( 1.-exp ( -1.*__Rd/ ( 3.1 ) ) ) * ( 1.-exp ( -1.*__Rd/ ( 3.1 ) ) ) * ( 1.- 0.3*exp ( -1.*__Rt*__Rt/42.25 ) );
        if ( _emod==1 ) { // KAPPA - EPSILON TURBULENCE MODEL
            diss[0]   = epsilon/kappa;                // -> it can be set implicitly in matrix
            diss[1]   = epsilon*__C20*f_exp/kappa;
            source[0] = prod_k*kCorr;
            source[1] = __C10*prod_k*epsilon/kappa;
            if ( _YapCorr==1 ) {
                double yap_term = YapTerm( KappaAndOmega[0], KappaAndOmega[1], dist );
                source[1] += yap_term;
            }
        } else { // KAPPA - OMEGA AND LOG FORMULATION TURBULENCE MODEL
            diss[0]   = __CMU*omega;
            diss[1]   = __CMU*omega* ( __C20*f_exp -1 );
            source[0] = prod_k*kCorr;
            source[1] = ( __C10-1. ) *prod_k*omega*wCorr/kappa;
            if ( _YapCorr==1 ) {
                double yap_term = YapTerm( KappaAndOmega[0], KappaAndOmega[1], dist );
                source[1] += yap_term;
            }
        }
    }//----------------------------------------------------------------------------------------------------------
    if ( _Wilcox==1 ) { // WILCOX TURBULENCE MODEL-----------------------------------------------------------------------
        // FIRST EQUATION SOURCE AND DISSIPATION
        if ( _numod==1 ) { // NUT EQUATION
            diss[0]   = __BETAN*omega;
            source[0] = 0.5*vel_sp*__AN/omega;
        } else { // KAPPA AND LOG(KAPPA) EQUATION
            diss[0]   = __BETAS*omega;
            source[0] = prod_k*kCorr;
        }
        // SECOND EQUATION SOURCE AND DISSIPATION
        diss[1]   = __BETAW*omega;
        source[1] = __AW * 0.5 * vel_sp*wCorr;
    }//----------------------------------------------------------------------------------------------------------
    return;
}



double TurbUtils::YapTerm( double & Kappa, double & Omega, double WallDist ){
    double yapterm = 0.;
    
    double kappa = ( 1.-_klog ) * Kappa + _klog*exp ( _klog*Kappa );
    double epsilon = Omega;
    double omega   = Omega;
    
    if ( _emod !=1 ){
        omega   =  ( 1.-_wlog ) * Omega + _wlog*exp ( _wlog*Omega );
        epsilon = __CMU * kappa * omega;
    }
    else{ 
        omega = epsilon / (__CMU * kappa);
    }
    
    double l_epsilon = pow ( __CMU,0.75 ) * 0.41*WallDist;
    double p1 = 0.83 * epsilon*epsilon/kappa;
    double p2 = ( kappa*sqrt ( kappa ) / ( epsilon * l_epsilon ) - 1 );
    double p3 = ( kappa*sqrt ( kappa ) / ( epsilon * l_epsilon ) );
    p3 *= p3;
    yapterm = p1*p2*p3;
    
    if ( _emod!=1 ){
        const double wCorr = ( _wlog/omega + ( 1.-_wlog ) );
        yapterm *= wCorr;
    }
    
    yapterm = max(yapterm, 0.);  
    
    return yapterm;
}


double TurbUtils::KaysPrt (double KappaAndOmega[], double dist){
    double nut = CalcMuTurb ( KappaAndOmega, dist );
    const double kays_prt = (0.85 + 0.7*__IPr/(nut + 1.e-10));
    return kays_prt;
}


double TurbUtils::CalcAlphaTurb ( double KappaAndOmega[], double TKappaAndOmega[], double dist ) {
    // REAL OMEGA AND OMEGAT
    double kappa  = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
    double omega  = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );
    double kappaT = ( 1.-_khlog ) * TKappaAndOmega[0] + _khlog*exp ( _khlog*TKappaAndOmega[0] );
    double omegaT = ( 1.-_whlog ) * TKappaAndOmega[1] + _whlog*exp ( _whlog*TKappaAndOmega[1] );

    double epsilon  = (_emod ==1) ?  KappaAndOmega[1]: omega*kappa * 0.09;
    double epsilonT = (_ehmod==1) ? TKappaAndOmega[1]: omegaT*kappaT * 0.09;

    kappa    = (kappa > 1.e-10)?    kappa:1.e-10;
    epsilon  = (epsilon > 1.e-10)?  epsilon:1.e-10;
    omega    = (omega > 1.e-10)?    omega:1.e-10;
    kappaT   = (kappaT > 1.e-10)?   kappaT:1.e-10;
    epsilonT = (epsilonT > 1.e-10)? epsilonT:1.e-10;
    omegaT   = (omegaT > 1.e-10)?   omegaT:1.e-10;

    if(_emod==1)  omega  = epsilon  / (__CMU * kappa);
    if(_ehmod==1) omegaT = epsilonT / (__CMU * kappaT);
    
    const double Rmax = 100.;
    const double R    = omega / omegaT;
    const double R_corr = (R > Rmax)? Rmax / R : 1.;
    
    
    // MU_TURB calculation
    // __fmu __fcorr __Ret __Rt __Rd
    double nut = CalcMuTurb ( KappaAndOmega, dist );
    double alphaT = 0.;
    
    if(_Kays == 1){
        const double InvKaysPrt = 1./(KaysPrt(KappaAndOmega, dist) + 1.e-10);
        alphaT = nut * InvKaysPrt;
    }
    else{

        nut = __Ret;   // nut without near wall corrections
        __rT      = R_corr* omega/omegaT;
        __F1t     = ( 1.-exp ( -__Rd/ ( sqrt ( __IPr ) *14. ) ) ) * ( 1.-exp ( -__Rd/14. ) );
        __F2at    = exp ( -4.e-6*__Rt*__Rt );
        __F2bt    = exp ( -2.e-5*__Rt*__Rt );
        
        const double InvKaysPrt = 1./(KaysPrt(KappaAndOmega, dist) + 1.e-10);    
        const double pr_inf = (1-_4PwithKays)*__Prdl_inf + _4PwithKays*InvKaysPrt;
        double  IPrdlT = ( 0.1/__CMU ) *__F1t* (
                                   __Prdl_inf                                                   /* Asymptotic contribution */
                                   + 2.*__rT/ ( __rT+0.3 ) *__F2at                              /* Contribution far from wall */
                                   + 1.3*__IPr*sqrt ( 2.*__rT ) / ( pow ( __Rt,0.75 ) ) *__F2bt /* Near Wall contribution */
                               );
        
        alphaT = fabs(nut * IPrdlT); //*_nu
    
    }
    return alphaT;
}


double TurbUtils::CalcAlphaTurb ( double KappaAndOmega[], double TKappaAndOmega[], double dist , double Nut) {
    // REAL OMEGA AND OMEGAT
    double kappa  = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
    double omega  = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );
    double kappaT = ( 1.-_khlog ) * TKappaAndOmega[0] + _khlog*exp ( _khlog*TKappaAndOmega[0] );
    double omegaT = ( 1.-_whlog ) * TKappaAndOmega[1] + _whlog*exp ( _whlog*TKappaAndOmega[1] );

    double epsilon  = (_emod ==1) ?  KappaAndOmega[1]: omega*kappa * 0.09;
    double epsilonT = (_ehmod==1) ? TKappaAndOmega[1]: omegaT*kappaT * 0.09;

    kappa    = (kappa > 1.e-10)?    kappa:1.e-10;
    epsilon  = (epsilon > 1.e-10)?  epsilon:1.e-10;
    omega    = (omega > 1.e-10)?    omega:1.e-10;
    kappaT   = (kappaT > 1.e-10)?   kappaT:1.e-10;
    epsilonT = (epsilonT > 1.e-10)? epsilonT:1.e-10;
    omegaT   = (omegaT > 1.e-10)?   omegaT:1.e-10;

    if(_emod==1)  omega  = epsilon  / (__CMU * kappa);
    if(_ehmod==1) omegaT = epsilonT / (__CMU * kappaT);
    
    const double Rmax = 100.;
    const double R    = omega / omegaT;
    const double R_corr = (R > Rmax)? Rmax / R : 1.;
    
    // MU_TURB calculation
    // __fmu __fcorr __Ret __Rt __Rd
    double nut = CalcMuTurb ( KappaAndOmega, dist );
    double alphaT = 0.;
    
    if(_Kays == 1){
        const double InvKaysPrt = 1./(KaysPrt(KappaAndOmega, dist) + 1.e-10);
        alphaT = Nut /(0.85 + 0.7*__IPr/(Nut + 1.e-10));
    }
    else{

        nut = __Ret;   // nut without near wall corrections
        __rT       = R_corr * omega/omegaT;
        __F1t     = ( 1.-exp ( -__Rd/ ( sqrt ( __IPr ) *14. ) ) ) * ( 1.-exp ( -__Rd/14. ) );
        __F2at    = exp ( -4.e-6*__Rt*__Rt );
        __F2bt    = exp ( -2.e-5*__Rt*__Rt );
        
        const double InvKaysPrt = 1./(KaysPrt(KappaAndOmega, dist) + 1.e-10);    
        const double pr_inf = (1-_4PwithKays)*__Prdl_inf + _4PwithKays*InvKaysPrt;
        double  IPrdlT = ( 0.1/__CMU ) *__F1t* (
                                   __Prdl_inf                                                   /* Asymptotic contribution */
                                   + 2.*__rT/ ( __rT+0.3 ) *__F2at                              /* Contribution far from wall */
                                   + 1.3*__IPr*sqrt ( 2.*__rT ) / ( pow ( __Rt,0.75 ) ) *__F2bt /* Near Wall contribution */
                               );
        
        alphaT = fabs(nut * IPrdlT); //*_nu
    
    }
    return alphaT;
}

void TurbUtils::CalcThermTurSourceAndDiss ( double KappaAndOmega[], // Dynamic Turbulence
        double TKappaAndOmega[],  // Thermal Turbulence
        double dist,              // Wall Distance
        double sp,                // Squared value of velocity derivatives
        double st,                // Modulus of temperature gradient
        double &alphaturb,        // Eddy thermal diffusivity
        double source[],          // Source terms
        double diss[],            // Dissipation terms
        double meccterm[] ) {     // Mechanical contribution

    double kappa  = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
    double omega  = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );
    double kappaT = ( 1.-_khlog ) * TKappaAndOmega[0] + _khlog*exp ( _khlog*TKappaAndOmega[0] );
    double omegaT = ( 1.-_whlog ) * TKappaAndOmega[1] + _whlog*exp ( _whlog*TKappaAndOmega[1] );

    double epsilon  = (_emod ==1) ?  KappaAndOmega[1]: omega*kappa * 0.09;
    double epsilonT = (_ehmod==1) ? TKappaAndOmega[1]: omegaT*kappaT * 0.09;

    double prod_k_corr = 1.;    
    kappa    = (kappa > 1.e-10)?    kappa:1.e-10;
    if(fabs(kappa)<1.e-8) prod_k_corr = 0.;
    
    epsilon  = (epsilon > 1.e-10)?  epsilon:1.e-10;
    
    kappaT   = (kappaT > 1.e-10)?   kappaT:1.e-10;
    epsilonT = (epsilonT > 1.e-10)? epsilonT:1.e-10;
    omegaT   = (omegaT > 1.e-10)?   omegaT:1.e-10;

    if(_emod==1)  omega  = epsilon  / (__CMU * kappa);
    if(_ehmod==1) omegaT = epsilonT / (__CMU * kappaT);  
    omega    = fabs(omega);
    
    const double Rmax = 100.;
    const double R    = omega / omegaT;
    const double R_corr = (R > Rmax)? Rmax / R : 1.;
    
    alphaturb = CalcAlphaTurb ( KappaAndOmega, TKappaAndOmega, dist );
    const double f_exp  = ( __CD2* ( 1.-0.3*exp ( -__Rt*__Rt/42.25 ) )-1. ) * ( 1. - exp ( -__Rd/5.7 ) ) * ( 1.-exp ( -__Rd/5.7 ) );

    double muturb = __fmu*__fcorr*kappa/omega;
    muturb = (muturb>1.e-10)? muturb:1.e-10;
    if ( _Durbin==1 ) {
        if ( muturb > kappa/ ( sqrt ( sp ) *_nu ) ) {
            muturb = kappa/ ( sqrt ( sp ) *_nu );
        }
    }
    const double prod_k = 0.5*sp*muturb;
    const double prod_kt = _nu*alphaturb*st;

    if ( _ehmod==0 ) {// THERMAL KAPPA AND OMEGA TURBULENCE MODEL

        const double TkCorr = ( ( 1.-_khlog ) + _khlog/kappaT );
        const double TwCorr = ( ( 1.-_whlog ) + _whlog/omegaT );

        // THERMAL SOURCE
        source[0] = prod_kt*TkCorr;
        source[1] = ( __CP1-1. ) *prod_kt*TwCorr*omegaT/kappaT;

        // THERMAL DISSIPATION
        diss[0]   = __CMU*omegaT;
        diss[1]   = ( __CD1-1. ) *__CMU*omegaT;

        // MECHANICAL DISSIPATION AND SOURCE
        meccterm[0]  = f_exp*__CMU*omega*R_corr;                      // DISSIPATION
        meccterm[1]  = prod_k_corr*__CP2*prod_k*omegaT*TwCorr/kappa;       // SOURCE
        
//         if ( _YapCorr==1 ) {
//             double yap_term = YapTerm( KappaAndOmega[0], KappaAndOmega[1], dist );
//             meccterm[1]  += yap_term;
//         }

    } else {// THERMAL KAPPA AND EPSILON TURBULENCE MODEL
        double kappaT = TKappaAndOmega[0];
        double epsilonT = TKappaAndOmega[1];

        // THERMAL SOURCE
        source[0] = prod_kt;
        source[1] = __CP1 *prod_kt*epsilonT/kappaT;

        // THERMAL DISSIPATION
        diss[0]   = epsilonT/kappaT;
        diss[1]   = __CD1 *__CMU*epsilonT/kappaT;

        // MECHANICAL DISSIPATION AND SOURCE
        meccterm[0]  = f_exp*__CMU*omega/kappaT;          // DISSIPATION
        meccterm[1]  = prod_k_corr*__CP2*prod_k*epsilonT/kappa;       // SOURCE
    }
    return;
}

void TurbUtils::CalcThermTurSourceAndDiss ( double KappaAndOmega[], // Dynamic Turbulence
        double TKappaAndOmega[],  // Thermal Turbulence
        double dist,              // Wall Distance
        double sp,                // Squared value of velocity derivatives
        double st,                // Modulus of temperature gradient
        double &alphaturb,        // Eddy thermal diffusivity
        double source[],          // Source terms
        double diss[],            // Dissipation terms
        double meccterm[],
        double Nut ) {     // Mechanical contribution

    double kappa  = ( 1.-_klog ) * KappaAndOmega[0] + _klog*exp ( _klog*KappaAndOmega[0] );
    double omega  = ( 1.-_wlog ) * KappaAndOmega[1] + _wlog*exp ( _wlog*KappaAndOmega[1] );
    double kappaT = ( 1.-_khlog ) * TKappaAndOmega[0] + _khlog*exp ( _khlog*TKappaAndOmega[0] );
    double omegaT = ( 1.-_whlog ) * TKappaAndOmega[1] + _whlog*exp ( _whlog*TKappaAndOmega[1] );

    double epsilon  = (_emod ==1) ?  KappaAndOmega[1]: omega*kappa * 0.09;
    double epsilonT = (_ehmod==1) ? TKappaAndOmega[1]: omegaT*kappaT * 0.09;

    double prod_k_corr = 1.;    
    kappa    = (kappa > 1.e-10)?    kappa:1.e-10;
    if(fabs(kappa)<1.e-8) prod_k_corr = 0.;
    
    epsilon  = (epsilon > 1.e-10)?  epsilon:1.e-10;
    
    kappaT   = (kappaT > 1.e-10)?   kappaT:1.e-10;
    epsilonT = (epsilonT > 1.e-10)? epsilonT:1.e-10;
    omegaT   = (omegaT > 1.e-10)?   omegaT:1.e-10;

    if(_emod==1)  omega  = epsilon  / (__CMU * kappa);
    if(_ehmod==1) omegaT = epsilonT / (__CMU * kappaT);  
    omega    = fabs(omega);
    const double Rmax = 100.;
    const double R    = omega / omegaT;
    const double R_corr = (R > Rmax)? Rmax / R : 1.;
    
    alphaturb = CalcAlphaTurb ( KappaAndOmega, TKappaAndOmega, dist, Nut );
    const double f_exp  = ( __CD2* ( 1.-0.3*exp ( -__Rt*__Rt/42.25 ) )-1. ) * ( 1. - exp ( -__Rd/5.7 ) ) * ( 1.-exp ( -__Rd/5.7 ) );

    double muturb = fabs(Nut) * _nu;

    const double prod_k = 0.5*sp*muturb;
    const double prod_kt = _nu*alphaturb*st;

    if ( _ehmod==0 ) {// THERMAL KAPPA AND OMEGA TURBULENCE MODEL

        const double TkCorr = ( ( 1.-_khlog ) + _khlog/kappaT );
        const double TwCorr = ( ( 1.-_whlog ) + _whlog/omegaT );

        // THERMAL SOURCE
        source[0] = prod_kt*TkCorr;
        source[1] = ( __CP1-1. ) *prod_kt*TwCorr*omegaT/kappaT;

        // THERMAL DISSIPATION
        diss[0]   = __CMU*omegaT;
        diss[1]   = ( __CD1-1. ) *__CMU*omegaT;

        // MECHANICAL DISSIPATION AND SOURCE
        meccterm[0]  = f_exp*__CMU*omega*R_corr;                           // DISSIPATION
        meccterm[1]  = prod_k_corr*__CP2*prod_k*omegaT*TwCorr/kappa;       // SOURCE
        
    } else {// THERMAL KAPPA AND EPSILON TURBULENCE MODEL
        double kappaT = TKappaAndOmega[0];
        double epsilonT = TKappaAndOmega[1];

        // THERMAL SOURCE
        source[0] = prod_kt;
        source[1] = __CP1 *prod_kt*epsilonT/kappaT;

        // THERMAL DISSIPATION
        diss[0]   = epsilonT/kappaT;
        diss[1]   = __CD1 *__CMU*epsilonT/kappaT;

        // MECHANICAL DISSIPATION AND SOURCE
        meccterm[0]  = f_exp*__CMU*omega/kappaT;                      // DISSIPATION
        meccterm[1]  = prod_k_corr*__CP2*prod_k*epsilonT/kappa;       // SOURCE
    }
    return;
}

void TurbUtils::DynTurInitValues ( double & kappa, double & omega, double WallDist, bool FlatProfile ) {

    if ( FlatProfile ) {
        const double Re_h     = _vmid*_diameter/_nu;
        const double In       = 0.16*pow ( Re_h,-0.125 );
        const double len      = 0.07*_diameter;
        const double k_in     = 1.5* ( In*_vmid ) * ( In*_vmid );
        const double e_in     = __CMU*k_in*sqrt ( k_in ) /len;
        const double w_in     = sqrt ( k_in ) /len;
        const double n_in    = k_in/w_in;

        if ( _numod==1 ) kappa = n_in;
        else          kappa = ( 1-_klog ) *k_in + _klog*log ( k_in );
        if ( _emod==1 ) omega = e_in;
        else         omega = ( 1-_wlog ) *w_in + _wlog*log ( w_in );
    } else {
        // FIRST MESH POINT ASSUMED TO BE AT y+ = 10.
        const double utau = ( _InputUtau<0 ) ? 0.5*_nu/_BoundWallDist : _InputUtau;
        DynTurInitValues ( kappa, omega, WallDist, utau );
    }

    return;
}

void TurbUtils::DynTurInitValues ( double & kappa, double & omega, double WallDist, double Utau ) {
    double yp     = WallDist*Utau/_nu;
//     if(yp<1.)     yp = 1.;
    double WD     = yp*_nu/Utau;

    double nut    = _nu/ ( 1./ ( 0.001093*yp*yp*yp ) + 1./ ( 0.41*yp ) ); // _nu/(1./(0.001093*WD*WD*WD) + 1./(0.41*WD));
    double w_lin  = 2.*_nu/ ( 0.09*WD*WD );
    double w_log  = Utau/ ( sqrt ( __CMU ) *0.41*WD );
//     double w      = (yp < 5.)? w_lin:w_log;
    double w      = sqrt ( w_lin*w_log );
    w = w_log;
    double k      = w * nut;

    double klin = Utau*Utau* ( 0.041*yp*yp );
    double klog = Utau*Utau/0.3;
    double k3   = Utau*Utau*pow ( yp,-0.2 ) *9;
    k = 1./ ( 1./klin + 1./k3 );

    double k_log2 = 1/ ( 5* ( yp+1.e-10 ) );
    k = 1./ ( 1./k + 1./k_log2 );
//     k = Utau*Utau/sqrt(__CMU);
//     if(k < exp(-11.)) k = exp(-11.);
//     w = w_log/3.;

    if ( _numod==1 ) kappa = nut;
    else  kappa = ( 1-_klog ) *k + _klog*log ( k );

    if ( _emod==1 ) omega = Utau*Utau*Utau/ ( 0.41*WD );
    else omega = ( 1-_wlog ) *w + _wlog*log ( w );

    return;
}

void TurbUtils::TherTurInitValues ( double& kappaT, double& omegaT, double WallDist, double NormFlux, double AvVel, double Diameter, bool FlatProfile ) {

    double kappa, omega, yp;
    DynTurInitValues ( kappa, omega, WallDist, FlatProfile );

    double real_k, real_w;
    real_k = ( 1-_klog ) *kappa + _klog*exp ( _klog*kappa );
    real_w = ( 1-_wlog ) *omega + _wlog*exp ( _wlog*omega );
    if ( FlatProfile ) {
        kappaT = real_k* ( _alpha/_nu );
        omegaT = real_w* ( _alpha/_nu );
        kappaT = real_k;
        omegaT = real_w;

    } else {
        // FIRST MESH POINT ASSUMED TO BE AT y+ = 10.
        const double utau  = 10.*_nu/_BoundWallDist;
        const double T_tau = NormFlux/utau;            // NormFlux needs to be _qs/_rhof*_cp0
        yp     = utau*WallDist/_nu;

        double wh1, wh2;
        wh2    = real_w*sqrt ( _alpha/_nu );
        wh1    = 2.*_alpha/ ( __CMU*WallDist*WallDist );
        omegaT     = sqrt ( wh1*wh2 );
        omegaT     = real_w*4;                                   // Assumed value of R = omega/omegaT = 0.25
        kappaT     = real_k*T_tau*T_tau/ ( utau*utau*_alpha/_nu );
        const double yPr = yp*_nu/_alpha;
        const double khlin = 0.0001*pow ( yPr,1.5 );
        const double khlog = 1./ ( yPr*yPr*yPr );
        kappaT     = khlin*khlog/ ( khlin+khlog );
    }

    if ( _ehmod==1 ) omegaT = kappaT * omegaT * __CMU;          // epsilon
    else omegaT = ( 1-_whlog ) * omegaT + _whlog * log ( omegaT ); // omega or log_omega
    kappaT = ( 1-_khlog ) * kappaT + _khlog * log ( kappaT );   // kappa or log_kappa
    return;
}




void TurbUtils::PrintStatus ( ) {

    ofstream output;
    output.open ( "SimulationSpecs.txt" );
    output << " MESSAGE PRINTED BY TURBUTILS CLASS "<<endl;
    output << endl << " TurbUtils filled: "<<_IsFilled<<endl;
    output << "---------------------------------------" <<endl;
    output << " Solved Equations: "<<endl;
    output << "\t Navier Stokes:        "<<_SolveNS<<endl;
    output << "\t Temperature:          "<<_SolveT<<endl;
    output << "\t Dynamical Turbulence: "<<_SolveTBK<<endl;
    output << "\t Thermal Turbulence:   "<<_SolveTTBK<<endl;
    output << "---------------------------------------" <<endl;
    output << " Realizability and Constrains:"<<endl;
    output << "\t Durbin Realizability Constrain: "<<_Durbin<<endl;
    output << "\t Yap Correction:                 "<<_YapCorr<<endl;
    output << "---------------------------------------" <<endl;
    output << " Geometry Parameters: "<<endl;
    output << "\t Fixed Wall Distance: "<<_BoundWallDist <<endl;
    output << "---------------------------------------" <<endl;
    output << " Physical Properties: "<<endl;
    output << "\t Kinematic Viscosity: "<<_nu<<endl;
    output << "\t Thermal diffusivity: "<<_alpha<<endl;

    return;
}

void TurbUtils::read_file() { // READING Tparameter.in ======================================
    std::string string_value;  // read double, string, dummay
    std::ostringstream file, file1, file2;
    std::vector<std::string> FILES;

    if ( __NMeshes>0 ) {
        file   <<getenv ( "APP_PATH" ) <<"/DATA/DATA"<<std::to_string ( __MeshID ) <<"/Turbulence.in";
        file1  <<getenv ( "APP_PATH" ) <<"/DATA/DATA"<<std::to_string ( __MeshID ) <<"/GeometrySettings.in";
        file2  <<getenv ( "APP_PATH" ) <<"/DATA/DATA"<<std::to_string ( __MeshID ) <<"/MaterialProperties.in";
    } else {              
        file   <<getenv ( "APP_PATH" ) <<"/DATA/Turbulence.in";
        file1  <<getenv ( "APP_PATH" ) <<"/DATA/GeometrySettings.in";
        file2  <<getenv ( "APP_PATH" ) <<"/DATA/MaterialProperties.in";
    }

    FILES.push_back ( file1.str() );
    FILES.push_back ( file.str() );
    FILES.push_back ( file2.str() );

    for ( int i=0; i<FILES.size(); i++ ) {

        std::ifstream fin;
        fin.open ( FILES[i].c_str() ); // stream file
        std::string buf="";
        if ( fin.is_open() ) {
            std::cout << "\nInit Reading = " << FILES[i] <<  std::endl;
        }
        if ( fin.is_open() ) {
            while ( buf != "/" ) {
                fin >> buf;    // find "/" file start
            }
            fin >> buf;
            while ( buf != "/" ) {
                if ( buf == "#" ) {
                    getline ( fin, buf ); // comment line
                } else {
                    fin >> string_value;
                    _FileMap[buf] = string_value;
                    std::cout<<buf<<"\t"<<string_value<<std::endl;
                }
                fin >> buf;
            }
        } else {
            std::cerr << "TurbUtils.read_file(): no parameter file found\t ->" <<FILES[i]<< std::endl;
            abort();
        }
        std::cout << "TurbUtils.read_file() End Reading file " <<  FILES[i]<< std::endl;
        fin.close();
    }

    print_par();
    return;
}// END READING ==========================================================================


void TurbUtils::print_par() {

    std::cout << "\033[038;5;196;1m \n----------------------------------------------\n  TURBUTILS MAP \n \033[0m"  << std::endl;
    for ( std::map<std::string, std::string >::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it ) {
        std::cout << it->first << " " << it->second << "\n";
    }
    std::cout << "\033[038;5;196;1m \n----------------------------------------------\n \033[0m"  << std::endl;

    return;
}


void TurbUtils::CalcWallFuncKappaAndOmega ( double KappaAndOmega[], int NodeOnBound, double WallDist, double utau ) {
    double k,w;
    double wd = WallDist + 1.e-10;
    double yplus = wd*utau/_nu;
    const double wlin  = 2.*_nu/ ( 0.09*wd*wd );
    const double wlog  = utau/ ( 0.41 * wd * sqrt ( 0.09 ) );

    const double w_der_lin = -4.*_nu/ ( 0.09*wd*wd*wd );
    const double w_der_log = -1.*utau/ ( sqrt ( 0.09 ) *0.41*wd*wd );

    w = ( yplus<10. ) ?  wlin: wlog;
    const double WDerMidPoint = ( yplus<10. ) ? w_der_lin:w_der_log;
    if ( fabs ( utau ) <1.e-9 )   w = wlin;
    if ( NodeOnBound == 0 )   w = w - WDerMidPoint*wd;
    double Mulin = _nu* ( 0.001093*yplus*yplus*yplus );
    double Mulog = _nu* ( 0.41*yplus );
    double Mu = 1/ ( 1/ ( Mulin ) + 1/ ( Mulog ) );
    double klin = utau*utau* ( 0.041*yplus*yplus );
    double klog = utau*utau/0.3;
    k = 1/ ( 1/klin + 1/klog );
    if ( _Nagano==1 ) {
        if ( _numod==1 )   KappaAndOmega[0] = Mu*_nu;
        else      KappaAndOmega[0] = ( 1.-_klog ) *k + _klog*log ( k );
        if ( _emod==1 ) KappaAndOmega[1] = __CMU*k*w;
        else      KappaAndOmega[1] = ( 1.-_wlog ) *w + _wlog*log ( w );
    } else if ( _Wilcox==1 ) {
        KappaAndOmega[0] = ( 1.-_klog ) *k + _klog*log ( k );
        KappaAndOmega[1] = ( 1.-_wlog ) *w + _wlog*log ( w );
    }
    return;
}
void TurbUtils::CalcWallFuncThermalKappaAndOmega ( double KappaAndOmega[], int NodeOnBound, double WallDist, double utau ) {
    double k,w;
    double yplus = WallDist*utau/_nu;
    const double wlin  = 2.*_nu/ ( 0.09*WallDist*WallDist );
    const double wlog  = utau/ ( 0.41 * WallDist * sqrt ( 0.09 ) );

    const double w_der_lin = -4.*_nu/ ( 0.09*WallDist*WallDist*WallDist );
    const double w_der_log = -1.*utau/ ( sqrt ( 0.09 ) *0.41*WallDist*WallDist );

    w = ( yplus<10. ) ?  wlin: wlog;
    const double WDerMidPoint = ( yplus<10. ) ? w_der_lin:w_der_log;
    if ( fabs ( utau ) <1.e-9 )   w = wlin;
    if ( NodeOnBound == 0 )   w = w - WDerMidPoint*WallDist;
    double Mulin = _nu* ( 0.001093*yplus*yplus*yplus );
    double Mulog = _nu* ( 0.41*yplus );
    double Mu = 1/ ( 1/ ( Mulin ) + 1/ ( Mulog ) );
    double klin = utau*utau* ( 0.041*yplus*yplus );
    double klog = utau*utau/0.3;
    k = 1/ ( 1/klin + 1/klog );

    w *= _alpha/_nu;
    k *= _alpha/_nu;
    if ( _Nagano==1 ) {
        if ( _numod==1 )   KappaAndOmega[0] = Mu*_nu;
        else      KappaAndOmega[0] = ( 1.-_klog ) *k + _klog*log ( k );
        if ( _emod==1 ) KappaAndOmega[1] = __CMU*k*w;
        else      KappaAndOmega[1] = ( 1.-_wlog ) *w + _wlog*log ( w );
    } else if ( _Wilcox==1 ) {
        KappaAndOmega[0] = ( 1.-_klog ) *k + _klog*log ( k );
        KappaAndOmega[1] = ( 1.-_wlog ) *w + _wlog*log ( w );
    }
    return;
}


void TurbUtils::DynTurNearWallValues ( double & kappa, double & omega, double WallDist, double Utau ) {

    const double y_plus = WallDist * Utau / _nu;

    // omega derivative on cell mid point
    const double w_der_lin = -4.*_nu/ ( 0.09*WallDist*WallDist*WallDist );
    const double w_der_log = -1.*Utau/ ( sqrt ( 0.09 ) *0.41*WallDist*WallDist );
    double OmegaNearWallDer = ( y_plus<5. ) ? w_der_lin:w_der_log;
    // omega modelled value on cell mid point
    const double wlin  = 2.*_nu/ ( 0.09*WallDist*WallDist );
    const double wlog  = Utau/ ( 0.41 * WallDist * sqrt ( 0.09 ) );
    double OmegaNearWall = ( y_plus<5 ) ? wlin:wlog;
    // omega linear reconstruction up to wall
    double OmegaWall = OmegaNearWall - WallDist * OmegaNearWallDer;

    // eddy viscosity on cell mid point
    double a = 0.41;
    double c = 0.001093;
    double den = ( a + c*y_plus*y_plus );

    double Mu = _nu/ ( 1./ ( c*y_plus*y_plus*y_plus ) + 1./ ( 0.41*y_plus ) );
    double MuDer = _nu*c*a * ( 3.*a*y_plus*y_plus + c*y_plus*y_plus*y_plus*y_plus ) / ( den*den );
    double Mu2Der = _nu*c*a *2*y_plus*c* ( a*y_plus*y_plus - 3.*c ) / ( den*den*den );

    // eddy viscosity reconstructed up to wall
    double MuWall = Mu - MuDer * y_plus + Mu2Der*y_plus*y_plus;
    MuWall = ( MuWall > 0 ) ? MuWall:1.e-20;

    // reconstruction of kappa wall value
    double KappaWall = MuWall * OmegaWall;

    if ( _numod==1 ) kappa = MuWall;
    else  kappa = ( 1-_klog ) *KappaWall + _klog*log ( KappaWall );

    if ( _emod==1 ) omega = __CMU * KappaWall * OmegaWall;
    else omega = ( 1-_wlog ) *OmegaNearWall + _wlog*log ( OmegaNearWall );
    return;
}

void TurbUtils::ThermTurNearWallValues ( double & kappa, double & omega, double WallDist, double Utau ) {

    const double y_plus = WallDist * Utau / _nu;

    // omega derivative on cell mid point
    const double w_der_lin = -4.*_alpha/ ( 0.09*WallDist*WallDist*WallDist );
    const double w_der_log = -1.*_alpha*Utau/ ( _nu*sqrt ( 0.09 ) *0.41*WallDist*WallDist );
    double OmegaNearWallDer = ( y_plus<5. ) ? w_der_lin:w_der_log;
    // omega modelled value on cell mid point
    const double wlin  = 2.*_nu/ ( 0.09*WallDist*WallDist );
    const double wlog  = Utau/ ( 0.41 * WallDist * sqrt ( 0.09 ) );
    double OmegaNearWall = ( y_plus<5 ) ? wlin:wlog;
    // omega linear reconstruction up to wall
    double OmegaWall = OmegaNearWall - WallDist * OmegaNearWallDer;

    double KappaWall = 1.e-8;

    kappa = ( 1-_khlog ) *KappaWall + _khlog*log ( KappaWall );

    if ( _ehmod==1 ) omega = __CMU * KappaWall * OmegaWall;
    else omega = ( 1-_whlog ) *OmegaNearWall + _whlog*log ( OmegaNearWall );
    return;
}
