#include "EquationsMap.h"
#include "EquationSystemsExtendedM.h"
#include "Equations_conf.h"
#include "MGSolverBase.h"
#include <vector>
#include <fstream>


#ifdef  NS_EQUATIONS
#include "MGSolverNS_proj.h"  // NS ================
#include "MGSolverNS_coup.h"  // NS ================
#include "MGSolverP.h"
// #if NS_EQUATIONS ==3
// #include "MGSolverNSunique.h"
// #endif
#endif // -----------------------------------------
#ifdef  NSA_EQUATIONS
#include "MGSolverNSA.h"  // NS ================
#endif // -----------------------------------------
#ifdef  TBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTBK.h"
#endif // -----------------------------------------
#ifdef  RANS_EQUATIONS  // Turbulence -------------
#include "Nagano_Log.h"
#include "Nagano_KW.h"
#include "Nagano_KE.h"
#include "Wilcox.h"
#include "WilcoxLog.h"
#endif // -----------------------------------------
#ifdef _TURBULENCE_
#include "MGSolverTURB.h"
#endif
#ifdef  TBKA_EQUATIONS  // Turbulence -------------
#include "MGSolverTBKA.h"
#endif // -----------------------------------------
#ifdef  T_EQUATIONS // Temperature ================
#include "MGSolverT.h"
#endif // -----------------------------------------
#ifdef  CTRL_EQUATIONS // Temperature ================
#include "MGSolverCTRL.h"
#endif // -----------------------------------------
#ifdef  TA_EQUATIONS // Temperature adjoint ================
#include "MGSolverTA.h"
#endif // -----------------------------------------
#ifdef  T_G_EQUATIONS // Temperature control ================
#include "MGSolverT_G.h"
#endif // -----------------------------------------
#ifdef  T_COUP_EQUATIONS // Temperature control coupled opt system ================
#include "MGSolverT_COUP.h"
#endif // -----------------------------------------
#ifdef  ALFA_EQUATIONS
#include "MGSolverALFA.h"  // two fluids ===========
#endif // -----------------------------------------
#ifdef TTBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTTBK.h"
#ifdef _TURBULENCE_
#include "MGSolverTURB.h"
#endif
#endif
#ifdef RANS_THERMAL_EQUATIONS  // Turbulence -------------
#include "Nagano_LogT.h"
#include "Nagano_KWT.h"
#include "Nagano_KET.h"
#endif
#ifdef SM_EQUATIONS  // Structural Mechanics ======
#include "MGSolverSM.h"
#include "MGSolverDS.h"
#endif
#ifdef FSI_EQUATIONS  // FSI =======================
#include "MGSolverFSI.h"
#include "MGSolverDS.h"
#ifdef FSIC_EQUATIONS
#include "MGSolverFSIC.h"
#endif
#endif
#ifdef FSIA_EQUATIONS  // FSI adjoint =======================
#include "MGSolverFSIA.h"
#endif
#ifdef COLOR_EQUATIONS
#include "MGSolverCOL.h"
#endif
#ifdef IMMERSED_BOUNDARY
#include "MGSolverIB.h"
#endif
#ifdef DS_EQUATIONS
#include "MGSolverDS.h"
#endif
#include "MGSolverDA.h"
#ifdef   TWO_PHASE
#include "MGSolverCC.h"
#endif



/// EquationsMap standard constructor
EquationsMap::EquationsMap() {
    Fill_FIELD_map();
    ReadEquationsToAdd();
    Fill_pbName();
    }


/// Function that handles #EquationsMap::_myproblemP vector filling
void EquationsMap:: Fill_pbName ( ) {
    _myproblemP.clear();

#ifdef _TURBULENCE_
    ReadTurbulenceInfo();
#endif


    int i = 0;
    std::cout << "------------------------------------------------------------------------------ \n";
    std::cout << " EquationMap::Fill_pbName() \n";

    for ( std::map<std::string, int>::iterator it = _EquationsToAdd.begin(); it != _EquationsToAdd.end(); ++it ) {
        if ( ( it->first ).compare ( 0, 3, "MG_" ) == 0 ) {
            if ( it->second != 0 ) {
                FIELDS ff =  _map_str2field[it->first];
                _myproblemP.push_back ( ff );
                std::cout << " Adding field " << it->first << " with Equation tab number " << _myproblemP[i] << "\n";
                i++;
#ifdef _TURBULENCE_
                if ( it->first == "MG_DynamicalTurbulence" ) {
                    _myproblemP.push_back ( DIST );
                    std::cout << " Adding wall distance field with Equation tab number " << _myproblemP[i] << "\n";
                    i++;
                    _myproblemP.push_back ( MU_T );
                    std::cout << " Adding dynamical turbulence field with Equation tab number " << _myproblemP[i] << "\n";
                    i++;
                    }

                if ( it->first == "MG_ThermalTurbulence" ) {
                    _myproblemP.push_back ( ALPHA_T );
                    std::cout << " Adding thermal turbulence with Equation tab number " << _myproblemP[i] << "\n";
                    i++;
                    }
#endif
                }
            }
        }

    std::cout << " ------------------------------------------------------------------------------ \n";

    return;
    }

// =============================================================================
/// This function adds equations to #FEMUS::_mg_equations_map
void EquationsMap::FillEquationMap (
    EquationSystemsExtendedM & EqMap
) { // ==========================================================================

    int n_equations = _myproblemP.size();

    // Initialization of problem classes

    for ( int iname = 0; iname < n_equations; iname++ ) {
        FIELDS problem = _myproblemP[iname];

        if ( _myproblemP[iname] == NS_F  || _myproblemP[iname] == NSX_F || _myproblemP[iname] == NSY_F || _myproblemP[iname] == NSZ_F ) {
            initNavierStokes ( EqMap );
            }

        if ( _myproblemP[iname] == FS_F ) {
            initFluidStructure ( EqMap );
            }

        if ( _myproblemP[iname] == SDS_F ) {
            initDisplacements ( EqMap );
            }

        if ( _myproblemP[iname] == T_F ) {
            initTemperature ( EqMap );
            }

        if ( _myproblemP[iname] == K_F   ||  _myproblemP[iname] == EW_F ) {
            initDynamicTurbulence ( EqMap );
            }

        if ( _myproblemP[iname] == KTT_F ||  _myproblemP[iname] == EWTT_F ) {
            initThermalTurbulence ( EqMap );
            }

        if ( _myproblemP[iname] == CO_F ) {
            initColor ( EqMap );
            }

        if ( _myproblemP[iname] == NSA_F  || _myproblemP[iname] == NSAX_F || _myproblemP[iname] == NSAY_F || _myproblemP[iname] == NSAZ_F ) {
            initAdjointNavierStokes ( EqMap );
            }

        if ( _myproblemP[iname] == FSA_F ) {
            initAdjointFSI ( EqMap );
            }

        if ( _myproblemP[iname] == TA_F ) {
            initAdjointTemperature ( EqMap );
            }

        if ( _myproblemP[iname] == KA_F   ||  _myproblemP[iname] == EWA_F ) {
            initAdjointDynamicTurbulence ( EqMap );
            }

        if ( _myproblemP[iname] == CTRL_F ) {
            initControlTemperature ( EqMap );
            }

        if ( _myproblemP[iname] == IB_F ) {
            initImmersedBoundary ( EqMap );
            }
          
      if ( _myproblemP[iname]== TCOUP_F ) {
          initCoupledTemperature ( EqMap );
          } 

        }
        

    // This Function calls the MGSolDA::init_ext_fields()
    for ( auto eqn = EqMap._equations.begin(); eqn != EqMap._equations.end(); eqn++ ) {
        MGSolBase * mgsol = eqn->second; // get the pointer
        mgsol -> setUpExtFieldData();
        setProblems ( mgsol );
        }

    return;
    }

// ================================== NS_EQUATIONS ================================
void EquationsMap::initNavierStokes ( EquationSystemsExtendedM & EqMap ) {

#ifdef NS_EQUATIONS
    const int NS_sol_type = _EquationsToAdd["MG_NavierStokes"] ;

    if ( NDOF_K > 0 ) {
        _nvars[0] = 1;    // konstant(0) approx
        _nvars[1] = 0;
        }

    if ( NS_sol_type == 2 ) { // SEGREGATED SOLVER
        _nvars[2] = 1; // quadratic(2) approx
        _nvars[0] = 0;
        _nvars[1] = 0; // linear(1) approx
        EqMap.AddSolver<MGSolNS_proj> ( "NS0X", NS_F,     _nvars[0], _nvars[1], _nvars[2], "u" );
        EqMap.AddSolver<MGSolNS_proj> ( "NS0Y", NS_F + 1, _nvars[0], _nvars[1], _nvars[2], "v" );
#if DIMENSION==3
        EqMap.AddSolver<MGSolNS_proj> ( "NS0Z", NS_F + 2, _nvars[0], _nvars[1], _nvars[2], "w" );
#endif
        }

    if ( NS_sol_type == 1 ) { // MONOLITHIC SOLVER
        _nvars[2] = DIMENSION; // quadratic(2) approx
        _nvars[0] = 0;
        _nvars[1] = 1; // linear(1) approx
        EqMap.AddSolver<MGSolNS_coup> ( "NS0", NS_F, _nvars[0], _nvars[1], _nvars[2], "u" );
        }

    if ( NS_sol_type == 2 ) { // PRESSURE FOR SEGREGATED SOLVER
        _nvars[0] = 0;
        _nvars[1] = 1;
        _nvars[2] = 0; // only  Linear(1) approx
        EqMap.AddSolver<MGSolP> ( "NS2P", NS_F, _nvars[0], _nvars[1], _nvars[2], "p" );
        }

// #if NS_EQUATIONS==3
//   _nvars[0] = _nvars[1] = 0;
//   _nvars[2] = 1;
//   EqMap.AddSolver<MGSolNS_1comp> ( "NS0", NS_F, _nvars[0], _nvars[1], _nvars[2], "u" );
// #endif

#endif
    return;
    }
// ================================== FSI_EQUATIONS ================================
void EquationsMap::initFluidStructure ( EquationSystemsExtendedM & EqMap ) {
#ifdef FSI_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = FSI_EQUATIONS % 2;
    _nvars[2] = ( ( FSI_EQUATIONS == 2 ) ? 1 : DIMENSION );

    if ( NDOF_K > 0 ) {
        _nvars[0] = 1;
        _nvars[1] = 0;
        }

#if FSI_EQUATIONS!=2
    EqMap.AddSolver<MGSolFSI> ( "FSI0", FS_F, _nvars[0], _nvars[1], _nvars[2], "u" );
#else
    EqMap.AddSolver<MGSolFSI> ( "FSI0X", FS_F, _nvars[0], _nvars[1], _nvars[2], "u" );
    EqMap.AddSolver<MGSolFSI> ( "FSI0Y", FS_F + 1, _nvars[0], _nvars[1], _nvars[2], "v" );
#if DIMENSION==3
    EqMap.AddSolver<MGSolFSI> ( "FSI0Z", FS_F + 2, _nvars[0], _nvars[1], _nvars[2], "z" );
#endif
#endif
#if (FSI_EQUATIONS%2==0)
    _nvars[0] = 0;
    _nvars[1] = 1;
    _nvars[2] = 0;
    EqMap.AddSolver<MGSolFSIP> ( "FSIP", FS_F + 3, _nvars[0], _nvars[1], _nvars[2], "p" );
#endif

// DISPLACEMENTS
    initDisplacements ( EqMap );
#endif

    return;
    }
// ================================== DS_EQUATIONS ================================
void EquationsMap::initDisplacements ( EquationSystemsExtendedM & EqMap ) {
#ifdef DS_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = 1;
    EqMap.AddSolver<MGSolDS> ( "SDSX", SDSX_F, _nvars[0], _nvars[1], _nvars[2], "dx" );
    EqMap.AddSolver<MGSolDS> ( "SDSY", SDSY_F, _nvars[0], _nvars[1], _nvars[2], "dy" );
#if (DIMENSION==3)
    EqMap.AddSolver<MGSolDS> ( "SDSZ", SDSZ_F, _nvars[0], _nvars[1], _nvars[2], "dz" );
#endif
#endif

    return;
    }

// ================================== T_EQUATIONS ================================
void EquationsMap::initTemperature ( EquationSystemsExtendedM & EqMap ) {
#ifdef T_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = 1; // only Quadratic[2] approx
    EqMap.AddSolver<MGSolT> ( "T", T_F, _nvars[0], _nvars[1], _nvars[2], "T" );
#endif
    return;
    }

// ================================== TBK_EQUATIONS ================================
void EquationsMap::initDynamicTurbulence ( EquationSystemsExtendedM & EqMap ) {
#ifdef TBK_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = ( TBK_EQUATIONS % 2 ) + 1;

#if ((TBK_EQUATIONS/2)==1)  // k-epsilon      
    EqMap.AddSolver<MGSolTBK> ( "K", K_F, _nvars[0], _nvars[1], _nvars[2], "kt" );
    EqMap.AddSolver<MGSolTBK> ( "K2", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "et" );
#endif

#if ((TBK_EQUATIONS/2)==2)  // k-omega  
    EqMap.AddSolver<MGSolTBK> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "kt" );
    EqMap.AddSolver<MGSolTBK> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "wt" );
#endif
#endif

#ifdef RANS_EQUATIONS
    _nvars[0] = 0;                   // Costant(0)
    _nvars[1] = 0;                   // Linear(1)
    _nvars[2] = ( RANS_EQUATIONS % 2 ) + 1; // Quadratic(2)


    const std::string TurbModel = ( _TurbulenceModel["RANS_dynamic"] != "" ) ? _TurbulenceModel["RANS_dynamic"] : "default";
    int Model = ::DynTurbModelMap.at ( TurbModel );

    switch ( Model ) {
        case nagano_log: {
            EqMap.AddSolver<MGSolNaganoLog> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "lk" );
            EqMap.AddSolver<MGSolNaganoLog> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "lw" );
            break;
            }

        case nagano_kw: {
            EqMap.AddSolver<MGSolNaganoKW> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "kappa" );
            EqMap.AddSolver<MGSolNaganoKW> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "omega" );
            break;
            }

        case ::nagano_ke: {
            EqMap.AddSolver<MGSolNaganoKE> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "kappa" );
            EqMap.AddSolver<MGSolNaganoKE> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "epsilon" );
            break;
            }

        case wilcox: {
            EqMap.AddSolver<MGSolWilcox> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "kappa" );
            EqMap.AddSolver<MGSolWilcox> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "omega" );
            break;
            }

        case ::wilcox_log: {
            EqMap.AddSolver<MGSolWilcoxLog> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "log_k" );
            EqMap.AddSolver<MGSolWilcoxLog> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "log_w" );
            break;
            }

        case default_dyn: {
            std::cerr << "==================================================================\n";
            std::cerr << " MGEquationSystem: settin default dyn turb model: nagano_log \n";
            std::cerr << "==================================================================\n";
            EqMap.AddSolver<MGSolNaganoLog> ( "K2K", K_F, _nvars[0], _nvars[1], _nvars[2], "lk" );
            EqMap.AddSolver<MGSolNaganoLog> ( "K1W", K_F + 1, _nvars[0], _nvars[1], _nvars[2], "lw" );
            break;
            }

        default: {
            std::cerr << "==================================================================\n";
            std::cerr << " MGEquationSystem: unknown turbulence model " << TurbModel << "\n";
            std::cerr << " change value inside Turbulence.in -> RANS_dynamic value    \n";
            std::cerr << "==================================================================\n";
            abort();
            break;

            }
        }

#endif    // end TBK model =============================================================    

#ifdef _TURBULENCE_
    _nvars[2] = 1;
    EqMap.AddSolver<MGSolTURB> ( "DIST", DIST, _nvars[0], _nvars[1], _nvars[2], "dist" );
    EqMap.AddSolver<MGSolTURB> ( "MU_T", MU_T, _nvars[0], _nvars[1], _nvars[2], "muT" );
#endif


    return;
    }

// ================================== TTBK_EQUATIONS ================================
void EquationsMap::initThermalTurbulence ( EquationSystemsExtendedM & EqMap ) {
#ifdef TTBK_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = ( TTBK_EQUATIONS % 2 ) + 1;

#if ((TTBK_EQUATIONS/2)==1)  // k-epsilon      
    EqMap.AddSolver<MGSolTTBK> ( "TK", KTT_F, _nvars[0], _nvars[1], _nvars[2], "kh" );
    EqMap.AddSolver<MGSolTTBK> ( "TK2", KTT_F + 1, _nvars[0], _nvars[1], _nvars[2], "eh" );
#endif

#if ((TTBK_EQUATIONS/2)==2)  // k-omega  
    EqMap.AddSolver<MGSolTTBK> ( "TK", KTT_F, _nvars[0], _nvars[1], _nvars[2], "kh" );
    EqMap.AddSolver<MGSolTTBK> ( "TK2", KTT_F + 1, _nvars[0], _nvars[1], _nvars[2], "wh" );
#endif

#ifdef _TURBULENCE_
    _nvars[2] = 1;
    EqMap.AddSolver<MGSolTURB> ( "ALPHA_T", ALPHA_T, _nvars[0], _nvars[1], _nvars[2], "alphaT" );
#endif

#endif
    return;
    }

// ================================== COLOR_EQUATIONS ================================
void EquationsMap::initColor ( EquationSystemsExtendedM & EqMap ) {
// COLOR EQUATION
#ifdef COLOR_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = 1; // only Quadratic[2] approx
    EqMap.AddSolver<MGSolCOL> ( "C", CO_F, _nvars[0], _nvars[1], _nvars[2], "c" );
    EqMap.AddSolver<MGSolCOL> ( "CK", CO_F + 1, _nvars[0], _nvars[1], _nvars[2], "k" );
#endif
    return;
    }


void EquationsMap::initAdjointFSI ( EquationSystemsExtendedM & EqMap ) {
#ifdef FSIA_EQUATIONS //    FLUID-STRUCTURE adjoint---  
// ====================================================================

    _nvars[0] = 0;
    _nvars[1] = FSIA_EQUATIONS % 2; // Costant(1)  Linear(0)
    _nvars[2] = ( ( FSIA_EQUATIONS == 2 ) ? 1 : DIMENSION ); // quadratic(2) Approximation

    if ( NDOF_K > 0 ) {
        _nvars[0] = 1;    // konstant(0) approx
        _nvars[1] = 0;
        }

#if FSIA_EQUATIONS!=2
    EqMap.AddSolver<MGSolFSIA> ( "FSIA0", FSA_F, _nvars[0], _nvars[1], _nvars[2], "ua" );
#else
    EqMap.AddSolver<MGSolFSIA> ( "FSIA0X", FSA_F, _nvars[0], _nvars[1], _nvars[2], "ua" );
    EqMap.AddSolver<MGSolFSIA> ( "FSIA0Y", FSA_F + 1, _nvars[0], _nvars[1], _nvars[2], "va" );
#if DIMENSION==3
    EqMap.AddSolver<MGSolFSIA> ( "FSIA0Z", FSA_F + 2, _nvars[0], _nvars[1], _nvars[2], "wa" );
#endif
#endif
#if (FSIA_EQUATIONS%2==0) // - FSI_EQUATIONS==0,2  projection -
    _nvars[0] = 0;
    _nvars[1] = 1;
    _nvars[2] = 0; // only  Linear(1) approx
    EqMap.AddSolver<MGSolFSIAP> ( "FSIAP", FSA_F + 3, _nvars[0], _nvars[1], _nvars[2], "pa" );
#endif

// #ifdef DSA_EQUATIONS // -------------  Displacement adjoint--------------------
//   nvars_in[0]=0; nvars_in[1]=0;   // Costant(1)  Linear(0)
//   nvars_in[2]=1;                  // quadratic(2) Approximation
//
//   MGSolDSA* mgsdsdxa=new MGSolDSA(*this,nvars_in,"SDSAX","dxa");  init_eqs(mgsdsdxa);
//   MGSolDSA* mgsdsdya=new MGSolDSA(*this,nvars_in,"SDSAY","dya");  init_eqs(mgsdsdya);
// #if (DIMENSION==3)
//   MGSolDSA* mgsdsdza=new MGSolDSA(*this,nvars_in,"SDSAZ","dza");  init_eqs(mgsdsdza);
// #endif
#endif // ----------------  end  Disp adjoint ---------------------------------  
    return;
    }
// ================================== NSA_EQUATIONS ================================
void EquationsMap::initAdjointNavierStokes ( EquationSystemsExtendedM & EqMap ) {
#ifdef NSA_EQUATIONS
    _nvars[2] = ( ( NSA_EQUATIONS == 2 ) ? 1 : DIMENSION ); // quadratic(2) approx
    _nvars[0] = 0;
    _nvars[1] = NSA_EQUATIONS % 2; // linear(1) approx

    if ( NDOF_K > 0 ) {
        _nvars[0] = 1;    // konstant(0) approx
        _nvars[1] = 0;
        }

#if NSA_EQUATIONS==2     // - NS_EQUATIONS==2 -
    EqMap.AddSolver<MGSolNSA> ( "NSA0X", NSA_F, _nvars[0], _nvars[1], _nvars[2], "ua" );
    EqMap.AddSolver<MGSolNSA> ( "NSA0Y", NSA_F, _nvars[0], _nvars[1], _nvars[2], "va" );
#if DIMENSION==3
    EqMap.AddSolver<MGSolNSA> ( "NSA0Z", NSA_F, _nvars[0], _nvars[1], _nvars[2], "wa" );
#endif
#endif
#if (NSA_EQUATIONS==0 || NSA_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
    EqMap.AddSolver<MGSolNSA> ( "NSA0", NSA_F, _nvars[0], _nvars[1], _nvars[2], "ua" );
#endif

#if (NSA_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
    _nvars[0] = 0;
    _nvars[1] = 1;
    _nvars[2] = 0; // only  Linear(1) approx
    EqMap.AddSolver<MGSolPA> ( "NSA2P", NSA_F, _nvars[0], _nvars[1], _nvars[2], "pa" );
#endif
#endif
    return;
    }
// ================================== TA_EQUATIONS ================================
void EquationsMap::initAdjointTemperature ( EquationSystemsExtendedM & EqMap ) {
#ifdef TA_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = 1; // only Quadratic[2] approx
    EqMap.AddSolver<MGSolTA> ( "TA", TA_F, _nvars[0], _nvars[1], _nvars[2], "TA" );
#endif
    return;
    }
// ================================== TBKA_EQUATIONS ================================
void EquationsMap::initAdjointDynamicTurbulence ( EquationSystemsExtendedM & EqMap ) {
#ifdef TBKA_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = ( TBKA_EQUATIONS % 2 ) + 1;

#if ((TBKA_EQUATIONS/2)==1)  // k-epsilon      
    EqMap.AddSolver<MGSolTBKA> ( "KA", KA_F, _nvars[0], _nvars[1], _nvars[2], "kta" );
    EqMap.AddSolver<MGSolTBKA> ( "K2A", KA_F + 1, _nvars[0], _nvars[1], _nvars[2], "eta" );
#endif

#if ((TBKA_EQUATIONS/2)==2)  // k-omega  
    EqMap.AddSolver<MGSolTBKA> ( "K2KA", KA_F, _nvars[0], _nvars[1], _nvars[2], "kta" );
    EqMap.AddSolver<MGSolTBKA> ( "K1WA", KA_F + 1, _nvars[0], _nvars[1], _nvars[2], "wta" );
#endif
#endif
    return;
    }
// ================================== CTRL_EQUATIONS ================================
void EquationsMap::initControlTemperature ( EquationSystemsExtendedM & EqMap ) {
#ifdef CTRL_EQUATIONS
    _nvars[0] = 0;
    _nvars[1] = 0; // Costant(1)  Linear(0)
    _nvars[2] = 1;              // quadratic(2) Approximation

    EqMap.AddSolver<MGSolCTRL> ( "CTRLX", CTRLX_F, _nvars[0], _nvars[1], _nvars[2], "cx" );

    if ( _EquationsToAdd["MG_ControlTemperature"] == 2 ) {
        EqMap.AddSolver<MGSolCTRL> ( "CTRLY", CTRLX_F + 1, _nvars[0], _nvars[1], _nvars[2], "cy" );

#if (DIMENSION==3)
        EqMap.AddSolver<MGSolCTRL> ( "CTRLZ", CTRLX_F + 2, _nvars[0], _nvars[1], _nvars[2], "cz" );
#endif
        }

#endif
    return;
    }
// ================================== T_COUPLED_EQUATIONS ================================
void EquationsMap::initCoupledTemperature ( EquationSystemsExtendedM & EqMap ) {
#ifdef T_COUP_EQUATIONS
    _nvars[0]=0;
    _nvars[1]=0;          // Costant(1)  Linear(0)
    _nvars[2]=3;          // quadratic(2) Approximation

    EqMap.AddSolver<MGSolTCOUP> ( "TCOUP0", TCOUP_F, _nvars[0],_nvars[1],_nvars[2],"T" );

#endif
    return;
}

void EquationsMap::initImmersedBoundary ( EquationSystemsExtendedM & EqMap ) {
#ifdef IMMERSED_BOUNDARY
    _nvars[0] = 0;
    _nvars[1] = 0;
    _nvars[2] = 1; // only Quadratic[2] approx
    EqMap.AddSolver<MGSolIB> ( "IB1", IB_F, _nvars[0], _nvars[1], _nvars[2], "Col" );
    EqMap.AddSolver<MGSolIB> ( "IB2", IB_F + 1, _nvars[0], _nvars[1], _nvars[2], "VolFrac" );
#endif
    return;
    }



/// This function fills #EquationsMap::_map_str2field map
void EquationsMap::Fill_FIELD_map (
) { // map equation   map_str2field

    _map_str2field["MG_NavierStokes"] = NS_F;
    _map_str2field["NS_F"] = NS_F;   // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["NSX_F"] = NSX_F;   // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["NSY_F"] = NSY_F;  // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["NSZ_F"] = NSZ_F;   // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["MG_FluidStructure"] = FS_F;
    _map_str2field["FS_F"] = FS_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["FSX_F"] = FSX_F;   // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["FSY_F"] = FSY_F;   // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["FSZ_F"] = FSZ_F;    // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["MG_StructuralMechanics"] = SM_F;
    _map_str2field["SM_F"]  = SM_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["SMX_F"] = SMX_F;    // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["SMY_F"] = SMY_F;    // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["SMZ_F"] = SMZ_F;    // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    _map_str2field["MG_Pressure"] = P_F;
    _map_str2field["P_F"]   = P_F;    // [3] -> Pressure (linear (1);NS_EQUATIONS==0 or 2)
    _map_str2field["MG_Temperature"] = T_F;
    _map_str2field["T_F"]   = T_F;    // [4] -> Temperature   (quadratic (2);T_EQUATIONS)
    _map_str2field["MG_DynamicalTurbulence"] = K_F;
    _map_str2field["K_F"]   = K_F;    // [5] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
    _map_str2field["EW_F"]  = EW_F;    // [6] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
    _map_str2field["MG_ThermalTurbulence"] = KTT_F;
    _map_str2field["KTT_F"] = KTT_F;    // [7] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
    _map_str2field["EWTT_F"] = EWTT_F;   // [8] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
    _map_str2field["MG_Displacement"] = SDS_F;
    _map_str2field["SDS_F"] = SDS_F;   // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
    _map_str2field["SDSX_F"] = SDSX_F;   // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
    _map_str2field["SDSY_F"] = SDSY_F;  // [10]-> Displacement (quadratic (2); DS_EQUATIONS)
    _map_str2field["SDSZ_F"] = SDSZ_F;  // [11]-> Displacement (quadratic (2); DS_EQUATIONS)
    _map_str2field["MG_DA"] = DA_F;
    _map_str2field["DA_F"]   = DA_F;  // [12]-> DA solver (quadratic (2); DA_EQUATIONS)
    // adjoint
    _map_str2field["MG_AdjointNavierStokes"] = NSA_F;
    _map_str2field["NSA_F"]  = NSA_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["NSAX_F"]  = NSAX_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["NSAY_F"]  = NSAY_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["NSAZ_F"]  = NSAZ_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["MG_AdjointFluidStructure"] = FSA_F;
    _map_str2field["FSA_F"]  = FSA_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["FSAX_F"]  = FSAX_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["FSAY_F"]  = FSAY_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["FSAZ_F"]  = FSAZ_F;  // [15]-> adjoint NS or FSI equations (Dimension)
    _map_str2field["MG_AdjointDA"] = DA_P;
    _map_str2field["DA_P"]   = DA_P;  // [13]-> DA solver (piecewise; DA_EQUATIONS)
    _map_str2field["MG_AdjointTemperature"] = TA_F;
    _map_str2field["TA_F"]   = TA_F;  // [14]-> Temp adjoint
    _map_str2field["MG_AdjointTurbulence"] = KA_F;
    _map_str2field["KA_F"]   = KA_F;  // [18] + 19 Adjoint turbulence
    _map_str2field["EWA_F"]   = EWA_F;  // [18] + 19 Adjoint turbulence
    //control
    _map_str2field["MG_ColorFunction"] = CO_F;
    _map_str2field["MG_Laplacian"] = CO_F;
    _map_str2field["CO_F"]   = CO_F;   // [16]-> Color function for FSI equations
    _map_str2field["MG_ControlTemperature"] = CTRL_F;
    _map_str2field["CTRL_F"]  = CTRL_F;   // [16]-> Color function for FSI equations
    _map_str2field["CTRLX_F"]  = CTRLX_F;   // [16]-> Color function for FSI equations
    _map_str2field["CTRLY_F"]  = CTRLY_F;   // [16]-> Color function for FSI equations
    _map_str2field["CTRLZ_F"]  = CTRLZ_F;   // [16]-> Color function for FSI
    _map_str2field["MG_ImmersedBoundary"]  = IB_F;    // [16]-> Color function for FSI
    _map_str2field["MG_CoupledTemperature"] = TCOUP_F;  // [29]-> Coupled Temperature optimality system
    _map_str2field["TCOUP_F"]   = TCOUP_F;    // [4] -> Temperature   (quadratic (2);T_EQUATIONS)
    return;
    }

/// Solution sharing between activated equations
void EquationsMap::setProblems ( MGSolBase *& ProbObj ) {
    _PieceEq = _LinearEq = _QuadEq = 0;

    std::cout << "\n------------------------------------------------\n";
    std::cout << "  Activating fields for problem " << ProbObj->_eqname;


    for ( std::map<std::string, int>::iterator eq = _EquationsToAdd.begin(); eq != _EquationsToAdd.end(); ++eq ) {
        std::string EqnName = eq->first;
        int EqnLabel = eq->second;

        if ( EqnLabel > 0 ) {

            // NAVIER STOKES SYSTEM
            if ( EqnName == "MG_NavierStokes" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, NS_F, "NS0", _QuadEq, coupled );

                    if ( coupled == 0 ) {
                        if ( NDOF_K == 1 ) {
                            ProbObj->ActivateScalar ( 0, P_F, "NS2P", _PieceEq );
                            }
                        else {
                            ProbObj->ActivateScalar ( 1, P_F, "NS2P", _LinearEq );
                            }
                        }
                    }

                if ( EqnLabel == 3 ) {
                    ProbObj->ActivateScalar ( 2, NS_F, "NS0", _QuadEq );
                    }
                }

            // ADJOINT NAVIER STOKES
            if ( EqnName == "MG_AdjointNavierStokes" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, NSA_F, "NSA0", _QuadEq, coupled );

                    if ( coupled == 0 ) {
                        ProbObj->ActivateScalar ( 1, P_F, "NSAP", _LinearEq );
                        }
                    }
                }

            // FLUID STRUCTURE
            if ( EqnName == "MG_FluidStructure" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, FS_F, "FSI0", _QuadEq, coupled );

                    if ( coupled == 0 ) {
                        ProbObj->ActivateScalar ( 1, P_F, "FSIP", _LinearEq );
                        }
                    }
                }

            // ADJOINT FLUID STRUCTURE
            if ( EqnName == "MG_AdjointFluidStructure" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, NSA_F, "FSIA0", _QuadEq, coupled );

                    if ( coupled == 0 ) {
                        ProbObj->ActivateScalar ( 1, P_F, "FSIAP", _LinearEq );
                        }
                    }
                }

            // STRUCTURAL MECHANICS
            if ( EqnName == "MG_StructuralMechanics" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, SM_F, "SM0", _QuadEq, coupled );
                    }
                }

            // SOLID DISPLACEMENTS
            if ( EqnName == "MG_Displacement" ) {
                if ( EqnLabel <= 2 ) {
                    int coupled = ( EqnLabel == 1 ) ? 1 : 0;
                    ProbObj->ActivateVectField ( 2, SDS_F, "SDS", _QuadEq, coupled );
                    }
                }

            // TEMPERATURE
            if ( EqnName == "MG_Temperature" ) {
                ProbObj->ActivateScalar ( 2, T_F, "T", _QuadEq );
                }

            // CONTROL TEMPERATURE
            if ( EqnName == "MG_ControlTemperature" ) {
                int vector;
                vector = ( _EquationsToAdd["MG_ControlTemperature"] == 2 ) ? 1 : 0;
                ProbObj->ActivateControl ( 2, CTRL_F, "CTRL", _QuadEq, vector, DIMENSION );
                }

            // ADJOINT TEMPERATURE
            if ( EqnName == "MG_AdjointTemperature" ) {
                ProbObj->ActivateScalar ( 2, TA_F, "TA", _QuadEq );
                }

            // DYNAMICAL TURBULENCE
            if ( EqnName == "MG_DynamicalTurbulence" ) {
                ProbObj->ActivateCoupled ( 2, K_F, "K2K", "K1W", _QuadEq );
#ifdef _TURBULENCE_
                ProbObj->ActivateScalar ( 2, DIST, "DIST", _QuadEq );
                ProbObj->ActivateScalar ( 2, MU_T, "MU_T", _QuadEq );
#endif
                }

            // THERMAL TURBULENCE
            if ( EqnName == "MG_ThermalTurbulence" ) {
                ProbObj->ActivateCoupled ( 2, KTT_F, "TK", "TK2", _QuadEq );
#ifdef _TURBULENCE_
                ProbObj->ActivateScalar ( 2, ALPHA_T, "ALPHA_T", _QuadEq );
#endif
                }

            // COLOR FUNCTION
            if ( EqnName == "MG_ColorFunction" ) {
                ProbObj->ActivateCoupled ( 2, CO_F, "C", "CK", _QuadEq );
                }

            // IMMERSED BOUNDARY
            if ( EqnName == "MG_ImmersedBoundary" ) {
                ProbObj->ActivateCoupled ( 2, IB_F, "IB1", "IB2", _QuadEq );
                }

            // ADJOINT TURBULENCE
            if ( EqnName == "MG_AdjointTurbulence" ) {
                ProbObj->ActivateCoupled ( 2, KA_F, "K2KA", "K1WA", _QuadEq );
                }

            // DA
            if ( EqnName == "MG_DA" ) {
                ProbObj->ActivateScalar ( 2, DA_F, "DA", _QuadEq );
                }
                
            // COUPLED OPTIMAL CONTROL TEMPERATURE
            if ( EqnName == "MG_CoupledTemperature" ) {
//                 if ( EqnLabel <=2 ) {
                    int coupled = 1;
                    ProbObj->ActivateVectField ( 2, TCOUP_F, "TCOUP0", _QuadEq, coupled );
//                     if ( coupled == 0 ) {
//                         ProbObj->ActivateScalar ( 1, P_F, "NSAP", _LinearEq );
//                         }
//                     }
                }
            }

        }

    return;
    }

void EquationsMap::ReadEquationsToAdd() {
    std::ostringstream filename;
    filename << getenv ( "APP_PATH" ) << "/DATA/Equations.in";
    std::ifstream fin;

    fin.open ( filename.str().c_str() ) ; // stream file
    std::string buf = "";

    int value;

    if ( fin.is_open() ) { // -------------------------------------------------------------
        while ( buf != "/" ) {
            fin >> buf;    // find "/" file start
            }

        fin >> buf;

        while ( buf != "/" ) {
            if ( buf == "#" ) {
                getline ( fin, buf );    // comment line
                }
            else {
                fin >> value;
                _EquationsToAdd.insert ( std::pair<std::string, int> ( buf, value ) );
                }

            fin >> buf;
            }
        } // --------------------------------------------------------------------------------------

    return;
    }

void  EquationsMap::ReadTurbulenceInfo() {
    std::ostringstream filename;
    filename << getenv ( "APP_PATH" ) << "/DATA/Turbulence.in";
    std::ifstream fin;

    fin.open ( filename.str().c_str() ) ; // stream file
    std::string buf = "";

    std::string value;

    if ( fin.is_open() ) { // -------------------------------------------------------------
        while ( buf != "/" ) {
            fin >> buf;    // find "/" file start
            }

        fin >> buf;

        while ( buf != "/" ) {
            if ( buf == "#" ) {
                getline ( fin, buf );    // comment line
                }
            else {
                fin >> value;
                _TurbulenceModel.insert ( std::pair<std::string, std::string> ( buf, value ) );
                }

            fin >> buf;
            }
        } // --------------------------------------------------------------------------------------

    return;

    }


