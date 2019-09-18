#ifndef __EquationMapFiller_h__
#define __EquationMapFiller_h__

#include <map>
#include <vector>

class EquationSystemsExtendedM;
class MGSolBase;

enum  FIELDS : int {
  MG_NavierStokes = 0,
  NS_F  = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  NSX_F = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  NSY_F = 1,    // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  NSZ_F = 2,    // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  MG_FluidStructure = 0,
  FS_F  = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  FSX_F = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  FSY_F = 1,    // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  FSZ_F = 2,    // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  MG_StructuralMechanics = 0,
  SM_F  = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  SMX_F = 0,    // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  SMY_F = 1,    // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  SMZ_F = 2,    // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
  MG_Pressure = 3,
  P_F   = 3,    // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
  MG_Temperature = 4,
  T_F   = 4,    // [4] -> Temperature   (quadratic (2),T_EQUATIONS)
  MG_DynamicalTurbulence = 5,
  K_F   = 5,    // [5] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
  EW_F  = 6,    // [6] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
  MG_ThermalTurbulence = 7,
  KTT_F = 7,    // [7] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
  EWTT_F = 8,   // [8] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
  MG_Boundary = 9,
  B_F = 9,   // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
  BX_F = 9,   // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
  BY_F = 10,  // [10]-> Displacement (quadratic (2), DS_EQUATIONS)
  BZ_F = 11,  // [11]-> Displacement (quadratic (2), DS_EQUATIONS)
  MG_Displacement = 9,
  SDS_F = 9,   // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
  SDSX_F = 9,   // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
  SDSY_F = 10,  // [10]-> Displacement (quadratic (2), DS_EQUATIONS)
  SDSZ_F = 11,  // [11]-> Displacement (quadratic (2), DS_EQUATIONS)
  MG_DA = 12,
  DA_F   = 12,  // [12]-> DA solver (quadratic (2), DA_EQUATIONS)
  // adjoint
  MG_AdjointNavierStoke = 13,
  NSA_F  = 13,  // [15]-> adjoint NS or FSI equations (Dimension)
  NSAX_F  = 13,  // [15]-> adjoint NS or FSI equations (Dimension)
  NSAY_F  = 14,  // [15]-> adjoint NS or FSI equations (Dimension)
  NSAZ_F  = 15,  // [15]-> adjoint NS or FSI equations (Dimension)
  MG_AdjointFluidStructure = 13,
  FSA_F  = 13,  // [15]-> adjoint NS or FSI equations (Dimension)
  FSAX_F  = 13,  // [15]-> adjoint NS or FSI equations (Dimension)
  FSAY_F  = 14,  // [15]-> adjoint NS or FSI equations (Dimension)
  FSAZ_F  = 15,  // [15]-> adjoint NS or FSI equations (Dimension)
  MG_AdjointDA = 16,
  DA_P   = 16,  // [13]-> DA solver (piecewise, DA_EQUATIONS)
  MG_AdjointTemperature = 17,
  TA_F   = 17,  // [14]-> Temp adjoint
  MG_AdjointTurbulence = 18,
  KA_F   = 18,  // [18] + 19 Adjoint turbulence
  EWA_F  = 19,
  //control
  MG_ColorFunction = 20,
  MG_Laplacian = 20,
  CO_F   = 20,   // [16]-> Color function for FSI equations
  MG_ControlTemperature = 21,
  CTRL_F  = 21,   // [16]-> Color function for FSI equations
  CTRLX_F  = 21,   // [16]-> Color function for FSI equations
  CTRLY_F  = 22,   // [16]-> Color function for FSI equations
  CTRLZ_F  = 23,   // [16]-> Color function for FSI equations
  // IMMERSED BOUNDARY
  IB_F = 24,
  IB_COL = 24,
  IB_VOL = 25,
  // TURBULENCE
  DIST = 26,
  MU_T = 27,
  ALPHA_T = 28,
  MG_CoupledTemperature = 29, // [29]-> Coupled Temperature optimality system 
  TCOUP_F=29,                  // [29]-> Coupled Temperature optimality system
  A_F = 30,                   // [30]-> 1D area solver 
  Q_F=31,                      // [31]-> 1D flow rate solver
  // ANISOTROPIC TURBULENCE
    MG_ReynoldsStressTensor = 32,
    TAU_F = 32,
    TAUXX_F=32,
    TAUXY_F=33,
    TAUYY_F=34,
    MG_ReynoldsHeatFlux =35,
    THF_F=35,
    THFX_F=35,
    THFY_F=36,
    THFZ_F=37
  };


// ==========================================================
/// Equation handler  
/*! This class handles the construction of equation objects. Equations to activate
 * are read from input parameter file Equations.in. After their constructor has been called,
 * their pointers are added #FEMUS::_mg_equations_map equation map. Finally all equations share
 * their solution by filling MGSolDA::_data_eq structure for each problem object */
class EquationsMap {

    
  public:
    /// Vector containing equation numbers of equations that must be initialized
    std::vector<FIELDS> _myproblemP;
    
    /// Map associating equation names to corresponding number coming from #FIELDS enumeration
    std::map<std::string, FIELDS> _map_str2field;
    
    /// Map associating names and numbers of equation to be solved
    std::map<std::string, int> _EquationsToAdd;
    
    /// Map for turbulence model choice
    std::map<std::string, std::string> _TurbulenceModel;

    /// Global amount of equations with piecewise solved variables
    int _PieceEq;
    
    /// Global amount of equations with linear solved variables
    int _LinearEq;
    
    /// Global amount of equations with quadratic solved variables
    int _QuadEq; 
 
    /// Finite element approximations for solved variables
    /*! Array containing number of piecewise (_nvars[0]), linear (_nvars[1]) and quadratic
     * (_nvars[2]) variables that solved for each initialized equation */    
    int _nvars[3];

    /*! Standard constructor of EquationsMap class. When the class object is initialized
     * functions     #EquationsMap::Fill_FIELD_map(), #EquationsMap::ReadEquationsToAdd()
     * and #EquationsMap::Fill_pbName() are called */
    EquationsMap();
    ~EquationsMap() {};
    
    /*! This function fills #EquationsMap::_EquationsToAdd map with the names of 
     * equations that will be activated. Equations are read from Equations.in 
     * input file */
    void ReadEquationsToAdd();
    
    /*! This function fills #EquationsMap::_TurbulenceModel map with info read from
     Turbulence.in input file*/
    void ReadTurbulenceInfo();

    /*! This function fills #EquationsMap::_map_str2field map: equation names are 
     * associated to corresponding number of #FIELDS enumeration */
    void Fill_FIELD_map ( ) ;

    /*! This function fills #EquationsMap::_myproblemP vector with equation numbers.
     * Equation numbers are read from #EquationsMap::_map_str2field map using entries
     * stored in #EquationsMaps::_EquationsToAdd map */
    void Fill_pbName ( );

    /*! This function calls the construction of equations objects and add their 
    * pointers to #FEMUS::_mg_equations_map. Activated equations share their 
    * solutions through #MGSolDA::_data_eq structure */
    void FillEquationMap ( EquationSystemsExtendedM & EqMap );
    
    /*! This functions fills #MGSolDA::_data_eq structure for each activated equation
     * in order to let equations share their solutions */
    void setProblems ( MGSolBase *& problem );
    
    /// Initializaiton of Navier Stokes system
    void initNavierStokes ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Fluid Structure Interaction system
    void initFluidStructure ( EquationSystemsExtendedM & EqMap );
    
     /// Initializaiton of boundary system
    void initBoundary ( EquationSystemsExtendedM & EqMap ) ;
    
    /// Initializaiton of solid displacement system
    void initDisplacements ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Temperature equation
    void initTemperature ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of dynamical turbulence system
    void initDynamicTurbulence ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of thermal turbulence
    void initThermalTurbulence ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of color-curvature system
    void initColor ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Adjoint Navier Stokes system
    void initAdjointNavierStokes ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Adjoint Fluid Structure Interaction system
    void initAdjointFSI ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Adjoint temperature equation
    void initAdjointTemperature ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Adjoint dynamical turbulence system
    void initAdjointDynamicTurbulence ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Control temperature system
    void initControlTemperature ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of Temperature control optimality system
    void initCoupledTemperature ( EquationSystemsExtendedM & EqMap );

    /// Initializaiton of Immersed Boundary system
    void initImmersedBoundary ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of one-dimensional Area system
    void initAreaMono ( EquationSystemsExtendedM & EqMap );
    
    /// Initializaiton of one-dimensional Flowrate system
    void initFlowrateMono ( EquationSystemsExtendedM & EqMap );
    
    /// Initialization of Reynolds Stress Tensor
    void initReynoldsStressTensor (EquationSystemsExtendedM & EqMap) ;
        
    /// Initialization of Reynolds Heat Flux
    void initReynoldsHeatFlux (EquationSystemsExtendedM & EqMap) ;

  };



#endif
