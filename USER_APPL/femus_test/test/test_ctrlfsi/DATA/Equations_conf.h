#ifndef _equationsconfw
#define _equationsconfw

// ==================================================
 //                   Fluid mechanics
 // ==================================================

 //  Navier-Stokes -----------------------------------
 //  NS projection -> 0 coupled -> 1  split -> 2
 //  -------------------------------------------------
// #define  WALL_FUNC_APP
//    #define  NS_EQUATIONS (1)
//   #define  NSA_EQUATIONS (1)
// #define  Nat_Conv   // For Natural Convection

 //  Turbulence --------------------------------------
 //  Length   ->0     Spalart-Allmaras -> 1
 //  KE  coupled -> 3  split -> 2
 //  KW  coupled -> 5  split -> 4
 // --------------------------------------------------
//   #define  TBK_EQUATIONS  (4)
//    #define  LOG_W
//   #define  LOG_K   
 // LOG_W is for logarithmic model of k-w

 // ==================================================
 //               ENERGY
 // ==================================================

 //  Temperature equation ---------------------------
 //  -------------------------------------------------
//     #define T_EQUATIONS  (1)

 //  Energy Turbulence   ----------------------------
 //  KE  coupled -> 3  split -> 2
 //  KW  coupled -> 5  split -> 4
 //  -------------------------------------------------
// #define  TTBK_EQUATIONS  (4)
// #define  LOG_WT
// #define  LOG_KT  
 // LOG_WT is for logarithmic model of kt-wt

 // ==================================================
 //          STRUCTURAL MECHANICS
 // ==================================================
 //  SM projection -> 0 coupled -> 1  split -> 2
 //  -------------------------------------------------
 //  #define SM_EQUATIONS


 // ==================================================
 //               NEUTRONICS
 // ==================================================
 //  NEUTRON DIFFUSION
 //  -------------------------------------------------

 //  #define NEU_EQUATIONS


 // ==================================================
 //                ELECTROSTATICS
 // ==================================================
 //  POTENTIAL EQUATION
 //  -------------------------------------------------

 //  #define V_EQUATIONS


 // ==================================================
 //             FLUID STRUCTURE  FSI
 // ==================================================

 //  FSI equations
 //  FSI projection -> 0 coupled -> 1  split -> 2
 //  -------------------------------------------------
   #define FSIA_EQUATIONS  (1)
   #define FSI_EQUATIONS (1)
//    #define COLOR_FLAG  (1)
   #define CTRL_EQUATIONS (2)
 // ==================================================
 //                   MHD
 // ==================================================

 // FSI equations
 // FSI projection -> 0 coupled -> 1  split -> 2
 //  -------------------------------------------------
 // #define MHD_EQUATIONS

 // ==================================================
 //                 TWO_PHASE
 // ==================================================
 //  --------------------------------------
 // #define TWOPH_EQUATIONS


//   #define DA1_EQUATIONS      // solver DA3D, test system
 // ================================================
 // =========== only if you know how to do =========
 // ================================================
// #define DA_EQUATIONS         // solver DA3D, test system

 #ifdef FSI_EQUATIONS
 /// Displacement
 #define  DS_EQUATIONS
 #if (FSI_EQUATIONS%2==0)       // projection method
  #define FSIP_EQUATIONS (1)     // need separated P equations
 #endif
 #endif

 #ifdef SM_EQUATION
 #if (SM_EQUATIONS%2==0)       // projection method
 #define SMP_EQUATIONS (1)     // need separated P equations
 #endif
 #endif
 // ================================================


#endif  // end file _equationsconf_
