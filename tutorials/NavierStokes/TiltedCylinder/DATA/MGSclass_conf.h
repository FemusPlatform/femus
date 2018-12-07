#include <Equations_conf.h>

//
//
/*===================================== DA Equation =======================================================*/
//
#ifdef DA_EQUATIONS
// ===============================->0.4268 ==============
#ifndef __mgsnsc_h__
#define __mgsnsc_h__



//  -----------------------------
// NAVIER-STOKES ADVECTION term
// --------------------------------
// A) Stokes flow   ADV 0. B)  Navier-Stokes ADV 1.
#define ADV 1.
//  Navier-Stokes  A) Picard iteration  ADV1 (0.) B) Newton ADV 1.
#define ADV1 0.
//  Navier-Stokes  stab  +0.5*(div u,u.v)
#define STAB 1.
//


#endif

#endif /* End DA Equation ---------------------------------------------------------------------------------*/
//
//
/*===================================== FSI Equation ======================================================*/
//
#ifdef FSI_EQUATIONS
// =============================================

#ifndef __mgsnscfsi_h__
#define __mgsnscfsi_h__


//
// MULTIGRID PARAMETERS
// -------------------------------
// Navier-Stokes solver type
#define SOLVER_FSI GMRESM  // options -> GMRESM BICGSTABM

// Pressure solver type (projection method only)
#define SOLVER_FSIP CGM  // options -> GMRESM CGM BICGSTABM


// ------------------------
//   SOLID MODEL
// ----------------------------
// geometric non-linearity
 #define NL_GEOM  (1)

// define penalty  (only with FSIP_EQUATIONS==1)
#define  PENALTY_FSI (100.)



 // ---------------------------
//   SOLID-FLUID REGIONS
// -----------------------------
//  #define SOLID 0
//  #define STIFF 10



//  --------------------------
// 3D NAVIER-STOKES ADVECTION term
// ---------------------------------
// A) Stokes flow    ADVPIC_SM 0. ADVNEW_SM=0
// B)  Navier-Stokes ADVPIC_SM 1. ADVNEW_SM={0,1}
#define ADVPIC_FSI 0.
//  Navier-Stokes  nonlinear iterations
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_SM=1, ADVNEW_SM=0
// B) Newton iteration  ADVPIC_SM=1, ADVNEW_SM=1
#define ADVNEW_FSI 0.


// -----------  stabilization ------------------
// --------------------------------------------
//  Navier-Stokes  stab=STAB_SM*0.5*(div u,u.v)
// #define STAB_FSI 0.
//  Navier-Stokes  compressibility=KOMP_SM*dp/dt
#define KOMP_FSI 1.e-20
//  Navier-Stokes  compressibility  upwind=Re+UP_WIND_SM*v^2
#define UP_WIND_FSI (1.)
//  Navier-Stokes  penalty =LAMBDA grad div u (solid)
// #define LAMBDA (2100)
// Navier-Stokes supg
// #define SUPG_FSI (0.)
// Navier-Stokes  c -> antisymmetric (0.5)
// #define ADV_ASYM 0.


// #define LQ (2)


// Crank-Nicolson first order 0. 2nd order 0.5
// #define CN_TIME 0.

//  boundary integral
// #define P_0 (0.)

// turbulent
// #define MU_LOW (1.e-12)
// #define MU_TOP (1.e+12)


// #define SQCMU (0.3)
// #define YPLUS (1.)
// #define ALPHA0 (1.)



// P solver ------------------------------------

// #define SOLVERT VANKATM ---------------------
// #define SOLVERFSIP GMRESM


// temperature lows -----------------------------
//  #define  CONST 1
// // #define densityT(x) (1.)
// // #define cipiT(x) (1.)
// // #define kappa(x)  (1.)
//
// // quadratic QL=2 linear QL=1 -------------------
// #define LQ (2)
// #define LQ_P (1)
// // boundary pressure --------------------------
// #define P_BD (1.)


#endif
#endif /* End FSI Equation --------------------------------------------------------------------------------*/
//
//
/*===================================== DS Equation =======================================================*/
//
#ifdef DS_EQUATIONS
// =============================================


#endif /* End DS Equation ---------------------------------------------------------------------------------*/
//
//
//
//
/*===================================== Thermal Turbulence ================================================*/
//
#ifdef TTBK_EQUATIONS

#ifndef  __mgsttbkconf_h__
#define  __mgsttbkconf_h__

// turbulence  ==========================
#define ADVE 1.

// #define SOLVERT VANKATM =======================
#define SOLVERTBK GMRESM

// temperature lows ================================
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)

// #define  LINEAR 1 quad (2)
//  #define LQ_TB  (2)

// #define LES (100.)

#define UP_WIND_TTK (1.)

//boundary
// #define SQCMU (0.3)
// #define KAPPA_VK (0.4)
//
/*-------------------------------------- For Nagano K-E ----------------------------------------------------*/
//
#if ((TTBK_EQUATIONS/2)==1)

#define MAX_TAUT (1.e+8)                 /* Limit for maximum thermal characteristic time */
#define MIN_TAUT (1.e-8)                 /* Limit for minimum thermal characteristic time */

#define SIGMA_EH (0.714285714286)
#define SIGMA_KH (0.714285714286)

#define CD1 (1.0)
#define CP1 (0.925)
#define CD2 (1.9)
#define CP2 (0.9)

#define CMU (0.09)

#endif /* End Nagano K-E -----------------------------------------------------------------------------------*/
//
/*-------------------------------------- The Thermal K-W turbulence models ---------------------------------*/
//
#if ((TTBK_EQUATIONS/2)==2)
//
#define MAX_TAUT (1.e+14)                /* Limit for maximum thermal characteristic time */
#define MIN_TAUT (1.e-14)                /* Limit for minimum thermal characteristic time */
#define CMU (0.09)
#define BETA (0.09)


/*-------------------------------------- For SST -----------------------------------------------------------*/
// #define SST (1)
#ifdef SST
#define SIGMA_WH (0.714285714286)
#define SIGMA_KH (0.714285714286)
#define CD1 (1.4)
#define CP1 (1.1)
#define CD2 (0.8)
#define CP2 (0.6)
#endif /* End SST ------------------------------------------------------------------------------------------*/
//
/*-------------------------------------- For Nagano K-w ----------------------------------------------------*/
//
// #define NAGANO (1)
#ifdef NAGANO
#define SIGMA_WH (0.714285714286)
#define SIGMA_KH (0.714285714286)
#define CD1 (0.1)
#define CP1 (0.025)
#define CD2 (1.9)
#define CP2 (0.9)
#endif /* End Nagano K-W -----------------------------------------------------------------------------------*/

#endif /* End K-W models -----------------------------------------------------------------------------------*/

#endif
#endif /* End TTBK_EQUATIONS ------------------------------------------------------------------------------*/
