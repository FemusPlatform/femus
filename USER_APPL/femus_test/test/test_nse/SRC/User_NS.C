// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS
// #if NS_EQUATIONS==1
// ======================================================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file


// config file --------------------------------------------------------------------------
#include "MGGeomEl.h"        // Geometrical element
#include "MGFE_conf.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ======================================================================================
/**     \addtogroup user_function User functions
 * @{
 */

/// \ingroup  user_function ic_read
// ======================================================================================
/// This function generates the initial conditions for the NS system:
void MGSolNS::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]
) {// ===================================================================================

  double ILref = 1./_lref;

  // =======================================================================================
  //    boundary conditions box file plane_ch_rotate.med ================================== 
  u_value[0]=0.; u_value[1]=0.;   u_value[DIMENSION]=0.1;

  // 2D CASE =========================================================================================

   double mod = 0.5;  
//    if(fabs(xp[0] - 0.025)<1.e-5 && fabs(xp[1] + 0.0375)<1.e-5 && fabs(xp[2] + 0.021651)<1.e-5) mod= 2.;
   
   u_value[0] =-0.5*mod;  u_value[1] =0.75*mod;  u_value[2] =0.433013*mod; 
   
        
  return;
}

/// \ingroup  user_function bc_read
// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolNS::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], 	// normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0  ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 ->  tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;
   bc_Neum[0] = outflow;
   bc_Neum[0] = _NS_parameter._map_NSgroup[bc_gam];

// ############################################
//    if(bc_gam==20)   bc_Neum[0]=44;
  // ########################################### 
   // The low and right bottom corners are not correct with groups only
//   if(xp[1]<-0.0683) bc_Neum[0]=velocity; // low bottom corner
//    if(xp[0]>0.0683) bc_Neum[0]=velocity; // right bottom corner
//     if(xp[1]>0.0683) bc_Neum[0]=44;
//  if(xp[0]<-0.0683) bc_Neum[0]=44; // low bottom corner
  return;
}









#endif  //ENDIF NS_EQUATIONS
