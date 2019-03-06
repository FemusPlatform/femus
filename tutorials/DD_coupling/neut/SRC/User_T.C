#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================

// class files --------------------------------------------------------------------------

#include "MGSolverT.h"       // Navier-Stokes class header file
#include "UserT.h"

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
/// This function generates the initial conditions for the energy equation:
/// \ingroup  user_function ic_read
// =================================================
void MGSolT::ic_read(int bc_gam,int bc_mat, double xp[],int iel,double u_value[]) {
// xp[]=(xp,yp) u_value[]=(u,v,p)
  u_value[0] =  573;//*(1.-xp[1]);
  return;
}
// ========================================================
/// This function  defines the boundary conditions for the energy equation:
/// \ingroup  user_function bc_read
void MGSolT::bc_read(int bc_gam,int bc_mat, double xp[],int bc_Neum[],int bc_flag[]) {
//   double ILref = 1./_lref;
   bc_Neum[0] = insulation;
//   bc_Neum[0] = _T_parameter._map_Tgroup[bc_gam];
  if(xp[2]<0.) bc_Neum[0] = Twall;
  
  
//   if(xp[1] < -0.005816 + BDRY_TOLL && xp[0] < BDRY_TOLL) bc_Neum[0] = Twall0;
} // end boundary conditions ==========================
/** @} */

#endif
