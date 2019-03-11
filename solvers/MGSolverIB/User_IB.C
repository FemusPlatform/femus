#include "Equations_conf.h"

// ============================================
#ifdef IMMERSED_BOUNDARY // 3D-2D Energy equation
// ============================================
#include "MGSolverIB.h"       // Navier-Stokes class header file
#include "UserIB.h"
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library


// ======================================================================================

/**     \addtogroup user_function User functions
 * @{
 */

// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
/// \ingroup  user_function ic_read
// =================================================
void MGSolIB::ic_read(int bc_gam,int bc_mat, double xp[],int iel,double u_value[]) {
//    boundary conditions box plane_ch_rotate.med ==================================
  u_value[0] = 1.;
}



// ========================================================
/// This function  defines the boundary conditions for the energy equation:
/// \ingroup  user_function bc_read
void MGSolIB::bc_read(int bc_gam,int bc_mat, double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================

bc_Neum[0]=1;
for(int k=0; k<_WallGroupsIDs.size(); k++){
  if(bc_gam == _WallGroupsIDs[k]){
       bc_Neum[0]=0;
  }
}
  return;
} // end boundary conditions ==========================
/** @} */

#endif
