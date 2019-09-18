#include "Equations_conf.h"

#ifdef _TURBULENCE_

#include "MGSolverTURB.h"  // Navier-Stokes class header file
#include "UserTURB.h"

/**     \addtogroup user_function User functions
 * @{
 */

// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
/// \ingroup  user_function ic_read
// =================================================
void MGSolTURB::ic_read(int bc_gam, int bc_mat, double xp[], int iel, double u_value[]) {
  u_value[0] = 0.;
  return;
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
/// \ingroup  user_function bc_read
void MGSolTURB::bc_read(int bc_gam, int bc_mat, double xp[], int bc_Neum[], int bc_flag[]) {
  // =================================================

  bc_Neum[0] = 1;

  for (int k = 0; k < _WallGroupsIDs.size(); k++)
    if (bc_gam == _WallGroupsIDs[k]) { bc_Neum[0] = 0; }

  return;
}  // end boundary conditions ==========================
  /** @} */

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
