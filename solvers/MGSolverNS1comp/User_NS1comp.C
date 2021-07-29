// =======================================================================================
//                 INITIAL AND BOUNDARY CONDITION FOR VELOCITY EQUATION
// =======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS


#include "MGSolverNS_1comp.h"
#include "MGSolverNS.h"
#include "MGUtils.h"
#include "MGMesh.h"

// INITIAL CONDITIONS

void MGSolNS_1comp::ic_read(
    int bc_gam,  /**<1 */
    int bc_mat,  /**<2 */
    double xp[], /**<3 */
    int iel,     /**<4 */
    double u_value[] /**< 5*/) {
  u_value[0] =0.;
  return;
}

// BOUNDARY CONDITIONS
// -------------------------------------------------------------------
void MGSolNS_1comp::bc_read(int bc_gam, int bc_mat, double xp[], int bc_Neum[], int bc_flag[]) {
  bc_Neum[0] = outflow;
  bc_Neum[0] = _NS_parameter._map_NSgroup[bc_gam];
  return;
}

#endif
