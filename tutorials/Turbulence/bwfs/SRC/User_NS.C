// =======================================================================================
//                 INITIAL AND BOUNDARY CONDITION FOR VELOCITY EQUATION
// =======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

#include "MGSolverNS.h"
#include "MGUtils.h"

// INITIAL CONDITIONS
// --------------------------------------------------------------------
void MGSolNS::ic_read(int bc_gam,  /**<1 */
                      int bc_mat,  /**<2 */
                      double xp[], /**<3 */
                      int iel,     /**<4 */
                      double u_value[] /**< 5*/) {

  double walldist = min(xp[0], 0.121 - xp[0]);
  double vert_vel = 0.;
  if (_NS_parameter._FlatProfile == 0 && xp[0] < 0.121 - BDRY_TOLL &&
      xp[1] < 1.e-10) {
    vert_vel = 0.002907 * _mgutils._TurbParameters->Musker(walldist, 0.002907);
  }

  double vel[3];

  vel[0] = 0.;
  vel[1] = vert_vel;
  vel[2] = 0.;

  if (_Coupled == 1) {
    u_value[0] = vel[0];
    u_value[1] = vel[1];
    u_value[DIMENSION - 1] = (_nNSdim == 2) ? vel[1] : vel[2];
    u_value[DIMENSION] = 0.0 * 0.01;
  } else if (_Coupled == 0) {
    u_value[0] = vel[_dir];
  }

  return;
}

// BOUNDARY CONDITIONS
// -------------------------------------------------------------------
void MGSolNS::bc_read(int bc_gam,
                      int bc_mat,
                      double xp[], int bc_Neum[], int bc_flag[]) {

  bc_Neum[0] = outflow;
  bc_Neum[0] = _NS_parameter._map_NSgroup[bc_gam];

  return;
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
