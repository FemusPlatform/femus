// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

#include "MGSolverNS.h"       // Navier-Stokes class header file
#include <fstream>
#include "MGUtils.h"

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

  double vel[3];
  double mod = 0.1;
  double notaL = 0.00001;

  double pi_2 = acos(0.);

  vel[0] = sin(pi_2*4/3)*mod +notaL;
  vel[1] = cos(pi_2*4/3)*mod +notaL;
  vel[2] = 0.;
  
if(_Coupled==1){
  u_value[0] = vel[0];
  u_value[1] = vel[1];
  u_value[DIMENSION-1] = (_nNSdim==2)? vel[1]:vel[2];
  u_value[DIMENSION] = 0.0  * 0.01;
}
else if (_Coupled==0){
  u_value[0] = vel[_dir];  
}
  
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
) {
  double ILref = 1./_lref;
  bc_Neum[0] = outflow;
  bc_Neum[0] = _NS_parameter._map_NSgroup[bc_gam];

  // The low and right bottom corners are not correct with groups only
  return;
}



#endif  //ENDIF NS_EQUATIONS
