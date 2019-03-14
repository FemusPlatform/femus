// ======================================================================================
//                  NAVIER-STOKES INITIAL AND BOUNDARY CONDITIONS
// ======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

#include "MGSolverNS.h"      
#include <fstream>
#include "MGUtils.h"


// INITIAL CONDITION ---------------------------------------------------------------------
void MGSolNS::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]
) {

  double vel[3];

  vel[0] = 0.;
  vel[1] = 0.;
  vel[2] = 0.;
  
if(_Coupled==1){
  u_value[0] = vel[0];
  u_value[1] = vel[1];
  u_value[DIMENSION-1] = (_nNSdim==2)? vel[1]:vel[2];
  u_value[DIMENSION] = 0.0;
}
else if (_Coupled==0){
  u_value[0] = vel[_dir];  
}
  
  return;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolNS::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          
  int bc_Neum[], 
  int bc_flag[]         
) {

  bc_Neum[0] = outflow;
  bc_Neum[0] = _NS_parameter._map_NSgroup[bc_gam];

  return;
}


#endif 
