// ======================================================================================
//                        PRESSURE INITIAL AND BOUNDARY CONDITIONS
// ======================================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include "MGSolverP.h"
#include <fstream>
#include "MGUtils.h"
#include "UserP.h"
#include "Pparameters.h"

// INITIAL CONDITION ---------------------------------------------------------------------
void MGSolP::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]) 
{
  u_value[0] = 0.;
  return;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolP::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int bc_Neum[],
  int bc_flag[]) 
{
  bc_Neum[0] = vel_fix;
  bc_Neum[0] = _P_parameter._map_NSgroup[bc_gam];
  return;
} 

#endif


