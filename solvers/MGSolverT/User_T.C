// =======================================================================================
//                   ENERGY EQUATION INITIAL AND BOUNDARY CONDITIONS
// =======================================================================================
#include "Equations_conf.h"

#ifdef T_EQUATIONS // 3D-2D Energy equation

#include "MGSolverT.h"       
#include "UserT.h"

// INITIAL CONDITION ---------------------------------------------------------------------
void MGSolT::ic_read(
    int bc_gam,
    int bc_mat,
    double xp[],
    int iel,
    double u_value[])
{
  u_value[0] =  573;
  return;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolT::bc_read(
    int bc_gam,
    int bc_mat,
    double xp[],
    int bc_Neum[],
    int bc_flag[]) 
{

  bc_Neum[0] = insulation;
  bc_Neum[0] = _T_parameter._map_Tgroup[bc_gam];
  
  return;
} 

#endif
