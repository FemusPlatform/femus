// =======================================================================================
//                 INITIAL AND BOUNDARY CONDITION FOR PRESSURE EQUATION
// =======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

#include "MGSolverP.h"
#include "UserP.h"

// INITIAL CONDITIONS --------------------------------------------------------------------
void MGSolP::ic_read (
    int bc_gam,
    int bc_mat,
    double xp[],
    int iel,
    double u_value[] )
{
    u_value[0] = ( 1.331 - xp[1] ) * 9.81 ;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolP::bc_read (
    int bc_gam,
    int bc_mat,
    double xp[],int bc_Neum[],int bc_flag[] )
{

    bc_Neum[0] = vel_fix;
    bc_Neum[0] = _P_parameter._map_NSgroup[bc_gam];

    return;
} 

#endif


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
