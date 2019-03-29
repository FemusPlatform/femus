// =======================================================================================
//  INITIAL AND BOUNDARY CONDITION FOR EDDY VISCOSITY AND THERMAL DIFFUSIVITY EQUATIONS
// =======================================================================================
#include "Equations_conf.h"

#ifdef _TURBULENCE_

#include "MGSolverTURB.h"       // Navier-Stokes class header file
#include "UserTURB.h"

// INITIAL CONDITIONS --------------------------------------------------------------------
void MGSolTURB::ic_read
( int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]
)
{

    u_value[0] = 0.;
    return;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolTURB::bc_read ( int bc_gam, int bc_mat, double xp[], int bc_Neum[], int bc_flag[] )
{
    
    bc_Neum[0] = 1;

    for ( int k = 0; k < _WallGroupsIDs.size(); k++ )
        if ( bc_gam == _WallGroupsIDs[k] ) {
            bc_Neum[0] = 0;
        }
        
    return;
}
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
