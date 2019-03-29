// =======================================================================================
//            INITIAL AND BOUNDARY CONDITION FOR DYNAMICAL TURBULENCE EQUATIONS
// =======================================================================================
#include "Equations_conf.h"

#ifdef RANS_EQUATIONS 

#include "MGSolverRANS.h"       
#include "User_RANS.h"
#include "MGUtils.h"
#include "MGMesh.h"

// INITIAL CONDITIONS --------------------------------------------------------------------
void MGSolRANS::ic_read (
     int bc_gam,
     int bc_mat,
     double xp[],
     int iel,
     double u_value[]
)
{

     double k_in, w_in, wall_dist, initial_value[2];
     wall_dist = _mgmesh._dist[iel];

     _mgutils._TurbParameters->DynTurInitValues ( k_in, w_in, wall_dist, _FlatProfile );

//      if ( xp[1]<0.121 ) {
//           double walldist = min ( xp[0], 0.121 - xp[0] ) + 1.e-5;
//           double utau = 0.002907;
//           _mgutils._TurbParameters->DynTurNearWallValues ( k_in, w_in, wall_dist, utau );
//      }

     initial_value[0] = k_in;
     initial_value[1] = w_in;

     // KAPPA
     if ( _dir == 0 ) {
          u_value[0] = initial_value[0];
     }

     // OMEGA
     if ( _dir == 1 ) {
          u_value[0] = initial_value[1];
     }
     
     return;
}



// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolRANS::bc_read (
     int bc_gam,
     int bc_mat,
     double xp[],
     int bc_Neum[],
     int bc_flag[]
)
{

     bc_Neum[0] = Kinsulation;

     if ( _dir == 1 ) {
          bc_Neum[0] = _RANS_parameter._map_DTWgroup[bc_gam];
     }

     if ( _dir == 0 ) {
          bc_Neum[0] = _RANS_parameter._map_DTKgroup[bc_gam];
     }
     return;
}

#endif
// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 
