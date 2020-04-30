// ======================================================================================
// --------------   NAVIER-STOKES system [FS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef FSI_EQUATIONS
// ======================================================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"  // Navier-Stokes class conf file
#include "MGSolverFSI.h"    // Navier-Stokes class header file

// config file --------------------------------------------------------------------------
#include "MGFE_conf.h"       // FEM approximation
#include "MGGeomEl.h"        // Geometrical element
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
#include "MGEquationsSystem.h"  // Equation map class
#include "MGMesh.h"             // Mesh class
#include "MGSystem.h"           // System class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>  // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"    // algebra dense matrices
#include "dense_vectorM.h"    // algebra dense vectors
#include "linear_solverM.h"   // algebra solvers
#include "numeric_vectorM.h"  // algebra numerical vectors
#include "sparse_matrixM.h"   // algebra sparse matrices
// ======================================================================================
/**     \addtogroup user_function User functions
 * @{
 */

/// \ingroup  user_function ic_read
// ======================================================================================
/// This function generates the initial conditions for the NS system:
void MGSolFSI::ic_read(
    int bc_gam, int bc_mat, double xp[], int iel,
    double
        u_value[]) {  // ===================================================================================
  // =======================================================================================
  //    boundary conditions box file channel_fsi1t.med ==================================
  // =======================================================================================
  //    boundary conditions box file plane_ch_rotate.med ==================================
  u_value[0] = 0.;
  u_value[1] = 0.;
  double u = 1.;
  //   if(xp[0]<0.0001)
  if (bc_mat == 2) u_value[0] = 1.5 * u * 4 / 0.1681 * xp[1] * (0.41 - xp[1]);

  double pref = _refvalue[DIMENSION];
  u_value[DIMENSION] = 15. - 3 * (xp[0] - 0.25);

  return;
}

/// \ingroup  user_function bc_read
// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolFSI::bc_read(
    int bc_gam, int bc_mat,
    double xp[],    // xp[] is the NON-DIMENSIONAL node coordinates
    int bc_Neum[],  // normal
    int bc_flag[]   // boundary condition flag
) {                 // ===================================
  bc_Neum[0] = wall;
  bc_flag[0] = 0;
  bc_Neum[0] = _FSI_parameter._map_FSIgroup[bc_gam];

  return;
}

#endif  // ENDIF FSI_EQUATIONS
