#include "Equations_conf.h"

// ============================================
#ifdef DS_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "MGGeomEl.h"
#include "MGSolverDS.h"

// configuration files -----------
#include "MGEquationsSystem.h"
#include "MGFE_conf.h"
#include "Printinfo_conf.h"

// local include -------------------
#include "MGSystem.h"

void MGSolDS::ic_read(
    int /*bc_gam*/, int /*mat_gam*/, double /*xp*/[], int iel,
    double u_value[]) {  // =======================================
  u_value[0] = 0.;
  return;
}

// ============================================================================
/// This function  defines the boundary conditions for the DS system:
void MGSolDS::bc_read(
    int bc_gam, int mat_gam,
    double xp[],    // xp[] is the NON-DIMENSIONAL node coordinates
    int bc_Neum[],  // normal
    int bc_flag[]   // boundary condition flag
) {                 // =========================================================================

  // channel_fsi1t.med ===============================================
  bc_Neum[0] = fix_disp0;
  bc_flag[0] = 0;
  bc_Neum[0] = _DS_parameter._map_DSgroup[bc_gam];

  return;
}
#endif
