#ifndef __user_RANS_thermal_h__
#define __user_RANS_thermal_h__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_THERMAL_EQUATIONS
// ==================================

// T boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
enum bound_condRANS_T {
  // ========================================================================
  TKwall0 = 0,        ///< 0 = Dirichlet homogeneus     (T=0)
  TKwall = 1,         ///< 1 = Dirichlet non homogeneus for wall boundaries
  TKinlet = 2,        ///< 2 = Dirichlet non homogeneus for inlet boundaries
  TKinsulation = 10,  ///< 10 = Neuman homogeneus         (dT.n=0)
  TKheat_flux = 11,   ///< 11 = Neuman nonhomogeneus      (dT.n=q_0)
  TKrobinT = 12,      ///< 12 = Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
  TKLext_conv = 13,   ///< 13 = Turbulence      (dT.n=tau_0^2=beta*ug(u))
  TKsimmetry = 20,    ///< 20 = Simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
  TKWallFuncGrad = 21
  // ========================================================================
};

#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
