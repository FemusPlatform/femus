#ifndef __userT_h__
#define __userT_h__

#include "Equations_conf.h"
// ===================================
#ifdef T_EQUATIONS
// ==================================

// T boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
enum bound_condT {
  // ========================================================================
  Twall0 = 0,       ///< 0= Dirichlet homogeneus     (T=0)
  Twall = 1,        ///< 1= Dirichlet nonhomogeneus  (T=T_0)
  insulation = 10,  ///< 10= Neuman homogeneus         (dT.n=0)
  heat_flux = 11,   ///< 11= Neuman nonhomogeneus      (dT.n=q_0)
  robinT = 12,      ///< 12= Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
  ext_conv = 13,    ///< 13= Turbulence      (dT.n=tau_0^2=beta*ug(u))
  simmetry = 20     ///< 20= Simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
                    // ========================================================================
};

#endif
#endif
