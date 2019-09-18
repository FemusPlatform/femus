#ifndef __userRANS_h__
#define __userRANS_h__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_EQUATIONS
// ==================================

// T boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
enum bound_condK {
  Kinlet = 2,
  Kwall0 = 0,  ///< 0= Dirichlet homogeneus     (T=0)
  Kwall = 1,   ///< 1= Dirichlet nonhomogeneus  (T=T_0)
  Kinit = 4,
  Kinsulation = 10,  ///< 10= Neuman homogeneus       (dT.n=0)
  Kheat_flux = 11,   ///< 11= Neuman nonhomogeneus    (dT.n=q_0)
  KrobinT = 12,  ///< 12= Robin                   \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
  KLext_conv = 13,  ///< 13= Turbulence              (dT.n=tau_0^2=beta*ug(u))
  KrobinTb = 16,   ///< 12= Robin                   \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
  Ksimmetry = 20,  ///< 20= Simmetry                \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
  KWallFuncGrad = 21
};
/*!  \defgroup RANS_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup RANS_param
// ============================================================================

#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
