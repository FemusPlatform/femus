#ifndef __userP_h__
#define __userP_h__
// ============================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
#include "Solvertype_enum.h"

#include <map>
#include <vector>

class MGUtils;

// ============================================================================
// P boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
// P boundary conditions===========================================================================
enum bound_cond_p {

  // ================================================================================================
  // outflowp0 = 0= Dirichlet homogeneus     \f$ p= 0   \f$    (p=0)
  // outflowp  = 4= Dirichlet nonhomogeneus  \f$ p= p_0 \f$   (p=p_0)
  //
  // vel_fix   =10= Neuman homogeneus or simmetry \f$ \nabla p \cdot \widehat{n}= 0 \f$  (dp.n=0)
  // interiorp =11= Neuman nonhomogeneus \f$ \nabla p \cdot \widehat{n}= \Delta p_0 \f$   (dp.n=dp0)
  // ================================================================================================
  // Dirichlet
  outflowp0 = 0,  ///< 0=Dirichlet homogeneus     \f$ p= 0   \f$    (p=0)
  outflowp = 4,   ///< 1=Dirichlet nonhomogeneus  \f$ p= p_0 \f$    (p=p_0)
                  // Neuman
  vel_fix = 10,   ///< 10= Neuman homogeneus or simmetry \f$ \nabla p \cdot \widehat{n}= 0 \f$ (dp.n=0)
  interiorp = 11  ///< 11= Neuman nonhomogeneus \f$ \nabla p \cdot \widehat{n}=\Delta p_0 \f$  (dp.n=dp0)
};

#endif
#endif  //    __userP_h__
        //    <-----------------------------------------------------------------------------------------
