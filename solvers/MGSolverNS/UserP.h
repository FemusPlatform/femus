
#ifndef __userP_h__
#define __userP_h__
// ============================================================================
#include "Equations_conf.h"
#if NS_EQUATIONS%2==0

// ============================================================================

/*!  \defgroup T_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup P_param
// ============================================================================
/// Class P_parameter for Pressure equation
/// (use only pressure-velocity splitting formulation) 
class P_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int AXISYM  /**< axisymmetry  */;
    // non linear -------------------------------------------------------------
    int    NL_ITER;/**< NON LIN IT */    int    NL_ITER0;/**< INIT NON LIN IT */
    double NL_TIME0; // initial time for  non linear regime

public:
    // constructor --------------------------------------------------------------
    P_param() {
        AXISYM=0;
        // non linear ------------------------------------------------------------
        NL_ITER=0;/**< NON LIN IT */ NL_ITER0=0;/**< INIT NON LIN IT */
        NL_TIME0=-.0001;
    }

    inline void set_AXISYM ( int val ) {
        AXISYM=val;
    }

};





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
    outflowp0 =0,///< 0=Dirichlet homogeneus     \f$ p= 0   \f$    (p=0)
    outflowp  =4,///< 1=Dirichlet nonhomogeneus  \f$ p= p_0 \f$    (p=p_0)
// Neuman
    vel_fix  =10,///< 10= Neuman homogeneus or simmetry \f$ \nabla p \cdot \widehat{n}= 0 \f$ (dp.n=0)
    interiorp=11 ///< 11= Neuman nonhomogeneus \f$ \nabla p \cdot \widehat{n}=\Delta p_0 \f$  (dp.n=dp0)

};


#endif



#endif   //    __userP_h__  <-----------------------------------------------------------------------------------------
