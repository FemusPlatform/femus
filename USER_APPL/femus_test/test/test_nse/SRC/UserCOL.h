
#ifndef __userCOL_h__
#define __userCOL_h__

#include "Equations_conf.h"
// ===================================
#ifdef COLOR_EQUATIONS
// ==================================

/*!  \defgroup COL_param   Class Table:  energy equation parameters (COL_param) */    
/// \ingroup COL_param
// ============================================================================
class COL_param
{//< This class defines the physical and numerical  energy equation parameters
public:
     int AXISYM  /**< axisymmetry  */;
    // stabilization NS -------------------------------------------------------
    double SUPG/**< supg */; double UPWIND;/**< Normal Upwind */ double UPWIND2;/**< Transv Upwind */
    // turbulence
    double LES/**< les */;  double PR_T/**< Prandtl */; double DIST_FIX/**< distance from wall */;
    // non linear -------------------------------------------------------------
    int    NL_ITER;/**< NON LIN IT */    int    NL_ITER0;/**< INIT NON LIN IT */
    double NL_TIME0; // initial time for  non linear regime
    // bc  -------------------------------------------------------------
    double H_EXT;/**< h convective in bc  */ double T_EXT;/**< T_ext convective in bc  */ 
    // time discretization
    double CRANK_NICK/**< implicit (1) explicit (0) Crank-Nicolson (0.5)*/; 
    int UNSTEADY;    /**< un(1)/steady(0) flag  */
public:    
    // constructor --------------------------------------------------------------
    COL_param() {
        
        AXISYM=0;
        // stabilization NS -----------------------------------------------------
        SUPG=1.0/**< SUPG */; UPWIND=.0/**< Normal Upwind */; UPWIND2=.0/**< Transv Upwind */;
        // turbulence
        LES=.0/**< LES */; PR_T=0.9/**< turb Prandtl */;DIST_FIX=1.e-3/**< distance from wall */;
        // non linear ------------------------------------------------------------
        NL_ITER=0;/**< NON LIN IT */ NL_ITER0=0;/**< INIT NON LIN IT */
        NL_TIME0=-.0001;
        // bc  -------------------------------------------------------------
        H_EXT=0.;/**< h convective in bc  */ T_EXT=0.;/**< T_ext convective in bc  */ 
        // time discretization
        CRANK_NICK=1.;  /**< implicit (1) explicit (0) Crank-Nicolson (0.5) */
        UNSTEADY=1;     /**< un(1)/steady(0) flag  */
    }
    inline void set_AXISYM(int val){AXISYM=val;}
    inline void set_PRT(double val)  {PR_T=val;}
    inline void set_HEXT(double val) {H_EXT=val;}
    inline void set_TEXT(double val) {T_EXT=val;}
    inline void set_UNSTEADY(double val) {UNSTEADY=val;}
};





// T boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */ 
/// \ingroup Boundary_conditions
enum bound_condC {
    // ========================================================================
    Twall0c=     0, ///< 0= Dirichlet homogeneus     (T=0)
    Twallc=      1, ///< 1= Dirichlet nonhomogeneus  (T=T_0)

    insulationc=10, ///< 10= Neuman homogeneus         (dT.n=0)
    heat_fluxc= 11, ///< 11= Neuman nonhomogeneus      (dT.n=q_0)
    robinTc=    12, ///< 12= Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
    ext_convc=  13, ///< 13= Turbulence      (dT.n=tau_0^2=beta*ug(u))
    
    simmetryc=  20  ///< 20= Simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
               // ========================================================================
       
};
#endif
#endif