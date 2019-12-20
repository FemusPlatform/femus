#ifndef __userNS_h__
#define __userNS_h__

#include "Equations_conf.h"
// ===================================
#ifdef NS_EQUATIONS
// ==================================


/*! \defgroup NS_Boundary_conditions     Enum Table: NS Boundary conditions  */
/// \ingroup NS_Boundary_conditions
// ========================================================================
/// Navier-Stokes boundary conditions
enum bound_cond {
    interior               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$
    nostress               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$
    outflow                =       11,     ///<  \f$  p=0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero pressure and stress along tangential direction
    pressure_outlet        =       28,     ///<  \f$  p=0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Zero pressure and velocity field along tangential direction
    outflow_p              =       31,     ///<  \f$  p=p_0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Fixed pressure and null stress along tangential direction
    pressure_inlet         =       38,     ///<  \f$  p=p_0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed pressure and zero velocity field along tangential direction
    penalty_turb           =       0,
    periodic_stress        =       44,
    accelerating_swirl     =       85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_swirl     =      -85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    swirl                  =       85,
    accelerating_stress    =       84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_stress    =      -84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    stress                 =       84,
    slip                   =       81,     ///<  \f$ {\bf u} \cdot {\bf n} = 0 \wedge ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero velocity field along normal direction and free stress along tangential directions
    wall                   =       88,     ///<  \f$ {\bf u} = 0 \f$: Zero velocity field, noslip condition
    velocity               =       99,     ///<  \f$ {u_n} = {u_n^*} \wedge {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Fixed velocity components along normal and tangential directions
    velocity_norm          =       98,     ///<  \f$  {u_n} = {u_n^*}  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed velocity along normal direction and zero velocity components along tangential direction
    velocity_tang          =       89      ///<  \f$  {u_n} = 0  \wedge  {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Zero velocity field along normal direction and fixed velocity components along tangential direction
};

#endif
#endif
