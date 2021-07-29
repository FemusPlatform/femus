#ifndef __userNS1c_h__
#define __userNS1c_h__

#include "Equations_conf.h"
// ===================================
#ifdef NS_EQUATIONS
// ==================================


/*! \defgroup NS_Boundary_conditions     Enum Table: NS Boundary conditions  */
/// \ingroup NS_Boundary_conditions
// ========================================================================
/// Navier-Stokes boundary conditions
enum bound_cond1c {
    interior1c               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$
    nostress1c               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$
    outflow1c                =       11,     ///<  \f$  p=0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero pressure and stress along tangential direction
    pressure_outlet1c        =       28,     ///<  \f$  p=0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Zero pressure and velocity field along tangential direction
    outflow_p1c              =       31,     ///<  \f$  p=p_0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Fixed pressure and null stress along tangential direction
    pressure_inlet1c         =       38,     ///<  \f$  p=p_0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed pressure and zero velocity field along tangential direction
    penalty_turb1c           =       44,
    accelerating_swirl1c     =       85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_swirl1c     =      -85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    swirl1c                  =       85,
    accelerating_stress1c    =       84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_stress1c    =      -84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    stress1c                 =       84,
    slip1c                   =       81,     ///<  \f$ {\bf u} \cdot {\bf n} = 0 \wedge ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero velocity field along normal direction and free stress along tangential directions
    wall1c                  =       88,     ///<  \f$ {\bf u} = 0 \f$: Zero velocity field, noslip condition
    velocity1c               =       99,     ///<  \f$ {u_n} = {u_n^*} \wedge {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Fixed velocity components along normal and tangential directions
    velocity_norm1c          =       98,     ///<  \f$  {u_n} = {u_n^*}  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed velocity along normal direction and zero velocity components along tangential direction
    velocity_tang1c          =       89      ///<  \f$  {u_n} = 0  \wedge  {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Zero velocity field along normal direction and fixed velocity components along tangential direction
};

#endif
#endif
