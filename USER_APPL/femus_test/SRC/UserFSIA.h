
#ifndef __userFSIA_h__
#define __userFSIA_h__

#include "Equations_conf.h"
// ===================================
#ifdef FSIA_EQUATIONS
// ==================================
#include "Solvertype_enum.h"

/*! \defgroup FSIA_Boundary_conditions     Enum Table: NS Boundary conditions  */
/// \ingroup FSIA_Boundary_conditions
// ========================================================================
/// Adjoint Navier-Stokes boundary conditions
enum bound_cond_fsia {
    interior_fsia               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$  
    nostress_fsia               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$  
    outflow_fsia               =       11,     ///<  \f$  p=0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero pressure and stress along tangential direction
    pressure_outlet_fsia        =       28,     ///<  \f$  p=0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Zero pressure and velocity field along tangential direction
    outflow_p_fsia              =       31,     ///<  \f$  p=p_0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Fixed pressure and null stress along tangential direction
    pressure_inlet_fsia         =       38,     ///<  \f$  p=p_0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed pressure and zero velocity field along tangential direction
    penalty_turb_fsia           =       44,
//     accelerating_swirl_fsia     =       85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
//     decelerating_swirl_fsia     =      -85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    swirl_fsia                  =      85,
//     accelerating_stress_fsia    =       84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
//     decelerating_stress_fsia    =      -84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    stress_fsia               =       84,
    slip_fsia                   =       81,     ///<  \f$ {\bf u} \cdot {\bf n} = 0 \wedge ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero velocity field along normal direction and free stress along tangential directions
    wall_fsia                   =       88,     ///<  \f$ {\bf u} = 0 \f$: Zero velocity field, noslip condition
    velocity_fsia               =       99,     ///<  \f$ {u_n} = {u_n^*} \wedge {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Fixed velocity components along normal and tangential directions
    velocity_norm_fsia          =       98,     ///<  \f$  {u_n} = {u_n^*}  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed velocity along normal direction and zero velocity components along tangential direction
    velocity_tang_fsia          =       89      ///<  \f$  {u_n} = 0  \wedge  {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Zero velocity field along normal direction and fixed velocity components along tangential direction
};

/*! \defgroup FSIA_param   Class Table:  Navier-Stokes equation parameters (NS_param) */
/// \ingroup FSIA_param
// ============================================================================
class FSIA_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int _SolveSteady;                    /// Flag for steady state solution - value: 0 or 1
    int _Supg;                           /// Flag for SUPG stabilization of advection term - value: 0 or 1
    int _Upwind;                         /// Flag for normal upwind stabilization - value: 0 or 1
    double _UnderRelaxation;             /// Flag for stabilization based on UnderRelaxation - value between 0 and 1
    int _ReactionNumberBased;            /// Flag for stabilization based on reaction number - value: 0 or 1
    int _FlatProfile;                    /// Flag for initial solution: flat profile or modulated as function of wall distance - value: 0 or 1
//     bool _SolveFSI;                         /// Flag for solution of temperature equation
    int _WallFunctionApproach;
    int _InterpolatedMuTurb;
    int _Compressible;
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_cond_fsia>  _map_FSIAgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_cond_fsia> _BoundMap;   /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    int _Les;                                      /// Flag for turbulence with LES
    int _MaxNonLinearIterations;                   /// Maximum number of non linear interations
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    double _Penalty_n;
    double    _Penalty_tg;
    double    _Tg1_stress;
    double    _Tg2_stress;
    FSIA_param() {
        _SolverType              =GMRESM;
        _SolveSteady             =0;
        _Supg                    =1;
        _Upwind                  =0.;
        _ReactionNumberBased     =0;
        _FlatProfile             =1;
        _UnderRelaxation         =0.;
        _Les                     =0.;
        _MaxNonLinearIterations  =1;
        _WallFunctionApproach    =0;
	_InterpolatedMuTurb      =0;
        _Penalty_n               =1e+5;
        _Penalty_tg              =0.;
	_Tg1_stress              =0.;
        _Tg2_stress              =0.;
        _Compressible            =1;
        _BoundMap["interior_fsia"]        = interior_fsia;
        _BoundMap["nostress_fsia"]        = nostress_fsia;
        _BoundMap["outflow_fsia"]         = outflow_fsia;
        _BoundMap["pressure_outlet_fsia"] = pressure_outlet_fsia;
        _BoundMap["outflow_p_fsia"]       = outflow_p_fsia;
        _BoundMap["pressure_inlet_fsia"]  = pressure_inlet_fsia;
        _BoundMap["slip_fsia"]            = slip_fsia;
        _BoundMap["wall_fsia"]            = wall_fsia;
        _BoundMap["penalty_turb_fsia"]    = penalty_turb_fsia;
        _BoundMap["velocity_fsia"]        = velocity_fsia;
        _BoundMap["velocity_norm_fsia"]   = velocity_norm_fsia;
        _BoundMap["velocity_tang_fsia"]   = velocity_tang_fsia; 
//         _BoundMap["accelerating_swirl_fsia"]           = accelerating_swirl_fsia;
//         _BoundMap["decelerating_swirl_fsia"]           = decelerating_swirl_fsia;
//         _BoundMap["accelerating_stress_fsia"]    = accelerating_stress_fsia;
//         _BoundMap["decelerating_stress_fsia"]    = decelerating_stress_fsia;
        _BoundMap["stress_fsia"]    = stress_fsia;
        _BoundMap["swirl_fsia"]    = swirl_fsia;
    }
    
   ~FSIA_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_FSIAgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

#endif
#endif
