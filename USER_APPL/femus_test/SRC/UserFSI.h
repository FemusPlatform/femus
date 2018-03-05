
#ifndef __userFSI_h__
#define __userFSI_h__

#include "Equations_conf.h"
// ===================================
#ifdef FSI_EQUATIONS
// ==================================
#include "Solvertype_enum.h"

/*! \defgroup FSI_Boundary_conditions     Enum Table: NS Boundary conditions  */
/// \ingroup FSI_Boundary_conditions
// ========================================================================
/// Navier-Stokes boundary conditions
enum bound_cond {
    interior               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$  
    nostress               =       11,     ///< 1 (n) \f$ {\sigma}_{nn}=0  \f$  + 1(t)   \f$   {\sigma}_{nt}=0   \f$  
    outflow                =       11,     ///<  \f$  p=0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero pressure and stress along tangential direction
    pressure_outlet        =       28,     ///<  \f$  p=0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Zero pressure and velocity field along tangential direction
    outflow_p              =       31,     ///<  \f$  p=p_0  \wedge  ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Fixed pressure and null stress along tangential direction
    pressure_inlet         =       38,     ///<  \f$  p=p_0  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed pressure and zero velocity field along tangential direction
    penalty_turb           =       44,
    accelerating_swirl     =       85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_swirl     =      -85,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t2} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t2}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    swirl                  =      85,
    accelerating_stress    =       84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Accelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    decelerating_stress    =      -84,     ///<  \f$ ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_{t1} } \propto {\bf u_{t1}} \f$: Decelerating stress for velocity component \f$ \bf{u_{t1}} \f$ velocity along \f$ \bf{t_1} \f$ (see #MGSolNS::TangType())
    stress                 =       84,
    slip                   =       81,     ///<  \f$ {\bf u} \cdot {\bf n} = 0 \wedge ( {\boldsymbol \tau} \cdot {\bf n} ) \cdot {\boldsymbol \varphi_t } = 0  \f$: Zero velocity field along normal direction and free stress along tangential directions
    wall                   =       88,     ///<  \f$ {\bf u} = 0 \f$: Zero velocity field, noslip condition
    velocity               =       99,     ///<  \f$ {u_n} = {u_n^*} \wedge {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Fixed velocity components along normal and tangential directions
    velocity_norm          =       98,     ///<  \f$  {u_n} = {u_n^*}  \wedge  u_{t1} = 0 \wedge u_{t2} = 0  \f$: Fixed velocity along normal direction and zero velocity components along tangential direction
    velocity_tang          =       89      ///<  \f$  {u_n} = 0  \wedge  {u_{t1}} = {u_{t1}^*} \wedge {u_{t2}} = {u_{t2}^*} \f$: Zero velocity field along normal direction and fixed velocity components along tangential direction
};

/*! \defgroup FSI_param   Class Table:  Navier-Stokes equation parameters (NS_param) */
/// \ingroup FSI_param
// ============================================================================
class FSI_param
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
    std::map<int,  bound_cond>  _map_FSIgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_cond> _BoundMap;   /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    int _Les;                                      /// Flag for turbulence with LES
    int _MaxNonLinearIterations;                   /// Maximum number of non linear interations
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    double _Penalty_n;
    double    _Penalty_tg;
    double    _Tg1_stress;
    double    _Tg2_stress;
    FSI_param() {
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
        _BoundMap["interior"]        = interior;
        _BoundMap["nostress"]        = nostress;
        _BoundMap["outflow"]         = outflow;
        _BoundMap["pressure_outlet"] = pressure_outlet;
        _BoundMap["outflow_p"]       = outflow_p;
        _BoundMap["pressure_inlet"]  = pressure_inlet;
        _BoundMap["slip"]            = slip;
        _BoundMap["wall"]            = wall;
        _BoundMap["penalty_turb"]    = penalty_turb;
        _BoundMap["velocity"]        = velocity;
        _BoundMap["velocity_norm"]   = velocity_norm;
        _BoundMap["velocity_tang"]   = velocity_tang; 
        _BoundMap["accelerating_swirl"]           = accelerating_swirl;
        _BoundMap["decelerating_swirl"]           = decelerating_swirl;
        _BoundMap["accelerating_stress"]    = accelerating_stress;
        _BoundMap["decelerating_stress"]    = decelerating_stress;
        _BoundMap["stress"]    = stress;
        _BoundMap["swirl"]    = swirl;
    }
    
   ~FSI_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_FSIgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

#endif
#endif
