#ifndef __NSparameters_h__
#define __NSparameters_h__

#include "Equations_conf.h"
// ===================================
#ifdef NS_EQUATIONS
// ==================================
#include "Solvertype_enum.h"
#include "UserNS.h"
#include <vector>
#include <map>

class MGUtils;
enum bound_cond;

/*! \defgroup NS_param   Class Table:  Navier-Stokes equation parameters (NS_param) */
/// \ingroup NS_param
// ============================================================================
class NS_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int _SolveSteady;                    /// Flag for steady state solution - value: 0 or 1
    int _Supg;                           /// Flag for SUPG stabilization of advection term - value: 0 or 1
    int _Upwind;                         /// Flag for normal upwind stabilization - value: 0 or 1
    double _UnderRelaxation;             /// Flag for stabilization based on UnderRelaxation - value between 0 and 1
    int _ReactionNumberBased;            /// Flag for stabilization based on reaction number - value: 0 or 1
    int _FlatProfile;                    /// Flag for initial solution: flat profile or modulated as function of wall distance - value: 0 or 1
    bool _SolveNS;                         /// Flag for solution of temperature equation
    int _WallFunctionApproach;
    int _InterpolatedMuTurb;
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_cond>  _map_NSgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_cond> _BoundMap;   /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    int _Les;                                      /// Flag for turbulence with LES
    int _MaxNonLinearIterations;                   /// Maximum number of non linear interations
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    double _Penalty_n;
    double    _Penalty_tg;
    double    _Tg1_stress;
    double    _Tg2_stress;
    double   _Threshold;
    int      _TimeDisc;
    int      _NumRestartSol;
    
    NS_param();

    ~NS_param()
    {
        _BoundMap.clear();
        _FileMap.clear();
        _map_NSgroup.clear();
        _BoundaryGroupsIDs.clear();
    };

    void  read_param ( MGUtils & mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

#endif
#endif
