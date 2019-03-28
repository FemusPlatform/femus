#ifndef __RANS_parameters_h__
#define __RANS_parameters_h__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_EQUATIONS
// ==================================
#include <string>
#include <map>
#include <vector>
#include "Solvertype_enum.h"

class MGUtils;

class RANS_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int _SolveSteady;
    int _Supg;
    std::string _ModifiedSupg;
    int _Scpg;
    double _Upwind;
    int _ReactionNumberBased;
    double _DynamicUnderRelaxation;
    int _MaxNonLinearIterations;
    int _FlatProfile;
    bool _Solve;
    int _InterpolatedMuTurb;
    int _WallFunctionApproach;
    int _FractionalStep;
    int _TimeDer;
    std::vector < int >_BoundaryGroupsIDs;
    std::map < int, bound_condK > _map_DTKgroup;
    std::map < int, bound_condK > _map_DTWgroup;
    std::map < std::string, bound_condK > _BoundMap;
    std::map < std::string, std::string > _FileMap;
    SolverTypeM _SolverType;	// for other solver types see Solvertype_enum.h
public:
    // constructor --------------------------------------------------------------
    RANS_param ();

    ~RANS_param ();

    void read_param ( MGUtils & mgutils );	/// This function sets all the T_param parameters
    void read_file ();		/// This function reads the parameters from Tproperties.in file
    void print_par ();		/// This function prints the parameters contained in Tproperties.in
};





#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
