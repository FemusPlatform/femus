#ifndef __Tparameters_h__
#define __Tparameters_h__

#include "Equations_conf.h"
// ===================================
#ifdef T_EQUATIONS
// ==================================
#include <string>
#include <map>
#include <vector>
#include "Solvertype_enum.h"
#include "User_T.h"

using namespace std;

class MGUtils;

/*!  \defgroup T_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup T_param
// ============================================================================
class T_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int _SolveSteady;                    /// Flag for steady state solution - value: 0 or 1
    int _Supg;                           /// Flag for SUPG stabilization of advection term - value: 0 or 1
    double   _Upwind;                    /// Flag for normal upwind stabilization - value: 0 <Up <1
    bool _Solve;                         /// Flag for solution of temperature equation
    double _Prt;
    int _InterpolatedAlphaTurb;
    int _NumRestartSol;
    int _TimeDer;
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_condT>  _map_Tgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_condT> _BoundMap;  /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h


public:
    // constructor --------------------------------------------------------------
    T_param(); 

    ~T_param(); 

    void  read_param ( MGUtils &mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

#endif
#endif
