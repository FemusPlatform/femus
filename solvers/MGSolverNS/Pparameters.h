#ifndef __Pparameters_h__
#define __Pparameters_h__
// ============================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
#include "Solvertype_enum.h"
#include "UserP.h"
#include <map>
#include <vector>

class MGUtils;

class P_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_cond_p>  _map_NSgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_cond_p> _BoundMap;   /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    int _TimeDisc;
    int _AssembleOnce;
    int _NodeIDrefPressure;
    int _NumRestartSol;
    
    P_param();
   ~P_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_NSgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils, int proc=0 ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};


#endif
#endif   //    __userP_h__  <-----------------------------------------------------------------------------------------
