#ifndef __RANS_thermal_parameters_h__
#define __RANS_thermal_parameters_h__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_THERMAL_EQUATIONS
// ==================================
#include <map>
#include <string>
#include <vector>
#include "Solvertype_enum.h"
#include "User_RANS_thermal.h"

class MGUtils;

/*!  \defgroup TK_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup TK_param
// ============================================================================

class RANS_thermal_param {
  //< This class defines the physical and numerical  energy equation parameters
 public:
  int _SolveSteady;
  int _Supg;
  std::string _ModifiedSupg;
  int _Scpg;
  int _Upwind;
  int _ReactionNumberBased;
  double _DynamicUnderRelaxation;
  int _MaxNonLinearIterations;
  int _FlatProfile;
  bool _Solve;
  int _InterpolatedAlphaTurb;
  int _WallFunctionApproach;
  int _FractionalStep;
  std::vector<int> _BoundaryGroupsIDs;
  std::map<int, bound_condRANS_T> _map_TTKgroup;
  std::map<int, bound_condRANS_T> _map_TTWgroup;
  std::map<std::string, bound_condRANS_T> _BoundMap;
  std::map<std::string, std::string> _FileMap;
  SolverTypeM _SolverType;  // for other solver types see Solvertype_enum.h
 public:
  // constructor --------------------------------------------------------------
  RANS_thermal_param();

  ~RANS_thermal_param() {
    _BoundMap.clear();
    _FileMap.clear();
    _map_TTKgroup.clear();
    _map_TTWgroup.clear();
    _BoundaryGroupsIDs.clear();
  };

  void read_param(MGUtils& mgutils);  /// This function sets all the T_param parameters
  void read_file();                   /// This function reads the parameters from Tproperties.in file
  void print_par();                   /// This function prints the parameters contained in Tproperties.in
};
#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
