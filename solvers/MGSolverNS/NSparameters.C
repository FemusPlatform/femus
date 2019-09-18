// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

#include <fstream>
#include "MGUtils.h"
#include "NSparameters.h"  // Navier-Stokes class header file

// ======================================================================
/// Boundary and parameter setting for NS equation from NSproperties.in
///   BOUNDARY CONDITIONS  in the file -> NSproperties.in :
///   In order to set the boundary condition follow the example below
///  NSgroup<Group Number>    <Boundary Condition for kh>,<Boundary Condition for wh>
///         NSgroup10       <Bound_cond>
///   ! A condition for each group defined in SimulationConfiguration.in file must be specified !
///   NS PARAMETERS  in the file -> NSproperties.in :
///   SolveSteady  1  If 1 (yes) then steady state equation is solved - without  time derivative
///   MaxNonLinearIterations    0  -> non linear
///   Supg                        1 Supg: standard supg formulation
///   Upwind                      0      -> (yes) or 0 (no)
///   Les                         0       ->  1(yes) or 0(no)
///   ReactionNumberBased         0      ->    (yes) or 0 (no)
///   DynamicUnderRelaxation      0.     ->    0 < UnderRelaxation < 1
///   FlatProfile                 1      -> on inlet
//  =====================================================================

NS_param::NS_param() {
  _SolverType = GMRESM;
  _SolveSteady = 0;
  _Supg = 1;
  _Upwind = 0.;
  _ReactionNumberBased = 0;
  _FlatProfile = 1;
  _UnderRelaxation = 0.;
  _Les = 0.;
  _MaxNonLinearIterations = 1;
  _WallFunctionApproach = 0;
  _InterpolatedMuTurb = 0;
  _Penalty_n = 1e+5;
  _Penalty_tg = 0.;
  _Tg1_stress = 0.;
  _Tg2_stress = 0.;
  _Threshold = 0.495;
  _TimeDisc = 1;
  _NumRestartSol = 1;
  _BoundMap["interior"] = interior;
  _BoundMap["nostress"] = nostress;
  _BoundMap["outflow"] = outflow;
  _BoundMap["pressure_outlet"] = pressure_outlet;
  _BoundMap["outflow_p"] = outflow_p;
  _BoundMap["pressure_inlet"] = pressure_inlet;
  _BoundMap["slip"] = slip;
  _BoundMap["wall"] = wall;
  _BoundMap["penalty_turb"] = penalty_turb;
  _BoundMap["velocity"] = velocity;
  _BoundMap["velocity_norm"] = velocity_norm;
  _BoundMap["velocity_tang"] = velocity_tang;
  _BoundMap["accelerating_swirl"] = accelerating_swirl;
  _BoundMap["decelerating_swirl"] = decelerating_swirl;
  _BoundMap["accelerating_stress"] = accelerating_stress;
  _BoundMap["decelerating_stress"] = decelerating_stress;
  _BoundMap["stress"] = stress;
  _BoundMap["swirl"] = swirl;
}

void NS_param::read_param(MGUtils& mgutils) {
  //  Reading parameter  -> NSproperties.in ------------------------------------------------
  read_file();

  // Boundary condition block ------------------------------------------------------------
  std::cout << "Navier-Stokes boundary condition block \n";
  // boundary group names  from
  std::string GroupString = mgutils._sim_config["BoundaryGroups"];
  int count = 0;
  int pos1 = 0;
  int Length = GroupString.length();

  while (count < Length) {
    if (GroupString.at(count) == ',') {
      std::string temps = GroupString.substr(pos1, count - pos1);
      _BoundaryGroupsIDs.push_back(stoi(temps));
      pos1 = count + 1;
    }

    count++;
  }

  // boundary group values from
  _BoundaryGroupsIDs.push_back(stod(GroupString.substr(pos1, Length - pos1)));

  for (int i = 0; i < _BoundaryGroupsIDs.size(); i++) {
    std::string BDcond = "NSgroup" + std::to_string(_BoundaryGroupsIDs[i]);
    _map_NSgroup[_BoundaryGroupsIDs[i]] = _BoundMap[_FileMap[BDcond]];
  }

  // ---------------------  end bc ------------------------------------
  // Navier-Stokes parameters block ------------------------------------------------------------
  std::cout << " NAVIER-STOKES PARAMETER map (_NS_parameter) in  NSproperties.in + UserNS.h :\n";

  if (_FileMap["SolveSteady"] != "") {
    _SolveSteady = stoi(_FileMap["SolveSteady"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._SolveSteady (in UserNS.h) \n";
  }

  if (_FileMap["MaxNonLinearIterations"] != "") {
    _MaxNonLinearIterations = stoi(_FileMap["MaxNonLinearIterations"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._MaxNonLinearIterations (in UserNS.h)\n";
  }

  if (_FileMap["DynamicUnderRelaxation"] != "") {
    _UnderRelaxation = stof(_FileMap["DynamicUnderRelaxation"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._DynamicUnderRelaxation (in UserNS.h)\n";
  }

  std::cout << " NSproperties.in: default value for _NS_parameter._SolverType  (in UserNS.h) \n";

  if (_FileMap["SolverType"] != "") { _SolverType = mgutils._SolverTypeMap[_FileMap["SolverType"]]; }

  if (_FileMap["Supg"] != "") {
    _Supg = stoi(_FileMap["Supg"]);
  } else {
    std::cout << " NSproperties.in:  default value for _NS_parameter.Supg (in UserNS.h)\n";
  }

  if (_FileMap["Upwind"] != "") {
    _Upwind = stof(_FileMap["Upwind"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._Upwind (in UserNS.h)\n";
  }

  if (_FileMap["Les"] != "") {
    _Les = stof(_FileMap["Les"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._Les (in UserNS.h)\n";
  }

  if (_FileMap["ReactionNumberBased"] != "") {
    _ReactionNumberBased = stoi(_FileMap["ReactionNumberBased"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._ReactionNumberBased  (in UserNS.h)\n";
  }

  if (_FileMap["FlatProfile"] != "") {
    _FlatProfile = stoi(_FileMap["FlatProfile"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._FlatProfile (in UserNS.h)\n";
  }

  if (_FileMap["InterpolatedMuTurb"] != "") {
    _InterpolatedMuTurb = stoi(_FileMap["InterpolatedMuTurb"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._InterpolatedMuTurb (in UserNS.h)\n";
  }

  if (_FileMap["WallFunctionApproach"] != "") {
    _WallFunctionApproach = stoi(_FileMap["WallFunctionApproach"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._WallFunctionApproach (in UserNS.h)\n";
  }

  if (_FileMap["Penalty_n"] != "") {
    _Penalty_n = stof(_FileMap["Penalty_n"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter.Penalty_n (in UserNS.h)\n";
  }

  if (_FileMap["Penalty_tg"] != "") {
    _Penalty_tg = stof(_FileMap["Penalty_tg"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter.Penalty_tg (in UserNS.h)\n";
  }

  if (_FileMap["Tg1_stress"] != "") {
    _Tg1_stress = stof(_FileMap["Tg1_stress"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter.Tg1_stress (in UserNS.h)\n";
  }

  if (_FileMap["Tg2_stress"] != "") {
    _Tg2_stress = stof(_FileMap["Tg2_stress"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter.Tg2_stress (in UserNS.h)\n";
  }

  if (_FileMap["Threshold"] != "") {
    _Threshold = stof(_FileMap["Threshold"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter.Threshold (in UserNS.h)\n";
  }
  if (_FileMap["TimeDisc"] != "") {
    _TimeDisc = stoi(_FileMap["TimeDisc"]);
  } else {
    std::cout << " NSproperties.in: default value for _NS_parameter._TimeDisc (in UserNS.h)\n";
  }
  if (_FileMap["NumRestartSol"] != "") { _NumRestartSol = stoi(_FileMap["NumRestartSol"]); }
  // ---------------------  end parameter NS ------------------------------------
  _FileMap.clear();
  return;
}

// ============================================================
/// This function reads NS parameter file (NSparameter.in)
void NS_param::read_file() {  //  ===========================================================

  //  getting file name --------------------------------------------------------------------
  std::ostringstream file;
  file << getenv("APP_PATH") << "/DATA/NSproperties.in";
  std::ifstream fin;
  fin.open(file.str().c_str());  // stream file
#ifdef PRINT_INFO

  if (fin.is_open()) { std::cout << "\nInit Reading = " << file.str() << std::endl; }

#endif

  //  reading param file -----------------------------------------------------------------
  if (fin.is_open()) {
    std::string string_value;
    std::string buf = "";  // read double, string, dummy

    while (buf != "/") {
      fin >> buf;  // find "/" file start
    }

    fin >> buf;

    while (buf != "/") {
      if (buf == "#") {
        getline(fin, buf);  // comment line
      } else {
        fin >> string_value;
        _FileMap[buf] = string_value;
      }

      fin >> buf;
    }
  } else {
    std::cerr << "NS_param.read_file(): no parameter file found" << std::endl;
    abort();
  }

// printing after reading
// ----------------------------------------------------------------------------------------
#ifdef PRINT_INFO
  std::cout << "\033[038;5;"
            << "\n  NAVIER-STOKES PARAMETER MAP \n \033[0m" << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }

  std::cout << "\033[038;5;" << NS_F + 50 << ";1m \
                \n----------------------------------------------\n\033[0m" << std::endl;
#endif

  fin.close();
  return;
}

#endif  // ENDIF NS_EQUATIONS
