// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include <fstream>
#include "MGUtils.h"
#include "Pparameters.h"
#include "UserP.h"

P_param::P_param() {
  _SolverType = GMRESM;
  _NumRestartSol = 1;
  _BoundMap["interior"] = interiorp;
  _BoundMap["nostress"] = interiorp;
  _BoundMap["outflow"] = interiorp;
  _BoundMap["pressure_outlet"] = outflowp;
  _BoundMap["outflow_p"] = outflowp;
  _BoundMap["pressure_inlet"] = outflowp0;
  _BoundMap["slip"] = vel_fix;
  _BoundMap["wall"] = vel_fix;
  _BoundMap["penalty_turb"] = vel_fix;
  _BoundMap["velocity"] = vel_fix;
  _BoundMap["velocity_norm"] = vel_fix;
  _BoundMap["velocity_tang"] = vel_fix;
  _BoundMap["accelerating_swirl"] = vel_fix;
  _BoundMap["decelerating_swirl"] = vel_fix;
  _BoundMap["accelerating_stress"] = vel_fix;
  _BoundMap["decelerating_stress"] = vel_fix;
  _BoundMap["stress"] = vel_fix;
  _BoundMap["swirl"] = vel_fix;
}

// ============================================================
/// This function reads NS parameter file (NSparameter.in)
void P_param::read_file() {  //  ===========================================================

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

void P_param::read_param(MGUtils& mgutils, int Proc) {
  //  Reading parameter  -> NSproperties.in ------------------------------------------------
  read_file();

  // Boundary condition block ------------------------------------------------------------
  std::cout << "Navier-Stokes - Pressure Equation boundary condition block \n";
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

  if (Proc == 0) {  // PRINTING ON MESSAGE.LOG
    std::cerr << "=============================================================================\n";
    std::cerr << " NAVIER-STOKES - PRESSURE EQUATION - BOUNDARY CONDITIONS \n";
    for (auto& x : _map_NSgroup) {
      std::cerr << "Group n. " << x.first << "  Pressure condition " << x.second << std::endl;
    }
  }

  if (_FileMap["TimeDisc"] != "") {
    _TimeDisc = stoi(_FileMap["TimeDisc"]);
  } else {
    _TimeDisc = 1;
  }
  if (_FileMap["AssembleOnce"] != "") {
    _AssembleOnce = stoi(_FileMap["AssembleOnce"]);
  } else {
    _AssembleOnce = 0;
  }

  if (_FileMap["NodeIDrefPressure"] != "") {
    _NodeIDrefPressure = stoi(_FileMap["NodeIDrefPressure"]);
  } else {
    _NodeIDrefPressure = 0;
  }
  if (_FileMap["NumRestartSol"] != "") { _NumRestartSol = stoi(_FileMap["NumRestartSol"]); }

  _FileMap.clear();
  return;
}

#endif
