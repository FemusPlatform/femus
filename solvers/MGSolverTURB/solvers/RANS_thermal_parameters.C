#include "Equations_conf.h"

#ifdef RANS_THERMAL_EQUATIONS

#include "RANS_thermal_parameters.h"
#include <stdio.h>
#include <fstream>
#include "MGUtils.h"

RANS_thermal_param::RANS_thermal_param() {
  _SolverType = BICGSTABM;  // BICGSTABM  GMRESM
  _SolveSteady = 0;
  _Supg = 1;
  _ModifiedSupg = "supgc";
  _Scpg = 0;
  _Upwind = 0;
  _ReactionNumberBased = 0;
  _DynamicUnderRelaxation = 0.;
  _MaxNonLinearIterations = 1;
  _FlatProfile = 1;
  _InterpolatedAlphaTurb = 0;
  _WallFunctionApproach = 0;
  _FractionalStep = 0;
  _BoundMap["TKwall0"] = TKwall0;
  _BoundMap["TKwall"] = TKwall;
  _BoundMap["TKinlet"] = TKinlet;
  _BoundMap["TKinsulation"] = TKinsulation;
  _BoundMap["TKheat_flux"] = TKheat_flux;
  _BoundMap["TKrobinT"] = TKrobinT;
  _BoundMap["TKLext_conv"] = TKLext_conv;
  _BoundMap["TKsimmetry"] = TKsimmetry;
  _BoundMap["TKWallFuncGrad"] = TKWallFuncGrad;
}

void RANS_thermal_param::read_param(MGUtils& mgutils) {
  read_file();
  print_par();

  // BOUNDARY CONDITION BLOCK ============================================================
  std::cout << "Thermal Turbulence boundary condition block \n";
  std::string GroupString = mgutils._sim_config["BoundaryGroups"];
  std::string temps;
  int Length = GroupString.length();
  int count, pos1;
  std::string BDcondK, BDcondW;
  count = pos1 = 0;

  while (count < Length) {
    if (GroupString.at(count) == ',') {
      temps = GroupString.substr(pos1, count - pos1);
      _BoundaryGroupsIDs.push_back(stoi(temps));
      pos1 = count + 1;
    }

    count++;
  }

  _BoundaryGroupsIDs.push_back(stod(GroupString.substr(pos1, Length - pos1)));
  const int numero = _BoundaryGroupsIDs.size();

  for (int i = 0; i < _BoundaryGroupsIDs.size(); i++) {
    BDcondK = "TTgroup" + std::to_string(_BoundaryGroupsIDs[i]) + "k";
    _map_TTKgroup[_BoundaryGroupsIDs[i]] = _BoundMap[_FileMap[BDcondK]];
    BDcondW = "TTgroup" + std::to_string(_BoundaryGroupsIDs[i]) + "w";
    _map_TTWgroup[_BoundaryGroupsIDs[i]] = _BoundMap[_FileMap[BDcondW]];
  }

  // END BOUNDARY CONDITION BLOCK ========================================================

  // SETTING OTHER PARAMETERS ------------------------------------------------------------
  std::cout << "Thermal Turbulence parameters block \n";
  _SolveSteady = stoi(_FileMap["SolveSteady"]);
  _Supg = stoi(_FileMap["Supg"]);
  _Upwind = stoi(_FileMap["Upwind"]);
  _ReactionNumberBased = stoi(_FileMap["ReactionNumberBased"]);
  _DynamicUnderRelaxation = stod(_FileMap["DynamicUnderRelaxation"]);
  _FlatProfile = stoi(_FileMap["FlatProfile"]);
  _ModifiedSupg = _FileMap["ModifiedSupg"];
  _Scpg = stoi(_FileMap["Scpg"]);
  _MaxNonLinearIterations = stoi(_FileMap["MaxNonLinearIterations"]);
  _InterpolatedAlphaTurb = stoi(_FileMap["InterpolatedAlphaTurb"]);
  _WallFunctionApproach = stoi(_FileMap["WallFunctionApproach"]);

  if (_FileMap["FractionalStep"] != "") { _FractionalStep = stoi(_FileMap["FractionalStep"]); }

  _FileMap.clear();
  return;
}

void RANS_thermal_param::read_file()  // READING Tparameter.in ======================================
{
  std::string string_value;
  std::string buf = "";  // read double, string, dummay
  std::ostringstream file, file2;
  file << getenv("APP_PATH") << "/DATA/RANS_thermal_properties.in";
  std::ifstream fin;
  fin.open(file.str().c_str());  // stream file

#ifdef PRINT_INFO

  if (fin.is_open()) { std::cerr << "\nInit Reading = " << file.str() << std::endl; }

#endif

  if (fin.is_open()) {
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
    std::cerr << "TK_param.read_file(): no parameter file found" << std::endl;
    abort();
  }

  std::cerr << "End Reading file " << file.str() << std::endl;
  fin.close();

  return;
}  // END READING ==========================================================================

void RANS_thermal_param::print_par() {
#ifdef PRINT_INFO
  std::cerr << "\033[038;5;" << KTT_F + 50
            << ";1m  \n----------------------------------------------\n THERMAL TURBULENCE MAP \n\033[0m"
            << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
    std::cerr << it->first << " " << it->second << "\n";
  }

  std::cerr << "\033[038;5;" << KTT_F + 50 << ";1m \n----------------------------------------------\n \033[0m"
            << std::endl;
#endif

  return;
}

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
