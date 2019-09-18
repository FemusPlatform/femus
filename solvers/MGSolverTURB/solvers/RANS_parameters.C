#include "Equations_conf.h"

// ============================================
#ifdef RANS_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "User_RANS.h"
#include "RANS_parameters.h"
#include "MGUtils.h"
#include <fstream>

RANS_param::RANS_param() {
  _SolverType = GMRESM;
  _SolveSteady = 0;
  _Supg = 1;
  _ModifiedSupg = "supgc";
  _Scpg = 0;
  _Upwind = 0;
  _ReactionNumberBased = 0;
  _DynamicUnderRelaxation = 0.;
  _MaxNonLinearIterations = 1;
  _FlatProfile = 1;
  _InterpolatedMuTurb = 0;
  _WallFunctionApproach = 0;
  _FractionalStep = 0;
  _BoundMap["Kwall0"] = Kwall0;
  _BoundMap["Kwall"] = Kwall;
  _BoundMap["Kinlet"] = Kinlet;
  _BoundMap["Kinsulation"] = Kinsulation;
  _BoundMap["Kheat_flux"] = Kheat_flux;
  _BoundMap["KrobinT"] = KrobinT;
  _BoundMap["KLext_conv"] = KLext_conv;
  _BoundMap["Ksimmetry"] = Ksimmetry;
  _BoundMap["Kinit"] = Kinit;
  _BoundMap["KWallFuncGrad"] = KWallFuncGrad;
}

RANS_param::~RANS_param() {
  _BoundMap.clear();
  _FileMap.clear();
  _map_DTKgroup.clear();
  _map_DTWgroup.clear();
  _BoundaryGroupsIDs.clear();
};

void RANS_param::read_param(MGUtils& mgutils) {
  read_file();
  print_par();

  // BOUNDARY CONDITION BLOCK ============================================================
  std::cout << "Turbulence boundary condition block \n";
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
    BDcondK = "DTgroup" + std::to_string(_BoundaryGroupsIDs[i]) + "k";
    _map_DTKgroup[_BoundaryGroupsIDs[i]] = _BoundMap[_FileMap[BDcondK]];
    BDcondW = "DTgroup" + std::to_string(_BoundaryGroupsIDs[i]) + "w";
    _map_DTWgroup[_BoundaryGroupsIDs[i]] = _BoundMap[_FileMap[BDcondW]];
  }

  // END BOUNDARY CONDITION BLOCK ========================================================

  // SETTING OTHER PARAMETERS ------------------------------------------------------------
  std::cout << "Turbulence parameters block \n";
  _SolveSteady = stoi(_FileMap["SolveSteady"]);
  _Supg = stoi(_FileMap["Supg"]);
  _Upwind = stod(_FileMap["Upwind"]);
  _ReactionNumberBased = stoi(_FileMap["ReactionNumberBased"]);
  _DynamicUnderRelaxation = stod(_FileMap["DynamicUnderRelaxation"]);
  _FlatProfile = stoi(_FileMap["FlatProfile"]);
  _ModifiedSupg = _FileMap["ModifiedSupg"];
  _Scpg = stoi(_FileMap["Scpg"]);
  _MaxNonLinearIterations = stoi(_FileMap["MaxNonLinearIterations"]);
  _InterpolatedMuTurb = stoi(_FileMap["InterpolatedMuTurb"]);
  _WallFunctionApproach = stoi(_FileMap["WallFunctionApproach"]);

  if (_FileMap["FractionalStep"] != "") { _FractionalStep = stoi(_FileMap["FractionalStep"]); }

  if (mgutils._SolverTypeMap[_FileMap["SolverType"]] != 0) {
    _SolverType = mgutils._SolverTypeMap[_FileMap["SolverType"]];
  }

  if (_FileMap["TimeDer"] != "") {
    _TimeDer = stoi(_FileMap["TimeDer"]);
  } else {
    _TimeDer = 1;
  }

  std::cerr << "RANS SOLVER TYPE  " << _SolverType << std::endl;
  _FileMap.clear();
  return;
}

void RANS_param::read_file()  // READING Tparameter.in ======================================
{
  std::string string_value;
  std::string buf = "";  // read double, string, dummay
  std::ostringstream file, file2;
  file << getenv("APP_PATH") << "/DATA/RANS_dynamical_properties.in";
  std::ifstream fin;
  fin.open(file.str().c_str());  // stream file

#ifdef PRINT_INFO

  if (fin.is_open()) { std::cout << "\nInit Reading = " << file.str() << std::endl; }

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
    std::cerr << "RANS_param.read_file(): no parameter file found" << std::endl;
    abort();
  }

  std::cout << "End Reading file " << file.str() << std::endl;
  fin.close();

  return;
}  // END READING ==========================================================================

void RANS_param::print_par() {
#ifdef PRINT_INFO
  std::cout
      << "\033[038;5;217;1m \n----------------------------------------------\n  TURBULENCE MAP \n\033[0m"
      << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }

  std::cout << "\033[038;5;217;1m \n----------------------------------------------\n \033[0m" << std::endl;
#endif

  return;
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
