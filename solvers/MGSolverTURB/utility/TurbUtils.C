#include "TurbUtils.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Utils_DynTurbModels_4p.h"
#include "Utils_DynTurbModels_Wilcox.h"
#include "Utils_ThermTurbModels_4p.h"
//
using namespace std;

TurbUtils::TurbUtils(int MeshId, int nMeshes) : __Proc(0), __MeshID(MeshId), __NMeshes(nMeshes) {
  FillModelMap();
  FillParameters();  // DEFAULT VALUES
}

TurbUtils::~TurbUtils() {}

void TurbUtils::FillParameters() {
  _BoundWallDist = 1.e-3;
  _nu = 1.e-4;
  _alpha = 1.e-4;
  _klog = 0.;
  _wlog = 0.;
  _khlog = 0.;
  _whlog = 0.;

  read_file();

  _SolveAlphaT = _SolveMuT = 0;

  std::cout << "\n ========================================================= \n"
            << "\033[38;5;118m \t \t SETTING THE TURBULENCE MODEL \033[0m \n";

  _BoundWallDist = stod(_FileMap["Wall_dist"]);

  const double rho = stod(_FileMap["rho0"]);
  const double mu = stod(_FileMap["mu0"]);
  const double kappa = stod(_FileMap["kappa0"]);
  const double cp = stod(_FileMap["cp0"]);

  _nu = mu / rho;
  _alpha = kappa / (rho * cp);
  __IPr = _alpha / _nu;

  //
  std::cerr << "==================================================== \n";
  std::cerr << " TURB_UTILS: available dynamic turbulence models: \n";

  for (auto it = ::DynTurbModelMap.cbegin(); it != ::DynTurbModelMap.cend(); ++it) {
    std::cerr << "     " << it->first << " " << it->second << "\n";
  }

  std::cerr << " TURB_UTILS: available thermal turbulence models: \n";

  for (auto it = ::ThermTurbModelMap.cbegin(); it != ::ThermTurbModelMap.cend(); ++it) {
    std::cerr << "     " << it->first << " " << it->second << "\n";
  }

  std::cerr << "    Chosen dynamic turb model: " << _FileMap["RANS_dynamic"] << std::endl;
  std::cerr << "    Chosen themral turb model: " << _FileMap["RANS_thermal"] << std::endl;

  std::string DynModel = (_FileMap["RANS_dynamic"] != "") ? _FileMap["RANS_dynamic"] : "default";
  int RANS_dynamic = ::DynTurbModelMap.at(DynModel);

  switch (RANS_dynamic) {
    case default_dyn:
      _Nagano = 1;
      _Wilcox = _klog = _wlog = _emod = _numod = 0;
      _DynModel = new Dyn4P_KW();
      break;

    case nagano_ke:
      _Nagano = _emod = 1;
      _Wilcox = _klog = _wlog = _numod = 0;
      _DynModel = new Dyn4P_KE();
      break;

    case nagano_kw:
      _Nagano = 1;
      _Wilcox = _klog = _wlog = _emod = _numod = 0;
      _DynModel = new Dyn4P_KW();
      break;

    case nagano_log:
      _Nagano = _klog = _wlog = 1;
      _Wilcox = _numod = _emod = 0;
      _DynModel = new Dyn4P_logKW();
      break;

    case wilcox:
      _Nagano = _klog = _wlog = _emod = _numod = 0;
      _Wilcox = 1;
      _DynModel = new DynWilcox_KW();
      break;

    case wilcox_log:
      _Nagano = _emod = _numod = 0;
      _Wilcox = _klog = _wlog = 1;
      _DynModel = new DynWilcox_logKW();
      break;

    default:
      std::cout << "\033[1;31m\n=====================================================\n"
                << "   TURB_UTILS: unknown dynamical turbulence model " << _FileMap["RANS_dynamic"]
                << "\n=====================================================\n\033[0m";
      abort();
      break;
  }

  std::string ThermModel = (_FileMap["RANS_thermal"] != "") ? _FileMap["RANS_thermal"] : "default";
  int RANS_thermal = ::ThermTurbModelMap.at(ThermModel);
  _Kays = 0;
  _4PwithKays = 0;

  if (_FileMap["UseKays"] != "") { _4PwithKays = stoi(_FileMap["UseKays"]); }

  switch (RANS_thermal) {
    case default_therm:
      _ehmod = _khlog = _whlog = 0;
      _ThermModel = new Therm4P_KHWH();
      break;

    case nagano_keT:
      _ehmod = 1;
      _khlog = _whlog = 0;
      _ThermModel = new Therm4P_KHEH();
      break;

    case nagano_kwT:
      _ehmod = _khlog = _whlog = 0;
      _ThermModel = new Therm4P_KHWH();
      break;

    case nagano_logT:
      _ehmod = 0;
      _khlog = _whlog = 1;
      _ThermModel = new Therm4P_logKHWH();
      break;

    case kays: _Kays = 1; break;

    default:
      std::cout << "\033[1;31m\n=====================================================\n"
                << "   TURB_UTILS: unknown thermal turbulence model " << _FileMap["RANS_thermal"]
                << "\n=====================================================\n\033[0m";
      abort();
      break;
  }

  _IsFilled = true;

  if (_FileMap["YapCorrection"] != "") {
    _YapCorr = stoi(_FileMap["YapCorrection"]);
  } else {
    _YapCorr = 0;
  }

  if (_FileMap["DurbinConstrain"] != "") {
    _Durbin = stoi(_FileMap["DurbinConstrain"]);
  } else {
    _Durbin = 0;
  }

  if (_FileMap["Park"] != "") {
    _Park = stoi(_FileMap["Park"]);
  } else {
    _Park = 1;
  }

  if (_FileMap["Diameter"] != "") {
    _diameter = stod(_FileMap["Diameter"]);
  } else {
    _diameter = 0.1;
  }

  if (_FileMap["AvVelocity"] != "") {
    _vmid = stod(_FileMap["AvVelocity"]);
  } else {
    _vmid = 0.1;
  }

  if (_FileMap["utau"] != "") {
    _InputUtau = stod(_FileMap["utau"]);
  } else {
    _InputUtau = -1;  // if negative then profiles will be calculated with default utau method
  }

  if (_FileMap["SolveMuT"] != "") { _SolveMuT = stoi(_FileMap["SolveMuT"]); }

  if (_FileMap["SolveAlphaT"] != "") { _SolveAlphaT = stoi(_FileMap["SolveAlphaT"]); }

  if (__Proc == 0) {
    std::cerr << " \n =====================================================\n";
    std::cerr << "  TURB UTILS PARAMETERS   \n";
    std::cerr << "   _YapCorr    " << _YapCorr << std::endl;
    std::cerr << "   _Durbin     " << _Durbin << std::endl;
    std::cerr << "   _Park       " << _Park << std::endl;
    std::cerr << "   _diameter   " << _diameter << std::endl;
    std::cerr << "   _vmid       " << _vmid << std::endl;
    std::cerr << "   _InputUtau  " << _InputUtau << std::endl;
    std::cerr << "   _Nagano      " << _Nagano << std::endl;
    std::cerr << "   _Wilcox      " << _Wilcox << std::endl;
    std::cerr << "   _klog        " << _klog << std::endl;
    std::cerr << "   _wlog        " << _wlog << std::endl;
    std::cerr << "   _emod        " << _emod << std::endl;
    std::cerr << "   _numod       " << _numod << std::endl;
    std::cerr << " =====================================================\n";
  }

  _DynModel->SetNu(_nu);
  _DynModel->SetYap(_YapCorr);
  _DynModel->SetPark(_Park);
  _ThermModel->SetNu(_nu);
  _ThermModel->SetYap(_YapCorr);
  _ThermModel->SetPark(_Park);
  _ThermModel->SetIPr(__IPr);

  _klim = _DynModel->_LowerKappa;
  _wlim = _DynModel->_LowerOmega;
  _elim = _DynModel->_LowerOmega;
  _khlim = _ThermModel->_LowerKappah;
  _whlim = _ThermModel->_LowerKappah;
  _ehlim = _ThermModel->_LowerOmegah;

  return;
}

void TurbUtils::FillModelMap() {
  _DynamicModel["nagano_k"] = nagano_k;
  _DynamicModel["nagano_logk"] = nagano_logk;
  _DynamicModel["nagano_w"] = nagano_w;
  _DynamicModel["nagano_logw"] = nagano_logw;
  _DynamicModel["nagano_e"] = nagano_e;
  _DynamicModel["wilcox_k"] = wilcox_k;
  _DynamicModel["wilcox_logk"] = wilcox_logk;
  _DynamicModel["wilcox_nut"] = wilcox_nut;
  _DynamicModel["wilcox_w"] = wilcox_w;
  _DynamicModel["wilcox_logw"] = wilcox_logw;

  _DynamicModel["nagano_ke"] = nagano_ke;
  _DynamicModel["nagano_kw"] = nagano_kw;
  _DynamicModel["nagano_log"] = nagano_log;
  _DynamicModel["wilcox"] = wilcox;
  _DynamicModel["wilcox_log"] = wilcox_log;

  _ThermalModel["natural_kh"] = natural_kh;
  _ThermalModel["logarithmic_kh"] = logarithmic_kh;
  _ThermalModel["natural_omegah"] = natural_omegah;
  _ThermalModel["logarithmic_omegah"] = logarithmic_omegah;
  //     _ThermalModel["logarithmic_omegah"] = logarithmic_omegah;
  return;
}

double TurbUtils::CalcUtau(double vel_bound, double dist) {
  double umusk, ulog, ulin, utau, ulold, diff;
  double beta, vk;

  double yp;

  beta = 5.2;
  vk = 0.41;

  ulog = 0.;
  ulold = 11.6 / vel_bound;

  // Calculation of utau through the linear relation
  ulin = sqrt(vel_bound * _nu / dist);
  yp = ulin * dist / _nu;

  // Calculation of utau through the logarithmic relation
  diff = 100.;

  if (yp > 5.) {
    while (diff > 1.e-6) {
      ulog = vel_bound * vk / (log(exp(beta * vk) * dist * ulold / _nu));
      diff = fabs(ulold - ulog);
      ulold = ulog;
    }

    yp = ulog * dist / _nu;
  } else {
    ulog = 0.;
  }

  // Calculation of utau through the musker relation

  if (yp > 5. && yp < 40.) {
    diff = 100.;
    int cont = 0;
    double umuskold = ulog;

    while (diff > 1.e-6) {
      umusk = vel_bound / Musker(dist, umuskold);
      diff = fabs(umuskold - umusk);
      umuskold = umusk;
      cont++;

      if (cont > 3000) {
        umusk = 0.;
        break;
      }
    }
  }

  if (yp > 5.) {
    utau = max(ulog, umusk);

    if (umusk > 2. * ulog) { utau = ulog; }
  } else {
    utau = ulin;
  }

  return utau;
}

double TurbUtils::Musker(double dist, double utau) {
  double yplus = dist * utau / _nu;
  double vel = 5.424 * atan((2. * yplus - 8.15) / 16.7) + 4.1693 * log(yplus + 10.6) -
               0.8686 * log(yplus * yplus - 8.15 * yplus + 86) - 3.52;
  return vel;
}

// DYNAMICAL TURBULENCE SECTION ==========================================================

double TurbUtils::CalcMuTurb(double KappaAndOmega[], double dist, double vel_sp) {
  _MuTurb = _DynModel->CalcMuT(KappaAndOmega, dist);
  return _MuTurb;
}

void TurbUtils::CalcDynTurSourceAndDiss(
    double KappaAndOmega[], double dist, double vel_sp, double& muturb, double source[2], double diss[2],
    double div_g) {
  _MuTurb = _DynModel->CalcMuT(KappaAndOmega, dist);
  _DynModel->CalcDynSourceTerms(KappaAndOmega, source, diss, muturb, vel_sp, dist);
  return;
}

double TurbUtils::YapTerm(double& Kappa, double& Omega, double WallDist) {
  double TurbValues[2] = {Kappa, Omega};
  double yapterm = _DynModel->YapTerm(TurbValues, WallDist);
  return yapterm;
}

void TurbUtils::DynTurInitValues(double& kappa, double& omega, double WallDist, bool FlatProfile) {
  double TurbValues[2];

  if (FlatProfile) {
    _DynModel->DynTurInitKappaAndOmegaValues(TurbValues, _vmid, _diameter);
  } else {
    const double utau = (_InputUtau < 0) ? 0.5 * _nu / _BoundWallDist : _InputUtau;
    _DynModel->DynTurNearKappaAndOmegaValues(TurbValues, WallDist, utau);
  }

  kappa = TurbValues[0];
  omega = TurbValues[1];
  return;
}

void TurbUtils::DynTurInitValues(double& kappa, double& omega, double WallDist, double Utau) {
  double TurbValues[2];
  _DynModel->DynTurNearKappaAndOmegaValues(TurbValues, WallDist, Utau);
  kappa = TurbValues[0];
  omega = TurbValues[1];
  return;
}

void TurbUtils::DynTurNearWallValues(double& kappa, double& omega, double WallDist, double Utau) {
  double TurbValues[2];
  _DynModel->DynTurNearKappaAndOmegaValues(TurbValues, WallDist, Utau);
  kappa = TurbValues[0];
  omega = TurbValues[1];
  return;
}

// =======================================================================================

double TurbUtils::KaysPrt(double KappaAndOmega[], double dist) {
  double nut = CalcMuTurb(KappaAndOmega, dist);
  const double kays_prt = (0.85 + 0.7 * __IPr / (nut + 1.e-10));
  return kays_prt;
}

double TurbUtils::CalcAlphaTurb(double KappaAndOmega[], double TKappaAndOmega[], double dist) {
  _DynModel->CalcKappaAndOmega(KappaAndOmega);
  double alphaT = _ThermModel->CalcAlphaTurb(_DynModel->_KappaAndOmega, TKappaAndOmega, dist);
  return alphaT;
}

double TurbUtils::CalcAlphaTurb(double KappaAndOmega[], double TKappaAndOmega[], double dist, double Nut) {
  _DynModel->CalcKappaAndOmega(KappaAndOmega);
  double alphaT = _ThermModel->CalcAlphaTurb(_DynModel->_KappaAndOmega, TKappaAndOmega, dist);
  return alphaT;
}

void TurbUtils::TherTurInitValues(
    double& kappaT, double& omegaT, double WallDist, double RefDeltaT, double AvVel, double Diameter,
    bool FlatProfile) {
  double ThermTurInitKappaAndOmegaValues[2];
  if (FlatProfile) {
    _ThermModel->ThermTurInitKappaAndOmegaValues(ThermTurInitKappaAndOmegaValues, AvVel, Diameter, RefDeltaT);
  } else {
    const double Utau = 10. * _nu / _BoundWallDist;
    _ThermModel->ThermTurNearKappaHAndOmegaHValues(ThermTurInitKappaAndOmegaValues, WallDist, Utau);
  }
  kappaT = ThermTurInitKappaAndOmegaValues[0];
  omegaT = ThermTurInitKappaAndOmegaValues[1];
  return;
}

void TurbUtils::ThermTurNearWallValues(double& kappa, double& omega, double WallDist, double Utau) {
  double ThermTurInitKappaAndOmegaValues[2];
  _ThermModel->ThermTurNearKappaHAndOmegaHValues(ThermTurInitKappaAndOmegaValues, WallDist, Utau);
  kappa = ThermTurInitKappaAndOmegaValues[0];
  omega = ThermTurInitKappaAndOmegaValues[1];

  return;
}

void TurbUtils::CalcThermTurSourceAndDiss(
    double KappaAndOmega[],   // Dynamic Turbulence
    double TKappaAndOmega[],  // Thermal Turbulence
    double dist,              // Wall Distance
    double sp,                // Squared value of velocity derivatives
    double st,                // Modulus of temperature gradient
    double& alphaturb,        // Eddy thermal diffusivity
    double source[],          // Source terms
    double diss[],            // Dissipation terms
    double meccterm[]) {      // Mechanical contribution

  _DynModel->CalcKappaAndOmega(KappaAndOmega);
  double mut = _DynModel->CalcMuT(KappaAndOmega, dist);
  _ThermModel->CalcThermSourceTerms(
      _DynModel->_KappaAndOmega, TKappaAndOmega, source, diss, meccterm, st, sp, mut, alphaturb, dist);

  return;
}

void TurbUtils::CalcThermTurSourceAndDiss(
    double KappaAndOmega[],   // Dynamic Turbulence
    double TKappaAndOmega[],  // Thermal Turbulence
    double dist,              // Wall Distance
    double sp,                // Squared value of velocity derivatives
    double st,                // Modulus of temperature gradient
    double& alphaturb,        // Eddy thermal diffusivity
    double source[],          // Source terms
    double diss[],            // Dissipation terms
    double meccterm[],
    double Nut) {  // Mechanical contribution

  _DynModel->CalcKappaAndOmega(KappaAndOmega);
  _ThermModel->CalcThermSourceTerms(
      _DynModel->_KappaAndOmega, TKappaAndOmega, source, diss, meccterm, st, sp, Nut, alphaturb, dist);

  return;
}

void TurbUtils::PrintStatus() {
  ofstream output;
  output.open("SimulationSpecs.txt");
  output << " MESSAGE PRINTED BY TURBUTILS CLASS " << endl;
  output << endl << " TurbUtils filled: " << _IsFilled << endl;
  output << "---------------------------------------" << endl;
  output << " Solved Equations: " << endl;
  output << "\t Navier Stokes:        " << _SolveNS << endl;
  output << "\t Temperature:          " << _SolveT << endl;
  output << "\t Dynamical Turbulence: " << _SolveTBK << endl;
  output << "\t Thermal Turbulence:   " << _SolveTTBK << endl;
  output << "---------------------------------------" << endl;
  output << " Realizability and Constrains:" << endl;
  output << "\t Durbin Realizability Constrain: " << _Durbin << endl;
  output << "\t Yap Correction:                 " << _YapCorr << endl;
  output << "---------------------------------------" << endl;
  output << " Geometry Parameters: " << endl;
  output << "\t Fixed Wall Distance: " << _BoundWallDist << endl;
  output << "---------------------------------------" << endl;
  output << " Physical Properties: " << endl;
  output << "\t Kinematic Viscosity: " << _nu << endl;
  output << "\t Thermal diffusivity: " << _alpha << endl;

  return;
}

void TurbUtils::read_file() {  // READING Tparameter.in ======================================
  std::string string_value;    // read double, string, dummay
  std::ostringstream file, file1, file2;
  std::vector<std::string> FILES;

  if (__NMeshes > 0) {
    file << getenv("APP_PATH") << "/DATA/DATA" << std::to_string(__MeshID) << "/Turbulence.in";
    file1 << getenv("APP_PATH") << "/DATA/DATA" << std::to_string(__MeshID) << "/GeometrySettings.in";
    file2 << getenv("APP_PATH") << "/DATA/DATA" << std::to_string(__MeshID) << "/MaterialProperties.in";
  } else {
    file << getenv("APP_PATH") << "/DATA/Turbulence.in";
    file1 << getenv("APP_PATH") << "/DATA/GeometrySettings.in";
    file2 << getenv("APP_PATH") << "/DATA/MaterialProperties.in";
  }

  FILES.push_back(file1.str());
  FILES.push_back(file.str());
  FILES.push_back(file2.str());

  for (int i = 0; i < FILES.size(); i++) {
    std::ifstream fin;
    fin.open(FILES[i].c_str());  // stream file
    std::string buf = "";

    if (fin.is_open()) { std::cout << "\nInit Reading = " << FILES[i] << std::endl; }

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
          std::cout << buf << "\t" << string_value << std::endl;
        }

        fin >> buf;
      }
    } else {
      std::cerr << "TurbUtils.read_file(): no parameter file found\t ->" << FILES[i] << std::endl;
      abort();
    }

    std::cout << "TurbUtils.read_file() End Reading file " << FILES[i] << std::endl;
    fin.close();
  }

  print_par();
  return;
}  // END READING ==========================================================================

void TurbUtils::print_par() {
  std::cout
      << "\033[038;5;196;1m \n----------------------------------------------\n  TURBUTILS MAP \n \033[0m"
      << std::endl;

  for (std::map<std::string, std::string>::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }

  std::cout << "\033[038;5;196;1m \n----------------------------------------------\n \033[0m" << std::endl;

  return;
}
