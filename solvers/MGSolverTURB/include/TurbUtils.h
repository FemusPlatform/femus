#include <sstream>
#include <vector>
#include <map>
#ifndef __TurbUtils__
#define __TurbUtils__

#include "TurbModels.h"

class DYNturModels;
class Therm4P;

class TurbUtils {
 private:
  double _klim, _wlim, _khlim, _whlim, _elim, _ehlim;
  double __IPr;

  double __CMU = 0.09;

  const int __Proc;
  int __MeshID;
  int __NMeshes;

 public:
  double _utau;
  double _kwall;
  double _wwall;

  double _BoundWallDist;
  double _MuTurb;
  double _AlphaTurb;
  double _nu;
  double _alpha;
  int _wlog, _klog, _khlog, _whlog;
  int _Nagano, _Wilcox, _Kays, _4PwithKays;
  int _numod, _emod, _ehmod;
  double _vmid, _diameter;
  double _DynUnderRel, _ThermUnderRel;
  double _InputUtau;
  std::map<int, int> _MgToMed;

  double _Tmean, _Fraction;

  //   FLAGS
  bool _IsFilled = false;
  bool _MuTurbInitialized = false;
  bool _SolveNS, _SolveT, _SolveTBK, _SolveTTBK;
  int _YapCorr, _Durbin, _Park;
  bool _IsWallDistSet = false;
  int _SolveMuT;
  int _SolveAlphaT;

  std::map<std::string, std::string> _FileMap;  /// String map containing Tproperties.in parameters

  std::map<std::string, DynTurbModel> _DynamicModel;
  std::map<std::string, ThermTurbModel> _ThermalModel;

  DYNturModels* _DynModel;
  Therm4P* _ThermModel;

 public:
  TurbUtils(int MeshId = 0, int nMeshes = 0);

  ~TurbUtils();
  void read_file();
  void print_par();

  double CalcUtau(double vel_bound, double dist);
  double CalcMuTurb(double KappaAndOmega[], double dist, double vel_sp = 1.);

  void CalcDynTurSourceAndDiss(
      double KappaAndOmega[], double dist, double vel_sp, double& muturb, double source[2], double diss[2],
      double div_g = 0.);

  double CalcAlphaTurb(double KappaAndOmega[], double TKappaAndOmega[], double dist, double Nut);

  double CalcAlphaTurb(double KappaAndOmega[], double TKappaAndOmega[], double dist);

  double KaysPrt(double KappaAndOmega[], double dist);

  void CalcThermTurSourceAndDiss(
      double KappaAndOmega[], double TKappaAndOmega[], double dist, double sp, double st, double& alphaturb,
      double source[2], double diss[2], double meccterm[2]);

  void CalcThermTurSourceAndDiss(
      double KappaAndOmega[], double TKappaAndOmega[], double dist, double sp, double st, double& alphaturb,
      double source[2], double diss[2], double meccterm[2], double Nut);

  void FillParameters();
  void FillModelMap();

  void PrintStatus();

  void DynTurNearWallValues(double& kappa, double& omega, double WallDist, double Utau);
  void ThermTurNearWallValues(double& kappa, double& omega, double WallDist, double Utau);

  double YapTerm(double& kappa, double& omega, double WallDist);

  void DynTurInitValues(double& kappa, double& omega, double WallDist, bool FlatProfile);
  void DynTurInitValues(double& kappa, double& omega, double WallDist, double Utau);
  void TherTurInitValues(
      double& kappaT, double& omegaT, double WallDist, double RefDeltaT, double AvVel, double Diameter,
      bool FlatProfile);

  inline bool CheckIfFilled() { return _IsFilled; }
  inline void SetAvVel(double vel) { _vmid = vel; }
  inline void SetDiameter(double diameter) { _diameter = diameter; }
  double Musker(double dist, double utau);

  inline double GetKlim() { return _klim; }
  inline double GetElim() { return _elim; }
  inline double GetWlim() { return _wlim; }
  inline double GetKHlim() { return _khlim; }
  inline double GetEHlim() { return _ehlim; }
  inline double GetWHlim() { return _whlim; }
};

#endif
