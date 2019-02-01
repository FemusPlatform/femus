#include <sstream>
#include <vector>
#include <map>
#ifndef __TurbUtils__
#define __TurbUtils__

#include "TurbModels.h"

class TurbUtils{
private:  
  // NAGANO-DERIVED TURBULENCE MODEL COEFFICIENTS
  const double __CMU = 0.09;
  const double __CP1 = 1.025;
  const double __CP2 = 0.9;
  const double __CD1 = 1.1;
  const double __CD2 = 1.9;
  const double __C20 = 1.9;
  const double __C10 = 1.5;
  const double __Prdl_inf = 0.75;
  
  // WILCOX-DERIVED TURBULENCE MODEL COEFFICIENTS
  const double __AW = 13./25.;
  const double __AN = 1. - __AW;
  const double __BETAS = 8./100.;
  const double __BETAW = 3./40.;
  const double __BETAN = __BETAS - __BETAW;
  
  double _klim, _wlim, _khlim, _whlim, _elim, _ehlim;
  
  // DYNAMICAL TURBULENCE FUNCTIONS
  double __Rt    ;
  double __Rd    ;
  double __fmu   ;
  double __fcorr ;
  double __Ret   ;
  
  // THERMAL TURBULENCE FUNCTIONS
  double __IPr      ;  
  double __rT       ;
  double __F1t      ;
  double __F2at     ;
  double __F2bt     ;    
 
  const int __Proc;
  int       _MeshID;
  
  int _LibToMed_3D[27]    = {7, 4, 5, 6, 3, 0, 1, 2, 19, 16, 17, 18, 11, 8, 9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};
  int _LibToMed_2D[9]    = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  int **_CanElemMap;

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
  std::map<int,int> _MgToMed;
  
  double _Tmean, _Fraction;
  
  
//   FLAGS 
  bool   _IsFilled = false;
  bool   _MuTurbInitialized=false;
  bool   _SolveNS, _SolveT, _SolveTBK, _SolveTTBK;
  int    _YapCorr, _Durbin, _Park;
  bool   _IsWallDistSet = false;
  
  
  std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
  
  std::map<std::string, DynTurbModel> _DynamicModel;  
  std::map<std::string, ThermTurbModel> _ThermalModel; 
  
public:
  TurbUtils();  

  ~TurbUtils();
  void read_file();
  void print_par();
  
  double CalcUtau(double vel_bound, 
		  double dist);
  double CalcMuTurb(double KappaAndOmega[], 
		    double dist,
		   double vel_sp = 1.);
    
  void CalcDynTurSourceAndDiss(double KappaAndOmega[],
			       double dist, 
			       double vel_sp ,
			       double &muturb, 
			       double source[2], 
			       double diss[2],
                   double div_g=0.
			      );
  
  
  double CalcAlphaTurb(double KappaAndOmega[],
		       double TKappaAndOmega[],
		       double dist,
               double Nut);
  
  double CalcAlphaTurb(double KappaAndOmega[],
		       double TKappaAndOmega[],
		       double dist);
  
  
  double KaysPrt (double KappaAndOmega[], double dist);
    
  void CalcThermTurSourceAndDiss(double KappaAndOmega[], 
				 double TKappaAndOmega[], 
				 double dist, 
				 double sp, 
				 double st, 
				 double &alphaturb, 
				 double source[2], 
				 double diss[2], 
				 double meccterm[2]);
  
    void CalcThermTurSourceAndDiss(double KappaAndOmega[], 
				 double TKappaAndOmega[], 
				 double dist, 
				 double sp, 
				 double st, 
				 double &alphaturb, 
				 double source[2], 
				 double diss[2], 
				 double meccterm[2],
                 double Nut);
  
  void FillParameters();
  void FillModelMap();
  
  void PrintStatus();
  
  void DynTurNearWallValues(double & kappa, double & omega, double WallDist, double Utau);
  void ThermTurNearWallValues ( double & kappa, double & omega, double WallDist, double Utau );
  
  double YapTerm( double & kappa, double & omega, double WallDist );
  
  void DynTurInitValues(double & kappa, double & omega, double WallDist, bool FlatProfile);
  void DynTurInitValues(double & kappa, double & omega, double WallDist, double Utau);
  void TherTurInitValues(double & kappaT, double & omegaT, double WallDist, double NormFlux, double AvVel, double Diameter, bool FlatProfile);
  void CalcWallFuncKappaAndOmega(double KappaAndOmega[], int NodeOnBound, double WallDist, double utau);
  void CalcWallFuncThermalKappaAndOmega(double KappaAndOmega[], int NodeOnBound, double WallDist, double utau);
  inline bool CheckIfFilled(){return _IsFilled;}
  inline void SetAvVel(double vel){ _vmid = vel;}
  inline void SetDiameter(double diameter){ _diameter = diameter;}
  double Musker (double dist, double utau) ;  

};


#endif
