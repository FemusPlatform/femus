#include <sstream>
#include <vector>
#include <map>
#ifndef __TurbUtils__
#define __TurbUtils__


namespace MEDCoupling {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}


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
 
  const int __Levels;
  const int __Proc;
  int       _MeshID;
  
  int _LibToMed_3D[27]    = {7, 4, 5, 6, 3, 0, 1, 2, 19, 16, 17, 18, 11, 8, 9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};
  int _LibToMed_2D[9]    = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  int **_CanElemMap;
//   int CanonicalElementMaps
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
  int _Nagano, _Wilcox;
  int _numod, _emod;
  double _vmid, _diameter;
  double _DynUnderRel, _ThermUnderRel;
  double _InputUtau;
  std::map<int,int> _MgToMed;
  
//   FLAGS 
  bool   _IsFilled = false;
  bool   _MuTurbInitialized=false;
  bool   _SolveNS, _SolveT, _SolveTBK, _SolveTTBK;
  int    _YapCorr, _Durbin, _Park;
  bool   _IsWallDistSet = false;
  
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> _NodeWallDist;
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> _NodeMap;
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> _MuTurbField;
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> _AlphaTurbField;
  

  inline void SetNodeWallDist(MEDCoupling::MEDCouplingFieldDouble * NodeWallDist) {_NodeWallDist.push_back(NodeWallDist);};
  inline MEDCoupling::MEDCouplingFieldDouble * GetMuTurbField(int Level){return _MuTurbField[Level];};
  
  void SetMuTurbFieldAtLevel(int Level, MEDCoupling::MEDCouplingFieldDouble * MuTurb);
  void SetAlphaTurbFieldAtLevel(int Level, MEDCoupling::MEDCouplingFieldDouble * AlphaTurb);
  void SetWallDistAtLevel(int Level, MEDCoupling::MEDCouplingFieldDouble * WallDist);
  
  std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
  
  enum DynTurbModel{
//     natural_k = 10,
//     logarithmic_k = 11,
//     natural_omega = 20,
//     logarithmic_omega = 21
    // NAGANO BASED MODELS
    nagano_ke   =1,
    nagano_kw   =2,
    nagano_log  =3,
    wilcox      =4,
    wilcox_log  =5,
    nagano_k    =10,
    nagano_logk =11,
    nagano_w    =10,
    nagano_logw =11,
    nagano_e    =12,
    // WILCOX BASED MODELS
    wilcox_k    =20,
    wilcox_logk =21,
    wilcox_nut  =22,
    wilcox_w    =20,
    wilcox_logw =21
  };
  enum ThermTurbModel{
    natural_kh = 10,
    logarithmic_kh = 11,
    natural_omegah = 20,
    logarithmic_omegah = 21
  };
  
  
  std::map<std::string, DynTurbModel> _DynamicModel;  
  std::map<std::string, ThermTurbModel> _ThermalModel; 
  
  
public:
  TurbUtils();  
  TurbUtils( double wall_dist, double TurbModel[], double nu, double alpha = 1.e-4);
  TurbUtils(int proc,
            int levels,
            std::vector<MEDCoupling::MEDCouplingFieldDouble *>NodeMap,
            bool DynTurb, 
            bool TherTurb);
  TurbUtils ( int proc,
            int levels,
            std::vector<MEDCoupling::MEDCouplingFieldDouble *>NodeMap,
            bool DynTurb,
            bool TherTurb,
	        int MeshID);
  ~TurbUtils();
  void read_file();
  void print_par();
  
  double CalcUtau(double vel_bound, 
		  double dist);
  double CalcMuTurb(double KappaAndOmega[], 
		    double dist,
		   double vel_sp = 1.);
  void CalcMuTurb(
    MEDCoupling::MEDCouplingFieldDouble * FirstDynVar,
    MEDCoupling::MEDCouplingFieldDouble * SecDynVar,
    int Level
  );
  void CalcDynTurSourceAndDiss(double KappaAndOmega[],
			       double dist, 
			       double vel_sp ,
			       double &muturb, 
			       double source[2], 
			       double diss[2],
			       double div_g =0
			      );
  double CalcAlphaTurb(double KappaAndOmega[],
		       double TKappaAndOmega[],
		       double dist);
  
  void CalcAlphaTurb(
    MEDCoupling::MEDCouplingFieldDouble * FirstDynVar,
    MEDCoupling::MEDCouplingFieldDouble * SecDynVar,
    MEDCoupling::MEDCouplingFieldDouble * FirstThermVar,
    MEDCoupling::MEDCouplingFieldDouble * SecThermVar,
    int Level);  
  
  void CalcThermTurSourceAndDiss(double KappaAndOmega[], 
				 double TKappaAndOmega[], 
				 double dist, 
				 double sp, 
				 double st, 
				 double &alphaturb, 
				 double source[2], 
				 double diss[2], 
				 double meccterm[2]);
  
  void FillParameters(double wall_dist, 
		      double TurbModel[],
		      double nu, 
		      double alpha=1.e-4);
  void FillParameters(double wall_dist, 
		      std::string val1, 
		      std::string val2, 
		      std::string val3, 
		      std::string val4, 
		      double nu=1.e-4, 
		      double alpha=1.e-4);
  void FillParameters(double wall_dist, 
		      std::vector<std::string> TurbModel,
		      std::vector<std::string> SolveEqs,
		      std::vector<std::string> Constrain,
		      double UnderRelaxation[],
		      double nu=1.e-4, 
		      double alpha=1.e-4
		      );
  void FillParameters();
  void FillModelMap();
  
  void PrintStatus(std::vector<std::string> TurbModel);
  
  void DynTurNearWallValues(double & kappa, double & omega, double WallDist, double Utau);
  
  
  void DynTurInitValues(double & kappa, double & omega, double WallDist, bool FlatProfile);
  void DynTurInitValues(double & kappa, double & omega, double WallDist, double Utau);
  void TherTurInitValues(double & kappaT, double & omegaT, double WallDist, double NormFlux, double AvVel, double Diameter, bool FlatProfile);
  void CalcWallFuncKappaAndOmega(double KappaAndOmega[], int NodeOnBound, double WallDist, double utau);
  void CalcWallFuncThermalKappaAndOmega(double KappaAndOmega[], int NodeOnBound, double WallDist, double utau);
  inline bool CheckIfFilled(){return _IsFilled;}
  inline void SetAvVel(double vel){ _vmid = vel;}
  inline void SetDiameter(double diameter){ _diameter = diameter;}
  double Musker (double dist, double utau) ;  
  
  void GetLevelElemMuTurb(int iel, int level, double MuTurb[]);
  void GetLevelElemAlphaTurb(int iel, int level, double AlphaTurb[]);
  void GetLevelElemNodeWallDist(int iel, int level, double NodeWallDist[]);
  
  MEDCoupling::MEDCouplingFieldDouble * BuildInitTurbField(int TurbVar);
};


#endif
