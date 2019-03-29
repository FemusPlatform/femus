#ifndef __Utils_DynTurModels_4p_h__
#define __Utils_DynTurModels_4p_h__

#include <math.h>
#include <algorithm>
#include "Utils_DynTurbModels.h"

using namespace std;



struct Dyn4P:DYNturModels
{

  const double __CMU;
  const double __C20;
  const double __C10;

  double _Rt, _Rd, _fmu, _fcorr, _Ret;
  
    Dyn4P ():__CMU (0.09), __C20 (1.9), __C10 (1.5)
  {
  };

  void CalcDynCoefficients (double ydist);

  double CalcMuT (double *TurbVariables, double ydist = 1.e-5);

  virtual void CalcKappaAndOmega (double *TurbVariables)
  {
  };

  virtual void CalcDynSourceTerms (double *TurbVariables,
				   double *SourceTerms, double *DissTerms,
				   double MuT, double VelGradMod,
				   double ydist = 1.e-5)
  {
  };
  virtual void ConvertKandWtoLocal(double *TurbValues){};

};

struct Dyn4P_logKW:Dyn4P
{
  void CalcDynSourceTerms (double *TurbVariables,
			   double *SourceTerms, double *DissTerms, double MuT,
			   double VelGradMod, double ydist = 1.e-5);
  void CalcKappaAndOmega (double *TurbVariables);
  void ConvertKandWtoLocal(double *TurbValues);
};

struct Dyn4P_KW:Dyn4P
{
  void CalcDynSourceTerms (double *TurbVariables,
			   double *SourceTerms, double *DissTerms, double MuT,
			   double VelGradMod, double ydist = 1.e-5);
  void CalcKappaAndOmega (double *TurbVariables);
  void ConvertKandWtoLocal(double *TurbValues);
};

struct Dyn4P_KE:Dyn4P
{
  void CalcDynSourceTerms (double *TurbVariables,
			   double *SourceTerms, double *DissTerms, double MuT,
			   double VelGradMod, double ydist = 1.e-5);
  void CalcKappaAndOmega (double *TurbVariables);
  void ConvertKandWtoLocal(double *TurbValues);
};


#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
