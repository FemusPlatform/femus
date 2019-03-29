#ifndef __Utils_DynTurModels_Wilcox_h__
#define __Utils_DynTurModels_Wilcox_h__

#include <math.h>
#include <algorithm>
#include "Utils_DynTurbModels.h"

using namespace std;


struct DynWilcox:DYNturModels
{

  const double __AW = 13. / 25.;
  const double __AN = 1. - __AW;
  const double __BETAS = 8. / 100.;
  const double __BETAW = 3. / 40.;
  const double __BETAN = __BETAS - __BETAW;

    DynWilcox ():__AW (13. / 25.), __AN (1. - 13. / 25.), __BETAS (8. / 100.),
    __BETAW (3. / 40.), __BETAN (8. / 100. - 3. / 40.)
  {
  };

  double CalcMuT (double *TurbVariables, double ydist);
  virtual void CalcDynSourceTerms (double *TurbVariables,
				   double *SourceTerms, double *DissTerms,
				   double MuT, double VelGradMod,
				   double ydist = 1.e-5)
  {
  };
  virtual void CalcKappaAndOmega (double *TurbVariables)
  {
  };
  virtual void ConvertKandWtoLocal(double *TurbValues){};

};


struct DynWilcox_KW:DynWilcox
{
  void CalcDynSourceTerms (double *TurbVariables,
			   double *SourceTerms, double *DissTerms,
			   double MuT, double VelGradMod, double ydist =
			   1.e-5);
  void CalcKappaAndOmega (double *TurbVariables);
  void ConvertKandWtoLocal(double *TurbValues);
};

struct DynWilcox_logKW:DynWilcox
{
  void CalcDynSourceTerms (double *TurbVariables,
			   double *SourceTerms, double *DissTerms,
			   double MuT, double VelGradMod, double ydist =
			   1.e-5);
  void CalcKappaAndOmega (double *TurbVariables);
  void ConvertKandWtoLocal(double *TurbValues);
};


#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
