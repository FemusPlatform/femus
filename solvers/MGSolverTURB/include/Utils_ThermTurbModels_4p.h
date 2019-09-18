#ifndef __Utils_ThermTurModels_4p_h__
#define __Utils_ThermTurModels_4p_h__

#include <math.h>
#include <algorithm>
#include "Utils_DynTurbModels_4p.h"

using namespace std;

struct Therm4P : Dyn4P {
  const double __CP1;
  const double __CP2;
  const double __CD1;
  const double __CD2;
  const double __Prdl_inf;

  double _LowerKappah, _LowerOmegah;
  double _rT, _F1t, _F2at, _F2bt, _IPr;
  double _KappaHAndOmegaH[2];

  Therm4P()
      : __CP1(1.025),
        __CP2(0.9),
        __CD1(1.1),
        __CD2(1.9),
        __Prdl_inf(1.3),
        _LowerKappah(1.e-7),
        _LowerOmegah(1.e-7){};

  /*! Function called to compute values of eddy thermal diffusivity \f$\alpha_t\f$.
   *  The function is inherited by derived #Therm4P_logKHWH, #Therm4P_KHWH and #Therm4P_KHEH
   *  structs
   */
  double CalcAlphaTurb(
      double KappaAndOmega[],
      /**< Input variables must be real \f$k\f$ and \f$\omega\f$ */
      double TKappaAndOmega[],
      /**< Input variables for thermal turbulence */
      double dist
      /**< Wall distance */
  );

  /*! Function called to compute values of real \f$k_\theta\f$ and real
   *  \f$\omega_\theta\f$. The function is specialized by #Therm4P struct
   *  children: #Therm4P_logKHWH, #Therm4P_KHWH and #Therm4P_KHEH
   */
  virtual void CalcKappaHAndOmegaH(double* TurbVariables){};
  virtual void ConvertKHandWHtoLocal(double* TurbVariables){};

  void inline SetIPr(double IPr) { _IPr = IPr; };

  virtual void CalcThermSourceTerms(
      double* KappaAndOmega,      /**< Input variables must be real \f$k\f$ and \f$\omega\f$ */
      double* ThermTurbVariables, /**< Input variables for thermal turbulence */
      double* SourceTerms, double* DissTerms, double* MeccTerms, double TempGradMod, double VelGradMod,
      double muT, double alphaT, double ydist = 1.e-5){};

  void ThermTurNearKappaHAndOmegaHValues(
      double* ThermTurInitKappaAndOmegaValues, double WallDist, double Utau);
  void ThermTurInitKappaAndOmegaValues(
      double* ThermTurInitKappaAndOmegaValues, double vel, double diameter, double RefDeltaT);
};

struct Therm4P_logKHWH : Therm4P {
  /*! Function called to compute values of real \f$k_\theta\f$ and real
   *  \f$\omega_\theta\f$.
   */
  void CalcKappaHAndOmegaH(double* TurbVariables);

  void CalcThermSourceTerms(
      double* KappaAndOmega,      /**< Input variables must be real \f$k\f$ and \f$\omega\f$ */
      double* ThermTurbVariables, /**< Input variables for thermal turbulence */
      double* SourceTerms, double* DissTerms, double* MeccTerms, double TempGradMod, double VelGradMod,
      double muT, double alphaT, double ydist = 1.e-5);

  void ConvertKHandWHtoLocal(double* TurbVariables);
};

struct Therm4P_KHWH : Therm4P {
  /*! Function called to compute values of real \f$k_\theta\f$ and real
   *  \f$\omega_\theta\f$.
   */
  void CalcKappaHAndOmegaH(double* TurbVariables);

  void CalcThermSourceTerms(
      double* KappaAndOmega,      /**< Input variables must be real \f$k\f$ and \f$\omega\f$ */
      double* ThermTurbVariables, /**< Input variables for thermal turbulence */
      double* SourceTerms, double* DissTerms, double* MeccTerms, double TempGradMod, double VelGradMod,
      double muT, double alphaT, double ydist = 1.e-5);

  void ConvertKHandWHtoLocal(double* TurbVariables);
};

struct Therm4P_KHEH : Therm4P {
  /*! Function called to compute values of real \f$k_\theta\f$ and real
   *  \f$\omega_\theta\f$.
   */
  void CalcKappaHAndOmegaH(double* TurbVariables);

  void CalcThermSourceTerms(
      double* KappaAndOmega,      /**< Input variables must be real \f$k\f$ and \f$\omega\f$ */
      double* ThermTurbVariables, /**< Input variables for thermal turbulence */
      double* SourceTerms, double* DissTerms, double* MeccTerms, double TempGradMod, double VelGradMod,
      double muT, double alphaT, double ydist = 1.e-5);

  void ConvertKHandWHtoLocal(double* TurbVariables);
};

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
