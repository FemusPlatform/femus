#include "Utils_DynTurbModels_4p.h"
#include <iostream>

// FATHER Dyn4P CALC COEFFICIENTS
void Dyn4P::CalcDynCoefficients(double ydist) {
  const double kappa = _KappaAndOmega[0];
  const double omega = _KappaAndOmega[1];

  _Ret = kappa / (_nu * omega);  //  viscosity ratio
  _Rt = _Ret / __CMU;            // turbulent Reynolds number
  _Rd = ydist * sqrt(kappa / sqrt(_Rt)) / _nu;
  _fmu = (1. - exp(-1. * _Rd / 14.)) * (1. - exp(-1. * _Rd / 14.));
  _fcorr = 1. + 5. / pow(_Rt, 0.75) * exp(-1. * _Rt * _Rt / 40000.);

  return;
}

// 4 PARAMETERS: LOG KAPPA AND LOG OMEGA MODEL -------------------------------------------

void Dyn4P_logKW::CalcKappaAndOmega(double* TurbVariables) {
  double kappa = exp(TurbVariables[0]), omega = exp(TurbVariables[1]);
  _KappaAndOmega[0] = (kappa > _LowerKappa) ? kappa : _LowerKappa;
  _KappaAndOmega[1] = (omega > _LowerOmega) ? omega : _LowerOmega;
  return;
}

void Dyn4P_logKW::ConvertKandWtoLocal(double* TurbValues) {
  TurbValues[0] = log(_KappaAndOmega[0]);
  TurbValues[1] = log(_KappaAndOmega[1]);

  return;
}

double Dyn4P::CalcMuT(double* TurbVariables, double ydist) {
  CalcKappaAndOmega(TurbVariables);
  CalcDynCoefficients(ydist);

  double MUT = _Ret * _fcorr * _fmu;
  return MUT;
}

void Dyn4P_logKW::CalcDynSourceTerms(
    double* TurbVariables, double* SourceTerms, double* DissTerms, double MuT, double VelGradMod,
    double ydist) {
  CalcKappaAndOmega(TurbVariables);
  CalcDynCoefficients(ydist);

  double prod_k = 0.5 * _nu * MuT * VelGradMod;
  if (_Park == 1) { prod_k = Park_KSource(prod_k, _KappaAndOmega[0], VelGradMod); }

  const double f_exp = (1. - exp(-1. * _Rd / (3.1))) * (1. - exp(-1. * _Rd / (3.1))) *
                       (1. - 0.3 * exp(-1. * _Rt * _Rt / 42.25));
  DissTerms[0] = __CMU * _KappaAndOmega[1];
  DissTerms[1] = __CMU * _KappaAndOmega[1] * (__C20 * f_exp - 1);
  SourceTerms[0] = prod_k / _KappaAndOmega[0];
  SourceTerms[1] = (__C10 - 1.) * prod_k / _KappaAndOmega[0];
  if (_YapCorr == 1) {
    double KappaAndEpsilon[2];
    KappaAndEpsilon[0] = _KappaAndOmega[0];
    KappaAndEpsilon[1] = __CMU * _KappaAndOmega[0] * _KappaAndOmega[1];

    // source term for epsilon equation
    double yap_term = YapTerm(KappaAndEpsilon, ydist);
    yap_term *= 1. / _KappaAndOmega[1];

    SourceTerms[1] += yap_term;
  }

  return;
}

// 4 PARAMETERS: KAPPA AND OMEGA MODEL ---------------------------------------------------

void Dyn4P_KW::CalcKappaAndOmega(double* TurbVariables) {
  _KappaAndOmega[0] = (TurbVariables[0] > _LowerKappa) ? TurbVariables[0] : _LowerKappa;
  _KappaAndOmega[1] = (TurbVariables[1] > _LowerOmega) ? TurbVariables[1] : _LowerOmega;
  return;
}

void Dyn4P_KW::ConvertKandWtoLocal(double* TurbValues) {
  TurbValues[0] = (_KappaAndOmega[0]);
  TurbValues[1] = (_KappaAndOmega[1]);

  return;
}

void Dyn4P_KW::CalcDynSourceTerms(
    double* TurbVariables, double* SourceTerms, double* DissTerms, double MuT, double VelGradMod,
    double ydist) {
  CalcKappaAndOmega(TurbVariables);
  CalcDynCoefficients(ydist);

  double prod_k = 0.5 * _nu * MuT * VelGradMod;
  if (_Park == 1) { prod_k = Park_KSource(prod_k, _KappaAndOmega[0], VelGradMod); }

  const double f_exp = (1. - exp(-1. * _Rd / (3.1))) * (1. - exp(-1. * _Rd / (3.1))) *
                       (1. - 0.3 * exp(-1. * _Rt * _Rt / 42.25));
  DissTerms[0] = __CMU * _KappaAndOmega[1];
  DissTerms[1] = __CMU * _KappaAndOmega[1] * (__C20 * f_exp - 1);
  SourceTerms[0] = prod_k;
  SourceTerms[1] = (__C10 - 1.) * prod_k * _KappaAndOmega[1] / _KappaAndOmega[0];
  if (_YapCorr == 1) {
    double KappaAndEpsilon[2];
    KappaAndEpsilon[0] = _KappaAndOmega[0];
    KappaAndEpsilon[1] = __CMU * _KappaAndOmega[0] * _KappaAndOmega[1];

    // source term for epsilon equation
    double yap_term = YapTerm(KappaAndEpsilon, ydist);
    SourceTerms[1] += yap_term;
  }

  return;
}

// 4 PARAMETERS: KAPPA AND EPSILON MODEL -------------------------------------------------

void Dyn4P_KE::CalcKappaAndOmega(double* TurbVariables) {
  double omega = (TurbVariables[1] / (TurbVariables[0] * __CMU));
  _KappaAndOmega[0] = (TurbVariables[0] > _LowerKappa) ? TurbVariables[0] : _LowerKappa;
  _KappaAndOmega[1] = (omega > _LowerOmega) ? omega : _LowerOmega;
  return;
}

void Dyn4P_KE::ConvertKandWtoLocal(double* TurbValues) {
  TurbValues[0] = (_KappaAndOmega[0]);
  TurbValues[1] = (_KappaAndOmega[1] * _KappaAndOmega[0] * __CMU);

  return;
}

void Dyn4P_KE::CalcDynSourceTerms(
    double* TurbVariables, double* SourceTerms, double* DissTerms, double MuT, double VelGradMod,
    double ydist) {
  CalcKappaAndOmega(TurbVariables);
  CalcDynCoefficients(ydist);

  double prod_k = 0.5 * _nu * MuT * VelGradMod;
  if (_Park == 1) { prod_k = Park_KSource(prod_k, _KappaAndOmega[0], VelGradMod); }

  const double f_exp = (1. - exp(-1. * _Rd / (3.1))) * (1. - exp(-1. * _Rd / (3.1))) *
                       (1. - 0.3 * exp(-1. * _Rt * _Rt / 42.25));
  DissTerms[0] = __CMU * _KappaAndOmega[1];
  DissTerms[1] = __CMU * _KappaAndOmega[1] * __C20 * f_exp;
  SourceTerms[0] = prod_k;
  SourceTerms[1] = __CMU * __C10 * prod_k * _KappaAndOmega[1];
  if (_YapCorr == 1) {
    double KappaAndEpsilon[2];
    KappaAndEpsilon[0] = _KappaAndOmega[0];
    KappaAndEpsilon[1] = __CMU * _KappaAndOmega[0] * _KappaAndOmega[1];

    // source term for epsilon equation
    double yap_term = YapTerm(KappaAndEpsilon, ydist);
    SourceTerms[1] += yap_term;
  }

  return;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
