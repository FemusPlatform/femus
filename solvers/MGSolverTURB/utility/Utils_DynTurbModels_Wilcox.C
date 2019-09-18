#include <iostream>
#include "Utils_DynTurbModels_Wilcox.h"

double DynWilcox::CalcMuT(double* TurbVariables, double ydist) {
  CalcKappaAndOmega(TurbVariables);
  double MuTurb = _KappaAndOmega[0] / (_nu * _KappaAndOmega[1]);
  return MuTurb;
}

void DynWilcox_KW::CalcKappaAndOmega(double* TurbVariables) {
  _KappaAndOmega[0] = (TurbVariables[0] > _LowerKappa) ? TurbVariables[0] : _LowerKappa;
  _KappaAndOmega[1] = (TurbVariables[1] > _LowerOmega) ? TurbVariables[1] : _LowerOmega;
  return;
}

void DynWilcox_KW::ConvertKandWtoLocal(double* TurbValues) {
  TurbValues[0] = (_KappaAndOmega[0]);
  TurbValues[1] = (_KappaAndOmega[1]);

  return;
}

void DynWilcox_KW::CalcDynSourceTerms(
    double* TurbVariables, double* SourceTerms, double* DissTerms, double MuT, double VelGradMod,
    double ydist) {
  CalcKappaAndOmega(TurbVariables);
  double prod_k = 0.5 * _nu * MuT * VelGradMod;
  if (_Park == 1) { prod_k = Park_KSource(prod_k, _KappaAndOmega[0], VelGradMod); }

  SourceTerms[0] = prod_k;
  SourceTerms[1] = __AW * prod_k * _KappaAndOmega[1] / _KappaAndOmega[0];
  DissTerms[0] = __BETAS * _KappaAndOmega[1];
  DissTerms[1] = __BETAW * _KappaAndOmega[1];

  return;
}

void DynWilcox_logKW::CalcKappaAndOmega(double* TurbVariables) {
  double kappa = exp(TurbVariables[0]), omega = exp(TurbVariables[1]);
  _KappaAndOmega[0] = (kappa > _LowerKappa) ? kappa : _LowerKappa;
  _KappaAndOmega[1] = (omega > _LowerOmega) ? omega : _LowerOmega;
  return;
}

void DynWilcox_logKW::ConvertKandWtoLocal(double* TurbValues) {
  TurbValues[0] = log(_KappaAndOmega[0]);
  TurbValues[1] = log(_KappaAndOmega[1]);

  return;
}

void DynWilcox_logKW::CalcDynSourceTerms(
    double* TurbVariables, double* SourceTerms, double* DissTerms, double MuT, double VelGradMod,
    double ydist) {
  CalcKappaAndOmega(TurbVariables);
  double prod_k = 0.5 * _nu * MuT * VelGradMod;
  if (_Park == 1) { prod_k = Park_KSource(prod_k, _KappaAndOmega[0], VelGradMod); }
  SourceTerms[0] = prod_k / _KappaAndOmega[0];
  SourceTerms[1] = __AW * prod_k / _KappaAndOmega[0];
  DissTerms[0] = __BETAS * _KappaAndOmega[1];
  DissTerms[1] = __BETAW * _KappaAndOmega[1];

  return;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
