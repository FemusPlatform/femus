#include <iostream>
#include "Utils_DynTurbModels.h"

double DYNturModels::YapTerm(double* KappaAndEpsilon, double ydist) {
  double YAP_source = 0.;

  double kappa = KappaAndEpsilon[0];
  double epsilon = KappaAndEpsilon[1];

  double l_epsilon = pow(0.09, 0.75) * 0.41 * ydist;
  double p1 = 0.83 * epsilon * epsilon / kappa;
  double p2 = (kappa * sqrt(kappa) / (epsilon * l_epsilon) - 1);
  double p3 = (kappa * sqrt(kappa) / (epsilon * l_epsilon));
  p3 *= p3;
  YAP_source = p1 * p2 * p3;
  YAP_source = max(YAP_source, 0.);

  return YAP_source;
}

void DYNturModels::DynTurNearKappaAndOmegaValues(double* TurbVariables, double WallDist, double Utau) {
  const double y_plus = WallDist * Utau / _nu;

  // OMEGA VALUES
  const double yplim = 2. * 0.41 / sqrt(0.09);
  const double wlin = 2. * _nu / (0.09 * WallDist * WallDist);
  const double wlog = Utau / (0.41 * WallDist * sqrt(0.09));
  _KappaAndOmega[1] = (y_plus > yplim) ? wlog : wlin;

  // MUTURB VALUES
  const double a = 0.41;
  const double c = 0.001093;
  double Mu = _nu / (1. / (c * y_plus * y_plus * y_plus) + 1. / (0.41 * y_plus));

  // KAPPA VALUES
  if (y_plus > 11.6) {
    const double ck = -0.416;
    const double bk = 8.366;
    _KappaAndOmega[0] = Utau * Utau * (ck / 0.41 * log(y_plus) + bk);
  } else {
    const double c = 11.0;
    const double cf = (1.0 / ((y_plus + c) * (y_plus + c)) + 2.0 * y_plus / (c * c * c) - 1. / (c * c));
    const double ceps2 = 1.9;
    _KappaAndOmega[0] = Utau * Utau * 2400. * cf / (ceps2 * ceps2);
  }

  //     kappa = omega * MU;
  ConvertKandWtoLocal(TurbVariables);
  return;
}

void DYNturModels::DynTurInitKappaAndOmegaValues(double* TurbVariables, double vel, double diameter) {
  const double Re_h = vel * diameter / _nu;
  const double In = 0.16 * pow(Re_h, -0.125);
  const double len = 0.07 * diameter;
  _KappaAndOmega[0] = 1.5 * (In * vel) * (In * vel);
  _KappaAndOmega[1] = sqrt(_KappaAndOmega[0]) / len;
  ConvertKandWtoLocal(TurbVariables);

  return;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
