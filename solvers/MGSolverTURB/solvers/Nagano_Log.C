#include "Equations_conf.h"

// ============================================
#ifdef RANS_EQUATIONS  // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGUtils.h"
#include "Nagano_Log.h"

MGSolNaganoLog::MGSolNaganoLog(
    MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
    const int nvars_in[],                    ///< KLQ number of variables
    std::string eqname_in,                   ///< equation name
    std::string varname_in                   ///< basic variable name
    )
    : MGSolRANS(
          mg_equations_map_in,  ///<  mg_equations_map_in pointer
          nvars_in,             ///< KLQ number of variables
          eqname_in,            ///< equation name
          varname_in            ///< basic variable name
          )                     // parameter  conductivity reference
{
  //  =========================================================================

  _var_names[0] = varname_in;

  if (!varname_in.compare("lk")) {
    _dir = 0;  // kappa
    _refvalue[0] = _uref * _uref;
  }

  if (!varname_in.compare("lw")) {
    _dir = 1;  // omega
    _refvalue[0] = _uref * _uref * _uref / _lref;
  }

  _InvSigma = 1. / 1.4;
  _ExplicitNearWallDer[0] = _ExplicitNearWallDer[1] = 1;

  double k_lower_lim = _mgutils._TurbParameters->GetKlim();
  double w_lower_lim = _mgutils._TurbParameters->GetWlim();

  _TurLowerLim[0] = log(k_lower_lim);
  _TurLowerLim[1] = log(w_lower_lim);

  return;
}

void MGSolNaganoLog::CalcAdvectiveAndDiffusiveTerms(int i, int j, int el_ndof2, double f_upwind) {
  _Adv = _Lap = _LapSupg = _LapMuTurb = _Cross[0] = _Cross[1] = _Log_Cross[0] = _Log_Cross[1] = 0.;

  for (int idim = 0; idim < _nTdim; idim++) {  // LOOP OVER SPACE COMPONENTS
    const double dphiidxg = _dphi_g[2][i + idim * el_ndof2];
    const double dphijdxg = _dphi_g[2][j + idim * el_ndof2];

    _Adv += _Vel_g[idim] * dphijdxg;       // advection
    _Lap += _nueff * dphijdxg * dphiidxg;  // diffusion
    _Cross[1] += 2. * _nueff * _T_dxg[0][idim] * dphijdxg;

    _Log_Cross[0] += _nueff * (_T_dxg[0][idim]) * dphijdxg;
    _Log_Cross[1] += _nueff * (_T_dxg[1][idim]) * dphijdxg;

    if (_UPWIND > 0.001) {
      _Lap += _UPWIND * f_upwind * _Vel_g[idim] * _Vel_g[idim] * dphijdxg * dphiidxg;  // normal upwind
    }

    if (_SUPG && _mgutils._TurbParameters->_SolveMuT == 1 && _WallElement != 1) {
      _LapMuTurb += dphijdxg * _MuTurbDxg[idim];
    }

    if (_SUPG) { _LapSupg += _nueff * _ddphi_g[2][j * _nTdim * _nTdim + idim * _nTdim + idim]; }
  }

  return;
}

void MGSolNaganoLog::CalcSourceAndDiss(int el_ndof2) {
  _y_dist = _ub_g[2][_FF_idx[DIST]];

  // function must be called using _kappa_g and not values from non linear iterations
  double mut = max(0., _ub_g[2][_FF_idx[MU_T]]);
  _mgutils._TurbParameters->CalcDynTurSourceAndDiss(
      _kappa_g, _y_dist, _sP, mut, _source, _diss);  // point wise source and diss
  _mu_turb = max(0., _ub_g[2][_FF_idx[MU_T]]);

  double kappa = exp(_kappa_g[0]);
  const double k_lim = _mgutils._TurbParameters->GetKlim();
  kappa = (kappa < k_lim) ? k_lim : kappa;
  double MuDurbin = kappa / (sqrt(_sP + 1.e-10) * _IRe);  // Durbin limit value

  // setting explicit and implicit source and dissipation
  _explicit_source[_dir] = _source[_dir];
  _explicit_diss[0] = _diss[0];
  _explicit_diss[1] = fabs(_diss[1]);
  _implicit_diss[_dir] = 0.;
  _implicit_source[_dir] = 0.;

  // Durbin correction -> mu turb and explicit source
  if (_mgutils._TurbParameters->_Durbin == 1) {
    _explicit_source[_dir] /= (_mu_turb + 1.e-10);
    _mu_turb = (MuDurbin < _mu_turb) ? MuDurbin : _mu_turb;
    _explicit_source[_dir] *= _mu_turb;
  }

  _nueff = _IRe * (1. + _InvSigma * _mu_turb);

  return;
}

void MGSolNaganoLog::VelocityForSUPG(
    double& mod2_vel, double vel_g[],
    double VEL[])  // NUMERICAL STABILIZATION - UPWIND AND SUPG ===============================
{
  mod2_vel = 0.;
  for (int i = 0; i < DIMENSION; i++) {
    VEL[i] = vel_g[i];
    if (_ModifiedSupg) {
      VEL[i] -= _nueff * ((1 - _dir) * _T_dxg[0][i] + _dir * (_T_dxg[1][i] + 2. * _T_dxg[0][i]));
    }
    mod2_vel += VEL[i] * VEL[i];
  }
  mod2_vel = sqrt(mod2_vel);

  return;
}  //==========================================================================

#endif
// kate: indent-mode cstyle; indent-width 5; replace-tabs on;
