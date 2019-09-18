#include "Equations_conf.h"

// ============================================
#ifdef RANS_THERMAL_EQUATIONS  // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "Nagano_KWT.h"
#include "Printinfo_conf.h"

#include <sstream>
#include "MGGeomEl.h"
#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"
#include "MGSystem.h"
#include "MGFE.h"
#include "MGUtils.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"
#include "parallelM.h"

MGSolNaganoKWT::MGSolNaganoKWT(
    MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
    const int nvars_in[],                    ///< KLQ number of variables
    std::string eqname_in,                   ///< equation name
    std::string varname_in                   ///< basic variable name
    )
    : MGSolRANS_thermal(
          mg_equations_map_in,  ///<  mg_equations_map_in pointer
          nvars_in,             ///< KLQ number of variables
          eqname_in,            ///< equation name
          varname_in            ///< basic variable name
          )                     // parameter  conductivity reference
{
  //  =========================================================================

  _var_names[0] = varname_in;

  if (!varname_in.compare("kh")) {
    _dir = 0;  // kappa
    _refvalue[0] = _uref * _uref;
  }

  if (!varname_in.compare("wh")) {
    _dir = 1;  // omega
    _refvalue[0] = _uref * _uref * _uref / _lref;
  }

  _ExplicitNearWallDer[0] = _ExplicitNearWallDer[1] = 0;

  double khlim = _mgutils._TurbParameters->GetKHlim();
  double whlim = _mgutils._TurbParameters->GetWHlim();

  _LowerLimit = (_dir == 0) ? khlim : whlim;

  return;
}

void MGSolNaganoKWT::CalcAdvectiveAndDiffusiveTerms(int i, int j, int el_ndof2, double f_upwind) {
  _Adv = _Lap = _LapSupg = _LapMuTurb = _Cross[0] = _Cross[1] = _Log_Cross[0] = _Log_Cross[1] = 0.;
  const double kappaT = _Tkappa_g[0];

  for (int idim = 0; idim < _nTKdim; idim++) {  // LOOP OVER SPACE COMPONENTS
    const double dphiidxg = _dphi_g[2][i + idim * el_ndof2];
    const double dphijdxg = _dphi_g[2][j + idim * el_ndof2];

    _Cross[1] += 2. * _alpha_eff * _KH_der[idim] * dphijdxg / kappaT;

    _Adv += _Vel_g[idim] * dphijdxg;           // advection
    _Lap += _alpha_eff * dphijdxg * dphiidxg;  // diffusion

    if (_UPWIND > 0.001) {
      _Lap += _UPWIND * f_upwind * _Vel_g[idim] * _Vel_g[idim] * dphijdxg * dphiidxg;  // normal upwind
    }

    if (_SUPG && _RANS_t_parameter._InterpolatedAlphaTurb == 1 && _WallElement != 1) {
      _LapMuTurb += dphijdxg * _AlphaTurbDxg[idim];
    }

    if (_SUPG) { _LapSupg += _alpha_eff * _ddphi_g[2][j * _nTKdim * _nTKdim + idim * _nTKdim + idim]; }
  }

  return;
}

void MGSolNaganoKWT::VelocityForSUPG(
    double& mod2_vel, double vel_g[],
    double VEL[])  // NUMERICAL STABILIZATION - UPWIND AND SUPG ===============================
{
  for (int i = 0; i < DIMENSION; i++) {
    VEL[i] = vel_g[i];

    if (_ModifiedSupg) { VEL[i] -= _alpha_eff * (_dir * 2. * _KH_der[i]); }

    mod2_vel += VEL[i] * VEL[i];
  }

  mod2_vel = sqrt(mod2_vel);

  return;
}  //==========================================================================

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
