#include "Equations_conf.h"  // <--- Equations configure

#ifdef RANS_THERMAL_EQUATIONS

#include "User_RANS_thermal.h"
#include "MGSolverRANS_thermal.h"
#include "MGUtils.h"
#include "MGMesh.h"

// INITIAL CONDITION ---------------------------------------------------------------------
void MGSolRANS_thermal::ic_read(int bc_gam, int bc_mat, double xp[], int iel, double u_value[]) {
  double kh_value, wh_value, kh, wh, vel, diameter, wall_dist, NormFlux;
  vel = 0.25456231405;  // stod(_mgutils.get_file("AvVelocity"));
  diameter = 0.03025;   // stod(_mgutils.get_file("Diameter"));
  wall_dist = _WallDist + _mgmesh._dist[iel];
  NormFlux = _qs / (_rhof * _cp0);
  double FluxCorr = 1.;

  if (_mgutils.get_name() == 1) {
    FluxCorr = 0.0041 * 0.2 * 0.04 / (4. * 0.01 * 10. * 3.1415 * (0.03025 * 0.03025 - 0.0041 * 0.0041));
  }

  _mgutils._TurbParameters->TherTurInitValues(
      kh, wh, wall_dist, FluxCorr * NormFlux, vel, diameter, _FlatProfile);
  kh = 0.0001;
  wh_value = (wh);
  kh_value = (kh);

  if (_dir == 0) { u_value[0] = kh_value; }

  if (_dir == 1) { u_value[0] = wh_value; }

  return;
}

// BOUNDARY CONDITIONS -------------------------------------------------------------------
void MGSolRANS_thermal::bc_read(int bc_gam, int bc_mat, double xp[], int bc_Neum[], int bc_flag[]) {
  double ILref = 1. / _lref;
  bc_Neum[0] = TKinsulation;

  if (_dir == 0) { bc_Neum[0] = _RANS_t_parameter._map_TTKgroup[bc_gam]; }

  if (_dir == 1) { bc_Neum[0] = _RANS_t_parameter._map_TTWgroup[bc_gam]; }

  return;
}

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
