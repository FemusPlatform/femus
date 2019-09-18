// =======================================================================================
//            INITIAL AND BOUNDARY CONDITION FOR THERMAL TURBULENCE EQUATIONS
// =======================================================================================
#include "Equations_conf.h"

#ifdef RANS_THERMAL_EQUATIONS

#include "MGMesh.h"
#include "MGSolverRANS_thermal.h"
#include "MGUtils.h"
#include "User_RANS_thermal.h"

// INITIAL CONDITION ---------------------------------------------------------------------
void MGSolRANS_thermal::ic_read(int bc_gam, int bc_mat, double xp[], int iel, double u_value[]) {
  double ref_DT;
  double kh_value, wh_value, kh, wh, vel, diameter, wall_dist;
  vel = _mgutils._TurbParameters->_vmid;
  diameter = _mgutils._TurbParameters->_diameter;
  wall_dist = _WallDist + _mgmesh._dist[iel];

  ref_DT = _qs * diameter / (_kappa0);

  _mgutils._TurbParameters->TherTurInitValues(kh, wh, wall_dist, ref_DT, vel, diameter, _FlatProfile);
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
