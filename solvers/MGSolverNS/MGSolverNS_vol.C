// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include "MGFE.h"        // Mesh class
#include "MGFE_conf.h"   // FEM approximation
#include "MGSolverNS.h"  // Navier-Stokes class header file

void MGSolNS::matrixrhsvol(int el_ndof[], const int mode, int el_conn[]) {
  double dphijdx_g2[DIMENSION], dphiidx_g2[DIMENSION];
  double rho = 1.;
  const int el_ngauss = _fe[2]->_NoGauss1[_nNSdim - 1];  // elem gauss points

  for (int qp = 0; qp < el_ngauss; qp++) {
    const double det2 = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);   // quadratic Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nNSdim - 1][qp];  // quadratic weight

    GetVolumeTestFunctions(qp);
    InterpolateSolutions(el_ndof);

    if (_AxiSym == 1) JxW_g2 *= _xyzg[0];

    _IReEff = _IRe;

    if (_FF_idx[MU_T] >= 0) { _IReEff += _IRe * _ub_g[2][_FF_idx[MU_T]]; }

    _f_upwind = CalcFUpwind(_u_1ts_g, _dphi_g[2], _IReEff, _nNSdim, el_ndof[2]);

    // MOMENTUM EQUATION
    MomentumEquation(JxW_g2, el_ndof, qp);

    // CONTINUITY EQUATION
    ContinuityEquation(JxW_g2, el_ndof);

  }  // end of the quadrature point qp-loop

  // ====================== end volume (element) =======================================
  return;
}

void MGSolNS::CalcTimeDer(double& TimeDer, TIME_DER_TYPE Type, int SpaceDir, int nPhi, int interpNodeID) {
  TimeDer = 0.;

  if (Type == ACTUAL_STEP) {
    const int TD_coeff = (_NS_parameter._TimeDisc == 2) ? 1.5 : 1.;

    if (_NS_parameter._SolveSteady == 0)
      TimeDer = TD_coeff * (_phi_g[2][nPhi] + _Phi_supg) * _phi_g[2][interpNodeID] / _dt;
  }

  if (Type == OLDER_STEPS) {
    double TD_1ts = (_NS_parameter._TimeDisc == 2) ? 2. : 1.;
    double TD_2ts = (_NS_parameter._TimeDisc == 2) ? 0.5 : 0.;

    if (_NS_parameter._SolveSteady == 0)
      TimeDer =
          (TD_1ts * _u_1ts_g[SpaceDir] - TD_2ts * _u_2ts_g[SpaceDir]) * (_phi_g[2][nPhi] + _Phi_supg) / (_dt);
  }

  return;
}

void MGSolNS::CalcAdv_Lap_LapSupg(
    double& Adv, double& Lap, double& LapSupg, int nPhi, int j, int el_ndof[], int dir) {
  double dphijdx_g2[DIMENSION];

  double ADV_1ts = (_NS_parameter._TimeDisc == 2) ? 2. : 1.;
  double ADV_2ts = (_NS_parameter._TimeDisc == 2) ? 1. : 0.;

  Adv = Lap = LapSupg = 0.;

  if (_AxiSym == 1 && dir == 0) {  // axysimmetry only --------------------------------------------------
    Lap = 2. * (_IReEff + _NS_parameter._Upwind * _f_upwind * _u_nl_g[0] * _u_nl_g[0]) * _phi_g[2][j] *
          (_phi_g[2][nPhi] + _Phi_supg) / (_xyzg[0] * _xyzg[0]);
  }  // --------------------------------------------------------------------------------

  for (int kdim = 0; kdim < _nNSdim; kdim++) {
    dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof[2]];
    double up = _NS_parameter._Upwind * _f_upwind * _u_1ts_g[kdim] * _u_1ts_g[kdim];
    Adv += (ADV_1ts * _u_1ts_g[kdim] - ADV_2ts * _u_2ts_g[kdim]) *
           dphijdx_g2[kdim];  // Adv_g +=_u_1ts_g[kdim]*dphijdx_g2[kdim]*phii_g;
    Lap += (_IReEff + up) * dphijdx_g2[kdim] * _dphi_g[2][nPhi + kdim * el_ndof[2]];
    if (_NS_parameter._Supg == 1) {
      LapSupg += (_IReEff)*_ddphi_g[2][j * _nNSdim * _nNSdim + kdim * _nNSdim + kdim];
    }
  }

  Adv *= (_phi_g[2][nPhi] + _Phi_supg);
  LapSupg* _Phi_supg;

  return;
}

void MGSolNS::CalcPresGrad(double& PresGrad, int nPhi, int interpNodeID, int spaceDir, int el_ndof[]) {
  PresGrad = 0.;

  if (_AxiSym == 1 && spaceDir == 0) {
    PresGrad = -1. * _phi_g[2][nPhi] * _phi_g[1][interpNodeID] / _xyzg[0];
  }

  PresGrad +=
      (-1. * _phi_g[1][interpNodeID] * _dphi_g[2][nPhi + spaceDir * el_ndof[2]] +
       _dphi_g[1][interpNodeID + spaceDir * el_ndof[1]] * (0. * _phi_g[2][nPhi] + _Phi_supg));
  return;
}

void MGSolNS::CalcVelDivergence(
    double& Div, int PresOrder, int nPhi, int InterpNodeID, int SpaceDir, int el_ndof[]) {
  Div = 0.;
  const double psii_g = _phi_g[PresOrder][nPhi];

  if (_AxiSym == 1 && _nNSdim == 2) {
    const double phij_g = _phi_g[2][InterpNodeID];
    Div += (1 - SpaceDir) * psii_g * phij_g / _xyzg[0];
  }

  Div += psii_g * _dphi_g[2][InterpNodeID + SpaceDir * el_ndof[2]];

  return;
}

void MGSolNS::GetVolumeTestFunctions(int qp) {
  _fe[0]->get_phi_gl_g(_nNSdim, qp, _phi_g[0]);
  _fe[1]->get_phi_gl_g(_nNSdim, qp, _phi_g[1]);  // linear shape funct
  _fe[2]->get_phi_gl_g(_nNSdim, qp, _phi_g[2]);  // quadratic shape function

  _fe[1]->get_dphi_gl_g(_nNSdim, qp, _InvJac2, _dphi_g[1]);  // global coord deriv
  _fe[2]->get_dphi_gl_g(_nNSdim, qp, _InvJac2, _dphi_g[2]);  // global coord deriv

  _fe[2]->get_ddphi_gl_g(_nNSdim, qp, _InvJac2, _ddphi_g[2]);  // global second derivatives

  return;
}

void MGSolNS::RowSetUp(int nPhi, int indx_eq, int qp, int el_ndof[]) {
  _Phi_supg = 0.;

  if (_bc_el[indx_eq] == 0 && _NS_parameter._Supg == 1) {
    for (int idim = 0; idim < _nNSdim; idim++) {  // Phi_supg = 0. on boundaries with Dirichlet bc
      _Phi_supg += _f_upwind * _u_1ts_g[idim] * _dphi_g[2][nPhi + idim * el_ndof[2]];  // f_upwind
    }
  }

  if (_bc_el[indx_eq] == -8 && qp == 0) {
    std::cout << "ATTENTION!!! el node " << nPhi << " row " << indx_eq << " bc_el " << _bc_el[indx_eq]
              << " bc_vol " << _bc_vol[nPhi] << " bc_bd " << _bc_bd[nPhi] << std::endl;
    _bc_el[indx_eq] = 1;
  }

  _BoundEquation = (_bc_el[indx_eq] <= -1) ? 1 : 0;
  _NormalEquation = 0;

  if (_BoundEquation && _bc_el[indx_eq] == -3) {
    _NormalEquation = 1;

    for (int dim = 0; dim < _nNSdim; dim++) { _ProjDir[2][dim] = _normal_pt[dim + nPhi * _nNSdim]; }
  }

  if (_BoundEquation && !_NormalEquation) {
    CalcTangDir(nPhi, el_ndof[2], _bc_el[indx_eq], qp, _ProjDir[0], _ProjDir[1], VEL_BASED);
  }

  _NComp =
      (_BoundEquation) ? _nNSdim : 1;  // Default: 1 vel comp for equation in row indx -> ivar = row_shift

  _ImmVal = 1;
  if (_ImmersedBoundary == 1)
    if (_data_eq[2].ub[(_FF_idx[IB_COL]) * NDOF_FEM + nPhi] > 0.9) _ImmVal = 0;

  return;
}

void MGSolNS::InterpolateSolutions(int el_ndof[]) {
  interp_el_sol(_xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof[2], _xyzg);
  interp_el_sol(
      _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof[2],
      _ub_g[2]);                                                        // field _ub_g[2][DIM]
  interp_el_sol(_u_2ts, 0, _nNSdim, _phi_g[2], el_ndof[2], _u_2ts_g);   // field _ub_g[2][DIM]
  interp_el_sol(_u_nl, 0, _nNSdim, _phi_g[2], el_ndof[2], _u_nl_g);     // field _ub_g[2][DIM]
  interp_el_gdx(_u_1ts, 0, _nNSdim, _dphi_g[2], el_ndof[2], _vel_gdx);  // derivatives  vel_gdx[DIM][DIM]
  interp_el_gdx(_Pressure, 0, 1, _dphi_g[1], el_ndof[1], _P_gdx);       // derivatives  vel_gdx[DIM][DIM]
  interp_el_gddx(_u_1ts, 0, _nNSdim, _ddphi_g[2], el_ndof[2], _vel_gddx);

  for (int idim = 0; idim < _nNSdim; idim++) { _u_1ts_g[idim] = _ub_g[2][_FF_idx[NS_F] + idim]; }

  return;
}

void MGSolNS::CalcVolume() {
  const int el_ngauss = _fe[2]->_NoGauss1[_nNSdim - 1];  // elem gauss points
  _ElemVolume = 0.;
  for (int qp = 0; qp < el_ngauss; qp++) {
    const double det2 = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);   // quadratic Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nNSdim - 1][qp];  // quadratic weight
    _ElemVolume += JxW_g2;
  }

  return;
}

#endif
