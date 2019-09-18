#include "Equations_conf.h"

// ============================================
#ifdef RANS_THERMAL_EQUATIONS  // 3D-2D Energy equation
// ============================================

#include "MGFE.h"                  // Mesh class
#include "MGSolverRANS_thermal.h"  // Navier-Stokes class header file

/// This function sets the _bc_el with boundary condition flags.
/// This function also assembles the surface integral to obtain the algebraic sytem.
/// It is called by MGSolTBK::GenMatRhs in MGSolverT.C
/// This function sets  the  functional defined (user write}
void MGSolRANS_thermal::bc_set(int sur_toply[], int el_ndof2, int elb_ndof2, int elb_ngauss) {
  double det = _fe[2]->JacSur(elb_ngauss - 1, _xxb_qnds, _InvJac2);  // jacobian
  double Ipenalty = 1.;                                              // Dirichlet bc flag
  double utau, normal[DIMENSION], yplus, vel_mod, x_m[DIMENSION];
  int sign;
  int bc_face = _bc_vol[sur_toply[NDOF_FEMB - 1]];

  {
    // CALCULATION OF UTAU AND UNIT NORMAL VECTOR
    double vel_norm, vel_bound;

    for (int idim = 0; idim < _nTKdim; idim++) {  // based on scalar product between v and each element side
      x_m[idim] = 0.;

      for (int d = 0; d < NDOF_FEM; d++) { x_m[idim] += _xx_qnds[idim * NDOF_FEM + d] / NDOF_FEM; }
    }

    _fe[2]->normal_g(_xxb_qnds, x_m, normal, sign);

    vel_mod = vel_norm = 0.;

    for (int dim = 0; dim < _nTKdim; dim++) {
      double u_dim = _data_eq[2].ub[(_FF_idx[NS_F] + dim) * NDOF_FEM + sur_toply[NDOF_FEMB - 1]];
      vel_mod += u_dim * u_dim;
      vel_norm += u_dim * normal[dim];
    }

    vel_bound = sqrt(vel_mod - vel_norm * vel_norm);
    utau = _mgutils._TurbParameters->CalcUtau(vel_bound, _WallDist);
    yplus = _WallDist * utau / _IRe;
    vel_mod = sqrt(vel_mod);
  }

  if (_RANS_t_parameter._WallFunctionApproach == 1 && bc_face != 10) {
    if (bc_face == 11 || bc_face == 12) { _WallElement = 0; }

    if (bc_face == 1) {
      _WallElement = 2;

      for (int dir = 0; dir < _nTKdim; dir++) { _NormMidCell[dir] += normal[dir]; }
    }

    if (bc_face == 21) {
      _WallElement = 5;

      for (int dir = 0; dir < _nTKdim; dir++) { _NormMidCell[dir] += normal[dir]; }
    }
  }

  if (_RANS_t_parameter._WallFunctionApproach == 0 || _WallElement == 0) {
    // DIRICHLET BOUNDARY CONDITIONS ============================================================
    for (int lb_node = 0; lb_node < elb_ndof2; lb_node++) {
      if (_bc_el[sur_toply[lb_node]] == 0) {
        int lv_node = sur_toply[lb_node];  // local vol index
        int bc_s = _bc_vol[lv_node] % 10;  // if bc_s == 0 then homogeneous dirichlet
        int bc_wall = bc_s % 2;
        int bc_inlet = (bc_s & 2) >> 1;

        double k_value = _data_eq[2].ub[(_FF_idx[KTT_F]) * NDOF_FEM + lv_node];
        double w_value = _data_eq[2].ub[(_FF_idx[KTT_F] + 1) * NDOF_FEM + lv_node];

        double wall_value[2];
        wall_value[0] = wall_value[1] = 1.;

        if (bc_wall) {  // DIRICHLET VALUES ON WALL BOUNDARIES
          double WallDistance = _WallDist;
          double kappa, omega;

          if (_RANS_t_parameter._WallFunctionApproach == 1) { WallDistance = _y_dist; }

          _mgutils._TurbParameters->ThermTurNearWallValues(kappa, omega, WallDistance, utau);
          wall_value[0] = kappa;
          wall_value[1] = omega;
        }

        double inlet_value[2];
        inlet_value[0] = inlet_value[1] = 1.;

        if (bc_inlet) {  // DIRICHLET VALUES ON INLET BOUNDARIES   ->
                         // http://support.esi-cfd.com/esi-users/turb_parameters/
          inlet_value[0] =
              _data_eq[2].ub[(_FF_idx[KTT_F]) * NDOF_FEM + lv_node];  // k_value;//(1.-_khlog)*k_in +
                                                                      // _khlog*log(k_in);
          inlet_value[1] =
              _data_eq[2].ub[(_FF_idx[KTT_F] + 1) * NDOF_FEM + lv_node];  // w_value;//(1.-_whlog)*w_in +
                                                                          // _whlog*log(w_in);
        }

        double Value = (bc_wall * wall_value[_dir] + bc_inlet * inlet_value[_dir]);
        double FemVal = (Value > 1.) ? 1 : Value;
        double KemVal = (Value > 1.) ? Value : 1;
        _FeM(lv_node) += Ipenalty * FemVal;
        _KeM(lv_node, lv_node) += Ipenalty / KemVal;
      }  // END DIRICHLET BOUNDARY CONDITION
    }    // END LOOP OVER BOUNDARY FACE NODES

    // NEUMANN BOUNDARY CONDITIONS ============================================================
    if (_bc_vol[sur_toply[NDOF_FEMB - 1]] / 10 > 0) {  // bc 10-30
      int bc_s = _bc_vol[sur_toply[NDOF_FEMB - 1]] % 10;
      int bc_alpha = (int)((bc_s & 2) >> 1);  // (1?) linear term
      int bc_beta = (int)(bc_s % 2);          // (?1)  constant term

      for (int qp = 0; qp < elb_ngauss; qp++) {                   // GAUSS LOOP
        double det = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);     // local coord _phi_g and jac
        double JxW_g2 = det * _fe[2]->_weight1[_nTKdim - 2][qp];  // weight
        _fe[2]->get_phi_gl_g(_nTKdim - 1, qp, _phi_g[2]);         // global coord _phi_g

        if (_AxiSym == 1) {  // axisymmetric  (index ->0)
          double xyg[DIMENSION];
          interp_el_bd_sol(_xx_qnds, sur_toply, elb_ndof2, 0, _nTKdim, _phi_g[2], elb_ndof2, xyg);
          JxW_g2 *= xyg[0];
        }

        for (int lsi_node = 0; lsi_node < elb_ndof2; lsi_node++) {  // local side loop (over the node face)
          const double phii_g = _phi_g[2][lsi_node];                // boundary test function
          const int lei_node = sur_toply[lsi_node];                 // local element index

          if (_bc_el[lei_node] != 0) {
            double dtxJxW_g = JxW_g2;
            double wall_der[2];
            wall_der[0] = wall_der[1] = 0.;

            if (bc_alpha) {  // NEAR WALL DERIVATIVES
              const double k_der = -2. * _alpha / (_y_dist);
              const double w_der = -2. * _alpha / (_y_dist);
              wall_der[1] = (yplus < 11.) ? -w_der : 0;
              wall_der[0] = (yplus < 11.) ? k_der : 0.;
            }

            const double FemBC = _ExplicitNearWallDer[_dir];
            const double KemBC = 1 - _ExplicitNearWallDer[_dir];

            if (FemBC == 1) { _FeM(lei_node) += dtxJxW_g * phii_g * (bc_alpha * FemBC * wall_der[_dir]); }

            // Assemblying Matrix ---------------------------------
            if (KemBC == 1) {
              for (int lsj_node = 0; lsj_node < elb_ndof2; lsj_node++) {
                _KeM(lei_node, sur_toply[lsj_node]) +=
                    dtxJxW_g * bc_alpha * KemBC * wall_der[_dir] * phii_g * _phi_g[2][lsj_node];
              }  // END LOOP OVER BOUNDARY NODES
            }    // END KEM NEUMANN BC
          }      // END IF BC_EL != 0
        }        // END LOOP OVER BOUNDARY TEST FUNCTIONS
      }          // END QUADRATURE LOOP - qp INDEX
    }  // END NEUMANN BOUNDARY CONDITIONS ===========================================================
  }

  return;
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
