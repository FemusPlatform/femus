// // ===============================================================
// --------------   NAVIER-STOKES system [DS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef DS_EQUATIONS

// ==============================================================
// DS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// DS_EQUATIONS==1 coupled    solver (u,v,w,p)
// DS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"  // Navier-Stokes class conf file
#include "MGSolverDS.h"     // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"  // FEM approximation

// local Femus class include -----------------------------------
#include "MGFE.h"  // Mesh class
#include "MeshExtended.h"

// Thermodinamical Properties
// ========================================================================================== constant
double _DSdensity(double T, double T_r) {
  return 1.;  // water (t K) T_ref
}
double _DSkviscosity(double T, double T_r) {
  return 1.;  // lead T K
}
//  ======================================================================================================================

// ==============================================================================================
void MGSolDS::matrixrhsvol_sol_ds(DenseMatrixM& KeM, DenseVectorM& FeM, int el_ndof[]) {
  for (int i = 0; i < el_ndof[2]; i++) {
    // set up row i
    for (int idim = 0; idim < _nDSdim; idim++) { _dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof[2]]; }
    //             double dtxJxW_g=JxW_g2*(_bc_el[i]);
    const double ds_old = _data_eq[2].ub[(_FF_idx[SDSX_F] + _dir) * NDOF_FEM + i];
    const double v_old = _data_eq[2].ub[(_FF_idx[FS_F] + _dir) * NDOF_FEM + i];

    // rhs  Assemblying ---------------------------
    FeM(i) = (ds_old + v_old * _dt);

    //       FeM(i) = 1.*_dir*(0.25-_xx_qnds[i])*(0.25-_xx_qnds[i]);
    // Matrix Assemblying ---------------------------
    KeM(i, i) = 1.;
  }  // ----------------------------------------
  // ====================== end volume (element) =======================================
  return;
}

// ==============================================================================================
void MGSolDS::matrixrhsvol_liq_ds(
    DenseMatrixM& KeM, DenseVectorM& FeM, int el_ndof[], int flag_group[], double u_old[],
    double myDisp[]) {  // ==================================  Volume
                        // ===============================================
  double l_gdx[DIMENSION * DIMENSION];
  double dist = 0.;
  double dist_x = 0.;
  double dist_y = 0.;
  const int el_ngauss = _fe[2]->_NoGauss1[_nDSdim - 1];  // elem gauss points
  for (int qp = 0; qp < el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2 = _fe[2]->Jac(qp, _xx_qnds, _InvJac2);               // Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nDSdim - 1][qp];        // weight
    _fe[2]->get_phi_gl_g(_nDSdim, qp, _phi_g[2]);                    // shape funct
    _fe[2]->get_dphi_gl_g(_nDSdim, qp, _InvJac2, _dphi_g[2]);        // global coord deriv
    interp_el_gdx(u_old, 0, _nDSdim, _dphi_g[2], el_ngauss, l_gdx);  // derivatives  vel_gdx[DIM][DIM]
    interp_el_sol(_xx_qnds, 0, _nDSdim, _phi_g[2], el_ngauss, _xyzg);

    //     double J=0.;
    //     J=(1.+l_gdx[0])*(1.+l_gdx[3])-l_gdx[1]*l_gdx[2];

    /// d) Local (element) assemblying energy equation
    // *********************** *******************************

    for (int i = 0; i < el_ndof[2]; i++) {
      for (int idim = 0; idim < _nDSdim; idim++) { _dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof[2]]; }
      int fl_int = ((fabs(flag_group[i]) > 999) ? 0 : 1);
      double dtxJxW_g = JxW_g2 * _bc_el[i] * fl_int;
      const double phii_g = _phi_g[2][i];
      dist_x = _xyzg[0] - (0.6 + myDisp[0]);
      dist_y = _xyzg[1] - (0.2 + myDisp[1]);
      dist = sqrt(dist_x * dist_x + dist_y * dist_y);
      //       cout<<dist<<endl;
      if (_bc_vol[i] == 3) dtxJxW_g = 0.;
      // Matrix Assemblying ---------------------------
      for (int j = 0; j < el_ndof[2]; j++) {
        double Lap = 0.;
        // set up test function
        const double phij_g = _phi_g[2][j];
        for (int idim = 0; idim < _nDSdim; idim++) {
          _dphijdx_g2[idim] = _dphi_g[2][j + idim * el_ndof[2]];
          Lap += 1. / (1. + 100000. * dist) * _dphijdx_g2[idim] * _dphiidx_g2[idim];  // Laplacian
        }
        // energy-equation
        KeM(i, j) += dtxJxW_g * (Lap);
      }
    }  // ----------------------------------------
  }    // end of the quadrature point qp-loop ***********************
  return;
}

#endif
