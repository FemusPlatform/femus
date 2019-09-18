#include "Equations_conf.h"

// ============================================
#ifdef IMMERSED_BOUNDARY  // 3D-2D Energy equation
// ============================================
#include "MGFE.h"        // Mesh class,double vel_g[]
#include "MGSolverIB.h"  // Navier-Stokes class header file

// ===============================================================================================
void MGSolIB::vol_integral(
    DenseMatrixM& KeM, DenseVectorM& FeM, const int el_ndof2, const int el_ngauss, double xx_qnds[],
    const int unsteady, const int mode, int el_conn[], double u_old[],
    int iel) {  // ==============================================================================================
  double xyz_g[DIMENSION], dc_gdx[DIMENSION];
  // --------------------------------------------------------------------------------------------------------------------
  /// c) gaussian integration loop (n_gauss)
  // ------------------------------------------------------------------------------------------------------------------

  double alpha = _mgutils._mat_prop["ColDiffusion"];
  if (_dir == 0) alpha = 0.;

  for (int qp = 0; qp < el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2 = _fe[2]->Jac(qp, xx_qnds, _InvJac2);         // Jacobian
    double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];  // weight
    _fe[2]->get_phi_gl_g(_nTdim, qp, _phi_g[2]);              // shape funct
    _fe[2]->get_dphi_gl_g(_nTdim, qp, _InvJac2, _dphi_g[2]);  // global coord deriv
    _fe[1]->get_dphi_gl_g(_nTdim, qp, _InvJac2, _dphi_g[1]);  // global coord deriv
#ifdef AXISYM                                                 // axisymmetric (index -> 0)
    interp_el_sol(_xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g);
    JxW_g[2] *= xyz_g[0];
#endif
    /// d) Local (element) assemblying energy equation
    // =====================================================================================================================
    for (int i = 0; i < el_ndof2; i++) {
      if (_dir == 1) {
        //                 if ( _y_dist > _mgutils._IBParameter->_LowerThreshold ) {
        if (_bc_el[i] == 1) {
          double dtxJxW_g = JxW_g2;  // area with bc and weight
          const double phii_g = _phi_g[2][i];
          // Matrix Assemblying
          // ------------------------------------------------------------------------------------------------
          FeM(i) += dtxJxW_g * (_y_dist)*phii_g;
          for (int j = 0; j < el_ndof2; j++) {
            double Lap = 0.;
            for (int idim = 0; idim < _nTdim; idim++)
              Lap += _dphi_g[2][j + idim * el_ndof2] * _dphi_g[2][i + idim * el_ndof2];
            const double phij_g = _phi_g[2][j];
            //                             FeM ( i )   += dtxJxW_g* ( ( 0.5*_dir + ( 1-_dir ) )
            //                             *u_old[j]*phij_g ) *phii_g;
            KeM(i, j) += dtxJxW_g * (                    // energy-equation matrix
                                        phij_g * phii_g  // time term
                                        + Lap * alpha    // diffusion term
                                    );
          }
        } else {
          FeM(i) = 0.;
          KeM(i, i) = 1.;
        }
        //                 } // ----------------------------------------
        //                 else {
        //                     FeM ( i ) = 0.;
        //                     KeM ( i,i ) = 1000.;
        //                 }
      } else {
        //                 if ( _bc_el[i]==1 ) {
        double dtxJxW_g = JxW_g2;  // area with bc and weight
        const double phii_g = _phi_g[2][i];
        for (int j = 0; j < el_ndof2; j++) {
          double Lap = 0.;
          const double phij_g = _phi_g[2][j];
          FeM(i) += dtxJxW_g * (u_old[j] * phij_g) * phii_g;
          KeM(i, j) += dtxJxW_g * (phij_g * phii_g);
        }
        //                 }
        //                 else {
        //                     FeM ( i ) = 0.;
        //                     KeM ( i,i ) = 1.;
        //                 }
      }
    }
  }  // end of the quadrature point qp-loop ***********************
  return;
}

#endif
