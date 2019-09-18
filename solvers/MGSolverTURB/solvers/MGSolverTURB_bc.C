#include "Equations_conf.h"

// ============================================
#ifdef _TURBULENCE_  // 3D-2D Energy equation
// ============================================

#include "MGSolverTURB.h"  // Navier-Stokes class header file
#include "MGUtils.h"
#include "MGFE.h"  // Mesh class

// =========================================================================================
/// This function sets the _bc_el with boundary condition flags.
/// This function also assembles the surface integral to obtain the algebraic sytem.
/// It is called by MGSolTURB::GenMatRhs in MGSolverT.C
/// This function sets  the  functional defined (user write}
void MGSolTURB::bc_set(
    DenseMatrixM& KeM, DenseVectorM& FeM, int sur_toply[], int el_ndof2, int elb_ndof2, int elb_ngauss,
    int sign_normal)  // ====================================================================================
{
  double xyz_g[DIMENSION];
  //   Dirichlet  boundary conditions ---------------------------------------------------->
  double det = _fe[2]->JacSur(elb_ngauss - 1, _xxb_qnds, _InvJac2);  // jacobian
  double Ipenalty = det / _dt;                                       // Dirichlet bc flag
  // local boundary loop   ---------------------------------------
  //   for(int lb_node=0; lb_node< elb_ndof2; lb_node++)  if(_bc_el[sur_toply[lb_node]] ==0) {
  //       int lv_node= sur_toply[lb_node]; // local vol index
  //       int bc_s=_bc_vol[lv_node]%10;
  //       // Assemblying  Matrix & rhs
  //       FeM(lv_node) += bc_s*Ipenalty*_data_eq[2].ub[_FF_idx[CO_F]*NDOF_FEM+lv_node];
  //       KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
  // //       std::cout << _data_eq[2].ub[_FF_idx[CO_F]*NDOF_FEM+lv_node]<<    "  ";
  //     }// lb_node -end  local boundary loop
  //   // < ------------------------------------------------- end Dirichlet  boundary conditions

  // ========================  Neumann  boundary conditions  ==============================
  if (_bc_el[sur_toply[NDOF_FEMB - 1]] == 0 && _dir != 0) {  // bc 10-30

    for (int qp = 0; qp < elb_ngauss; qp++) {
      // quad/linear  [2]=quad [1]=linear------------------------------------
      double det = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);    // local coord _phi_g and jac
      double JxW_g2 = det * _fe[2]->_weight1[_nTdim - 2][qp];  // weight
      _fe[2]->get_phi_gl_g(_nTdim - 1, qp, _phi_g[2]);         // global coord _phi_g

      interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof2, 0, _nTdim, _phi_g[2], elb_ndof2, xyz_g);
#ifdef AXISYM  // axisymmetric  (index ->0)
               //          interp_el_bd_sol(_xx_qnds,sur_toply,elb_ndof2,0, _nTdim,_phi_g[2],elb_ndof2,xyz_g);
      JxW_g2 *= xyz_g[0];
#endif

      for (int lsi_node = 0; lsi_node < elb_ndof2; lsi_node++) {  // local side loop (over the node face)
        // set up row i
        const double phii_g = _phi_g[2][lsi_node];  // boundary test function
        const int lei_node = sur_toply[lsi_node];   // local element index
        double dtxJxW_g = JxW_g2;                   // Neumann bc flag and Jac

        // Assemblying Matrix ---------------------------------
        for (int lsj_node = 0; lsj_node < elb_ndof2; lsj_node++) {
          KeM(lei_node, sur_toply[lsj_node]) += dtxJxW_g * _DiffusionCoefficient * phii_g *
                                                _phi_g[2][lsj_node] * 3. /
                                                _y_dist;  // Robin bc  (k*dt/dn = h*(-T))
        }                                                 // end lsj_node
      }                                                   // lsi_node
    }  // end of the quadrature point qp-loop **********************
  }    // ======================== end Neumann  boundary conditions  ===========================

  return;
}

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
