#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverT.h"       // Navier-Stokes class header file
#include "MGFE.h"            // Mesh class


void  MGSolT::bc_set (
    int sur_toply[],
    int el_ndof2,int elb_ndof2,int elb_ngauss,
    int sign_normal
) { // ====================================================================================
    double xyz_g[DIMENSION];
    //   Dirichlet  boundary conditions ---------------------------------------------------->
    double det= _fe[2]->JacSur ( elb_ngauss-1,_xxb_qnds, _InvJac2 ); // jacobian
    double Ipenalty=det;                               // Dirichlet bc flag
    // local boundary loop   ---------------------------------------
    for ( int lb_node=0; lb_node< elb_ndof2; lb_node++ )
        if ( _bc_el[sur_toply[lb_node]] ==0 ) {
            int lv_node= sur_toply[lb_node]; // local vol index
            int bc_s=_bc_vol[lv_node]%10;
            // Assemblying  Matrix & rhs
            _FeM ( lv_node ) += bc_s*Ipenalty*_T_1ts[lv_node];
            _KeM ( lv_node,lv_node ) += Ipenalty; //  Dirichlet bc
        }// lb_node -end  local boundary loop

    // ========================  Neumann  boundary conditions  ==============================
    if ( _bc_vol[sur_toply[NDOF_FEMB-1]]/10>0 ) { // bc 10-30

        // flag setup: +++++++++++++++++++++++++++++++++++++++
        //   (Neumann DT.n=bc_alpha*T+bc_beta*value)
        double alpha_eff=1.;
        const int  bc_s=_bc_vol[sur_toply[NDOF_FEMB-1]]%10;
        const int  bc_alpha = ( int ) ( ( bc_s&2 ) >>1 ); // (1?) linear term
        const int  bc_beta = ( int ) ( bc_s%2 );  // (?1)  constant term
        // Non homogenous Neumann boundary conditions  ***********************************
        // gaussian integration loop (n_gauss)
        // -----------------------------------------------
        for ( int qp=0; qp<  elb_ngauss; qp++ ) {
            // quad/linear  [2]=quad [1]=linear------------------------------------
            double det  = _fe[2]->JacSur ( qp,_xxb_qnds,_InvJac2 ); // local coord _phi_g and jac
            double JxW_g2=det*_fe[2]->_weight1[ _nTdim-2][qp];// weight
            _fe[2]->get_phi_gl_g ( _nTdim-1,qp,_phi_g[2] ); // global coord _phi_g

            // axisymmetric  (index ->0)
            if ( _Axisym==1 )  {
                interp_el_bd_sol ( _xx_qnds,sur_toply,el_ndof2,0, _nTdim,_phi_g[2],elb_ndof2,xyz_g );
                JxW_g2  *=xyz_g[0];
            }

            for ( int lsi_node=0; lsi_node< elb_ndof2; lsi_node++ ) { // local side loop (over the node face)
                // set up row i
                const double phii_g=_phi_g[2][lsi_node]; // boundary test function
                const int lei_node= sur_toply[lsi_node]; // local element index
                double dtxJxW_g=JxW_g2*_bc_el[lei_node]; // Neumann bc flag and Jac
                // Assemblying rhs ----------------------------
                _FeM ( lei_node ) += dtxJxW_g*bc_beta* (
                                        + ( 1-bc_alpha ) *_qs*phii_g   // Neumann heat flux bc
                                    );
                // Assemblying Matrix ---------------------------------
                for ( int lsj_node=0; lsj_node< elb_ndof2;  lsj_node++ ) {
                    _KeM ( lei_node,sur_toply[lsj_node] ) +=
                        dtxJxW_g*bc_alpha*alpha_eff*phii_g*_phi_g[2][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
                }// end lsj_node
            }// lsi_node

        } // end of the quadrature point qp-loop **********************

    } // ======================== end Neumann  boundary conditions  ===========================
    return;
}


#endif
