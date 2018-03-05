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
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverDS.h"       // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class


// Thermodinamical Properties ==========================================================================================
// constant
// double  _CTYLdensity(double T,double T_r) {return 1.;} // water (t K) T_ref
// double  _CTRLkviscosity(double T,double T_r) {return 1.;} // lead T K
// water
//     double  _DSdensity(double T,double T_r){return (1. - (T-3.9863)*(T-3.9863)*(T+288.9414)/(508929.2*(T+68.12963)))/
//         (1. - (T_r-3.9863)*(T_r-3.9863)*(T_r+288.9414)/(508929.2*(T_r+68.12963)));} // water (t C) T_ref
//      double  _DSkviscosity(double T){return 1.;} //  T K
//  // Lead
//     double  _DSdensity(double T,double T_r){return (11441-1.2795*T)/(11441-1.2795*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return return 4.55e-4*exp(1069/T);} // lead T K
//       // LBE
//     double  _DSdensity(double T,double T_r){return (11065-1.293*T)/(11065-1.293*T_r);} // water (t K) T_ref
//     double  _NSviscosity(double T){return 4.94e-4*exp(754.1/T);} // lead T K
//
//  ======================================================================================================================

// // v/*oid MGSolDS::get_el_field_data(
//   int iel, int Level,
//  int el_conn [], int offset,int el_ndof[],int ndof_lev,
//   double u_old[],  double u_oold[],   double u_nl[],
//   double p_proj[], double dp_proj[]
// ) {
//
//   int el_ndofp=el_ndof[1];
//
//   // element nodes coordinates ----------------------------------------------------------------------------------------
//   for(int idim=0; idim<DIMENSION; idim++) for(int d=0; d< NDOF_FEM; d++) {
//       _data_eq[2].ub[idim*NDOF_FEM+d]=_xx_qnds[idim*NDOF_FEM+d];
//     }
//   // external fields (from constant 0 to quadratic 2) -------------------------------------------------------------------
//   for(int deg=0; deg<3; deg++)
//     for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
//       _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
//                                            el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
//     }
//
//   //  internal field data (NS) -------------------------------------------------------------------------------------------
//
// #if DS_EQUATIONS==0 // pressure as external field (projection)// --------------------------------------------------------------
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_nonl_sol(0,_nDSdim,el_ndof[2],el_conn, offset,0,u_nl);
//   _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0,1,el_ndofp,el_conn, offset,0,p_proj);
//   _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndofp,el_conn, offset,0,dp_proj);
//
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_sol(0,_nDSdim,el_ndof[2],el_conn, offset,0,u_old);
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_oldsol(0,_nDSdim,el_ndof[2],el_conn, offset,0,u_oold);
// #endif
//
// #if DS_EQUATIONS==1 //   coupled ----------------------------------------------------------------------------------------
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_sol(_nDSdim,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub);    // pressure
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_sol(0,_nDSdim,el_ndof[2],el_conn, offset,0,u_old);     // old vel
//   _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F]]->get_el_nonl_sol(0,_nDSdim,el_ndof[2],el_conn, offset,0,u_nl);   //  non linear vel
// #endif
//
// #if DS_EQUATIONS==2 // pressure as external field (splitting) -----------------------------------------------------------
//   for(int kdim=0; kdim< _nDSdim; kdim++) {
//     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F+kdim]]->get_el_nonl_sol(0,1,el_ndof[2],el_conn,offset,kdim,u_nl);
//     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F+kdim]]->get_el_sol(0,1,el_ndof[2],el_conn, offset,kdim,u_old);
//     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DS_F+kdim]]->get_el_oldsol(0,1,el_ndof[2],el_conn, offset,kdim,u_oold);
//   }
//   _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0,1,el_ndofp,el_conn, offset,0,p_proj);       // old pressure
//   _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndofp,el_conn, offset,0,dp_proj);   // dp pressure
// #endif
//
//   return;
//
//
// }*/




// ==============================================================================================
void MGSolDS::matrixrhsvol_sol_ds(
  DenseMatrixM &KeM,
  DenseVectorM &FeM,
  int el_ndof[]
) {// ==================================  Volume ===============================================

  double u_old; double v_old;
  double ddt=0.001*_dt;
  // ======================================================================
  // Volume =============================================================
  // ======================================================================

  /// c) gaussian integration loop (n_gauss)
  // --------------------------------------------
  const  int  el_ngauss = _fe[2]->_NoGauss1[ _nDSdim-1];                //elem gauss points
  for(int qp=0; qp< el_ngauss; qp++) {


    // shape functions at gaussian points -----------------------------------
    double  det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);     // Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[_nDSdim-1][qp];       // weight
    _fe[2]->get_phi_gl_g(_nDSdim,qp,_phi_g[2]);               // shape funct
    _fe[2]->get_dphi_gl_g(_nDSdim,qp,_InvJac2,_dphi_g[2]); // global coord deriv
    for(int ideg=0; ideg<2; ideg++) if(_nvars[ideg]>0) { _fe[ideg]->get_phi_gl_g(_nDSdim,qp,_phi_g[ideg]); }    // shape funct
    double  det1      = _fe[1]->Jac(qp,_xx_qnds,_InvJac2);     // Jacobian
    _fe[1]->get_phi_gl_g(_nDSdim,qp,_phi_g[1]);               // shape funct
    _fe[1]->get_dphi_gl_g(_nDSdim,qp,_InvJac2,_dphi_g[1]); // global coord deriv
    interp_el_sol(_xx_qnds,0,_nDSdim,_phi_g[2],el_ndof[2],_xyz_g);
#ifdef AXISYM
    JxW_g2  *=_xyz_g[0];   // std::cout<<_ub_g[2][0]<<"   "<<std::endl;
#endif
    // quadratic fields ---------------------------------------------------------------------
    interp_el_sol(_data_eq[2].ub, _FF_idx[SDSX_F], DIMENSION, _phi_g[2],el_ndof[2],_ub_g[2]);
    interp_el_sol(_data_eq[2].ub,_FF_idx[FS_F], DIMENSION, _phi_g[2],el_ndof[2],_ub_g[2]);
  

    for(int i=0; i<el_ndof[2]; i++)     {
      // set up row i
//       const double phii_g=_phi_g[2][i];
      for(int idim=0; idim< _nDSdim; idim++) { _dphiidx_g2[idim]=_dphi_g[2][i+idim*el_ndof[2]]; }
//       double dtxJxW_g=JxW_g2*(_bc_vol[i]%2)*_dt;
      double dtxJxW_g=JxW_g2*(_bc_el[i]);
      const double   ds_old=_data_eq[2].ub[(_FF_idx[SDSX_F]+_dir)*NDOF_FEM+i];    
      const double   v_old=_data_eq[2].ub[(_FF_idx[FS_F]+_dir)*NDOF_FEM+i];

      // rhs  Assemblying ---------------------------
//       FeM(i)+=dtxJxW_g*(v_old*_dt+ ds_old  )*phii_g;
       FeM(i)+=ds_old+  
               v_old*_dt;
//        cout<<ds_old<<"  "<<v_old<<endl;
      // Matrix Assemblying ---------------------------
//       KeM(i,i) += JxW_g2*(1-(_bc_el[i]));
       KeM(i,i) += 1;
//       for(int j=0; j<el_ndof[2]; j++) {
// #ifdef ALE_MOVING
//         const double phij_g=_phi_g[2][j];// set up test function
//         double Lap=0.;
//         for(int kdim=0; kdim< _nDSdim; kdim++) {
// //        _dphijdx_g2[kdim]=_dphi_g[2][j+kdim*el_ndof[2]]; Lap +=_dphijdx_g2[kdim]*_dphiidx_g2[kdim]; // Laplacian
//         }
//         // matrix assemblying
//         KeM(i,j) +=dtxJxW_g*(
//                      phij_g*phii_g           // time displ
//                    ) ;
// #endif //        END ALE MOVING MESH
//       }
    } // ----------------------------------------

  }
// ====================== end volume (element) =======================================
  return;
}


// ==============================================================================================
void MGSolDS::matrixrhsvol_liq_ds(
  DenseMatrixM &KeM,
  DenseVectorM &FeM,
  int el_ndof[],
  int flag_group[]
) {// ==================================  Volume ===============================================

  double u_old; double v_old;
  const  int  el_ngauss = _fe[2]->_NoGauss1[ _nDSdim-1];                //elem gauss points
  const int ddt=1000.;
  for(int qp=0; qp< el_ngauss; qp++) {

    // shape functions at gaussian points -----------------------------------
    double det2      = _fe[2]->Jac(qp,_xx_qnds,_InvJac2);     // Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[_nDSdim-1][qp];       // weight
    _fe[2]->get_phi_gl_g(_nDSdim,qp,_phi_g[2]);               // shape funct
    _fe[2]->get_dphi_gl_g(_nDSdim,qp,_InvJac2,_dphi_g[2]); // global coord deriv

    /// d) Local (element) assemblying energy equation
    // *********************** *******************************

    for(int i=0; i<el_ndof[2]; i++)     {

      // set up row i
//       const double phii_g=_phi_g[2][i];

      for(int idim=0; idim<_nDSdim; idim++) { _dphiidx_g2[idim]=_dphi_g[2][i+idim*el_ndof[2]]; }

//  fl_int=1;  if(interface==1) fl_int=0;
//           int fl_int= ((_bc_vol[i]>3.5)? 0:1);
      int fl_int= ((fabs(flag_group[i])>999)? 0:1);
//        if(interface==1) std::cout<< _bc_vol[i]<<" \n";
      double dtxJxW_g=JxW_g2*_bc_el[i]*fl_int;
     u_old=_data_eq[2].ub[(_FF_idx[SDSX_F]+_dir)*NDOF_FEM+i];
//         FeM(i)=dtxJxW_g*(u_old*phii_g);
      // Matrix Assemblying ---------------------------
      for(int j=0; j<el_ndof[2]; j++) {



// #ifdef ALE_MOVING

        double Lap=0.;
        // set up test function
        const double phij_g=_phi_g[2][j];
        for(int idim=0; idim< _nDSdim; idim++) {
          _dphijdx_g2[idim]=_dphi_g[2][j+idim*el_ndof[2]];
          Lap +=_dphijdx_g2[idim]*_dphiidx_g2[idim]; // Laplacian
        }
        // energy-equation
        KeM(i,j) +=dtxJxW_g*(
//                      phij_g*phii_g+
                     Lap
                   );

// #endif //   END ALE MOVING MESH


      }
    } // ----------------------------------------
  } // end of the quadrature point qp-loop ***********************


  return;
}




// // ==============================================================================================
//   void MGSolDS::set_bc_matrix ( DenseMatrixM &KeM, DenseVectorM &FeM,int dir_maxnormal,
//                          int sur_toply[], int el_ndof[],int elb_ndof[], int elb_ngauss,  double u_old[], double normal[],double Ipenalty
// ) {// ==================================  Volume ===============================================
//
//   // Dirichlet boundary conditions  ***********************************
//           if ((_bc_vol[sur_toply[NDOF_FEMB-1]]%2) ==0) {
//
//             //[NDOF_FEMB-1] is the midpoint of a quadaratic FEM element (HEX9 (2D) or HEX27 (3D))
//             int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]];     // b cond
//             double det= _fe[2]->JacSur(elb_ngauss-1,_xxb_qnds,_InvJac2);// jacobian
//             double Ipenalty=det/_dt;                               // Dirichlet bc flag
//             // local boundary loop   ---------------------------------------
//             for (int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {
//
//               int lv_node= sur_toply[lb_node]; // local vol index
//               // flag setup (\int bc_var*T+bc_var*val)
// //            int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
//               int  bc_var = (int)(bc_s%2);       // (?1) variable
//               // Assemblying  Matrix & rhs
//               int fl_int= ((_bc_vol[lv_node]>3.5)? 0:1);
// //      if (interface==1 ) std::cout<< "\n "<< _bc_vol[lv_node];
// //               if (mode == 1)
//                 FeM(lv_node) += fl_int*bc_var*Ipenalty*(_data_eq[2].ub[_FSI_idx*NDOF_FEM+lv_node]*_dt +  _data_eq[2].ub[_DS_idx*NDOF_FEM+lv_node]);
//               KeM(lv_node,lv_node) += fl_int*Ipenalty;  //  Dirichlet bc
//             }// lb_node -end  local boundary loop -------------------------
//           } // end if Dirichlet  boundary conditions
//           // **********************************************************************
//
//
//
//   return;
// }


#endif


