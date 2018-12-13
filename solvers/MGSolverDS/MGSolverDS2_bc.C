// ===============================================================
// --------------   DS system [SDSX_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef DS_EQUATIONS
// ==============================================================
// DS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// DS_EQUATIONS==1 coupled    solver (u,v,w,p)
// DS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// local Femus class include -----------------------------------
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverDS.h"       // Navier-Stokes class header file
#include "MGFE.h"          // Mesh class



void MGSolDS::set_bc_matrix(
  DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
  int dir_maxnormal,                          ///<  normal dir
  int sur_toply[],                            ///< boundary topology map
  int el_ndof[],                              ///< number of volume dofs
  int elb_ndof[],                             ///< number of boundary dofs
  int elb_ngauss,                             ///<  number of surface gaussian points
  double normal[],                            ///< normal
  const int iaxis                              ///< axisymmetric
) {

  double eff_stress=1.; double pressure;
  double vel_bound[1]; double xyz_g[DIMENSION];
//   const int elb_dof0=(NDOF_K>1)?elb_ndof[_pres_order]*NDOF_K:elb_ndof[_pres_order];
  double det2= _fe[2]->JacSur(elb_ngauss-1,_xxb_qnds, _InvJac2);// jacobian
  double Ipenalty=det2;                            // Dirichlet bc flagdouble xyz_bg[DIMENSION];


  //  ============================================================================================
  //  =================== Dirichlet boundary condition (Dirichlet=0; Neumann=1) ====================
  //  ============================================================================================
  // enum bound_cond
  // velocity_in0=      1,   velocity_tg0=      2,   wall=              3,
  // velocity_in=       5,   velocity_tg=       6,   velocity=          7,
  //  ============================================================================================

  for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {

    int bc_dn=(_bc_vol[sur_toply[lbnode]]/10);
    if(bc_dn==0) {  // dirichlet -------------------------------------------------------------------------------------------------------------------

      int bc_s=_bc_vol[sur_toply[lbnode]]%10;
      //  (0;??1) normal              (0;?1?) tg  // (0;1??) nonhomoge neous
      int bc_normal=(bc_s%2); int bc_tg=((bc_s&2)>>1); int bc_rhs=((bc_s&4)>>2);
     
        const int  indx_sol=sur_toply[lbnode]+_dir*el_ndof[2];//ivar=0;
        const int  indx_row=sur_toply[lbnode];//ivar=0;

        if(_dir == dir_maxnormal)  {  // normal ---------------------------------------------------------------------------------------------------
       
          _bc_el[sur_toply[lbnode]]=0;
          //          Dirichlet (homogenoues and non  homogenoues) Assemblying
          const double nii=Ipenalty*normal[_dir]/**normal[ivar+_dir]*/;
          KeM(indx_row,indx_row) = nii; //diagonal ->  u_n=n(u*nx+v*ny+w*nz)
          FeM(indx_row) = bc_normal*bc_rhs*nii*_data_eq[2].ub[indx_sol+(_FF_idx[SDSX_F]+_dir)*el_ndof[2]]; // only normal non-homogeneous
          for(int   jvar=_dir+1; jvar< _nDSdim+_dir; jvar++)    {  // non-diagonal part  ->  u_n=n(u*nx+v*ny+w*nz)
            int      jvar1=jvar%DIMENSION;
            const double nij= Ipenalty*normal[jvar1]/**normal[ivar+_dir]*/;
            const double valj=_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
            FeM(indx_row) += (bc_rhs*bc_normal-1)*nij*valj; //  KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += aij;
          }// jvar loop

        } //  end  Dirichlet normal --------------------------------------------------------------------------------------------------------------------------

        else  {  // Non normal direction  ---------------------------------------------------------------------------------------------------------------

          _bc_el[sur_toply[lbnode]]=0; // storage bc for vol part

          //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
          // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
          //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
          // diagonal part
          const double nii= normal[_dir]*normal[_dir]*fabs(bc_tg-bc_normal);
          KeM(indx_row,indx_row) =  Ipenalty*(1-nii);
          FeM(indx_row) =  Ipenalty*bc_rhs*bc_tg*(1-nii)*_data_eq[2].ub[indx_sol+(_FF_idx[SDSX_F]+_dir)*el_ndof[2]];  
          // non-diagonal part
          for(int   jvar=_dir+1; jvar< _nDSdim+_dir; jvar++)    {
            int      jvar1=jvar%_nDSdim;
            const double nij=Ipenalty*normal[jvar1]*normal[_dir]*fabs(bc_tg-bc_normal);
            const double valj=_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
            FeM(indx_row) -= (bc_rhs*bc_tg-1)*nij*valj;
//             FeM(indx_row) += aij*valj*fabs(bc_tg-bc_normal);

//             KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += -aij;
          }// jvar loop
        } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------

      }  // dc
 

    // ==============================================================================================================
    // simmetry for vector field
    // simmx=            21,   simmy=            22,   simmz=            23,
    // simmxy=           24,   simmxz=           25,   simmyz=           26 <-  not used
    // ===============================================================================================================
    if(bc_dn==2) {  // simmetry--------------------------------------------------------------------------------------------
      int bc_s=_bc_vol[sur_toply[lbnode]]%10;

      if((bc_s%3) <4) { // simmx=21,24 (x-plane),simmy=22,25(y-plane);simmz=23,26(z-plane);
        int sim_dir= (bc_s%3)-1 ;
        if( sim_dir == dir_maxnormal && _dir==dir_maxnormal)  {  // normal ---------------------------------------------------------------------------------------------------
//           const int  indx_sol=sur_toply[lbnode]+( sim_dir+_dir)*el_ndof[2];//ivar=0;
          const int  indx_row=sur_toply[lbnode];//ivar=0;
          _bc_el[sur_toply[lbnode]]=0;
          //          Dirichlet (homogenoues and non  homogenoues) Assemblying
          KeM(indx_row,indx_row) = 1.;FeM(indx_row) = 0.; //diagonal ->  u_n=n(u*nx+v*ny+w*nz)
          
        } //  end  ivar == dir_maxnormal------------------------------------------------------------------------------------------------------------
      }

      
      // -----------------------------------------------------------
    }// dc  end dc==0
  } // lnopde


  //  ============================================================================================
  //  ====================================== Neumann  ============================================
  // outflow=          10,   interior=         11,   pressure_outlet=  12,
  // wall_turb=        13,   outflow_p=        14,   pressure_inlet=   16,
  //  ============================================================================================
  int ivar=0;
//   for(int  ivar=0; ivar< _nvars[2]; ivar++)
  {
    if(_bc_vol[sur_toply[NDOF_FEMB-1]]/10 ==1) {  /*bc_dn=1; }*/
//     if(bc_dn==1) {
      // Neumann  -------------------------------------------------------------------------------------------------------------------
      int bc_s=_bc_vol[sur_toply[NDOF_FEMB-1]]%10;
      int bc_normal=(bc_s%2);  int bc_tg=((bc_s&2)>>1); int bc_rhs=((bc_s&4)>>2); // (0;?1?) tg (0;??1) normal (0;1??) nonhomogeneous
   
//       if(bc_s !=3)    for(int k=0; k<elb_dof0 ; k++) { _bc_el[_pres_order*sur_toply[k]+ _nDSdim*NDOF_FEM]= 0  ; }// pressure bc
      int flag_normal=(ivar+_dir == dir_maxnormal)?1:0 ;
      if(bc_s == 3) {
        flag_normal= !flag_normal; // turb condition  bc_tg=1;bc_normal=1;
      }
      if(flag_normal)  {  // normal ---------------------------------------------------------------------------------------------------
        double   det;
        double sign=1.; if(normal[dir_maxnormal]<0) {sign=-1.;}
        for(int  qp=0; qp< elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)
          // quad/linear  [2]=quad [1]=linear------------------------------------
          det     = _fe[2]->JacSur(qp,_xxb_qnds, _InvJac2);   // local coord _phi_g and jac
          double  JxW_g=det*_fe[2]->_weight1[_nDSdim-2][qp]; // weight
          _fe[2]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[2]);   // global coord _phi_g
          _fe[1]->get_phi_gl_g(_nDSdim-1,qp,_phi_g[1]);   // global coord _phi_g
          interp_el_bd_sol(_xx_qnds,sur_toply,el_ndof[2],0,_nDSdim,_phi_g[2],elb_ndof[2],_xyz_g);

          if(iaxis==1) { JxW_g  *=_xyz_g[0]; } // axisymmetric  (index ->0)

          // Assemblying NS equation
          for(int i=0; i< elb_ndof[2]; i++) {  // Assemblying NS equation
            // set up row i
            const double phii_g=_phi_g[2][i];
            const int   indx_var=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
            const int   indx_row=sur_toply[i]+(ivar)*el_ndof[2];// volume dof index
            // boundary flag
            int bc_v1=(int)(_bc_vol[sur_toply[i]]/10);
            //                 int bc_s1=(int)_bc_bd[indx_row];
            //                  bc_rhs   =((bc_s1&4)>>2); // (1??) -> nonhomogeneous
            //                  bc_tg    =((bc_s1&2)>>1); // (?1?) -> tg
            //                  bc_normal=(bc_s1%2);      // (??1) -> normal
            double dtJxW_g=JxW_g*_bc_el[sur_toply[i]+(ivar+_dir)*NDOF_FEM];

            // wall_turb boundary conditions
            double un=0;  for(int  kdim=0; kdim<_nDSdim; kdim++) { un +=_data_eq[2].ub[sur_toply[i]+(_FF_idx[SDSX_F]+kdim)*NDOF_FEM];}
//               u_old[sur_toply[i]+kdim*NDOF_FEM]*normal[kdim]; }
            double ut=0;  for(int  jvar=0; jvar<_nDSdim; jvar++) {
              const double tmp=_data_eq[2].ub[sur_toply[i]+(_FF_idx[SDSX_F]+jvar)*NDOF_FEM]-normal[jvar]*un; ut +=tmp*tmp;
            }
            vel_bound[0] =sqrt(ut);    double u_tau= eval_var2(vel_bound);
#ifdef TBK_EQUATIONS
            eff_stress =u_tau*u_tau/fabs(vel_bound[0]+1.e-20);
#endif
            for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad ------------------------------------------------------------------------------------------
              const double phij_g= _phi_g[2][j];
//               if(j<elb_ndof[1]) {
//                 double pressure_ctrl=1*_data_eq[_pres_order].ub[_pres_order*sur_toply[j]];
//                 FeM(indx_row)  += -1.*sign*dtJxW_g*phii_g*_phi_g[1][j]*pressure_ctrl;//_refvalue[3]);
//               }
              KeM(indx_row,indx_row) +=   dtJxW_g*bc_tg*bc_normal*eff_stress*phij_g*phii_g;
              for(int  jvar=ivar+_dir+1; jvar< _nDSdim+_dir; jvar++)    {  // u
//                 KeM(indx_row,sur_toply[j]+(jvar%DIMENSION)*el_ndof[2]) +=    dtJxW_g*eff_stress*phij_g*phii_g*bc_normal*bc_tg;
              }// jvar-loop
            }//  j-loop ---------------------------------------------------------------------------------------------------------
          }// i-loop
        }// end gaussian  integration
      }//  end  normal direction
      //  -------------------------------------------------------------------------------------------------------------------
      else  {  // Non normal direction  -------------------------------------------------------------------------------
        for(int lbnode=0; lbnode< elb_ndof[2]; lbnode++) {  // Assemblying NS equation
          // set up row                int  iflag=0;
          bc_s=(_bc_vol[sur_toply[lbnode]]%10)%4;
          bc_tg=((bc_s&2)>>1);   // (0;?1?) tg
          bc_normal=(bc_s%2);    // (0;??1) normal
          bc_rhs=((bc_s&4)>>2);  // (0;1??) nonhomogeneous
          if(_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM] != 0) {
            const int  indx_sol=sur_toply[lbnode]+(ivar+_dir)*el_ndof[2];//ivar=0;
            const int  indx_row=sur_toply[lbnode]+ivar*el_ndof[2];//ivar=0;
            _bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]= ((bc_tg+bc_normal) == 0) ? 1:0 ;  // storage bc for vol part
            //          Dirichlet (homogenoues and non  homogenoues) Assemblying Matrix
            // u=u_n+u_t;   u_n=n(u*n) ; u_t=u-(u.n)n  flags=(?,bc_normal,bc_tg,bc_rhs)
            //  (?00?) single comp bc (?01?) Normal Dirichlet bc ; (?10?) tg Dirichlet bc
            // diagonal part
            KeM(indx_row,indx_row) += ((bc_tg+bc_normal)%2)*Ipenalty;
            double aii=Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*((bc_tg+bc_normal));
            KeM(indx_row,sur_toply[lbnode]+ivar*el_ndof[2]) += -Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(bc_tg+bc_normal);
            FeM(indx_row) +=  bc_rhs*bc_tg*Ipenalty*_data_eq[2].ub[_FF_idx[SDSX_F]*NDOF_FEM+indx_sol]
                              -bc_tg*bc_rhs*aii*_data_eq[2].ub[_FF_idx[SDSX_F]*NDOF_FEM+indx_sol];   // Dirichlet tang  homogeneous only
            // non-diagonal part
            for(int   jvar=ivar+_dir+1; jvar< _nDSdim+ivar+_dir; jvar++)    {
              int      jvar1=jvar%DIMENSION;
              double aij=_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]*Ipenalty*normal[jvar1]*normal[ivar+_dir]*(bc_tg+bc_normal);
              double valj=_bc_el[sur_toply[lbnode]+ivar*NDOF_FEM]*_data_eq[2].ub[(_FF_idx[SDSX_F]+jvar1)*NDOF_FEM+sur_toply[lbnode]];
              FeM(indx_row) += -bc_rhs*bc_tg*aij*valj;
//               KeM(indx_row,sur_toply[lbnode]+jvar1*el_ndof[2]) += -aij;
              FeM(indx_row) +=aij*valj;
            }// jvar loop
            //                     bc_rhs*Ipenalty*u_tau[i]   + bc_yplus*Ipenalty*y1_plus[i];
          } // end non-diagonal part ---------------------------------------------------------------------------------------------------------------------------
        }
      }
    }  //  ivar |||||
  }  //  bc_dn |||||
  // ========================================= Dirichlet boundary conditions ==============================
  return;
}


#endif
