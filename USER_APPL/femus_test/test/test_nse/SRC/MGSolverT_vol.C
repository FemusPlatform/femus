#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================



#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverT.h"       // Navier-Stokes class header file

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MGFE.h"          // Mesh class,double vel_g[]



// ===============================================================================================
void MGSolT::vol_integral(
  DenseMatrixM &KeM, DenseVectorM &FeM,
  const int el_ndof2,const int el_ngauss,
  double xx_qnds[],
  const int unsteady,const int mode,const int iaxisim
) { // ==============================================================================================

  double vel_g[DIMENSION];
  double T_der[DIMENSION],T_secder[DIMENSION*DIMENSION];
  double xyz_g[DIMENSION];
  double rhocp=1.;
  
   double dt=_dt;  if(dt>1) dt=1. ; if(dt<1.e-10) dt=1.e-10;
   
// --------------------------------------------------------------------------------------------------------------------
  /// c) gaussian integration loop (n_gauss)
  // ------------------------------------------------------------------------------------------------------------------
  for(int qp=0; qp< el_ngauss; qp++) {
    // shape functions at gaussian points -----------------------------------
    double det2      = _fe[2]->Jac(qp,xx_qnds,_InvJac2);     // Jacobian
    double JxW_g2 =det2*_fe[2]->_weight1[_nTdim-1][qp]*dt;       // weight
    _fe[2]->get_phi_gl_g(_nTdim,qp,_phi_g[2]);               // shape funct
    _fe[2]->get_dphi_gl_g(_nTdim,qp,_InvJac2,_dphi_g[2]); // global coord deriv
    _fe[2]->get_ddphi_gl_g(_nTdim,qp,_InvJac2,_ddphi_g[2]); // local second deriv

    //  fields --------------------------------------------------------------------------------------------------------
    interp_el_sol(_xx_qnds,0,_nTdim,_phi_g[2],el_ndof2,xyz_g);
    interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2]); // quadratic
    //  derivatives
    interp_el_gdx(_data_eq[2].ub,_FF_idx[T_F],1,_dphi_g[2],el_ndof2,T_der);
    interp_el_gddx(_data_eq[2].ub,_FF_idx[T_F],1,_ddphi_g[2],el_ndof2,T_secder);
    
   // axisymmetric (index -> 0)
   if(iaxisim==1) JxW_g2  *=xyz_g[0];

    // Temperature field -> [T_F] -> (quad, _indx_eqs[NS_F])---------------------------------------->
    double alpha_eff = _IPrdl*_IRe;
// #ifndef CONST
//     rhocp *=densityT(uold_g[0])*cipiT(uold_g[0]);// nondimensional density (rho/rhoref)
//     alpha_eff *=kappaT(uold_g[0]);// nondimensional viscosity (mu/muref)
// #endif
    // Velocity field -> [NS_F] -> (quad, _indx_eqs[NS_F]) ----------------------------------------->
    if(_FF_idx[NS_F]>-1) {
      for(int idim=0; idim<  _nTdim; idim++) { vel_g[idim] =_ub_g[2][_FF_idx[NS_F]+idim]; } // velocity field
    } else {
      vel_g[_nTdim-1] = 0.;    vel_g[0] = .0;  vel_g[1] = 1.;
    }
    double mod2_vel=1.e-20; for(int idim=0; idim<  _nTdim; idim++) { mod2_vel +=vel_g[idim]*vel_g[idim]; }
    mod2_vel =sqrt(mod2_vel); // velocity modulus
    
    // h_eff=2/sum|s.dN| and f_upwind ---------------------------------------------------------------------------------
    double h_eff=1.e-21;
    for(int i=0; i<el_ndof2; i++) {
      double hh=1.e-20; for(int idim=0; idim< _nTdim; idim++) hh += vel_g[idim]*_dphi_g[2][i+idim*el_ndof2]/mod2_vel;
      h_eff += fabs(hh);
    }
    h_eff=2./h_eff; if(h_eff<1.e-10) {h_eff=1. ; std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";   }

    double Pe_h=0.5*mod2_vel*h_eff/(_IRe*_IPrdl);
    double a_opt=(1./tanh(Pe_h)-1./Pe_h); if(a_opt >1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
    double f_upwind=rhocp*0.5*a_opt*h_eff/(mod2_vel);   // upwind
    
    double T_g=_ub_g[2][_FF_idx[T_F]]+273.15;
    rhocp =1.;//(3.284+0.01617*T_g-2.305e-6*T_g*T_g)/_kappa0;
    alpha_eff = _IPrdl*_IRe;//*(3.284+0.01617*T_g-2.305e-6*T_g*T_g)/_kappa0;
    
    if(_FF_idx[K_F]>=0) {// turbulence _mu_turb evaluation
     double kappa_g=fabs(_ub_g[2][_FF_idx[K_F]]); double omega_g=fabs(_ub_g[2][_FF_idx[K_F]+1])+1e-20;
     double mu_turb= kappa_g/omega_g;   //   _mu_turb=eval_var1(val_tbg);
      if(mu_turb<1.e-4) mu_turb=1.e-4;
        _IPrdl_turb=1./_T_parameter._Prt;//        _IPrdl_turb=2.;    _IPrdl_turb=0.85+1/(1+ mu_turb*145.7*_muf0/kappa);
      alpha_eff += mu_turb*_IRe*_IPrdl_turb;
    }
// alpha_eff =

    /// d) Local (element) assemblying energy equation
    // =====================================================================================================================
    for(int i=0; i<el_ndof2; i++)     {
 
      double dtxJxW_g=JxW_g2*_bc_el[i];           // area with bc and weight
      double Phi_supg=0.,Lap_expl=0.,Lap_supg=0.,Adv_expl=0.;  // supg, Lap explicit , Lap explicit supg
      for(int idim=0; idim< _nTdim; idim++) {
        const double  dphiidxg=_dphi_g[2][i+idim*el_ndof2];
        Adv_expl += vel_g[idim]*T_der[idim];
        Lap_supg += alpha_eff*T_secder[idim*_nTdim+idim];       // explicit Laplacian supg
        Lap_expl += alpha_eff*T_der[idim]* dphiidxg;            // explicit Laplacian
        Phi_supg += _T_parameter._Supg*_bc_el[i]* f_upwind*vel_g[idim]* dphiidxg; // phii_g+
      }

      const double phii_g=_phi_g[2][i];

      // Rhs Assemblying  ---------------------------------------------------------------------------------------------------
      if(mode == 1) { // rhs
        FeM(i) += dtxJxW_g*(
                    (unsteady*rhocp*_ub_g[2][_FF_idx[T_F]])*(phii_g+Phi_supg)/_dt  // time
                     +_qheat*(phii_g+Phi_supg)                             // source term
                  );
      }
      
 
      // Matrix Assemblying ------------------------------------------------------------------------------------------------
      for(int j=0; j<el_ndof2; j++) {

        double Adv=0.,Lap=0.,Lap_supgi=0.;
        for(int idim=0; idim<  _nTdim; idim++) {
          const double  dphiidxg=_dphi_g[2][i+idim*el_ndof2];
          const double  dphijdxg=_dphi_g[2][j+idim*el_ndof2];
          Adv +=vel_g[idim]*dphijdxg;                                                                // advection
          Lap += alpha_eff*dphijdxg*dphiidxg;                                                        // diffusion
          Lap += _T_parameter._Upwind*f_upwind*vel_g[idim]*vel_g[idim]*dphijdxg*dphiidxg;             // normal upwind 
          
          
          for(int jdim=idim+1; jdim< idim+ _nTdim; jdim++)  Lap += 0.*
            f_upwind*vel_g[(jdim%_nTdim)]*vel_g[idim]*_dphi_g[2][j+(jdim%_nTdim)*el_ndof2]*dphiidxg; // skew upwind Laplacian   _T_parameter.UPWIND2
          Lap_supgi += alpha_eff*_ddphi_g[2][j*_nTdim*_nTdim+idim*_nTdim+idim];
        }

        const double phij_g= _phi_g[2][j];
        KeM(i,j) +=dtxJxW_g*( // energy-equation matrix
                     unsteady*rhocp*phij_g*(phii_g +Phi_supg)/_dt // time term
                     + rhocp*Adv*(phii_g+Phi_supg)     //advection term
                     + Lap                             //diffusion term
                     - Lap_supgi*Phi_supg              //diff supg term
                   );
//         std::cout << "Matrix \n"  <<KeM  << "FeM\n"  <<FeM;
      }
    } // ----------------------------------------
  } // end of the quadrature point qp-loop ***********************

  return;
}

#endif
