// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// local Femus class include -----------------------------------
// class files --------------------------------------------------
// #include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS_1comp.h"       // Navier-Stokes class header file
#include "MGFE.h"          // Mesh class
#include <iomanip>      // std::setprecision

void MGSolNS_1comp::set_bc_matrix (
    DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
    int dir_maxnormal,                          ///<  normal dir
    int sur_toply[],                            ///< boundary topology map
    int el_ndof[],                              ///< number of volume dofs
    int elb_ndof[],                             ///< number of boundary dofs
    int elb_ngauss,                             ///<  number of surface gaussian points
    double normal[],                            ///< normal
    double u_old[],
    int el_conn[]
) {

    const int bd_face = abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 100;
    const int norm_face = ( abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 10000 ) / 1000 - 1 ;
    double U_old[NDOF_FEM];
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol ( 0, 0, 1, el_ndof[2], el_conn, _offset, 0, U_old );

    for ( int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++ ) {
        int gnode = sur_toply[lbnode];
        const int bd_node = abs ( _bc_bd[gnode] ) % 100;
        if ( bd_node==88 ) {
            _bc_el[gnode] = 1;
            KeM ( gnode, gnode ) = 1.;
            FeM ( gnode ) = 0.;
        }
        if ( bd_node==98 || bd_node==99 ) {
            _bc_el[gnode] = 1;
            KeM ( gnode, gnode ) = 1.;
            FeM ( gnode ) = U_old[gnode];
        }
    }
    double VelOnGauss = U_old[NDOF_FEM-1];
    double u_tau = 0.;
    double y_plus = 0.;
    double eff_stress = 0.;
    
    
    for ( int  qp=0; qp< elb_ngauss; qp++ ) { // START LOOP OVER GAUSS POINTS ===============================================
        double det   = _fe[2]->JacSur ( qp,_xxb_qnds, _InvJac2 ); // local coord _phi_g and jac
        double JxW_g = det*_fe[2]->_weight1[_nNSdim-2][qp]/**dt*/;  // weight
        _fe[2]->get_phi_gl_g ( _nNSdim-1,qp,_phi_g[2] );    // global coord _phi_g
        _fe[1]->get_phi_gl_g ( _nNSdim-1,qp,_phi_g[1] );    // global coord _phi_g

        if ( _AxiSym==1 ) {
            interp_el_bd_sol ( _xx_qnds,sur_toply,el_ndof[2],0,_nNSdim,_phi_g[2],elb_ndof[2],_xyzg );
            JxW_g  *=_xyzg[0];    // axisymmetric  (index ->0)
        }


        const double WDist =_mgutils._geometry["Wall_dist"] +_y_bcout;// _mgutils._TurbParameters->_BoundWallDist;
        if ( _FF_idx[K_F]>-1 ) {
            // Calculation of effective stress through u_tau
            u_tau      = _mgutils._TurbParameters->CalcUtau ( VelOnGauss + 1.e-10,  WDist );
            y_plus     = _mgutils._geometry["Wall_dist"]*u_tau/_IRe;
            eff_stress = ( fabs ( VelOnGauss ) < 1.e-10 ) ? _IRe/ ( _mgutils._geometry["Wall_dist"] ) : -1.*u_tau*u_tau;
	    
// 	    eff_stress = -_IRe/ ( WDist );
            if(bd_face== 84 && qp==0) std::cout<<"utau "<<u_tau<<"  yplus "<<y_plus<<"  wd "<<_mgutils._geometry["Wall_dist"]<<"  vel "<<VelOnGauss<<std::endl;
        }

    if (bd_face==84){ 
        for ( int i=0; i< elb_ndof[2]; i++ ) { // LOOP OVER TEST FUNCTIONS ===================================================

            const int bc_var_check2=_bc_bd[sur_toply[i]];  // total  bc_var
            const int bc_var_npt= bc_var_check2/10000;                              // recurrent  points in one element
            const int sign = bc_var_check2/abs ( bc_var_check2 );
            const int bc_var_check= abs ( _bc_bd[sur_toply[i]] ) %10000;                        // bc_var
            const int bc_var_normal= bc_var_check/1000-1;                                 //
            const int bd_node = abs ( _bc_bd[sur_toply[i]] ) %100;
            const int bc_varn_flag=bd_node/10;
            const int bc_vartg_flag=bd_node%10;

            double phii_g  = _phi_g[2][i];

            if ( bd_face==bd_node ) {

                    int indx     = sur_toply[i];
                    int NComp = ( _bc_el[indx] < 0 ) ? _nNSdim:1;
                    int bc_el = ( _bc_el[indx]<=0 ) ? 1 :0;
                    const double dtJxW_g =bc_el* JxW_g;

                    bool NormBC = ( _bc_el[indx] == -3 ) ? true : false;

                    // Stress component along tangential directions (-1,-2) ---------------------------------------------------
                    const int  bc_rhs= ( ( bc_vartg_flag%2 ) == 1 ) ? 1:0;

                    double stress = 0;
		    double exp_cont = 0.;
		    
		    FeM ( indx ) += - (exp_cont)*u_tau*u_tau*dtJxW_g*_phi_g[2][i];
// 		    
                    for ( int j=0; j<elb_ndof[2]; j++ ) {
		        eff_stress = ( fabs ( U_old[sur_toply[j]] ) < 1.e-5 ) ? -1.*_IRe/ ( WDist ) : -1.*u_tau*u_tau/U_old[sur_toply[j]];
// 			eff_stress = -1.*_IRe/ ( WDist );
			stress = eff_stress*dtJxW_g*_phi_g[2][i];
                        KeM ( indx,sur_toply[j] ) += - (1.-exp_cont)*stress*_phi_g[2][j];
                    }// END LOOP OVER MATRIX ROWS RELATIVE TO TEST FUNCTION I

            }// END GAUSSIAN INTEGRATION ========================================================================================
	}
    }
    
 }
    
            return;
        }

#endif






