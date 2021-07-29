// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// #if NS_EQUATIONS==3
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// class files --------------------------------------------------
// #include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS_1comp.h"       // Navier-Stokes class header file
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class

void MGSolNS_1comp::get_el_field_data (
    int iel, int Level,
    int el_conn [], int offset, int el_ndof[], int ndof_lev,
    double u_old[],  double u_oold[],   double u_nl[]
) {
    int el_ndofp = el_ndof[1];
    for ( int deg = 0; deg < 3; deg++ )
        for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
            _data_eq[deg].mg_eqs[eq]->get_el_sol (0, 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                                   el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
        }
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F ]]->get_el_nonl_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, u_nl );
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F ]]->get_el_sol ( 0, 0, 1, el_ndof[2], el_conn, offset, 0, u_old );
    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F ]]->get_el_sol (1, 0, 1, el_ndof[2], el_conn, offset, 0, u_oold );

    return;
}

// ==============================================================================================
void MGSolNS_1comp::matrixrhsvol (
    DenseMatrixM &KeM, DenseVectorM &FeM,
    int el_ndof[],
    double u_old[],
    double u_nl[],
    const int unsteady_flag,
    const int mode,
    int el_conn[]
) {
    double dphijdx_g2[DIMENSION];
    double dphiidx_g2[DIMENSION];
    double vel_gdx[DIMENSION];
    const int el_ngauss = _fe[2]->_NoGauss1[ _nNSdim - 1];              //elem gauss points
    const int el_ndof2 = el_ndof[2];
    
    double utau = 0.;
    double WallOnG[1];

    for ( int qp = 0; qp <  el_ngauss; qp++ ) {

        // shape functions at gaussian points (qp) --------------------------------------------------------------------------------
        // quadratic continuous (2)  (velocity)
        const double det2      = _fe[2]->Jac ( qp, _xx_qnds, _InvJac2 ); // quadratic Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[ _nNSdim - 1][qp]; // quadratic weight
        _fe[2]->get_phi_gl_g ( _nNSdim, qp, _phi_g[2] );                // quadratic shape function
        _fe[2]->get_dphi_gl_g ( _nNSdim, qp, _InvJac2, _dphi_g[2] );    // global coord deriv

        // discontinuous (0) (disc pressure)
        if ( _nvars[0] > 0 ) _fe[0]->get_phi_gl_g ( _nNSdim, qp, _phi_g[0] ); // piecewise shape function

        // interpolation fields at gaussian points (qp) ---------------------------------------------------------------------------
        interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], el_ndof2, _ub_g[2] ); // field _ub_g[2][DIM]
        interp_el_gdx ( u_old, 0, 1, _dphi_g[2], el_ndof2, vel_gdx ); // derivatives  vel_gdx[DIM][DIM]


        
        if ( _AxiSym == 1 ) {
            interp_el_sol ( _xx_qnds, 0, _nNSdim, _phi_g[2], el_ndof2, _xyzg );
            JxW_g2  *= _xyzg[0];
        }
        // Velocity, Reynolds and upwind  --------------------------------------------------------------------------
        _sP = 1.e-20;
        for ( int jdim = 0; jdim <  _nNSdim; jdim++ ) _sP += 2.* ( vel_gdx[ jdim ] * vel_gdx[ jdim] );
        _mu_turb = 0.;
        // -------------------- Turbulence [K_F] -> (quad,_indx_eqs[K_F]) -------------------------------------
	
        if ( _FF_idx[MU_T]>=0 ){  // turbulence _mu_turb evaluation
            _mu_turb = _IRe*_ub_g[2][_FF_idx[MU_T]];
        }
        
        double IRe_eff = _IRe+_mu_turb;

        for ( int i = 0; i < el_ndof2; i++ ) { // LOOP OVER TEST FUNCTIONS -------------------------------------------------
            // set up test function phii_g, derivatives dphiidx_g2[dim], and supg Phi_supg
            int indx_eq     = i ;
            double phii_g = _phi_g[2][i];
            for ( int idim = 0; idim < _nNSdim; idim++ ) dphiidx_g2[idim] = _dphi_g[2][i + idim * el_ndof2];

            // zero line for Dirichlet  bc -----------------------------------------------------------
            if ( _bc_el[indx_eq] <=0 ) {
                const int ivar = 0;
                // Regularization supg term ------------------------------------------------------------------------------------
                if ( _NS_parameter._SolveSteady == 0 ) 
                    FeM ( i )  +=  JxW_g2  * _ub_g[2][_FF_idx[NS_F]] * phii_g / _dt ;
                
                FeM ( i )  +=  JxW_g2  * _IFr * _dirg[1] * phii_g ;

                for ( int j = 0; j < el_ndof2; j++ ) {
                    const double phij_g = _phi_g[2][j];
                    double Lap_g = 0.;
                    if ( _AxiSym == 1 ) { // axysimmetry only --------------------------------------------------
                        Lap_g = 2.* ( IRe_eff ) * phij_g * ( phii_g ) / ( _xyzg[0] * _xyzg[0] );
                    } // --------------------------------------------------------------------------------
                    for ( int kdim = 0; kdim < _nNSdim; kdim++ ) {
                        dphijdx_g2[kdim] = _dphi_g[2][j + kdim * el_ndof2];
                        Lap_g     += IRe_eff * dphijdx_g2[kdim] * dphiidx_g2[kdim];
                    }
                    //--------------------------------- Diagonal blocks [1-5-9] ------------------------------------------------------
                    const double TimeDerivative = JxW_g2 * phii_g * phij_g / _dt ;
                    const double KeM_ij = JxW_g2 * Lap_g;
                    if ( _NS_parameter._SolveSteady == 0 ) KeM ( i, j ) += TimeDerivative;
                    KeM ( i, j ) += KeM_ij;
                } // end A element matrix quad -quad (end loop on j)---------------------------------------------
            } // end loop ivar
        } // end loop i
    }
    return;
}


#endif
// #endif









