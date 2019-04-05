#include "Equations_conf.h"

// ============================================
#ifdef RANS_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverRANS.h"       // Navier-Stokes class header file
#include "MGFE.h"          // Mesh class



/// This function sets the _bc_el with boundary condition flags.
/// This function also assembles the surface integral to obtain the algebraic sytem.
/// It is called by MGSolRANS::GenMatRhs in MGSolverT.C
/// This function sets  the  functional defined (user write}
void  MGSolRANS::bc_set (
    int sur_toply[],
    int el_ndof2, int elb_ndof2, int elb_ngauss, int el_conn[]
)
{

    double det = _fe[2]->JacSur ( elb_ngauss - 1, _xxb_qnds, _InvJac2 ); // jacobian
    double Ipenalty = 1.;                             // Dirichlet bc flag
    double utau, normal[DIMENSION], yplus, vel_mod, x_m[DIMENSION];
    int sign;
    const int FaceBD = _bc_vol[sur_toply[NDOF_FEMB - 1]] % 100;

    {
        // CALCULATION OF UTAU AND UNIT NORMAL VECTOR
        double vel_norm, vel_bound;

        for ( int idim = 0; idim < _nTdim; idim++ ) { // based on scalar product between v and each element side
            x_m[idim] = 0.;

            for ( int d = 0; d < NDOF_FEM; d++ ) {
                x_m[idim] += _xx_qnds[idim * NDOF_FEM + d] / NDOF_FEM;
            }
        }

        _fe[2]->normal_g ( _xxb_qnds, x_m, normal, sign );

        vel_mod = vel_norm = 0.;

        for ( int dim = 0; dim < _nTdim; dim ++ ) {
            vel_mod  += _data_eq[2].ub[ ( _FF_idx[NS_F] + dim ) * NDOF_FEM + NDOF_FEM - 1] * _data_eq[2].ub[ ( _FF_idx[NS_F] + dim ) * NDOF_FEM + NDOF_FEM - 1];
            vel_norm += _data_eq[2].ub[ ( _FF_idx[NS_F] + dim ) * NDOF_FEM + NDOF_FEM - 1] * normal[dim];
        }

        vel_bound = sqrt ( vel_mod - vel_norm * vel_norm );

        double WallDist = _WallDist;

        if ( _RANS_parameter._WallFunctionApproach == 1 ) {
            WallDist = _y_dist;
        }

        utau  = _mgutils._TurbParameters->CalcUtau ( vel_bound, _y_dist );
        yplus = WallDist * utau / _IRe;
        vel_mod = sqrt ( vel_mod );
    }

    _WallElement = 0;
    if ( _RANS_parameter._WallFunctionApproach == 1 ) {
        switch ( FaceBD ) {
        case 1:
        case 12:
            _WallElement = 1;
            break;
        default:
            break;
//             _WallElement = 0;
        }
    }

    if ( _RANS_parameter._WallFunctionApproach == 0 || _WallElement == 0 ) {
        // DIRICHLET BOUNDARY CONDITIONS ============================================================
        for ( int lb_node = 0; lb_node < elb_ndof2; lb_node++ ) {

            const int NodeBC = _bc_vol[sur_toply[lb_node]] % 100;

            if ( _bc_el[sur_toply[lb_node]] == 0 && NodeBC == FaceBD ) { // ONLY IF NODE BC IS THE SAME OF BOUNDARY MID POINT BC
                int lv_node = sur_toply[lb_node]; // local vol index
                int bc_s = NodeBC % 10; // if bc_s == 0 then homogeneous dirichlet

                int bc_wall, bc_inlet, bc_init;
                switch ( bc_s ) {
                case 1:
                    bc_wall = 1;
                    bc_inlet = bc_init = 0;
                    break;
                case 2:
                    bc_inlet = 1;
                    bc_wall = bc_init = 0;
                    break;
                case 4:
                    bc_init = 1;
                    bc_wall = bc_inlet = 0;
                    break;
                default:
                    bc_wall = 1;
                    bc_inlet = bc_init = 0;
                    break;
                }

                double wall_value[2];
                wall_value[0] = wall_value[1] = 1.;
                double init_val[2];
                init_val[0] = init_val[1] = 1.;

                init_val[0] = _data_eq[2].ub[_FF_idx[K_F] * NDOF_FEM + lv_node];
                init_val[1] = _data_eq[2].ub[ ( _FF_idx[K_F]+1 ) * NDOF_FEM + lv_node];

                if ( bc_wall == 1 ) { // DIRICHLET VALUES ON WALL BOUNDARIES
                    double WallDistance = _WallDist;
                    double kappa, omega;

                    if ( _RANS_parameter._WallFunctionApproach == 1 ) {
                        WallDistance = _y_dist;
                    }

                    _mgutils._TurbParameters->DynTurNearWallValues ( kappa, omega, _y_dist, utau );
                    
                    wall_value[0] = ( kappa );
                    wall_value[1] = ( omega );
                }

                double inlet_value[2];
                inlet_value[0] = inlet_value[1] = 1.;

                if ( bc_inlet == 1 ) { // DIRICHLET VALUES ON INLET BOUNDARIES   ->    http://support.esi-cfd.com/esi-users/turb_parameters/
                    double k_in, w_in;
                    _mgutils._TurbParameters->DynTurInitValues ( k_in, w_in, _WallDist, true );
                    inlet_value[1] = w_in ;
                    inlet_value[0] = k_in ;
                }

                _FeM ( lv_node )          = Ipenalty * ( bc_wall * wall_value[_dir] + bc_inlet * inlet_value[_dir] + bc_init * init_val[_dir] );
                _KeM ( lv_node, lv_node ) = Ipenalty;

            }// END DIRICHLET BOUNDARY CONDITION
        }// END LOOP OVER BOUNDARY FACE NODES


        // ======================== end Dirichlet  boundary conditions  ===========================

        // NEUMANN BOUNDARY CONDITIONS ============================================================
        if ( FaceBD / 10 > 0 ) { // bc 10-30
            //   (Neumann DT.n=bc_alpha*T+bc_beta*value)

            double alpha_eff = 100.;
            int  bc_s     = FaceBD % 10;
            int  bc_alpha = ( int ) ( ( bc_s & 2 ) >> 1 ); // (1?) linear term
            int  bc_beta  = ( int ) ( bc_s % 2 ); // (?1)  constant term
            // Non homogenous Neumann boundary conditions  ***********************************

            // GAUSS LOOP
            for ( int qp = 0; qp <  elb_ngauss; qp++ ) {
                // quad/linear  [2]=quad [1]=linear------------------------------------
                double det  = _fe[2]->JacSur ( qp, _xxb_qnds, _InvJac2 ); // local coord _phi_g and jac
                double JxW_g2 = det * _fe[2]->_weight1[ _nTdim - 2][qp]; // weight
                _fe[2]->get_phi_gl_g ( _nTdim - 1, qp, _phi_g[2] );     // global coord _phi_g

                if ( _AxiSym == 1 ) { // axisymmetric  (index ->0)
                    double xyg[DIMENSION];
                    interp_el_bd_sol ( _xx_qnds, sur_toply, elb_ndof2, 0, _nTdim, _phi_g[2], elb_ndof2, xyg );
                    JxW_g2  *= xyg[0];
                }

                for ( int lsi_node = 0; lsi_node < elb_ndof2; lsi_node++ ) { // local side loop (over the node face)
                    // set up row i
                    const double phii_g = _phi_g[2][lsi_node]; // boundary test function
                    const int lei_node = sur_toply[lsi_node]; // local element index

                    if ( _bc_el[lei_node] == 1 ) {
                        const int NodeBD = _bc_vol[lei_node] % 100;

                        if ( NodeBD == FaceBD ) {
                            double wall_der[2];
                            wall_der[0] = wall_der[1] = 0.;

                            if ( bc_alpha == 1 ) {
                                double WallDist = _y_dist;
                                int numod = _mgutils._TurbParameters->_numod;
                                const double k_der = -2. * _IRe / ( WallDist );
                                const double w_der =  2. * _IRe / ( WallDist );
                                const double w_der_log = _IRe / ( WallDist );
                                wall_der[1] = ( yplus < 2. ) ? w_der : w_der_log;
                                wall_der[0] = ( yplus < 15. ) ? k_der : 0.;
                            }

                            int nsides = _NodeOnNwallSides[lei_node];

                            const double FemBC = _ExplicitNearWallDer[_dir];
                            const double KemBC = 1 - _ExplicitNearWallDer[_dir];

                            double AggWallDerVal =  JxW_g2 * phii_g  * ( bc_alpha *  wall_der[_dir] * sign );
                            
                            _FeM ( lei_node ) += FemBC * AggWallDerVal ;

                            // Assemblying Matrix ---------------------------------
                            for ( int lsj_node = 0; lsj_node < elb_ndof2;  lsj_node++ ) {
                                _KeM ( lei_node, sur_toply[lsj_node] ) -= KemBC * AggWallDerVal * _phi_g[2][lsj_node];
                            }// END LOOP OVER BOUNDARY NODES
                        }
                    }// END LOOP OVER BOUNDARY TEST FUNCTIONS
                }
            }// END QUADRATURE LOOP - qp INDEX
        }// END NEUMANN BOUNDARY CONDITIONS ===========================================================
    }

    return;
}




#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
