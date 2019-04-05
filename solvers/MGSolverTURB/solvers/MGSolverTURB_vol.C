#include "Equations_conf.h"

// ============================================
#ifdef _TURBULENCE_ // 3D-2D Energy equation
// ============================================

#include "MGSolverTURB.h"       // Navier-Stokes class header file
#include "MGUtils.h"
#include "TurbUtils.h"
#include "MGMesh.h"
#include "MGFE.h"          // Mesh class,double vel_g[]


// ===============================================================================================
void MGSolTURB::vol_integral (
    DenseMatrixM & KeM,
    DenseVectorM & FeM,
    const int el_ndof2,
    const int el_ngauss,
    double xx_qnds[],
    const int unsteady,
    const int mode,
    int el_conn[]
)   // ==============================================================================================
{

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------
    double alpha = 1.e-10;
    double wall_dist[1], dynTur[2], thermTur[2];

    int WallElement = 0;

    for ( int i = 0; i < el_ndof2; i++ )
        if ( _bc_el[i] == 0 ) {
            WallElement = 1;
        }

    for ( int qp = 0; qp < el_ngauss; qp++ ) {
        // shape functions at gaussian points -----------------------------------
        double det2      = _fe[2]->Jac ( qp, xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];  // weight
        _fe[2]->get_phi_gl_g ( _nTdim, qp, _phi_g[2] );          // shape funct
        _fe[2]->get_dphi_gl_g ( _nTdim, qp, _InvJac2, _dphi_g[2] ); // global coord deriv
        _fe[2]->get_ddphi_gl_g ( _nTdim, qp, _InvJac2, _ddphi_g[2] ); // local second deriv

        //  fields --------------------------------------------------------------------------------------------------------
        interp_el_sol ( _KW, 0, 2, _phi_g[2], el_ndof2, dynTur ); // quadratic
        interp_el_sol ( _KHWH, 0, 2, _phi_g[2], el_ndof2, thermTur ); // quadratic
        interp_el_sol ( _DIST, 0, 1, _phi_g[2], el_ndof2, wall_dist ); // quadratic


#ifdef AXISYM   // axisymmetric (index -> 0)
        double xyz_g[DIMENSION];
        interp_el_sol ( _xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g );
        JxW_g[2]  *= xyz_g[0];
#endif
        double utau = 0.;
        double yplus = 0.;

        if ( WallElement == 1 ) { // CALCULATION OF UTAU AND UNIT NORMAL VECTOR
            double vel_bound;
            vel_bound = _data_eq[2].ub[ ( _FF_idx[NS_F] ) * NDOF_FEM + NDOF_FEM - 1];
            double WallDist = _DIST[NDOF_FEM - 1];
            utau  = _mgutils._TurbParameters->CalcUtau ( vel_bound, WallDist );
            yplus = WallDist * utau / _IRe;
        }

        /// d) Local (element) assemblying energy equation
        // =====================================================================================================================
        for ( int i = 0; i < el_ndof2; i++ )     {

            if ( _bc_el[i] == 1 ) {
                double dtxJxW_g = JxW_g2;         // area with bc and weight
                const double phii_g = _phi_g[2][i];
                double source = 0;
                double mu_turb    = _mgutils._TurbParameters->CalcMuTurb ( dynTur, wall_dist[0] );
                double alpha_turb = _mgutils._TurbParameters->CalcAlphaTurb ( dynTur, thermTur, wall_dist[0] );

                if ( _dir == 0 ) {
                    source = _y_dist + _mgutils._geometry["Wall_dist"];
                }

                if ( _dir == 1 ) {
                    source = mu_turb;
                }

                if ( _dir == 2 ) {
                    source = alpha_turb;
                }

                FeM ( i )   += dtxJxW_g * phii_g * ( source );
                _res += dtxJxW_g * phii_g * ( source );

                // Matrix Assemblying ------------------------------------------------------------------------------------------------
                for ( int j = 0; j < el_ndof2; j++ ) {
                    double Lap = 0.;

                    for ( int idim = 0; idim <  _nTdim; idim++ ) {
                        Lap += _dphi_g[2][j + idim * el_ndof2] * _dphi_g[2][i + idim * el_ndof2];
                    }

                    const double phij_g = _phi_g[2][j];
                    KeM ( i, j ) += dtxJxW_g * ( // energy-equation matrix
                                        phij_g * phii_g                          // time term
                                        + Lap * _DiffusionCoefficient//diffusion term
                                    );
                    _res += -1.*dtxJxW_g * ( // energy-equation matrix
                                phij_g * phii_g                          // time term
                                + Lap * _DiffusionCoefficient//diffusion term
                            ) * _data_eq[2].ub[ ( _FF_idx[DIST] + _dir ) * NDOF_FEM + j];
                }
            } // ----------------------------------------
            else {
                FeM ( i ) = 0.;
                KeM ( i, i ) = 1.;
            }

            if ( _dir == 0 &&  _bc_el[i] == 0 ) {
                FeM ( i ) = _mgutils._geometry["Wall_dist"];

                for ( int j = 0; j < el_ndof2; j++ ) {
                    KeM ( i, j ) = 0.;
                }

                KeM ( i, i ) = 1.;
            }
        }
    } // end of the quadrature point qp-loop ***********************

    return;
}


void MGSolTURB::compute_y_plus (
    const int el_ndof2,
    const int el_ngauss,
    double xx_qnds[],
    const int unsteady,
    const int mode,
    int el_conn[],
    int iel
)   // ==============================================================================================
{
    double grad[DIMENSION] = {0., 0.};
    grad[DIMENSION -1] = 0.;
    double area = 0.;
    for ( int qp = 0; qp < el_ngauss; qp++ ) {
        double det2      = _fe[2]->Jac ( qp, xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2    = det2 * _fe[2]->_weight1[_nTdim - 1][qp];   // weight
        area += JxW_g2;

        _fe[2]->get_dphi_gl_g ( _nTdim, qp, _InvJac2, _dphi_g[2] ); // global coord deriv

        double der_dist[DIMENSION] = {0., 0.};
        der_dist[DIMENSION -1] = 0.;


        for ( int n = 0; n < el_ndof2; n++ ) {
            for ( int dim = 0; dim<DIMENSION; dim++ ) {
                der_dist[dim] += JxW_g2 * _data_eq[2].ub[_FF_idx[DIST]*el_ndof2 + n] * _dphi_g[2][n + dim*el_ndof2];
            }
        }
        for ( int dim = 0; dim<DIMENSION; dim++ ) {
            grad[dim] += JxW_g2 * der_dist[dim];
        }
    }
    double norm[DIMENSION];
    double mod = 0;
    for ( int dim = 0; dim<DIMENSION; dim++ ) {
        mod += grad[dim]*grad[dim];
    }
    mod = sqrt ( mod ) + 1.e-10;
    for ( int dim = 0; dim<DIMENSION; dim++ ) {
        norm[dim] = -grad[dim]/mod;
    }


    // calculation of utau from mid cell point values
    double vel_mod =0., vel_norm = 0., vel_bound=0.;

    for ( int dim = 0; dim < _nTdim; dim ++ ) {
        double vel_d = _data_eq[2].ub[ ( _FF_idx[NS_F] + dim ) * NDOF_FEM + NDOF_FEM - 1];
        vel_mod  += vel_d * vel_d;
        vel_norm += vel_d * norm[dim];
    }

    vel_bound = sqrt ( vel_mod - vel_norm * vel_norm ) + 1.e-10;

    double WallDist = _data_eq[2].ub[_FF_idx[DIST]*el_ndof2 + el_ndof2 -1];
    double utau  = _mgutils._TurbParameters->CalcUtau ( vel_bound, WallDist );
    
    double yplus = WallDist * utau / _IRe;
    
    _mgmesh._yplus[iel] = yplus;
            
    return;
}

// ===============================================================================================
void MGSolTURB::rhs_integral (
    DenseVectorM & FeM,
    const int el_ndof2, const int el_ngauss,
    double xx_qnds[],
    const int unsteady, const int mode, int el_conn[]
)   // ==============================================================================================
{

    // --------------------------------------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------------------------------------
    double wall_dist[1], dynTur[2], thermTur[2];

    for ( int qp = 0; qp < el_ngauss; qp++ ) {
        // shape functions at gaussian points -----------------------------------
        double det2      = _fe[2]->Jac ( qp, xx_qnds, _InvJac2 ); // Jacobian
        double JxW_g2 = det2 * _fe[2]->_weight1[_nTdim - 1][qp];  // weight
        _fe[2]->get_phi_gl_g ( _nTdim, qp, _phi_g[2] );          // shape funct

        //  fields --------------------------------------------------------------------------------------------------------
        interp_el_sol ( _KW, 0, 2, _phi_g[2], el_ndof2, dynTur ); // quadratic
        interp_el_sol ( _KHWH, 0, 2, _phi_g[2], el_ndof2, thermTur ); // quadratic
        interp_el_sol ( _DIST, 0, 1, _phi_g[2], el_ndof2, wall_dist ); // quadratic

#ifdef AXISYM   // axisymmetric (index -> 0)
        double xyz_g[DIMENSION];
        interp_el_sol ( _xx_qnds, 0, _nTdim, _phi_g[2], el_ndof2, xyz_g );
        JxW_g[2]  *= xyz_g[0];
#endif

        /// d) Local (element) assemblying energy equation
        // =====================================================================================================================
        for ( int i = 0; i < el_ndof2; i++ )     {

            if ( _bc_el[i] == 1 ) {
                double source = 0;
                double mu_turb    = _mgutils._TurbParameters->CalcMuTurb ( dynTur, wall_dist[0] );
                double alpha_turb = _mgutils._TurbParameters->CalcAlphaTurb ( dynTur, thermTur, wall_dist[0] );

                if ( _dir == 0 ) {
                    source = _y_dist;
                }

                if ( _dir == 1 ) {
                    source = mu_turb;
                }

                if ( _dir == 2 ) {
                    source = alpha_turb;
                }

                FeM ( i )   += JxW_g2 * _phi_g[2][i] * source ;

            } // ----------------------------------------
            else {
                FeM ( i ) = 0.;
            }
        }
    } // end of the quadrature point qp-loop ***********************

    return;
}


#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
