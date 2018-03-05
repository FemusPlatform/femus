#include "Equations_conf.h"

// ============================================
#ifdef COLOR_EQUATIONS // 3D-2D Energy equation
// ============================================
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverCOL.h"       // Navier-Stokes class header file

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MGFE.h"          // Mesh class,double vel_g[]


// ===============================================================================================
void MGSolCOL::vol_integral (
     DenseMatrixM &KeM, DenseVectorM &FeM,
     const int el_ndof2,const int el_ngauss,
     double xx_qnds[],
     const int unsteady,const int mode
)   // ==============================================================================================
{
     double xyz_g[DIMENSION];
// --------------------------------------------------------------------------------------------------------------------
     /// c) gaussian integration loop (n_gauss)
     // ------------------------------------------------------------------------------------------------------------------
     for ( int qp=0; qp< el_ngauss; qp++ ) {
          // shape functions at gaussian points -----------------------------------
          double det2      = _fe[2]->Jac ( qp,xx_qnds,_InvJac2 );  // Jacobian
          double JxW_g2 =det2*_fe[2]->_weight1[_nTdim-1][qp];       // weight
          _fe[2]->get_phi_gl_g ( _nTdim,qp,_phi_g[2] );            // shape funct
          _fe[2]->get_dphi_gl_g ( _nTdim,qp,_InvJac2,_dphi_g[2] ); // global coord deriv
          _fe[2]->get_ddphi_gl_g ( _nTdim,qp,_InvJac2,_ddphi_g[2] ); // local second deriv

          //  fields --------------------------------------------------------------------------------------------------------
          interp_el_sol ( _xx_qnds,0,_nTdim,_phi_g[2],el_ndof2,xyz_g );
          interp_el_sol ( _data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof2,_ub_g[2] ); // quadratic
          //  derivatives

#ifdef AXISYM   // axisymmetric (index -> 0)
          JxW_g[2]  *=xyz_g[0];
#endif

          double num_visc = _ub_g[2][_FF_idx[CO_F]+1]*_IRe;

          /// d) Local (element) assemblying energy equation
          // =====================================================================================================================
          for ( int i=0; i<el_ndof2; i++ )     {

               if ( _bc_el[i]==1 ) {
                    double dtxJxW_g=JxW_g2;           // area with bc and weight
                    const double phii_g=_phi_g[2][i];

                    // Matrix Assemblying ------------------------------------------------------------------------------------------------
                    for ( int j=0; j<el_ndof2; j++ ) {
                         double Lap=0.;
                         for ( int idim=0; idim<  _nTdim; idim++ ) {
                              const double  dphiidxg=_dphi_g[2][i+idim*el_ndof2];
                              const double  dphijdxg=_dphi_g[2][j+idim*el_ndof2];
                              Lap += dphijdxg*dphiidxg;                                                        // diffusion
                         }

                         const double phij_g= _phi_g[2][j];
                         const double SourceVal = ( 1-_dir ) *_y_dist + _dir*_data_eq[2].ub[j+ ( _FF_idx[CO_F]+1 ) *el_ndof2];

                         FeM ( i )   += dtxJxW_g* ( SourceVal ) *phii_g*phij_g;

                         KeM ( i,j ) += dtxJxW_g* ( // energy-equation matrix
                                            phij_g*phii_g// time term
                                            +Lap* ( 1.e-8* ( 1-_dir ) + 1.e-6*_dir ) //diffusion term
                                       );
                    }
               } // ----------------------------------------
               else {
                    FeM ( i ) = 0.;
                    KeM ( i,i ) = 1.;
               }
          }
     } // end of the quadrature point qp-loop ***********************

     return;
}

#endif
