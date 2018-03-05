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

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
#include "MGSolverP.h"       // Navier-Stokes class header file

// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif



// ================================================================


#if (NS_EQUATIONS%2==0)
#include "Solvertype_enum.h"
// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
// ==================================================================
/// This function constructs the 3d-2D MGSolP class I


MGSolP::MGSolP(
    MGEquationsSystem &mg_equations_map_in,
    const int             nvars_in[],    ///
    std::string     eqname_in,    ///< name equation (in)
    std::string     varname_in    ///< name variable (in)
):  MGSolDA(mg_equations_map_in, nvars_in, eqname_in, varname_in),
//
/*===================================================================================================*/
/*                                    A) Reading parameters                                          */
/*===================================================================================================*/
//
    _offset(_mgmesh._NoNodes[-1]),             // mesh nodes (top level)
    _dt(stod(_mgutils._sim_config["dt"])),                        // parameter  dt
    _uref(_mgutils._mat_prop["Uref"]),         // parameter  u reference
    _lref(_mgutils._mat_prop["Lref"]),         // parameter  l reference
    _rhof(_mgutils._mat_prop["rho0"]),         // parameter density
    _muf(_mgutils._mat_prop["mu0"])           // parameter viscosity
{
    // parameter viscosity

    _nPdim = DIMENSION;

    // class equation ---------------------------------------------------------------------
    for (int k_index = 0; k_index < 30; k_index++) {
        _FF_idx[k_index] = -1;
    }
    /*===================================================================================================*/
    /*                                    B) Setting class variables                                     */
    /*===================================================================================================*/
    _var_names[0] = "p";
    _refvalue[0] = _rhof * _uref * _uref; // class variable names
    /*===================================================================================================*/
    /*                                    C) Setting solver type                                         */
    /*===================================================================================================*/
    for (int  l = 0; l < _NoLevels; l++) {
        _solver[l]->set_solver_type(GMRESM);
    }
    /*===================================================================================================*/
    /*                                    D) Setting no _nNSdimensional parameters                       */
    /*===================================================================================================*/

//   _SolveP = true;

    return;
}
//
//  ====================================================
//  This function assembles the matrix and the rhs:
//  ====================================================
//
void  MGSolP::GenMatRhs(
    const double /* time*/, // time  <-
    const int  Level,  // Level <-
    const  int mode    // mode  <- (1=rhs) (0=only matrix)
)    // ===============================================
{

    /*===================================================================================================*/
    /*                                    A) Set up                                                      */
    /*===================================================================================================*/
    /*----------------------------------- Geometry  -----------------------------------------------------*/
    int axis = (int)(_mgutils._geometry["Axisym"]);
    const int   offset = _mgmesh._NoNodes[_NoLevels - 1];				// mesh nodes
//   double      xx_qnds[DIMENSION*NDOF_FEM];           				// element node coords
    //   int         el_conn[NDOF_FEM];                    			// element connectivity
    const int  el_sides = _mgmesh._GeomEl._n_sides[0];              	       // element nodes
    int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];             	        // element connectivity
    int        el_neigh[NDOF_FEM];                                 	        // element connectivity
    double     uvw_dxg[DIMENSION];
//   double     u_old_p[NDOF_FEM];double u_nl_p[NDOF_FEM];

    double p_old[NDOF_FEM];
    double dp[NDOF_FEM];

    /*----------------------------------- Gauss integration  --------------------------------------------*/

    double det[3], JxW_g[3], InvJac[3][DIMENSION * DIMENSION];           		// Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION]; 		// global derivatives at g point
    const int  el_ngauss = _fe[2]->_NoGauss1[ _nPdim - 1];                		//elem gauss points
    double u_nl[DIMENSION * NDOF_FEM];
    double u_old[DIMENSION * NDOF_FEM];

    double vel_gdx[DIMENSION * DIMENSION];

    // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
    int el_ndof[3];
    el_ndof[0] = 1;                // number of el dofs
    int el_mat_nrows = 0;                           // number of mat rows (dofs)
    for (int ideg = 1; ideg < 3; ideg++) {
        el_ndof[ideg] = _fe[ideg]->_NoShape[ _nPdim - 1];
        el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
    };
    int el_mat_ncols = el_mat_nrows;                     //square matrix
    std::vector<int > el_dof_indices(el_mat_ncols);      // element dof vector
//
    /*---------------------------- Fields -> Navier-Stokes  [NS_F]  -------------------------------------*/
//
//   int idx_ns= _data_eq[2].tab_eqs[NS_F];// Navier-Stokes ---------------------------------
//   _NS_idx=(idx_ns>=0)? _data_eq[2].indx_ub[idx_ns]:-1;
    for (int k = 0; k < 30; k++) { // coupling  basic system fields
        const int idx = _data_eq[2].tab_eqs[k];
        _FF_idx[k] = (idx >= 0) ? _data_eq[2].indx_ub[idx] : -1;
    }
//
    /* --------------------------- Element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --*/
//
    A[Level]->zero();
    if (mode == 1) {
        b[Level]->zero();    // global matrix+rhs
    }
    DenseMatrixM KeM;
    DenseVectorM FeM;                // local  matrix+rhs
    KeM.resize(el_mat_nrows, el_mat_ncols);
    FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs
//
    int ndof_lev = 0;
    for (int pr = 0; pr < _mgmesh._iproc; pr++) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
        ndof_lev += delta;
    }
//
    const int  nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
    const int  nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element
    for (int iel = 0; iel < (nel_e - nel_b); iel++) {

        // ===================================================================================================
        //                          B) Element  Loop over the volume (n_elem)
        // ===================================================================================================

        // ------------------------ Set to zero matrix and rhs and center -----------------------------------
        KeM.zero();
        FeM.zero();
        // ------------------------ Geometry and element fields ---------------------------------------------
        // ------------------------ Element Connectivity (el_conn) and coordinates (xx_qnds) ----------------
        _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, _xx_qnds);
        _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);
        get_el_dof_bc(Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);
        // element nodes coordinates
        for (int idim = 0; idim < _nPdim; idim++) {
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + idim]]->get_el_oldsol(0, 1, el_ndof[2], el_conn, offset, idim, u_old);
        }
        _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_sol(0, 1, el_ndof[1], el_conn, offset, 0, p_old);
        _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0, 1, el_ndof[1], el_conn, offset, 0, dp);

        for (int  j = 0; j < el_ndof[1]; j++) {
            _bc_el[j] = _bc_vol[j] / 10;
        }
        // ===================================================================================================
        //                          C) Gaussian integration loop (n_gauss)
        // ===================================================================================================

        for (int  qp = 0; qp < el_ngauss; qp++) {
            // ------------------------- Shape functions at gaussian points --------------------------------------
            for (int  ideg = 1; ideg < 3; ideg++) {  				// linear-quadratic  [1]=linear [2]=quad
                det[ideg]      = _fe[ideg]->Jac(qp, _xx_qnds, InvJac[ideg]);    	// Jacobian
                JxW_g[ideg] = det[ideg] * _fe[ideg]->_weight1[ _nPdim - 1][qp];      	// weight
                _fe[ideg]->get_phi_gl_g(_nPdim, qp, _phi_g[ideg]);                	// shape funct
                _fe[ideg]->get_dphi_gl_g(_nPdim, qp, InvJac[ideg], _dphi_g[ideg]); 	// global coord deriv
            }

            interp_el_gdx(u_old, 0, _nPdim, _dphi_g[2], el_ndof[2], vel_gdx); // derivatives  vel_gdx[DIM][DIM]
            interp_el_sol(_xx_qnds, 0, _nPdim, _phi_g[2], el_ndof[2], _xyzg);

            if (axis == 1) {
                JxW_g[2]  *= _xyzg[0];     //  printf("  %f \n",_ub_g[2][0]);
                JxW_g[1]  *= _xyzg[0];
            }

            // ===================================================================================================
            //                          D) Local (element) assemblying pressure equation
            // ===================================================================================================
            for (int  i = 0; i < el_ndof[1]; i++)     {
                const double phii_g = _phi_g[1][i];
                for (int  idim = 0; idim <  _nPdim; idim++) {
                    dphiidx_g[1][idim] = _dphi_g[1][i + idim * el_ndof[1]];
                }
                // =================================================
                //  enum bound_cond_p {outflowp0=0,outflowp=4,vel_fix=10,interiorp=11}
                // =================================================
                if (_bc_el[i] == 1) {
                    double dtxJxW_g = JxW_g[2];
                    // ------------------------- Rhs Assemblying  -------------------------------------------------------
                    double div = 0.;
                    for (int  idim = 0; idim < _nPdim; idim++) {
                        div += vel_gdx[idim * _nPdim + idim];    // velocity divergence
                    }
                    FeM(i) -= dtxJxW_g * div * phii_g / _dt;    // velocity divergence
                    
                    if (axis == 1)
                        for (int  j = 0; j < el_ndof[2]; j++) {
                            FeM(i) += dtxJxW_g * (u_old[j]) * _phi_g[2][j] * phii_g / (_xyzg[0]*_dt) ;    // axisym
                        }

                    // ------------------------- Matrix Assemblying  ----------------------------

                    for (int  j = 0; j < el_ndof[1]; j++) {
                        const double phij_g = _phi_g[1][j];
                        double Lap = 0.;
                        for (int  idim = 0; idim < _nPdim; idim++) {
                            if (axis == 1) {
                                Lap += (1 - (idim)) * phii_g * phij_g / (_xyzg[0] * _xyzg[0]);    // axysimmetry
                            }
                            dphijdx_g[1][idim] = _dphi_g[1][j + idim * el_ndof[1]];
                            Lap += dphijdx_g[1][idim] * dphiidx_g[1][idim]; // Laplacian
                        }
                        // Pressure matrix assembling ---------------------------------------------
                        KeM(i, j) += dtxJxW_g * Lap;
                    } // ---------------------------------------------
                } else if (_bc_el[i] == 0) {
                    if (_bc_vol[i] == 4)  {
                        int press_ind = _bc_vol[i] % 3;
                        FeM(i) = press_ind * (p_old[i]);// * _dt;
                    }
                    KeM(i, i) = 1.;
                }
            }
        } // end of the quadrature point qp-loop +++++++++++++++++++++++++

        // ===================================================================================================
        //                          E) Global assemblying pressure equation
        // ===================================================================================================
//    std::cout << "\n renormilized KeM" << KeM << "\n FeM" << FeM;
        A[Level]->add_matrix(KeM, el_dof_indices);                 // global matrix
        if (mode == 1) {
            b[Level]->add_vector(FeM, el_dof_indices);    // global rhs
        }

    } // end of element loop
    // clean and close
    el_dof_indices.clear();
    A[Level]->close();
    if (mode == 1) {
        b[Level]->close();
    }

#ifdef PRINT_INFO
    std::cout << " Matrix Assembled(P)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

    return;
} /******************************************************************************************************/
//
//

// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolP::MGTimeStep(
    const double time,  // time
    const int /*iter*/  // Number of max inter
)
{
    /* ========================================================================================= */
    /*              A) Set up the time step                                                      */
    /* ========================================================================================= */
    std::cout  << std::endl << "\033[038;5;" << 155 << ";1m "
               << "--------------------------------------------------- \n\t"
               <<  _eqname.c_str()
               << " solution of problem " << _mgutils.get_name()
               << "\n ---------------------------------------------------\n \033[0m";
    /* ========================================================================================= */
    /*              B) Assemblying of the Matrix-Rhs                                             */
    /* ========================================================================================= */
#if PRINT_TIME==1
    std::clock_t start_time = std::clock();
#endif
    GenMatRhs(time, _NoLevels - 1, 1);                                            // matrix and rhs
    for (int Level = 0 ; Level < _NoLevels - 1; Level++) {
        GenMatRhs(time, Level, 0);    // matrix
    }
#if PRINT_TIME==1
    std::clock_t end_time = std::clock();
    std::cout << "  Assembly time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
    /* ========================================================================================= */
    /*               C) Solution of the linear MGsystem (MGSolP::MGSolve)                        */
    /* ========================================================================================= */
    MGSolve(1.e-6, 50);
#if PRINT_TIME==1
    end_time = std::clock();
    std::cout << " Assembly+solution time -----> =" << double(end_time - start_time) / CLOCKS_PER_SEC
              << "s " << std::endl;
#endif
    /* ========================================================================================= */
    /*               D) Update of the old solution at the top Level		                     */
    /* ========================================================================================= */
//     x[_NoLevels - 1]->localize(*x_oold[_NoLevels - 1]); // dp
//     x_old[_NoLevels - 1]->add(1. / _dt, *x_oold[_NoLevels - 1]); // p


      x[_NoLevels - 1]->localize(*x_old[_NoLevels - 1]); // dp
    return;
} /******************************************************************************************************/
//
//

// /// =========================================================================================
// /// This function controls the assembly and the solution of the P_equation system:
// /// =========================================================================================
// void MGSolP::MGTimeStep_nl_setup(
//   const double time,  // time
//   const int /*iter*/  // Number of max inter
// ) {
// // =========================================================================================
//
//   return;
// } /******************************************************************************************************/
// //
// //

/// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
/// =========================================================================================

// void MGSolP::MGTimeStep_nl_sol_up(
//   const double time,  // time
//   const int /*iter*/  // Number of max inter
// ) {
//
//   /* ========================================================================================= */
//   /*              A) Set up the time step                                                      */
//   /* ========================================================================================= */
// //
//   std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
//
// //
//   /* ========================================================================================= */
//   /*              B) Assemblying of the Matrix-Rhs                                             */
//   /* ========================================================================================= */
// //
// #if PRINT_TIME==1
//   std::clock_t start_time=std::clock();
// #endif
//   GenMatRhs(time,_NoLevels-1,1);                                               // matrix and rhs
//   for(int Level = 0 ; Level < _NoLevels-1; Level++) { GenMatRhs(time,Level,0); } // matrix
// #if PRINT_TIME==1
//   std::clock_t end_time=std::clock();
//   std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
// #endif
// //
//   /* ========================================================================================= */
//   /*               C) Solution of the linear MGsystem (MGSolP::MGSolve)                        */
//   /* ========================================================================================= */
// //
//   MGSolve(1.e-6,50);
// #if PRINT_TIME==1
//   end_time=std::clock();
//   std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
//             << "s "<< std::endl;
// #endif
// //
//   /* ========================================================================================= */
//   /*               D) Update of the old solution at the top Level		                     */
//   /* ========================================================================================= */
// //
//   x[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);  // dp
//   x_old[_NoLevels-1]->add(1.,*x_oold[_NoLevels-1]);// p
//   return;
// }// =======================================================================================
//
// /// ======================================================
// /// This function controls the time step operations:
// /// ======================================================
// int MGSolP::MGTimeStep_nl_iter(const double time, int) {
//
//   return 0;
// } /******************************************************************************************************/
// //
// //


#endif  //ENDIF NS_EQUATIONS%2==0
// // *************************************************

#endif  //ENDIF NS_EQUATIONS
// #endif  // NS_equation is personal
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 

