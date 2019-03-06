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
#if NS_EQUATIONS==2
// local Femus class include -----------------------------------
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file
#include "MGFE.h"          // Mesh class
#include <iomanip>      // std::setprecision

void MGSolNS::set_bc_matrix(
    DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
    int dir_maxnormal,                          ///<  normal dir
    int sur_toply[],                            ///< boundary topology map
    int el_ndof[],                              ///< number of volume dofs
    int elb_ndof[],                             ///< number of boundary dofs
    int elb_ngauss,                             ///<  number of surface gaussian points
    double normal[],                            ///< normal
    double u_old[],
    int el_conn[]
)
{
    double ProjDir[3][DIMENSION];
    double ProjDir_Dirichlet[3][DIMENSION];
    for (int dir = 0; dir < _nNSdim; dir++) {
        ProjDir_Dirichlet[2][dir] = normal[dir];
        ProjDir[2][dir] = normal[dir];
    }

    const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;
    const int norm_face = (abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 10000) / 1000 - 1 ;

    // TANGENT VECTORS BASED ON BOUNDARY GEOMETRY
    double dir_tang1[DIMENSION], dir_tang2[DIMENSION];
    CalcTangDir(ProjDir_Dirichlet[0], ProjDir_Dirichlet[1], ProjDir_Dirichlet[2], GEOM_BASED);  // for dirichlet bc
    CalcTangDir(ProjDir[0], ProjDir[1], ProjDir[2], VEL_BASED);// for neumann bc
    
    SetBCFlags(normal, sur_toply, el_ndof[2], el_conn);
    double u_oold[DIMENSION * NDOF_FEM];
    for (int kdim = 0; kdim < _nNSdim; kdim++) {
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F + kdim]]->get_el_oldsol(0, 1, el_ndof[2], el_conn, _offset, kdim, u_oold);
    }

    if (bd_face == 88 || bd_face == 84) {
        if (_FF_idx[K_F] >= 0) {
            if (_NS_parameter._WallFunctionApproach == 1) {
                _WallElement = 1;
            }
            for (int dir = 0; dir < _nNSdim; dir++) {
                _normal_pt[(NDOF_FEM - 1) *_nNSdim + dir] += normal[dir];
            }
        }
    }
        
    // ...........................................................................................
    //     DIRICHLET BOUNDARY CONDITION
    // ...........................................................................................

    for (int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {

        const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;
        const int bc_var_check2 = abs(_bc_bd[sur_toply[lbnode]]);    // total  bc_var
        const int bc_var_check = bc_var_check2 % 10000;         // bc_var
        const int bc_var_normal = bc_var_check / 1000 - 1;      //
        const int bc_var = bc_var_check % 1000;
        const int bc_n_flag = bc_var / 10;
        const int bc_tg_flag = bc_var % 10;
        const int norm_node = (abs(_bc_bd[sur_toply[lbnode]]) % 10000) / 1000 - 1 ;

        if (_nNSdim == 3) {
            CorrectBCFlags(ProjDir_Dirichlet[0], ProjDir_Dirichlet[1], bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
        }

        if (bd_face == bd_node) {
            const int  indx_row = sur_toply[lbnode] + _dir * el_ndof[2]; // solution indx
            const int  node = sur_toply[lbnode];

            bool NormBC = (abs(_bc_el[indx_row]) == 3) ? true : false;
            bool DirichletBC = (_bc_el[indx_row] > 0) ? true : false;

            double u_old_node[3];
            for (int jvar = 0; jvar < _nNSdim; jvar++) {
                u_old_node[jvar] =  u_old[node + jvar * NDOF_FEM];
            }

            if (DirichletBC && norm_face == norm_node/*!_AlreadyWrittenDirBC[indx_row]*/) {
                int  bc_bc  = NormBC ? bc_n_flag : bc_tg_flag;
                const int  bc_rhs = ((bc_bc % 2) == 1) ? 1 : 0; // bc_rhs -> rhs bc
                bool OrDir = (bc_bc > 7) ? true : false;
                const double vart = bc_rhs * u_old[node + _dir * el_ndof[2]];

                double prod = 0., u_norm[DIMENSION];
                for (int i = 0; i < _nNSdim; i++) {
                    prod += u_old_node[i] * normal[i];
                }
                for (int i = 0; i < _nNSdim; i++) {
                    u_norm[i] = prod * normal[i];
                }

                double alpha_beta = 1.;
                if (bd_node == 88 || bd_node == 66) {
                    KeM(node, node) = 1.;
                    FeM(node)      = 0.;
                } else if (OrDir) {
                    //  IMPOSITION OF VELOCITY VECTOR ALONG GENERAL DIRECTION f
                    //  \vec(U)_f = (\vec(U)_{old} \dot \hat(f)) \hat(f)
                    //  f can be normal or tangential direction
                    if (bd_node == 98) {
                        KeM(node, node) = 1.;
                        FeM(node)       = u_norm[_dir];
                    } else {             
                        int direction = -2*(abs(_bc_el[indx_row])/3) + 2;                     
                        double scal_prod = 0.; double vel_proj[DIMENSION];
                        for (int i = 0; i < _nNSdim; i++) {
                                scal_prod += u_old_node[i] * ProjDir_Dirichlet[direction][i];
                            }
                            for (int i = 0; i < _nNSdim; i++) {
                                vel_proj[i] = scal_prod * ProjDir_Dirichlet[direction][i];
                            }
                        KeM(node, node)  = 1.;
                        FeM(node)        = vel_proj[_dir];
                    }
                } else {
                    KeM(node, node) = alpha_beta;
                    FeM(node) = vart * alpha_beta;
                }// beta_0 u -> diagonal
            }// END IF DIRICHLET BC
        }// END IF BD_NODE == BD_FACE
    }// END LOOP OVER BOUNDARY NODES - LBNODE

    // ...........................................................................................
    //     INTEGRATION OVER BOUNDARY ELEMENT - STRESS CALCULATION
    // ...........................................................................................
    for (int  qp = 0; qp < elb_ngauss; qp++) { // START LOOP OVER GAUSS POINTS ===============================================
        double det   = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);   // local coord _phi_g and jac
        double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp]; // weight
        
        _fe[2]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[2]);   // global coord _phi_g
        _fe[1]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[1]);   // global coord _phi_g

        if (_AxiSym == 1) {
            interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg);
            JxW_g  *= _xyzg[0];   // axisymmetric  (index ->0)
        }
        const double dtJxW_g = JxW_g;        
        
        double VelOnGauss[DIMENSION];
        interp_el_bd_sol(u_old, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], VelOnGauss);

        for (int i = 0; i < elb_ndof[2]; i++) { // LOOP OVER TEST FUNCTIONS ===================================================

            const int bc_var_check  = abs(_bc_bd[sur_toply[i]]) % 10000;                         // bc_var
            const int bc_var_normal = bc_var_check / 1000 - 1;                            //
            const int bd_node       = abs(_bc_bd[sur_toply[i]]) % 100;
            const int indx          = sur_toply[i] + _dir * el_ndof[2];
            const int NComp         = (_bc_el[indx] < 0) ? _nNSdim : 1;
            const int bc_flag       = (_bc_el[indx] == -3) ? bd_node / 10 : bd_node % 10;  // normal or tangential bc flags
            const int pen_cond      = (bd_node == 44) ? 0 : 1;
            const double lambda1    = (1 - pen_cond) * _NS_parameter._Penalty_tg;
            const double beta       = lambda1 * _NS_parameter._Penalty_n;

            if (bd_face == bd_node) {
                // Direction of maximum vector component of tangent 1 and 2 -> alignment of equations
                if (_nNSdim == 3) {
                    CorrectBCFlags(ProjDir[0], ProjDir[1], bc_var_normal, sur_toply[i], el_ndof[2], 1);
                }
                if (bc_flag < 6 && bd_node != 11) { // only for neumann bc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // NORMAL VECTOR FOR THE PROJECTION OF VOLUME PART EQUATION
                    // _normal_pt is used to project the equation (volume part) using the same normal used in the boundary
                    if (qp == 0) for (int dir = 0; dir < _nNSdim; dir++) {
                            _normal_pt[ sur_toply[i]*_nNSdim + dir] += normal[dir];
                        }
                    for (int  ivarN = _dir; ivarN < _dir + NComp; ivarN++) { // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
                        const int ivar = ivarN % DIMENSION;
                        // Stress component along tangential directions (-1,-2) ---------------------------------------------------
                        if (_bc_el[indx] == -1 || _bc_el[indx] == -2) {
                            if (bc_flag != 1) {
//                                    _                              _
//                                   /   =>   ->    -->             /    ->   -->
//                                  /   (T *  n) *  phi_t ds  =    /  a (u *  phi_t) ds
//                                -s                             -s
                                double alpha_beta = ProjDir[abs(_bc_el[indx]) - 1][ivar];
                                if (bc_flag == 5) {
                                    alpha_beta = ProjDir[0][ivar];
                                }
                                double tau = (_bc_el[indx] == -1) ? (_NS_parameter._Tg1_stress) : _NS_parameter._Tg2_stress ;    //utau*utau/vel_bound;
                                double stress = alpha_beta * (pen_cond * tau - lambda1) * dtJxW_g * _phi_g[2][i];

                                if (ivarN == _dir)
                                    for (int j = 0; j < elb_ndof[2]; j++) {
                                        KeM(sur_toply[i], sur_toply[j]) -= stress * _phi_g[2][j];
                                    }
                                else {
                                    FeM(sur_toply[i]) += stress * VelOnGauss[ivar];
                                }
                            }
                        }// End loop if equation is for tangential direction-----------------------------------------------------
                        else if (_bc_el[indx] == -3) {
                            // normal with penalty (u dot n) n dot phi
//                                    _
//                                   /   ->   ->  ->   -->
//                                  /   (u *  n)  n *  phi_n ds = 0
//                                -s
                            if (bc_flag == 4) {
                                double alpha_beta = 1.;
                                alpha_beta = normal[_dir];
                                double norm_stress = dtJxW_g * alpha_beta * beta * _phi_g[2][i] * normal[_dir];
                                if (ivarN == _dir)
                                    for (int j = 0; j < elb_ndof[2]; j++) {
                                        KeM(sur_toply[i], sur_toply[j]) -= norm_stress * normal[ivar] * _phi_g[2][j];
                                    }
                                else {
                                    FeM(sur_toply[i]) += norm_stress * normal[ivar] * VelOnGauss[ivar];
                                }
                            } //  ---------------------------------------------------------------------------------------------------
                        }// contribution along normal direction
                    }// END LOOP IF BC FLAG IS FOR INTEGRATION
                }// END LOOP OVER VARIABLES FOR LINEAR COMBINATION OF EQUATIONS ++++++++++++++++++++++++++++++++++++++++++++++
            }// END IF BC_NODE == BC_FACE
        }// END LOOP OVER BOUNDARY TEST FUNCTIONS =========================================================================
    }// END GAUSSIAN INTEGRATION ========================================================================================
    return;
}

void MGSolNS::set_bc_matrix_corr_step(
    DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
    int dir_maxnormal,                          ///<  normal dir
    int sur_toply[],                            ///< boundary topology map
    int el_ndof[],                              ///< number of volume dofs
    int elb_ndof[],                             ///< number of boundary dofs
    int elb_ngauss,                             ///<  number of surface gaussian points
    double normal[],                            ///< normal
    double u_nl[],
    double pressure[],
    double pressure_old[],
    int el_conn[]
)
{
    double dir_tang1[DIMENSION], dir_tang2[DIMENSION];
    CalcTangDir(dir_tang1, dir_tang2, normal, GEOM_BASED);
    SetBCFlags(normal, sur_toply, el_ndof[2], el_conn);
    const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;
    // ...........................................................................................
    //     DIRICHLET BOUNDARY CONDITION
    // ...........................................................................................
    for (int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {
        const int bc_var_check2 = abs(_bc_bd[sur_toply[lbnode]]);    // total  bc_var
        const int bc_var_npt = bc_var_check2 / 10000;                // recurrent  points in one element
        const int bc_var_check = bc_var_check2 % 10000;              // bc_var
        const int bc_var_normal = bc_var_check / 1000 - 1;           //
        const int bc_var = bc_var_check % 1000;
        const int bc_n_flag = bc_var / 10;
        const int bc_tg_flag = bc_var % 10;
        if (_nNSdim == 3) {
            CorrectBCFlags(dir_tang1, dir_tang2, bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
        }
        const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;

        if (bd_face == bd_node) {
            const int indx_row = sur_toply[lbnode] + _dir * el_ndof[2]; // solution indx
            const int node = sur_toply[lbnode];
            if (_bc_el[indx_row] > 0) {
                double value = _data_eq[2].ub[(_FF_idx[NS_F] + _dir) * NDOF_FEM + node];
                KeM(node, node) = 1.;
                FeM(node)       = value;
            } else if (_bc_el[indx_row] < 0 && _bc_el[indx_row] > -3) {
                for (int dir = 0; dir < _nNSdim; dir++) {
                    _normal_pt[ node * _nNSdim + dir] = normal[dir];
                }
            } else if (_bc_el[indx_row] == -3 || bd_node == outflow) {// _bc_el == -3 or 0
                for (int  qp = 0; qp < elb_ngauss; qp++) { // START LOOP OVER GAUSS POINTS ===============================================
                    double det   = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);  // local coord _phi_g and jac
                    double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp]/**dt*/; // weight
                    _fe[2]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[2]);   // global coord _phi_g
                    _fe[1]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[1]);   // global coord _phi_g
                    if (_AxiSym == 1) {
                        interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg);
                        JxW_g  *= _xyzg[0];   // axisymmetric  (index ->0)
                    }
                    const int bc_flag = bd_node / 10 ;
                    if (qp == 0) for (int dir = 0; dir < _nNSdim; dir++) {
                            _normal_pt[ sur_toply[lbnode]*_nNSdim + dir] = normal[dir];
                        }
                    int ncomp = (_bc_el[indx_row] == -3) ? _nNSdim : 1;
                    for (int  ivarN = _dir; ivarN < _dir + ncomp; ivarN++) { // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
                        int ivar = ivarN % _nNSdim;
                        double alpha_beta = (_bc_el[indx_row] == -3) ? normal[ivar] : 1;
                        // surface integral p phi dot n   --> pressure_inlet - pressure_outlet - outflow_p ------------------------------------------------
                        for (int j = 0; j < elb_ndof[1]; j++) {
                            double det = JxW_g * alpha_beta * normal[ivar];
                            double Pressure = _data_eq[1].ub[(_FF_idx[NS_F]) * NDOF_FEM + sur_toply[j]];
//                             FeM(sur_toply[lbnode])  -= det * _phi_g[2][lbnode] * _phi_g[1][j] * (pressure[sur_toply[j]] /*- pressure_old[sur_toply[j]]*/);//pressure[sur_toply[j]];
                        }
                    }// END LOOP IF BC FLAG IS FOR INTEGRATION
                }// END GAUSSIAN INTEGRATION ========================================================================================
            }
        }
    }
    return;
}


#endif
#endif






