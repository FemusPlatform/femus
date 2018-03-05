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
    const double delta = _Wall_dist;
    double pressure;
    double vel_bound[1];
    double xyz_g[DIMENSION];
    const int elb_dof0 = (NDOF_K > 1) ? elb_ndof[_pres_order] * NDOF_K : elb_ndof[_pres_order];
    double det2 = _fe[2]->JacSur(elb_ngauss - 1, _xxb_qnds, _InvJac2); // jacobian
    double Ipenalty = det2;                          // Dirichlet bc flagdouble xyz_bg[DIMENSION];
    Ipenalty = 1.e0;

    const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;
    const int norm_face = (abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 10000) / 1000 - 1 ;

    // TANGENT VECTORS BASED ON BOUNDARY GEOMETRY
    double dir_tang1[DIMENSION], dir_tang2[DIMENSION];
    CalcTangDir(dir_tang1, dir_tang2, normal, 0);
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
        const int bc_var_npt = bc_var_check2 / 10000;           // recurrent  points in one element
        const int bc_var_check = bc_var_check2 % 10000;         // bc_var
        const int bc_var_normal = bc_var_check / 1000 - 1;      //
        const int bc_var = bc_var_check % 1000;
        const int bc_n_flag = bc_var / 10;
        const int bc_tg_flag = bc_var % 10;
        const int norm_node = (abs(_bc_bd[sur_toply[lbnode]]) % 10000) / 1000 - 1 ;

        if (_nNSdim == 3) {
            CorrectBCFlags(dir_tang1, dir_tang2, bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
        }

        if (bd_face == bd_node) {
//       for(int row_shift=0; row_shift< _nvars[2]; row_shift++) {  // LOOP OVER SPACE DIRECTIONS =======================================

            const int  indx_row = sur_toply[lbnode] + _dir * el_ndof[2]; // solution indx
            const int  node = sur_toply[lbnode];
            int a = 0;
            if (el_conn[node] == 38 && _dir == 0) {
                a = 1;
            }
            if (el_conn[node] == 38 && _dir == 1) {
                a = 1;
            }
            if (el_conn[node] == 38 && _dir == 2) {
                a = 1;
            }

            bool NormBC = (abs(_bc_el[indx_row]) == 3) ? true : false;
            bool DirichletBC = (_bc_el[indx_row] > 0) ? true : false;

            double u_old_node[3], u_oold_node[3];
            for (int jvar = 0; jvar < _nNSdim; jvar++) {
                u_old_node[jvar] =  _data_eq[2].ub[(_FF_idx[NS_F] + jvar) * NDOF_FEM + node];
                u_oold_node[jvar] =  u_oold[jvar * NDOF_FEM + node];
            }

            if (DirichletBC && norm_face == norm_node/*!_AlreadyWrittenDirBC[indx_row]*/) {
//           _AlreadyWrittenDirBC[indx_row] = true;
                int  bc_bc  = NormBC ? bc_n_flag : bc_tg_flag;
                const int  bc_rhs = ((bc_bc % 2) == 1) ? 1 : 0; // bc_rhs -> rhs bc
                bool OrDir = (bc_bc > 7) ? true : false;
                const double vart = bc_rhs * _data_eq[2].ub[indx_row + _FF_idx[NS_F] * el_ndof[2]];

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
                        for (int   jvar = _dir ; jvar < _dir + _nNSdim; jvar++) {
                            const int comp = jvar % _nNSdim;
                            alpha_beta = -normal[comp];
                            if (!NormBC) {
                                alpha_beta = (_bc_el[indx_row] == 1) ? dir_tang1[comp] : dir_tang2[comp];
                            }
                            const int ImpCon = (jvar > _dir) ? 0 : 1;
                            if (ImpCon == 1) {// implicit contribution
                                KeM(node, node) = alpha_beta;
                            }
                            if (ImpCon == 0) {// explicit contribution
                                FeM(node)       -= u_old_node[comp] * alpha_beta ;
                            }
                            FeM(node) += u_old_node[comp] * bc_rhs * alpha_beta;
                        }
                    }
                    if (a == 1) {
                        std::cout << KeM(node, node) << "   " << FeM(node) << std::endl;
                    }
                } else {
                    KeM(node, node) = alpha_beta;
                    FeM(node) = vart * alpha_beta;
                }// beta_0 u -> diagonal
            }
//       }// END DIRICHLET BLOCK ====================================================================================
        }
    }

    // ...........................................................................................
    //     INTEGRATION OVER BOUNDARY ELEMENT - STRESS CALCULATION
    // ...........................................................................................

    double dt = _dt;
    if (dt > 1 || _NS_parameter._SolveSteady == 1) {
        dt = 1. ;
    }
    if (dt < 1.e-10) {
        dt = 1.e-10;
    }
    double tang1[DIMENSION], tang2[DIMENSION];
    double VelOnGauss[DIMENSION], u_tau, y_plus, eff_stress, VelOnBound, veltg;

    const int method = (_nNSdim == 2) ? 0 : 1;
    CalcTangDir(tang1, tang2, normal, method);

    for (int  qp = 0; qp < elb_ngauss; qp++) { // START LOOP OVER GAUSS POINTS ===============================================
        double det   = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);   // local coord _phi_g and jac
        double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp]/**dt*/; // weight
        _fe[2]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[2]);   // global coord _phi_g
        _fe[1]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[1]);   // global coord _phi_g

        if (_AxiSym == 1) {
            interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg);
            JxW_g  *= _xyzg[0];   // axisymmetric  (index ->0)
        }

        interp_el_bd_sol(u_old, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], VelOnGauss);
//     {
//       // Calculation of effective stress through u_tau
//       double unorm=0;  double utot=0;
//       for(int  kdim=0; kdim<_nNSdim; kdim++) {
//         unorm +=VelOnGauss[kdim]*normal[kdim];
//         utot += VelOnGauss[kdim]*VelOnGauss[kdim];
//       }
//       const double WDist =_y_bcout;// _mgutils._TurbParameters->_BoundWallDist;
//       veltg  = sqrt(utot - unorm*unorm);
//       u_tau      = _mgutils._TurbParameters->CalcUtau(veltg + 1.e-10,WDist);
//       y_plus     = WDist*u_tau/_IRe;
//       eff_stress = (fabs(veltg) < 1.e-10)? _IRe/(WDist):u_tau*u_tau;
//       VelOnBound = sqrt(utot);
//     }

        for (int i = 0; i < elb_ndof[2]; i++) { // LOOP OVER TEST FUNCTIONS ===================================================

            const int bc_var_check2 = _bc_bd[sur_toply[i]]; // total  bc_var
            const int bc_var_npt = bc_var_check2 / 10000;                           // recurrent  points in one element
            const int sign = bc_var_check2 / abs(bc_var_check2);
            const int bc_var_check = abs(_bc_bd[sur_toply[i]]) % 10000;                         // bc_var
            const int bc_var_normal = bc_var_check / 1000 - 1;                            //
            const int bd_node = abs(_bc_bd[sur_toply[i]]) % 100;
            const int bc_varn_flag = bd_node / 10;
            const int bc_vartg_flag = bd_node % 10;

            double phii_g  = _phi_g[2][i];

            if (bd_face == bd_node) {

                const int pen_cond = (bd_node == 44) ? 0 : 1;
                double lambda1 = (1 - pen_cond) * _NS_parameter._Penalty_tg;
                double beta = lambda1 * _NS_parameter._Penalty_n;
                // Direction of maximum vector component of tangent 1 and 2 -> alignment of equations
                if (_nNSdim == 3) {
                    CorrectBCFlags(tang1, tang2, bc_var_normal, sur_toply[i], el_ndof[2], 1);
                }

//                     for ( int  row_shift=0; row_shift< _nNSdim; row_shift++ ) { // LOOP OVER ROWS: i + row_shift*el_ndof2
                int indx     = sur_toply[i] + _dir * el_ndof[2];
                int NComp = (_bc_el[indx] < 0) ? _nNSdim : 1;
                int bc_el = (_bc_el[indx] <= 0) ? 1 : 0;
                const double dtJxW_g = bc_el * JxW_g;

                bool NormBC = (_bc_el[indx] == -3) ? true : false;
                const int bc_flag = (NormBC) ? bd_node / 10 : bd_node % 10;

                if (bc_flag < 6 && bd_node != 11) { // only for neumann bc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // NORMAL VECTOR FOR THE PROJECTION OF VOLUME PART EQUATION
                    // _normal_pt is used to project the equation (volume part) using the same normal used in the boundary
                    if (qp == 0) for (int dir = 0; dir < _nNSdim; dir++) {
                            _normal_pt[ sur_toply[i]*_nNSdim + dir] += normal[dir];
                        }

                    for (int  ivarN = _dir; ivarN < _dir + NComp; ivarN++) { // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
                        const int ivar = ivarN % DIMENSION;
                        const int ImpContrib = (ivarN > _dir) ? 0 : 1;
                        // Stress component along tangential directions (-1,-2) ---------------------------------------------------
                        if (_bc_el[indx] == -1 || _bc_el[indx] == -2) {
                            if (bc_flag != 1) {
//                                    _                              _
//                                   /   =>   ->    -->             /    ->   -->
//                                  /   (T *  n) *  phi_t ds  =    /  a (u *  phi_t) ds
//                                -s                             -s
                                double alpha_beta = (_bc_el[indx] == -1) ? tang1[ivar] : tang2[ivar];
                                if (bc_flag == 5) {
                                    alpha_beta = tang1[ivar];
                                }
                                double tau = (_bc_el[indx] == -1) ? (_NS_parameter._Tg1_stress) : _NS_parameter._Tg2_stress ;    //utau*utau/vel_bound;
                                double stress = alpha_beta * (pen_cond * tau - lambda1) * dtJxW_g * _phi_g[2][i];
                                if (ImpContrib == 1)
                                    for (int j = 0; j < elb_ndof[2]; j++) {
                                        KeM(sur_toply[i], sur_toply[j]) += -1.*stress * _phi_g[2][j];
                                    }
                                if (ImpContrib == 0) {
                                    FeM(sur_toply[i]) += stress * VelOnGauss[ivar];
                                }
                            }
                        }// End loop if equation is for tangential direction-----------------------------------------------------
                        else {
                            // Stress component along normal direction(-3,0) -----------------------------------------------------
                            double alpha_beta = 1.;
                            if (_bc_el[indx] == -3) {
                                alpha_beta = normal[_dir];
                                // normal with penalty (u dot n) n dot phi
//                                    _
//                                   /   ->   ->  ->   -->
//                                  /   (u *  n)  n *  phi_n ds = 0
//                                -s
                                if (bc_varn_flag == 4) {
                                    double norm_stress = dtJxW_g * alpha_beta * beta * _phi_g[2][i] * normal[_dir];
                                    if (ImpContrib == 1)
                                        for (int j = 0; j < elb_ndof[2]; j++) {
                                            KeM(sur_toply[i], sur_toply[j]) -= norm_stress * normal[ivar] * _phi_g[2][j];
                                        }
                                    if (ImpContrib == 0) {
                                        FeM(sur_toply[i]) += norm_stress * normal[ivar] * VelOnGauss[ivar];
                                    }
                                }

                            } //  ---------------------------------------------------------------------------------------------------
                            // surface integral p phi dot n   --> pressure_inlet - pressure_outlet - outflow_p ------------------------------------------------
//                                         if ( bc_varn_flag > 1 && bc_varn_flag < 4 ) {
//                                              for ( int j=0; j< elb_dof0; j++ ) {
//                                                   double det = dtJxW_g*alpha_beta*normal[ivar];
//                                                   KeM ( indx, sur_toply[j]+ _nNSdim*el_ndof[2] ) += det*_phi_g[2][i]*_phi_g[1][j];
// //                       FeM(indx)  -= det*_phi_g[2][i]*_phi_g[1][j]*_data_eq[_pres_order].ub[(_FF_idx[NS_F]) *NDOF_FEM+_pres_order*sur_toply[j]];
//                                              }
//                                         } //  ---------------------------------------------------------------------------------------------------
                        }// contribution along normal direction

                    }// END LOOP IF BC FLAG IS FOR INTEGRATION
                }// END LOOP OVER VARIABLES FOR LINEAR COMBINATION OF EQUATIONS ++++++++++++++++++++++++++++++++++++++++++++++
//                     }// END LOOP OVER MATRIX ROWS RELATIVE TO TEST FUNCTION I
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
    double u_old[],
    int el_conn[]
)
{
    double dir_tang1[DIMENSION], dir_tang2[DIMENSION];
    CalcTangDir(dir_tang1, dir_tang2, normal, 0);
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
        if (_nNSdim == 3) {
            CorrectBCFlags(dir_tang1, dir_tang2, bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
        }
        const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;
        if (bd_face == bd_node && dir_maxnormal == bc_var_normal) {
            const int  indx_row = sur_toply[lbnode] + _dir * el_ndof[2]; // solution indx
            const int  node = sur_toply[lbnode];
            if (_bc_el[indx_row] > 0) {
//               int a =1;
                KeM(node, node) = 1.;
                FeM(node)       = _data_eq[2].ub[(_FF_idx[NS_F] + _dir) * NDOF_FEM + node];
            } else if (_bc_el[indx_row] < 0) {
                for (int dir = 0; dir < _nNSdim; dir++) {
                    _normal_pt[ node*_nNSdim + dir] = normal[dir];
                }
            }
        }
    }

    return;
}


// =======================================================================================================
void MGSolNS::SetBCFlags(
    double normal[], // normal unit
    int sur_toply[], // surface connectivity (respect to volume)
    int el_ndof,     // number of dofs
    int el_conn[]    // volume connectivity
)  // -----------------------------------------------------------------------------------------------------
{
    // bc conditions
    // Dirichlet >5  (surface point setting ( no volume integration))
    // 6-> set zero component 7-> set component       ------------------------> _bc_el[] = 3,0,-3; _bc_el[] =1,2,-1,-2;
    // 8-> zero along dir (normal or tg) 9-> value along dir  ----------------> _bc_el[] = 3,0,-3; _bc_el[] =1,2,-1,-2;
    //  Neumann (volume integration (+surface integration))   ---------------->
    // 1-> no surface integration                             ---------------->
    // 2-> pressure integration (normal only)  3-> pressure setting    ------->
    // 4-> tau=a.velocity                      4-> stau=a.(velocity-velocity_0)------>
    const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100; // element normal-tg bc flag
    const int norm_face = (abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 10000) / 1000 - 1 ; // element max-normal-dir bc flag

    for (int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++) { // loop on surface nodes
        const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;   // point normal-tg  bc flag
        const int bc_var_check = abs(_bc_bd[sur_toply[lbnode]]) % 10000;  // bc_var
        int bc_var_normal = bc_var_check / 1000 - 1; // point max-normal-dir  bc flag
        const int bc_n_flag = bd_node / 10;         // normal bc flag
        const int bc_tg_flag = bd_node % 10;         // tg bc flag

        bool Norm_Dirichlet_Bound_Cond = (bc_n_flag > 5) ? true : false;
        bool Tang_Dirichlet_Bound_Cond = (bc_tg_flag > 5) ? true : false;
        int NormDir = bc_var_normal;

        // BC are imposed only if node bc equal to face bc
        if (bd_face == bd_node) {
            int NormRow = sur_toply[lbnode] + NormDir * el_ndof;

            // Setting bc flags for normal direction --------------------------
            if (_bc_el[NormRow] == _BdFlagId) {   // only _BdFlagId=-8 (first time)
                if (Norm_Dirichlet_Bound_Cond) {
                    _bc_el[NormRow] = 3;
                } else { // Neumann bc (surface integration)
                    if (bc_n_flag == 1 || bd_node == 31) {
                        _bc_el[NormRow] = 0;
                    }   else {
                        _bc_el[NormRow] = -3;
                    }
                }
            }//  norm  ----------------------------------------------------------
            // Setting bc flags for tangential directions -------------------
            int dir_tg0 = 0, NeuTg = 0;
            for (int kdir = NormDir + 1; kdir < NormDir + _nNSdim; kdir ++) {
                int idir = kdir % _nNSdim;
                int tg_dirRow = sur_toply[lbnode] + idir * el_ndof;
                if (_bc_el[tg_dirRow] == _BdFlagId) { // only _BdFlagId=-8 (first time)
                    if (Tang_Dirichlet_Bound_Cond) {
                        dir_tg0 ++;     // Dirichlet bc
                        _bc_el[tg_dirRow] = dir_tg0;
                    } else { // Neumann bc (surface integration)
                        if (bd_node == 11 || bd_node == 31) {
                            _bc_el[tg_dirRow] = 0;    // for outflow bc we don't project the equation along normal-tangential directions
                        } else {
                            NeuTg++;
                            _bc_el[tg_dirRow] = -NeuTg;
                        }
                    }
                }
            }//  Tan ----------------------------------------------------------
        }// if (bd_face==bd_node)
    }

    return;
}

// ===========================================================================
void MGSolNS::CorrectBCFlags(
    double Tang1[],  // tan in dir 1 (flow align)
    double Tang2[],  // tan in dir 2 (transverse flow)
    int MaxNormal,   // dir max normal
    int ElementNode, // node id
    int ElDof,       // elemet dofs
    int BcType       // 0 dirichlet 1 neumann
)
{
    int tg1 = (MaxNormal + 1) % DIMENSION;
    int tg2 = (MaxNormal + 2) % DIMENSION;

    int tang1_maxdir = (fabs(Tang1[tg1]) > fabs(Tang1[tg2])) ? tg1 : tg2;
    int tang2_maxdir = (fabs(Tang2[tg1]) > fabs(Tang2[tg2])) ? tg1 : tg2;

    if (tang1_maxdir == tang2_maxdir) {
        int last_dir = (tang1_maxdir == tg1) ? tg2 : tg1;
        if (fabs(Tang1[tang1_maxdir]) > fabs(Tang2[tang2_maxdir])) {
            tang2_maxdir = last_dir;
        } else {
            tang1_maxdir = last_dir;
        }
    }

    if (BcType == 1) { // neumann
        if (_bc_el[ElementNode + (tg1) *ElDof]*_bc_el[ElementNode +  tg2 * ElDof] == 2
                && _bc_el[ElementNode +  tg1 * ElDof] + _bc_el[ElementNode +  tg2 * ElDof] == -3) {
            _bc_el[ElementNode + (tang1_maxdir) *ElDof] = -1;
            _bc_el[ElementNode + (tang2_maxdir) *ElDof] = -2;
        }
    } else if (BcType == 0) { // dirichlet
        if (_bc_el[ElementNode + (tg1) *ElDof]*_bc_el[ElementNode +  tg2 * ElDof] == 2
                && _bc_el[ElementNode +  tg1 * ElDof] + _bc_el[ElementNode +  tg2 * ElDof] == 3) {
            _bc_el[ElementNode + (tang1_maxdir) *ElDof] = 1;
            _bc_el[ElementNode + (tang2_maxdir) *ElDof] = 2;
        }
    }

    return;
}



#endif




#endif





