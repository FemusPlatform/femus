// ===============================================================
// --------------   NAVIER-STOKES system [FS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef FSI_EQUATIONS
// ==============================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#if FSI_EQUATIONS == 1
// local Femus class include -----------------------------------
// class files --------------------------------------------------
#include <iomanip>          // std::setprecision
#include "MGFE.h"           // Mesh class
#include "MGSclass_conf.h"  // Navier-Stokes class conf file
#include "MGSolverFSI.h"    // Navier-Stokes class header file

void MGSolFSI::set_bc_matrix(
    DenseMatrixM& KeM, DenseVectorM& FeM,  ///< local Matrix and rhs
    int dir_maxnormal,                     ///<  normal dir
    int sur_toply[],                       ///< boundary topology map
    int el_ndof[],                         ///< number of volume dofs
    int elb_ndof[],                        ///< number of boundary dofs
    int elb_ngauss,                        ///<  number of surface gaussian points
    double normal[],                       ///< normal
    double u_old[], int el_conn[]) {
  double xyz_g[DIMENSION];
  const int elb_dof0 = (NDOF_K > 1) ? elb_ndof[_pres_order] * NDOF_K : elb_ndof[_pres_order];
  double det2 = _fe[2]->JacSur(elb_ngauss - 1, _xxb_qnds, _InvJac2);  // jacobian
  double Ipenalty = det2;  // Dirichlet bc flagdouble xyz_bg[DIMENSION];
  Ipenalty = 1.e0;

  for (int ivar = 0; ivar < DIMENSION; ivar++) {
    xyz_g[ivar] = 0;
    for (int i = 0; i < NDOF_FEMB; i++) { xyz_g[ivar] += _xxb_qnds[i + ivar * NDOF_FEMB] / NDOF_FEMB; }
  }

  const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;

  // TANGENT VECTORS BASED ON BOUNDARY GEOMETRY
  double dir_tang1[DIMENSION], dir_tang2[DIMENSION];
  CalcTangDir(dir_tang1, dir_tang2, normal, 0);
  SetBCFlags(normal, sur_toply, el_ndof[2], el_conn);

  // ...........................................................................................
  //     DIRICHLET BOUNDARY CONDITION
  // ...........................................................................................

  for (int lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {
    const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;
    const int bc_var_check2 = abs(_bc_bd[sur_toply[lbnode]]);  // total  bc_var
    const int bc_var_check = bc_var_check2 % 10000;            // bc_var
    const int bc_var_normal = bc_var_check / 1000 - 1;         //
    const int bc_var = bc_var_check % 1000;
    const int bc_n_flag = bc_var / 10;
    const int bc_tg_flag = bc_var % 10;

    if (_nNSdim == 3) {
      CorrectBCFlags(dir_tang1, dir_tang2, bc_var_normal, sur_toply[lbnode], el_ndof[2], 0);
    }

    if (bd_face == bd_node) {
      for (int row_shift = 0; row_shift < _nvars[2];
           row_shift++) {  // LOOP OVER SPACE DIRECTIONS =======================================
        const int indx_row = sur_toply[lbnode] + row_shift * el_ndof[2];  // solution indx

        bool NormBC = (abs(_bc_el[indx_row]) == 3) ? true : false;
        bool DirichletBC = (_bc_el[indx_row] > 0) ? true : false;

        if (DirichletBC /*!_AlreadyWrittenDirBC[indx_row]*/) {
          //           _AlreadyWrittenDirBC[indx_row] = true;
          int bc_bc = NormBC ? bc_n_flag : bc_tg_flag;
          const int bc_rhs = ((bc_bc % 2) == 1) ? 1 : 0;  // bc_rhs -> rhs bc
          bool OrDir = (bc_bc > 7) ? true : false;
          const double vart = bc_rhs * _data_eq[2].ub[indx_row + _FF_idx[FS_F] * el_ndof[2]];

          double alpha_beta = 1.;
          if (OrDir) {
            if (bd_node == 88 || bd_node == 66) {
              KeM(indx_row, indx_row) = 1.;
              FeM(indx_row) = 0.;
              //                             if (xyz_g[0]>0.199&&xyz_g[0]<0.301 && xyz_g[1]>0.4)
              //                                     FeM(indx_row) += -20*_dx_old [0][0]*row_shift;
              //                             if (xyz_g[0]>0.199&&xyz_g[0]<0.301 && xyz_g[1]<0.01)
              //                                     FeM(indx_row) +=  20*_dx_old [0][0]*row_shift;
            } else {
              for (int jvar = 0; jvar < _nNSdim; jvar++) {  //  -beta_n*n*( u.n)
                const int indj = sur_toply[lbnode] + jvar * el_ndof[2];
                alpha_beta = -normal[jvar];
                if (!NormBC) { alpha_beta = (_bc_el[indx_row] == 1) ? dir_tang1[jvar] : dir_tang2[jvar]; }
                KeM(indx_row, indj) += alpha_beta;
                FeM(indx_row) += _data_eq[2].ub[(_FF_idx[FS_F]) * NDOF_FEM + indj] * bc_rhs * alpha_beta;
              }
            }
          } else {
            KeM(indx_row, indx_row) = alpha_beta;
            FeM(indx_row) = vart * alpha_beta;
          }  // beta_0 u -> diagonal
        }
      }  // END DIRICHLET BLOCK
         // ====================================================================================
    }
  }

  // ...........................................................................................
  //     PRESSURE BOUNDARY CONDITION
  // ...........................................................................................

  const int bc_side_check2 = _bc_bd[sur_toply[elb_ndof[2] - 1]];  // total  bc_var
  const int bc_side_n_flag = ((bc_side_check2 % 10000) % 1000) / 10;

  if (bc_side_n_flag < 4 && bc_side_n_flag > 1) {
    const int rhs_p = bc_side_n_flag % 2;
    for (int i = 0; i < elb_ndof[1]; i++) {
      const int indx = _pres_order * sur_toply[i] + _nNSdim * NDOF_FEM;  // volume dof index
      _bc_el[indx] = 0;                                                  // SET TO 1 -> CONTROL IF RIGHT
      KeM(indx, indx) += Ipenalty;
      FeM(indx) += rhs_p * Ipenalty * _data_eq[_pres_order].ub[_pres_order * sur_toply[i]];
    }
  }

  // ...........................................................................................
  //     INTEGRATION OVER BOUNDARY ELEMENT - STRESS CALCULATION
  // ...........................................................................................

  double dt = _dt;
  if (dt > 1) { dt = 1.; }
  if (dt < 1.e-10) { dt = 1.e-10; }
  double tang1[DIMENSION], tang2[DIMENSION];
  CalcTangDir(tang1, tang2, normal, 0);

  if (bd_face == 44) {
    std::cout << "\n ------------------ \n";
    double mod1 = 0., mod2 = 0.;
    for (int dir = 0; dir < 3; dir++) {
      std::cout << dir << "  " << tang1[dir] << "    " << tang2[dir] << std::endl;
      mod1 += tang1[dir] * tang1[dir];
      mod2 += tang2[dir] * tang2[dir];
    }
    std::cout << mod1 << "  " << mod2 << std::endl;
  }

  for (int qp = 0; qp < elb_ngauss;
       qp++) {  // START LOOP OVER GAUSS POINTS ===============================================
    double det = _fe[2]->JacSur(qp, _xxb_qnds, _InvJac2);         // local coord _phi_g and jac
    double JxW_g = det * _fe[2]->_weight1[_nNSdim - 2][qp] * dt;  // weight
    _fe[2]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[2]);             // global coord _phi_g
    _fe[1]->get_phi_gl_g(_nNSdim - 1, qp, _phi_g[1]);             // global coord _phi_g

    if (_AxiSym == 1) {
      interp_el_bd_sol(_xx_qnds, sur_toply, el_ndof[2], 0, _nNSdim, _phi_g[2], elb_ndof[2], _xyzg);
      JxW_g *= _xyzg[0];  // axisymmetric  (index ->0)
    }

    //         interp_el_bd_sol(u_old,sur_toply,el_ndof[2],0,_nNSdim,_phi_g[2],elb_ndof[2],VelOnGauss);

    for (int i = 0; i < elb_ndof[2];
         i++) {  // LOOP OVER TEST FUNCTIONS ===================================================

      const int bc_var_check = abs(_bc_bd[sur_toply[i]]) % 10000;  // bc_var
      const int bc_var_normal = bc_var_check / 1000 - 1;           //
      const int bd_node = abs(_bc_bd[sur_toply[i]]) % 100;
      const int bc_varn_flag = bd_node / 10;
      const int bc_vartg_flag = bd_node % 10;

      if (bd_face == bd_node) {
        const int pen_cond = (bd_node == 44) ? 0 : 1;
        double lambda1 = (1 - pen_cond) * _FSI_parameter._Penalty_tg;
        double beta = lambda1 * _FSI_parameter._Penalty_n;
        // Direction of maximum vector component of tangent 1 and 2 -> alignment of equations
        if (_nNSdim == 3) { CorrectBCFlags(tang1, tang2, bc_var_normal, sur_toply[i], el_ndof[2], 1); }

        for (int row_shift = 0; row_shift < _nvars[2];
             row_shift++) {  // LOOP OVER ROWS: i + row_shift*el_ndof2
          int indx = sur_toply[i] + row_shift * el_ndof[2];
          int NComp = (_bc_el[indx] < 0) ? _nNSdim : 1;
          int bc_el = (_bc_el[indx] <= 0) ? 1 : 0;
          const double dtJxW_g = bc_el * JxW_g;

          bool NormBC = (_bc_el[indx] == -3) ? true : false;
          const int bc_flag = (NormBC) ? bd_node / 10 : bd_node % 10;

          if (bc_flag < 6 &&
              bd_node !=
                  11) {  // only for neumann bc
                         // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // NORMAL VECTOR FOR THE PROJECTION OF VOLUME PART EQUATION
            // _normal_pt is used to project the equation (volume part) using the same normal used in the
            // boundary
            if (qp == 0)
              for (int dir = 0; dir < _nNSdim; dir++) {
                _normal_pt[sur_toply[i] * _nNSdim + dir] += normal[dir];
              }

            for (int ivarN = row_shift; ivarN < row_shift + NComp;
                 ivarN++) {  // LOOP OVER VELOCITY COMPONENTS FOR EQUATION WRITTEN IN ROW indx
              const int ivar = ivarN % DIMENSION;
              // Stress component along tangential directions (-1,-2)
              // ---------------------------------------------------
              if (_bc_el[indx] == -1 || _bc_el[indx] == -2) {
                if (bc_vartg_flag != 1) {
                  double alpha_beta = (_bc_el[indx] == -1) ? tang1[ivar] : tang2[ivar];
                  if (bc_vartg_flag == 5) { alpha_beta = tang1[ivar]; }
                  double tau = (_bc_el[indx] == -1) ? _FSI_parameter._Tg1_stress
                                                    : _FSI_parameter._Tg2_stress;  // utau*utau/vel_bound;
                  double stress = alpha_beta * (pen_cond * tau - lambda1) * dtJxW_g * _phi_g[2][i];
                  for (int j = 0; j < elb_ndof[2]; j++) {
                    KeM(indx, sur_toply[j] + ivar * el_ndof[2]) += -1. * stress * _phi_g[2][j];
                  }  // END LOOP OVER MATRIX COMPONENTS - IMPLICIT STRESS
                }
              }  // End loop if equation is for tangential
                 // direction-----------------------------------------------------
              else {
                // Stress component along normal direction(-3,0)
                // -----------------------------------------------------
                double alpha_beta = 1.;
                if (_bc_el[indx] == -3) {
                  alpha_beta = normal[ivar];
                  // normal with penalty (u dot n) n dot phi
                  if (bc_varn_flag == 4)
                    for (int j = 0; j < elb_ndof[2]; j++)
                      for (int jvar = 0; jvar < _nNSdim; jvar++) {
                        KeM(indx, sur_toply[j] + jvar * el_ndof[2]) -= alpha_beta * normal[ivar] *
                                                                       normal[jvar] * beta * dtJxW_g *
                                                                       _phi_g[2][j] * _phi_g[2][i];
                      }
                }  //  ---------------------------------------------------------------------------------------------------
                // surface integral p phi dot n   --> pressure_inlet - pressure_outlet - outflow_p
                // ------------------------------------------------
                if (bc_varn_flag > 1 && bc_varn_flag < 4) {
                  for (int j = 0; j < elb_dof0; j++) {
                    double det = dtJxW_g * alpha_beta * normal[ivar];
                    KeM(indx, sur_toply[j] + _nNSdim * el_ndof[2]) += det * _phi_g[2][i] * _phi_g[1][j];
                    //                       FeM(indx)  -=
                    //                       det*_phi_g[2][i]*_phi_g[1][j]*_data_eq[_pres_order].ub[(_FF_idx[FS_F])
                    //                       *NDOF_FEM+_pres_order*sur_toply[j]];
                  }
                }  //  ---------------------------------------------------------------------------------------------------
              }  // contribution along normal direction

            }  // END LOOP IF BC FLAG IS FOR INTEGRATION
          }    // END LOOP OVER VARIABLES FOR LINEAR COMBINATION OF EQUATIONS
               // ++++++++++++++++++++++++++++++++++++++++++++++
        }      // END LOOP OVER MATRIX ROWS RELATIVE TO TEST FUNCTION I
      }        // END IF BC_NODE == BC_FACE
    }          // END LOOP OVER BOUNDARY TEST FUNCTIONS
               // =========================================================================
  }            // END GAUSSIAN INTEGRATION
               // ========================================================================================

  return;
}

// =======================================================================================================
void MGSolFSI::SetBCFlags(
    double normal[],  // normal unit
    int sur_toply[],  // surface connectivity (respect to volume)
    int el_ndof,      // number of dofs
    int el_conn[]     // volume connectivity
) {  // -----------------------------------------------------------------------------------------------------
  // bc conditions
  // Dirichlet >5  (surface point setting ( no volume integration))
  // 6-> set zero component 7-> set component       ------------------------> _bc_el[] = 3,0,-3; _bc_el[]
  // =1,2,-1,-2; 8-> zero along dir (normal or tg) 9-> value along dir  ----------------> _bc_el[] = 3,0,-3;
  // _bc_el[] =1,2,-1,-2;
  //  Neumann (volume integration (+surface integration))   ---------------->
  // 1-> no surface integration                             ---------------->
  // 2-> pressure integration (normal only)  3-> pressure setting    ------->
  // 4-> tau=a.velocity                      4-> stau=a.(velocity-velocity_0)------>
  const int bd_face = abs(_bc_bd[sur_toply[NDOF_FEMB - 1]]) % 100;  // element normal-tg bc flag

  for (int lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {                // loop on surface nodes
    const int bd_node = abs(_bc_bd[sur_toply[lbnode]]) % 100;         // point normal-tg  bc flag
    const int bc_var_check = abs(_bc_bd[sur_toply[lbnode]]) % 10000;  // bc_var
    int bc_var_normal = bc_var_check / 1000 - 1;                      // point max-normal-dir  bc flag
    const int bc_n_flag = bd_node / 10;                               // normal bc flag
    const int bc_tg_flag = bd_node % 10;                              // tg bc flag

    bool Norm_Dirichlet_Bound_Cond = (bc_n_flag > 5) ? true : false;
    bool Tang_Dirichlet_Bound_Cond = (bc_tg_flag > 5) ? true : false;
    int NormDir = bc_var_normal;

    // BC are imposed only if node bc equal to face bc
    if (bd_face == bd_node) {
      int NormRow = sur_toply[lbnode] + NormDir * el_ndof;

      // Setting bc flags for normal direction --------------------------
      if (_bc_el[NormRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
        if (Norm_Dirichlet_Bound_Cond) {
          _bc_el[NormRow] = 3;
        } else {  // Neumann bc (surface integration)
          if (bc_n_flag == 1 || bd_node == 31) {
            _bc_el[NormRow] = 0;
          } else {
            _bc_el[NormRow] = -3;
          }
        }
      }  //  norm  ----------------------------------------------------------
      // Setting bc flags for tangential directions -------------------
      int dir_tg0 = 0, NeuTg = 0;
      for (int kdir = NormDir + 1; kdir < NormDir + _nvars[2]; kdir++) {
        int idir = kdir % _nNSdim;
        int tg_dirRow = sur_toply[lbnode] + idir * el_ndof;
        if (_bc_el[tg_dirRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
          if (Tang_Dirichlet_Bound_Cond) {
            dir_tg0++;  // Dirichlet bc
            _bc_el[tg_dirRow] = dir_tg0;
          } else {  // Neumann bc (surface integration)
            if (bd_node == 11 || bd_node == 31) {
              _bc_el[tg_dirRow] =
                  0;  // for outflow bc we don't project the equation along normal-tangential directions
            } else {
              NeuTg++;
              _bc_el[tg_dirRow] = -NeuTg;
            }
          }
        }
      }  //  Tan ----------------------------------------------------------
    }    // if (bd_face==bd_node)
    else if (bd_node == 88 || bd_node == 66) {
      int NormRow = sur_toply[lbnode] + NormDir * el_ndof;

      // Setting bc flags for normal direction --------------------------
      if (_bc_el[NormRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
        _bc_el[NormRow] = 3;
      }  //  norm  ----------------------------------------------------------
      // Setting bc flags for tangential directions -------------------
      int dir_tg0 = 0;
      for (int kdir = NormDir + 1; kdir < NormDir + _nvars[2]; kdir++) {
        int idir = kdir % _nNSdim;
        int tg_dirRow = sur_toply[lbnode] + idir * el_ndof;
        if (_bc_el[tg_dirRow] == _BdFlagId) {  // only _BdFlagId=-8 (first time)
          dir_tg0++;                           // Dirichlet bc
          _bc_el[tg_dirRow] = dir_tg0;
        }
      }  //  Tan ----------------------------------------------------------
    }
  }

  return;
}

// ===========================================================================
void MGSolFSI::CorrectBCFlags(
    double Tang1[],   // tan in dir 1 (flow align)
    double Tang2[],   // tan in dir 2 (transverse flow)
    int MaxNormal,    // dir max normal
    int ElementNode,  // node id
    int ElDof,        // elemet dofs
    int BcType        // 0 dirichlet 1 neumann
) {
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

  if (BcType == 1) {  // neumann
    if (_bc_el[ElementNode + (tg1)*ElDof] * _bc_el[ElementNode + tg2 * ElDof] == 2 &&
        _bc_el[ElementNode + tg1 * ElDof] + _bc_el[ElementNode + tg2 * ElDof] == -3) {
      _bc_el[ElementNode + (tang1_maxdir)*ElDof] = -1;
      _bc_el[ElementNode + (tang2_maxdir)*ElDof] = -2;
    }
  } else if (BcType == 0) {  // dirichlet
    if (_bc_el[ElementNode + (tg1)*ElDof] * _bc_el[ElementNode + tg2 * ElDof] == 2 &&
        _bc_el[ElementNode + tg1 * ElDof] + _bc_el[ElementNode + tg2 * ElDof] == 3) {
      _bc_el[ElementNode + (tang1_maxdir)*ElDof] = 1;
      _bc_el[ElementNode + (tang2_maxdir)*ElDof] = 2;
    }
  }

  return;
}

#endif

#endif
