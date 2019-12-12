// #include "Equations_tab.h"
// #include "Equations_conf.h"
// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ----------
// #include "Printinfo_conf.h"

// class files ----------------
#include "MGSolverDA.h"
// #include "MGSclass_conf.h"

// local includes --------------
#include "EquationSystemsExtendedM.h"
// #include "MGFE_conf.h"
#include "MGGeomEl.h"
// #include "MGSystem.h"
#include "MeshExtended.h"
// #include "dense_set.h"

// #include "MGEquationsSystem.h"
#include "MGFE.h"
// #include "MGFEMap.h"
#include "MGGraph.h"
#include "MGUtils.h"

// #ifdef HAVE_MED
// #include "InterfaceFunctionM.h"
// #include "MEDCouplingFieldDouble.hxx"
// #include "MEDCouplingUMesh.hxx"
// #include "MEDLoader.hxx"
// #endif

// algebric includes -----------
// #include "dense_matrixM.h"
// #include "dense_vectorM.h"
// #include "numeric_vectorM.h"
// #include "sparse_MmatrixM.h"
#include "sparse_matrixM.h"

// ==========================================================
//               MGSolDA functions
// ==========================================================

///  ====================================================
/// This function assembles the matrix and the rhs:
///  ====================================================
void MGSolDA::GenMatRhs(
    const double /*time*/,  // time  <-
    const int Level,        // Level <-
    const int mode          // mode  <- (1=rhs+matrix) (0=only matrix)
) {                         // ===============================================

  /// Set up
  // geometry and bc---------------------------------------------------------------------------------
  const int ndim = DIMENSION;                          // dimension
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];    // element sides
  int el_conn[NDOF_FEM];                               // elb_conn[NDOF_FEMB];   // element connectivity
  int el_neigh[NDOF_FEM];                              // bd element connectivity
  int elb_conn[NDOF_FEMB];                             // elb_conn[NDOF_FEMB];   // element connectivity
  int sur_toply[NDOF_FEMB];                            // boundary topology
  double xx_qnds[DIMENSION * NDOF_FEM];                // element node coords
  double xxb_qnds[DIMENSION * NDOF_FEMB];              // boundary of the element node coords
  int _bc_vol[NDOF_FEM * 10];                          // element  b.cond flags (Neu or Dir)
  int _bc_bd[NDOF_FEM * 10];                           // element  b.cond flags (different possibilities)

  // gauss integration  -----------------------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim - 1];          // quadratic element gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION - 2];    // quadratic bd elem gauss points
  double det[3], JxW_g[3], InvJac[3][DIMENSION * DIMENSION];  // determinant of Jacobean and det*gauss_weigth
  double dphijdx_g[3][DIMENSION];                             // global derivatives at gauss point
  double dphiidx_g[3][DIMENSION];                             // global derivatives at gauss point
  //   double _dt =0.1;

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int elb_ndof[3];
  elb_ndof[0] = NDOF_BK;
  elb_ndof[1] = NDOF_PB;
  elb_ndof[2] = NDOF_FEMB;  // number of element boundary dofs
  int el_mat_nrows = 0;     // number of matrix rows (dofs)

  for (int ideg = 0; ideg < 3; ideg++) { el_mat_nrows += _nvars[ideg] * _el_dof[ideg]; }

  int el_mat_ncols = el_mat_nrows;                // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);  // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  int DA_f = _data_eq[2].indx_ub[_data_eq[2].tab_eqs[DA_F]];  // Quadratic
  double _ub_g[3][14];                                        // values of external fields

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix+rhs) ---------------------------
  A[Level]->zero();

  if (mode == 1) {
    b[Level]->zero();  // global matrix A and rhs b
  }

  DenseMatrixM KeM;
  DenseVectorM FeM;  // local  matrix KeM and rhs FeM
  KeM.resize(el_mat_nrows, el_mat_ncols);
  FeM.resize(el_mat_nrows);  // resize local matrix and rhs

  // number of total elements for level
  int ndof_lev = 0;

  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
    ndof_lev += delta;
  }

  // test pass
  // #ifdef HAVE_COUPLING
  //     MGSystemExtended * ext_es= static_cast<MGSystemExtended *>(&_mgphys);
  //      double test_pass= ext_es->_nComp;
  // #endif

  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];      // stop element

  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    KeM.zero();
    FeM.zero();  // set to zero matrix and rhs

    // geometry and element  fields ------------------------------------

    _mgmesh.get_el_nod_conn(
        0, Level, iel, el_conn, xx_qnds);  // Element Connectivity (el_conn) and coordinates (xx_qnds)
    _mgmesh.get_el_neighbor(el_sides, 0, Level, iel, el_neigh);  // Neighbors of the element

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd);

    // grid and gaussian points
    for (int idim = 0; idim < DIMENSION; idim++) {
      for (int d = 0; d < NDOF_FEM; d++) {
        _data_eq[2].ub[idim * NDOF_FEM + d] = xx_qnds[idim * NDOF_FEM + d];  // element nodes xxg (DIM)
      }

      // element grid distance
    }

    // element field values
    //         for (int deg=0; deg<3; deg++) {
    for (int eq = 0; eq < _data_eq[2].n_eqs; eq++) {
      _data_eq[2].mg_eqs[eq]->get_el_sol(
          0, 0, _data_eq[2].indx_ub[eq + 1] - _data_eq[2].indx_ub[eq], _el_dof[2], el_conn, offset,
          _data_eq[2].indx_ub[eq], _data_eq[2].ub);
    }

    //         }
    // linear field
    _data_eq[1].mg_eqs[0]->get_el_sol(
        0, _nvars[2], _data_eq[1].indx_ub[0], _el_dof[1], el_conn, offset, 0, _data_eq[1].ub);
    //  external node quantities -------------------------------------

    _data_eq[0].mg_eqs[0]->get_el_sol_piece(
        _nvars[2] + _nvars[1], _data_eq[0].indx_ub[0], _el_dof[0], iel + ndof_lev, offset, 0, _data_eq[0].ub);

    // ======================================================================
    // Volume =============================================================
    // ======================================================================
    // ---------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // --------------------------------------------
    for (int qp = 0; qp < el_ngauss; qp++) {
      // shape functions at gaussian points -----------------------------------
      for (int ideg = 1; ideg < 3; ideg++) {                          // linear-quadratic  [1]=linear [2]=quad
        det[ideg] = _fe[ideg]->Jac(qp, xx_qnds, InvJac[ideg]);        // Jacobian
        JxW_g[ideg] = det[ideg] * _fe[ideg]->_weight1[ndim - 1][qp];  // weight
        _fe[ideg]->get_phi_gl_g(ndim, qp, _phi_g[ideg]);              // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim, qp, InvJac[ideg], _dphi_g[ideg]);  // global coord deriv
      }

      JxW_g[0] = JxW_g[2];
      _fe[0]->get_phi_gl_g(ndim, qp, _phi_g[0]);  // shape function piecewise

      //  fields -----------------------------------------------------------
      interp_el_sol(
          _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2], _el_dof[2],
          _ub_g[2]);                                                                              // quadratic
      interp_el_sol(_data_eq[1].ub, 0, _data_eq[1].indx_ub[0], _phi_g[1], _el_dof[1], _ub_g[1]);  // linear
      interp_el_sol(_data_eq[0].ub, 0, _data_eq[0].indx_ub[0], _phi_g[0], _el_dof[0], _ub_g[0]);  // constant

      /// d) Local (element) assemblying DA equation
      // *********************** *******************************

      for (int i = 0; i < _el_dof[2]; i++) {  //  --- QUADRATIC ---
        // set up row i
        const double phii_g = _phi_g[2][i];

        for (int idim = 0; idim < ndim; idim++) { dphiidx_g[2][idim] = _dphi_g[2][i + idim * _el_dof[2]]; }

        //         int gl_node=el_conn[i];

        for (int ivar = 0; ivar < _nvars[2]; ivar++) {
          double dtxJxW_g = JxW_g[2] * _bc_vol[i + ivar * _el_dof[2]];
          int index = i + ivar * _el_dof[2];

          // Rhs Assemblying  ----------------------------------------------------------
          if (mode == 1) {
            // rhs

            FeM(index) += dtxJxW_g * (_ub_g[2][DA_f + ivar] * phii_g / _dt  // time
                                      + 2. * phii_g                         // heat source
                                     );
          }

          // Matrix Assemblying ---------------------------
          for (int j = 0; j < _el_dof[2]; j++) {
            double phij_g = _phi_g[2][j];
            double Lap = 0.;

            for (int idim = 0; idim < ndim; idim++) {
              dphijdx_g[2][idim] = _dphi_g[2][j + idim * _el_dof[2]];
              Lap += dphijdx_g[2][idim] * dphiidx_g[2][idim];  // Laplacian
            }

            // energy-equation
            KeM(index, j + ivar * _el_dof[2]) += dtxJxW_g * (phii_g * phij_g / _dt  // time term
                                                             + Lap                  // diff
                                                            );
          }
        }  // quadratic variable cycle
      }    // ----------------------------------------END QUADRATIC

      for (int i = 0; i < _el_dof[1]; i++) {  //    --- LINEAR ---
        // set up row i
        const double phii_g = _phi_g[1][i];

        for (int idim = 0; idim < ndim; idim++) { dphiidx_g[1][idim] = _dphi_g[1][i + idim * _el_dof[1]]; }

        for (int ivar = 0; ivar < _nvars[1]; ivar++) {
          double dtxJxW_g = JxW_g[1] * _bc_vol[i + (ivar + _nvars[2]) * NDOF_FEM];
          int index = i + ivar * _el_dof[1] + _el_dof[2] * _nvars[2];

          // Rhs Assemblying  ----------------------------------------------------------
          if (mode == 1) {
            // rhs

            FeM(index) += dtxJxW_g * (_ub_g[1][ivar] * phii_g / _dt  // time
                                      + 1. * phii_g                  // heat source
                                     );
          }

          // Matrix Assemblying ---------------------------
          for (int j = 0; j < _el_dof[1]; j++) {
            double phij_g = _phi_g[1][j];
            int jndex = j + _el_dof[2] * _nvars[2] + ivar * _el_dof[1];
            double Lap = 0.;

            for (int idim = 0; idim < ndim; idim++) {
              dphijdx_g[1][idim] = _dphi_g[1][j + idim * _el_dof[1]];
              Lap += dphijdx_g[1][idim] * dphiidx_g[1][idim];  // Laplacian
            }

            // energy-equation
            KeM(index, jndex) += dtxJxW_g * (phii_g * phij_g / _dt  // time term
                                             + Lap                  // diff
                                            );
          }
        }  //// linear variable cycle
      }    // ----------------------------------END LINEAR

      for (int i = 0; i < _el_dof[0]; i++) {  //    --- Piecewise ---
        // set up row i
        const double phii_g = _phi_g[0][i];

        for (int ivar = 0; ivar < _nvars[0]; ivar++) {
          double dtxJxW_g = JxW_g[0] * _bc_vol[i + (ivar + _nvars[1] + _nvars[2]) * NDOF_FEM];
          int index = i + ivar * _el_dof[0] + _el_dof[2] * _nvars[2] + _el_dof[1] * _nvars[1];

          // Rhs Assemblying  ----------------------------------------------------------
          if (mode == 1) {
            // rhs
            FeM(index) += dtxJxW_g * (
                                         //                        _ub_g[2][0]*phii_g
                                         _ub_g[0][ivar] * phii_g / _dt              // time
                                         + _data_eq[2].ub[0] * phii_g * (1 - ivar)  // source x
                                         + _data_eq[2].ub[0] * _data_eq[2].ub[1] * _data_eq[2].ub[2] *
                                               phii_g * ivar  // source x*y*z
                                     );
          }

          // Matrix Assemblying ---------------------------
          for (int j = 0; j < _el_dof[0]; j++) {
            double phij_g = _phi_g[0][j];
            int jndex = j + ivar * _el_dof[0] + _el_dof[2] * _nvars[2] + _el_dof[1] * _nvars[1];

            KeM(index, jndex) += dtxJxW_g * (phii_g * phij_g / _dt);  // time term
          }
        }  ////  variable cycle
      }    // ----------------------------------END piecewise

    }  // end of the quadrature point qp-loop ***********************

    // ======================================================================
    // ====================== boundary ======================================
    // ======================================================================

    for (int iside = 0; iside < el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int idof = 0; idof < NDOF_FEMB; idof++) {  // idof -> boundary node
          sur_toply[idof] =
              _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside];  // use map to find global node
          int idofb = sur_toply[idof];                              // idofb -> element node
          elb_conn[idof] = el_conn[idofb];                          // connectivity vector

          for (int idim = 0; idim < DIMENSION; idim++) {
            xxb_qnds[idim * NDOF_FEMB + idof] = xx_qnds[idim * NDOF_FEM + idofb];  // get boundary coordinates
            _data_eq[2].ub[idim * NDOF_FEM + idof] = xxb_qnds[idim * NDOF_FEMB + idof];
          }
        }

        for (int iql = 2; iql < 3; iql++)  // --- QUADRATIC ---
          for (int ivar = 0; ivar < _nvars[iql]; ivar++) {
            // Dirichlet boundary conditions  ***********************************
            if (_bc_vol[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM] == 0) {
              //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))
              int bc_s = (int)_bc_bd[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM];  // b cond
              double det = _fe[iql]->JacSur(elb_ngauss - 1, xxb_qnds, InvJac[2]);  // jacobian
              double Ipenalty = det / _dt;                                         // Dirichlet bc flag

              // local boundary loop   ---------------------------------------
              for (int lb_node = 0; lb_node < elb_ndof[iql]; lb_node++) {
                int index = sur_toply[lb_node] + ivar * NDOF_FEM;  // local vol index
                int gl_node_bd = elb_conn[lb_node];
                // flag setup (\int bc_var*T+bc_var*val)
                int bc_val = (int)((bc_s & 2) >> 1);  // (1?)non homogeneous
                int bc_var = (int)(bc_s % 2);         // (?1) variable

                // Assemblying  Matrix & rhs
                if (mode == 1) {
                  FeM(index) += bc_val * Ipenalty * _data_eq[iql].ub[DA_f * NDOF_FEM + index];
                  FeM(index) += (1 - bc_val) * bc_var * Ipenalty * (0.);
                  FeM(index) += Ipenalty * (1.);
                }

                KeM(index, index) += Ipenalty;  //  Dirichlet bc
              }                                 // lb_node -end  local boundary loop -------------------------
            }                                   // end if Dirichlet  boundary conditions

            // **********************************************************************

            else if (_bc_vol[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM] != 0) {
              // Non homogenous Neumann boundary conditions  ***********************************

              // gaussian integration loop (n_gauss)
              // -----------------------------------------------
              for (int qp = 0; qp < elb_ngauss; qp++) {
                // quad/linear  [2]=quad [1]=linear------------------------------------
                det[iql] = _fe[iql]->JacSur(qp, xxb_qnds, InvJac[2]);      // local coord _phi_g and jac
                JxW_g[iql] = det[iql] * _fe[iql]->_weight1[ndim - 2][qp];  // weight
                _fe[iql]->get_phi_gl_g(ndim - 1, qp, _phi_g[iql]);         // global coord _phi_g

#ifdef AXISYM  // axisymmetric  (index ->0)
                interp_el_sol(_data_eq[iql].ub, 0, DIMENSION, _phi_g[iql], elb_ndof[iql], _ub_g[iql]);
                JxW_g[iql] *= _ub_g[iql][0];

#endif

                // local side loop (over the node face)
                for (int lsi_node = 0; lsi_node < elb_ndof[iql]; lsi_node++) {
                  // set up row i
                  const double phii_g = _phi_g[iql][lsi_node];        // boundary test function
                  int index = sur_toply[lsi_node] + ivar * NDOF_FEM;  // local element index
                  int bc_s = (int)_bc_bd[index];
                  int bc_v = (int)_bc_vol[index];
                  double dtxJxW_g = JxW_g[iql] * bc_v;  // Neumann bc flag and Jac

                  // flag setup: +++++++++++++++++++++++++++++++++++++++
                  //  \int (bc_var*T+bc_val*val)ds
                  int bc_val = (int)((bc_s & 2) >> 1);  // (1?) non-homogeneous
                  int bc_var = (int)(bc_s % 2);         // (?1) tau=A*variable
                  // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                  // Assemblying rhs ----------------------------
                  if (mode == 1) { FeM(index) += bc_val * (1 - bc_var) * dtxJxW_g * phii_g * 1.; }

                  // Assemblying Matrix ---------------------------------
                  for (int lsj_node = 0; lsj_node < elb_ndof[iql]; lsj_node++) {
                    int jndex = sur_toply[lsj_node] + ivar * NDOF_FEM;
                    KeM(index, jndex) +=
                        dtxJxW_g * bc_var * phii_g * _phi_g[iql][lsj_node];  // Robin bc  (k*dt/dn = h*(-T))
                  }  // end j  ---------------------------------------

                }  // i   +++++++++++++++++++++++++++++++
              }    // end of the quadrature point qp-loop **********************

            }  // Neumann non homog

          }  // ---END QUADRATIC ---

        //        } //end if side

        for (int iql = 1; iql < 2; iql++) {  // --- LINEAR ---
          for (int ivar = 0; ivar < _nvars[iql]; ivar++) {
            // Dirichlet boundary conditions  ***********************************
            //         int dn_flag=0;
            //            for (int lb_node=0; lb_node< elb_ndof[iql]; lb_node++)
            //       if (_bc_vol[sur_toply[lb_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM] ==1) dn_flag++;

            //       if (dn_flag< 2 ){ // Dirichlet  boundary conditions <-- ONLY THIS
            //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))

            double det = _fe[iql]->JacSur(0, xxb_qnds, InvJac[2]);  // jacobian
            double Ipenalty = det / _dt;                            // Dirichlet bc flag

            // local boundary loop   ---------------------------------------
            for (int lb_node = 0; lb_node < elb_ndof[iql]; lb_node++) {
              if (_bc_vol[sur_toply[lb_node] + (ivar + _nvars[2]) * NDOF_FEM] == 0) {
                int index_bd = sur_toply[lb_node] + (ivar + _nvars[2]) * NDOF_FEM;  // local vol index
                int index_fem = sur_toply[lb_node] + ivar * _el_dof[iql] + _nvars[2] * NDOF_FEM;
                // flag setup (\int bc_var*T+bc_var*val)
                int bc_s = (int)_bc_bd[index_bd];  // b cond
                //                             int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                int bc_var = (int)(bc_s % 2);  // (?1) variable

                // Assemblying  Matrix & rhs
                if (mode == 1) {
                  FeM(index_fem) +=
                      (1 - bc_var) * Ipenalty * _data_eq[iql].ub[ivar * NDOF_FEM + sur_toply[lb_node]];
                  FeM(index_fem) += bc_var * Ipenalty * (5.);
                }

                KeM(index_fem, index_fem) += Ipenalty;  //  Dirichlet bc
              }                                         // end if on bc_vol
            }  // lb_node -end  local boundary loop -------------------------

            //          } // end if Dirichlet  boundary conditions
            // **********************************************************************

            //         else if (_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {
            //           if (dn_flag>0){
            // Non homogenous Neumann boundary conditions  ***********************************

            // gaussian integration loop (n_gauss)
            // -----------------------------------------------
            //           for (int qp=0; qp<  elb_ngauss; qp++) {
            //
            //             // quad/linear  [2]=quad [1]=linear------------------------------------
            //             det[iql]  = _fe[iql]->JacSur(qp,xxb_qnds);    // local coord _phi_g and jac
            //             JxW_g[iql]=det[iql]*_fe[iql]->_weight1[ndim-2][qp];// weight
            //             _fe[iql]->get_phi_gl_g(ndim-1,qp,_phi_g[iql]);   // global coord _phi_g
            //
            //
            // #ifdef AXISYM   // axisymmetric  (index ->0)
            //             interp_el_sol(_data_eq[iql].ub,0,DIMENSION,_phi_g[iql],elb_ndof[iql],_ub_g[iql]);
            //             JxW_g[iql]  *=_ub_g[2][0];
            //
            // #endif
            //             // local side loop (over the node face)
            //             for (int lsi_node=0; lsi_node< elb_ndof[iql]; lsi_node++) {
            //
            //               // set up row i
            //               const double phii_g=_phi_g[iql][lsi_node]; // boundary test function
            //               int index= sur_toply[lsi_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM; // local
            //               element index int bc_s=(int)_bc_bd[index]; int bc_v=(int)_bc_vol[index]; double
            //               dtxJxW_g=JxW_g[iql]*bc_v; // Neumann bc flag and Jac
            //
            //               // flag setup: +++++++++++++++++++++++++++++++++++++++
            //               //  \int (bc_var*T+bc_val*val)ds
            //               int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
            //               int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
            //               // ++++++++++++++++++++++++++++++++++++++++++++++++++++
            //
            //               // Assemblying rhs ----------------------------
            //               if (mode == 1){
            //             FeM(index) += bc_val*(1-bc_var)*dtxJxW_g*phii_g*1.;
            //        }
            //
            //               // Assemblying Matrix ---------------------------------
            //               for (int lsj_node=0; lsj_node< elb_ndof[iql];  lsj_node++) {
            //    int jndex=sur_toply[lsj_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM;
            //                  KeM(index,jndex) += dtxJxW_g*bc_var*phii_g*_phi_g[iql][lsj_node]; // Robin bc
            //                  (k*dt/dn = h*(-T))
            //               }// end j  ---------------------------------------
            //
            //             }// i   +++++++++++++++++++++++++++++++
            //           } // end of the quadrature point qp-loop **********************
            //         } // Neumann non homog
          }  // end for variables
        }    // ---END LINEAR ---

        for (int iql = 0; iql < 1; iql++) {  // --- PIECEWISE ---
          for (int ivar = 0; ivar < _nvars[iql]; ivar++) {
            // Dirichlet boundary conditions  ***********************************
            double det = _fe[2]->JacSur(elb_ngauss - 1, xxb_qnds, InvJac[2]);  // jacobian
            double Ipenalty = det / _dt;                                       // Dirichlet bc flag

            if (_bc_vol[(ivar + _nvars[1] + _nvars[2]) * NDOF_FEM] == 0) {
              int index =
                  ivar * _el_dof[iql] + _nvars[1] * _el_dof[1] + _nvars[2] * _el_dof[2];  // local vol index

              // Assemblying  Matrix & rhs
              if (mode == 1) {
                FeM(index) += Ipenalty * 0.;  //_data_eq[iql].ub[ivar*NDOF_FEM];
              }

              KeM(index, index) += Ipenalty;  //  Dirichlet bc
            }                                 // end if on bc_vol
          }                                   // end loop variables
        }                                     // -- END PIECEWISE --

      }  // end if side

    }  // ======================  end for boundary =======================================

    //     std::cout << KeM << "\n";
    //     std::cout << FeM << "\n";
    /// e) Global assemblying energy equation
    A[Level]->add_matrix(KeM, el_dof_indices);  // global matrix

    if (mode == 1) {
      b[Level]->add_vector(FeM, el_dof_indices);  // global rhs
    }

  }  // end of element loop

  // clean
  el_dof_indices.clear();
  A[Level]->close();

  //      A[Level]->print();
  if (mode == 1) { b[Level]->close(); }

#ifdef PRINT_INFO
  std::cout << " Matrix Assembled(DA)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

  return;
}
