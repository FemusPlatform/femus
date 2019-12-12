// lib include ------------------
// #include <limits>
// #include <sstream>

// configure file -------------
// #include "MGFE_conf.h"
// #include "Printinfo_conf.h"
// #include "Solverlib_conf.h"
// class include
#include "MGSolverBase.h"
#ifdef TWO_PHASE
#include "MGSolverCC.h"
#endif
// local inlude -----------------
#include "MGMesh.h"
// #include "MGSystem.h"
// #include "MGUtils.h"

// #include "MGEquationsSystem.h"
#include "MGFEMap.h"

// alg include-------------------
#include "linear_solverM.h"
#include "numeric_vectorM.h"
// #include "sparse_MmatrixM.h"
// #include "sparse_matrixM.h"

// #include <mpi.h>  //for MPI_COMM_WORLD
// #include "petsc_macroM.h"

// ==================================================
//           SOLUTION FUNCTIONS
//======================================================
// =================================================================================
/// This function interpolates a vector field over the fem element
// It can be used for sol and deriv (not for 2deriv)
void MGSolBase::interp_el(
    const double uold[],    // node values <-
    const int ivar0,        // init variable  <-
    const int nvars,        // # of variables  <-
    const double phi[],     // shape functions  <-
    const int n_shape,      // # of shape functions  <-
    double u_int[],         // interpolated function ->
    int dim,                // deriv or
    const int sur_tpgly[],  // surface nodes topology <-
    const int el_ndof       // surface nodes topology <-
    ) const {               // =======================================
  for (int ivar = 0; ivar < nvars; ivar++) {
    for (int jdim = 0; jdim < dim; jdim++) u_int[ivar * dim + jdim] = 0.;  // set zero u_int
    for (int j = 0; j < n_shape; j++) {
      for (int jdim = 0; jdim < DIMENSION; jdim++)
        u_int[ivar * dim + jdim] += phi[j + jdim * n_shape] * uold[j + (ivar + ivar0) * NDOF_FEM];
    }
  }
  return;
}
// =================================================================================
/// This function interpolates a vector field over the fem element
void MGSolBase::interp_el_sol(
    const double uold_b[],  // node values <-
    const int ivar0,        // init variable  <-
    const int nvars,        // # of variables  <-
    const double phi[],     // shape functions  <-
    const int n_shape,      // # of shape functions  <-
    double uold[]           // interpolated function ->
    ) const {               // =======================================
  for (int ivar = ivar0; ivar < ivar0 + nvars; ivar++) {
    uold[ivar] = 0.;
    for (int eln = 0; eln < n_shape; eln++) { uold[ivar] += phi[eln] * uold_b[eln + ivar * n_shape]; }
  }
  return;
}

// =================================================================================
/// This function interpolates a vector field and its derivative over the fem element
void MGSolBase::interp_el_gdx(
    double uold_b[],      // node values <-
    const int ivar0,      // init variable  <-
    const int nvars,      // # of variables  <-
    const double dphi[],  // derivatives of the shape functions  <-
    const int n_shape,    // # of shape functions  <-
    double uold_dx[]      // interpolated derivatives ->
    ) const {             // =================================================================
  // All data vectors are NDOF_FEM long
  // variable loop
  //   uold_dx=[dujdx dujdy dujdz dvjdx dvjdy dvjdz dwjdx dwjdy dwjdz]
  for (int ivar = 0; ivar < nvars; ivar++) {
    for (int jdim = 0; jdim < DIMENSION; jdim++) uold_dx[ivar * DIMENSION + jdim] = 0.;  // set zero
    // interpolation with shape functions
    for (int j = 0; j < n_shape; j++) {
      for (int jdim = 0; jdim < DIMENSION; jdim++)
        uold_dx[ivar * DIMENSION + jdim] += dphi[j + jdim * n_shape] * uold_b[j + (ivar + ivar0) * NDOF_FEM];
    }
  }
  return;
}
// =================================================================================
/// This function interpolates a vector field and its second derivatives over the fem element
void MGSolBase::interp_el_gddx(
    double uold_b[],      // node values <-
    const int ivar0,      // init variable  <-
    const int nvars,      // # of variables  <-
    const double dphi[],  // derivatives of the shape functions  <-
    const int n_shape,    // # of shape functions  <-
    double uold_dx[]      // interpolated derivatives ->
    ) const {             // =================================================================
  // All data vectors are NDOF_FEM long
  // uold_dx=[ d2udxdx d2udxdy  d2udxdz  d2udydx d2udydy  d2udydz  d2udzdx d2udzdy  d2udzdz
  // d2vdxdx d2vdxdy  d2vdxdz  d2vdydx d2vdydy  d2vdydz  d2vdzdx d2vdzdy  d2vdzdz
  // d2wdxdx d2wdxdy  d2wdxdz  d2wdydx d2wdydy  d2wdydz  d2wdzdx d2wdzdy  d2wdzdz]
  // variable loop
  for (int ivar = 0; ivar < nvars; ivar++) {
    for (int jdim = 0; jdim < DIMENSION * DIMENSION; jdim++)
      uold_dx[ivar * DIMENSION * DIMENSION + jdim] = 0.;  // set zero
    // interpolation with shape functions
    for (int eln = 0; eln < n_shape; eln++) {
      for (int jdim = 0; jdim < DIMENSION * DIMENSION; jdim++)
        uold_dx[ivar * DIMENSION * DIMENSION + jdim] +=
            dphi[eln * DIMENSION * DIMENSION + jdim] * uold_b[eln + (ivar + ivar0) * NDOF_FEM];
    }
  }
  return;
}
// =================================================================================
/// This function interpolates a vector field over the fem element
void MGSolBase::interp_el_bd_sol(
    const double uold_b[],  // node values <-
    const int sur_tpgly[],  // surface nodes topology <-
    const int el_ndof,      // surface nodes topology <-
    const int ivar0,        // init variable  <-
    const int nvars,        // # of variables  <-
    const double phi[],     // shape functions  <-
    const int n_shape,      // # of shape functions  <-
    double uold[]           // interpolated function ->
    ) const {               // =======================================
  for (int ivar = 0; ivar < nvars; ivar++) {
    uold[ivar] = 0.;
    for (int eln = 0; eln < n_shape; eln++)
      uold[ivar] += phi[eln] * uold_b[sur_tpgly[eln] + (ivar + ivar0) * el_ndof];
  }

  return;
}

// =================================================================================
/// This function interpolates a vector field and its derivative over a boundary
void MGSolBase::interp_el_bd_gdx(
    const double uold_b[],  // node values <-
    const int sur_tpgly[],  // surface nodes topology <-
    const int el_ndof,      // surface nodes topology <-
    const int ivar0,        // init variable  <-
    const int nvars,        // # of variables  <-
    const double dphi[],    // derivatives of the shape functions  <-
    const int n_shape,      // # of shape functions  <-
    double uold_dx[]        // interpolated derivatives ->
    ) const {               // =================================================================
  // variable loop
  int dim_b = DIMENSION;
  for (int ivar = 0; ivar < nvars; ivar++) {  // loop over variables (u,v,w,...)
    for (int jdim = 0; jdim < dim_b; jdim++) uold_dx[ivar * dim_b + jdim] = 0.;  // set zero
    // interpolation with shape functions
    for (int eln = 0; eln < n_shape; eln++) {
      for (int jdim = 0; jdim < dim_b; jdim++)  // loop over directions (x,y,z,...)
        uold_dx[ivar * dim_b + jdim] +=
            dphi[eln + jdim * n_shape] * uold_b[sur_tpgly[eln] + (ivar + ivar0) * el_ndof];
    }
  }
  return;
}

// ========================================================================
//  SET FUNCTIONS
// ========================================================================
void MGSolBase::set_x_aux(int i, int k, double val) {
  x_aux[i]->set(k, val);
  return;
}  // set the field
void MGSolBase::set_d_aux(int i, int k, double val) {
  d_aux[i]->set(k, val);
  return;
}
void MGSolBase::set_sol(int i, int k, double val) {
  x_old[i][_NoLevels - 1]->set(k, val);
  return;
}
#ifdef TWO_PHASE_LIB
void MGSolBase::set_mgcc(MGSolCC& cc) { _msolcc = &cc; }
#endif

// ========================================================================
//  GET FUNCTIONS
// ========================================================================
double MGSolBase::get_x_aux(int i, int k) { return ((*x_aux[i])(k)); }
double MGSolBase::get_d_aux(int i, int k) { return ((*d_aux[i])(k)); }
double MGSolBase::get_sol(int i, int k) { return ((*x_old[i][_NoLevels - 1])(k)); }

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void MGSolBase::get_el_sol(
    const int i_step,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double uold[]         // element node values ->
    ) const {             // ==============================================================
  for (int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim

      const int kdof_top =
          _node_dof[_NoLevels - 1][el_conn[id] + (ivar + ivar0) * offset];  // dof from top level
      uold[id + (ivar + kvar0) * NDOF_FEM] = ((*x_old[i_step][_NoLevels - 1])(kdof_top));  // element sol
    }  // end quadratic ------------------------------------------------
  }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void MGSolBase::get_el_sol_piece(
    const int ivar0,   // initial variable  <-
    const int nvars,   // # of variables to get  <-
    const int el_nds,  // # of element nodes for this variable  <-
    const int iel,     // connectivity <-
    const int offset,  // offset for connectivity <-
    const int kvar0,   // offset  variable for  uold <-
    double uold[]      // element node values ->
    ) const {          // ==============================================================
  for (int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivar is like idim
      const int kdof_top =
          _node_dof[_NoLevels - 1][iel * el_nds + id + (ivar + ivar0) * offset];      // dof from top level
      uold[id + (ivar + kvar0) * NDOF_FEM] = ((*x_old[0][_NoLevels - 1])(kdof_top));  // element sol
    }  // end quadratic ------------------------------------------------
  }
  return;
}

// ==========================================================================================
/// This function gets  the solution  vector at the nodes of  an element.
void MGSolBase::get_el_nonl_sol(
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double uold[]         // element node values ->
    ) const {             // ==============================================================
  for (int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivarq is like idim
      const int kdof_top =
          _node_dof[_NoLevels - 1][el_conn[id] + (ivar + ivar0) * offset];  // dof from top level
      double val = ((*x_nonl[_NoLevels - 1])(kdof_top));
      uold[id + (kvar0 + ivar) * NDOF_FEM] = val;  // element sol
    }  // end quadratic ------------------------------------------------
  }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void MGSolBase::get_el_d_aux(
    const int istep,
    const int ivar0,      // initial variable  <-
    const int nvars,      // # of variables to get  <-
    const int el_nds,     // # of element nodes for this variable  <-
    const int el_conn[],  // connectivity <-
    const int offset,     // offset for connectivity <-
    const int kvar0,      // offset  variable for  uold <-
    double uold[]         // element node values ->
    ) const {             // ==============================================================
  for (int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for (int ivar = 0; ivar < nvars; ivar++) {  // ivarq is like idim
      const int kdof_top =
          _node_dof[_NoLevels - 1][el_conn[id] + (ivar + ivar0) * offset];               // dof from top level
      uold[id + (kvar0 + ivar) * NDOF_FEM] = ((d_aux[istep][_NoLevels - 1])(kdof_top));  // element sol
    }  // end quadratic ------------------------------------------------
    //     disp_old[_NoLevels-1]->print();
  }
  return;
}

/// ======================================================
/// This function change current dt
/// ======================================================

double MGSolBase::GetValue(int flag) {}

void MGSolBase::SetValue(double value) {}

void MGSolBase::SetValueVector(std::vector<double> value) {}

// ==========================================================================
/// This function copies values from "vec_from" to "vec_to" vectors.
/// Available vectors are ordered as: x->0, x_old->1, x_oold->2, x_ooold->3,
/// x_nonl->4, disp->5, disp_old->6,  disp_oold->7
// ==========================================================================
void MGSolBase::set_cp_vector(const int& vec_from, const int& vec_to) {
  for (int Level = 0; Level <= _NoLevels - 1; Level++) {
    const int offset = _mgmesh._NoNodes[Level];  // fine level # of nodes

    for (int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++) {
      int el_nds = NDOF_FEM;

      if (ivar >= _nvars[2]) {
        el_nds = NDOF_P;  // quad and linear
      }

      for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++)
        for (int iel = 0; iel < _mgmesh._off_el[0][iproc * _NoLevels + Level + 1] -
                                    _mgmesh._off_el[0][iproc * _NoLevels + Level];
             iel++) {
          int elem_gidx = (iel + _mgmesh._off_el[0][iproc * _NoLevels + Level]) * NDOF_FEM;

          for (int i = 0; i < el_nds; i++) {            // linear and quad
            int k = _mgmesh._el_map[0][elem_gidx + i];  // the global node
            double value = 0.;

            switch (vec_from) {
              case 1: value = (*x_old[0][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 2: value = (*x_old[1][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 3: value = (*x_old[2][Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 4: value = (*x_nonl[Level])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 5: value = (*d_aux[0])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 6: value = (*d_aux[1])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 7: value = (*d_aux[2])(_node_dof[_NoLevels - 1][k + ivar * offset]); break;
              case 8:
                value = (*d_aux[3])(_node_dof[_NoLevels - 1][k + ivar * offset]);
                break;
                //                case 9: value = (*x_oooold[Level])(_node_dof[_NoLevels - 1][k + ivar *
                //                offset]); break;
              default: std::cout << "Incorrect vec_from number in set_uoold function" << endl; break;
            }

            switch (vec_to) {
              case 1: x_old[0][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 2: x_old[1][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 3: x_old[2][Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 4: x_nonl[Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value); break;
              case 5:
                d_aux[0]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 6:
                d_aux[1]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 7:
                d_aux[2]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
              case 8:
                d_aux[3]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
                break;
                //
                //               case 9:
                //                 x_oooold[Level]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);
                //                 // set the field break;

              default: cout << "Incorrect vec_to number in set_uoold function" << endl; break;
            }
          }
        }
    }

    if (vec_from == 1 && vec_to == 3) x[Level]->localize(*x_old[1][Level]);  // for backward compatibility
  }

  int ndof_lev = 0;
  for (int pr = 0; pr < _mgmesh._iproc; pr++) {
    int delta =
        _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
    ndof_lev += delta;
  }
  return;
}

// ====================================================================
/// This function gets element dof indices using local-to-global map
// ====================================================================
void MGSolBase::get_el_dof_indices(
    const int Level,                                 /// Level                                (in)
    const int iel,                                   /// ELement number                       (in)
    const int el_conn[],                             /// ELement connectivity                 (in)
    const int el_dof[],                              /// Quadratic[2] Linear[1] Const[0] dofs (in)
    const int offset,                                /// Offset                               (in)
    std::map<int, std::vector<int>>& el_dof_indices  /// Global dof of iel (out)
    ) const {
  for (int id = 0; id < NDOF_FEM; id++) {
    // quadratic -------------------------------------------------
    if (id < el_dof[2])
      for (int ivar = 0; ivar < _nvars[2]; ivar++) {  // ivarq is like idim
        const int indx_loc_ql = id + ivar * el_dof[2];
        const int indx_glob = el_conn[id] + ivar * offset;

        el_dof_indices[iel][indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
      }  // end quadratic ------------------------------------------------

    //     // linear -----------------------------
    if (id < el_dof[1])
      for (int ivar = 0; ivar < _nvars[1]; ivar++) {  // ivarq is like idim
        const int indx_loc_ql = id + ivar * el_dof[1] + _nvars[2] * el_dof[2];
        const int indx_glob = el_conn[id] + (ivar + _nvars[2]) * offset;

        el_dof_indices[iel][indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
      }  // end quadratic ------------------------------------------------

    //     // piecewise -----------------------------
    if (id < el_dof[0])
      for (int ivar = 0; ivar < _nvars[0]; ivar++) {  // ivarq is like idim
        const int indx_loc_ql = id + ivar * el_dof[0] + _nvars[2] * el_dof[2] + _nvars[1] * el_dof[1];
        const int indx_glob = id + iel * el_dof[0] + (ivar + _nvars[2] + _nvars[1]) * offset;

        el_dof_indices[iel][indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
      }  // end piecewise ------------------------------------------------
  }
  return;
}
