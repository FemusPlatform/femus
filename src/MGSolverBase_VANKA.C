// lib include ------------------

// class include
#include "MGSolverBase.h"
#ifdef TWO_PHASE
#include "MGSolverCC.h"
#endif
// local inlude -----------------
#include "MGMesh.h"
#include "MGUtils.h"
#include "MGEquationsSystem.h"

// #define VANKA (0)

// ========================================================================
#ifdef VANKA
// ========================================================================
// ****************** HERE STARTS VANKA SECTION ***************************
// ========================================================================
#include "MeshExtended.h"
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "petsc_matrixM.h"
#include "petsc_vectorM.h"
//
// ====================================================================
/// This function does one multigrid step
double MGSolBase::Vanka_test(int Level) {
  // ====================================================================
  std::pair<int, double> rest(0, 0.);
  if (Level == 0) {  // coarse level ----------------------------------
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], 1.e-6, 40);
#ifdef PRINT_CONV
    std::cout << " Coarse res " << rest.second << " " << rest.first << std::endl;
#endif
  }  // --------------------------------------------------------------
  ;
  return 0;
}

// ====================================================================
/// This function does one multigrid step with Vanka
double MGSolBase::MGStep_Vanka(
    int Level,            // Level
    double Eps1,          // Tolerance
    int MaxIter,          // n iterations
    const int Gamma,      // Control V W cycle
    const int Nc_pre,     // n pre-smoothing cycles
    const int Nc_coarse,  // n coarse cycles
    const int Nc_post     // n post-smoothing cycles
) {
  // ====================================================================
  std::pair<int, double> rest(0, 0.);
  if (Level == 0) {  // coarse level ----------------------------------
    //       b[Level]->close(); double bNorm=b[Level]->l2_norm();
    //       x[Level]->close(); double xNorm=x[Level]->l2_norm();
    //       A[Level]->close(); double aNorm=A[Level]->linfty_norm();
    // coarse solution
    //    std::pair<int,double> rest1(0.,0);
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, 40);
    // coarse residual
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV
    std::cout << " Coarse res " << 40 << " " << res[Level]->l2_norm() << std::endl;
#endif

  }       // --------------------------------------------------------------
  else {  // fine levels

  // presmoothing (Nu1) ---------------------------------
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t start_time = std::clock();
#endif  //  TC +++++++++++++++
    int Nc_pre1 = Nc_pre / 4;
    if (Level < _NoLevels - 1) Nc_pre1 *= 2;
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, Nc_pre1);
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Pre Lev " << Level << " res " << res[Level]->l2_norm();
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
    // presmoothing residual
    //         res[Level]->resid(*b[Level],*x[Level],*A[Level]);
    // --------- end presmoothing (Nc_pre) ------------------------

    // restriction
    b[Level - 1]->matrix_mult(*res[Level], *Rst[Level - 1]);

    //  solving of system of equations for the residual on the coarser grid
    x[Level - 1]->close();
    x[Level - 1]->zero();
    for (int g = 1; g <= Gamma; g++)
      MGStep_Vanka(Level - 1, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);

    // interpolation of the solution from the coarser grid (projection)
    res[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);
    // adding the coarse solution
    x[Level]->add(*res[Level]);  //  *x[Level] +=*res[Level];

    // postsmooting (Nc_post) --------------------------------------------
#ifdef PRINT_TIME               //  TC +++++++++++++++
    start_time = std::clock();  //   initial set
#endif                          //  TC +++++++++++++++
    //  postsmooting
    int Nc_post1 = Nc_post / 4;
    if (Level < _NoLevels - 1) Nc_post1 *= 2;
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, Nc_post1);
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Post Lev " << Level << " res " << res[Level]->l2_norm();
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
    //  postsmooting residual
    //         res[Level]->resid(*b[Level],*x[Level],*A[Level]);
    // ----------------  end postsmoothing ---------------------------
  }
  // end cycle -------------------------------------
  res[Level]->close();
  double norm2 = res[Level]->l2_norm();
  //    std::cout<< " True res l2norm (not prec) " << norm2<< std::endl;
  return norm2;
}

// ========================================================
/// This function solves a Vanka step  Ax=b

void MGSolBase::Vanka_solve(
    int Level,                    /// Level
    SparseMatrixM& matrix_in,     /// Matrix A
    NumericVectorM& solution_in,  /// Previous solution x
    NumericVectorM& rhs_in,       /// Rhs b
    const double tol,             /// Tolerance
    const int m_its               /// n smoothing cycles
) {
  // Get el_dof to know the linear system size
  int el_dof[3];
  el_dof[0] = ((_nvars[0] > 0) ? NDOF_K : 0);    // Lagrange piecewise constant variables
  el_dof[1] = ((_nvars[1] > 0) ? NDOF_P : 0);    // Lagrange piecewise linear variables
  el_dof[2] = ((_nvars[2] > 0) ? NDOF_FEM : 0);  // Lagrange piecewise linear variables
  const int system_size = el_dof[2] * _nvars[2] + el_dof[1] * _nvars[1] + el_dof[0] * _nvars[0];  // ndofs
  // elem
  std::map<int, std::vector<int>> el_dof_indices;  // element dof vector
  std::map<int, std::vector<int>> subdom_el_map;   // subdom map vector

  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];

  // Now we use the mat group to create the Vanka blocks
  MeshExtended* ext_mesh = dynamic_cast<MeshExtended*>(&_mgmesh);
  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    int mat = ext_mesh->_mat_id[iel + nel_b];
    if (mat % 4 == 1) { subdom_el_map[0].push_back(iel); }
    if (mat % 4 == 2) {
      subdom_el_map[0].push_back(iel);
      subdom_el_map[1].push_back(iel);
    }
    if (mat % 4 == 3) {
      subdom_el_map[1].push_back(iel);
      subdom_el_map[2].push_back(iel);
    }
    if (mat % 4 == 0) { subdom_el_map[2].push_back(iel); }
  }
  int size;  // linear system size
  DenseMatrixM Mat;
  DenseVectorM sol, loc_rhs;

  // Petsc variables init
  PetscErrorCode ierr;
  PetscInt ncols;
  const PetscInt* cols;
  const PetscScalar* vals;

  // Make sure the data passed in are really of Petsc types
  PetscMatrixM* matrix = libmeshM_cast_ptr<PetscMatrixM*>(&matrix_in);
  PetscVectorM* solution = libmeshM_cast_ptr<PetscVectorM*>(&solution_in);
  PetscVectorM* rhs = libmeshM_cast_ptr<PetscVectorM*>(&rhs_in);

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();
  solution->close();
  rhs->close();

  /// a) Set up
  // geometry -----------------------------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[Level];  // mesh nodes
  int el_conn[NDOF_FEM];                       // element connectivity
  double xx_qnds[DIMENSION * NDOF_FEM];
  const int total_offset = _mgmesh._NoNodes[_NoLevels - 1];

  vector<double> x_loc[offset];
  solution->localize(*x_loc);
  vector<double> b_loc[offset];
  rhs->localize(*b_loc);

  //====== Get el_dof_indices map for all elements at this level===============================
  for (int iel = 0; iel < (nel_e - nel_b); iel++) { el_dof_indices[iel].resize(system_size); }

  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, xx_qnds);  // get connectivity and coord
    get_el_dof_indices(Level, iel, el_conn, el_dof, total_offset, el_dof_indices);
  }

  //======End Get el_dof_indices map for all elements at this level===============================

  int nsm = m_its;
  for (int ismooth = 0; ismooth <= nsm; ismooth++) {                       // smooth loop
    for (int n_subdom = 0; n_subdom < subdom_el_map.size(); n_subdom++) {  // loop over subdomains
      std::vector<int> subdom_indx;  // subdom restricted indx with no repeated dof
      subdom_indx.resize(system_size * subdom_el_map[n_subdom].size());          // max number
      for (int loc_el = 0; loc_el < subdom_el_map[n_subdom].size(); loc_el++) {  // el in subdom
        //               std::cout<<"subdom "<<n_subdom<<" iel "<<subdom_el_map[n_subdom][loc_el]<<std::endl;
        for (int i = 0; i < system_size; i++) {
          subdom_indx[i + loc_el * system_size] = el_dof_indices[subdom_el_map[n_subdom][loc_el]][i];
        }
      }
      std::sort(subdom_indx.begin(), subdom_indx.end());  // removing duplicates
      auto last = std::unique(subdom_indx.begin(), subdom_indx.end());
      subdom_indx.erase(last, subdom_indx.end());
      size = subdom_indx.size();
      PetscScalar matr[size * size];  // Resize with actual linear system size
      PetscInt idx[size];
      Mat.resize(size, size);
      sol.resize(size);
      loc_rhs.resize(size);
      for (int i = 0; i < size; i++) { idx[i] = subdom_indx[i]; }

      ierr = MatGetValues(matrix->mat(), size, idx, size, idx, matr);  // extract values from global matrix

      for (int inode = 0; inode < size; inode++) {  // loop over ndof
        //         std::cout << el_conn[inode] << " " ;
        ierr = MatGetRow(matrix->mat(), idx[inode], &ncols, &cols, &vals);

        loc_rhs(inode) = (*b_loc)[idx[inode]];  // RHS from assembly
        for (int i = 0; i < size; i++) {
          loc_rhs(inode) += matr[i + size * inode] * (*x_loc)[idx[i]];  // sum diagonal values
        }
        for (int i = 0; i < ncols; i++) {
          loc_rhs(inode) -= vals[i] * (*x_loc)[cols[i]];  // subtract all row values
        }

        ierr = MatRestoreRow(matrix->mat(), idx[inode], &ncols, &cols, &vals);

      }  // end inode

      Mat.zero();
      for (int j = 0; j < size; j++)
        for (int i = 0; i < size; i++) { Mat(j, i) = matr[i + j * size]; }

      Mat.lu_solve(loc_rhs, sol);  // Solve Mat*sol=loc_rhs with libmesh lu decomposition

      for (int i = 0; i < size; i++) { (*x_loc)[idx[i]] = sol(i); }

      if (ismooth == nsm) {  // We save the solution only after the nsm-th smoothing cycle
        for (int i = 0; i < size; i++) { (*x[Level]).set(idx[i], sol(i)); }
      }
    }  // end subdom loop
  }    // end ismooth loop
  return;
}
#endif
