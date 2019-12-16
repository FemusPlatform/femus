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

#include "linear_solverM.h"

// ========================================================================
#ifdef VANKA
// ========================================================================
// ****************** HERE STARTS VANKA SECTION ***************************
// ========================================================================
#include "dense_vectorM.h"
#include "dense_matrixM.h"
#include "petsc_vectorM.h"
#include "petsc_matrixM.h"
#include "MeshExtended.h"
//
// ====================================================================
/// This function does one multigrid step
double MGSolBase::Vanka_test(int Level) {
  // ====================================================================
  std::pair<int, double> rest(0, 0.);
  if (Level == 0) {  // coarse level ----------------------------------
                     //         std::vector<int> eldofs;
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
    // coarse solution
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, MaxIter, true);
    // coarse residual
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV
    std::cout << " Coarse res " << MaxIter << " " << res[Level]->l2_norm() << std::endl;
#endif

  }       // --------------------------------------------------------------
  else {  // fine levels

    // presmoothing (Nu1) ---------------------------------
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t start_time = std::clock();
#endif  //  TC +++++++++++++++
    int Nc_pre1 = Nc_pre;
    if (Level < _NoLevels - 1) { Nc_pre1 *= 2; }
    bool direct = /*(Level<5)?false:*/ true;
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, Nc_pre1, direct);
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Pre Lev " << Level << " res " << res[Level]->l2_norm();
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
    // --------- end presmoothing (Nc_pre) ------------------------

    // restriction
    b[Level - 1]->matrix_mult(*res[Level], *Rst[Level - 1]);

    //  solving of system of equations for the residual on the coarser grid
    x[Level - 1]->close();
    x[Level - 1]->zero();
    for (int g = 1; g <= Gamma; g++) {
      MGStep_Vanka(Level - 1, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post);
    }

    // interpolation of the solution from the coarser grid (projection)
    res[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);
    // adding the coarse solution
    x[Level]->add(*res[Level]);  //  *x[Level] +=*res[Level];

    // postsmooting (Nc_post) --------------------------------------------
#ifdef PRINT_TIME               //  TC +++++++++++++++
    start_time = std::clock();  //   initial set
#endif                          //  TC +++++++++++++++
    //  postsmooting
    int Nc_post1 = Nc_post;
    if (Level < _NoLevels - 1) { Nc_post1 *= 2; }
    Vanka_solve(Level, *A[Level], *x[Level], *b[Level], Eps1, Nc_post1, direct);
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Post Lev " << Level << " res " << res[Level]->l2_norm();
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
        // ----------------  end postsmoothing ---------------------------
  }
  // end cycle -------------------------------------
  res[Level]->close();
  double norm2 = res[Level]->l2_norm();
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
    const int m_its,              /// n smoothing cycles
    bool isdirect) {
  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1];
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc];

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

  vector<double> x_loc[offset];
  solution->localize(*x_loc);
  vector<double> b_loc[offset];
  rhs->localize(*b_loc);

  int nsm = m_its / 2;
  for (int ismooth = 0; ismooth <= nsm; ismooth++) {             // smooth loop
    for (int iel = 0; iel < (nel_e - nel_b); iel++) {            // loop over subdomains
      _mgmesh.get_el_nod_conn(0, Level, iel, el_conn, xx_qnds);  // get connectivity and coord
      std::vector<int> subdom_indx;  // subdom restricted indx with no repeated dof
      for (int idof = 0; idof < NDOF_P; idof++) {
        ierr = MatGetRow(matrix->mat(), el_conn[idof], &ncols, &cols, NULL);
        for (int icol = 0; icol < ncols; icol++) subdom_indx.push_back(cols[icol]);
        ierr = MatRestoreRow(matrix->mat(), el_conn[idof], &ncols, &cols, &vals);
      }

      std::sort(subdom_indx.begin(), subdom_indx.end());  // removing duplicates
      auto last = std::unique(subdom_indx.begin(), subdom_indx.end());
      subdom_indx.erase(last, subdom_indx.end());
      size = subdom_indx.size();

      PetscInt idx[size];
      for (int i = 0; i < size; i++) { idx[i] = subdom_indx[i]; }
      if (isdirect == true) {
        //                 cout<<"diretto"<<Level<<endl;
        PetscScalar matr[size * size];  // Resize with actual linear system size
        Mat.resize(size, size);
        sol.resize(size);
        loc_rhs.resize(size);
        ierr = MatGetValues(matrix->mat(), size, idx, size, idx, matr);  // extract values from global matrix

        for (int inode = 0; inode < size; inode++) {  // loop over ndof
          //               std::cout << el_conn[inode] << " " ;
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
      }
      if (isdirect == false) {
        _solver[Level]->restrict_solve_to(nullptr, SUBSET_DONT_TOUCH);
        _solver[Level]->restrict_solve_to(&subdom_indx, SUBSET_DONT_TOUCH);
        const int clearing = 0;
        std::pair<int, double> rest(0, 0.);
        rest = _solver[Level]->solve(matrix_in, solution_in, rhs_in, tol, m_its, clearing);
        //                 cout<<size<<endl;
        for (int i = 0; i < size; i++) { (*x_loc)[idx[i]] = solution_in(idx[i]); }

        if (ismooth == nsm) {  // We save the solution only after the nsm-th smoothing cycle
          for (int i = 0; i < size; i++) { (*x[Level]).set(idx[i], (*x_loc)[idx[i]]); }
        }
      }
    }  // end subdom loop
  }    // end ismooth loop
  return;
}
#endif
