#include "Precondtype_enum.h"
#include "Solverlib_conf.h"
#include "linear_solverM.h"

// C++ includes
#include <memory>
// Local Includes
// #include "auto_ptr.h"

#include "laspack_linear_solverM.h"
#include "petsc_linear_solverM.h"
// #include "trilinos_aztec_linear_solver.h"
// #include "auto_ptr.h"
#include "preconditionerM.h"

//------------------------------------------------------------------
// LinearSolver members

// =============================================================
std::unique_ptr<LinearSolverM> LinearSolverM::build(
    const ParallelM::Communicator& comm, const SolverPackageM solver_package) {
  // Build the appropriate solver
  switch (solver_package) {
#ifdef HAVE_LASPACKM
    case LASPACK_SOLVERSM: {
      std::unique_ptr<LinearSolverM> ap(new LaspackLinearSolverM(comm));
      return ap;
    }
#endif
#ifdef HAVE_PETSCM
    case PETSC_SOLVERSM: {
      std::unique_ptr<LinearSolverM> ap(new PetscLinearSolverM(comm));
      return ap;
    }
#endif
#ifdef HAVE_TRILINOS
    case TRILINOS_SOLVERSM: {
      std::unique_ptr<LinearSolverM> ap(new AztecLinearSolverM(comm));
      return ap;
    }
#endif

    default: std::cerr << "ERROR:  Unrecognized solver package: " << solver_package << std::endl; abort();
  }
  std::cerr << "ERROR:  Unrecognized solver package: " << solver_package << std::endl;
  abort();
  //   std::unique_ptr<LinearSolverM > ap(NULL);
  //   return ap;
}

// ============================================================
PreconditionerTypeM LinearSolverM::preconditioner_type() const {
  if (_preconditioner) return _preconditioner->type();
  return _preconditioner_type;
}

// ===========================================================
void LinearSolverM::set_preconditioner_type(const PreconditionerTypeM pct) {
  if (_preconditioner)
    _preconditioner->set_type(pct);
  else
    _preconditioner_type = pct;
}

// =============================================================
void LinearSolverM::attach_preconditioner(PreconditionerM* preconditioner) {
  if (this->_is_initialized) {
    std::cerr << "Preconditioner must be attached before the solver is initialized!" << std::endl;
    abort();
  }
  _preconditioner_type = SHELL_PRECONDM;
  _preconditioner = preconditioner;
}
// ======================================
void LinearSolverM::reuse_preconditioner(bool reuse_flag) { same_preconditioner = reuse_flag; }
// =================================================
void LinearSolverM::restrict_solve_to(
    const std::vector<int>* const dofs, const SubsetSolveModeM /*subset_solve_mode*/) {
  if (dofs != NULL) { std::cout << "libmesh_not_implemented"; }
}

// ===========================================
std::pair<int, double> LinearSolverM::adjoint_solve(
    SparseMatrixM& mat, NumericVectorM& sol, NumericVectorM& rhs, const double tol, const int n_iter) {
  // Log how long the linear solve takes.
  //     START_LOG("adjoint_solve()", "LinearSolver");

  // Take the discrete adjoint
  mat.close();
  mat.get_transpose(mat);
  // Call the solve function for the relevant linear algebra library and
  // solve the transpose matrix
  const std::pair<int, double> totalrval = this->solve(mat, sol, rhs, tol, n_iter, 0);
  // Now transpose back and restore the original matrix
  // by taking the discrete adjoint
  mat.get_transpose(mat);

  // Stop logging the nonlinear solve
  //     STOP_LOG("adjoint_solve()", "LinearSolver");
  return totalrval;
}
//------------------------------------------------------------------
// Explicit instantiations
// template class LinearSolver<Number>;
