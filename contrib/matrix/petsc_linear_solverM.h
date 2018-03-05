#ifndef __petsc_linear_solverMaa_h__
#define __petsc_linear_solverMaa_h__

#include "Solverlib_conf.h"

#ifdef HAVE_PETSCM

#include <mpi.h> //MPI_Comm_world

// Local includes
#include "linear_solverM.h"
#include "petsc_vectorM.h"
#include "petsc_matrixM.h"
#include "petsc_macroM.h"

// Petsc include files. 
EXTERN_C_FOR_PETSC_BEGIN
// #if PETSC_VERSION_LESS_THAN(2,2,0)
// #  include <petscsles.h>
// #else
#  include <petscksp.h>
// #endif
EXTERN_C_FOR_PETSC_END

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed for preconditioning
// 
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
// #if PETSC_VERSION_LESS_THAN(2,2,1)
//   typedef int PetscErrorCode;
//   typedef int PetscInt;
// #endif
#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  /// This function is called by PETSc to initialize the preconditioner ctx
  PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx);
  /// This function is called by PETSc to acctually apply the preconditioner ctx 
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y);
#else
  PetscErrorCode __libmesh_petsc_preconditioner_setup (PC);
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC, Vec x, Vec y);
#endif
} // end extern "C"

		
// ==========================================
/// This class provides an interface to PETSc iterative solvers 
class PetscLinearSolverM : public LinearSolverM
{// ============================================

   private:
     // data ---------------------------------
// #if PETSC_VERSION_LESS_THAN(2,2,0)  // SLES removed from >= PETSc 2.2.0
//   SLES _sles;///< Linear solver context   
// #endif
  PC _pc;  ///< Preconditioner context  
  KSP _ksp;///< Krylov subspace context 
  
  /// PETSc index set containing the dofs on which to solve (\p NULL  means solve on all dofs).
  IS _restrict_solve_to_is;
   /// PETSc index set, complement to \p _restrict_solve_to_is.  This
   /// will be created on demand by the method \p _create_complement_is().
  IS _restrict_solve_to_is_complement;
  
    /// Internal method that returns the local size of \p _restrict_solve_to_is.
  size_t _restrict_solve_to_is_local_size(void)const;
  
   /// Creates \p _restrict_solve_to_is_complement to contain all
   /// indices that are local in \p vec_in, except those that are
   ///contained in \p _restrict_solve_to_is.
  void _create_complement_is (const NumericVectorM &vec_in);
  
   /// If restrict-solve-to-subset mode is active, this member decides
   /// what happens with the dofs outside the subset.
  SubsetSolveModeM _subset_solve_mode;
  
public:
  // Constructor --------------------------------------
  ///  Constructor. Initializes Petsc data structures 
  PetscLinearSolverM (const ParallelM::Communicator &comm);
  /// Destructor.  
  ~PetscLinearSolverM ();
  /// Release all memory and clear data structures.
  void clear ();
  /// Initialize data structures if not done so already.
  void init ();
  /// Initialize data structures if not done so already plus much more
  void init (PetscMatrixM* matrix);
  
  // Solvers ------------------------------------------------------
  // ========================================================
//     std::pair< int, double> MGsolve(SparseMatrixM& matrix_in,
//     SparseMatrixM& precond_in, NumericVectorM& solution_in, NumericVectorM& rhs_in,  NumericVectorM& res_in,
//     SparseMMatrixM& prol_in, SparseMMatrixM& rest_in, const int nlevels, const int nup_its, const int ndown_its, 
//    const int mcoarse_its, const double tol);
  
  // ======================================================
    /// After calling this method, all successive solves will be
    ///  restricted to the given set of dofs, which must contain local
    ///  dofs on each processor only and not contain any duplicates.  This
   /// mode can be disabled by calling this method with \p dofs being a \p NULL pointer.
    virtual void restrict_solve_to (
      const std::vector<int>* const dofs,
      const SubsetSolveModeM subset_solve_mode=SUBSET_ZERO);
    
  /// Call the Petsc solver.  This function calls the method below, using the
  /// same matrix for the system and preconditioner matrices.    
  std::pair< int, double>   solve (
    SparseMatrixM  &matrix_in, 
    NumericVectorM &solution_in, 
    NumericVectorM &rhs_in,
    const double tol,
    const  int m_its,
    const int clearing  
)  {
    return this->solve(matrix_in, matrix_in, solution_in, rhs_in, tol, m_its, clearing);
  }

  /// This method allows you to call a linear solver while specifying
  /// the matrix to use as the (left) preconditioning matrix.  Note
  /// that the linear solver will not compute a preconditioner in this
  /// case, and will instead premultiply by the matrix you provide.
  /// In PETSc, this is accomplished by calling  PCSetType(_pc, PCMAT);
  /// before invoking KSPSolve().  Note: this functionality is not implemented
  /// in the LinearSolver class since there is not a built-in analog to this method for LasPack 
  std::pair< int, double>   solve (
         SparseMatrixM  &/*matrix*/,
	 SparseMatrixM  &/*preconditioner*/,
	 NumericVectorM &/*solution*/,
	 NumericVectorM &/*rhs*/,
	 const double /*tol*/,const  int /*Niter*/, const int clearing);  
  
  // =======================================================
  /// Call the Petsc solver.  It calls the method below, using the
  ///  same matrix for the system and preconditioner matrices.
  std::pair<int, double> adjoint_solve (
    SparseMatrixM  &matrix_in,
	 NumericVectorM &solution_in,
	 NumericVectorM &rhs_in,
	 const double tol,
         const  int m_its);
// This function solves a system whose matrix is a shell matrix.
//   std::pair< int, double>
//     solve (const ShellMatrix<T>& shell_matrix,
// 	   NumericVectorM& solution_in,
// 	   NumericVectorM& rhs_in,
// 	   const double tol,
// 	   const  int m_its);
  
//  /** This function solves a system whose matrix is a shell matrix, but
//   * a sparse matrix is used as preconditioning matrix, this allowing
//   * other preconditioners than JACOBI.
//   */
//   virtual std::pair< int, double>
//     solve (const ShellMatrix<T>& shell_matrix,
// 	   const SparseMatrixM& precond_matrix,
// 	   NumericVectorM& solution_in,
// 	   NumericVectorM& rhs_in,
// 	   const double tol,
// 	   const  int m_its);
// void MGSolve(std::vector<SparseMatrixM *>  &/*matrix_in*/,   // System Matrix
//	              std::vector<NumericVectorM *> &/*solution_in*/,// Solution vector
//		      std::vector<NumericVectorM *> &/*rhs_in*/,     // RHS vector
//		      Matrix */*P*/, Matrix */*R*/,  
//		      const double /*tol*/,const  int /*m_its*/){
// std::cout << "not implemented"; abort(); 
//}	    
// Returns -----------------------------------	      
  /// Returns the raw PETSc preconditioner context pointer.  This allows
  /// you to specify the PCShellSetApply() and PCShellSetSetUp() functions if you desire.
  PC pc() { this->init(); return _pc; }
  /// Returns the raw PETSc ksp context pointer.  This is useful if
  /// you are for example setting a custom convergence test with KSPSetConvergenceTest().
  KSP ksp() { this->init(); return _ksp; }
  /// Fills the input vector with residual norms from the latest iterative solve.
  void get_residual_history(std::vector<double>& hist);
  /// Returns just the initial residual for the solve just
  /// completed with this interface.  Use this method instead
  /// of the one above if you just want the starting residual and not the entire history.
  double get_initial_residual();

  // Print ----------------------------------------
  /// Prints a useful message about why the latest linear solve con(di)verged.
  virtual void print_converged_reason();
  
  // Setting --------------------------------------------
private:
  ///  Set the user-specified solver stored in \p _solver_type
  void set_petsc_solver_type ();
 // /// Internal function if shell matrix mode is used.
 // static PetscErrorCode _petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest);
 // /// Internal function if shell matrix mode is used.
 // static PetscErrorCode _petsc_shell_matrix_get_diagonal(Mat mat, Vec dest);
};


/*----------------------- functions ----------------------------------*/
// =================================================
inline PetscLinearSolverM::PetscLinearSolverM (
  const ParallelM::Communicator &comm
) :   LinearSolverM(comm),
    _restrict_solve_to_is(NULL),
    _restrict_solve_to_is_complement(NULL),
    _subset_solve_mode(SUBSET_ZERO) {
      
    if (this->n_processors() == 1)    this->_preconditioner_type = ILU_PRECONDM;
    else    this->_preconditioner_type = BLOCK_JACOBI_PRECONDM;
  
//   int i; MPI_Comm_size(MPI_COMM_WORLD,&i);  //TODO
// //  if (i == 1) this->_preconditioner_type = LU_PRECONDM;
//   if (i == 1) this->_preconditioner_type = ILU_PRECONDM;
// //  if (i == 1) this->_preconditioner_type = CHOLESKY_PRECONDM;
//   else  this->_preconditioner_type = BLOCK_JACOBI_PRECONDM;
// //  this->_preconditioner_type = ILU_PRECONDM;
}
// =============================================
inline PetscLinearSolverM::~PetscLinearSolverM (){this->clear ();}


// =============================
inline size_t PetscLinearSolverM::
_restrict_solve_to_is_local_size(void)const{
  assert(_restrict_solve_to_is);
  PetscInt s;
  int ierr = ISGetLocalSize(_restrict_solve_to_is,&s); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<size_t>(s);
}

// // ==========================================
inline
void PetscLinearSolverM::_create_complement_is (
  const NumericVectorM &  vec_in
  ){
  assert(_restrict_solve_to_is);
  if(_restrict_solve_to_is_complement==NULL)    {
      int ierr = ISComplement(_restrict_solve_to_is,
			      vec_in.first_local_index(),
			      vec_in.last_local_index(),
			      &_restrict_solve_to_is_complement);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}
#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
