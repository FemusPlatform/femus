#ifndef __linear_solverM_h__
#define __linear_solverM_h__


#include "Printinfo_conf.h"
#include "Solverlib_conf.h"

#include "SolverPackage_enum.h"
#include "Precondtype_enum.h"
#include "Solvertype_enum.h"
#include "parallel_objectM.h"
#include "enum_subset_solve_modeM.h"
// #include "VankaSmoothertype_enum.h"

// #include "enum_solver_package.h"
// #include "enum_solver_type.h"
// #include "enum_preconditioner_type.h"

// #include "reference_counted_object.h"
// #include "libmesh.h"
// C++ includes
#include <cstddef>
#include <vector>
#include <memory>

// forward declarations
// template <typename T> class AutoPtr;  //
class SparseMatrixM;
class NumericVectorM;
// template <typename T> class ShellMatrix;  //
class PreconditionerM;
// #include "matrix.h"  //TODO does it need this? oh, it's from LASPACK!
class ParallelObjectM;
// ================================================
// This class provides a uniform interface for linear solvers.  This base
// class is overloaded to provide linear solvers from different packages
// like PETSC or LASPACK.
// ===============================================
class LinearSolverM:
#ifdef LM_REFCOUNT
    public ReferenceCountedObject<LinearSolverM > ,
#endif
    public ParallelObjectM
{
    // ================================
    // DATA
    // ==============================
protected:
    /// Enum stating which type of iterative solver to use.
    SolverTypeM _solver_type;

    /// Flag indicating if the data structures have been initialized.
    bool _is_initialized;

    /// Enum statitng with type of preconditioner to use.
    PreconditionerTypeM _preconditioner_type;

    //   VankaSmootherTypeM _vanka_smoother_type;

    /// Holds the Preconditioner object to be used for the linear solves.
    PreconditionerM * _preconditioner;

    /// Boolean flag to indicate whether we want to use an identical preconditioner to the previous solve.
    bool same_preconditioner;
public:



    // =================================
    // CONSTR/DESTR
    // ================================
    ///  Constructor. Initializes Solver data structure
    LinearSolverM (
        const ParallelM::Communicator &comm_in   ///< parallel communicator <-
    );
    /// Destructor.
    virtual ~LinearSolverM();
    /// Release all memory and clear data structures.
    virtual void clear() {}
    /// Initialize data structures if not done so already.
    virtual void init() = 0;

    /// Builds a \p LinearSolverM using the linear solver in \p solver_package
    static std::unique_ptr<LinearSolverM > build(
        const ParallelM::Communicator  &comm_in,
        const SolverPackageM solver_package =LSOLVER
    );

    virtual void reuse_preconditioner(bool );

    // =================================
    // SETTING FUNCTIONS
    // ================================
    /// Sets the type of solver to use.
    void set_solver_type (const SolverTypeM st)  {
        _solver_type = st;
    }

    /// Sets the type of preconditioner to use.
    void set_preconditioner_type (const PreconditionerTypeM pct);

    /// Sets the type of vanka smoother to use.
//   void set_vanka_smoother_type (const VankaSmootherTypeM vsm) {_vanka_smoother_type = vsm; }

    /// Attaches a Preconditioner object to be used
    void attach_preconditioner(PreconditionerM * preconditioner);

    // =================================
    // RETURN FUNCTIONS
    // ================================
    /// ** @returns true if the data structures are
    bool get_same_preconditioner();
    /// @returns true if the data structures are
    bool initialized () const {  return _is_initialized;}
    /// Returns the type of solver to use.
    SolverTypeM solver_type () const { return _solver_type; }
    /// Returns the type of preconditioner to use.
    PreconditionerTypeM preconditioner_type () const;


    /// Returns the type of Vanka smoother to use.
//   VankaSmootherTypeM vanka_smoother_type () const { return _vanka_smoother_type; }

    // =================================
    // SOLVE FUNCTIONS
    // ================================

//   //sandro
//   /// This function calls the solver
//   virtual std::pair< int, double> solve(
//     SparseMatrixM&/*Sys Mtrx*/,NumericVectorM&/*sol_vec*/, NumericVectorM&/*rhs*/,
//     const double /*tol*/,const  int/*NIter*/) = 0;
//
//   ///   *This function calls the solver with preconditioner
//   virtual std::pair< int, double> solve(
//     SparseMatrixM&/*Sys Mtrx*/,SparseMatrixM&/*Prec Mtrx*/,
//     NumericVectorM&/*sol_vec*/,NumericVectorM&/*rhs*/,
//     const double /*tol*/,const  int/*NIter*/) = 0;

    /// Function to solve the adjoint system. Note that this method
    /// will compute the preconditioner from the system matrix. This is not a pure virtual
    // function and is defined linear_solver.C
    virtual std::pair< int, double> adjoint_solve (      SparseMatrixM&,  // System Matrix
            NumericVectorM&, // Solution vector
            NumericVectorM&, // RHS vector
            const double,      // Stopping tolerance
            const int); // N. Iterations

    /// This function calls the solver
    virtual std::pair< int, double> solve (      SparseMatrixM&,  // System Matrix
            NumericVectorM&, // Solution vector
            NumericVectorM&, // RHS vector
            const double,      // Stopping tolerance
            const  int,
            const int
                                          ) = 0; // N. Iterations

    ///   *This function calls the solver with preconditioner
    virtual std::pair< int, double> solve (      SparseMatrixM&,  // System Matrix
            SparseMatrixM&,  // Preconditioning Matrix
            NumericVectorM&, // Solution vector
            NumericVectorM&, // RHS vector
            const double,      // Stopping tolerance
            const  int,
            const int
                                          ) = 0; // N. Iteration

    /// This function calls the solver "_solver_type" preconditioned with
    /// the "_preconditioner_type" preconditioner.  The preconditioning
    /// matrix is used if it is provided, or the system matrix is used if
    /// precond_matrix is null
    std::pair<int, double> solve (SparseMatrixM& matrix,
                                  SparseMatrixM* precond_matrix,
                                  NumericVectorM&, // Solution vector
                                  NumericVectorM&, // RHS vector
                                  const double,      // Stopping tolerance
                                  const  int,
                                  const int
                                 ); // N. Iterations


    /// After calling this method, all successive solves will be
    /// restricted to the given set of dofs, which must contain local
    /// dofs on each processor only and not contain any duplicates.  This
    /// mode can be disabled by calling this method with \p dofs being a \p NULL pointer.
    virtual void restrict_solve_to (
        const std::vector< int>* const dofs,
        const SubsetSolveModeM subset_solve_mode=SUBSET_ZERO);

/// This function calls the solver
//  virtual void MGSolve(std::vector<SparseMatrixM *>  &matrix_in,   // System Matrix
// 	              std::vector<NumericVectorM *> &solution_in,// Solution vector
// 		      std::vector<NumericVectorM *> &rhs_in,     // RHS vector
// 		      Matrix *P, Matrix *R,
// 		      const double tol,
// 		      const  int m_its) =0;

    /// Prints a useful message about why the latest linear solve con(di)verged.
    virtual void print_converged_reason() = 0;
};


// =============================================
inline LinearSolverM::LinearSolverM(
    const ParallelM::Communicator &comm_in
) :
    ParallelObjectM(comm_in),
    _solver_type(GMRESM),                //     _solver_type(PREONLYM),
    _is_initialized(false),
    _preconditioner_type(ILU_PRECONDM), // _preconditioner_type(IDENTITY_PRECONDM),
    //     _vanka_smoother_type(NO_SMOOTHING_VANKA),
    _preconditioner(NULL),
    same_preconditioner(false)
    {}
// ========================================================
inline LinearSolverM::~LinearSolverM() {
    this->clear();
}

// ==============================
inline bool LinearSolverM::get_same_preconditioner() {
    return same_preconditioner;
}

// =====================================================
inline
std::pair<int, double>
LinearSolverM::solve (SparseMatrixM&   mat,
                      SparseMatrixM*   pc_mat,
                      NumericVectorM&  sol,
                      NumericVectorM&  rhs,
                      const double       tol,
                      const int n_iter,
                      const int clearing
                     ) {
    if (pc_mat)  return this->solve(mat, *pc_mat, sol, rhs, tol, n_iter,clearing);
    else   return this->solve(mat, sol, rhs, tol, n_iter,clearing);
}

#endif // #ifdef __solver_h__


