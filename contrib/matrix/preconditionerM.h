#ifndef __preconditionerM_h__
#define __preconditionerM_h__

// Local includes
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"
#include "SolverPackage_enum.h"
#include "Precondtype_enum.h"
#include "parallel_objectM.h"

// C++ includes
#include <memory>
#include <iostream>


// forward declarations
class SparseMatrixM;
class NumericVectorM;


// ========================================================================
/// This class provides a uniform interface for preconditioners.  This base
/// class is overloaded to provide linear solvers from different packages
///
/// Preconditioner Data:
///   SparseMatrixM * _matrix;                   P=_matrix
///   PreconditionerTypeM _preconditioner_type;  type
///   bool _is_initialized;                      flag (if P)
///
///  P Constructor is empty:
///    ParallelM::Communicator &comm      parallel communicator
///
///  P init:
///    not defined
///
///  build by using the linear solver package -> P
///     ParallelM::Communicator &comm ,           parallel communicator
///     SolverPackageM solver_package = LSOLVER   linear solver
///
///  P by user:
///        set_type   type of preconditioner
///        setup      your preconditioning matrix   
///
/// In the below comments P is the matrix to be preconditioned with Apply()
/// performing the equivalent of the matrix vector product P^-1 x.  This
/// can also be thought of as (usually approximately) solving for Py=x.
// ========================================================================
class PreconditionerM:
#ifdef LM_REFCOUNT
    public ReferenceCountedObject<PreconditionerM>,
#endif
    public ParallelObjectM {
protected:

    /// The matrix P... ie the matrix to be preconditioned.
    ///  This is often the actual system matrix of a linear sytem.
    SparseMatrixM * _matrix;

    /// Enum statitng with type of preconditioner to use.
    PreconditionerTypeM _preconditioner_type;

    /// Flag indicating if the data structures have been initialized.
    bool _is_initialized;

public:
    //  Constructor-Destructor -------------------------------------
    ///  Constructor
    PreconditionerM (
        const ParallelM::Communicator &comm  ///< parallel communicator
    );
    /// Destructor
    virtual ~PreconditionerM (){this->clear ();}

    /// Builds a \p PreconditionerM using the linear solver package
    /// specified by \p solver_package
    static PreconditionerM * build(
        const ParallelM::Communicator &comm CAN_DEFAULT_TO_COMMWORLD, ///< parallel communicator <-
        const SolverPackageM solver_package = LSOLVER  ///< linear solver <-
    );

    /// Release all memory and clear data structures.
    virtual void clear () {}

    /// Initialize data structures if not done so already.
    virtual void init () {};

    //  Set functions  -----------------------------

    /// Sets the matrix P to be preconditioned.
    void set_matrix(
      SparseMatrixM & mat              ///< sparse matrix <-
    ) { //If the matrix is changing then we (probably) need to reinitialize.
    _is_initialized = false; _matrix = &mat;
}
    /// Sets the type of preconditioner to use.
    void set_type (
      const PreconditionerTypeM pct  ///< preconditioner type <-
    ) {//If the preconditioner type changes we (probably) need to reinitialize.
    _is_initialized = false;    _preconditioner_type = pct;
}
/// If you need to fill in your preconditioning matrix.
    virtual void setup () {}

    //  Return functions  -----------------------------
    /// @returns true if the data structures are initialized
    bool initialized () const { return _is_initialized; }
    /// @returns th preconditioner type
    PreconditionerTypeM type () const {return _preconditioner_type;}

//  Operation functions  -----------------------------

    /// Computes the preconditioned vector "y" based on input "x".
    /// Usually by solving Py=x to get the action of P^-1 x.
    virtual void apply(
        const NumericVectorM & x,     ///< rhs  <-
        NumericVectorM & y            ///< solution <-
    ) = 0;

};




/*----------------------- inline functions ----------------------------------*/
// =============================================
inline PreconditionerM::PreconditionerM (
    const ParallelM::Communicator &comm_in  ///< communicator <-
) :
    ParallelObjectM(comm_in),
    _matrix(NULL),
    _preconditioner_type (ILU_PRECONDM),
    _is_initialized      (false)
{
}

// // ========================================================
// inline void PreconditionerM::set_matrix(
//   SparseMatrixM & mat           ///< matrix <-
// ) {
//     //If the matrix is changing then we (probably) need to reinitialize.
//     _is_initialized = false;
//     _matrix = &mat;
// }
// // ==========================================================
// inline void PreconditionerM::set_type (const PreconditionerTypeM pct) {
//     //If the preconditioner type changes we (probably) need to reinitialize.
//     _is_initialized = false;
//     _preconditioner_type = pct;
// }


#endif // #ifdef __preconditioner_h__
