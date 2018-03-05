#ifndef __petsc_preconditionerM_h__
#define __petsc_preconditionerM_h__

#include "Solverlib_conf.h"


#ifdef HAVE_PETSCM

#include "Precondtype_enum.h"

// Local includes
#include "preconditionerM.h"
// #include "libmesh_common.h"
// // // #include "enum_solver_package.h"   //TODO
// // // #include "enum_preconditioner_type.h"  //TODO
#include "parallel_objectM.h"

// Petsc includes
#include "petscpc.h"

// forward declarations
class SparseMatrixM;
class NumericVectorM;

/// This class provides an interface to  preconditioners  from Petsc.
class PetscPreconditionerM : public PreconditionerM{
  
  protected:
  PC _pc; ///<Preconditioner context
  Mat _mat;///< Petsc Matrix 
  
public:
  // Constructor
  ///  Constructor. Initializes PetscPreconditioner data structures 
  PetscPreconditionerM (
    const ParallelM::Communicator &comm
  ); 
  /// Destructor.
  virtual ~PetscPreconditionerM (){this->clear ();  }  
  /// Release all memory and clear data structures. 
  virtual void clear () {}
  /// Initialize data structures if not done so already.
  virtual void init ();
  // Return
  /// Returns the actual Petsc PC struct.     
  PC pc() { return _pc; }
  
  // Compute 
   /// Computes the preconditioned vector "y" based on input "x".
  /// Usually by solving Py=x to get the action of P^-1 x.
  virtual void apply(const NumericVectorM & x, NumericVectorM & y);
  /// Tells PETSC to use the user-specified preconditioner  
  static void set_petsc_preconditioner_type 
             (const PreconditionerTypeM & preconditioner_type, PC & pc);
private:
  /**
   * Some PETSc preconditioners (ILU, LU) don't work in parallel.  This function
   * is called from set_petsc_preconditioner_type() to set additional options
   * for those so-called sub-preconditioners.  This method ends up being static
   * so that it can be called from set_petsc_preconditioner_type().  Not sure
   * why set_petsc_preconditioner_type() needs to be static though...
   */
// #if PETSC_VERSION_LESS_THAN(3,0,0)
// //  In Petsc 2.3.3, PCType was #define'd as const char*
//   static void set_petsc_subpreconditioner_type(PCType type, PC& pc);
// #else
//  In later versions, PCType is #define'd as char*, so we need the const
  static void set_petsc_subpreconditioner_type(const PCType type, PC& pc);
// #endif	     
	     
};




/*----------------------- inline functions ----------------------------------*/
/*
inline PetscPreconditionerM::PetscPreconditionerM () :  PreconditionerM(){}
inline PetscPreconditionerM::~PetscPreconditionerM ()*/
// ===============================================
inline
PetscPreconditionerM::PetscPreconditionerM (
  const ParallelM::Communicator &comm) :
  PreconditionerM(comm),
  _pc(PETSC_NULL)
{}
#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_linear_solver_h__
