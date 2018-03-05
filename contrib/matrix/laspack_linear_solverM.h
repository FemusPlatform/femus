#ifndef __laspack_linear_solverM_h__
#define __laspack_linear_solverM_h__

#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM
//#if defined(LIBMESH_HAVE_LASPACK) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
// -----------------------------------
// C++ includes 
// -----------------------------------
// Local includes 
#include "linear_solverM.h"
#include "laspack_vectorM.h"
#include "laspack_matrixM.h"
#include <itersolv.h>         // laspack 
#include <rtc.h>              // laspack 
#include <errhandl.h>         // laspack 
// -----------------------------------

// ===============================================
// This class provides an interface to Laspack
// iterative solvers that is compatible with the \p libMesh
// \p LinearSolver<>
// =============================================

class LaspackLinearSolverM : public LinearSolverM
{
  private:
  /// Preconditioner type
  PrecondProcType _precond_type;
  
 public:
  ///  Constructor. Initializes Laspack data structures
  LaspackLinearSolverM ();
  /// Destructor.
  ~LaspackLinearSolverM ();
  /// Release all memory and clear data structures.
  void clear ();

  /// Initialize data structures if not done so already.
  void init ();
  
// ===================================
// SOLVE
// ===================================
  /// Call the Laspack solver    
  std::pair<unsigned int, Real> 
    solve (SparseMatrixM  &matrix,NumericVectorM &solution,NumericVectorM &rhs,
	   const double tol, const unsigned int m_its);
  /// Call the Laspack solver   
  std::pair<unsigned int, Real> 
    solve (SparseMatrixM  &matrix,SparseMatrixM  &pc,NumericVectorM &solution,NumericVectorM &rhs,
	   const double tol, const unsigned int m_its);
	   
  void MGStep(int NoLevels, int Level, int Gamma);   
	   
  /// This function calls the solver	       
  void MGSolve(std::vector<SparseMatrixM*>  &matrix_in, std::vector<NumericVectorM*> &solution_in,
		      std::vector<NumericVectorM*> &rhs_in,  Matrix *P, Matrix *R,   
		      const double tol,  const unsigned int m_its);
  /// Prints a useful message about why the latest linear solve con(di)verged.
  virtual void print_converged_reason();
  
 private:
  ///  * Tells LASPACK to use the user-specified preconditioner stored in \p _preconditioner_type
  void set_laspack_preconditioner_type ();

};


/*----------------------- functions ----------------------------------*/

inline LaspackLinearSolverM::LaspackLinearSolverM () : _precond_type (ILUPrecond){}
// ==================================================================
inline LaspackLinearSolverM::~LaspackLinearSolverM (){  this->clear ();}
// ====================================================================
inline std::pair<unsigned int, Real> LaspackLinearSolverM::solve 
                            (SparseMatrixM&, SparseMatrixM&,NumericVectorM&,NumericVectorM&,
			     const double,const unsigned int){
  std::cerr << "ERROR: LASPACK does not support a user-supplied preconditioner!" << std::endl; abort();
  std::pair<unsigned int, Real> p;
  return p;
}

#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifndef __laspack_linear_solver_h__
