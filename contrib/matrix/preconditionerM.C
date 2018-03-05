
// Standard include
#include <cstdlib>  //for abort

// Local config
#include "Solverlib_conf.h"

// Class Includes
#include "preconditionerM.h"

// Local Includes
#include "petsc_preconditionerM.h"
//here, only the Petsc Preconditioner is settled here, why not also the Laspack one?

//------------------------------------------------------------------
// Preconditioner members

PreconditionerM *PreconditionerM::build(
  const ParallelM::Communicator &comm,   ///< parallel communicator
  const SolverPackageM solver_package    ///< solver package
){
  // Build the appropriate solver
  switch (solver_package){

#ifdef _HAVE_LASPACK
    case LASPACK_SOLVERS:
      {AutoPtr<PreconditionerM > ap(new LaspackPreconditionerM(comm));return ap;}
#endif
#ifdef HAVE_PETSCM
    case PETSC_SOLVERSM:
      {return new PetscPreconditionerM(comm);}
#endif
#ifdef LIBMESH_HAVE_TRILINOS
    case TRILINOS_SOLVERS:
    {AutoPtr<PreconditionerM > ap(new AztecPreconditionerM(comm));return ap;}
#endif
    default:
      std::cerr << "ERROR:  Unrecognized solver package: "
		<< solver_package
		<< std::endl;
      abort();
    }
  return NULL;    
}

