#include "Solverlib_conf.h"

#ifdef HAVE_PETSCM

// C++ includes

// Local Includes
// #include "auto_ptr.h"
#include "petsc_preconditionerM.h"
#include "petsc_macroM.h"
#include "petsc_matrixM.h"
#include "petsc_vectorM.h"

#include <mpi.h>  //for MPI_COMM_WORLD

// #include "libmesh_common.h"

// =============================================================================
void PetscPreconditionerM::apply(const NumericVectorM & x, NumericVectorM & y){
  PetscVectorM & x_pvec = dynamic_cast<PetscVectorM&>(const_cast<NumericVectorM&>(x));
  PetscVectorM & y_pvec = dynamic_cast<PetscVectorM&>(const_cast<NumericVectorM&>(y));
  Vec x_vec = x_pvec.vec();  Vec y_vec = y_pvec.vec();
  PCApply(_pc,x_vec,y_vec);
}
// =================================================
void PetscPreconditionerM::init () {
  if(!this->_matrix) {
    std::cerr << "ERROR: No matrix set for PetscPreconditioner, but init() called" << std::endl;
    abort();
  }
  //Clear the preconditioner in case it has been created in the past
  if(!this->_is_initialized)  {
     // Should probably use PCReset(), but it's not working at the moment so we'll destroy instead
    if (_pc){ int ierr = PCDestroy(&_pc); CHKERRABORT(MPI_COMM_WORLD,ierr); }
 
    //Create the preconditioning object
    int ierr=PCCreate(this->comm().get(),&_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    //Set the PCType
    set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
// #ifdef LIBMESH_HAVE_PETSC_HYPRE
//     if(this->_preconditioner_type == AMG_PRECOND)
//       PCHYPRESetType(this->_pc, "boomerang");
// #endif
    PetscMatrixM * pmatrix = libmeshM_cast_ptr<PetscMatrixM*, SparseMatrixM >(this->_matrix);
    _mat = pmatrix->mat();
  }
  int ierr=PCSetOperators(_pc,_mat,_mat
#ifdef OLD_PETSC
  ,SAME_NONZERO_PATTERN
#endif
  ); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set the PCType.  Note: this used to be done *before* the call to
  // PCSetOperators(), and only when !_is_initialized, but
  // 1.) Some preconditioners (those employing sub-preconditioners,
  // for example) have to call PCSetUp(), and can only do this after
  // the operators have been set.
  // 2.) It should be safe to call set_petsc_preconditioner_type()
  // multiple times.
  set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
  
  
  this->_is_initialized = true;
}

// =====================================================
void PetscPreconditionerM::set_petsc_preconditioner_type
                (const PreconditionerTypeM & preconditioner_type, PC & pc){
  int ierr = 0;
   // get the communicator from the PETSc object
  ParallelM::communicator comm;
  PetscObjectGetComm((PetscObject)pc, &comm);
  ParallelM::Communicator communicator(comm);
  
  switch (preconditioner_type)  {
    
  case IDENTITY_PRECONDM:
    ierr = PCSetType (pc, (char*) PCNONE);      CHKERRABORT(MPI_COMM_WORLD,ierr); break;
	
  case CHOLESKY_PRECONDM:
    ierr = PCSetType (pc, (char*) PCCHOLESKY);  CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case ICC_PRECONDM:
    ierr = PCSetType (pc, (char*) PCICC);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;
  
  case ILU_PRECONDM:
    {
//       int nprocs; MPI_Comm_size(MPI_COMM_WORLD,&nprocs);  //TODO
      // In serial, just set the ILU preconditioner type
       if (communicator.size())
         {
          ierr = PCSetType (pc, (char*) PCILU);
          CHKERRABORT(MPI_COMM_WORLD,ierr);
         }
       else
         {
          // But PETSc has no truly parallel ILU, instead you have to set
          // an actual parallel preconditioner (e.g. block Jacobi) and then
          // assign ILU sub-preconditioners.
          ierr = PCSetType (pc, (char*) PCBJACOBI);
          CHKERRABORT(MPI_COMM_WORLD,ierr);

          // Set ILU as the sub preconditioner type
          set_petsc_subpreconditioner_type(PCILU, pc);
         }
      break;
    }

  case LU_PRECONDM:
    {
      
//       int nprocs; MPI_Comm_size(MPI_COMM_WORLD,&nprocs);  //TODO
      // In serial, just set the LU preconditioner type
         if (communicator.size())
         {
          ierr = PCSetType (pc, (char*) PCLU);
          CHKERRABORT(MPI_COMM_WORLD,ierr);
         }
       else
         {
          // But PETSc has no truly parallel LU, instead you have to set
          // an actual parallel preconditioner (e.g. block Jacobi) and then
          // assign LU sub-preconditioners.
          ierr = PCSetType (pc, (char*) PCBJACOBI);
          CHKERRABORT(MPI_COMM_WORLD,ierr);

          // Set LU as the sub preconditioner type
          set_petsc_subpreconditioner_type(PCLU, pc);
         }
      break;
    }
    
//   case SLU_PRECONDM:
//    ierr = PCSetType (pc, (char*) PCLU);   CHKERRABORT(MPI_COMM_WORLD,ierr); 
//    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);  CHKERRABORT(MPI_COMM_WORLD,ierr);   break;    //here we set the SuperLU_dist solver package
// 
// 
//   case MLU_PRECONDM:
//    ierr = PCSetType (pc, (char*) PCLU);                    CHKERRABORT(MPI_COMM_WORLD,ierr); 
//    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);  CHKERRABORT(MPI_COMM_WORLD,ierr);                    //here we set the MUMPS parallel direct solver package
//    break;
//    
//    
//   case MCC_PRECONDM:
//    ierr = PCSetType (pc, (char*) PCCHOLESKY);                    CHKERRABORT(MPI_COMM_WORLD,ierr); 
//    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);  CHKERRABORT(MPI_COMM_WORLD,ierr);                    //here we set the MUMPS parallel direct solver package
//    break; 
   
   
   
//   case ILU_PRECONDM:
//     ierr = PCSetType (pc, (char*) PCILU);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;
// 
//   case LU_PRECONDM:
//     ierr = PCSetType (pc, (char*) PCLU);        CHKERRABORT(MPI_COMM_WORLD,ierr); break;
      
  case ASM_PRECONDM:
    ierr = PCSetType (pc, (char*) PCASM);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case JACOBI_PRECONDM:
    ierr = PCSetType (pc, (char*) PCJACOBI);    CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case BLOCK_JACOBI_PRECONDM:
    ierr = PCSetType (pc, (char*) PCBJACOBI);   CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case SOR_PRECONDM:
    ierr = PCSetType (pc, (char*) PCSOR);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case EISENSTAT_PRECONDM:
    ierr = PCSetType (pc, (char*) PCEISENSTAT); CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case AMG_PRECONDM:
    ierr = PCSetType (pc, (char*) PCHYPRE);     CHKERRABORT(MPI_COMM_WORLD,ierr); break;
    
/*  case MG_PRECONDM:
    ierr = PCSetType (pc, (char*) PCMG);        CHKERRABORT(MPI_COMM_WORLD,ierr); break;*/ 

#if !(PETSC_VERSION_LESS_THAN(2,1,2))
    // Only available for PETSC >= 2.1.2      
  case USER_PRECONDM:
    ierr = PCSetType (pc, (char*) PCMAT);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;
#endif

  case SHELL_PRECONDM:
    ierr = PCSetType (pc, (char*) PCSHELL);     CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  default:
    std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
              << preconditioner_type << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }

  //Let the commandline override stuff
  if( preconditioner_type != AMG_PRECONDM /*&& preconditioner_type != MG_PRECONDM*/)   PCSetFromOptions(pc);  //!!!!!!
}


// template <typename T>
// #if PETSC_VERSION_LESS_THAN(3,0,0)
//  void PetscPreconditioner<T>::set_petsc_subpreconditioner_type(PCType type, PC& pc)
// #else
 void PetscPreconditionerM::set_petsc_subpreconditioner_type(const PCType type, PC& pc)
// #endif
{
  // For catching PETSc error return codes
  int ierr = 0;
  // get the communicator from the PETSc object
  ParallelM::communicator comm;
  PetscObjectGetComm((PetscObject)pc, &comm);
  ParallelM::Communicator communicator(comm);
  // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
  // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
  //
  // "Object is in wrong state!"
  // "Matrix must be set first."
  //
  // error messages...
  ierr = PCSetUp(pc);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // To store array of local KSP contexts on this processor
   KSP* subksps;
// 
//   // the number of blocks on this processor
   int n_local;
// 
//   // The global number of the first block on this processor.
//   // This is not used, so we just pass PETSC_NULL instead.
//   // int first_local;
// 
//   // Fill array of local KSP contexts
   ierr = PCBJacobiGetSubKSP(pc, &n_local, PETSC_NULL, &subksps);
   CHKERRABORT(comm,ierr);
// 
//   // Loop over sub-ksp objects, set ILU preconditioner
  for (int i=0; i<n_local; ++i)
    {
      // Get pointer to sub KSP object's PC
      PC subpc;
      ierr = KSPGetPC(subksps[i], &subpc);
      CHKERRABORT(comm,ierr);

      // Set requested type on the sub PC
      ierr = PCSetType(subpc, type);
      CHKERRABORT(comm,ierr);
    }
}


//------------------------------------------------------------------

#endif // #ifdef LIBMESH_HAVE_PETSC
