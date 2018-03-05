#include "Solverlib_conf.h"

#ifdef HAVE_PETSCM

// Local Includes
#include "petsc_macroM.h"
// #include "libmesh_logging.h"
#include "petsc_linear_solverM.h"
// #include "shell_matrix.h"
#include "petsc_preconditionerM.h"
#include "petsc_vectorM.h"
#include "petsc_matrixM.h"
#include "petsc_MmatrixM.h"

// C++ includes
#include <string.h>

// ==========================================================
extern "C" {
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  // ------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(void * ctx) {
    PreconditionerM * preconditioner = static_cast<PreconditionerM*>(ctx);
       if(!preconditioner->initialized()){
      std::cerr<<"Preconditioner not initialized!  Make sure you call init() before solve!"<<std::endl;
      exit(0);
    }
    preconditioner->setup();
    
//     preconditioner->init();
    return 0;
  }
// ------------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)  {
    PreconditionerM * preconditioner = static_cast<PreconditionerM*>(ctx);
     PetscVectorM x_vec(x, preconditioner->comm());
    PetscVectorM y_vec(y, preconditioner->comm());
//     PetscVectorM x_vec(x);
//     PetscVectorM y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#else
// ----------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(PC pc) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    PreconditionerM * preconditioner = static_cast<PreconditionerM*>(ctx);
    if(!preconditioner->initialized()){
      std::cerr<<"Preconditioner not initialized!  Make sure you call init() before solve!"<<std::endl;
      exit(0);
    }
    preconditioner->setup();
//     preconditioner->init();
    return 0;
  }
// --------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
    PreconditionerM * preconditioner = static_cast<PreconditionerM*>(ctx);
    PetscVectorM x_vec(x, preconditioner->comm());
    PetscVectorM y_vec(y, preconditioner->comm());
//     PetscVectorM x_vec(x); PetscVectorM y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#endif
} // end extern "C
// ================================================


// /*extern "C"
// {
// #if PETSC_VERSION_LESS_THAN(2,2,1)
//   typedef int PetscErrorCode;
//   typedef int PetscInt;
// #endif
// 
// 
// #if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
//   PetscErrorCode __libmesh_petsc_preconditioner_setup (void * ctx)
//   {
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
//     preconditioner->init();
// 
//     return 0;
//   }
//   
// 
//   PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)
//   {
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
// 
//     PetscVector<Number> x_vec(x);
//     PetscVector<Number> y_vec(y);
// 
//     preconditioner->apply(x_vec,y_vec);
// 
//     return 0;
//   }
// #else
//   PetscErrorCode __libmesh_petsc_preconditioner_setup (PC pc)
//   {
//     void *ctx;
//     PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
//     preconditioner->init();
// 
//     return 0;
//   }
// 
//   PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y)
//   {
//     void *ctx;
//     PetscErrorCode ierr = PCShellGetContext(pc,&ctx);CHKERRQ(ierr);
//     Preconditioner<Number> * preconditioner = static_cast<Preconditioner<Number>*>(ctx);
// 
//     PetscVector<Number> x_vec(x);
//     PetscVector<Number> y_vec(y);
// 
//     preconditioner->apply(x_vec,y_vec);
// 
//     return 0;
//   }
// #endif
// } // end extern "C"*/

// ==============================================
// ----------------------- functions ------
// ==============================================

void PetscLinearSolverM::clear() {
  if (this->initialized())    {
      /* If we were restricted to some subset, this restriction must
	 be removed and the subset index set destroyed.  */
      if(_restrict_solve_to_is!=NULL)	{
	  PetscErrorCode ierr = ISDestroy(&_restrict_solve_to_is);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  _restrict_solve_to_is = NULL;
	}

      if(_restrict_solve_to_is_complement!=NULL){
	  PetscErrorCode ierr = ISDestroy(&_restrict_solve_to_is_complement);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  _restrict_solve_to_is_complement = NULL;
	}
      this->_is_initialized = false;
      PetscErrorCode ierr=0;
#if PETSC_VERSION_LESS_THAN(2,2,0)
  // 2.1.x & earlier style
  ierr = SLESDestroy(_sles);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
  // 2.2.0 & newer style
  ierr = LibMeshKSPDestroy(&_ksp);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
      // Mimic PETSc default solver and preconditioner
      this->_solver_type           = GMRESM;
      if(!this->_preconditioner)      {
        if (this->n_processors() == 1) this->_preconditioner_type = ILU_PRECONDM;
        else   this->_preconditioner_type = BLOCK_JACOBI_PRECONDM;
      }
    }
}


// ==============================================================
void PetscLinearSolverM::init() {
  // Initialize the data structures if not done so already.
  if (!this->initialized()) {
    this->_is_initialized = true;  int ierr=0;

    // Create the linear solver context
    ierr = KSPCreate(this->comm().get(), &_ksp);   CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
//     ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    // NOT NECESSARY!!!!
    //ierr = PCSetFromOptions (_pc);
    //CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
     KSPType ksp_type;
     ierr = KSPGetType (_ksp, &ksp_type);   CHKERRABORT(MPI_COMM_WORLD,ierr);

      if (strcmp(ksp_type, "preonly"))        {
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE); CHKERRABORT(MPI_COMM_WORLD,ierr);
        }

    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    PetscPreconditionerM::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

    //If there is a preconditioner object we need to set the internal setup and apply routines
    if (this->_preconditioner) {
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}

// ========================================================
void PetscLinearSolverM::init(PetscMatrixM* matrix) {
  // Initialize the data structures if not done so already.
  if (!this->initialized())    {
    this->_is_initialized = true;   int ierr=0;
// #if PETSC_VERSION_LESS_THAN(2,2,0)  // 2.1.x & earlier style
//     // Create the linear solver context
//     ierr = SLESCreate(MPI_COMM_WORLD, &_sles);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Create the Krylov subspace & preconditioner contexts
//     ierr = SLESGetKSP(_sles, &_ksp);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = SLESGetPC(_sles, &_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Have the Krylov subspace method use our good initial guess rather than 0
//     ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     // Set user-specified  solver and preconditioner types
//     this->set_petsc_solver_type();
//     // Set the options from user-input
//     // Set runtime options, e.g.,
//     //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//     //  These options will override those specified above as long as
//     //  SLESSetFromOptions() is called _after_ any other customization
//     //  routines.
//     ierr = SLESSetFromOptions(_sles);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//
// #else // 2.2.0 & newer style
    // Create the linear solver context
    ierr = KSPCreate(this->comm().get(), &_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    //ierr = PCCreate (MPI_COMM_WORLD, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set operators. The input matrix works as the preconditioning matrix
    ierr = KSPSetOperators(_ksp, matrix->mat(), matrix->mat()
#ifdef OLD_PETSC
    , SAME_NONZERO_PATTERN
#endif      
    );
    CHKERRABORT(MPI_COMM_WORLD,ierr);  
    
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  KSPSetFromOptions() is called _after_ any other customization
      //  routines.
    ierr = KSPSetFromOptions (_ksp);CHKERRABORT(MPI_COMM_WORLD,ierr);  
    
    
     KSPType ksp_type;
     ierr = KSPGetType (_ksp, &ksp_type);CHKERRABORT(MPI_COMM_WORLD,ierr);  

      if (strcmp(ksp_type, "preonly")){
          ierr = KSPSetInitialGuessNonzero (_ksp, PETSC_TRUE);CHKERRABORT(MPI_COMM_WORLD,ierr);  
        }

    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
//     ierr = KSPSetFromOptions(_ksp);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    //ierr = PCSetFromOptions (_pc);CHKERRABORT(MPI_COMM_WORLD,ierr);

// #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    



    PetscPreconditionerM::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
    if (this->_preconditioner) {
      this->_preconditioner->set_matrix(*matrix);
      this->_preconditioner->init();
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}



// ============================================
void PetscLinearSolverM::restrict_solve_to (
  const std::vector<int>* const dofs,					 
  const SubsetSolveModeM subset_solve_mode){
  PetscErrorCode ierr=0;

  /* The preconditioner (in particular if a default preconditioner)
     will have to be reset.  We call this->clear() to do that.  This
     call will also remove and free any previous subset that may have
     been set before.  */
  this->clear();
  _subset_solve_mode = subset_solve_mode;

  if(dofs!=NULL)    {
      PetscInt* petsc_dofs = NULL;
      ierr = PetscMalloc(dofs->size()*sizeof(PetscInt), &petsc_dofs);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      for(size_t i=0; i<dofs->size(); i++){ petsc_dofs[i] = (*dofs)[i];}

      ierr = ISCreateGeneral(this->comm().get(),dofs->size(),petsc_dofs,PETSC_OWN_POINTER,&_restrict_solve_to_is);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
}

// ========================================================
std::pair< int, double> PetscLinearSolverM::solve(
  SparseMatrixM&  matrix_in,
  SparseMatrixM&  precond_in,  
  NumericVectorM& solution_in,  
  NumericVectorM& rhs_in,
    const double tol,   
    const  int m_its,
    const int clearing
                                                 ) {

//   START_LOG("solve()", "PetscLinearSolverM");
  // Make sure the data passed in are really of Petsc types
  PetscMatrixM* matrix   = dynamic_cast<PetscMatrixM*>(&matrix_in);
  PetscMatrixM* precond  = dynamic_cast<PetscMatrixM*>(&precond_in);
  PetscVectorM* solution = dynamic_cast<PetscVectorM*>(&solution_in);
  PetscVectorM* rhs      = dynamic_cast<PetscVectorM*>(&rhs_in);
  
//   Mat precond_ = precond->mat();
 if(clearing==1) {    this->clear(); std::cout << " Clearing the precond "<<std::endl;
} 
 this->init(matrix);  //as before
  
//   this->init();

  int ierr=0;  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;
  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();  precond->close();  solution->close(); rhs->close();

  
  Mat submat = NULL;  Mat subprecond = NULL;
  Vec subrhs = NULL;  Vec subsolution = NULL;
  VecScatter scatter = NULL;
  PetscMatrixM* subprecond_matrix = NULL;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if(_restrict_solve_to_is!=NULL)   {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      
      ierr = MatGetSubMatrix(matrix->mat(), _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
  
      
      
        /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO){
	  _create_complement_is(rhs_in);
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;	  Mat submat1 = NULL;	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0); CHKERRABORT(MPI_COMM_WORLD,ierr);
  
	    ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScatterDestroy(&scatter1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecDestroy(&subvec1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = MatDestroy(&submat1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond
#ifdef OLD_PETSC
			     ,this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN
#endif
			    );
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner){
	  subprecond_matrix = new PetscMatrixM(subprecond,this->comm());
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else    {

      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat()
#ifdef OLD_PETSC
			     ,this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN
#endif
	
      );
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner)      {
	this->_preconditioner->set_matrix(matrix_in);
        this->_preconditioner->init();
      }
    }

  
  
  
  
  
  
  
  //   // If matrix != precond, then this means we have specified a
//   // special preconditioner, so reset preconditioner type to PCMAT.
//   if (matrix != precond)
//     {
//       this->_preconditioner_type = USER_PRECOND;
//       this->set_petsc_preconditioner_type ();
//     }
  
//   if (this->_preconditioner) { 
//     this->_preconditioner->set_matrix(matrix_in); 
//   }
  
  
  // 2.2.1 & newer style
  // Set operators. The input matrix works as the preconditioning matrix
  
  
//     if( this->_preconditioner_type==MCC_PRECONDM || this->_preconditioner_type==ICC_PRECONDM ) {
// //       /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
//       ierr = MatSetOption(matrix->mat(),MAT_SYMMETRIC,PETSC_TRUE); CHKERRABORT(MPI_COMM_WORLD,ierr); 
//       ierr = MatSetOption(matrix->mat(),MAT_SPD,PETSC_TRUE);       CHKERRABORT(MPI_COMM_WORLD,ierr); /* set MUMPS id%SYM=1 */
//      }
  
//   if( this->_preconditioner_type==MLU_PRECONDM) {
//       PetscInt ival,icntl;
//       
//       ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, DIFFERENT_NONZERO_PATTERN);  CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       ierr = PCFactorSetUpMatSolverPackage(_pc);           CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//       ierr = PCFactorGetMatrix(_pc,&precond_);             CHKERRABORT(MPI_COMM_WORLD,ierr);
//       icntl=7; ival = 2;
//       ierr = MatMumpsSetIcntl(precond_,icntl,ival);        CHKERRABORT(MPI_COMM_WORLD,ierr);
//   } 
  
//   else if( this->_preconditioner_type==MCC_PRECONDM) {
//       PetscInt ival,icntl;
//       
//       ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, DIFFERENT_NONZERO_PATTERN);  CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//       ierr = PCFactorSetUpMatSolverPackage(_pc);           CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//       ierr = PCFactorGetMatrix(_pc,&precond_);             CHKERRABORT(MPI_COMM_WORLD,ierr);
//       icntl=7; ival = 2;
//       ierr = MatMumpsSetIcntl(precond_,icntl,ival);        CHKERRABORT(MPI_COMM_WORLD,ierr);
//       
//   } else {
  
//     if (!this->same_preconditioner)  {
//       ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, SAME_NONZERO_PATTERN);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
//      } else {
//        ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, SAME_PRECONDITIONER);
//        CHKERRABORT(MPI_COMM_WORLD,ierr);
//      }
//   }
  
  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT,PETSC_DEFAULT, max_its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  
  // Solve the linear system
  
//        PetscLogEvent USER_EVENT;
//      PetscLogDouble user_event_flops;
//      PetscLogEventRegister("User event",0,&USER_EVENT);
//      PetscLogEventBegin(USER_EVENT,0,0,0,0);
       
   // Solve the linear system
  if(_restrict_solve_to_is!=NULL){
      ierr = KSPSolve (_ksp, subrhs, subsolution);  CHKERRABORT(MPI_COMM_WORLD,ierr);  }
  else {ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());  CHKERRABORT(MPI_COMM_WORLD,ierr);}
//   ierr = KSPSolve(_ksp, rhs->vec(), solution->vec()); CHKERRABORT(MPI_COMM_WORLD,ierr);
//         PetscLogFlops(user_event_flops);
//      PetscLogEventEnd(USER_EVENT,0,0,0,0);

     // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber(_ksp, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm(_ksp, &final_resid); CHKERRABORT(MPI_COMM_WORLD,ierr);

   if(_restrict_solve_to_is!=NULL)    {
      switch(_subset_solve_mode)	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec()); CHKERRABORT(MPI_COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterDestroy(&scatter);      CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner)	{
	  /* Before we delete subprecond_matrix, we should give the
	     _preconditioner a different matrix.  */
	  this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
	  delete subprecond_matrix;
	  subprecond_matrix = NULL;
	}

      ierr = VecDestroy(&subsolution);    CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&subrhs);         CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&submat);         CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&subprecond);     CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
//   STOP_LOG("solve()", "PetscLinearSolverM");
  return std::make_pair(its, final_resid);
}

// ============================================
std::pair<int, double> PetscLinearSolverM::adjoint_solve (
  SparseMatrixM&  matrix_in,
				     NumericVectorM& solution_in,
				     NumericVectorM& rhs_in,
				     const double tol,
				     const int m_its){
//   START_LOG("solve()", "PetscLinearSolver");

  // Make sure the data passed in are really of Petsc types
  PetscMatrixM* matrix   = dynamic_cast<PetscMatrixM*>(&matrix_in);
  // Note that the matrix and precond matrix are the same
  PetscMatrixM* precond  = dynamic_cast<PetscMatrixM*>(&matrix_in);
  PetscVectorM* solution = dynamic_cast<PetscVectorM*>(&solution_in);
  PetscVectorM* rhs      = dynamic_cast<PetscVectorM*>(&rhs_in);

  this->init (matrix);

  PetscErrorCode ierr=0;
  PetscInt its=0, max_its = static_cast<PetscInt>(m_its);
  PetscReal final_resid=0.;

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close ();
  precond->close ();
  solution->close ();
  rhs->close ();

//   // If matrix != precond, then this means we have specified a
//   // special preconditioner, so reset preconditioner type to PCMAT.
//   if (matrix != precond)
//     {
//       this->_preconditioner_type = USER_PRECOND;
//       this->set_petsc_preconditioner_type ();
//     }

// 2.1.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)

  if(_restrict_solve_to_is!=NULL)    {    std::cout <<"    libmesh_not_implemented";      }

  // Based on http://wolfgang.math.tamu.edu/svn/public/deal.II/branches/MATH676/2008/deal.II/lac/source/petsc_solver.cc, http://tccc.iesl.forth.gr/AMS_EPEAEK/Elements/doc/in_html/petsc/SLES/index.html

  SLES sles;
  ierr = SLESCreate (this->comm().get(), &sles);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = SLESSetOperators (sles, matrix->mat(), precond->mat(), this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  KSP ksp;
  ierr = SLESGetKSP (sles, &ksp);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = SLESSetUp (sles, rhs->vec(), solution->vec());
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // See http://tccc.iesl.forth.gr/AMS_EPEAEK/Elements/doc/in_html/petsc/KSP/KSPSolveTrans.html#KSPSolveTrans
  ierr = SLESSolveTrans (ksp, &its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

// 2.2.0
#elif PETSC_VERSION_LESS_THAN(2,2,1)

  if(_restrict_solve_to_is!=NULL)    { std::cout <<"    libmesh_not_implemented";   }

  // Set operators. The input matrix works as the preconditioning matrix
  // This was commented earlier but it looks like KSPSetOperators is supported
  // after PETSc 2.2.0
  ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat(),
			 this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(MPI_COMM_WORLD,ierr);


  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  // Convergence is detected at iteration k if
  // ||r_k||_2 < max(rtol*||b||_2 , abstol)
  // where r_k is the residual vector and b is the right-hand side.  Note that
  // it is the *maximum* of the two values, the larger of which will almost
  // always be rtol*||b||_2.
  ierr = KSPSetTolerances (_ksp,
			   tol,           // rtol   = relative decrease in residual  (1.e-5)
			   PETSC_DEFAULT, // abstol = absolute convergence tolerance (1.e-50)
 			   PETSC_DEFAULT, // dtol   = divergence tolerance           (1.e+5)
			   max_its);
         CHKERRABORT(MPI_COMM_WORLD,ierr);


  // Set the solution vector to use
  ierr = KSPSetSolution (_ksp, solution->vec());
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Set the RHS vector to use
  ierr = KSPSetRhs (_ksp, rhs->vec());
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Solve the linear system
  ierr = KSPSolveTranspose (_ksp);
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);
         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);
         CHKERRABORT(MPI_COMM_WORLD,ierr);

// 2.2.1 & newer style
#else

  Mat submat = NULL;
  Mat subprecond = NULL;
  Vec subrhs = NULL;
  Vec subsolution = NULL;
  VecScatter scatter = NULL;
  PetscMatrixM* subprecond_matrix = NULL;

  // Set operators.  Also restrict rhs and solution vector to
  // subdomain if neccessary.
  if(_restrict_solve_to_is!=NULL)    {
      size_t is_local_size = this->_restrict_solve_to_is_local_size();

      ierr = VecCreate(this->comm().get(),&subrhs);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetSizes(subrhs,is_local_size,PETSC_DECIDE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subrhs);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecCreate(this->comm().get(),&subsolution);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetSizes(subsolution,is_local_size,PETSC_DECIDE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetFromOptions(subsolution);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is, subrhs,NULL, &scatter);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,rhs->vec(),subrhs,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      ierr = VecScatterBegin(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,solution->vec(),subsolution,INSERT_VALUES,SCATTER_FORWARD);
      CHKERRABORT(MPI_COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     PETSC_DECIDE,MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
      ierr = MatGetSubMatrix(matrix->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&submat);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatGetSubMatrix(precond->mat(),
			     _restrict_solve_to_is,_restrict_solve_to_is,
			     MAT_INITIAL_MATRIX,&subprecond);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif

      /* Since removing columns of the matrix changes the equation
	 system, we will now change the right hand side to compensate
	 for this.  Note that this is not necessary if \p SUBSET_ZERO
	 has been selected.  */
      if(_subset_solve_mode!=SUBSET_ZERO)
	{
	  _create_complement_is(rhs_in);

	  
	  size_t is_complement_local_size = rhs_in.local_size()-is_local_size;

	  Vec subvec1 = NULL;	  Mat submat1 = NULL;
	  VecScatter scatter1 = NULL;

	  ierr = VecCreate(this->comm().get(),&subvec1);  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetSizes(subvec1,is_complement_local_size,PETSC_DECIDE);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecSetFromOptions(subvec1);  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScatterCreate(rhs->vec(),_restrict_solve_to_is_complement, subvec1,NULL, &scatter1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScatterBegin(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecScatterEnd(scatter1,_subset_solve_mode==SUBSET_COPY_RHS ? rhs->vec() : solution->vec(),subvec1,INSERT_VALUES,SCATTER_FORWARD);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScale(subvec1,-1.0);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);

#if PETSC_VERSION_LESS_THAN(3,1,0)
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 PETSC_DECIDE,MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else
	  ierr = MatGetSubMatrix(matrix->mat(),
				 _restrict_solve_to_is,_restrict_solve_to_is_complement,
				 MAT_INITIAL_MATRIX,&submat1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif

	  ierr = MatMultAdd(submat1,subvec1,subrhs,subrhs);	  CHKERRABORT(MPI_COMM_WORLD,ierr);

	  ierr = VecScatterDestroy(&scatter1);
	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = VecDestroy(&subvec1);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  ierr = MatDestroy(&submat1);	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	}

      ierr = KSPSetOperators(_ksp, submat, subprecond
	#ifdef OLD_PETSC
			     ,this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN
#endif
      );
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner)	{
	  subprecond_matrix = new PetscMatrixM(subprecond,
						      this->comm());
	  this->_preconditioner->set_matrix(*subprecond_matrix);
          this->_preconditioner->init();
	}
    }
  else    {
      ierr = KSPSetOperators(_ksp, matrix->mat(), precond->mat()
      #ifdef OLD_PETSC
			     ,this->same_preconditioner ? SAME_PRECONDITIONER : DIFFERENT_NONZERO_PATTERN
#endif
			    );
      CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner)      {
	this->_preconditioner->set_matrix(matrix_in);
        this->_preconditioner->init();
      }
    }

  // Set the tolerances for the iterative solver.  Use the user-supplied
  // tolerance for the relative residual & leave the others at default values.
  ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
 			   PETSC_DEFAULT, max_its);    CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Solve the linear system
  if(_restrict_solve_to_is!=NULL){
      ierr = KSPSolveTranspose (_ksp, subrhs, subsolution);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  else    {
      ierr = KSPSolveTranspose (_ksp, rhs->vec(), solution->vec());
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber (_ksp, &its);         CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm (_ksp, &final_resid);         CHKERRABORT(MPI_COMM_WORLD,ierr);

  if(_restrict_solve_to_is!=NULL)    {
      switch(_subset_solve_mode)	{
	case SUBSET_ZERO:
	  ierr = VecZeroEntries(solution->vec());CHKERRABORT(MPI_COMM_WORLD,ierr);
	  break;

	case SUBSET_COPY_RHS:
	  ierr = VecCopy(rhs->vec(),solution->vec());	  CHKERRABORT(MPI_COMM_WORLD,ierr);
	  break;

	case SUBSET_DONT_TOUCH:
	  /* Nothing to do here.  */
	  break;

	}
      ierr = VecScatterBegin(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterEnd(scatter,subsolution,solution->vec(),INSERT_VALUES,SCATTER_REVERSE);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecScatterDestroy(&scatter);   CHKERRABORT(MPI_COMM_WORLD,ierr);

      if(this->_preconditioner)	{
	  /* Before we delete subprecond_matrix, we should give the
	     _preconditioner a different matrix.  */
	  this->_preconditioner->set_matrix(matrix_in);
          this->_preconditioner->init();
	  delete subprecond_matrix;
	  subprecond_matrix = NULL;
	}

      ierr = VecDestroy(&subsolution); CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecDestroy(&subrhs);      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&submat);      CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = MatDestroy(&subprecond);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

#endif

//   STOP_LOG("solve()", "PetscLinearSolver");
  // return the # of its. and the final residual norm.
  return std::make_pair(its, final_resid);
}
// std::pair< int, double>
// PetscLinearSolverM::solve (const ShellMatrix<T>& shell_matrix,
// 			     NumericVector<T>& solution_in,
// 			     NumericVector<T>& rhs_in,
// 			     const double tol,
// 			     const  int m_its)
// {
//
// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   // FIXME[JWP]: There will be a bunch of unused variable warnings
//   // for older PETScs here.
//   std::cout << "This method has been developed with PETSc 2.3.1.  "
// 	    << "No one has made it backwards compatible with older "
// 	    << "versions of PETSc so far; however, it might work "
// 	    << "without any change with some older version." << std::endl;
//   libmesh_error();
//   return std::make_pair(0,0.0);
//
// #else
//
//   START_LOG("solve()", "PetscLinearSolver");
//
//   // Make sure the data passed in are really of Petsc types
//   PetscVectorM* solution = dynamic_cast<PetscVectorM*>(&solution_in);
//   PetscVectorM* rhs      = dynamic_cast<PetscVectorM*>(&rhs_in);
//
//   this->init ();
//
//   int ierr=0;
//   int its=0, max_its = static_cast<int>(m_its);
//   PetscReal final_resid=0.;
//
//   // Close the matrices and vectors in case this wasn't already done.
//   solution->close ();
//   rhs->close ();
//
//   // Prepare the matrix.
//   Mat mat;
//   ierr = MatCreateShell(MPI_COMM_WORLD,
// 			rhs_in.local_size(),
// 			solution_in.local_size(),
// 			rhs_in.size(),
// 			solution_in.size(),
// 			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
// 			&mat);
//   /* Note that the const_cast above is only necessary because PETSc
//      does not accept a const void*.  Inside the member function
//      _petsc_shell_matrix() below, the pointer is casted back to a
//      const ShellMatrix<T>*.  */
//
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
//   ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set operators. The input matrix works as the preconditioning matrix
//   ierr = KSPSetOperators(_ksp, mat, mat,
// 			 SAME_NONZERO_PATTERN);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set the tolerances for the iterative solver.  Use the user-supplied
//   // tolerance for the relative residual & leave the others at default values.
//   ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
//  			   PETSC_DEFAULT, max_its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Solve the linear system
//   ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the number of iterations required for convergence
//   ierr = KSPGetIterationNumber (_ksp, &its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the norm of the final residual to return to the user.
//   ierr = KSPGetResidualNorm (_ksp, &final_resid);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Destroy the matrix.
//   ierr = MatDestroy(mat);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   STOP_LOG("solve()", "PetscLinearSolver");
//   // return the # of its. and the final residual norm.
//   return std::make_pair(its, final_resid);
//
// #endif
//
// }


// std::pair< int, double>
// PetscLinearSolverM::solve (const ShellMatrix<T>& shell_matrix,
// 			     const SparseMatrix<T>& precond_matrix,
// 			     NumericVector<T> &solution_in,
// 			     NumericVector<T> &rhs_in,
// 			     const double tol,
// 			     const  int m_its)
// {
//
// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   // FIXME[JWP]: There will be a bunch of unused variable warnings
//   // for older PETScs here.
//   std::cout << "This method has been developed with PETSc 2.3.1.  "
// 	    << "No one has made it backwards compatible with older "
// 	    << "versions of PETSc so far; however, it might work "
// 	    << "without any change with some older version." << std::endl;
//   libmesh_error();
//   return std::make_pair(0,0.0);
//
// #else
//
//   START_LOG("solve()", "PetscLinearSolver");
//
//   // Make sure the data passed in are really of Petsc types
//   const PetscMatrixM* precond  = dynamic_cast<const PetscMatrixM*>(&precond_matrix);
//   PetscVectorM* solution = dynamic_cast<PetscVectorM*>(&solution_in);
//   PetscVectorM* rhs      = dynamic_cast<PetscVectorM*>(&rhs_in);
//
//   this->init ();
//
//   int ierr=0;
//   int its=0, max_its = static_cast<int>(m_its);
//   PetscReal final_resid=0.;
//
//   // Close the matrices and vectors in case this wasn't already done.
//   solution->close ();
//   rhs->close ();
//
//   // Prepare the matrix.
//   Mat mat;
//   ierr = MatCreateShell(MPI_COMM_WORLD,
// 			rhs_in.local_size(),
// 			solution_in.local_size(),
// 			rhs_in.size(),
// 			solution_in.size(),
// 			const_cast<void*>(static_cast<const void*>(&shell_matrix)),
// 			&mat);
//   /* Note that the const_cast above is only necessary because PETSc
//      does not accept a const void*.  Inside the member function
//      _petsc_shell_matrix() below, the pointer is casted back to a
//      const ShellMatrix<T>*.  */
//
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = MatShellSetOperation(mat,MATOP_MULT,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_mult));
//   ierr = MatShellSetOperation(mat,MATOP_GET_DIAGONAL,reinterpret_cast<void(*)(void)>(_petsc_shell_matrix_get_diagonal));
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Set operators. The input matrix works as the preconditioning matrix
//   ierr = KSPSetOperators(_ksp, mat, const_cast<PetscMatrixM*>(precond)->mat(),
// 			 DIFFERENT_NONZERO_PATTERN);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   if(this->_preconditioner)
//     this->_preconditioner->set_matrix(const_cast<SparseMatrix<Number>&>(precond_matrix));
//
//   // Set the tolerances for the iterative solver.  Use the user-supplied
//   // tolerance for the relative residual & leave the others at default values.
//   ierr = KSPSetTolerances (_ksp, tol, PETSC_DEFAULT,
//  			   PETSC_DEFAULT, max_its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Solve the linear system
//   ierr = KSPSolve (_ksp, rhs->vec(), solution->vec());
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the number of iterations required for convergence
//   ierr = KSPGetIterationNumber (_ksp, &its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Get the norm of the final residual to return to the user.
//   ierr = KSPGetResidualNorm (_ksp, &final_resid);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Destroy the matrix.
//   ierr = MatDestroy(mat);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   STOP_LOG("solve()", "PetscLinearSolver");
//   // return the # of its. and the final residual norm.
//   return std::make_pair(its, final_resid);
//
// #endif
//
// }


// =========================================================================
void PetscLinearSolverM::get_residual_history(std::vector<double>& hist) {
  int ierr = 0;  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p; ierr = KSPGetResidualHistory(_ksp, &p, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check for early return
  if (its == 0) return;
  // Create space to store the result
  hist.resize(its);
  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i) { hist[i] = *p; p++; }
}

// ======================================================
double PetscLinearSolverM::get_initial_residual() {
  int ierr = 0;  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p; ierr = KSPGetResidualHistory(_ksp, &p, &its); CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check no residual history
  if (its == 0) {std::cerr << "No iterations have been performed, returning 0." << std::endl; return 0.; }
  // Otherwise, return the value pointed to by p.
  return *p;
}

// =================================================
void PetscLinearSolverM::set_petsc_solver_type() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CGM:
    ierr = KSPSetType(_ksp, (char*) KSPCG);         CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CRM:
    ierr = KSPSetType(_ksp, (char*) KSPCR);         CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CGSM:
    ierr = KSPSetType(_ksp, (char*) KSPCGS);        CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case BICGM:
    ierr = KSPSetType(_ksp, (char*) KSPBICG);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case TCQMRM:
    ierr = KSPSetType(_ksp, (char*) KSPTCQMR);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case TFQMRM:
    ierr = KSPSetType(_ksp, (char*) KSPTFQMR);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case LSQRM:
    ierr = KSPSetType(_ksp, (char*) KSPLSQR);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case BICGSTABM:
    ierr = KSPSetType(_ksp, (char*) KSPBCGS);       CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case MINRESM:
    ierr = KSPSetType(_ksp, (char*) KSPMINRES);     CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case GMRESM:
    ierr = KSPSetType(_ksp, (char*) KSPGMRES);      CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case RICHARDSONM:
    ierr = KSPSetType(_ksp, (char*) KSPRICHARDSON); CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  case CHEBYSHEVM:
    ierr = KSPSetType(_ksp, (char*) KSPCHEBYSHEV);  CHKERRABORT(MPI_COMM_WORLD,ierr); return;
//   case PREONLYM:  
//     ierr = KSPSetType(_ksp, (char*) KSPPREONLY);    CHKERRABORT(MPI_COMM_WORLD,ierr); return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}

// =======================================================
void PetscLinearSolverM::print_converged_reason() {
// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   std::cout << "This method is currently not supported "
//             << "(but may work!) for Petsc 2.3.0 and earlier." << std::endl;
// #else
  KSPConvergedReason reason;  KSPGetConvergedReason(_ksp, &reason);
 std::cout << "Linear solver convergence/divergence reason: " << KSPConvergedReasons[reason] << std::endl;
// #endif
}

// PetscErrorCode PetscLinearSolverM::_petsc_shell_matrix_mult(Mat mat, Vec arg, Vec dest)
// {
//   /* Get the matrix context.  */
//   int ierr=0;
//   void* ctx;
//   ierr = MatShellGetContext(mat,&ctx);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   /* Get user shell matrix object.  */
//   const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
//
//   /* Make \p NumericVector instances around the vectors.  */
//   PetscVectorM arg_global(arg);
//   PetscVectorM dest_global(dest);
//
//   /* Call the user function.  */
//   shell_matrix.vector_mult(dest_global,arg_global);
//
//   return ierr;
// }

// =====================================================
// PetscErrorCode PetscLinearSolverM::_petsc_shell_matrix_get_diagonal(Mat mat, Vec dest){
//   /* Get the matrix context.  */
//   int ierr=0;
//   void* ctx;
//   ierr = MatShellGetContext(mat,&ctx);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   /* Get user shell matrix object.  */
//   const ShellMatrix<T>& shell_matrix = *static_cast<const ShellMatrix<T>*>(ctx);
//
//   /* Make \p NumericVector instances around the vector.  */
//   PetscVectorM dest_global(dest);
//
//   /* Call the user function.  */
//   shell_matrix.get_diagonal(dest_global);
//
//   return ierr;
// }



//------------------------------------------------------------------
// Explicit instantiations
// template class PetscLinearSolver<Number>;



#endif // #ifdef LIBMESH_HAVE_PETSC
