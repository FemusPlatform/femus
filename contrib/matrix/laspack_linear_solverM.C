#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM
// -----------------------------------
// C++ includes
// -----------------------------------
// Local Includes
#include "laspack_linear_solverM.h" // Numeric
#include "sparse_matrixM.h"      // Numeric
#include "laspack_matrixM.h"     // Numeric
#include "numeric_vectorM.h"     // Numeric
#include "laspack_vectorM.h"     // Numeric
#include "linear_solverM.h"      // Numeric
#include "mlsolv.h"              // Laspack
#include "factor.h"              // Laspack
// -----------------------------------
// #include "libmesh_common.h"     
// Numeric
/*----------------------- functions ----------------------------------*/

void LaspackLinearSolverM::clear (){
  if (this->initialized())    {
      this->_is_initialized = false;
      this->_solver_type         = GMRESM;
      this->_preconditioner_type = ILU_PRECONDM;
    }
}

// ============================================
void LaspackLinearSolverM::init (){
  // Initialize the data structures if not done so already.
  if (!this->initialized())      this->_is_initialized = true;
 // SetRTCAuxProc (print_iter_accuracy);
}


// ==================================================
std::pair<unsigned int, Real> 
LaspackLinearSolverM::solve (SparseMatrixM &matrix_in,
			       NumericVectorM &solution_in,
			       NumericVectorM &rhs_in,
			       const double tol,
			       const unsigned int m_its)
{
//   START_LOG("solve()", "LaspackLinearSolver");
  this->init ();

  // Make sure the data passed in are really in Laspack types
  LaspackMatrixM* matrix   = static_cast<LaspackMatrixM*>(&matrix_in);
  LaspackVectorM* solution = static_cast<LaspackVectorM*>(&solution_in);
  LaspackVectorM* rhs      = static_cast<LaspackVectorM*>(&rhs_in);  

  // Zero-out the solution to prevent the solver from exiting in 0
  // iterations (?)
  //TODO:[BSK] Why does Laspack do this?  Comment out this and try ex13...
//   solution->zero();
  matrix->close ();  solution->close ();  rhs->close ();
  // Set the preconditioner type
  this->set_laspack_preconditioner_type ();

  // Set the solver tolerance
  SetRTCAccuracy (tol);

  // Solve the linear system
  switch (this->_solver_type)  {
    case CGM:  {// Conjugate-Gradient
	CGIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break; }
    case CGNM:  {  // Conjugate-Gradient Normalized
	CGNIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break; }   
    case CGSM: {// Conjugate-Gradient Squared
	CGSIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break;  }
    case BICGM:  {// Bi-Conjugate Gradient
	BiCGIter (&matrix->_QMat,&solution->_vec,&rhs->_vec, m_its, _precond_type, 1.);
	break;  }  
    case BICGSTABM: {// Bi-Conjugate Gradient Stabilized
	BiCGSTABIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break;  }
    case QMRM:  { // Quasi-Minimum Residual
	QMRIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break; }
    case SSORM:      {// Symmetric over-relaxation
	SSORIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break;  }
    case JACOBIM:{// Jacobi Relaxation
	JacobiIter (&matrix->_QMat, &solution->_vec,&rhs->_vec, m_its,_precond_type,1.);
	break;   }
    case GMRESM:{// Generalized Minimum Residual
	SetGMRESRestart (30);
	GMRESIter (&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
	break;    }
	#ifdef HAVE_LASPACKM
     case VANKATM:{// Vanka Temperature
       VankaIter(&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
 	break;    }
     case VANKANSM:{// Vanka NavierStokes
// //       std::cout << "DDDDDDDDDDDDDDD" << std::endl;
      VankaNSIter(&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
 	break;    }
     case LUMPM:{// Vanka NavierStokes
//        std::cout << "DDDDDDDDDDDDDDD" << std::endl;
       LumpIter(&matrix->_QMat,&solution->_vec,&rhs->_vec,m_its,_precond_type,1.);
 	break;    }
 	#endif
      // Unknown solver, use GMRES
    default:      {
	std::cerr << "ERROR:  Unsupported LASPACK Solver: "
		  << this->_solver_type      << std::endl
		  << "Continuing with GMRES" << std::endl;
	
	this->_solver_type = GMRESM;
	return this->solve (*matrix,*solution,*rhs,tol,m_its);
      }
    }
  // Check for an error
  if (LASResult() != LASOK){
      std::cerr << "ERROR:  LASPACK Error: " << std::endl;
      WriteLASErrDescr(stdout);
      abort();
   }
//   STOP_LOG("solve()", "LaspackLinearSolver");
  // Get the convergence step # and residual 
  return std::make_pair(GetLastNoIter(), GetLastAccuracy());
}

/* put out accuracy after each multigrid iteration */
void IterStatusM(int Iter,double rNorm,double bNorm,IterIdType IterId) {
  if (IterId == MGIterId || IterId == MGPCGIterId || IterId == BPXPCGIterId) {
    printf("%3d. iteration ... accuracy = ", Iter);
    if (!_LPIsZeroReal(bNorm))  printf("%11.4e\n", rNorm / bNorm);
    else     printf("    ---\n");
  }
}
// -------------------------------------------------------
// solution of discret problem with multigrid solver
void LaspackLinearSolverM::MGSolve(std::vector<SparseMatrixM*> &matrix_in,
			       std::vector<NumericVectorM*> &solution_in,
			       std::vector<NumericVectorM*> &rhs_in,
			       Matrix *P, Matrix *R,  
			       const double /*tol*/,
			       const unsigned int/* m_its*/) {
  uint   MaxIter  = 15; 
  double Eps      = 1.e-6;
  int Gamma       = 2; // ciclo V W
  int NoMGIter    =1;  //  
//   int RestrType   =1;
  int  Nu1        =       8; // Number of post-smoothing iterations  
  int Nu2=8;  // Relaxation parameter for smoothing: 
  double Omega=0.98;
  
  IterIdType MLSolverId= MGIterId;
 _LPBoolean  NestedMG  =_LPFalse; 
  IterProcType SmoothProc=GMRESIter;// Smoother
  PrecondProcType PrecondProc=ILUPrecond;// Preconditioner

  // Coarse Grid 
  int NuC     =40; 
  double OmegaC=0.98;
  IterProcType SolvProc      =GMRESIter;// Smoother coarse level
  PrecondProcType PrecondProcC=ILUPrecond;// Preconditioner coarse level
  
  
  int NoIter; /* number of performed iterations */
  int Level;
  double AccBeg, AccEnd; /* reached accuracy (of residuum) */
  double ContrRatePerMGIter, ContrRatePerSec;
  double ClockFactor, CPUTime;
  clock_t BegClock, EndClock;
  uint NoLevels=matrix_in.size();
   this->init ();

  // Make sure the data passed in are really in Laspack types
  QMatrix *LA=new QMatrix[NoLevels];
  QVector *uwx=new QVector[NoLevels];
  QVector *frb=new QVector[NoLevels];
  for(uint level=0;level<NoLevels;level++){
     LaspackMatrixM *tmp = dynamic_cast<LaspackMatrixM *>(matrix_in[level]);
      LA[level] = tmp->_QMat;
     uwx[level]= (dynamic_cast<LaspackVectorM *>(solution_in[level]))->_vec;
     frb[level]= (dynamic_cast<LaspackVectorM *>(rhs_in[level]))->_vec;
  }
  
  /* setting of RTC parameters */
  SetRTCAccuracy(Eps); SetRTCAuxProc(IterStatusM);

  // new !!!!!!!!!!!!!!!!!!!!!!
//      V_SetAllCmp(&uw[NoLevels - 1], 0.0);
//      RungeKutta4(&L[NoLevels - 1], &uw[NoLevels - 1],&fr[NoLevels - 1],50,
// 	   PrecondProcC,5.);
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  exit(0);

  // Eigenvalues *************************************
  /* estimate extremal eigenvalues */
//   if (SmoothProc == ChebyshevIter || SolvProc == ChebyshevIter) {
//     /* start of time counting */
//     ClockFactor = 1.0 / (double)CLOCKS_PER_SEC;   BegClock = clock();
// 
//     /* initialization random-number generator */
//     srand(1);
//     printf(" Estimating extremal eigenvalues : ");
//     if (SmoothProc == ChebyshevIter) {
//       for (Level = NoLevels - 1; Level > 0; Level--) {
//         GetMinEigenval(&L[Level], PrecondProc, Omega);
//         GetMaxEigenval(&L[Level], PrecondProc, Omega);
//       }
//     }
//     if (SolvProc == ChebyshevIter) {
//       GetMinEigenval(&L[0], PrecondProcC, OmegaC);
//       GetMaxEigenval(&L[0], PrecondProcC, OmegaC);
//     }
// 
//     /* end of time counting and geting out the CPU time */
//     EndClock = clock();
//     CPUTime = (double)(EndClock - BegClock) * ClockFactor;
//     printf(" CPU time: %7.2f s\n", CPUTime);
//   }
  // Preparation preconditioners ***************************
  /* factorize matrices */
  if (PrecondProc == ILUPrecond || PrecondProcC == ILUPrecond) {
    /* start of time counting */
    ClockFactor = 1.0 / (double)CLOCKS_PER_SEC; BegClock = clock();

    printf(" Incomplete factorization of matrices : ");
    if (PrecondProc == ILUPrecond    && SmoothProc != JacobiIter
        && SmoothProc != SORForwIter && SmoothProc != SORBackwIter
        && SmoothProc != SSORIter) {
      for (Level = NoLevels - 1; Level > 0; Level--)
         ILUFactor(&LA[Level]);
    }
    if (PrecondProcC == ILUPrecond    && SolvProc != JacobiIter
        && SolvProc != SORForwIter    && SolvProc != SORBackwIter
        && SolvProc != SSORIter) {
        ILUFactor(&LA[0]);
    }

    /* end of time counting and geting out the CPU time */
    EndClock = clock();
    CPUTime = (double)(EndClock - BegClock) * ClockFactor;
    printf(" CPU time: %7.2f s\n", CPUTime);
  }
  // end preparation preconditioners ********************



  /* solving of system of equations by nested multigrid method */
  if (MLSolverId == MGIterId && NestedMG) {
    /* initialisation of vector of unknowns */
    V_SetAllCmp(&uwx[0],0.0);

    /* approximation of solution by nested multigrid method */
    printf("\n Doing nested multigrid iterations ...\n");
    NestedMGIter(NoLevels,LA, uwx, frb, R, P, Gamma,
                 SmoothProc, Nu1, Nu2, PrecondProc, Omega,SolvProc,
                 NuC,  PrecondProcC ,OmegaC);
    AccBeg = GetLastAccuracy();
  }
  else {
    /* initialisation of vector of unknowns */
    // restart *****************************
     V_SetAllCmp(&uwx[NoLevels-1],0.0);
    // *************************************
    AccBeg = 1.0;
  }

  /* start of time counting */
  ClockFactor = 1.0 / (double)CLOCKS_PER_SEC;
  BegClock = clock();

  /* solving of system of equations */
  switch (MLSolverId) {
    case MGIterId:
      printf("\n Doing multigrid iterations \n");
      MGIter(NoLevels, LA, uwx, frb, R, P, MaxIter, Gamma,
             SmoothProc, Nu1, Nu2, PrecondProc, Omega,
             SolvProc, NuC, PrecondProcC, OmegaC);
// 	   MGIter(NoLevels, L, uw, fr, R, P, MaxIter, Gamma,
//              SmoothProc, Nu1, Nu2, PrecondProc, Omega,
//              SolvProc, NuC, PrecondProcC, OmegaC);
      break;
    case MGPCGIterId:
      printf(" Doing multigrid preconditioned CG iterations \n");
      MGPCGIter(NoLevels,  LA, uwx, frb, R, P, MaxIter, NoMGIter, Gamma,
                SmoothProc, Nu1, Nu2, PrecondProc, Omega,
                SolvProc, NuC, PrecondProcC, OmegaC);
      break;
    case BPXPCGIterId:
      printf(" Doing BPX preconditioned CG iterations \n");
      BPXPCGIter(NoLevels,  LA, uwx, frb, R, P, MaxIter,
                 SmoothProc, Nu1, PrecondProc, Omega,
                 SolvProc, NuC, PrecondProcC, OmegaC);
      break;
    default:
      break;
  }

  /* end of time counting and geting out of CPU time */
  EndClock = clock();
  CPUTime = (double)(EndClock - BegClock) * ClockFactor;
  printf("  CPU time: %7.2f s\n", CPUTime);

  AccEnd = GetLastAccuracy();
  NoIter = GetLastNoIter();

  /* computing of middle contraction rates */
  if (NoIter > 0) ContrRatePerMGIter = pow(AccEnd / AccBeg, 1.0 / (double)NoIter);
  else
    ContrRatePerMGIter = 0.0;
  if (CPUTime > DBL_EPSILON) ContrRatePerSec = pow(AccEnd / AccBeg, 1.0 / CPUTime);
  else
    ContrRatePerSec = 0.0;
// #ifdef PRINT_INFO
  printf(" Middle contraction rate :");
  printf(" Referred to one iteration: %10.3e ;", ContrRatePerMGIter);
  printf(" to 1 s CPU time:  %10.3e\n", ContrRatePerSec);
// #endif
  //   for(int Level=1;Level<V_GetDim(&uw[NoLevels - 1])-1;Level++) printf(" %g ",V_GetCmp(&uw[NoLevels - 1],Level));
  /* LASPack error messages */
  if (LASResult() != LASOK) {
    printf("LASPack error: ");
    WriteLASErrDescr(stderr);
  }
  fprintf(stderr, "\n");
  
  
//    for(uint level=0;level<NoLevels;level++){
//      Q_Destr(&LA[level]);
//      V_Destr(&uwx[level]);
//      V_Destr(&frb[level]);
//   }
//   delete [] LA;delete [] uwx;delete [] frb;
//   std::pair<uint,double> rest=_solver->solve(*A[NoLevels-1],*x[NoLevels-1],*b[NoLevels-1],1.e-6,40);

  //   GMRESIter(&A[NoLevels-1]->_QMat,&x[NoLevels-1]->_vec,&->_vec,40,ILUPrecond,1.);
//   std::cout << "computed solution \n iter= " << rest.first  <<"\n Res= " << rest.second <<"\n";
  
}


// void LaspackLinearSolverM::MGStep(int NoLevels, int Level, int Gamma)
// /* one multigrid iteration */
// {
//     int CoarseMGIter; /* multi grid iteration counter for coarser grid */
// 
//     if (Level == 0) {
//         /* solving of system of equations for the residual on the coarsest grid */
//         (*SolvProc)(&A[Level], &x[Level], &b[Level], NuC, PrecondProcC, OmegaC);
//     } else {
//         /* pre-smoothing - Nu1 iterations */
//         (*SmoothProc)(&A[Level], &x[Level], &b[Level], Nu1, PrecondProc, Omega);
//         /* restiction of the residual to the coarser grid */
//         Asgn_VV(&b[Level - 1], Mul_MV(&R[Level - 1],
// 	    Sub_VV(&b[Level], Mul_QV(&A[Level], &x[Level]))));
//         /* initialisation of vector of unknowns on the coarser grid */
//         V_SetAllCmp(&x[Level - 1], 0.0);
//         /* solving of system of equations for the residual on the coarser grid */
//         for (CoarseMGIter = 1; CoarseMGIter <= Gamma; CoarseMGIter++)
//             MGStep(NoLevels, A, x, b, R, P, Level - 1, Gamma,
// 		   SmoothProc, Nu1, Nu2, PrecondProc, Omega,
//                    SolvProc, NuC, PrecondProcC, OmegaC);
//         /* interpolation of the solution from the coarser grid */
// 	if (P != NULL)
//             AddAsgn_VV(&x[Level], Mul_MV(&P[Level], &x[Level - 1]));
// 	else
//             AddAsgn_VV(&x[Level], Mul_MV(Transp_M(&R[Level - 1]), &x[Level - 1]));
//         /* post-smoothing - Nu2 iterations */
//         (*SmoothProc)(&A[Level], &x[Level], &b[Level], Nu2, PrecondProc, Omega);
//     }
// 
//     return(&x[Level]);
// }


// ================================================================
void LaspackLinearSolverM::set_laspack_preconditioner_type (){
  switch (this->_preconditioner_type)  {
    case IDENTITY_PRECONDM: _precond_type = NULL; return;
    case ILU_PRECONDM:      _precond_type = ILUPrecond; return;
    case JACOBI_PRECONDM:   _precond_type = JacobiPrecond; return;
    case SSOR_PRECONDM:     _precond_type = SSORPrecond; return;
    default:
      std::cerr << "ERROR:  Unsupported LASPACK Preconditioner: "
		<< this->_preconditioner_type << std::endl
		<< "Continuing with ILU"      << std::endl;
      this->_preconditioner_type = ILU_PRECONDM;
      this->set_laspack_preconditioner_type();      
    }
}
// ==========================================================
void LaspackLinearSolverM::print_converged_reason(){
  std::cout << "print_converged_reason() is not supported" << std::endl;
}



//------------------------------------------------------------------
// Explicit instantiations
// template class LaspackLinearSolver<Number>;
 

#endif // #ifdef LIBMESH_HAVE_LASPACK
