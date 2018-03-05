#include "Solverlib_conf.h"
#ifndef HAVE_LASPACKM

#include <iostream>
// this class conf
#include "MGFemusInit.h"
#include "Printinfo_conf.h"
// MPI -------------------------
#if defined HAVE_MPI
# include <mpi.h>
#endif
// PETSC ----------------------
#ifdef HAVE_PETSCM
# include "petsc_macroM.h"
EXTERN_C_FOR_PETSC_BEGIN
# include <petsc.h>
# include <petscerror.h>
EXTERN_C_FOR_PETSC_END
#endif


// ============================================================================
// MGFemusInit data -----------------------------------------------------------
// ============================================================================
//  ParallelM::Communicator _comm;   ///< parallel communicator
//  bool _is_initialized;            ///< initialized flag
//  bool _is_initialized_mpi;        ///< mpi flag 
//  bool _is_initialized_petsc;      ///< petsc flag
//  int _size;                       ///< n processors
//  int _rank;                       ///< id processor


// =============================================================================
// MGFemusInit functions--------------------------------------------------------
// =============================================================================
// =============================================================================
/// This function initializes the libraries if it is parallel:
/// flags to initialize are: _is_initialized (femus),
/// _is_initialized_mpi (MPI), _is_initialized_petsc (PETSC)
#ifndef HAVE_MPI
MGFemusInit::MGFemusInit(int argc, const char* const* argv)
#else
MGFemusInit::MGFemusInit(
  int & argc,            // integer program input
  char** & argv,         // char program input
  MPI_Comm comm_world_in // communicator for MPI direct
)
#endif
{// ============================================================================

  // setting default
  _is_initialized_mpi=false;          // MPI initialization (a)
  _is_initialized_petsc=false;        // PETSC initialization (b)
  _is_initialized=false;              // femus iniitalization (c)
  assert (_is_initialized==false);
  
  // a) check  MPI initialization  ******************************************************
#ifdef HAVE_MPI
  int flag; MPI_Initialized (&flag);
  if (!flag) {
    MPI_Init (&argc, const_cast<char***>(&argv));
    _is_initialized_mpi = true;
  }
  // Duplicate the input communicator for internal use
  // And get a Parallel::Communicator copy too, to use
  // as a default for that API
  this->_comm = comm_world_in;
  _size = static_cast< int>(this->comm().size());
  assert (_size >= 0);
  _rank = static_cast<unsigned int>(this->comm().rank());
  assert (_rank >= 0);
  // Let's be sure we properly initialize on every processor at once:
  parallel_onlyM(this->comm());
#endif
  
 // b) check  PETSC initialization ******************************************************
#ifdef HAVE_PETSCM
  int ierr=0;
  PETSC_COMM_WORLD = comm_world_in;
  // Check whether the calling program has already initialized
  // PETSc, and avoid duplicate Initialize/Finalize
  PetscBool petsc_already_initialized;
  ierr = PetscInitialized(&petsc_already_initialized); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  if (petsc_already_initialized != PETSC_TRUE) _is_initialized_petsc= true;
  ierr = PetscInitialize (&argc, &argv, NULL, NULL); CHKERRABORT(PETSC_COMM_WORLD,ierr);

#ifdef  PRINT_PROC  // ------- info ------------------	      
  std::cout << " MGFemusInit: PETSC_COMM_WORLD initialized \n";
#endif // ------ end info ----------------------------
#endif

#ifdef  PRINT_PROC  // ------- info ------------------	
  std::cout << " MGFemusInit: MPI_COMM_WORLD initialized from proc " << _rank << "\n";
#endif // ------ end info ----------------------------
  
  // c) Femus initialization ************************************************************
  // redirect rank!=0 to null POINTER
  if ( _rank != 0 ) {std::cout.rdbuf(NULL);}
  _is_initialized= true; // we go!
  
  return;

}
// // =======================================================2014year
// /// This function initializes the libraries if it is parallel
// MGFemusInit::~MGFemusInit() { // ========================
//
// // 1 proc= nothing to do
// #ifdef HAVE_PETSCM
//       PetscFinalize();
// #ifdef  PRINT_PROC  // ------- info ------------------
//  std::cout << " ~MGFemusInit(): PETSC_COMM_WORLD ends \n";
// #endif // ------ end info ----------------------------
// #endif
//
// #ifdef HAVE_MPI
// //       MPI_Comm_free (MPI_COMM_WORLD);
//        MPI_Finalize();
// #ifdef  PRINT_PROC  // ------- info ------------------
//   std::cout << " ~MGFemusInit(): MPI_COMM_WORLD ends \n";
// #endif // ------ end info ----------------------------
// #endif
//
//  return;
// }
// =======================================================
/// This function initializes the libraries if it is parallel
MGFemusInit::~MGFemusInit()   // ========================
{
  // Let's be sure we properly close on every processor at once:
  parallel_onlyM(this->comm());
#if defined(HAVE_MPI)
  // We may be here in only one process,
  // because an uncaught libmesh_error() exception
  // called the LibMeshInit destructor.
  //
  // If that's the case, we need to MPI_Abort(),
  // not just wait for other processes that
  // might never get to MPI_Finalize()
  if (_is_initialized_mpi &&  std::uncaught_exception()) {
    std::cerr << "Uncaught exception - aborting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    abort();
  }
#endif
#if defined(HAVE_PETSCM)
  if (_is_initialized_petsc)  PetscFinalize();
#endif
#if defined(HAVE_MPI)
  // Allow the user to bypass MPI finalization
  this->_comm.clear();
//        ParallelM::Communicator_World.clear();
  if (_is_initialized_mpi) MPI_Finalize();
#endif
  // Set the initialized() flag to false
  _is_initialized = false;
  return;
}
#endif // end  HAVE_LASPACKM

