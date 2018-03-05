#ifndef _femus_init_
#define _femus_init_

#include "Solverlib_conf.h"
#include "parallelM.h"
#include "parallel_objectM.h"

#ifndef HAVE_LASPACKM
//use this class only with petsc-mpi
#ifdef HAVE_MPI
#include <mpi.h>  //For MPI_COMM_WORLD
#endif
// ====================================================
/// Initializes the MPI_COMM_WORLD environment variable for MPI_Comm
/// and PETSC for algebraic solutions.
// Note:
// Initialize the library for use, with the command line options
// provided.  This will e.g. call PetscInitialize if PETSC is
// available.  You must create a LibMeshInit object before using any
// of the library functionality.  This method may take an optional
// parameter to use a user-specified MPI communicator.
// ====================================================

class MGFemusInit {

private:
  // ========================================
  //  data -----------------------------------------------
  // ========================================
  ParallelM::Communicator _comm;   ///< parallel communicator 
  bool _is_initialized;            ///< initialized flag      
  bool _is_initialized_mpi;        ///< mpi flag
  bool _is_initialized_petsc;      ///< petsc flag
  int _size;                       ///< n processors
  int _rank;                       ///< id processor

public:
  // ============================================
  // Constructor-Destructor
  // ============================================
 
#ifdef HAVE_MPI
  // Constructor-Desctructor ------------------------------------------------------------
  /// Constructor
  MGFemusInit(
    int & argc,char** & argv,              ///< inline commands (no use in)
    MPI_Comm COMM_WORLD_IN=MPI_COMM_WORLD  ///< MPI communicator(in)
  );
#else
  MGFemusInit(
    int argc, const char* const* argv      ///< inline commands (no use in)
  );
#endif
  /// Destructor
  virtual ~MGFemusInit();
  
// return functions ---------------------------------------------------------------------
  const ParallelM::Communicator& comm() {return _comm; }

};

#endif // end no HAVE_LASPACKM -------------
#endif // end _femus_init ----------------------

