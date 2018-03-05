#ifndef __solverlibconf_h__
#define __solverlibconf_h__


// ======== SOLVER LIBRARY ======
//   #define HAVE_LASPACKM
    #define HAVE_PETSCM 
    #define HAVE_MPI   //what if I want to use Petsc without MPI?

    
// Med library
#define    HAVE_MED
    
#define LSOLVER LASPACK_SOLVERSM  
#ifdef HAVE_PETSCM
#undef LSOLVER
 #define LSOLVER  PETSC_SOLVERSM
#endif



//******PETSC VERSION ***88
#ifdef HAVE_PETSCM


/* PETSc's major version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MAJOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MAJOR  3 
#endif

/* PETSc's minor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MINOR  2
#endif

/* PETSc's subminor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR  0 
#endif

#endif



// #define HDF5_VERSIONM 1812

#define HDF5_VERSIONM 1810

// #define HDF5_VERSIONM 188



#endif