#ifndef _domainconf
#define _domainconf
//Dimension of the problem
// DIMENSION = 2 or 3 (D)
#define DIMENSION  (3)
#define NUM_MESH (1)
#define  BDRY_TOLL  1.e-12 //tolerance for setting the BCs
// libemesh no this 
#define MATBC_INTERFACE       // boundary conditions and material
// #define HAVE_GROUP (1) // see in the gencase file only for med
#endif
