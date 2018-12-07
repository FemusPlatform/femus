// *****************************************************************
//      MGGeomElconf
// *****************************************************************

#ifndef _quadtri_cfg
#define _quadtri_cfg

// ======= GEOMETRIC ELEMENT TYPE ====
//Geometric element,
// only quadratic,
// volume and boundary

//here no structure, but the numbers

#include "Domain_conf.h"


//***********
#define ELTYPE 27        //quadrilateral
// #define ELTYPE 10     //triangular 
//***********

#define NUMGEOMELS 1  //=========number of geometric element types


//***********
#if ELTYPE == 27

#if DIMENSION==2 
  #define  NNDS        (9)  //Lagrange Quad9
  #define  NNDS_B      (3)  //Lagrange Edge3
  #define  MED_EL_TYPE INTERP_KERNEL::NORM_QUAD9
  #define  MED_EL_BDTYPE INTERP_KERNEL::NORM_SEG3
#endif
#if DIMENSION==3 
  #define  NNDS       (27)
  #define  NNDS_B      (9)  //Lagrange 
  #define  MED_EL_TYPE INTERP_KERNEL::NORM_HEXA27
   #define  MED_EL_BDTYPE INTERP_KERNEL::NORM_QUAD9
#endif

#endif

#if ELTYPE == 10

#if DIMENSION==2 // ----------------------
#define  NNDS        (6)  //Lagrange Tri6
#define  NNDS_B      (3)  //Lagrange Edge3
  #define  MED_EL_TYPE INTERP_KERNEL::NORM_TRI6 
  #define  MED_EL_BDTYPE INTERP_KERNEL::NORM_SEG3
#endif

#if DIMENSION==3  // ------------------------ 
#define  NNDS       (10)
#define  NNDS_B      (6)
#define  MED_EL_TYPE INTERP_KERNEL::NORM_TETRA10
#define  MED_EL_BDTYPE INTERP_KERNEL::NORM_TRI6
#endif

#endif

 
// *****************************************************************
//      MGFEconf
// *****************************************************************

/*
enum FECouple {
 Q2Q1=0, Q2Q0, Q2S0, Q2P1 
}; 

static const FECouple myfemc = Q2Q1;      
#define  Q2Q   0 */



#define  Q2Q1 (1)
#define  NDOF_K      (0)
#define  NDOF_BK      (0)

// #define  Q2Q0 (1)
// #define  NDOF_K      (1)
// #define  NDOF_BK      (1)
// 
// #define  Q2P1 (1)
// #define  NDOF_K      (3)
// #define  NDOF_bK      (3)

#if DIMENSION==2 // ===============================

// Quad9-Quad8-Quad4
#if ELTYPE == 27

#define  NDOF_FEM    (9)
#define  NDOF_FEMB   (3)
#define  NDOF_P      (4)
#define  NDOF_PB     (2)
// #define  NDOF_K      (myfemc)
#define  NSUBDOM    (4)
#endif
// Tri6-Tri3
#if ELTYPE == 10

#define  NDOF_FEM    (6)
#define  NDOF_FEMB   (3)
#define  NDOF_P      (3)
#define  NDOF_PB     (2)
#define  NSUBDOM    (3)
#endif

// 2d axisymmetric case -----------------
//#define AXISYMX
#endif


#if DIMENSION==3 // ========================

// Hex27-Hex20-Hex8
#if ELTYPE==27

#define  NDOF_FEM   (27)
#define  NDOF_FEMB   (9)
#define  NDOF_P      (8)
#define  NDOF_PB     (4)
// #define  NDOF_K      (1+myfemc)
#define  NSUBDOM    (8)
#endif
// Tetra10-Tetra4
#if ELTYPE==10

#define  NDOF_FEM   (10)
#define  NDOF_FEMB   (6)
#define  NDOF_P      (4)
#define  NDOF_PB     (3)
#define  NSUBDOM    (4)
#endif

#endif

#endif


