// ======================================================================================
// --------------   NAVIER-STOKES system [FS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"

#ifdef FSI_EQUATIONS
// #if FSI_EQUATIONS==1
// ======================================================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverFSI.h"       // Navier-Stokes class header file


// config file --------------------------------------------------------------------------
#include "MGGeomEl.h"        // Geometrical element
#include "MGFE_conf.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ======================================================================================
/**     \addtogroup user_function User functions
 * @{
 */




/// \ingroup  user_function ic_read
// ======================================================================================
/// This function generates the initial conditions for the NS system:
void MGSolFSI::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]
) {// ===================================================================================

  double ILref = 1./_lref;

  // =======================================================================================
  //    boundary conditions box file plane_ch_rotate.med ==================================
   u_value[0] = 0.;
  u_value[1] = 0.;
//    if(xp[1]<-0.2499 ) u_value[1] = 1.;
//   if(bc_gam==11)u_value[1] = 2.;
  if(xp[1]<-0.24999 &&  bc_gam<20 )u_value[1] = .5;
//   u_value[1] = 2.;
//   if (xp[0] > .1999 && xp[1]>0.6999){ u_value[0] =2.; u_value[1] = 0.;}
//    if (xp[1] > .899999){ u_value[0] =2.; u_value[1] = 0.;}
//   if (xp[1] < 0.1) u_value[1] =0.1*1.5*(1.-xp[0])*(1.+xp[0])/(0.25);
  // cyl
//     if (xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =-0.25e+3*(xp[0]-LXE*ILref)*(LXE*ILref+xp[0])/_uref;//0.0165; 
    // annulus
//     if (xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =1.; 
  
    double pref = _refvalue[DIMENSION];
//     u_value[DIMENSION] = 15.;
      u_value[DIMENSION] = 10.-3*(xp[1]-0.25);

  
 
  return;
}

/// \ingroup  user_function bc_read
// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolFSI::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], 	// normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0  ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 ->  tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;
   bc_Neum[0] =wall; bc_flag[0]=0;
   bc_Neum[0] = _FSI_parameter._map_FSIgroup[bc_gam];
   
   // channel_fsi1t.med ===============================================
  
//  if (bc_gam == 1000 && xp[1] < 0. ) {bc_Neum[0]=wall;} //inlet
//  if (bc_gam == 21) {bc_Neum[0]=wall;}     //solid bottom
//  if (bc_gam == 20) {bc_Neum[0]=outflow;}        //solid right
//  if (bc_gam == 15) {bc_Neum[0]=wall;}        //liquid left
//  if (bc_gam == 11)  {bc_Neum[0]=velocity;} //liquid bottom
//  if (bc_gam == 13) {bc_Neum[0]=pressure_inlet;} //liquid top
//  if (bc_gam == 1000 && xp[1] > 0. ) {bc_Neum[0]=wall;} //top interface
//  if (bc_gam == 23) {bc_Neum[0]=wall;}        //solid top
//  if (xp[0]>0.2999) {bc_Neum[0]=outflow;}     //solid right
  if (xp[0] > 0.2999 && xp[1] < -0.2499) {bc_Neum[0]=wall;} //right vertex
  if (xp[0] > 0.2999 && xp[1] >  0.2499) {bc_Neum[0]=wall;}
   
   
  
  return;
}









#endif  //ENDIF FSI_EQUATIONS
