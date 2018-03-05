#include "Equations_conf.h"



// ============================================
#ifdef DS_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverDS.h"
#include "MGGeomEl.h"

// configuration files -----------
#include "MGFE_conf.h"
#include "Printinfo_conf.h"
#include "MGEquationsSystem.h"

// class local configuration -------


// local include -------------------
// #include "MGMesh.h"
 #include "MGSystem.h"
// #include "numeric_vectorM.h"
// #include "dense_vectorM.h"
// #include "sparse_matrixM.h"
// #include "dense_matrixM.h"
// #include "linear_solverM.h"
// #include "parallelM.h"





// ========================================
/// This function  defines the boundary conditions for the NS system:
// void MGSolDS::bc_intern_read(
//   int bc_gam,
//   int mat_gam,
//   double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
//   int bc_Neum[], // normal
//   int bc_flag[]  // boundary condition flag
// ) {// ===================================
// //    if (xp[0]<(LXE-LXB)/2-BDRY_TOLL){ bc_flag[0]=1;bc_Neum[0]=1; }
// //    if (xp[0]>(LXE-LXB)/2+BDRY_TOLL){ bc_flag[0]=1;bc_Neum[0]=3 ;}
// //   if (fabs(xp[0]-(LXE-LXB)/2 )<BDRY_TOLL) { bc_flag[0]=1;bc_Neum[0]=5;}
// //    if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
// //    if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
// //    if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
//   return;
// }

void MGSolDS::ic_read(
  int /*bc_gam*/,
  int /*mat_gam*/,
  double /*xp*/[],
  int iel,
  double u_value[]
) {// =======================================
   u_value[0] = 0.;
  return;
}

// ============================================================================
/// This function  defines the boundary conditions for the DS system:
void MGSolDS::bc_read(
  int bc_gam,
  int mat_gam,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]         // boundary condition flag
) {// =========================================================================
///   enum bound_condDS --------------------------------------------------------
///   Diriclet  (point constraint < 10) ----------------------------------------
///     simm=0, disp_in0=1,disp_tg0=2,fix_disp0_fix=3,
///     disp_in=5,disp_tg=6,fix_disp0_disp=7,
///    Neumann (integratral constraint >9) -------------------------------------
///     free_disp=10, free_dispn=12,interior=11,slip=13,
// ----------------------------------------------------------------------------
  
//   double ILref = 1./_lref;
#if DIMENSION==2  // ----- 2D boundary conditions ----------
  

  
//   if (bc_gam< 20){bc_Neum[0]== fix_disp0_fix; }
//   else if (bc_gam<30) {bc_Neum[0]= free_disp;}
//   else if (bc_gam>900)   {bc_Neum[0]=fix_disp0_fix;}
  
  // channel_fsi1t.med ===============================================
    bc_Neum[0]= fix_disp0;bc_flag[0]=0;
    bc_Neum[0] = _DS_parameter._map_DSgroup[bc_gam];
//  if (bc_gam == 11) {bc_Neum[0]= fix_disp0;} //liquid inlet
//  if (bc_gam == 1000 && xp[1] < 0. ) {bc_Neum[0]= fix_disp0;} //solid inlet
//  if (bc_gam == 21) {bc_Neum[0]= fix_disp0;} //inlet
//  if (bc_gam == 20) {bc_Neum[0]=free_disp;}        // right free displacement
//  if (bc_gam == 15) {bc_Neum[0]=fix_disp0;}        //left fix_disp0
//  if (bc_gam == 13) {bc_Neum[0]=fix_disp0;}     //liquid outlet
//  if (bc_gam == 1000 && xp[1] > 0. ) {bc_Neum[0]=fix_disp0;} //outlet
//  if (bc_gam == 23) {bc_Neum[0]=fix_disp0;}        //solid fix_disp0
 if (xp[1]>0.2499) {bc_Neum[0]=fix_disp0;} 
 if (xp[1]<-0.2499) {bc_Neum[0]=fix_disp0;} 
 // end  channel_fsi1t.med=============================================== 
  
//   // T_FSI.med ===============================================
// if (bc_gam == 11) {bc_Neum[0]=fix_disp0;}  //bottom
// if (bc_gam == 12) {bc_Neum[0]=fix_disp0;}     //upper walls
// if (bc_gam == 13) {bc_Neum[0]=fix_disp0;}     //left wall
// if (bc_gam == 14) {bc_Neum[0]=fix_disp0;}  // outlet 
// if (bc_gam == 15) {bc_Neum[0]=fix_disp0;}   //right wall
// if (bc_gam == 1000) {bc_Neum[0]=free_disp;}  //interface
// // end  T-FSI.med===============================================

#endif  
// ============================================================================
 #if DIMENSION==3   // ===================3D ==================================
 
    bc_Neum[0]= fix_disp0;bc_flag[0]=0;
    bc_Neum[0] = _DS_parameter._map_DSgroup[bc_gam];
 
 
 #endif    // ================================================================
// ============================================================================
  return;
 }
#endif
