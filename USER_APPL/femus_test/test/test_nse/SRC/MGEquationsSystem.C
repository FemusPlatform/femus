// std libraries -----------------
#include <iomanip>
#include <sstream>

// class
#include "MGEquationsSystem.h"

// conf files ------------------------------------------
#include "Equations_conf.h"   //choose the EQNS to solve
#include "Equations_tab.h"   //choose the EQNS to solve
#include "Printinfo_conf.h"
#include "MGFE_conf.h"
#include "Domain_conf.h"

// local inlcudes ----------------------
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGMesh.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "numeric_vectorM.h"

// classes included in the map ---------


#ifdef  NS_EQUATIONS
#include "MGSolverNS.h"  // NS ================
#if NS_EQUATIONS==2
#include "MGSolverP.h"
#endif
#endif // -----------------------------------------
#ifdef  NSA_EQUATIONS
#include "MGSolverNSA.h"  // NS ================
#endif // -----------------------------------------
#ifdef  TBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTBK.h"
#endif // -----------------------------------------
#ifdef  TBKA_EQUATIONS  // Turbulence -------------
#include "MGSolverTBKA.h"
#endif // -----------------------------------------
#ifdef  T_EQUATIONS // Temperature ================
#include "MGSolverT.h"
#endif // -----------------------------------------
#ifdef  CTRL_EQUATIONS // Temperature ================
#include "MGSolverCTRL.h"
#endif // -----------------------------------------
#ifdef  T_ADJ_EQUATIONS // Temperature adjoint ================
#include "MGSolverT_ADJ.h"
#endif // -----------------------------------------
#ifdef  T_G_EQUATIONS // Temperature control ================
#include "MGSolverT_G.h"
#endif // -----------------------------------------
#ifdef  T_COUP_EQUATIONS // Temperature control ================
#include "MGSolverT_COUP.h"
#endif // -----------------------------------------
#ifdef  ALFA_EQUATIONS
#include "MGSolverALFA.h"  // two fluids ===========
#endif // -----------------------------------------
#ifdef TTBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTTBK.h"
#endif
#ifdef SM_EQUATIONS  // Structural Mechanics ======
#include "MGSolverSM.h"
#include "MGSolverDS.h"
#endif
#ifdef FSI_EQUATIONS  // FSI =======================
#include "MGSolverFSI.h"
#include "MGSolverDS.h"
#ifdef FSIC_EQUATIONS
#include "MGSolverFSIC.h"
#endif
#endif
#ifdef FSIA_EQUATIONS  // FSI adjoint =======================
#include "MGSolverFSIA.h"
#endif
#ifdef COLOR_EQUATIONS
#include "MGSolverCOL.h"
#endif

#include "MGSolverDA.h"
#ifdef   TWO_PHASE
#include "MGSolverCC.h"
#endif



// =====================================================
// Table external fields
// ----------------------------------------------------


// ============================================
/// This function set the various MGSystems
void MGEquationsSystem::init(const std::vector<FIELDS> & pbName)  {

  int nvars_in[3]; nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=0;
  int n_equations=pbName.size();
//  for(int iname=0;iname<n_equations;iname++){


#ifdef DA_EQUATIONS // 
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== DA_F) {
      //number of variables of [0] piecewise [1] linear and [2] quadratic order
      nvars_in[0]=1;  nvars_in[1]=3; nvars_in[2]=1;
      MGSolDA   *mgsDA=new MGSolDA(*this,nvars_in);  // class def
      set_eqs(mgsDA);                                // set class -> equation_map
    }
#endif   // end DA_EQUATION  (Test system)


// ====================================================================
#ifdef NS_EQUATIONS // ---------  Navier-Stokes -----------
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== NS_F || pbName[iname]==NSX_F ||
        pbName[iname]==NSY_F || pbName[iname]==NSZ_F) {

      nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION);          // quadratic(2) approx
      nvars_in[0]=0;   nvars_in[1]=NS_EQUATIONS%2; // linear(1) approx
      if(NDOF_K>0) {   nvars_in[0]=1;   nvars_in[1]=0;}  // konstant(0) approx

// if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
#if NS_EQUATIONS==2     // - NS_EQUATIONS==2 -
      MGSolNS* mgs= new MGSolNS(*this,nvars_in,"NS0X","u");  set_eqs(mgs);
      MGSolNS* mgs2=new MGSolNS(*this,nvars_in,"NS0Y","v");  set_eqs(mgs2);
#if DIMENSION==3
      MGSolNS* mgs3=new MGSolNS(*this,nvars_in,"NS0Z","w");  set_eqs(mgs3);
#endif
#endif
#if (NS_EQUATIONS==0 || NS_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
      MGSolNS* mgs=new MGSolNS(*this,nvars_in,"NS0");  set_eqs(mgs);
//        MGSolFSI* mgs=new MGSolFSI(*this,nvars_in,"FSI0");
      set_num_eqs(mgs->_eqname,NS_F);
#endif

#if (NS_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=0;// only  Linear(1) approx
      MGSolP   *mgsP=new MGSolP(*this,nvars_in); 
      set_eqs(mgsP);
       set_num_eqs(mgsP->_eqname,P_F);
#endif
    }
#endif // -------------  end  Navier-Stokes ---------------

// ====================================================================
#ifdef NSA_EQUATIONS // ---------  Navier-Stokes Adjoint -----------
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== FS_F) {
      nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION);          // quadratic(2) approx
      nvars_in[0]=0;   nvars_in[1]=NS_EQUATIONS%2; // linear(1) approx
      if(NDOF_K>0) {   nvars_in[0]=1;   nvars_in[1]=0;}  // konstant(0) approx
//       nvars_in[1]=1; nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
// if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
#if NSA_EQUATIONS==2     // - NS_EQUATIONS==2 -
      MGSolNSA* mgs= new MGSolNSA(*this,nvars_in,"NSA0X","ua");   set_eqs(mgs);
      MGSolNSA* mgs2=new MGSolNSA(*this,nvars_in,"NSA0Y","va");   set_eqs(mgs2);
#if DIMENSION==3
      MGSolNSA* mgs3=new MGSolNSA(*this,nvars_in,"NSA0Z","wa");   set_eqs(mgs3);
#endif
#endif
#if (NSA_EQUATIONS==0 || NSA_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -

      MGSolNSA* mgs=new MGSolNSA(*this,nvars_in);  set_eqs(mgs);
#endif

    }
#endif // -------------  end  Navier-Stokes Adjoint ---------------



#ifdef ALFA_EQUATIONS
  nvars_in[0]=0; nvars_in[1]=1; nvars_in[2]=0;
  MGSolALFA   *mgsALFA=new MGSolALFA(*this,nvars_in, "ALFA", "a");     // class def
  set_eqs(mgsALFA);                                 // set class -> equation_map
#endif

// ====================================================================
// Energy-Equation
// ====================================================================
#ifdef T_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== T_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;  // only Quadratic[2] approx
      MGSolT   *mgsT=new MGSolT(*this,nvars_in); set_eqs(mgsT);     // class def
      set_num_eqs(mgsT->_eqname,T_F);
    }
#endif
#ifdef CTRL_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== CTR_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;  // only Quadratic[2] approx
      MGSolCTRL   *mgsCTRL=new MGSolCTRL(*this,nvars_in); set_eqs(mgsCTRL);     // class def
      set_num_eqs(mgsCTRL->_eqname,CTR_F);
    }
#endif
// ====================================================================
// Energy-Equation Boundary control
// ====================================================================

#ifdef T_ADJ_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== TA_F) {

      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;
      MGSolT_ADJ   *mgsT_adj=new MGSolT_ADJ(*this,nvars_in); set_eqs(mgsT_adj);     // class def
    }
#endif
#ifdef T_G_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== FS_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;
      MGSolT_G   *mgsT_g=new MGSolT_G(*this,nvars_in); set_eqs(mgsT_g);     // class def
    }
#endif

#ifdef T_COUP_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== TA_F) {
#if T_COUP_EQUATIONS == 1
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=2;
#else
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=1;
#endif
      MGSolT_COUP   *mgsT_coup=new MGSolT_COUP(*this,nvars_in); set_eqs(mgsT_coup);     // class def
    }
#endif

// ====================================================================
#ifdef SM_EQUATIONS // ---------  STRUCTURAL MECHANICS -----------
// ====================================================================
  for(int iname=0; iname<n_equations; iname++)   if(pbName[iname]== SM_F) {
      nvars_in[0]=0; nvars_in[1]=SM_EQUATIONS%2;   // Costant(1)  Linear(0)
      nvars_in[2]=((SM_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation


#if (SM_EQUATIONS==2)     // - SM_EQUATIONS==2 --------------------------
      MGSolSM* mgsm= new MGSolSM(*this,nvars_in,"SM0X","u");   set_eqs(mgms);
      MGSolSM* mgsm2=new MGSolSM(*this,nvars_in,"SM0Y","v");   set_eqs(mgsm2);
#if (DIMENSION==3)
      MGSolSM* mgsm3=new MGSolSM(*this,nvars_in,"SM0Z","w"); set_eqs(mgsm3);
#endif
#else                   // - SM_EQUATIONS==0,1 -----------------------
      MGSolSM* mgsm=new MGSolSM(*this,nvars_in);  set_eqs(mgsm);
#endif
    }
#endif // =================  end SM ===================================
// ====================================================================

// ====================================================================
// ====================================================================
// ====================================================================
#ifdef FSI_EQUATIONS //    FLUID-STRUCTURE ---  FSI
// ====================================================================
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== FS_F) {
//   nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION);          // quadratic(2) approx
//   nvars_in[0]=0;   nvars_in[1]=NS_EQUATIONS%2; // linear(1) approx
//   if(NDOF_K>0){   nvars_in[0]=1;   nvars_in[1]=0;}   // konstant(0) approx
      nvars_in[0]=0; nvars_in[1]=FSI_EQUATIONS%2;   // Costant(1)  Linear(0)
      nvars_in[2]=((FSI_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
      if(NDOF_K>0) {   nvars_in[0]=1;   nvars_in[1]=0;}  // konstant(0) approx
#if FSI_EQUATIONS!=2     // - NS_EQUATIONS==2 -
      MGSolFSI* mgs=new MGSolFSI(*this,nvars_in,"FSI0");
      set_eqs(mgs);
      set_num_eqs(mgs->_eqname,FS_F);

#else
      MGSolFSI* mgs= new MGSolFSI(*this,nvars_in,"FSI0X","u");  set_eqs(mgs);
      MGSolFSI* mgs2=new MGSolFSI(*this,nvars_in,"FSI0Y","v");  set_eqs(mgs2);
#if DIMENSION==3
      MGSolFSI* mgs3=new MGSolFSI(*this,nvars_in,"FSI0Z","w");  set_eqs(mgs3);
#endif
#endif
#if (FSI_EQUATIONS%2==0) // - FSI_EQUATIONS==0,2  projection -
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=0;// only  Linear(1) approx
      MGSolFSIP   *mgsP=new MGSolFSIP(*this,nvars_in,"FSIP","p");  set_eqs(mgsP);
#endif
    }
#ifdef DS_EQUATIONS // -------------  Displacement --------------------
  nvars_in[0]=0; nvars_in[1]=0;   // Costant(1)  Linear(0)
  nvars_in[2]=1;                  // quadratic(2) Approximation

  MGSolDS* mgsdsdx=new MGSolDS(*this,nvars_in,"SDSX","dx"); 
  set_eqs(mgsdsdx);  set_num_eqs(mgsdsdx->_eqname,SDSX_F);
  MGSolDS* mgsdsdy=new MGSolDS(*this,nvars_in,"SDSY","dy");  
  set_eqs(mgsdsdy);  set_num_eqs(mgsdsdy->_eqname,SDSY_F);
#if (DIMENSION==3)
  MGSolDS* mgsdsdz=new MGSolDS(*this,nvars_in,"SDSZ","dz"); 
  set_eqs(mgsdsdz);  set_num_eqs(mgsdsdz->_eqname,SDSZ_F);
#endif


#endif // ----------------  end  Disp ---------------------------------  

//   #ifdef COLOR_EQUATIONS //   COLORS
// // ====================================================================
// //    if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
//   nvars_in[0]=1; nvars_in[1]=0;   // Costant(1)  Linear(0)
//   nvars_in[2]=1; // quadratic(2) Approximation
//
//   MGSolCOL* mgscol= new MGSolCOL(*this,nvars_in);   set_eqs(mgscol);
// #endif // =================  end COLORS ==================================
  int a=1;


  // ====================================================================
#ifdef FSIA_EQUATIONS //    FLUID-STRUCTURE adjoint---  FSI
// ====================================================================
//   std::cout << " pointer2 FSI "  << _num_equations["FSI0"];
//   std::cout << " pointer2 FSIA "  << _num_equations["FSIA0"];

  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== FSA_F) {

//       nvars_in[0]=0; nvars_in[1]=FSIA_EQUATIONS%2;   // Costant(1)  Linear(0)
//       nvars_in[2]=((FSIA_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
      nvars_in[0]=0; nvars_in[1]=FSIA_EQUATIONS%2;   // Costant(1)  Linear(0)
      nvars_in[2]=((FSIA_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
      if(NDOF_K>0) {   nvars_in[0]=1;   nvars_in[1]=0;}  // konstant(0) approx

#if FSIA_EQUATIONS!=2
      MGSolFSIA* mgsa=new MGSolFSIA(*this,nvars_in,"FSIA0"); set_eqs(mgsa); set_num_eqs(mgsa->_eqname,FSA_F);


#else
      MGSolFSIA* mgsa= new MGSolFSIA(*this,nvars_in,"FSIA0X","ua");  set_eqs(mgsa);
      MGSolFSIA* mgsa2=new MGSolFSIA(*this,nvars_in,"FSIA0Y","va");  set_eqs(mgsa2);
#if DIMENSION==3
      MGSolFSIA* mgsa3=new MGSolFSIA(*this,nvars_in,"FSIA0Z","wa");  set_eqs(mgsa3);
#endif
#endif
#if (FSIA_EQUATIONS%2==0) // - FSI_EQUATIONS==0,2  projection -
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=0;// only  Linear(1) approx
      MGSolFSIAP   *mgsPa=new MGSolFSIAP(*this,nvars_in,"FSIAP","pa");  set_eqs(mgsPa);
#endif
    }
// #ifdef DSA_EQUATIONS // -------------  Displacement adjoint--------------------
//   nvars_in[0]=0; nvars_in[1]=0;   // Costant(1)  Linear(0)
//   nvars_in[2]=1;                  // quadratic(2) Approximation
//
//   MGSolDSA* mgsdsdxa=new MGSolDSA(*this,nvars_in,"SDSAX","dxa");  set_eqs(mgsdsdxa);
//   MGSolDSA* mgsdsdya=new MGSolDSA(*this,nvars_in,"SDSAY","dya");  set_eqs(mgsdsdya);
// #if (DIMENSION==3)
//   MGSolDSA* mgsdsdza=new MGSolDSA(*this,nvars_in,"SDSAZ","dza");  set_eqs(mgsdsdza);
// #endif
#endif // ----------------  end  Disp adjoint ---------------------------------  
#endif // =================  end FSI adjoint ==================================
// ====================================================================

// #endif // =================  end FSI ==================================

#ifdef COLOR_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== CO_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;  // only Quadratic[2] approx
      MGSolCOL   *mgsCOL=new MGSolCOL(*this,nvars_in,"C","c"); set_eqs(mgsCOL);     // class def
      set_num_eqs(mgsCOL->_eqname,CO_F);
      MGSolCOL   *mgsKur=new MGSolCOL(*this,nvars_in,"CK","k"); set_eqs(mgsKur);     // class def
      set_num_eqs(mgsKur->_eqname,CO_F+1);
    }
#endif
//   std::cout << " pointer "  << _equations[FSI0];

//   double r=0.;


// ========================================================================================
// ==================  Turbulence NS   ====================================================
#ifdef TBK_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== K_F ||  pbName[iname]==EW_F) {
      nvars_in[0]=0;                     // Costant(0)
      nvars_in[1]=0;                     // Linear(1)
      nvars_in[2]=(TBK_EQUATIONS%2)+1;   // Quadratic(2)
// ---------------  Turbulence K model -> 1 ----------------------------------------------
#if ((TBK_EQUATIONS/2)==0)
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in); // class def
      set_eqs(mgsTBK);                                 // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++

#if ((TBK_EQUATIONS/2)==1)    // Turbulence K-E model -> 1 +++++++++++++++++++++++++++++++
#if (TBK_EQUATIONS%2==0)                        // splitting [2]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in,"K","kt");  // class def (kt)
      set_eqs(mgsTBK);                                           // set class -> equation_map
      set_num_eqs(mgsTBK->_eqname,K_F);
      MGSolTBK   *mgsTBW=new MGSolTBK(*this,nvars_in,"K2","et"); // class def (et)
      set_eqs(mgsTBW);                                           // set class -> equation_map
      set_num_eqs(mgsTBW->_eqname,K_F+1);
#else                                           // coupled [3]  (kt,et)
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in);           // class def  (kt,et)
      set_eqs(mgsTBK);                                           // set class -> equation_map
      set_num_eqs(mgsTBK->_eqname,K_F);
#endif
#endif               // end K-e model ++++++++++++++++++++++

#if ((TBK_EQUATIONS/2)==2)    // Turbulence KW model -> 1 +++++++++++++++++++++++++++++++ 
#if (TBK_EQUATIONS%2==0)                         // splitting  [4]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in,"K2K","kt");  // class def
      set_eqs(mgsTBK);                                           // set class -> equation_map
      set_num_eqs(mgsTBK->_eqname,K_F);
      MGSolTBK   *mgsTBW=new MGSolTBK(*this,nvars_in,"K1W","wt"); // class def
      set_eqs(mgsTBW);                                           // set class -> equation_map
      set_num_eqs(mgsTBW->_eqname,K_F+1);
#else                                            // coupled [5]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in);         // class def   (kt,et)
      set_eqs(mgsTBK);                                         // set class -> equation_map
       set_num_eqs(mgsTBK->_eqname,K_F);
#endif                                                     // end splitting-coupled 
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif    // end TBK model =============================================================

// ========================================================================================
// ==================  Turbulence adjoint NS   ====================================================
#ifdef TBKA_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== KTT_F) {
      nvars_in[0]=0;                     // Costant(0)
      nvars_in[1]=0;                     // Linear(1)
      nvars_in[2]=(TBKA_EQUATIONS%2)+1;   // Quadratic(2)
// ---------------  Turbulence K model -> 1 ----------------------------------------------
#if ((TBKA_EQUATIONS/2)==0)
      MGSolTBKA   *mgsTBKA=new MGSolTBKA(*this,nvars_in); // class def
      set_eqs(mgsTBKA);                                 // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++
#if ((TBKA_EQUATIONS/2)==2)    // Turbulence KW model -> 1// coupled [5] +++++++++++++++++++++++++++++++ 
      MGSolTBKA   *mgsTBKA=new MGSolTBKA(*this,nvars_in);         // class def   (kt,et)
      set_eqs(mgsTBKA);                                         // set class -> equation_map
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif    // end TBK adjoint model =============================================================



// ============================================================================================
// ================     Turbulence energy  ====================================================
#ifdef TTBK_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== KTT_F || pbName[iname]== EWTT_F) {
      nvars_in[0]=0;                     // Costant   -> 0
      nvars_in[1]=0;                     // Linear    -> 1
      nvars_in[2]=(TTBK_EQUATIONS%2)+1;  // Quadratic -> 2
// --------------------- Turbulence K equation (kh)  ------------------------------------------
#if ((TTBK_EQUATIONS/2)==0)    // Turbulence K model -> 1 +++++++++++++++++++++++++++++++++++++
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);          // class def   (kt)
      set_eqs(mgsTTBK);                                            // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++
// --------------------- Turbulence K-e equation (kh,eh)---------------------------------------
#if ((TTBK_EQUATIONS/2)==1)          // Turbulence K-E model -> 1 +++++++++++++++++++++++++++++
#if (TTBK_EQUATIONS%2==0)                             // splitting [2]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in,"TK","kh");  // class def   (kt,et)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
      MGSolTTBK   *mgsTTBW=new MGSolTTBK(*this,nvars_in,"TK2","eh"); // class def   (kt,et)
      set_eqs(mgsTTBW);                                              // set class -> equation_map
#else                                                  // coupled [3]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);            // class def   (kt,et)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
#endif
#endif                             // end K-e model +++++++++++++++++++++++++++++++++++++++++++
// -----------------------Turbulence K-w equation (kh,wh)--------------------------------------
#if ((TTBK_EQUATIONS/2)==2)    // Turbulence KW model -> 1 ++++++++++++++++++++++++++++++++++++
#if (TTBK_EQUATIONS%2==0)          // splitting [4]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in,"TK","kh");  // class def   (kh)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
      MGSolTTBK   *mgsTTBW=new MGSolTTBK(*this,nvars_in,"TK2","wh"); // class def   (wh)
      set_eqs(mgsTTBW);// set class -> equation_map
#else                                                 // coupled [5]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);             // class def   (kh,wh)
      set_eqs(mgsTTBK);                                               // set class -> equation_map
#endif                   // end splitting-coupled 
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif  // ========================end Turbulence energy======================================


  // ===================================  external ==============================================
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase* mgsol = eqn->second;// get the pointer
    mgsol -> set_ext_fields(pbName);          // init ext fields
  }
  // ================================================================================
  return;
}
#ifdef   TWO_PHASE
void  MGEquationsSystem::set_mgcc(MGSolCC & cc) {
   // Reading operators
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase* mgsol = eqn->second;// get the pointer
    mgsol -> set_mgcc(cc);          // set mgcc
  }
}
#endif


// =================================================================
// ==========================================================================================
/// This function performes all the MGSystem time step routines for control problem
void MGEquationsSystem::eqnmap_timestep_loop_control(
  const double time,             // real time
  const int delta_t_step_in     // integer time
) {


#ifdef  TBKA_EQUATIONS   //NS + Turbulence distributed control
  MGSolBase* mgsolns  =  get_eqs("NS0");
#if NS_EQUATIONS%2==0
  MGSolBase* mgsolnsp =  get_eqs("NSP");
#endif
  MGSolBase* mgsolk   =  get_eqs("K");
  MGSolBase* mgsolnsa =  get_eqs("NSA0");
  MGSolBase* mgsolka  =  get_eqs("KA");

  const int NoLevels=mgsolns->_NoLevels;

  mgsolnsa->x_old[NoLevels-1]->localize(*mgsolnsa ->x_oold[NoLevels-1]);

  const int max_iter=1000; const int max_iter_adj=1000; const int max_func_iter=50; const double toll_conv=1.e-6;

// //   mgsolns->reset_dt(); mgsolk->reset_dt();
  for(int isolns=0; isolns<max_iter; isolns++)  {
    (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->disp[NoLevels-1]);
    (mgsolk->x_old[NoLevels-1])->localize(*mgsolk ->disp[NoLevels-1]);
    mgsolns -> MGTimeStep(time,delta_t_step_in);
#if NS_EQUATIONS%2==0
    mgsolnsp -> MGTimeStep(time,delta_t_step_in);
#endif
    mgsolk -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_ns=fabs(((mgsolns->disp[NoLevels-1]->l2_norm())-(mgsolns ->x_old[NoLevels-1]->l2_norm()))/
                                  (mgsolns ->x_old[NoLevels-1]->l2_norm()));
    const double err_norm_k=fabs(((mgsolk->disp[NoLevels-1]->l2_norm())-(mgsolk ->x_old[NoLevels-1]->l2_norm()))/
                                 (mgsolk ->x_old[NoLevels-1]->l2_norm()));
    if(err_norm_ns< toll_conv && err_norm_k< toll_conv) {
      std::cout << "\nSteady state NS-K found, iteration " << isolns+1 << " \n \n";
      break;
    } else if(mgsolns->x_old[NoLevels-1]->l2_norm() > 1.e+10 || mgsolk->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
      std::cout << "\nSystem NS-K broken, aborting" << '\n';
      abort();
      break;
    }
// //      /*if ((isolns+1)%(2000)==0)*/ {mgsolns->set_dt(1.02); mgsolk->set_dt(1.02);}
    if(isolns % 100 == 0) { std::cout << "NS: "<< isolns << ' ' << err_norm_ns << " and " << err_norm_k << '\n'; }
    if(isolns == max_iter-1) { std::cout << "\nMaximum iteration reached in NS!!! Error norm is " << err_norm_ns << " and " << err_norm_k << '\n'; }
  }

  double func0= mgsolns->MGFunctional(0,0.); //compute functional with _eta   std::cout << "\n Functional is " << func0 << endl;

// // // //   mgsolnsa->reset_dt(); mgsolka->reset_dt();
  for(int isolad=0; isolad<max_iter_adj; isolad++)  {
    (mgsolnsa->x_old[NoLevels-1])->localize(*mgsolnsa ->disp[NoLevels-1]);
    (mgsolka->x_old[NoLevels-1])->localize(*mgsolka ->disp[NoLevels-1]);
    mgsolka -> MGTimeStep(time,delta_t_step_in);
    mgsolnsa -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_nsa=fabs(((mgsolnsa ->disp[NoLevels-1]->l2_norm())-(mgsolnsa ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolnsa ->x_old[NoLevels-1]->l2_norm()));
//     const double err_norm_nsa=1.;
    const double err_norm_ka=fabs(((mgsolka ->disp[NoLevels-1]->l2_norm())-(mgsolka ->x_old[NoLevels-1]->l2_norm()))/
                                  (mgsolka ->x_old[NoLevels-1]->l2_norm()));
    if(err_norm_nsa < toll_conv && err_norm_ka < toll_conv) {
      std::cout << "\nSteady state NSA-KA found , iteration " << isolad+1 << '\n';
      break;
    } else if(mgsolnsa->x_old[NoLevels-1]->l2_norm() > 1.e+10 || mgsolka->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
      std::cout << "\nSystem NSA-KA broken, aborting" << '\n';
      abort();
      break;
    }
//     if ((isolad+1)%(2000)==0) {mgsolnsa->set_dt(2.); mgsolka->set_dt(2.);}
    if(isolad % 100 == 0) { std::cout << "Adjoint: "<< isolad << " " << err_norm_nsa << " " << err_norm_ka << '\n'; }
    if(isolad == max_iter_adj-1) { std::cout << "\nMaximum iteration reached in NSA!!! Error norm is "   << err_norm_nsa << " " << err_norm_ka << '\n'; }
  }


  (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_nonl[NoLevels-1]); //store the good value of NS

  for(int isolfunc=0; isolfunc<max_func_iter; isolfunc++)  {
//     mgsolns->reset_dt(); mgsolk->reset_dt();
    for(int isolns=0; isolns<max_iter; isolns++)  {
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns->disp[NoLevels-1]);
      (mgsolk->x_old[NoLevels-1])->localize(*mgsolk->disp[NoLevels-1]);
      mgsolns -> MGTimeStep(time,delta_t_step_in);
#if NS_EQUATIONS%2==0
      mgsolnsp -> MGTimeStep(time,delta_t_step_in);
#endif
      mgsolk -> MGTimeStep(time,delta_t_step_in);
      const double err_norm_ns=fabs(((mgsolns->disp[NoLevels-1]->l2_norm())-(mgsolns ->x_old[NoLevels-1]->l2_norm()))/
                                    (mgsolns ->x_old[NoLevels-1]->l2_norm()));
      const double err_norm_k=fabs(((mgsolk->disp[NoLevels-1]->l2_norm())-(mgsolk ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolk ->x_old[NoLevels-1]->l2_norm()));
      if(err_norm_ns< 1.e-6 && err_norm_k< 1.e-6) {
        std::cout << "\n Steady state NS-K found, iteration " << isolns+1 << '\n';
        break;
      } else if(mgsolns->x_old[NoLevels-1]->l2_norm() > 1.e+10 || mgsolk->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
        std::cout << "\n System NS-K broken, let's try with smaller _eta!" << '\n';
        isolfunc--;
        break;
      }
//       if ((isolns+1)%(2000)==0) {mgsolns->set_dt(2.); mgsolk->set_dt(2.);}
      if(isolns % 100 == 0) { std::cout << "NS: "<< isolns << ' ' << err_norm_ns << " and " << err_norm_k << '\n'; }
      if(isolns == max_iter-1) { std::cout << "\n Maximum iteration reached in NS!!! Error norm is " << err_norm_ns << " and " << err_norm_k << '\n'; }
    }
    double func1= mgsolns->MGFunctional(0,0);
    std::cout <<"Iter : " << isolfunc+1<< " old functional is " << func0 << " and new " << func1 << "\n";
    if((fabs(func0-func1)/func0) < toll_conv) {
      std::cout << "\n Convergence of the optimal control problem reached for equal functionals!! \n";
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_nonl[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      mgsolns->MGFunctional(1,0.);   //set _eta=0;
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_nonl[NoLevels-1]);
      break;
    } else if(func0<func1) {
      mgsolns->MGFunctional(1,0.66);  // set _eta=0.5*_eta
      (mgsolns->x_nonl[NoLevels-1])->localize(*mgsolns ->x_old[NoLevels-1]); //reset the good value of NS
    } else if(func0>func1) {
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_nonl[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      mgsolns->MGFunctional(1,1.5);   //set _eta=2.*_eta;
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_nonl[NoLevels-1]);
      break;
    }
    if(isolfunc==max_func_iter-1) {
      std::cout << "\n Convergence of the optimal control problem reached for maximum iteration!! \n";
// // //       mgsolns->MGFunctional(1,0.);   //set _eta=0;
      (mgsolns->x_nonl[NoLevels-1])->localize(*mgsolns ->x_old[NoLevels-1]); //reset the good value of NS
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_nonl[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      break;
    }
  }

  std::cout << endl;
#endif // end ifdef TBKA_EQUATIONS NS + Turbulence distributed control
// -------------------------------------------------------------------------------------------



//   // loop for time steps
//   for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
//     MGSolBase* mgsol = eqn->second;
//     mgsol -> MGTimeStep(time,delta_t_step_in);
//   }

#ifdef  FSIA_EQUATIONS   //FSI distributed control

  // FSI equations
#if FSI_EQUATIONS != 2  // couple
  MGSolBase* mgsolfsi    =  get_eqs("FSI0");
  const int NoLevels=mgsolfsi->_NoLevels;
#else
  MGSolBase* mgsolfsix    =  get_eqs("FSI0X"); // uncoupled
  MGSolBase* mgsolfsiy    =  get_eqs("FSI0Y");
  const int NoLevels=mgsolfsix->_NoLevels;
#if DIMENSION==3
  MGSolBase* mgsolfsiz    =  get_eqs("FSI0Z");
#endif
#endif
#if FSI_EQUATIONS % 2 == 0
  MGSolBase* mgsolfsip    =  get_eqs("FSIP");
#endif

  // displacement equation
  MGSolBase* mgsoldsx    =  get_eqs("SDSX");
  MGSolBase* mgsoldsy    =  get_eqs("SDSY");
#if DIMENSION==3
  MGSolBase* mgsoldsz    =  get_eqs("SDSZ");
#endif

  // adjoint equation
  MGSolBase* mgsolfsia   =  get_eqs("FSIA0");
//   MGSolBase* mgsoldsxa   =  get_eqs("SDSAX");
//   MGSolBase* mgsoldsya   =  get_eqs("SDSAY");
//   MGSolBase* mgsolcol    =  get_eqs("CO_F");
#if DIMENSION==3
//    MGSolBase* mgsoldsza   =  get_eqs("SDSAZ");
#endif


//   mgsolcol-> MGTimeStep(time,delta_t_step_in);




  // define flags ==================================================================================
  const int max_iter=50000; const int max_iter_adj=50000; const int max_func_iter=30; const double toll_conv=1.e-12;
#define FSI_PROB
//   #define FSI_ADJ
//   #define FSI_FUNCTIONAL
  // =================================================================================================



// #if FSI_EQUATIONS != 2
// //   mgsolfsi ->MGFunctional(1,0.); //set eta to zero, so use just x_oold of dsxa
//      mgsolfsi ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsi ->MGFunctional(1,0.5); //set eta to 0.5
// #else
//   //   mgsolfsi ->MGFunctional(1,0.); //set eta to zero, so use just x_oold of dsxa
//      mgsolfsix ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsix ->MGFunctional(1,0.5); //set eta to 0.5
//       mgsolfsiy ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsiy ->MGFunctional(1,0.5); //set eta to 0.5
// #if DIMENSION==3
//       mgsolfsiz ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsiz ->MGFunctional(1,0.5); //set eta to 0.5
// #endif
// #endif

  // ================================================================================
  // FSI DIRECT PROBLEM
  // ================================================================================
#ifdef FSI_PROB
  for(int isolns=0; isolns<max_iter; isolns++)  {
#if FSI_EQUATIONS != 2
    (mgsolfsi->x_old[NoLevels-1])->localize(*mgsolfsi ->disp[NoLevels-1]);
    mgsolfsi -> MGTimeStep(time,delta_t_step_in);
#else
    (mgsolfsix->x_old[NoLevels-1])->localize(*mgsolfsix ->disp[NoLevels-1]);
    (mgsolfsiy->x_old[NoLevels-1])->localize(*mgsolfsiy ->disp[NoLevels-1]);
    mgsolfsix -> MGTimeStep(time,delta_t_step_in);
    mgsolfsiy -> MGTimeStep(time,delta_t_step_in);
#if DIMENSION==3
    (mgsolfsiz->x_old[NoLevels-1])->localize(*mgsolfsiz ->disp[NoLevels-1]);
    mgsolfsiz -> MGTimeStep(time,delta_t_step_in);
#endif
#endif
#if FSI_EQUATIONS % 2 == 0
    mgsolfsip -> MGTimeStep(time,delta_t_step_in);
#endif
    mgsoldsx -> MGTimeStep(time,delta_t_step_in);
    mgsoldsy -> MGTimeStep(time,delta_t_step_in);
#if DIMENSION==3
    mgsoldsz -> MGTimeStep(time,delta_t_step_in);
#endif
#if FSI_EQUATIONS != 2
    const double err_norm_fsi=fabs(((mgsolfsi->disp[NoLevels-1]->l2_norm())-(mgsolfsi ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolfsi ->x_old[NoLevels-1]->l2_norm()));
#else
    const double err_norm_fsi=fabs(((mgsolfsix->disp[NoLevels-1]->l2_norm())-(mgsolfsix ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolfsix ->x_old[NoLevels-1]->l2_norm()))+
                              fabs(((mgsolfsiy->disp[NoLevels-1]->l2_norm())-(mgsolfsiy ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolfsiy ->x_old[NoLevels-1]->l2_norm()))
#if DIMENSION==3
                              +fabs(((mgsolfsiz->disp[NoLevels-1]->l2_norm())-(mgsolfsiz ->x_old[NoLevels-1]->l2_norm()))/
                                    (mgsolfsiz ->x_old[NoLevels-1]->l2_norm()))
#endif
                              ;
#endif
    if(err_norm_fsi< toll_conv) {
      std::cout << "Steady state FSI found, iteration " << isolns+1 << "\n \n";
      break;
    }
//      else if (mgsolfsi->x_old[NoLevels-1]->l2_norm() > 1.e+10){
//        std::cout << "\nSystem FSI broken, aborting" << '\n';
//        abort();
//         break;
//      }
    if(isolns % 100 == 0) { std::cout << "FSI: "<< isolns << ' ' << err_norm_fsi << '\n'; }
    if(isolns == max_iter-1) { std::cout << "\nMaximum iteration reached in FSI!!! Error norm is " << err_norm_fsi << '\n'; }
  }
#endif // FSI_PROB

  // ================================================================================
  // END FSI DIRECT PROBLEM
  // ================================================================================
  double ctrl=0.;
#if FSI_EQUATIONS != 2
  const double func0= mgsolfsi->MGFunctional(0,ctrl); //compute functional
#else
  const double func0= mgsolfsix->MGFunctional(0,ctrl); //compute functional
#endif

//    std::cout << "\nOld Functional is " << func0 << endl<< endl;
//double func3= mgsolfsi->MGFunctional(0,0.); //compute functional
//   std::cout <<"\n new functional is " << func3;

#ifdef FSI_ADJ
  for(int isolad=0; isolad<max_iter_adj; isolad++)  {
    (mgsolfsia->x_old[NoLevels-1])->localize(*mgsolfsia ->disp[NoLevels-1]);
    mgsolfsia -> MGTimeStep(time,delta_t_step_in);
//     mgsoldsxa -> MGTimeStep(time,delta_t_step_in);
//     mgsoldsya -> MGTimeStep(time,delta_t_step_in);
//     #if DIMENSION==3
//      mgsoldsza -> MGTimeStep(time,delta_t_step_in);
//      #endif
    const double err_norm_fsia=fabs(((mgsolfsia ->disp[NoLevels-1]->l2_norm())-(mgsolfsia ->x_old[NoLevels-1]->l2_norm()))/
                                    (mgsolfsia ->x_old[NoLevels-1]->l2_norm()));
    if(err_norm_fsia < toll_conv) {
      std::cout << "Steady state FSIA found , iteration " << isolad+1 << " \n\n";
      break;
    } else if(mgsolfsia->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
      std::cout << "\nSystem FSIA broken, aborting" << '\n';
      abort();
      break;
    }
    if(isolad % 100 == 0) { std::cout << "Adjoint: "<< isolad << " " << err_norm_fsia << '\n'; }
    if(isolad == max_iter_adj-1) { std::cout << "\nMaximum iteration reached in FSIA!!! Error norm is "   << err_norm_fsia << '\n'; }
  }
#endif // FSI_ADJ
//


#ifdef FSI_FUNCTIONAL

#if FSI_EQUATIONS != 2
  (mgsolfsi->x_old[NoLevels-1])->localize(*mgsolfsi ->x_user[NoLevels-1]); //store the good value of NS
#else
  (mgsolfsix->x_old[NoLevels-1])->localize(*mgsolfsix ->x_user[NoLevels-1]);
  (mgsolfsiy->x_old[NoLevels-1])->localize(*mgsolfsiy ->x_user[NoLevels-1]);
#if DIMENSION==3
  (mgsolfsiz->x_old[NoLevels-1])->localize(*mgsolfsiz ->x_user[NoLevels-1]);
#endif
#endif

#if FSI_EQUATIONS % 2 == 0
  (mgsolfsip->x_old[NoLevels-1]) -> localize(*mgsolfsip ->x_user[NoLevels-1]);
#endif
  (mgsoldsx->x_old[NoLevels-1])->localize(*mgsoldsx ->x_user[NoLevels-1]); //store the good value of NS
  (mgsoldsy->x_old[NoLevels-1])->localize(*mgsoldsy ->x_user[NoLevels-1]); //store the good value of NS
#if DIMENSION==3
  (mgsoldsz->x_old[NoLevels-1])->localize(*mgsoldsz ->x_user[NoLevels-1]); //reset the good value of NS
#endif

// #if FSI_EQUATIONS != 2
// //   mgsolfsi ->MGFunctional(1,0.); //set eta to zero, so use just x_oold of dsxa
//      mgsolfsi ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsi ->MGFunctional(1,0.5); //set eta to 0.5
// #else
//   //   mgsolfsi ->MGFunctional(1,0.); //set eta to zero, so use just x_oold of dsxa
//      mgsolfsix ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsix ->MGFunctional(1,0.5); //set eta to 0.5
//       mgsolfsiy ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsiy ->MGFunctional(1,0.5); //set eta to 0.5
// #if DIMENSION==3
//       mgsolfsiz ->MGFunctional(2,1.); //set eta to 1
//       mgsolfsiz ->MGFunctional(1,0.5); //set eta to 0.5
// #endif
// #endif

  for(int isolfunc=0; isolfunc<max_func_iter; isolfunc++)  {
//     mgsolns->reset_dt(); mgsolk->reset_dt();
    for(int isolns=0; isolns<max_iter; isolns++)  {
#if FSI_EQUATIONS != 2
      (mgsolfsi->x_old[NoLevels-1])->localize(*mgsolfsi ->disp[NoLevels-1]);
      mgsolfsi -> MGTimeStep(time,delta_t_step_in);
#else
      (mgsolfsix->x_old[NoLevels-1])->localize(*mgsolfsix ->disp[NoLevels-1]);
      (mgsolfsiy->x_old[NoLevels-1])->localize(*mgsolfsiy ->disp[NoLevels-1]);
      mgsolfsix -> MGTimeStep(time,delta_t_step_in);
      mgsolfsiy -> MGTimeStep(time,delta_t_step_in);
#if DIMENSION==3
      (mgsolfsiz->x_old[NoLevels-1])->localize(*mgsolfsiz ->disp[NoLevels-1]);
      mgsolfsiz -> MGTimeStep(time,delta_t_step_in);
#endif
#endif
#if FSI_EQUATIONS % 2 == 0
      mgsolfsip -> MGTimeStep(time,delta_t_step_in);
#endif
      mgsoldsx -> MGTimeStep(time,delta_t_step_in);
      mgsoldsy -> MGTimeStep(time,delta_t_step_in);
#if DIMENSION==3
      mgsoldsz -> MGTimeStep(time,delta_t_step_in);
#endif
#if FSI_EQUATIONS != 2
      const double err_norm_fsi=fabs(((mgsolfsi->disp[NoLevels-1]->l2_norm())-(mgsolfsi ->x_old[NoLevels-1]->l2_norm()))/
                                     (mgsolfsi ->x_old[NoLevels-1]->l2_norm()));
      const double norm_fsi=mgsolfsi ->x_old[NoLevels-1]->l2_norm();
#else
      const double err_norm_fsi=fabs(((mgsolfsix->disp[NoLevels-1]->l2_norm())-(mgsolfsix ->x_old[NoLevels-1]->l2_norm()))/
                                     (mgsolfsix ->x_old[NoLevels-1]->l2_norm()))+
                                fabs(((mgsolfsiy->disp[NoLevels-1]->l2_norm())-(mgsolfsiy ->x_old[NoLevels-1]->l2_norm()))/
                                     (mgsolfsiy ->x_old[NoLevels-1]->l2_norm()))
#if DIMENSION==3
                                +fabs(((mgsolfsiz->disp[NoLevels-1]->l2_norm())-(mgsolfsiz ->x_old[NoLevels-1]->l2_norm()))/
                                      (mgsolfsiz ->x_old[NoLevels-1]->l2_norm()))
#endif
                                ;
      const double norm_fsi=mgsolfsix ->x_old[NoLevels-1]->l2_norm();
#endif
      if(err_norm_fsi< toll_conv) {
        std::cout << "Steady state FSI found, iteration " << isolns+1 << " \n";
        break;
      }
      if(norm_fsi > 1.e+10) {
        std::cout << "\nSystem FSI broken, retry with smaller eta, aborting" << '\n';
        abort();
      }
      if(isolns % 100 == 0) { std::cout << "FSI: "<< isolns << ' ' << err_norm_fsi << '\n'; }
      if(isolns == max_iter-1) { std::cout << "\nMaximum iteration reached in FSI!!! Error norm is " << err_norm_fsi << '\n'; }
    }
#if FSI_EQUATIONS != 2
    double func1= mgsolfsi->MGFunctional(0,0);
#else
    double func1= mgsolfsix->MGFunctional(0,0); //compute functional
#endif
    std::cout <<"Iter : " << isolfunc+1<< " old functional is " << func0 << " and new " << func1 << "\n";
    if((fabs(func0-func1)/func0) < toll_conv) {
      std::cout << "\n Convergence of the optimal control problem reached for equal functionals!! \n";
//       /* compute and store actual control */
// #if FSI_EQUATIONS !=2
//       const double actual_eta=mgsolfsi->MGFunctional(5,0.);
// #else
//       const double actual_eta=mgsolfsix->MGFunctional(5,0.);
// #endif
//       mgsoldsxa->MGFunctional(0,actual_eta);
//       mgsoldsya->MGFunctional(1,actual_eta);
//       #if DIMENSION==3
//       mgsoldsza -> MGFunctional(2,actual_eta);
//       #endif
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_oold[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_oold[NoLevels-1]);
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_old[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_old[NoLevels-1]);
//       #if DIMENSION==3
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza ->x_old[NoLevels-1]);
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza->x_oold[NoLevels-1]);
//        #endif
//       /* end - compute and store actual control */
      /* compute and store actual control */
#if FSI_EQUATIONS !=2
      const double actual_eta=mgsolfsi->MGFunctional(5,0.);
#else
      const double actual_eta=mgsolfsix->MGFunctional(5,0.);
#endif
      mgsolfsia->MGFunctional(0,actual_eta);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_oold[NoLevels-1]);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_old[NoLevels-1]);
      /* end - compute and store actual control */
      break;
    } else if(func0<func1) {
#if FSI_EQUATIONS != 2
      mgsolfsi ->MGFunctional(1,0.66); // multiply eta with 0.66
#else
      mgsolfsix ->MGFunctional(1,0.66); // multiply eta with 0.66
      mgsolfsiy ->MGFunctional(1,0.66); // multiply eta with 0.66
#if DIMENSION==3
      mgsolfsiz ->MGFunctional(1,0.66); // multiply eta with 0.66
#endif
#endif
#if FSI_EQUATIONS != 2
      (mgsolfsi->x_user[NoLevels-1])->localize(*mgsolfsi ->x_old[NoLevels-1]); //reset the good value of FSI
#else
      (mgsolfsix->x_user[NoLevels-1])->localize(*mgsolfsix ->x_old[NoLevels-1]);
      (mgsolfsiy->x_user[NoLevels-1])->localize(*mgsolfsiy ->x_old[NoLevels-1]);
#if DIMENSION==3
      (mgsolfsiz->x_user[NoLevels-1])->localize(*mgsolfsiz ->x_old[NoLevels-1]);
#endif
#endif
#if FSI_EQUATIONS % 2 == 0
      (mgsolfsip->x_user[NoLevels-1]) -> localize(*mgsolfsip ->x_old[NoLevels-1]);
#endif
      mgsoldsx -> MGUndo_disp();
      mgsoldsy -> MGUndo_disp();
#if DIMENSION==3
      mgsoldsz -> MGUndo_disp();
#endif
      (mgsoldsx->x_user[NoLevels-1])->localize(*mgsoldsx ->x_old[NoLevels-1]); //reset the good value of NS
      (mgsoldsy->x_user[NoLevels-1])->localize(*mgsoldsy ->x_old[NoLevels-1]); //reset the good value of NS
#if DIMENSION==3
      (mgsoldsz->x_user[NoLevels-1])->localize(*mgsoldsz ->x_old[NoLevels-1]); //reset the good value of NS
#endif
    } else if(func0>func1) {
//       /* compute and store actual control */
// #if FSI_EQUATIONS !=2
//       const double actual_eta=mgsolfsi->MGFunctional(5,0.);
// #else
//       const double actual_eta=mgsolfsix->MGFunctional(5,0.);
// #endif
//       mgsoldsxa->MGFunctional(0,actual_eta);
//       mgsoldsya->MGFunctional(1,actual_eta);
//       #if DIMENSION==3
//       mgsoldsza -> MGFunctional(2,actual_eta);
//       #endif
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_oold[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_oold[NoLevels-1]);
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_old[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_old[NoLevels-1]);
//       #if DIMENSION==3
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza ->x_old[NoLevels-1]);
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza->x_oold[NoLevels-1]);
//        #endif
//       /* end - compute and store actual control */
      /* compute and store actual control */
#if FSI_EQUATIONS !=2
      const double actual_eta=mgsolfsi->MGFunctional(5,0.);
#else
      const double actual_eta=mgsolfsix->MGFunctional(5,0.);
#endif
      mgsolfsia->MGFunctional(0,actual_eta);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_oold[NoLevels-1]);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_old[NoLevels-1]);
      /* end - compute and store actual control */
      break;
    }
    if(isolfunc==max_func_iter-1) {
      std::cout << "\n Convergence of the optimal control problem reached for maximum iteration!! \n";
//       /* compute and store actual control */
// #if FSI_EQUATIONS !=2
//       const double actual_eta=mgsolfsi->MGFunctional(5,0.);
// #else
//       const double actual_eta=mgsolfsix->MGFunctional(5,0.);
// #endif
//       mgsoldsxa->MGFunctional(0,actual_eta);
//       mgsoldsya->MGFunctional(1,actual_eta);
//       #if DIMENSION==3
//       mgsoldsza -> MGFunctional(2,actual_eta);
//       #endif
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_oold[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_oold[NoLevels-1]);
//       (mgsoldsxa->x_user[NoLevels-1])->localize(*mgsoldsxa ->x_old[NoLevels-1]);
//       (mgsoldsya->x_user[NoLevels-1])->localize(*mgsoldsya ->x_old[NoLevels-1]);
//       #if DIMENSION==3
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza ->x_old[NoLevels-1]);
//       (mgsoldsza->x_user[NoLevels-1])->localize(*mgsoldsza->x_oold[NoLevels-1]);
//        #endif
//       /* end - compute and store actual control */
      /* compute and store actual control */
#if FSI_EQUATIONS !=2
      const double actual_eta=mgsolfsi->MGFunctional(5,0.);
#else
      const double actual_eta=mgsolfsix->MGFunctional(5,0.);
#endif
      mgsolfsia->MGFunctional(0,actual_eta);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_oold[NoLevels-1]);
      (mgsolfsia->x_user[NoLevels-1])->localize(*mgsolfsia ->x_old[NoLevels-1]);
      /* end - compute and store actual control */
      break;
    }
  }
#endif // FSI_FUNCTIONAL

  std::cout << endl;
#endif // end ifdef FSIA_EQUATIONS
// -------------------------------------------------------------------------------------------




// -------------------------------------------------------------------------------------------
#ifdef  T_ADJ_EQUATIONS   //Temperature boundary control


  MGSolBase* mgsolt  =  get_eqs("T");
  MGSolBase* mgsolta =  get_eqs("T_ad");
  MGSolBase* mgsoltg =  get_eqs("T_g");


  mgsolt->MGFunctional(1,0.); // set _eta=0

  const int NoLevels=mgsolt->_NoLevels;

  const int max_iter=10000; const int max_iter_adj=10000; const int max_func_iter=100; const double toll_conv=1.e-7;

  (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]);

  for(int isolt=0; isolt<max_iter; isolt++)  {
    (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->disp[NoLevels-1]);
    mgsolt -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_t=fabs(((mgsolt->disp[NoLevels-1]->l2_norm())-(mgsolt ->x_old[NoLevels-1]->l2_norm()))/
                                 (mgsolt ->x_old[NoLevels-1]->l2_norm()));
    if(err_norm_t< toll_conv) {
      std::cout << "\nSteady state T found, iteration " << isolt+1 << " \n \n";
      break;
    } else if(mgsolt->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
      std::cout << "\nSystem T broken, aborting" << '\n';
      abort();
      break;
    }
    if(isolt % 100 == 0) { std::cout << "T: "<< isolt << ' ' << err_norm_t <<'\n'; }
    if(isolt == max_iter-1) { std::cout << "\nMaximum iteration reached in T!!! Error norm is " << err_norm_t << '\n'; }
  }

  double func0= mgsolt->MGFunctional(0,0.); //compute functional with _eta   std::cout << "\n Functional is " << func0 << endl;

  for(int isolad=0; isolad<max_iter_adj; isolad++)  {
    (mgsolta->x_old[NoLevels-1])->localize(*mgsolta ->disp[NoLevels-1]);
    (mgsoltg->x_old[NoLevels-1])->localize(*mgsoltg ->disp[NoLevels-1]);
    mgsolta -> MGTimeStep(time,delta_t_step_in);
    mgsoltg -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_ta=fabs(((mgsolta ->disp[NoLevels-1]->l2_norm())-(mgsolta ->x_old[NoLevels-1]->l2_norm()))/
                                  (mgsolta ->x_old[NoLevels-1]->l2_norm()));
    const double err_norm_tg=fabs(((mgsoltg ->disp[NoLevels-1]->l2_norm())-(mgsoltg ->x_old[NoLevels-1]->l2_norm()))/
                                  (mgsoltg ->x_old[NoLevels-1]->l2_norm()));
    if(err_norm_ta < toll_conv) {
      std::cout << "\nSteady state TA-TG found , iteration " << isolad+1 << "\n \n";
      break;
    } else if(mgsolta->x_old[NoLevels-1]->l2_norm() > 1.e+10 || mgsoltg->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
      std::cout << "\nSystem TA-TG broken, aborting" << '\n';
      abort();
      break;
    }
    if(isolad % 100 == 0) { std::cout << "Adjoint: "<< isolad << " " << err_norm_ta << " " << err_norm_tg << '\n'; }
    if(isolad == max_iter_adj-1) { std::cout << "\nMaximum iteration reached in TA-TG!!! Error norm is "   << err_norm_ta << " " << err_norm_tg << '\n'; }
  }

  (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]); //store the good value of NS

  mgsolt->MGFunctional(2,0.); // set _eta=1
  mgsolt->MGFunctional(1,.1); // set _eta

  for(int isolfunc=0; isolfunc<max_func_iter; isolfunc++)  {
    for(int isolt=0; isolt<max_iter; isolt++)  {
      (mgsolt->x_old[NoLevels-1])->localize(*mgsolt->disp[NoLevels-1]);
      mgsolt -> MGTimeStep(time,delta_t_step_in);
      const double err_norm_t=fabs(((mgsolt->disp[NoLevels-1]->l2_norm())-(mgsolt ->x_old[NoLevels-1]->l2_norm()))/
                                   (mgsolt ->x_old[NoLevels-1]->l2_norm()));
      if(err_norm_t < toll_conv) {
        std::cout << "\n Steady state T found, iteration " << isolt+1 << '\n';
        break;
      } else if(mgsolt->x_old[NoLevels-1]->l2_norm() > 1.e+10) {
        std::cout << "\n System T broken, let's try with smaller _eta!" << '\n';
        isolfunc--;
        break;
      }
      if(isolt % 100 == 0) { std::cout << "T: "<< isolt << ' ' << err_norm_t << '\n'; }
      if(isolt == max_iter-1) { std::cout << "\n Maximum iteration reached in T!!! Error norm is " << err_norm_t << '\n'; }
    }
    double func1= mgsolt->MGFunctional(0,0);
    std::cout <<"Iter : " << isolfunc+1<< " old functional is " << func0 << " and new " << func1 << "\n";
    if((fabs(func0-func1)/func0) < toll_conv) {
      std::cout << "\n Convergence of the optimal control problem reached for equal functionals!! \n";
      (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]); //save the good value of NS
      break;
    } else if(func0<func1) {
      mgsolt->MGFunctional(1,0.66);  // set _eta=0.66*_eta
      (mgsolt->x_oold[NoLevels-1])->localize(*mgsolt ->x_old[NoLevels-1]); //reset the good value of NS
    } else if(func0>func1) {
      (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]); //save the good value of NS
      break;
    }
    if(isolfunc==max_func_iter-1) {
      std::cout << "\n Convergence of the optimal control problem reached for maximum iteration!! \n";
      (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]); //save the good value of NS
      break;
    }
  }

  std::cout << endl;
#endif // end ifdef T_ADJ_EQUATIONS  Temperature boundary control



  return;
}

// =================================================================







 
