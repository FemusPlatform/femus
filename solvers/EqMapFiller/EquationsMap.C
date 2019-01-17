#include "EquationsMap.h"
#include "FEMUS.h"
#include "EquationSystemsExtendedM.h"
#include "Equations_conf.h"
#include <vector>
#include <MGUtils.h>


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
#ifdef  TA_EQUATIONS // Temperature adjoint ================
#include "MGSolverTA.h"
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
#ifdef IMMERSED_BOUNDARY
#include "MGSolverIB.h"
#endif

#include "MGSolverDA.h"
#ifdef   TWO_PHASE
#include "MGSolverCC.h"
#endif


// ====================================================================
/// This is a constructor  of  EquationsMap (with mgutils class )
/// from Simulation_Configuration.in file
/// For example:  EquationsMap eqmapclass(mgutils[0]); 
/// where  std::vector<MGUtils*> mgutils;
// ---------------------------------------------------------------------
EquationsMap::EquationsMap(
  MGUtils & mgutils_in  ///< mgutils class 
) {
  Fill_FIELD_map(); // fill _map_str2field  
  mgutils_in.FillFieldsVector(*this,_myproblemP); // fill _myproblemP
 
  return;
}

// ====================================================================
// This function fills the problem vector _myproblemP from a std::vector<std::string> pbname  
void EquationsMap::Fill_pbName(std::vector<std::string> pbname) {

  _myproblemP.clear();
  int pbsize=pbname.size();
  for(int i=0; i<pbsize ; i++) { FIELDS ff=  _map_str2field[pbname[i]]; _myproblemP.push_back(ff); }

  // print ---------------------------------------------------------------------
  std::cout << "\n ========= myproblemP vector ======================= \n\n";
  for(int ii=0; ii< pbsize; ii++)  std::cout << "myproblemP["<< ii << "]="<< _myproblemP[ii]<<"\n";
  std::cout << "\n =============================== \n";



  return;
}


// =============================================================================
/// This function fills the FEMUS problem with equations (listed in _pbName)
void EquationsMap::FillEquationMap(
  FEMUS& FemusProblem                ///< FemusProblem
//   const std::vector<FIELDS>& pbName   ///< list of equations (from Simulation_Configuration.in)
) { // ==========================================================================

  EquationSystemsExtendedM & EqMap = FemusProblem.get_MGExtSystem();
  int n_equations=_myproblemP.size();

  for(int iname=0; iname<n_equations; iname++) {
    FIELDS problem = _myproblemP[iname];

    if(_myproblemP[iname]== NS_F  || _myproblemP[iname]==NSX_F || _myproblemP[iname]==NSY_F || _myproblemP[iname]==NSZ_F)         setNavierStokes(EqMap);
    if(_myproblemP[iname]== FS_F)  setFluidStructure(EqMap);
    if(_myproblemP[iname]== SDS_F) setDisplacements(EqMap);
    if(_myproblemP[iname]== T_F)  setTemperature(EqMap);
    if(_myproblemP[iname]== K_F   ||  _myproblemP[iname]==EW_F)   setDynamicTurbulence(EqMap);
    if(_myproblemP[iname]== KTT_F ||  _myproblemP[iname]==EWTT_F) setThermalTurbulence(EqMap);
    if(_myproblemP[iname]== CO_F)  setColor(EqMap);
    if(_myproblemP[iname]== NSA_F  || _myproblemP[iname]==NSAX_F || _myproblemP[iname]==NSAY_F || _myproblemP[iname]==NSAZ_F)         setAdjointNavierStokes(EqMap);
    if(_myproblemP[iname]== FSA_F) setAdjointFSI(EqMap);
    if(_myproblemP[iname]== TA_F)  setAdjointTemperature(EqMap);
    if(_myproblemP[iname]== KA_F   ||  _myproblemP[iname]==EWA_F)   setAdjointDynamicTurbulence(EqMap);
    if(_myproblemP[iname]== CTRL_F) setControlTemperature(EqMap);
  }

  return;
}


// ================================== NS_EQUATIONS ================================
void EquationsMap::setNavierStokes(EquationSystemsExtendedM & EqMap) {
#ifdef NS_EQUATIONS
  _nvars[2]=((NS_EQUATIONS==2)?1:DIMENSION);          // quadratic(2) approx
  _nvars[0]=0;   _nvars[1]=NS_EQUATIONS%2; // linear(1) approx
  if(NDOF_K>0) {   _nvars[0]=1;   _nvars[1]=0;}  // konstant(0) approx

#if NS_EQUATIONS==2     // - NS_EQUATIONS==2 -
  EqMap.AddSolver<MGSolNS> ("NS0X", NS_F, _nvars[0],_nvars[1],_nvars[2],"u");
  EqMap.AddSolver<MGSolNS> ("NS0Y", NS_F, _nvars[0],_nvars[1],_nvars[2],"v");
#if DIMENSION==3
  EqMap.AddSolver<MGSolNS> ("NS0Z", NS_F, _nvars[0],_nvars[1],_nvars[2],"w");
#endif
#endif
#if (NS_EQUATIONS==0 || NS_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
  EqMap.AddSolver<MGSolNS> ("NS0", NS_F, _nvars[0],_nvars[1],_nvars[2],"u");
#endif

#if (NS_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
  _nvars[0]=0;  _nvars[1]=1; _nvars[2]=0;// only  Linear(1) approx
  EqMap.AddSolver<MGSolP> ("NS2P", NS_F, _nvars[0],_nvars[1],_nvars[2],"p");
#endif
#endif
  return;
}
// ================================== FSI_EQUATIONS ================================
void EquationsMap::setFluidStructure(EquationSystemsExtendedM& EqMap) {
#ifdef FSI_EQUATIONS
  _nvars[0]=0; _nvars[1]=FSI_EQUATIONS%2;
  _nvars[2]=((FSI_EQUATIONS==2)?1:DIMENSION);
  if(NDOF_K>0) {   _nvars[0]=1;   _nvars[1]=0;}
#if FSI_EQUATIONS!=2
  EqMap.AddSolver<MGSolFSI> ("FSI0", FS_F, _nvars[0],_nvars[1],_nvars[2],"u");
#else
  EqMap.AddSolver<MGSolFSI> ("FSI0X", FS_F, _nvars[0],_nvars[1],_nvars[2],"u");
  EqMap.AddSolver<MGSolFSI> ("FSI0Y", FS_F + 1, _nvars[0],_nvars[1],_nvars[2],"v");
#if DIMENSION==3
  EqMap.AddSolver<MGSolFSI> ("FSI0Z", FS_F + 2, _nvars[0],_nvars[1],_nvars[2],"z");
#endif
#endif
#if (FSI_EQUATIONS%2==0)
  _nvars[0]=0;  _nvars[1]=1; _nvars[2]=0;
  EqMap.AddSolver<MGSolFSIP> ("FSIP", FS_F + 3, _nvars[0],_nvars[1],_nvars[2],"p");
#endif

// DISPLACEMENTS
  setDisplacements(EqMap);
#endif

  return;
}
// ================================== DS_EQUATIONS ================================
void EquationsMap::setDisplacements(EquationSystemsExtendedM& EqMap) {
#ifdef DS_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0;  _nvars[2]=1;
  EqMap.AddSolver<MGSolDS> ("SDSX", SDSX_F , _nvars[0],_nvars[1],_nvars[2],"dx");
  EqMap.AddSolver<MGSolDS> ("SDSY", SDSY_F , _nvars[0],_nvars[1],_nvars[2],"dy");
#if (DIMENSION==3)
  EqMap.AddSolver<MGSolDS> ("SDSZ", SDSZ_F , _nvars[0],_nvars[1],_nvars[2],"dz");
#endif
#endif

  return;
}

// ================================== T_EQUATIONS ================================
void EquationsMap::setTemperature(EquationSystemsExtendedM& EqMap) {
#ifdef T_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0; _nvars[2]=1;  // only Quadratic[2] approx
  EqMap.AddSolver<MGSolT> ("T", T_F, _nvars[0],_nvars[1],_nvars[2],"T");
#endif
  return;
}

// ================================== TBK_EQUATIONS ================================
void EquationsMap::setDynamicTurbulence(EquationSystemsExtendedM& EqMap) {
#ifdef TBK_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0;  _nvars[2]=(TBK_EQUATIONS%2)+1;

#if ((TBK_EQUATIONS/2)==1)  // k-epsilon      
  EqMap.AddSolver<MGSolTBK> ("K", K_F, _nvars[0],_nvars[1],_nvars[2],"kt");
  EqMap.AddSolver<MGSolTBK> ("K2", K_F + 1, _nvars[0],_nvars[1],_nvars[2],"et");
#endif

#if ((TBK_EQUATIONS/2)==2)  // k-omega  
  EqMap.AddSolver<MGSolTBK> ("K2K", K_F, _nvars[0],_nvars[1],_nvars[2],"kt");
  EqMap.AddSolver<MGSolTBK> ("K1W", K_F + 1, _nvars[0],_nvars[1],_nvars[2],"wt");
#endif
#endif
  return;
}

// ================================== TTBK_EQUATIONS ================================
void EquationsMap::setThermalTurbulence(EquationSystemsExtendedM& EqMap) {
#ifdef TTBK_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0;  _nvars[2]=(TTBK_EQUATIONS%2)+1;

#if ((TTBK_EQUATIONS/2)==1)  // k-epsilon      
  EqMap.AddSolver<MGSolTTBK> ("TK", KTT_F, _nvars[0],_nvars[1],_nvars[2],"kh");
  EqMap.AddSolver<MGSolTTBK> ("TK2", KTT_F + 1, _nvars[0],_nvars[1],_nvars[2],"eh");
#endif

#if ((TTBK_EQUATIONS/2)==2)  // k-omega  
  EqMap.AddSolver<MGSolTTBK> ("TK", KTT_F, _nvars[0],_nvars[1],_nvars[2],"kh");
  EqMap.AddSolver<MGSolTTBK> ("TK2", KTT_F + 1, _nvars[0],_nvars[1],_nvars[2],"wh");
#endif
#endif
  return;
}

// ================================== COLOR_EQUATIONS ================================
void EquationsMap::setColor(EquationSystemsExtendedM& EqMap) {
// COLOR EQUATION
#ifdef COLOR_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0; _nvars[2]=1;  // only Quadratic[2] approx
  EqMap.AddSolver<MGSolCOL> ("C", CO_F, _nvars[0],_nvars[1],_nvars[2],"c");
  EqMap.AddSolver<MGSolCOL> ("CK", CO_F + 1, _nvars[0],_nvars[1],_nvars[2],"k");
#endif
  return;
}


void EquationsMap::setAdjointFSI(EquationSystemsExtendedM& EqMap) {
#ifdef FSIA_EQUATIONS //    FLUID-STRUCTURE adjoint---  
// ====================================================================

      _nvars[0]=0; _nvars[1]=FSIA_EQUATIONS%2;   // Costant(1)  Linear(0)
      _nvars[2]=((FSIA_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
      if(NDOF_K>0) {   _nvars[0]=1;   _nvars[1]=0;}  // konstant(0) approx

#if FSIA_EQUATIONS!=2
      EqMap.AddSolver<MGSolFSIA> ("FSIA0", FSA_F, _nvars[0],_nvars[1],_nvars[2],"ua");
#else
      EqMap.AddSolver<MGSolFSIA> ("FSIA0X", FSA_F, _nvars[0],_nvars[1],_nvars[2],"ua");
      EqMap.AddSolver<MGSolFSIA> ("FSIA0Y", FSA_F + 1, _nvars[0],_nvars[1],_nvars[2],"va");
#if DIMENSION==3
      EqMap.AddSolver<MGSolFSIA> ("FSIA0Z", FSA_F + 2, _nvars[0],_nvars[1],_nvars[2],"wa");
#endif
#endif
#if (FSIA_EQUATIONS%2==0) // - FSI_EQUATIONS==0,2  projection -
      _nvars[0]=0;  _nvars[1]=1; _nvars[2]=0;// only  Linear(1) approx
      EqMap.AddSolver<MGSolFSIAP> ("FSIAP", FSA_F + 3, _nvars[0],_nvars[1],_nvars[2],"pa");
#endif
    
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
      return;
}
// ================================== NSA_EQUATIONS ================================
void EquationsMap::setAdjointNavierStokes(EquationSystemsExtendedM & EqMap) {
#ifdef NSA_EQUATIONS
  _nvars[2]=((NSA_EQUATIONS==2)?1:DIMENSION);          // quadratic(2) approx
  _nvars[0]=0;   _nvars[1]=NSA_EQUATIONS%2; // linear(1) approx
  if(NDOF_K>0) {   _nvars[0]=1;   _nvars[1]=0;}  // konstant(0) approx

#if NSA_EQUATIONS==2     // - NS_EQUATIONS==2 -
  EqMap.AddSolver<MGSolNSA> ("NSA0X", NSA_F, _nvars[0],_nvars[1],_nvars[2],"ua");
  EqMap.AddSolver<MGSolNSA> ("NSA0Y", NSA_F, _nvars[0],_nvars[1],_nvars[2],"va");
#if DIMENSION==3
  EqMap.AddSolver<MGSolNSA> ("NSA0Z", NSA_F, _nvars[0],_nvars[1],_nvars[2],"wa");
#endif
#endif
#if (NSA_EQUATIONS==0 || NSA_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
  EqMap.AddSolver<MGSolNSA> ("NSA0", NSA_F, _nvars[0],_nvars[1],_nvars[2],"ua");
#endif

#if (NSA_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
  _nvars[0]=0;  _nvars[1]=1; _nvars[2]=0;// only  Linear(1) approx
  EqMap.AddSolver<MGSolPA> ("NSA2P", NSA_F, _nvars[0],_nvars[1],_nvars[2],"pa");
#endif
#endif
  return;
}
// ================================== TA_EQUATIONS ================================
void EquationsMap::setAdjointTemperature(EquationSystemsExtendedM& EqMap) {
#ifdef TA_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0; _nvars[2]=1;  // only Quadratic[2] approx
  EqMap.AddSolver<MGSolTA> ("TA", TA_F, _nvars[0],_nvars[1],_nvars[2],"TA");
#endif
  return;
}
// ================================== TBKA_EQUATIONS ================================
void EquationsMap::setAdjointDynamicTurbulence(EquationSystemsExtendedM& EqMap) {
#ifdef TBKA_EQUATIONS
  _nvars[0]=0;  _nvars[1]=0;  _nvars[2]=(TBKA_EQUATIONS%2)+1;

#if ((TBKA_EQUATIONS/2)==1)  // k-epsilon      
  EqMap.AddSolver<MGSolTBKA> ("KA", KA_F, _nvars[0],_nvars[1],_nvars[2],"ka");
  EqMap.AddSolver<MGSolTBKA> ("K2A", KA_F + 1, _nvars[0],_nvars[1],_nvars[2],"ea");
#endif

#if ((TBKA_EQUATIONS/2)==2)  // k-omega  
  EqMap.AddSolver<MGSolTBKA> ("K2KA", KA_F, _nvars[0],_nvars[1],_nvars[2],"ka");
  EqMap.AddSolver<MGSolTBKA> ("K1WA", KA_F + 1, _nvars[0],_nvars[1],_nvars[2],"wa");
#endif
#endif
  return;
}
// ================================== CTRL_EQUATIONS ================================
void EquationsMap::setControlTemperature(EquationSystemsExtendedM& EqMap) {
#ifdef CTRL_EQUATIONS
  _nvars[0]=0; _nvars[1]=0;   // Costant(1)  Linear(0)
  _nvars[2]=1;                // quadratic(2) Approximation

  EqMap.AddSolver<MGSolCTRL> ("CTRLX", CTRLX_F, _nvars[0],_nvars[1],_nvars[2],"cx");
  #if CTRL_EQUATIONS!=1 
  EqMap.AddSolver<MGSolCTRL> ("CTRLY", CTRLX_F + 1, _nvars[0],_nvars[1],_nvars[2],"cy");
  
   #if (DIMENSION==3)
  EqMap.AddSolver<MGSolCTRL> ("CTRLZ", CTRLX_F + 2, _nvars[0],_nvars[1],_nvars[2],"cz");
   #endif
    
  #endif
#endif
    return;
}

