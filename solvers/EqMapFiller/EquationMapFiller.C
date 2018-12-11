#include "EquationMapFiller.h" 
#include "FEMUS.h"
#include "EquationSystemsExtendedM.h"
#include "Equations_tab.h"
#include "Equations_conf.h"


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


void EquationMapFiller::FillEquationMap(FEMUS& FemusProblem, const std::vector<FIELDS>& pbName)
{
    
  EquationSystemsExtendedM & EqMap = FemusProblem.get_MGExtSystem();
  int n_equations=pbName.size();
  
  for(int iname=0; iname<n_equations; iname++){
     FIELDS problem = pbName[iname];
     
     if(pbName[iname]== NS_F || pbName[iname]==NSX_F || pbName[iname]==NSY_F || pbName[iname]==NSZ_F) 
        setNavierStokes(EqMap);
     if(pbName[iname]== T_F) 
        setTemperature(EqMap);
     if(pbName[iname]== K_F ||  pbName[iname]==EW_F) 
        setDynamicTurbulence(EqMap);
     if(pbName[iname]== KTT_F ||  pbName[iname]==EWTT_F) 
        setThermalTurbulence(EqMap);
     if(pbName[iname]== CO_F) 
        setColor(EqMap);
     if(pbName[iname]== FS_F) 
        setFluidStructure(EqMap);
     if(pbName[iname]== SDS_F) 
        setDisplacements(EqMap);
  }

    return;
}

void EquationMapFiller::setNavierStokes(EquationSystemsExtendedM & EqMap)
{
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


void EquationMapFiller::setTemperature(EquationSystemsExtendedM& EqMap)
{
#ifdef T_EQUATIONS
    _nvars[0]=0;  _nvars[1]=0; _nvars[2]=1;  // only Quadratic[2] approx
    EqMap.AddSolver<MGSolT> ("T", T_F, _nvars[0],_nvars[1],_nvars[2],"T");
#endif
    return;
}


void EquationMapFiller::setDynamicTurbulence(EquationSystemsExtendedM& EqMap)
{
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


void EquationMapFiller::setThermalTurbulence(EquationSystemsExtendedM& EqMap)
{
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


void EquationMapFiller::setColor(EquationSystemsExtendedM& EqMap)
{
// COLOR EQUATION
#ifdef COLOR_EQUATIONS
      _nvars[0]=0;  _nvars[1]=0; _nvars[2]=1;  // only Quadratic[2] approx
      EqMap.AddSolver<MGSolCOL> ("C", CO_F, _nvars[0],_nvars[1],_nvars[2],"c");    
      EqMap.AddSolver<MGSolCOL> ("CK", CO_F + 1, _nvars[0],_nvars[1],_nvars[2],"k");    
#endif
   return;
}


void EquationMapFiller::setFluidStructure(EquationSystemsExtendedM& EqMap)
{
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

void EquationMapFiller::setDisplacements(EquationSystemsExtendedM& EqMap)
{
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

