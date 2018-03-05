// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"
#include  "Solverlib_conf.h"  // Solver library options 
// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGMesh.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"
#include "FEMUS.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#endif

#define HAVE_CONTROL
// #define HAVE_BFGS
/// Set up
// =======================================
// Main program
// =======================================

int main(int argc, char **argv) {
  argc = argc ;
  argv=argv;  // no unused warning

  std::cout<<" ============ MGUtils ===================================== \n";
  std::cout<<" =========================================================== \n";
  // setting MGUtils -> parameters and file name ------------------------
  std::vector<MGUtils *> mgutils;
  std::string MeshPosAndName[NUM_MESH];

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    mgutils.push_back(new MGUtils(i_mesh+1));
    MeshPosAndName[i_mesh] = mgutils[i_mesh]->_mesh_dir + mgutils[i_mesh]->_interface_mesh;
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "<< MeshPosAndName[i_mesh] <<"\n ";
  }
  std::cout<<" ============ end loop mesh ================================ \n";
  std::cout<<" =========================================================== \n";

  int    n_steps = stoi(mgutils[0]->_sim_config["nsteps"]);
  double      dt = stod(mgutils[0]->_sim_config["dt"]);
  int print_step = stoi(mgutils[0]->_sim_config["printstep"]);
  int    itime_0 = stoi(mgutils[0]->_sim_config["itime"]);
  int    rstart_0= stoi(mgutils[0]->_sim_config["restart"]);
  double time    = 0.;
  double ctrl_alpha=stof(mgutils[0]->_sim_config["ctrl_alpha"]);
  double ctrl_beta=stof(mgutils[0]->_sim_config["ctrl_beta"]);
  int levels = mgutils[0]->_geometry["nolevels"];

  // CONSTRUCTION OF FINITE ELEMENT CLASS -----------------------------------
  MGFEMap *mgfemap;  mgfemap=new MGFEMap();
  MGFE *dfe_q;  dfe_q=new MGFE(2,ELTYPE); dfe_q->init_qua();
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE);  dfe_l->init_lin();
  MGFE *dfe_k;  dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();

  mgfemap->set_FE(dfe_q);    // quadratic fem
  mgfemap->set_FE(dfe_l);    // linear fem
  mgfemap->set_FE(dfe_k);    // piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();

  // system problem =========================================================
  // CONSTRUCTION OF FEMUS PROBLEM ------------------------------------------
  std::vector<FIELDS> myproblemP;
  mgutils[0]->FillFieldsVector(myproblemP);
  FEMUS P;

  // std::cout<<" --------------------- INIT PARAM ----------------------\n";
  P.init_param(*mgutils[0]);                      // init parameter
  // std::cout<<" --------------------- INIT FEM ----------------------\n";
  P.init_fem(*mggeomel,*mgfemap);                 // init fem
  // std::cout<<" --------------------- SET MESH ----------------------\n";
  P.setMesh();                                // set MGmesh
  // std::cout<<" --------------------- SET SYSTEM ----------------------\n";
  P.setSystem(myproblemP);                        // set system
  // std::cout<<" --------------------- SET CTRL DOMAIN ----------------------\n";
  int i_proc=P.get_proc();

  // INITIALIZATION OF EQUATIONS TO SOLVE -----------------------------------

  // TIME LOOP --------------------------------------------------------------
//      for ( int itime=itime_0+1; itime<= itime_0 + n_steps; itime++ ) {
//           P.solve_onestep ( itime_0,itime,print_step,time,dt,0,11 ); // solving state
//              P.solve_onestep ( itime_0,itime,print_step,time,dt,12,19 ); // solving adjoint
//                 P.solve_onestep ( itime_0,itime,print_step,time,dt,20,22 ); // solving control
//           time +=dt;
//      }

  // Control set up =================================================
  P.setCtrlDomain(0.2,0.4,-0.05,0.05,0,0);          // set control region
  double  ctrl_old =1.; double  ctrl=ctrl_old;//     // initial control and updated control
//   double  func_grad=0.; double  v_adj=0.;
  int     num_its=0;         //number of nonlinear iter
  double  func_old=0.;  double  func =0.;     // old functional -------------
  P.set_uooold(0,1.e-4,0.,0,23);   //0 x_old->x_ooold all fields
  func =P.System_functional(SDSX_F,1,ctrl);
  func +=P.System_functional(FS_F,1,ctrl);
  func_old=(rstart_0==0)?1.e+15:func; // ------------------------------------

  // std::cout<<" --------------------- SOLVE SETUP ----------------------\n";
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)
#ifndef HAVE_CONTROL  //--------------------no control-----------------------------------
  for(int it=rstart_0+1; it<=rstart_0+n_steps; it++) {
    P.solve_onestep(itime_0,it,print_step,time,dt,0,11);  // Solve P problem
//           P.solve_onestep(itime_0,it,print_step,time,dt,18,20);  // Solve P problem
    time+=dt;
  }   // end time loop
  func=P.System_functional(SDSX_F,1,ctrl);
  if(i_proc==0) {
    std::ofstream funct;
    funct.open("AndamentoFunzionale.dat", std::ios_base::app);
    funct <<  " " << func_old << "    " << func << endl;
    funct.close();
  }
#else  //----------------------end no control-------------------------------------
//  -------------------CONTROL------------------------------------------------------
int it_tot=0;
  for(int it=rstart_0+1; it<=rstart_0+n_steps; it++) {   //control loop
    it_tot=it;
    num_its++;
    time =dt*it;

    for(int itime=itime_0; itime<= itime_0 +10; itime++) {  // Solve steady state FSI problem (0-11 equations_tab.h) ---------
      P.solve_onestep(itime_0,it,print_step,time,dt,0,11);
      time +=dt;
    }   // end time loop -------------------------------------------------------------

    func=P.System_functional(SDSX_F,1,ctrl); //displacement functional contribution
    func+=P.System_functional(FS_F,1,ctrl); //regularization term contribution
    std::cout << " Old functional is -> "<< func_old << " and new is -> "<< func <<  " \n";

    //Equal functionals: optimal control reached ====================================================================
    if(fabs(((func-func_old)/func_old))<1.e-3) {
      std::cout<<" \033[38;5;46m ********* Convergence reached for equal functionals \033[0m"<< std::endl;   
      break;
    }//===============================================================================================================

    //Func decreasing: control step successful =======================================================================
    if(func <  func_old  && func>1.e-20) {
      if(i_proc==0) {    //Print functional in a file
        std::ofstream funct;
        funct.open("AndamentoFunzionale.dat", std::ios_base::app);
        funct<<it<<"\t"<<func_old<<"\t"<<func<<"\t"<</*func_old-func_grad<<"\t" <<func_grad<<*/"\t"<< ctrl <<"\t"/*<<beta<<"\t"<<eta*/<<endl;
        funct.close();
      }
      func_old=func;
#ifdef HAVE_BFGS
      P.set_uooold(8,1.e-4,0.,0,9);   //"old" force in disp_oold for BFGS method
      P.set_uooold(4,1.e-4,0.,20,23);   //"old" control equation in disp_old for BFGS method
#endif
      P.set_uooold(0,1.e-4,0.,0,11);    //store new FSI fields in uooold
      P.set_uooold(2,1.e-4,0.,9,11);   // store new SDSX fields (disp) in disp_old
      P.set_uooold(2,1.e-4,0.,0,9);   // store new control(force) fields (disp) in disp_old

       std::cout << " \033[38;5;46m   Func decreasing: "<< func << "  ctrl is: \033[0m" <<  ctrl_old << std::endl ;
      
      for(int itime=itime_0; itime<= itime_0 +5; itime++) {   // number of iterations to reach stationary state
        P.solve_onestep(itime_0,it,print_step,time,dt,13,15);    // solve new adjoint fields (13-15)
        time +=dt;
      }   // end time loop
      for(int itime=itime_0; itime<= itime_0 +0; itime++) {   // number of iterations to reach stationary state
        P.solve_onestep(itime_0,it,print_step,time,dt,20,23);    // solve new control equaiton
        time +=dt;
      }   // end time loop
      P.set_uooold(0,1.e-4,0.,13,23); //store new fields adjoint and control in uooold
      #ifdef HAVE_BFGS
        ctrl=ctrl_old;   
      #else  
        ctrl/=0.7*0.7*0.7*0.7;
      #endif
      P.System_functional(SDSX_F,1,ctrl); //passes ctrl value to DS solver
      P.System_functional(FS_F,1,ctrl);   //passes ctrl value to FSI solver
      // -------------------------------------------------------------------------
    }
    //Func increasing: control step unsuccesful ====================================================================
    else {
      P.set_uooold(1,1.e-4,0.,0,11);   // get back fields from uooold
      P.set_uooold(3,1.e-4,0.,9,11);   // get back fields (disp) from disp_old  SDSX
      it--;           // it not successful
      ctrl *=.7;      //reducing control -------------------------------------------------------------------
      P.System_functional(SDSX_F,1,ctrl);   //passes ctrl value to DS solver
      P.System_functional(FS_F,1,ctrl);     //passes ctrl value to FSI solver
      std::cout << "\033[38;5;196m  Func increasing: from "<< func_old << "  to: " <<  func << "  so new ctrl:  \033[0m" <<  ctrl << std::endl ;
      // --------------------------------------------------------------------------------------------------
    }
    //==============================================================================================================
    if(fabs(ctrl)<1.e-7)   break; 

  }   // end control iterative  loop
      P.set_uooold(7,1.e-4,0.,0,9);   // moves control fields from disp_old -->x_ooold
      P.set_uooold(1,1.e-4,0.,0,9);   // moves control fields from x_ooold -->x
      P.dummy_step(it_tot,it_tot+1,it_tot,ctrl,ctrl);
//   //  ========CONTROL END==============================================================
  if(i_proc==0) {
    std::ofstream funct;
    funct.open("AndamentoFunzionale.dat", std::ios_base::app);funct <<"#"<<num_its<<endl;funct.close();
  }

  
#endif


  // CLEANING MEMORY --------------------------------------------------------
  P.terminate();
  mgutils.clear();
  delete dfe_q;
  delete dfe_l;
  delete dfe_k;  // delete   fem
  delete mggeomel;
  delete mgfemap;
  return 0;
}
