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
  double time    = 0.;

  int levels = mgutils[0]->_geometry["nolevels"];

  // CONSTRUCTION OF FINITE ELEMENT CLASS -----------------------------------
  MGFEMap *mgfemap;  mgfemap=new MGFEMap();  MGFE *dfe_q;
  dfe_q=new MGFE(2,ELTYPE);  dfe_q->init_qua();
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE);  dfe_l->init_lin();
  MGFE *dfe_k;  dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();

  mgfemap->set_FE(dfe_q);    // quadratic fem
  mgfemap->set_FE(dfe_l);    // linear fem
  mgfemap->set_FE(dfe_k);    // piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();

  
  
  
  
  // system problem =========================================================
  std::vector<FIELDS> myproblemP;
  mgutils[0]->FillFieldsVector(myproblemP);
  // CONSTRUCTION OF FEMUS PROBLEM ------------------------------------------
  FEMUS P;
  std::cout<<" --------------------- INIT PARAM ----------------------\n";
  P.init_param(*mgutils[0]);                      // init parameter
  std::cout<<" --------------------- INIT FEM ----------------------\n";
  P.init_fem(*mggeomel,*mgfemap);                 // init fem
  std::cout<<" --------------------- SET MESH ----------------------\n";
  P.setMesh();                                // set MGmesh
  std::cout<<" --------------------- SET SYSTEM ----------------------\n";
  P.setSystem(myproblemP);                        // set system
  // INITIALIZATION OF EQUATIONS TO SOLVE -----------------------------------
  std::cout<<" --------------------- SOLVE SETUP ----------------------\n";
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)
  
  // TIME LOOP --------------------------------------------------------------
  for(int itime=itime_0+1; itime<= itime_0 + n_steps; itime++) {
    P.solve_and_update(itime_0,itime,print_step,time,dt,0,2);  // solving u star
    P.solve_onestep(itime_0,itime,print_step,time,dt,3,3);     // solving pressure
    P.solve_and_update(itime_0,itime,print_step,time,dt,0,2);  // solving correction step
    time +=dt;
  }
  // CLEANING MEMORY --------------------------------------------------------
  P.terminate();
  mgutils.clear();
  delete dfe_q;  delete dfe_l;  delete dfe_k; // delete   fem
  delete mggeomel;
  delete mgfemap;
  return 0;
}
