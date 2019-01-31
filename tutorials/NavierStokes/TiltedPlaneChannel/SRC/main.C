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
#include "InterfaceProjection.h"

#include "EquationsMap.h"
// #include "MGSolverNS.h"


int main ( int argc, char** argv ) {
    argc = argc ;
    argv=argv;  // no unused warning


    std::vector<MGUtils*> mgutils;
    std::string MeshPosAndName[NUM_MESH];

    for ( int i_mesh=0; i_mesh< NUM_MESH; i_mesh++ ) {
        mgutils.push_back ( new MGUtils ( i_mesh+1 ) );
        MeshPosAndName[i_mesh] = mgutils[i_mesh]->_mesh_dir + mgutils[i_mesh]->_interface_mesh;
        std::cout<<" \n P mesh file "<< i_mesh+1 << "= "<< MeshPosAndName[i_mesh] <<"\n ";
    }

    int    n_steps = stoi ( mgutils[0]->_sim_config["nsteps"] );
    double      dt = stod ( mgutils[0]->_sim_config["dt"] );
    int print_step = stoi ( mgutils[0]->_sim_config["printstep"] );
    int    itime_0 = stoi ( mgutils[0]->_sim_config["itime"] );
    double time    = 0.;

    // CONSTRUCTION OF FEMUS PROBLEM ------------------------------------------
    FEMUS P(*mgutils[0]);  //  parameter list <- mgutils[0]

    // INITIALIZATION OF EQUATIONS TO SOLVE -----------------------------------
    P.solve_setup ( itime_0,time );                 // initial time loop (t=0)

    for ( int itime=itime_0 +1; itime<= itime_0 + n_steps; itime++ ) {
        time +=dt;
        P.solve_onestep ( itime_0,itime,print_step,time,dt,0,19 ); // solving P
    }

    // CLEANING MEMORY --------------------------------------------------------
    P.terminate();
    mgutils.clear();

    return 0;
}

