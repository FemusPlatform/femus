// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"  // Solver library options 

#include "FEMUS.h"
#include "TurbUtils.h"
#include "MGUtils.h"



void PrintTime ( double begin, double actual, int timestep );

int main ( int argc, char ** argv )
{
    argc = argc ;
    argv = argv; // no unused warning

    TurbUtils TurbParameter = TurbUtils ( 0, 0 );

    std::vector<MGUtils *> mgutils;
    std::string MeshPosAndName[NUM_MESH];

    for ( int i_mesh = 0; i_mesh < NUM_MESH; i_mesh++ ) {
        mgutils.push_back ( new MGUtils ( i_mesh + 1, &TurbParameter ) );
        MeshPosAndName[i_mesh] = mgutils[i_mesh]->_mesh_dir + mgutils[i_mesh]->_interface_mesh;
        std::cout << " \n P mesh file " << i_mesh + 1 << "= " << MeshPosAndName[i_mesh] << "\n ";
    }

    int    n_steps = stoi ( mgutils[0]->_sim_config["nsteps"] );
    double      dt = stod ( mgutils[0]->_sim_config["dt"] );
    int print_step = stoi ( mgutils[0]->_sim_config["printstep"] );
    int    itime_0 = stoi ( mgutils[0]->_sim_config["itime"] );
    double time    = 0.;

    // system problem =========================================================
    FEMUS P ( *mgutils[0] ); //  parameter list <- mgutils[0]

    // INITIALIZATION OF EQUATIONS TO SOLVE -----------------------------------
    P.solve_setup ( itime_0, time );                // initial time loop (t=0)

    // calculation of wall dist - galerkin interpolation of piecewise field
    P.solve_onestep ( itime_0, itime_0 + 1, print_step, time, dt, DIST, DIST ); // solving P

    std::clock_t begin_time = std::clock();

    for ( int itime = itime_0 + 1; itime <= itime_0 + n_steps; itime++ ) {
        time += dt;
        P.solve_onestep ( itime_0, itime, print_step, time, dt, MU_T, ALPHA_T ); // solving P
        P.solve_onestep ( itime_0, itime, print_step, time, dt, 0, 4 ); // solving P
        P.solve_onestep ( itime_0, itime, print_step, time, dt, 5, 19 ); // solving P


        std::clock_t par_time = std::clock();
        PrintTime ( begin_time, par_time, itime );

    }

    // CLEANING MEMORY --------------------------------------------------------
    P.terminate();
    mgutils.clear();

    return 0;
}



void PrintTime ( double begin_time, double par_time, int itime )
{
    int clock_per_min = ( int ) 60 * CLOCKS_PER_SEC;

    double secondi = double ( par_time - begin_time ) / CLOCKS_PER_SEC;
    double minuti  = secondi / 60.;
    double ore     = minuti /
                     60.;

    std::cout << " Step n. "
              << itime;

    if ( minuti < 1. ) {
        std::cout << " Computational time until now " << secondi << " s  \n";
    } else if ( minuti < 60. && minuti > 0.9999 ) {
        std::cout << " Computational time until now " << ( int ) minuti << " m and " << secondi - ( int ) minuti * 60 << " s  \n";
    } else {
        std::cout << " Computational time until now " << ( int ) ore << " h and " << ( int ) ( minuti - 60 * ( int ) ore ) << " m and " << secondi - ( int ) minuti * 60 << " s  \n";
    }

}


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
