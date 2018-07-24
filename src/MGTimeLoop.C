// ==================================================================
//                  Class MGTimeLoop
// ==================================================================
// lib include
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>

// config include
// #include "FilePaths_conf.h"
#include "Printinfo_conf.h"
#include "Equations_tab.h"

// local include
#include "MGTimeLoop.h"

// class include
#include "MGUtils.h"
#include "EquationSystemsExtendedM.h"
#include "MGMesh.h"

// =============================================================================
/// Constructor
MGTimeLoop::MGTimeLoop(
  MGUtils & mgutils_in,  
  EquationSystemsExtendedM & mgeqmap_in
):
  _mgutils(mgutils_in),_mgeqmap(mgeqmap_in)
{// ===========================================================================
  
}

// ============================================================================
///This function does the  start or the restart
void MGTimeLoop::transient_setup(
    const int & restart,         ///< restart iteration (0=no restart)  (in)
    int& t_in,             ///< initial time iteration                  (in)
    double&  time                ///< running time                      (in)
)  {// ========================================================================

  if(restart != 0)  { _mgeqmap.read_soln(restart);}
  else   {
    _mgeqmap.print_soln(restart);  //print initial step sol.0.h5 and sol.0.xmf
    if(_mgeqmap._mgmesh._iproc == 0) _mgeqmap._mgmesh.write_c(restart);
  }
  // print --------------------------------------
  _mgeqmap.print_case(restart);   //print caseN.xmf&h5 = IC + BC flags
  _mgeqmap.print_case_mat_h5(restart);// print caseN.xmf&h5 = mat extended
   _mgeqmap.print_case_bc_h5(restart);// print caseN.xmf&h5 = mat extended
  this->transient_print_xmf(restart); //print timeN.xmf

  // restart -------------------------------------------
  if(restart) {
    double restartime = stod(_mgutils._sim_config["restartime"]);   // restart or not
    
    t_in=restart; time +=restartime; // set initial time
    std::cout << "\n *+*+* MGTimeLoop::transient_setup: RESTART  "
              << restart << " time= " << restartime << std::endl;
  }
  return;
}

  void MGTimeLoop::set_uooold(
    const int& flag,               ///< initial time iteration          (in)
    const double&  toll,                  ///< running time                    (in)
    const double&  time,                  ///< running time                    (in)
    const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  ) { // ===================================================================================
 _mgeqmap.set_uooold(flag, toll,time,eq_min,eq_max); // solve one step
 
 return;
  }



// ===========================================================================
/// This function controls the transient loop
void MGTimeLoop::steady(
  const int & nmax_step,  ///< number max of steps
  const double & toll,  ///< tolerance
    const int & i_step,               ///< running time iteration      (in)
    const int & print_step,           ///< print every                 (in)
  const  double    &  dt,                   ///< step time                   (in)  
   const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
   const int  &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
)  { // ===================================================================================

  // A Soving the system ******************************************************************
   std::cout<<"\n  ** Solving nonlinear iteration  step "<< i_step <<" ***"<< std::endl;
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t  start_time=std::clock();
#endif   // -------------------------------------------------------------------             
  _mgeqmap.eqnmap_steady_loop(nmax_step, toll,dt,eq_min,eq_max); // solve one step
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step ********************************************
  // print solution in xdmf format
  if(_mgeqmap._mgmesh._iproc==0){
//    if(((t_step-t_in)/print_step)*print_step == (t_step-t_in) && t_step-t_in>0) {
    _mgeqmap.print_soln(i_step);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(i_step); // print mesh.N.h5
//   }
  }
//   time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}

// ============================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_loop(
    const int& t_in,               ///< initial time iteration          (in)
    double&  time,                  ///< running time                    (in)
  const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)  { //========================================================================

  //  parameters
  int    n_steps = stoi(_mgutils._sim_config["nsteps"]);
  double      dt = stod(_mgutils._sim_config["dt"]);
  int print_step = stoi(_mgutils._sim_config["printstep"]);
  
  // transient loop
  for(int t_step= t_in; t_step< t_in + n_steps; t_step++) {
    transient_onestep(t_in,t_step,print_step,time,dt, eq_min,eq_max);  ///< step time  
  }   // end time loop
  return;
}


// ===========================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_onestep(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt,                   ///< step time                   (in) 
  const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)  { // ===================================================================================

  // A Soving the system ******************************************************************
   std::cout<<"\n  ** Solving time step "<<t_step<< ", time = "<<time<<" ***"<< std::endl;
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t  start_time=std::clock();
#endif   // -------------------------------------------------------------------             
  _mgeqmap.eqnmap_timestep_loop(time, t_step-t_in-1,eq_min,eq_max); // solve one step
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step ********************************************
  // print solution in xdmf format
  if(_mgeqmap._mgmesh._iproc==0){
  if(((t_step-t_in)/print_step)*print_step == (t_step-t_in) && t_step-t_in>0) {
    _mgeqmap.print_soln(t_step);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(t_step); // print mesh.N.h5
  }
  }
//   time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}

// ===========================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_solve_and_update(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt,                   ///< step time                   (in) 
  const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)  { // ===================================================================================

  // A Soving the system ******************************************************************
   std::cout<<"\n  ** Solving time step "<<t_step<< ", time = "<<time<<" ***"<< std::endl;
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t  start_time=std::clock();
#endif   // -------------------------------------------------------------------             
  _mgeqmap.eqnmap_timestep_loop_and_update(time, t_step-t_in-1,eq_min,eq_max); // solve one step
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step ********************************************
  // print solution in xdmf format
  if(_mgeqmap._mgmesh._iproc==0){
  if(((t_step-t_in)/print_step)*print_step == (t_step-t_in) && t_step-t_in>0) {
    _mgeqmap.print_soln(t_step);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(t_step); // print mesh.N.h5
  }
  }
//   time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}
// ======================================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_control_onestep(
    const int  &nmax_step,             ///< number max of steps         (in)
    const int  & it,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt,                   ///< step time                   (in)  
    const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)  { // =================================================================================

  // A Soving the system ****************************************************************
#if PRINT_TIME==1 // only for cpu time check --------------------------------------------
  std::clock_t  start_time=std::clock();
#endif   // ----------------------------------------------------------------------------- 
  std::cout<<"\n **Solving control step "<<it<<", time= "<<time<<"***"<<std::endl;
            
  _mgeqmap.eqnmap_timestep_loop_control(nmax_step,it, dt,eq_min,eq_max); // solve one step control

#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step *****************************************
  // print solution in xdmf format
//   if(((t_step-t_in)/print_step)*print_step == (t_step-t_in) && t_step-t_in>0) {
    _mgeqmap.print_soln(it);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(it); // print mesh.N.h5
//   }
//   time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}

// ============================================================================
/// Xdmf transient  print
 /// This function  prints a file time.****.xmf in RESU
  void  MGTimeLoop::transient_print_xmf(
    int t_init                          ///< initial time iteration     (in)
)  {// ========================================================================

  // time parameters
   
    const int ndigits   = stoi(_mgutils._sim_config["ndigits"]);
  
//   const int ndigits     = _mgutils.get_par("ndigits");
  const int print_step   =  (int)(stoi(_mgutils._sim_config["printstep"]) +0.5);
  const int nsteps       =  (int)(stoi(_mgutils._sim_config["nsteps"])+0.5);
// 	for(int imesh=0;imesh<NUM_MESH_MAIN;imesh++){
  // dir names
  std::string basetime = _mgutils.get_file("BASETIME");
  std::string  basesol = _mgutils.get_file("BASESOL");
  std::string aux_xdmf = _mgutils.get_file("AUX_XDMF");
  // file

  std::ostringstream Name;
  Name << _mgutils._inout_dir << basetime << "."
       << setw(ndigits) << setfill('0') <<  t_init << ".xmf";
  std::ofstream out(Name.str().c_str());

  int nprt=1;
  std::string gname[3];  gname[0]=basesol;

// #ifdef TWO_PHASE
//   nprt=3;        gname[1]="int";        gname[2]="ccf";
// #endif
  //   Mesh -----------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM "
      <<  "\"" << _mgutils._contrib_dir << aux_xdmf << "\"" << "[]>\n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  for(int kp=0; kp< nprt; kp++)    {
    out << "<Grid Name=\""<< gname[kp].c_str() <<"\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
    // time loop for grid sequence
    for(int it=t_init; it<=t_init+nsteps; it=it+print_step) /*if(it%print_step ==0)*/   {
        out << "<xi:include href=\""
            << basesol << "."
            << setw(ndigits) << setfill('0') << it <<  ".xmf" << "\""
            << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< kp+1 <<"])\" >\n";
        out << "<xi:fallback />\n";
        out << " </xi:include>\n";
      }
    out << "</Grid> \n";
  }
  // Grid Collection end
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();
//       }
//       else {
// 	        std::ostringstream Name;
//         Name << _mgutils._inout_dir << basetime << "."
//              << setw(ndigits) << setfill('0') <<  t_init << ".xmf";
//         std::ofstream out(Name.str().c_str());
//
//         int nprt=1;
//         std::string gname[3];  gname[0]=basesol;
//
// #ifdef TWO_PHASE
//         nprt=3;        gname[1]="int";        gname[2]="ccf";
// #endif
//         //   Mesh -----------------------------------
//         out << "<?xml version=\"1.0\" ?> \n";
//         out << "<!DOCTYPE Xdmf SYSTEM "
//         <<  "\"" << _mgutils._contrib_dir << aux_xdmf << "\"" << "[]>\n";
//         out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
//         out << "<Domain> \n";
//         for (int kp=0;kp< nprt; kp++)    {
//           out << "<Grid Name=\""<< gname[kp].c_str() <<"\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
//           // time loop for grid sequence
//           for (int it=t_init+1;it<=t_init+nsteps;it++) if (it%print_step ==0)   {
//               out << "<xi:include href=\""
//               << basesol<< "."
//               << setw(ndigits) << setfill('0') << it <<  ".xmf" << "\""
//               << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< kp+1 <<"])\" >\n";
//               out << "<xi:fallback />\n";
//               out << " </xi:include>\n";
//             }
//           out << "</Grid> \n";
//         }
//         // Grid Collection end
//         out << "</Domain> \n";
//         out << "</Xdmf> \n";
//         out.close();
//       }
  return;
}

// ===========================================================================
/// This function controls the transient loop
void MGTimeLoop::dummy_step(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt                   ///< step time                   (in)  
)  { // =====================================================================

  // A Soving the system ****************************************************
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t  start_time=std::clock();
#endif   // ------------------------------------------------------------------- 
  std::cout << "\n  ** Dumming time step " << t_step
            << ", time = " << time << " ***" << std::endl;
            
//   _mgeqmap.eqnmap_timestep_loop(time, t_step-t_in-1); // solve one step
            
            
//   if(stoi(_mgutils._sim_config["FluidStructure"])!=0)     _mgeqmap.movemesh();  

#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step ****************************
  // print solution in xdmf format
  if(_mgeqmap._mgmesh._iproc==0){
//   if(((t_step-t_in)/print_step)*print_step == (t_step-t_in) && t_step-t_in>0)
  {
    _mgeqmap.print_soln(t_step);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(t_step); // print mesh.N.h5
  }
  }
//   time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}

