// std libraries -----------------
#include <math.h>
#include <iomanip>
#include <sstream>

// class
#include "MGEquationsSystem.h"

// conf files ------------------------------------------
#include "Equations_conf.h"   //choose the EQNS to solve
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



// ====================================================
/// This function constructs all the MGSystems
MGEquationsSystem::MGEquationsSystem(
  MGUtils& mgutils_in,// MGUtils pointer
//   MGSystem& mgphys_in,// MGSystem pointer
  MGMesh& mgmesh_in,  // MGMesh pointer
  MGFEMap& mgfemap_in, // MGFEMap pointer
  int np_data,
  int ncell_data
):MGSystem(mgutils_in,mgmesh_in,np_data,ncell_data),
//   _mgutils(mgutils_in),
//   _mgphys(mgphys_in),
//   _mgmesh(mgmesh_in),
  _mgfemap(mgfemap_in)   {// ====================================
}


// ====================================================
MGEquationsSystem::~MGEquationsSystem() {
  clean(); //deallocates the map of equations

}

// ====================================================
/// This function destroys all the MGSystems
void MGEquationsSystem::clean() {
  /*
   for(MGEquationsSystem::iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
     delete eqn->second;
   } */
  _equations.clear();
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
// ====================================================
/// This sets dof initial and boundary conditions and sets the operators
void  MGEquationsSystem::setDofBcOpIc() {

  // Reading operators
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase* mgsol = eqn->second;// get the pointer
    mgsol -> MGDofBcOp();          // init dof, GenBc, ReadOperators
    mgsol -> GenIc();              // initial solution
  }
  return;
}



// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_steady_loop(
   const int & nmax_step,  ///< number max of steps
  const double & toll,  ///< tolerance
  const double delta_t_step_in,  //   (in)
  const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int     &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
) {
  // Loop for time steps
  int NoLevels=0;
  double norm_new=1.e-20; double norm_old=1.e-20; double diff_norm=1.e-20;
  double diff_norm_old=10000.;
 
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
   
    if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
      NoLevels=mgsol->_NoLevels;
      norm_old += mgsol ->x_old[NoLevels-1]->l2_norm();
//       mgsol->x_ooold[NoLevels-1]=mgsol->x_old[NoLevels-1];
    }
  }
  double time_step =delta_t_step_in;
  double time=0;

  for(int istep=1; istep<= nmax_step; istep++) {
    //;*(istep);
//     std::cout<<"\n  *** Solving steady iteration n ="<<istep<<" ** Time step= "<<time_step<< " ***"<< std::endl;


    // equation loop -----------------------------------------------------------------------
    norm_new=1.e-20;
    for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
      MGSolBase* mgsol = eqn->second;
//       if(_num_equations[eqn->first] <flag_state) {
         if(_num_equations[eqn->first] >=eq_min && _num_equations[eqn->first] <= eq_max ) {
         NoLevels=mgsol->_NoLevels;
        mgsol ->set_dt(time_step);
        mgsol -> MGTimeStep(time,delta_t_step_in);
        norm_new += mgsol ->x_old[NoLevels-1]->l2_norm();
      }
    }
    diff_norm=fabs(norm_old-norm_new);
    // ---------------------------------------------------------------------------------------
    std::cout<<"\n step "<<istep <<": old norm="<< norm_old<< "; new norm="<< norm_new<<
    "; err  ="<<  diff_norm/(norm_new)    << std::endl;
    if(diff_norm/norm_old < toll) {  // diff_norm/norm_old < toll
      std::cout<<"*** Steady state found on  n ="<<istep<<" ** Time step= "<<time_step<< " ***"<< std::endl;
      break;
    } else {            //   diff_norm/norm_old > toll
      std::cout<<"\n  *** Steady state NOT found: relative norm difference is " << diff_norm;
      if(diff_norm< diff_norm_old) {  // diff_norm< diff_norm_old -> step ok ----------------------
        time_step  *=1.25;   std::cout<<", increasing  step to  "<< time_step<<std::endl;
      } // --------------
      else {
        time_step  *=0.95; std::cout<<", reduce step  to  "<< time_step << std::endl;
      }
      diff_norm_old= diff_norm;  norm_old =norm_new;
//       for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
//         MGSolBase* mgsol = eqn->second;
// //         if(_num_equations[eqn->first] <flag_state) {
//            if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
//           mgsol->x_ooold[NoLevels-1]=mgsol->x_old[NoLevels-1];
//           norm_old += mgsol ->x_old[NoLevels-1]->l2_norm();
//           
//         }
//       }
//       std::cout<<"\n  step end  "<<   diff_norm_old << "  "<< diff_norm <<" "<< norm_old<< "  "<< norm_new<< std::endl;
      time += time_step;

    } // end else  diff_norm/norm_old < toll
  }
  return;
}


// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::set_uooold(
   const int & flag,  ///<  0 xold-> x_ooold   1 x_ooold-> xold
  const double & toll,  ///< tolerance
  const double delta_t_step_in,  //   (in)
  const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
  const int     &  eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
) {
  // loop for time steps
  int NoLevels=0;
  
      for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
        MGSolBase* mgsol = eqn->second;
          NoLevels=mgsol->_NoLevels;
           if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) {
             if (flag==1) mgsol->set_xooold2x();
              else  mgsol->set_vector(flag);
            }
        }
  return;
}


// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_timestep_loop(
  const double time,             // real time
  const int delta_t_step_in,     // integer time
    const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) {
  // loop for time steps
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
     if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) 
    mgsol -> MGTimeStep(time,delta_t_step_in);
  }
  return;
}

// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_timestep_loop_and_update(
  const double time,             // real time
  const int delta_t_step_in,     // integer time
  const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
  const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) {
  // loop for time steps
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
     if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) 
    mgsol -> MGTimeStep_no_up(time,delta_t_step_in);
  }
    for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
     if(_num_equations[eqn->first] >= eq_min && _num_equations[eqn->first] <= eq_max ) 
    mgsol -> MGUpdateStep();
  }
  return;
}
// ==========================================================================================
/// This function sets the controlled domain
void MGEquationsSystem::eqnmap_ctrl_domain(
     const double xMin,
     const double xMax,
     const double yMin,
     const double yMax,
     const double zMin,
     const double zMax
) {
  // loop for time steps
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
    mgsol->set_ctrl_dom(xMin,xMax,yMin,yMax,zMin,zMax);
  }
  return;
}


/// This function prints xdmf and hdf5 file
void MGEquationsSystem::print_soln(const int t_step // time step
                                  ) {

  const int    iproc =_mgmesh._iproc;
  if(iproc==0) {  // print only one processor

    print_soln_h5(t_step);                   // print sol h5
    int n_lines=0,n_cells=0;
// #ifdef TWO_PHASE  // --------- cc ----------------
//     print_h5CC(file,t_flag,n_lines,n_lines); // print CC
// #endif // ----------------- cc --------------------
    print_soln_xmf(t_step,n_lines,n_cells);  // print xdmf file
  }

  return;
}

// =================================================================
/// This function prints the attributes into the corresponding hdf5 file
void MGEquationsSystem::print_soln_h5(const int t_flag // time flag
                                     ) {

  const int NoLevels = (int)(_mgutils._geometry["nolevels"]);
  const int ndigits  = stoi(_mgutils._sim_config["ndigits"]);

  // file  ---------------------------------------------
  // file name
  std::ostringstream filename;
  filename << _mgutils._inout_dir << _mgutils.get_file("BASESOL") << "." << setw(ndigits) << setfill('0') << t_flag << ".h5";
  // open file for hf5 storage
  hid_t   file= H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  H5Fclose(file);

  // print all systems ---------------------------
  MGEquationsSystem::const_iterator pos=_equations.begin();
  MGEquationsSystem::const_iterator pos_e=_equations.end();
  for(; pos!=pos_e; pos++)    {
    MGSolBase *mgsol=pos->second;
    mgsol->print_u(filename.str(),NoLevels-1);
#ifdef CTRL_EQUATIONS
       mgsol->print_weight_med("mesh_sol.med",NoLevels-1);
//     mgsol->print_u_med("mesh_sol.med",NoLevels-1);    
#endif
  }
// printing cell system data (MGSystem-> printdata)
  for(int idata=0; idata<_n_data[0]+_n_data[1]; idata++) {
    std::ostringstream dir_name;
    dir_name << "DATA" << idata;
    print_data_view(filename.str(),idata,dir_name.str());
  }
  return;
}

// ===================================================================
/// It prints the attributes in  Xdmf format for one time step
void MGEquationsSystem::print_soln_xmf(const int t_step, int /*n_lines*/,int /*n_cells*/) {

  const int NoLevels = (int)(_mgutils._geometry["nolevels"]);
  const int ndigits  = stoi(_mgutils._sim_config["ndigits"]);

  //  Mesh ----------------------
  const MGMesh& mgmesh=_mgmesh;
  int n_nodes    = mgmesh._NoNodes[NoLevels-1];
  int n_elements    = mgmesh._NoElements[0][NoLevels-1];

  // get parameters
  std::string inout_dir  = _mgutils._inout_dir;
  std::string basesol    = _mgutils.get_file("BASESOL");
  std::string basemesh   = _mgutils.get_file("BASEMESH");
  std::string contrib_dir= _mgutils.get_file("CONTRIB_DIR");
//   std::string aux_xdmf   = _mgutils.get_file("AUX_XDMF");
  std::string connlin    = _mgutils.get_file("CONNLIN");
  std::string   basecase = _mgutils.get_file("BASECASE");
  const double dt = stod(_mgutils._sim_config["dt"]);

  // files
//   std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
//   conn_file /* <<femus_dir  <<  "/" << appl_dir << "/"  << myapp << "/" << input_dir */ << basemesh;
  std::ostringstream topol_file; // topology file file (mesh_conn_lin.h5)
  topol_file  << basemesh   << connlin <<".h5";
//   conn_file << ".h5";
  std::ostringstream coord_time_file; // connectivity file (mesh_conn_lin.h5)
  coord_time_file /*<<femus_dir  << "/"<< output_dir << "/"*/ << basemesh << "." << std::setw(ndigits) << std::setfill('0') << t_step << ".h5";

  std::ostringstream filename; //  solution file xmf
  filename << inout_dir  << basesol << "." << setw(ndigits) << setfill('0') << t_step;
  std::ostringstream casefilename; //  solution file xmf
  casefilename << basecase << "." << setw(ndigits) << setfill('0') << 0 <<".h5";

  std::ostringstream attr_file;//  solution file h5
  /* attr_file << filename.str()<< ".h5"; */
  attr_file << basesol << "." << setw(ndigits) << setfill('0') << t_step<< ".h5";
  filename << ".xmf";

  // solution file  xmf
  std::ofstream out(filename.str().c_str());
  std::cout << " Solution written to= " << filename.str().c_str() << std::endl;
//  ++++++++++++ Header ++++++++++++++
  out << "<?xml version=\"1.0\" ?> \n";
//   out << "<!DOCTYPE Xdmf SYSTEM ";
//   out <<  "\"" << aux_xdmf << "\" > \n";
  out << "<Xdmf> \n" << "<Domain> \n"<< "<Grid Name=\"Mesh\"> \n";
  // time
  const double restartime = (stoi(_mgutils._sim_config["restart"]) != 0) ?  stod(_mgutils._sim_config["restartime"]) : 0.0;
  const int t_in =stoi(_mgutils._sim_config["itime"]);
  out << "<Time Value =\"" << restartime + (t_step-t_in)*dt << "\" /> \n";
  // +++++ Topology ++++++++++
  _mgmesh.print_xmf_topology(out,topol_file.str(),NoLevels-1,0);
  // +++++++  Geometry +++++++++++++++++
  _mgmesh.print_xmf_geometry(out,coord_time_file/*conn_file*/.str(),NoLevels-1,0);
  // ++++  Attributes ++++++++++++
  MGEquationsSystem::const_iterator pos1   = _equations.begin();
  MGEquationsSystem::const_iterator pos1_e = _equations.end();
  for(; pos1!=pos1_e; pos1++)   {
    MGSolBase *mgsol=pos1->second;
    mgsol->print_xml_attrib(out,n_nodes,n_elements*NSUBDOM,attr_file.str());
  }
//   print_xml_attrib(out,n_elements,n_nodes,attr_file.str());
//   printining cell attributes
  if(_n_data[0]+_n_data[1]>0) { print_xml_attrib(out,n_elements,n_nodes,attr_file.str()); }

#ifdef FSI_EQUATIONS
  print_xml_mat(out,n_nodes,n_elements*NSUBDOM,casefilename.str());
#endif
  // #ifdef TWO_PHASE
//   // print of CC
//   print_xmfCC(out,t_step,n_lines,n_cells);
// #endif
  out << "</Grid>\n" << "</Domain> \n" << "</Xdmf> \n";
  out.close();
  return;
}

// ========================================================================
/// This function read the solution for all the system (restart)
void MGEquationsSystem::read_soln(const int t_step) {
// --------------------------------------------------------------
  const int ndigits  = stoi(_mgutils._sim_config["ndigits"]);
//   const int ndigits  = _mgutils.get_par("ndigits");
  const  int restart_lev_flag = stoi(_mgutils._sim_config["restart_lev"]);

  // open file -----------------------------
  std::ostringstream namefile;
  namefile << _mgutils._inout_dir << _mgutils.get_file("BASESOL") << "."
           << setw(ndigits) << setfill('0') << t_step << ".xmf";

#ifdef PRINT_INFO // --------  info ------------------ 
  std::cout << "\n MGEquationsSystem::read_soln: Reading time  from "
            << namefile.str().c_str();
#endif  // -------------------------------------------
  std::ifstream in ; in.open(namefile.str().c_str());  //associate the file stream with the name of the file
  if(!in.is_open()) { std::cout << " MGCase: restart .xmf file not found "  << std::endl; abort(); }

  // reading time from xmf file --------------
  std::string buf="";  while(buf != "<Time") { in >> buf; }
  in >> buf >> buf;  buf=buf.substr(2,buf.size()-3);
//create an istringstream from a string
  std::istringstream buffer(buf); double restart_time;  buffer >> restart_time;

  //add parameter to system
//   _mgutils.set_par("restartime",restart_time);
  _mgutils.set_sim_par("restartime",std::to_string(restart_time));
  // ---------------------------------------------------
  // reading data from  sol.N.h5
  // ---------------------------------------------------
  // file name -----------------------------------------
  namefile.str("");  //empty string
  namefile << _mgutils._inout_dir << _mgutils.get_file("BASESOL")
           << "." << setw(ndigits) << setfill('0') << t_step << ".h5";

#ifdef PRINT_INFO  // --------------- info ---------------
  std::cout << "\n MGEquationsSystem::read_soln: Reading from file "
            << namefile.str().c_str() << std::endl;
#endif // ---------------------------------------------
  // loop reading over the variables ---------------------
  for(MGEquationsSystem::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase *mgsol=eqn->second;
    mgsol->read_u(namefile.str(),restart_lev_flag);
  } //  loop --------------------------------------------------------

// #ifdef TWO_PHASE
//   readCC(t_step);
// //   _mgsys.get_mgcc()->read_fine ( t_init,_mgsys.get_mgcc()->_Nlev-1 );
// //   readCC ( t_step );
// #endif

  return;
}


// ========================================================================
/// This function print data from single class equation to the mesh system  on vect_data
void MGEquationsSystem::print_mesh_data(double vect_data[]) {
// --------------------------------------------------------------

  // loop reading/printing over the equation ---------------------
  int count=0;
//   eqn=_equations.begin();
  for(MGEquationsSystem::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    if(_mgutils.get_file("MESHNUMBER")=="mesh1") {
      if(count==1) {
        MGSolBase *mgsol=eqn->second;
        mgsol->print_ext_data(vect_data); // compute from a system/mesh to vect_data
      }
      count++;
    }
    if(_mgutils.get_file("MESHNUMBER")=="mesh2") {
      if(count==0) {
        MGSolBase *mgsol=eqn->second;
        mgsol->print_ext_data(vect_data); // compute from a system/mesh to vect_data
      }
      count++;
    }
  } //  loop --------------------------------------------------------

  return;
}
// ====================================================

// ====================================================
/// This function prints initial and boundary data in xdmf+hdf5 format
void MGEquationsSystem::print_case(const int t_init) {

  const int    iproc =_mgmesh._iproc;

  if(iproc==0) {  // print only one processor
    print_case_h5(t_init); // ic+bc print format h5

    int n_lines=0,n_cells=0; // for VOF
// #ifdef TWO_PHASE
//     print_h5CC(hid_t file,t_init,&n_lines,&n_cells);
// #endif
    print_case_xmf(t_init,n_lines,n_cells); // xml format
  }
  return;
}





// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void MGEquationsSystem::print_case_h5(const int t_init) {

  const int NoLevels = (int)(_mgutils._geometry["nolevels"]);
  const int ndigits  = stoi(_mgutils._sim_config["ndigits"]);
  std::string output_dir = _mgutils._inout_dir;
  std::string   basecase = _mgutils.get_file("BASECASE");
  //  Mesh ---- ---------------------------------------------
  const MGMesh& mgmesh = _mgmesh;

  // file ---------------------------------------
  std::ostringstream filename; // file name
  filename << output_dir << basecase << "." << setw(ndigits) << setfill('0') << t_init << ".h5";

  hid_t   file = H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  mgmesh.print_subdom_hf5(filename.str())  ; // PID (n processor)

 if(stoi(_mgutils._sim_config["DynamicalTurbulence"])!=0)  mgmesh.print_dist_hf5(filename.str(),mgmesh._dist,"DIST")  ; // PID (n processor)
 
  H5Fclose(file);

  // loop over all systems ---------------------------
  MGEquationsSystem::const_iterator pos   = _equations.begin();
  MGEquationsSystem::const_iterator pos_e = _equations.end();
  for(; pos!=pos_e; pos++) {
    MGSolBase *mgsol=pos->second;
    mgsol->print_u(filename.str(),NoLevels-1);   // initial solution
    mgsol->print_bc(filename.str(),NoLevels-1); // boundary condition
  }

  return;
}



// ====================================================================
/// It prints the Xdmf file to read the initial and boundary conditions
void MGEquationsSystem::print_case_xmf(const int t_init,const int /*n_lines*/,
                                       const int /*n_cells*/) {

  // file ----------------------------------------------
  std::string inout_dir  = _mgutils._inout_dir;
  std::string   basecase = _mgutils.get_file("BASECASE");
  std::string   basemesh = _mgutils.get_file("BASEMESH");
  std::string  contrib_dir = _mgutils.get_file("CONTRIB_DIR");
  std::string  aux_xdmf = _mgutils.get_file("AUX_XDMF");
  std::string  connlin = _mgutils.get_file("CONNLIN");
  const int NoLevels = (int)(_mgutils._geometry["nolevels"]);
  const int ndigits  = stoi(_mgutils._sim_config["ndigits"]);

  // files
  std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
  conn_file << basemesh;
  std::ostringstream topol_file; // topology file file (mesh_conn_lin.h5)
  topol_file << conn_file.str() << connlin <<".h5";
  conn_file << ".h5";
  std::ostringstream filename; //  solution file xmf
  filename <<  inout_dir  << basecase << "." << setw(ndigits) << setfill('0') <<  t_init<< ".xmf";

  //  Mesh ---- ---------------------------------------------
  const MGMesh&    mgmesh = _mgmesh;
  int n_nodes    = mgmesh._NoNodes[NoLevels-1];
  int n_elements = mgmesh._NoElements[0][NoLevels -1];
  std::string var_name[3];    std::string var_type[3];

  //   File xdmf -------------------------------
  std::ofstream out(filename.str().c_str());
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM ";
  out <<  "\"" <<  aux_xdmf << "\" \n";
//    out << " [ <!ENTITY HeavyData \"\"> ] ";
  out << ">\n";
  out << "<Xdmf> \n" << "<Domain> \n" << "<Grid Name=\"Mesh\"> \n";
  // +++++ Topology ++++++++++
  _mgmesh.print_xmf_topology(out,topol_file.str(),NoLevels-1,0);
  // +++++++  Geometry +++++++++++++++++
  _mgmesh.print_xmf_geometry(out,conn_file.str(),NoLevels-1,0);

  // ++++  Attributes for each system +++++++++++++++++++
  MGEquationsSystem::const_iterator pos1=_equations.begin();
  MGEquationsSystem::const_iterator pos1_e=_equations.end();
  for(; pos1!=pos1_e; pos1++)   {
    MGSolBase *mgsol=pos1->second;
    for(int ivar=0; ivar<mgsol->_nvars[1]+mgsol->_nvars[2]; ivar++)     {

      // Volume and boundary conditions
      var_name[0] = mgsol->_var_names[ivar];
      var_name[1] =var_name[0]+"bd";
      var_name[2] =var_name[0]+"vl";
      var_type[0] ="Float";
      var_type[1] ="Float";
      var_type[2] ="Float";
      for(int ibvar=0; ibvar<3; ibvar++) {
        out << "<Attribute Name=\""<< var_name[ibvar] <<"\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        out << "<DataItem  DataType=\""<< var_type[ibvar].c_str()
            << "\" Precision=\"8\" Dimensions=\"" << n_nodes << "  "
            << 1 << "\" Format=\"HDF\">  \n";
        out <</* femus_dir << "/" << output_dir <<*/ basecase << "."
            << setw(ndigits) << setfill('0') << t_init << ".h5"
            << ":" << var_name[ibvar].c_str() << "\n";
        out << "</DataItem>\n" << "</Attribute>\n";
      }
    }
    for(int ivar=0; ivar<mgsol->_nvars[0]; ivar++)     {
      // Volume and boundary conditions
      var_name[0] = mgsol->_var_names[mgsol->_nvars[1]+mgsol->_nvars[2]+ivar];
      var_name[1] =var_name[0]+"bd";
      var_name[2] =var_name[0]+"vl";
      var_type[0] ="Float";
      var_type[1] ="Int";
      var_type[2] ="Int";
      for(int ibvar=0; ibvar<1; ibvar++) { // only initial condition ->var_name[0]
        out << "<Attribute Name=\""<< var_name[ibvar] <<"\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        out << "<DataItem  DataType=\""<< var_type[ibvar].c_str()
            << "\" Precision=\"8\" Dimensions=\"" << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
            << 1 << "\" Format=\"HDF\">  \n";
        out << basecase << "."
            << setw(ndigits) << setfill('0') << t_init << ".h5"
            << ":" << var_name[ibvar].c_str() << "\n";
        out << "</DataItem>\n" << "</Attribute>\n";
      }
    }
  }

  // ----------------------------------------------------------
  out << "<Attribute Name=\"PID\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
      << 1 << "\" Format=\"HDF\">  \n";
  out << basecase << "."
      << setw(ndigits) << setfill('0')<< t_init << ".h5" << ":PID\n";
  out << "</DataItem>\n" << "</Attribute>\n";
 if(stoi(_mgutils._sim_config["DynamicalTurbulence"])!=0){
  // ----------------------------------------------------------
  out << "<Attribute Name=\"DIST\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
      << 1 << "\" Format=\"HDF\">  \n";
  out << basecase << "."
      << setw(ndigits) << setfill('0')<< t_init << ".h5" << ":DIST\n";
  out << "</DataItem>\n" << "</Attribute>\n";
 }
// #ifdef TWO_PHASE
//   // print of CC
//   print_xmfCC(out, t_init,n_lines,n_cells);
// #endif
  out << "</Grid>\n"<< "</Domain> \n" << "</Xdmf> \n";
  out.close();

  return;
}

// ============================================================================
 double  MGEquationsSystem::System_functional(
  const int  & ff,                 ///< initial time iteration
  double parameter,
  double    &  control                   ///< step time
){
   double value=0.;
    for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
    if(_num_equations[eqn->first] == ff) {
      const int NoLevels=mgsol->_NoLevels;
      value +=mgsol->MGFunctional(parameter,control);
      std::cout << "  functional " <<_num_equations[eqn->first] <<"\n -----";
    }
  }
   
  return value;
 }
// #ifdef TWO_PHASE
// // ==================================================
// ///readCC
// void MGCase::readCC(const int flag_print) {
//
//   std::ostringstream filename;
//   filename <<"output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
//   hid_t  file_id = H5Fopen(filename.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//
//   hid_t dtset = H5Dopen(file_id, "/cres", H5P_DEFAULT);
//   hid_t filespace = H5Dget_space(dtset);
//   hsize_t dims[2]; H5Sget_simple_extent_dims(filespace, dims, NULL);
//
//   double *ccf = (new double[dims[0]]); double *xf = (new double[4*dims[0]]);
//   double *yf = (new double[4*dims[0]]); double *zf = (new double[4*dims[0]]);
//
//   hid_t status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ccf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X1", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X2", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,yf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X3", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,zf);
//   H5Dclose(dtset);
//
//   _mgsys.get_mgcc()->read_fine(dims[0],ccf,xf,yf,zf);
//
//   delete []ccf; delete []xf; delete []yf; delete []zf;
//   return;
// }
//
//
// /// It prints the format Xdmf for one time step of CC
// void MGCase::print_xmfCC(std::ofstream &out, const int t_step,
//                          const int n_lines, const int n_cells) {
//   //  Mesh ---- ---------------------------------------------
//   const MGMesh &mgmesh=_mgsys.get_mesh();
// //  int n_nodes=mgmesh._NoNodes[_cslevel];
//   int n_elements=mgmesh._NoElements[_cslevel];
//   const int nvrt=4* (DIMENSION-1);
// #if DIMENSION==2
//   std::string btype="Polyline";
//   std::string mtype="Quadrilateral";
// #else
//   std::string btype="Triangle";
//   std::string mtype="Hexahedron";
// #endif
// // time parameters
//   const double dt=_mgsys.get_par("dt");
//
//   // ============
//   // coarse color function
//   // ============
//   out << "<Attribute Name=\"CC\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
//   out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements*nvrt << "  " << 1 <<
//       "\" Format=\"HDF\">  \n";
//   out << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":CC" << "\n";
//   out << "</DataItem>\n";
//   out << "</Attribute>\n";
//   // Grid Collection end
//   out << "</Grid> \n";
//
//   // ==================
//   //  Interface
//   // =================
//   out << "<Grid Name=\"Interf\"> \n";
//   out << "<Time Value =\"" << t_step*dt << "\" /> \n";
//   // +++++ Topology ++++++++++
//   out << "<Topology Type=\""<< btype << "\"   Dimensions=\""<<  n_lines <<    "\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_lines <<"  "<< DIMENSION <<
//       "\" Format=\"HDF\">  \n";
//   out  << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intconn" << "\n";
//   out << "</DataStructure> \n" << "</Topology> \n";
//   // +++++++  Geometry +++++++++++++++++
//   out << "<Geometry Type=\"X_Y_Z\"> \n";
//   for(int ix=1; ix<4; ix++) {
//     out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_lines*DIMENSION << "  " << 1 <<
//         "\" Format=\"HDF\">  \n";
//     out <<  "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intX" <<ix<< "\n";
//     out << "</DataStructure> \n";
//   }
//   out << " </Geometry>\n";
//   // Grid Collection end
//   out << "</Grid> \n";
//   // =========================
//   // cc on fine grid
//   // =========================
//   out << "<Grid Name=\"ccf\"> \n";
//   out << "<Time Value =\"" << t_step*dt << "\" /> \n";
//   // +++++ Topology ++++++++++
//   out << "<Topology Type=\""<< mtype << "\"   Dimensions=\""<<  n_cells <<    "\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_cells <<"  "<< nvrt <<
//       "\" Format=\"HDF\">  \n";
//   out  << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 <<":" << "conn" << "\n";
//   out << "</DataStructure> \n" << "</Topology> \n";
//   // +++++++  Geometry +++++++++++++++++
//   out << "<Geometry Type=\"X_Y_Z\"> \n";
//   for(int ix=1; ix<4; ix++) {
//     out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_cells*nvrt << "  " << 1 <<
//         "\" Format=\"HDF\">  \n";
//     out <<  "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":X" <<ix<<"\n";
//     out << "</DataStructure> \n";
//   }
//   out << " </Geometry>\n";
//   // ccf
//   out  << "<Attribute Name=\"ccf\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
//   out  << "<DataItem Dimensions=\""<< n_cells << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\"> \n ";
//   out << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "ccf" << "\n";
//   out << "</DataItem>\n" << "</Attribute>\n";
//
//   return;
// }
//
//
// // =========================================================
// void MGCase::print_h5CC(hid_t file,const int flag_print,
//                         int *n_l_out,int *n_c_out)  {
//   //  Mesh ---- ---------------------------------------------
//   const MGMesh &mgmesh=_mgsys.get_mesh();
// //  int n_nodes=mgmesh._NoNodes[_cslevel];
//   const int n_elements=mgmesh._NoElements[_cslevel];
//
//   int nvrt= (DIMENSION-1) *4;
//   // ==========
//   //  function cc
//   // ===========
//   double *ucc;  ucc=new double[n_elements*(DIMENSION-1) *4];
//   hsize_t dimsf[2];  dimsf[0] = n_elements* (DIMENSION-1) *4;  dimsf[1] = 1;
//   MGSolCC *mgcc=_mgsys.get_mgcc();
//   mgcc->Write_xmf(ucc);
//   std::string name="CC";  print_Dhdf5(file,name,dimsf,ucc);
//   delete []ucc;
//   // ==============
//   // interface
//   // =============
//   const int n_lines=mgcc->get_nel();
//   // line coordinates
//   double *xcoord;
//   xcoord=new double[n_lines*DIMENSION];
//   double *ycoord;
//   ycoord=new double[n_lines*DIMENSION];
//   double *zcoord;
//   zcoord=new double[n_lines*DIMENSION];
//   int *tcon;
//   tcon=new int[n_lines*DIMENSION];
//   mgcc->get_int(n_lines,tcon,xcoord,ycoord,zcoord);
//   // coordinate datasets --------------------
//   dimsf[0] = n_lines*DIMENSION;  dimsf[1] = 1;
//   std::string name="intX1";
//   print_Dhdf5(file,name,dimsf,xcoord);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file, "intX1",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, xcoord );
// //   H5Dclose ( dataset );
//   name="intX2";
//   print_Dhdf5(file,name,dimsf,ycoord);
// //   dataset = H5Dcreate ( file, "intX2",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, ycoord );
// //   H5Dclose ( dataset );
//   name="intX3";
//   print_Dhdf5(file,name,dimsf,zcoord);
// //   dataset = H5Dcreate ( file, "intX3",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, zcoord );
// //   H5Dclose ( dataset );
// //   H5Sclose ( dataspace );
//   // line connectivity ------------------------
//   dimsf[0] = n_lines;  dimsf[1] = DIMENSION;
//   name="intconn";    print_Ihdf5(file,name,dimsf,tcon);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file, "intconn", H5T_NATIVE_INT, dataspace,
// //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT,tcon );
// //   H5Dclose ( dataset );
//   // Close/release resources.
// //   H5Sclose ( dataspace );
//   delete []tcon;  delete []xcoord;
//   delete []ycoord;  delete []zcoord;
//
//   // ==========
//   //  ccf
//   // ==========
//   std::ostringstream filename;
//   filename << "output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
//   hid_t file1= H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
//
//   // ccf
//   const int n_cells=mgcc->get_nef();
//   // coordinates
//   xcoord=new double[n_cells*nvrt];
//   ycoord=new double[n_cells*nvrt];
//   zcoord=new double[n_cells*nvrt];
//   double *ccc = (new double[n_cells]) ;
//   double *cres = (new double[n_cells]) ;
//   tcon=new int[n_cells*nvrt];
//   mgcc->Writefine_xmf(ccc,cres,xcoord,ycoord,zcoord);
//   // coordinate datasets --------------
//   dimsf[0] = n_cells*nvrt;  dimsf[1] = 1;
//
//   name="X1";  print_Dhdf5(file1,name,dimsf,xcoord);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file1, "X1",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, xcoord );
// //   H5Dclose ( dataset );
//   name="X2";  print_Dhdf5(file1,name,dimsf,ycoord);
// //   dataset = H5Dcreate ( file1, "X2",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, ycoord );
// //   H5Dclose ( dataset );
//   name="X3";  print_Dhdf5(file1,name,dimsf,zcoord);
// //   dataset = H5Dcreate ( file1, "X3",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, zcoord );
// //   H5Dclose ( dataset );
// //   H5Sclose ( dataspace );
//   // ccf --------------------------------
//   dimsf[0] = n_cells;  dimsf[1] = 1;
//   name="ccf";  print_Dhdf5(file1,name,dimsf,ccc);
//   name="cres";  print_Dhdf5(file1,name,dimsf,cres);
//   // line connectivity ----------------
//   for(int it=0; it<n_cells*nvrt; it++) tcon[it]=it;
//   dimsf[0] = n_cells;  dimsf[1] = nvrt;
//   name="conn";  print_Dhdf5(file1,name,dimsf,tcon);
//
//   H5Fclose(file1);
//   delete []tcon;  delete []ccc;
//   delete []cres;  delete []xcoord;
//   delete []ycoord;  delete []zcoord;
//
//   *n_l_out = n_lines;  *n_c_out = n_cells;
//   return;
// }
//
// #endif

void MGEquationsSystem::get_eqs_names(std::vector<string> &FieldsNames){
  
  for(MGEquationsSystem::iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++){
    FieldsNames.push_back(eqn->first);
//     std::cout<<" MGEquationsSystem::get_eqs_names "<< eqn->first <<std::endl;
  }
  return ;
}

// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void MGEquationsSystem::movemesh() {
  std::cout<<"Now we try to move mesh"<<std::endl;
  MGSolBase* mgsoldsx    =  get_eqs("SDSX");
  MGSolBase* mgsoldsy    =  get_eqs("SDSY");
  #if DIMENSION==3
   MGSolBase* mgsoldsz    =  get_eqs("SDSZ");
#endif
//   (mgsolt->x_old[NoLevels-1])->localize(*mgsolt ->x_oold[NoLevels-1]);
  
      const int flag_moving_mesh = (int)(_mgutils._geometry["moving_mesh"]);
          
      const int NoLevels = (int)(_mgutils._geometry["nolevels"]);
      
    const int n_nodes=_mgmesh._NoNodes[NoLevels-1];
    const int n_elem=_mgmesh._NoElements[0][NoLevels-1];
      
    if ( flag_moving_mesh ) {
        /// E) mesh update

//         mgsoldsx->x_old[NoLevels-1]->localize(*x_new);
        const int n_nodes=_mgmesh._NoNodes[NoLevels-1];
        int offsetp=0*n_nodes;
        
                for ( int inode=0; inode<n_nodes; inode++ ) {
                  offsetp=0*n_nodes;
            double disp= (*mgsoldsx->x_old[NoLevels-1])(inode)-(*mgsoldsx->x_oold[NoLevels-1])(inode);
//             _mgmesh._xyz[inode+offsetp] += disp; // cerroni
            _mgmesh._xyz[inode+offsetp]=  _mgmesh._xyzo[inode+offsetp]+(*mgsoldsx->x_old[NoLevels-1])(inode);
            _mgmesh._dxdydz[inode+offsetp] = disp;
            offsetp=1*n_nodes;
            disp= (*mgsoldsy->x_old[NoLevels-1])(inode)-(*mgsoldsy->x_oold[NoLevels-1])(inode);
//             _mgmesh._xyz[inode+offsetp] += disp; cerroni
            _mgmesh._dxdydz[inode+offsetp] = disp;
              _mgmesh._xyz[inode+offsetp]=  _mgmesh._xyzo[inode+offsetp]+(*mgsoldsy->x_old[NoLevels-1])(inode);
            
#if DIMENSION==3
             offsetp=2*n_nodes;
                      disp= (*mgsoldsz->x_old[NoLevels-1])(inode)-(*mgsoldsz->x_oold[NoLevels-1])(inode);
//             _mgmesh._xyz[inode+offsetp] += disp; cerroni
            _mgmesh._dxdydz[inode+offsetp] = disp;
              _mgmesh._xyz[inode+offsetp]=  _mgmesh._xyzo[inode+offsetp]+(*mgsoldsz->x_old[NoLevels-1])(inode);
#endif
        }
        
//         mgsoldsx->x_old[NoLevels-1]->localize(*mgsoldsx ->x_oold[NoLevels-1]);
//         mgsoldsy->x_old[NoLevels-1]->localize(*mgsoldsy ->x_oold[NoLevels-1]);
//         #if DIMENSION==3
//         mgsoldsz->x_old[NoLevels-1]->localize(*mgsoldsz ->x_oold[NoLevels-1]);
// #endif

    }
}

// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
