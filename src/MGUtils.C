// Class includes ---------------
#include "MGUtils.h"
#include "Solverlib_conf.h"  // petsc conf
#include "Printinfo_conf.h"  // petsc conf

// std libraries ----------------
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

// ===============================
//  MGUtils Class functions
// ===============================
// ===================================================
/// Constructor
MGUtils::MGUtils() :_ProbID(-1){
  read();
  StandardBuild();
}
MGUtils::MGUtils(int a) :_ProbID(a){

  std::ostringstream loc_par_file;
  loc_par_file << "/DATA/param_files_msh"<< a <<".in";
  const std::string &name=loc_par_file.str();
  read(name);
  StandardBuild();
  
  std::string mesh_nameP = get_file("F_MESH_READ"); // name mesh
  std::ostringstream filenameP;
  int posP = mesh_nameP.find(".");  // position of "live" in str
  filenameP <<   mesh_nameP.substr(0,posP)  << "_MedToMg.med" ;
  
  _interface_mesh = filenameP.str();
}

void MGUtils::StandardBuild(){
  _user_dir = _femus_dir +  "/USER_APPL/";
  _app_dir= _user_dir + _myapp_name;
  _inout_dir = _app_dir + "/RESU/";
  _mesh_dir = _user_dir + _mgfiles.find("MESH_DIR")->second ;
//     _mesh_dir = _user_dir + "/MESH/";
  _data_dir = _app_dir + "/DATA/";
  _fem_dir = _femus_dir + "/fem/";
  _contrib_dir=_femus_dir + "/contrib/";
  read_par();       // read parameters
#ifdef PRINT_INFO
  print();
  print_par();      // print in console parameters
#endif 
}

// ==============================================
// // //  OLD READ_PAR FUNCTION
/// This function reads all the  parameters
// void MGUtils::read_par() {
//
//     // file names ------------------
//
//     std::string    baseparutils = get_file("BASEPARUTILS");
//
//     std::ostringstream filename;
//     filename << _data_dir << baseparutils;
//  std::cout << "Init Reading" << filename.str() <<  std::endl;
//     std::ifstream fin;
//     fin.open(filename.str().c_str());
//     std::string buf="";
//     double value;
// #ifdef PRINT_INFO
//     std::cout << "Init Reading" << filename.str() <<  std::endl;
// #endif
//     if (fin.is_open()) {
//         while (buf != "/" )  fin >> buf;
//         fin >> buf;
//         while (buf != "/" ) {
//
//             if (buf == "#") getline(fin, buf); // comment line
//             else {
//                 fin >> value;    // set new parameter
//                 set_par(buf,value);
//             }
//             fin >> buf;// std::cerr <<buf.c_str() << "\n ";
//         }
//     }
//     else {
//         std::cerr << "MGUtils::read_par: no parameter file found" << std::endl;
//         abort();
//     }
// #ifdef PRINT_INFO
//     std::cout << "End Reading file " <<  filename.str() << std::endl;
// #endif
//     return;
// }

void MGUtils::read_par() { // READ PARAMETER FILES AND FILL RELATIVE MAPS ===========================

// file names ------------------
  std::vector<std::string> maps; maps.resize(4); // PARAMETER files
  maps[0] = get_file("BASEPARUTILS");   maps[1] = get_file("GEOM_PAR"); maps[2] = get_file("SIM_CONFIG"); maps[3] = get_file("MAT_PROP");

  double double_value;  std::string string_value; std::string buf=""; // read double, string, dummay

  for(int i=1; i<4; i++) { // CYCLE ON FILES TO READ -------------------------------------------
    std::ostringstream filename;   filename << _data_dir << maps[i]; // file name
    std::ifstream fin;  fin.open(filename.str().c_str()); // stream file
    buf="";
#ifdef PRINT_INFO
    if(fin.is_open()) {  std::cout << "Init Reading = " << filename.str() <<  std::endl; }
#endif
    if(fin.is_open()) {  // -------------------------------------------------------------
      while(buf != "/") { fin >> buf; }  // find "/" file start
      fin >> buf;
      while(buf != "/") {
        if(buf == "#") { getline(fin, buf); }  // comment line
        else {
          if(i==0) {fin >> double_value; set_par(buf,double_value); }
          else if(i==1) {  fin >> double_value;   set_geom_par(buf,double_value); }
          else if(i==2) {  fin >> string_value;   set_sim_par(buf,string_value); }
          else {           fin >> double_value;   set_mat_par(buf,double_value); }
        }
        fin >> buf;// std::cerr <<buf.c_str() << "\n ";
      }
    } // --------------------------------------------------------------------------------------
    else {
      std::cerr << "MGUtils::read_par: no parameter file found" << std::endl;   abort();
    }
#ifdef PRINT_INFO
    std::cout << "End Reading file " <<  filename.str() << std::endl;
#endif   
    fin.close();
  }// END CYCLE ON FILES TO READ ----------------------------------------------------------
  maps.clear();
  

// _SolverTypeMap["CGNM"]             =  CGNM              ;
// _SolverTypeMap["CGSM"]             =  CGSM             ;
// _SolverTypeMap["CRM"]              =  CRM               ;
// _SolverTypeMap["QMRM"]             =  QMRM              ;
// _SolverTypeMap["TCQMRM"]           =  TCQMRM            ;
// _SolverTypeMap["TFQMRM"]           =  TFQMRM            ;
// _SolverTypeMap["BICGM"]            =  BICGM             ;
// _SolverTypeMap["BICGSTABM"]        =  BICGSTABM         ;
// _SolverTypeMap["MINRESM"]          =  MINRESM           ;
// _SolverTypeMap["GMRESM"]           =  GMRESM            ;
// _SolverTypeMap["VANKATM"]          =  VANKATM           ;
// _SolverTypeMap["VANKANSM"]         =  VANKANSM          ;
// _SolverTypeMap["LSQRM"]            =  LSQRM             ;
// _SolverTypeMap["JACOBIM"]          =  JACOBIM           ;
// _SolverTypeMap["SOR_FORWARDM"]     =  SOR_FORWARDM      ;
// _SolverTypeMap["SOR_BACKWARDM"]    =  SOR_BACKWARDM     ;
// _SolverTypeMap["SSORM"]            =  SSORM         ;
// _SolverTypeMap["RICHARDSONM"]      =  RICHARDSONM       ;
// _SolverTypeMap["CHEBYSHEVM"]       =  CHEBYSHEVM;   
// _SolverTypeMap["LUMPM"]            =  LUMPM        ;
// _SolverTypeMap["INVALID_SOLVERM"]  =  INVALID_SOLVERM    ;
  
  return;
}// END READ_PAR FUNCTION ================================================================

// =================================================
/// This function reads the file names
void MGUtils::read(const std::string &name_file_in) {

  // read femus dir from shell -----------------------
  _femus_dir = getenv("FEMUS_DIR");
  if(_femus_dir == "") {
    std::cout << " Set FEMUS_DIR in your shell environment" << std::endl;
    abort();
  }

  _myapp_name = getenv("FM_MYAPP");
  if(_myapp_name == "") {
    std::cout << " Set MYAPP in your shell environment" << std::endl;
    abort();
  }

#ifdef PRINT_INFO
  std::cout << " femus_dir is " << _femus_dir << std::endl;
  std::cout << " myapp_name is " << _myapp_name << std::endl;
#endif

  // read param_files -----------------------------
  std::ostringstream filename;
//   filename << _femus_dir << "/" << "config/param_files.in";

  filename << _femus_dir <<  "/USER_APPL/"  <<  _myapp_name  <<"/" <<
           name_file_in;

//   "/config/param_files.in";
#ifdef PRINT_INFO
  std::cout << filename.str().c_str() <<" ";
#endif

  std::ifstream fin(filename.str().c_str());
  std::string buf="";
  std::string value;

  if(fin != NULL) {
    while(!fin.eof()) {
      fin >> buf;
      if(buf == "#") { fin.ignore(200,'\n'); }
      else { fin >> value;  set(buf,value);} // set new parameter
    }
  } else {
    std::cerr << " MGFiles::read: no parameter file found" << std::endl;
    abort();
  }
  // cleaning and check
  fin.close();
//   check_dirs();

  return;
}

// ============================================================
/// This file prints the file name map
void MGUtils::print() {
  std::cout << "\n ================================================ ";
  std::cout << "\n Class MGFiles: "<< (int)_mgfiles.size()
            << " file names: " << std::endl;
  std::map<std::string,std::string>::const_iterator pos=_mgfiles.begin();
  std::map<std::string,std::string>::const_iterator pos_e=_mgfiles.end();
  for(; pos!=pos_e; pos++) {
    std::string  name=pos->first; // get the string
    std::string value=pos->second;// get the name
    std::cout  <<" "<< std::left << std::setw(15)
               << name << " = " << value << std::endl;
  }
  return;
}

// ============================================================
/// This function prints all the  parameters to stdout
void MGUtils::print_par() const {
  std::cout << "\n ============================= ";
  std::cout << "\n   MGUtils: " << (int)_param_utils.size() << " parameters: \n";

  std::map<std::string,double>::const_iterator   pos=_param_utils.begin();
  std::map<std::string,double>::const_iterator pos_e=_param_utils.end();
  for(; pos!=pos_e; pos++) {
    std::string name=pos->first; // get name
    double value=pos->second;    // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = "
              << std::setprecision(12) << value << std::endl;
  }

  std::cout << " ---------------- GEOMETRY PARAMETERS ---------------\n  " ;

  std::map<std::string,double>::const_iterator   pos1=_geometry.begin();
  std::map<std::string,double>::const_iterator pos_e1=_geometry.end();
  for(; pos1!=pos_e1; pos1++) {
    std::string name=pos1->first; // get name
    double value=pos1->second;    // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = "
              << std::setprecision(12) << value << std::endl;
  }

  std::cout << " ---------------- SIMULATION PARAMETERS ---------------\n  " ;

  std::map<std::string,std::string>::const_iterator   pos2=_sim_config.begin();
  std::map<std::string,std::string>::const_iterator pos_e2=_sim_config.end();
  for(; pos2!=pos_e2; pos2++) {
    std::string name=pos2->first; // get name
    std::string value=pos2->second;    // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = "
              << std::setprecision(12) << value << std::endl;
  }

  return;
}

void MGUtils::FillFieldsVector(std::vector<FIELDS> &myproblemP) {
  if(stoi(_sim_config["DA"])!=0) { myproblemP.push_back(DA_F) ; }
  
  if(stoi(_sim_config["NavierStokes"])!=0) { myproblemP.push_back(NS_F) ; }
  if(stoi(_sim_config["FluidStructure"])!=0) { myproblemP.push_back(NS_F) ; }
  if(stoi(_sim_config["StructuralMechanics"])!=0) { myproblemP.push_back(NS_F) ; }
  
  if(stoi(_sim_config["AdjointNavierStokes"])!=0) { myproblemP.push_back(NSA_F) ; }
  if(stoi(_sim_config["AdjointFluidStructure"])!=0) { myproblemP.push_back(NSA_F) ; }
  
  if(stoi(_sim_config["Temperature"])!=0) { myproblemP.push_back(T_F) ; }
  if(stoi(_sim_config["DynamicalTurbulence"])!=0) { myproblemP.push_back(K_F) ; }
  if(stoi(_sim_config["ThermalTurbulence"])!=0) { myproblemP.push_back(KTT_F) ; }
  
  if(stoi(_sim_config["Displacement"])!=0) { myproblemP.push_back(SDS_F) ; }
  
  if(stoi(_sim_config["Laplacian"])!=0) { myproblemP.push_back(TA_F) ; }
  if(stoi(_sim_config["AdjointTemperature"])!=0) { myproblemP.push_back(TA_F) ; }
  
  if(stoi(_sim_config["AdjointTurbulence"])!=0) { myproblemP.push_back(KA_F) ; }
  
  if(stoi(_sim_config["ColorFunction"])!=0) { myproblemP.push_back(CO_F) ; }
  if(stoi(_sim_config["Control"])!=0) { myproblemP.push_back(CTRL_F) ; } 
  
  return;
}


#if HDF5_VERSIONM == 188
// =============================================================
hid_t MGUtils::read_Dhdf5(hid_t file,const std::string &name,double data[]) {

  hid_t  dataset = H5Dopen(file,name.c_str());
  hid_t status=H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

hid_t MGUtils::print_Dhdf5(hid_t file,const std::string &name, hsize_t dimsf[],double data[]) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

/// Print int data into dhdf5 file

hid_t MGUtils::print_Ihdf5(hid_t file,const std::string &name, hsize_t dimsf[],int data[]) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t MGUtils::read_Ihdf5(hid_t file,const std::string &name,int data[]) {

  hid_t  dataset = H5Dopen(file,name.c_str());
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

#else
// =============================================================
hid_t MGUtils::read_Dhdf5(hid_t file,const std::string &name,double data[]) {

  hid_t  dataset = H5Dopen(file,name.c_str()
#if HDF5_VERSIONM != 1808
                           , H5P_DEFAULT
#endif
                          );
  hid_t status=H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

hid_t MGUtils::print_Dhdf5(hid_t file,const std::string &name, hsize_t dimsf[],double data[]) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
                            dataspace,
#if HDF5_VERSIONM != 1808
                            H5P_DEFAULT, H5P_DEFAULT,
#endif
                            H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

/// Print int data into dhdf5 file

hid_t MGUtils::print_Ihdf5(hid_t file,const std::string &name, hsize_t dimsf[],int data[]) {

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace,
#if HDF5_VERSIONM != 1808
                            H5P_DEFAULT, H5P_DEFAULT,
#endif
                            H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t MGUtils::read_Ihdf5(hid_t file,const std::string &name,int data[]) {

  hid_t  dataset = H5Dopen(file,name.c_str()
#if HDF5_VERSIONM != 1808
                           , H5P_DEFAULT
#endif
                          );
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}
#endif
