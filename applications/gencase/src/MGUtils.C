// Class includes ---------------
#include "MGUtils.h"
#include "Printinfo_conf.h"  // petsc conf
#include "Solverlib_conf.h"  // petsc conf
// std libraries ----------------
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// ===============================
//  MGUtils Class functions
// ===============================
// ===================================================
/// Constructor
MGUtils::MGUtils() {
  _app_dir = getenv("APP_PATH");
  _inout_dir = _app_dir + "/RESU/";
  _data_dir = _app_dir + "/DATA/";
  _contrib_dir = _femus_dir + "/contrib/";
  read();
  _mesh_dir = _app_dir + "/" + _mgfiles.find("MESH_DIR")->second;
  read_par();  // read parameters
#ifdef PRINT_INFO
  print();
  print_par();  // print in console parameters
#endif
}
MGUtils::MGUtils(int a) {
  std::ostringstream loc_par_file;
  loc_par_file << "/DATA/param_files_msh" << a << ".in";
  const std::string& name = loc_par_file.str();

  _app_dir = getenv("APP_PATH");
  _inout_dir = _app_dir + "/RESU/";
  _data_dir = _app_dir + "/DATA/";
  if (NUM_MESH > 1) { _data_dir = _app_dir + "/DATA/DATA" + std::to_string(a) + "/"; }
  _contrib_dir = _femus_dir + "/contrib/";

  read(name);
  _mesh_dir = _app_dir + "/" + _mgfiles.find("MESH_DIR")->second;
  read_par();  // read parameters
#ifdef PRINT_INFO
  print();
  print_par();  // print in console parameters
#endif
}

void MGUtils::read_par() {  // READ PARAMETER FILES AND FILL RELATIVE MAPS ===========================

  // file names ------------------
  std::vector<std::string> maps;
  maps.resize(4);  // PARAMETER files
  maps[0] = get_file("BASEPARUTILS");
  maps[1] = get_file("GEOM_PAR");
  maps[2] = get_file("SIM_CONFIG");
  maps[3] = get_file("MAT_PROP");

  double double_value;
  std::string string_value;
  std::string buf = "";  // read double, string, dummay

  for (int i = 1; i < 4; i++) {  // CYCLE ON FILES TO READ -------------------------------------------
    std::ostringstream filename;
    filename << _data_dir << maps[i];  // file name
    std::ifstream fin;
    fin.open(filename.str().c_str());  // stream file
    buf = "";
#ifdef PRINT_INFO
    if (fin.is_open()) { std::cout << "Init Reading = " << filename.str() << std::endl; }
#endif
    if (fin.is_open()) {                  // -------------------------------------------------------------
      while (buf != "/") { fin >> buf; }  // find "/" file start
      fin >> buf;
      while (buf != "/") {
        if (buf == "#") {
          getline(fin, buf);
        }  // comment line
        else {
          if (i == 0) {
            fin >> double_value;
            set_par(buf, double_value);
          } else if (i == 1) {
            fin >> double_value;
            set_geom_par(buf, double_value);
          } else if (i == 2) {
            fin >> string_value;
            set_sim_par(buf, string_value);
          } else {
            fin >> double_value;
            set_mat_par(buf, double_value);
          }
        }
        fin >> buf;  // std::cerr <<buf.c_str() << "\n ";
      }
    }  // --------------------------------------------------------------------------------------
    else {
      std::cerr << "MGUtils::read_par: no parameter file found" << std::endl;
      abort();
    }
#ifdef PRINT_INFO
    std::cout << "End Reading file " << filename.str() << std::endl;
#endif

    fin.close();
  }  // END CYCLE ON FILES TO READ ----------------------------------------------------------
  maps.clear();
  return;
}  // END READ_PAR FUNCTION ================================================================

// =================================================
/// This function reads the file names
void MGUtils::read(const std::string& name_file_in) {
  // read femus dir from shell -----------------------
  _femus_dir = getenv("FEMUS_DIR");
  if (_femus_dir == "") {
    std::cout << " Set FEMUS_DIR in your shell environment" << std::endl;
    abort();
  }

  _myapp_name = getenv("FM_MYAPP");
  if (_myapp_name == "") {
    std::cout << " Set MYAPP in your shell environment" << std::endl;
    abort();
  }

#ifdef PRINT_INFO
  std::cout << " femus_dir is " << _femus_dir << std::endl;
  std::cout << " myapp_name is " << _myapp_name << std::endl;
#endif

  // read param_files -----------------------------
  std::ostringstream filename;

  filename << _app_dir << "/" << name_file_in;

#ifdef PRINT_INFO
  std::cout << filename.str().c_str() << " ";
#endif

  std::ifstream fin(filename.str().c_str());
  std::string buf = "";
  std::string value;

  while (!fin.eof()) {
    fin >> buf;
    if (buf == "#") {
      fin.ignore(200, '\n');
    } else {
      fin >> value;
      set(buf, value);
    }  // set new parameter
  }
  fin.close();

  return;
}

// ============================================================
/// This file prints the file name map
void MGUtils::print() {
  std::cout << "\n ================================================ ";
  std::cout << "\n Class MGFiles: " << (int)_mgfiles.size() << " file names: " << std::endl;
  std::map<std::string, std::string>::const_iterator pos = _mgfiles.begin();
  std::map<std::string, std::string>::const_iterator pos_e = _mgfiles.end();
  for (; pos != pos_e; pos++) {
    std::string name = pos->first;    // get the string
    std::string value = pos->second;  // get the name
    std::cout << " " << std::left << std::setw(15) << name << " = " << value << std::endl;
  }
  return;
}

// ============================================================
/// This function prints all the  parameters to stdout
void MGUtils::print_par() const {
  std::cout << "\n ============================= ";
  std::cout << "\n   MGUtils: " << (int)_param_utils.size() << " parameters: \n";

  std::map<std::string, double>::const_iterator pos = _param_utils.begin();
  std::map<std::string, double>::const_iterator pos_e = _param_utils.end();
  for (; pos != pos_e; pos++) {
    std::string name = pos->first;  // get name
    double value = pos->second;     // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }

  std::cout << " ---------------- GEOMETRY PARAMETERS ---------------\n  ";

  std::map<std::string, double>::const_iterator pos1 = _geometry.begin();
  std::map<std::string, double>::const_iterator pos_e1 = _geometry.end();
  for (; pos1 != pos_e1; pos1++) {
    std::string name = pos1->first;  // get name
    double value = pos1->second;     // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }

  std::cout << " ---------------- SIMULATION PARAMETERS ---------------\n  ";

  std::map<std::string, std::string>::const_iterator pos2 = _sim_config.begin();
  std::map<std::string, std::string>::const_iterator pos_e2 = _sim_config.end();
  for (; pos2 != pos_e2; pos2++) {
    std::string name = pos2->first;    // get name
    std::string value = pos2->second;  // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " << std::setprecision(12) << value
              << std::endl;
  }

  return;
}

#if HDF5_VERSIONM == 188
// =============================================================
hid_t MGUtils::read_Dhdf5(hid_t file, const std::string& name, double data[]) {
  hid_t dataset = H5Dopen(file, name.c_str());
  hid_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}

hid_t MGUtils::print_Dhdf5(hid_t file, const std::string& name, hsize_t dimsf[], double data[]) {
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

/// Print int data into dhdf5 file

hid_t MGUtils::print_Ihdf5(hid_t file, const std::string& name, hsize_t dimsf[], int data[]) {
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(file, name.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t MGUtils::read_Ihdf5(hid_t file, const std::string& name, int data[]) {
  hid_t dataset = H5Dopen(file, name.c_str());
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}

#else
// =============================================================
hid_t MGUtils::read_Dhdf5(hid_t file, const std::string& name, double data[]) {
  hid_t dataset = H5Dopen(
      file, name.c_str()
#if HDF5_VERSIONM != 1808
                ,
      H5P_DEFAULT
#endif
  );
  hid_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}

hid_t MGUtils::print_Dhdf5(hid_t file, const std::string& name, hsize_t dimsf[], double data[]) {
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
      file, name.c_str(), H5T_NATIVE_DOUBLE, dataspace,
#if HDF5_VERSIONM != 1808
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

/// Print int data into dhdf5 file

hid_t MGUtils::print_Ihdf5(hid_t file, const std::string& name, hsize_t dimsf[], int data[]) {
  hid_t dataspace = H5Screate_simple(2, dimsf, NULL);
  hid_t dataset = H5Dcreate(
      file, name.c_str(), H5T_NATIVE_INT, dataspace,
#if HDF5_VERSIONM != 1808
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ===========================================================================
hid_t MGUtils::read_Ihdf5(hid_t file, const std::string& name, int data[]) {
  hid_t dataset = H5Dopen(
      file, name.c_str()
#if HDF5_VERSIONM != 1808
                ,
      H5P_DEFAULT
#endif
  );
  hid_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  return status;
}
#endif
