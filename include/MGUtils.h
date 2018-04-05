#ifndef __mgutils_h__
#define __mgutils_h__

// configure files --------------------
#include "Domain_conf.h"
#include "Equations_tab.h"
// std libraries ----------------------
#include <iostream>
#include <cstdlib>
#include <map>
#include <string>
#include "hdf5.h"
// #ifdef TBK_EQUATIONS
#include "TurbUtils.h"
// #include "IbUtils.h"
// #endif
// Forwarding class -------------------
// class MGFiles;

class IbUtils;
// class TurbUtils;
// =========================================
//              MGUtils
// =========================================
/// This class contains information about parameters and directory path
class MGUtils  {

protected:
  // data ---------------------------------------
  std::map<std::string, double>       _param_utils; ///< Parameters map  
  int                                 _name;        ///< Problem name
  int                           _ProbID;
public:
  std::map<std::string, double>       _geometry;    ///< Geometry properties: bounding box for libmesh mesh generation, axisym, wall_dist
  std::map<std::string, double>       _mat_prop;    ///< Material properties: density, viscosity and other physical prop.
  std::map<std::string, std::string>  _mgfiles;     ///< Files map
  std::map<std::string, std::string>  _sim_config;  ///< Simulation properties: dt, print step, n. of steps, etc.  
    
  double _factor = 1.;  
  double _press  = 0.;
  ///@{ \name USEFUL DIRECTORY NAMES
  std::string                       _user_dir;      ///< User directory
  std::string                       _data_dir;      ///< Data directory
  std::string                       _app_dir;       ///< Application directory
  std::string                       _femus_dir;     ///< Femus directory
  std::string                       _myapp_name;    ///< Application name
  std::string                       _inout_dir;     ///< RESU directory
  std::string                       _mesh_dir;      ///< MESH directory
  std::string                       _fem_dir;       ///< FEM directory
  std::string                       _contrib_dir;   ///< FEM directory
  std::string                       _interface_mesh;
//------------------------------------------------------------------------
// #ifdef TBK_EQUATIONS
 TurbUtils                        *_TurbParameters = NULL;
 IbUtils                          *_IBParameter = NULL;
// #endif
  ///@}
  
//   MGFiles&           _files;     ///< MGFiles class pointer

  // Constructor-Destructor --------------------
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGUtils(/*MGFiles& mgfiles_in*/);  ///< Constructor -> read and print parameters 
  MGUtils(int a);  ///< Constructor -> read and print parameters 
  MGUtils(int a, TurbUtils TurbParameters);
  
  void StandardBuild();
  void StandardBuild(int a);
  /// Destructor
  ~MGUtils() {clean();}
  void clean(){
    _param_utils.clear();
    _mgfiles.clear();
    if(_TurbParameters != NULL) {
      _TurbParameters->~TurbUtils();
      delete _TurbParameters;
    }
//     if(_IBParameter != NULL) {
//       _IBParameter->~IbUtils();
//       delete _IBParameter;
//     }
  }
   ///@}
  //--------------------------------------------------------------------------------------
  ///@{ \name GET/SET
  inline double get_par(const std::string & name) const {
    if(_param_utils.end() == _param_utils.find(name))  {std::cout << "MGUtils: Get par failed ";abort();}
    return _param_utils.find(name)->second;
  }                               ///< Return a parameter from the corresponding map
  inline std::string get_file(const std::string & name) const {
    return _mgfiles.find(name)->second; 
  }                               ///< Return file name from the corresponding map
  
  inline void   set(const std::string & name, std::string & value)  {
    _mgfiles.insert(make_pair(name,value));
  }                               ///< Set file name into the corresponding map
    inline void   set_name(const int value)  {
    _name=value;
  }
    inline int  get_name()  {
    return _name;
  }
  
  /// Set a parameter in the corresponding map
  inline void   set_par(const std::string & name, double value)  {
    _param_utils[name] = value;
  }
  inline void   set_geom_par(const std::string & name, double value)  {
    _geometry[name] = value;
  }
  inline void   set_mat_par(const std::string & name, double value)  {
    _mat_prop[name] = value;
  }
  inline void   set_sim_par(const std::string & name, std::string value)  {
    _sim_config[name] = value;
  }
  
// #ifdef TBK_EQUATIONS
  inline void set_turbulence_info(TurbUtils *Parameters){
    _TurbParameters =  Parameters;
  }
  inline void set_IB_info(IbUtils *Parameters){
    _IBParameter =  Parameters;
  }
// #endif
  ///@}
  //-----------------------------------------------------------------------------------------
  ///@{ \name READ-PRINT (IN CONSTRUCTOR)
  void read_par();                                            ///< Read parameters from file
  void read(const std::string & name="/DATA/param_files.in"); ///< Read file names from file
  void print();                                               ///< Print in console the file names
  
  /// Print in console the parameters read
  void print_par() const;   
  ///@}
  
  //hdf5 ------------------------------------
  ///@{ \name READ-PRINT HDF5
  ///< \param  <file>   HDF5 file name (hid_t format)    
  ///< \param  <name>   Directory name inside the HDF5 file 
  ///< \param  <dimsf>  Dimension of the vector in the directory <name>
  ///< \param  <data>   Vector where put or get data
  
  //-----------------------------------------------------------------------------------------
  hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double data[]); ///< Print vector of double
  hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int data[]);    ///< Print vector of integer
  hid_t read_Dhdf5(hid_t file,const std::string & name,double data[]);                   ///< Read vector of double
  
  /// Read vector of integer
  hid_t read_Ihdf5(hid_t file,const std::string & name,int data[]);                   
 ///@}
 void cross(const double* a,const double* b, double* res) const;
 
 void FillFieldsVector(std::vector<FIELDS> &myproblemP);
};

inline void MGUtils::cross(const double* a,const double* b, double* res) const {
//a,b,res are 3D vectors
//clean then fill
//   for (uint i=0; i<3; i++) res[i]=0.;
  for (int i=0; i<3; i++) res[i] = (a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3]);
  return;
}
#endif
