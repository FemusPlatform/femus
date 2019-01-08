#ifndef __mg_system2__
#define __mg_system2__

// std libraries ---------------
#include <map>
#include <string>
#include <vector>

// configure files ----------------
// #include "Equations_tab.h"
// #include "MGSclass_conf.h"
// #include "MGMesh.h"
// Forward class -------------
class MGUtils;
class MGMesh;

// ========================================
//            MGSystem class
// ========================================
/// This class contains information about the physical properties and system 
class MGSystem {
public: 
// protected:
  // data ----------------------------------
  // mesh
  MGUtils&   _mgutils;  ///< MGUtils pointer
  MGMesh&    _mgmesh;   ///< MGMesh  pointer
  // parameters
  std::map<std::string,double> parameters;   ///< parameter map
public:  
  // cell or pt data field
  int       _n_data[2];       ///< n cell and point data  (0=point 1=cell)
  double*   _sys_data[6];   ///< pt or cell data field

// // Only for Femlcore app
//   int _flag_probe;
//   int _node_glob_probe[1000];
// only for Femlcore app


  // Nondimensional groups and ref values -----------

// #ifdef T_EQUATIONS
//   double _Pr;    ///< Prandl Number
// #endif
//  #ifdef MHD_EQUATIONS
//   double _S;     ///< Coupling Magnetic coeff
//   double _Rem;   ///< Magnetic Reynolds number
//   double _Hm;    ///< Hartmann Number
//  #endif
// #ifdef TWO_PHASE
//   double _We;    ///< Weber Number
// #endif
//-------------------------------------------------------------------------
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  MGSystem(MGUtils& mgutils_in,MGMesh& mgmesh_in,int n_pt_in=0,int n_cell_in=0); 
  ~MGSystem(){clean();}   
  void clean(){for(int i=0;i<_n_data[0]+_n_data[1];i++) delete []_sys_data[i];} 
  ///@}
  // Temperature dependence --------------
#ifndef PARAM_CONST
  ///@{ \name NONDIMENSIONAL PARAMETERS
  inline double adensity(double /*Temp*/) const {return 1.;}   ///< Density (rho/rhoref)
  inline double aviscosity(double /*Temp*/) const {return 1.;} ///< Viscosity (mu/muref)
  inline double adensityT(double /*Temp*/) const {return 1.;}   ///< Density
  inline double akappaT(double /*Temp*/) const {return 1.;}     ///< Conductivity
  
  /// Specific heat
  inline double acipiT(double /*Temp*/) const {return 1.;}      
  ///@}
#endif
  // -------------------------------------------------------------
  // PARAMETERS  
  // --------------------------------------------------------------
  // --------------------------------------------------------------
  ///@{ \name RETURN FUNCTIONS 
  inline double get_par(const std::string & name) const {
    return parameters.find(name)->second;  } ///< Return map value
  ///@}
  // -------------------------------------------------------------
  ///@{ \name SET FUNCTIONS
   
 		    
  
  /// Set a value in the map of parameters
  inline void   set_par(const std::string & name, double value)  {
    parameters.insert(make_pair(name,value));}  
  ///@}
  //-------------------------------------------------------------
  ///@{ \name READ-PRINT    
  void read_par();        ///< read parameter file
  
  /// print parameter map
  void print_par() const;   
  ///@}
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  ///@{ \name SYSTEM DATA  
  void read_data(const std::string & file_name);  ///<  Read data
  void clear_data();                              ///<  Clear data
  
  virtual void  print_data_view(
  const std::string & filename, // filename
  int i_choice,
  std::string dir_name
) const;                                         ///<
   void print_data(const std::string & file_name) const;
  void init_data(const int
 mode = 0);
  void print_xml_attrib(std::ofstream &out,int nelements,int nodes,std::string file_name) const; 
  void print_xml_mat(std::ofstream &out,int nelements,int nodes,std::string file_name) const; 
  /// Point data field
  double F_ext(const double xx[],const int
 ivarq) const;
  /// Return point data
    inline double get_sys_data(const int inode,const int iflag)const {
    return _sys_data[iflag][inode];
  }
  ///@}

  
private:
   // set nondimensional groups and ref values
//   void set_nondimgroups();

};

#endif
