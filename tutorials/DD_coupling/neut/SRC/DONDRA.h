#ifndef __DONDRA_HXX__
#define __DONDRA_HXX__

#include <map>
#include "Cle2000.hxx"
#include "Solverlib_conf.h"

using namespace boost;
using namespace std;

using namespace ganlib;
// #ifdef HAVE_MED
namespace MEDCoupling {
class MEDLoader;
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}  // namespace MEDCoupling
class InterfaceFunctionDD;
// #endif

class MeshC;
//     struct MeshC;
class MGUtils;
// #include "Communication.hxx"

// ===============================================================
//
// ===============================================================
class DONDRA {
  // ====================================================================================
  // Data interface Donjon neutronic code
 private:
  int _iproc;         ///< proc id
  double _power_;     ///< thermal reactor power in MW
  int _n_fuel_array;  ///< reactor tot # of fuel array
  int _n_array;       ///< reactor tot # of  array
  double _Lref;       ///< reference Length

  ClcmPtr _Fmap;   ///< Fuel map
  ClcmPtr _Matex;  ///< assembly reactor matrix
                   //   ClcmPtr _MacroP;
  ClcmPtr _Cpo;    ///< cpo
                   //   ClcmPtr _CpRefl;   ///< cpo
  ClcmPtr _Track;  ///< Track

#ifdef HAVE_MED
  const MEDCoupling::MEDCouplingUMesh* _med_mesh;        ///< MED-mesh
  std::map<int, InterfaceFunctionDD*> _interfaceFunMap;  ///< MG-Med interface map
#endif

 public:
  LifoPtr ipLifo2;
  Cle2000Ptr PowComponent;
  MeshC* _mesh_c;  ///< Donjon-cartesian mesh

  // ====================================================================================
  // Constructor-Destructor and init functions ------------------------------------------
  /// This constructor  creates a new DONDRA (Donjon) object
  DONDRA();
  // Constructor-Destructor and init functions ------------------------------------------
  /// This constructor  creates a new DONDRA (Donjon) object
  DONDRA(int iproc);
  /// This function initializes the Donjon code
  void initialize(
      double power,  ///< thermal reactor power in MW
      double Lref    ///<  reference Length
  );

  // Set ----------------------------------------------------------------------------
  /// This function sets up the DONDRA object
  void set_up(std::string filename_in  // name init power  (ex IniPowCompo)
  );
#ifdef HAVE_MED
  void set_mesh_med(const MEDCoupling::MEDCouplingUMesh& mg_mesh_med_in  ///< MEDMesh (in)
  ) {
    _med_mesh = &mg_mesh_med_in;
    return;
  }
#endif

  // Returns ----------------------------------------------------------------------------
  /// This function gets the number of fuel assemblies
  int get_n_fuel_array() { return _n_fuel_array; }
  /// This function gets the number of assemblies (fuel + dummy)
  int get_n_array() { return _n_array; }

  // Solve -------------------------------------------------------------
  /// This function executes the DONDRA object
  void solve_onestep();

  // Print --------------------------------------------------------------------------------
  // ======================================================================================
  // ======================================================================================
  /// This function prints in hdf format on linear cartesian grids (HEX27->HEX8)
  /// the following power quantites: n_cmp=1 -> power; n_cmp=2 power+burnup
  void print_power_med2hdf5(
      std::string file_name,                     ///< file name (where to print)
      MEDCoupling::MEDCouplingFieldDouble& phi,  ///< mesh field
      int n_cmp,                                 ///< stampa n_cmp=1->power;n_cmp=2 power+burnup
      std::string name_cmp[]                     ///< dataset name vector  ("power","burnup")
  );
  void print_power_xmf(
      const int t_step, MGUtils& mg_utils, std::string file_name,
      int n_cmp,              ///< stampa n_cmp=1->power;n_cmp=2 power+burnup
      std::string name_cmp[]  ///< dataset name vector  ("power","burnup")
  );
  void print_time_xmf(const int t_in, const int t_end);
  /// 8)
  void print_med(std::string namefile, const MEDCoupling::MEDCouplingFieldDouble& field);

#ifdef HAVE_MED

  // Interface functions -----------------------------------------------------------------
  // ======================================================================================
  // DONDRA:: Interface management Table (see below for details)
  // ======================================================================================
  // 1. void                                        init_interface
  //    (const int interface_name,const int interface_id,const std::string & medfile_name)
  // 2. InterfaceFunctionC * get_interface_fun(int id)
  // 3. MEDCoupling::MEDCouplingFieldDouble *        getValuesPW_elem
  //                                                    (int id,const char *variable_name)
  // 4. MEDCoupling::MEDCouplingFieldDouble *        getValuesPWBUP_elem
  //                                                    (int id,const char *variable_name)
  // 5. void   setBC
  //            (int name,int  n_cmp,const MEDCoupling::MEDCouplingFieldDouble *field)
  // 6. void set_values
  //        (int id_boundary_name,  std::string mgsystem_name, int n_cmp,int first_cmp)

  // ======================================================================================
  /// 1. This function gets the Group mesh from the med file (medfile_name) and sets the id
  /// and the mesh to the interfaces function through the index of the interface-functions
  void init_interface(
      const int interface_name,        ///< interface name            (in)
      const int interface_id,          ///< interface id              (in)
      const std::string& medfile_name  ///< medfile name              (in)
  );
  // ======================================================================================
  /// 2. This function retrun the pointer of the InterfaceFunctionC from the id name
  InterfaceFunctionDD* get_interface_fun(int id  ///< interface identity
  ) {
    return (_interfaceFunMap.find(id) == _interfaceFunMap.end() ? NULL : _interfaceFunMap[id]);
  }
  // ======================================================================================
  /// 3) This function gets the value of power in each element: variable -> "BUND-PW"
  MEDCoupling::MEDCouplingFieldDouble* getValuesPW_elem(
      int id,                    ///< int boundary identity   (in)
      const char* variable_name  ///< system name             (in)
  );
  // ======================================================================================
  /// 4) This function gets the value of power and burnup in each element:
  /// variable -> "BUND-PW"   variable -> "BUND-BUP"
  MEDCoupling::MEDCouplingFieldDouble* getValuesPWBUP_elem(
      int id,                    ///< int boundary identity   (in)
      const char* variable_name  ///< system name             (in)
  );

  // ======================================================================================
  /// 5) This function sets the field "field" to the interface with name "name"
  void setBC(
      int name,  //< interface name        (in)
                 //   int  n_cmp,                                     //< number of componenets (in)
      const MEDCoupling::MEDCouplingFieldDouble* field  //< Med-field             (in)
  );
  // ====================================================================================
  /// 6)
  void set_tc_values(
      const MEDCoupling::MEDCouplingFieldDouble& bdy,
      int id_boundary_name,       ///< identity interface name (in)
      std::string mgsystem_name,  ///< system name             (in)
      int itime

  );
  // ======================================================================================
  /// 7) This function sets the interface field with n_comp componenets starting from
  /// the first_cmp to the interface with name "id_boundary_name"
  ///
  void set_values(
      int id_boundary_name,       ///< identity interface name           (in)
      std::string mgsystem_name,  ///< system name                       (in)
      int n_cmp,                  ///<  variable system                  (in)
      int first_cmp               ///< from variable system              (in)

  );
  // ================================================================================

#endif

  // ---------------------------------------------------------------------------------------
};  // end class DONDRA (interface donjon-dragon)

#endif
