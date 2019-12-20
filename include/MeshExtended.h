#ifndef __MeshExtendedM__
#define __MeshExtendedM__

#include "MGMesh.h"
#include <vector>
#include "Solverlib_conf.h" 

// #ifdef HAVE_MED
namespace MEDCoupling {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
class DataArrayInt;
class DataArrayDouble;
}
// #endif
class MeshExtended : public MGMesh {
public:
   std::vector<int>  _bc_id;
   std::vector<int> _mat_id;
  // Constructor-Destructor
  MeshExtended(        const ParallelObjectM &comm,
                       MGUtils& mgutils,
                       MGGeomEl& geomel
                      );
  
  virtual ~MeshExtended();

  void read_bc_id(int Level);
  void read_mat_id(int Level);
  
  void  B_dist(const  int Level); ///< Calculate distance from a wall
  void print_mat_hf5(std::string filename,std::string dir_name) const;
  void print_bc_hf5(std::string filename,std::string dir_name) const;
#ifdef HAVE_MED
  void print_med(int Level,std::string filename);
#endif
  void restart_cell_array(double array[], std::string FieldName);
};

#endif
