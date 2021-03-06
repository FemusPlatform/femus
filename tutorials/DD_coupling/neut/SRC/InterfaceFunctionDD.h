#ifndef __INTERFACE_FUNCTION__
#define __INTERFACE_FUNCTION__

#include "Solverlib_conf.h"

#ifdef HAVE_MED
#include <iostream>
#include <map>
#include <vector>

namespace MEDCoupling {
class MEDCouplingFieldDouble;
class MEDCouplingUMesh;
}  // namespace MEDCoupling

class MeshC;
//

// class MeshExtended;

class InterfaceFunctionDD {
 protected:
  // data interface function ==================================================
  //   const MeshExtended                 * _mesh_mg;     ///< MGMesh
  const MEDCoupling::MEDCouplingUMesh* _support_med;  ///< MEDCoupling mesh
  MEDCoupling::MEDCouplingFieldDouble* _field;        ///< field

  int _n;         ///< number of nodes/elements in the interface
  int* _map_med;  ///< node/element map index->MGMesh  (size: _n)
  int* _map_mg;   ///< node/element map index->MEDCoupling (size: _n)

 public:
  // Constructors- Destructors --------------------------------------
  InterfaceFunctionDD()
      :                      ///< Constructor ============================
                             //      _mesh_mg(NULL),       // Femus-mesh
        _support_med(NULL),  // med-sub/mesh
        _field(NULL),        // function
        _map_med(NULL),
        _map_mg(NULL) {}
  ~InterfaceFunctionDD();  ///< Destructor =========================

  void eval(         ///< Evaluation function ====================
      int node_med,  ///< index node med <-             (in)
      int n_cmp,     ///< index componenent         (in)
      double val[]   ///< node field vector         (out)
  );

  // set functions ------------------------------------------------------------
  // ============================================================================
  void set_mesh_interface(    ///< Setting the interface-meshes
      const bool nodes_only,  ///< nodes only flag
                              //     const MeshExtended * mesh,                          ///< MGMesh
      const MEDCoupling::MEDCouplingUMesh* support  ///< MedCouplng <-
  ) {}

  void set_analytic_field(  ///< Setting the field (analytic)=
      const char* code,     // symbolic expression  <-
      int nComp);           // n of components      <-
  void set_analytic_field_elem(
      const char* symbolic_eq,  // symbolic function
      int nComp                 // number of componenents
  );

  void set_field(                                     ///< Setting the field (MEDfield)=================
      const MEDCoupling::MEDCouplingFieldDouble* f);  // field  <-
  void set_field_source(                              ///< Setting the field (MEDfield)=================
      const MEDCoupling::MEDCouplingFieldDouble* f);  // field  <-

  int get_n() { return _n; }  // ===============
  int* get_map_mg() {
    if (_map_mg == NULL) return NULL;
    return _map_mg;
  }  // ===============
  int* get_map_med() {
    if (_map_med == NULL) return NULL;
    return _map_med;
  }  // ===============
  MEDCoupling::MEDCouplingFieldDouble* getField(const char* vName);
  const MEDCoupling::MEDCouplingUMesh* getSupport() { return _support_med; }

  // ===========================================================================
  // This function prints the interface function
  void printOn(
      std::ostream& out,  ///< ostream file
      int id              ///< interface name
      ) const;
  // ============================================================================
  /// This function computes the mesh interface
  /// to set in the interface-function. The field in the
  /// interface-function is set by set_mesh_femus_interface
  void set_mesh_BFinterface_nodeID(
      //   const MeshExtended * mesh,                ///< Femus-mesh        (in)
      const MEDCoupling::MEDCouplingUMesh* support,  ///< med-mesh          (in)
      const int interface_id,                        ///< inrface identity  (in)
      const int order_cmp                            ///< order pt (1 or 2) (in)
  );
  // ============================================================================
  /// This function computes the mesh interface
  /// to set in the interface-function. The field in the
  /// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_nodeID(
      //   const MeshExtended * mesh,                ///< Femus-mesh        (in)
      const MEDCoupling::MEDCouplingUMesh* support,  ///< med-mesh          (in)
      const int interface_id,                        ///< inrface identity  (in)
      const int order_cmp                            ///< order pt (1 or 2) (in)
  );
  // ============================================================================
  /// This function computes the mesh interface
  /// to set in the interface-function. The field in the
  /// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_elemID(
      const MeshC& mesh_c,                           ///< Femus-mesh        (in)
      const MEDCoupling::MEDCouplingUMesh* support,  ///< med-mesh          (in)
      const int interface_id                         ///< inrface identity  (in)
  );
  // ============================================================================
};

#endif
#endif
