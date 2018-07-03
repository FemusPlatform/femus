#ifndef _EQNEXTENDED__
#define _EQNEXTENDED__

#include <map>
#include "MGEquationsSystem.h"
#include "Solverlib_conf.h"

#ifdef HAVE_MED
namespace MEDCoupling
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}
class InterfaceFunctionM;
#endif

class MeshExtended;
class BoundInterp;

// ==========================================================
// ==========================================================

class EquationSystemsExtendedM:public MGEquationsSystem
{


protected:
  const MeshExtended *_mg_mesh;	///< MultiGrid mesh
#ifdef HAVE_MED
  const MEDCoupling::MEDCouplingUMesh * _med_mesh;	///< MED mesh
  std::map < int, InterfaceFunctionM * >_interfaceFunMap;	///< Map where interface pointes are stored
#endif

public:

  EquationSystemsExtendedM (MGUtils & mgutils_in,	///<  MGUtils class       (in)
			    MeshExtended & mgmesh_in,	///<  MeshExtended class  (in)
			    MGFEMap & mgfemap_in,	///<  MGFEMap class       (in)
			    int npoint_data,	///<  point data          (in)
			    int ncell_data	///<  data cell           (in)
   );

  virtual ~ EquationSystemsExtendedM ();	///< Desctructor function


  void set_mesh_mg (MeshExtended & mg_mesh_femus_in)
  {
    _mg_mesh = &mg_mesh_femus_in;
    return;
  }
  const MeshExtended *getMeshMG () const
  {
    return _mg_mesh;
  }
// ==================  Mesh ===================================================
// ============================================================================
  // Print mat into case (extended only)
  void print_case_mat_h5 (const int t_init);
  void print_case_bc_h5 (const int t_init);

#ifdef HAVE_MED
  const MEDCoupling::MEDCouplingUMesh * getMeshMED () const
  {
    return _med_mesh;
  }

  void set_mesh_med (MEDCoupling::MEDCouplingUMesh & mg_mesh_med_in)
  {
    _med_mesh = &mg_mesh_med_in;
    return;
  }

  //! Function returning the pointer to interface with name <id>
  InterfaceFunctionM *get_interface_fun (int id	///< interface identity
    )
  {
    return _interfaceFunMap.find (id) ==
      _interfaceFunMap.end ()? NULL : _interfaceFunMap[id];
  }

// ============================================================================
  /// This function erases the interface with identity id
  void erase_interface_fun (int id	///< interface identity
    )
  {
    _interfaceFunMap.erase (id);
    return;
  }
//=============================================================================
  ///  This functon add an interface-function (med-mesh,femus-mesh, field)
  ///  in the interface-function boundary map (_interfaceFunMap):
  ///  the field in the interface-function is not assigned here
  void add_interface_fun (const int interface_name, const int bd_name_int,	///<[in] boundary id (int)     (in)
			  const MEDCoupling::MEDCouplingUMesh * b,	///<[in] med-submesh           (in)
			  const bool on_nodes,	///<[in] values on nodes (true)
			  const int order_cmp = 2	///<[in] order component (2=quad;lin=1)
    );

  //! This function adds an interface to interface map
  void add_interface_fun (const int interface_name,	///[in] Interface name
			  InterfaceFunctionM * fun	///[in] Pointer to interface
    );

  //! This function writes a field starting from an analytical expression
  void setBC (int InterfaceName,	 /**<[in] Interface ID */
	      int NumOfComp,		 /**<[in] Number of field components to write */
	      const char *AnalyticExpr	 /**<[in] String containing the analytical expression */
    );

  //! This function writes a field starting from a MED numerical field
  void setBC (int InterfaceName,					/**<[in] Interface ID */
	      int NumOfComp,						/**<[in] Number of field components to write */
	      const MEDCoupling::MEDCouplingFieldDouble * NumericField	/**<[in] MED field containing the values to write */
    );

  //! Set a MED field in system solution
  /*!
   * This function takes a MED field, stored into an interface, and write its
   * values into #MGSolBase::x_old numeric vector of the selected equation
   */
  void write_Boundary_value (int id_name,	///< [in] Identity interface name 
			     std::string mgsystem_name,	///< [in] System name          
			     int n_cmp,	///< [in] Number of components to write 
			     int first_cmp = 0	///< [in] First component   
    );

  //! Set a MED field as external field
  /*!
   * This function takes a MED field and sets it as an external field in #MGSolBase::_ExtField
   * The MED field values are not written into #MGSolBase::x_old numeric vector
   */
  void write_Boundary_value (std::string mgsystem_name,	///[in] Equation name
			     MEDCoupling::MEDCouplingFieldDouble * bcField	///[in] MED field to write
    );
  MEDCoupling::MEDCouplingFieldDouble * GetField(const std::string &systemName);
  //! Get values on mesh nodes
  /*!
   *  The function returns a MED field containing the solution of equation with name <system_name>.
   *  The values are relative to mesh nodes of interface with name <id>
   */
  MEDCoupling::MEDCouplingFieldDouble * getValuesOnBoundary_nodes (int id,	///< [in] Interface ID
								  const char *system_name,	///< [in] Equation name
								  int n_cmp,	///< [in] Number of components
								  int first_cmp = 0	///< [in] First component
    );

  //! Get values on mesh cells
  /*!
   *  The function returns a MED field containing the solution of equation with name <system_name>.
   *  The values are relative to mesh cells of interface with name <id>
   */
  MEDCoupling::MEDCouplingFieldDouble * getValuesOnBoundary_elem (int id,	///< [in] Interface ID
								 const char *system_name,	///< [in] Equation name
								 int n_cmp,	///< [in] Number of components
								 int first_cmp = 0	///< [in] First component
    );


  //! Solution of actual proc
  /*!
   *  This function returns a MED field containing the solution of equation with name <system_name>.
   *  The solution values are relative to the nodes belonging to proc mesh
   */
  MEDCoupling::MEDCouplingFieldDouble * getProcValues (MEDCoupling::MEDCouplingFieldDouble * NodeMap,	///<[in] MED field containing the map of proc nodes
						      int InterfaceId,	///<[in] Interface ID
						      const char *system_name,	///<[in] System name           
						      int n_cmp,	///<[in] Number of components to get
						      int first_cmp = 0,	///<[in] First component
                              int Level = 0
    );


  const MEDCoupling::MEDCouplingUMesh * getUMeshCoupling (int name);

  inline void getEqsNames (std::vector < string > &FieldNames)
  {
    get_eqs_names (FieldNames);
  };

  const MEDCoupling::MEDCouplingUMesh * getUMeshCoupling_orig (int name);
#endif


};				//end class 

#endif
