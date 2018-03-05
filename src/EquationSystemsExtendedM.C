
// std libraries -----------------
#include <iomanip>
#include <sstream>
#include <set>
#include <cstring>

#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"
#include "numeric_vectorM.h"
#include "MGFE_conf.h"
#include "MGUtils.h"

#ifdef HAVE_MED
#include "InterfaceFunctionM.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

// using namespace libMesh;
using namespace ParaMEDMEM;
#endif
// ========================================================================
/// Constructor
EquationSystemsExtendedM::EquationSystemsExtendedM(
  MGUtils& mgutils_in,
  MeshExtended& mgmesh_in,
  MGFEMap& mgfemap_in,
  int npoint_data,
  int ncell_data
) :
  MGEquationsSystem(mgutils_in,mgmesh_in,mgfemap_in,npoint_data,ncell_data),
  _mg_mesh(&mgmesh_in) {  // mesh class in

}

// ========================================================================
EquationSystemsExtendedM::~EquationSystemsExtendedM() {
#ifdef HAVE_MED
  std::map<int, InterfaceFunctionM *>::iterator it = _interfaceFunMap.begin();
  for(; it != _interfaceFunMap.end(); it++) {
    if(it->second) { delete it->second; }
    it->second = NULL;
  }
#endif
  return;
}

// ============================================================================
void EquationSystemsExtendedM::print_case_mat_h5(const int t_init) {
  if(_mg_mesh->_iproc==0) {
    // file ---------------------------------------
    const int NoLevels  = (int)(_mgutils._geometry["nolevels"]);
    const int ndigits   = stoi(_mgutils._sim_config["ndigits"]);
    std::string output_dir = _mgutils._inout_dir;
    std::string   basecase = _mgutils.get_file("BASECASE");
    std::ostringstream filename; // file name
    filename << output_dir << basecase << "." << setw(ndigits) << setfill('0') << t_init << ".h5";
    hid_t   file =H5Fopen(filename.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
    // print mat
    _mg_mesh->print_mat_hf5(filename.str(),"COLOR");

    H5Fclose(file);
  }
  return;
}

// ============================================================================
void EquationSystemsExtendedM::print_case_bc_h5(const int t_init) {


  if(_mg_mesh->_iproc==0) {
    // file ---------------------------------------
//     const int NoLevels  = _mgutils.get_par("nolevels");
//     const int ndigits   = _mgutils.get_par("ndigits");
    
     const int NoLevels  = (int)(_mgutils._geometry["nolevels"]);
    const int ndigits   = stoi(_mgutils._sim_config["ndigits"]);
    std::string output_dir = _mgutils._inout_dir;
    std::string   basecase = _mgutils.get_file("BASECASE");
    std::ostringstream filename; // file name
    filename << output_dir << basecase << "." << setw(ndigits) << setfill('0') << t_init << ".h5";
    hid_t   file =H5Fopen(filename.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
    // print mat
    _mg_mesh->print_bc_hf5(filename.str(),"BC");

    H5Fclose(file);
  }
  return;
}

#ifdef HAVE_MED
// =======================================================================
void EquationSystemsExtendedM::setBC(
  int name,
  int  n_cmp,
  const MEDCouplingFieldDouble *field
) {// =======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);  // interface-function
  if(fct == NULL) { return; }
  fct->set_field(field);
#ifdef PRINT_MED
//   fct->printOn(std::cout,name);
#endif

  return;
}

// =========================================================================
void EquationSystemsExtendedM::setBC(
  int name,
  int n_cmp,
  const char *s
) { // ======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);
  if(fct == NULL) { return; }
  if(name>9) {
    fct->set_analytic_field(s, n_cmp);
  } else {
    fct->set_analytic_field_elem(s, n_cmp);
  }

#ifdef  PRINT_MED
//  fct->printOn(std::cout,name);
#endif

  return;
}


// ============================================================================
// ===============  Interface function routines ===============================

// ============================================================================
// This functon add an interface-function (med-mesh,femus-mesh, field)
//  in the interface-function boundary map (_interfaceFunMap):
//  the field in the interface-function is not assigned here
void EquationSystemsExtendedM::add_interface_fun(       ///<map position(return)
  const int  interface_name,
  const int  interface_id,    ///< boundary id (int)  (in)
  const ParaMEDMEM::MEDCouplingUMesh * b,  ///< med-submesh        (in)
//   const int /*from_cmp*/,               ///< initial id component
//   const int /*n_cmp*/,                  ///< n components
  const bool on_nodes,                     ///< values on nodes (true)
  const int order_cmp                      ///< order component (2=quad;lin=1)
) { // ========================================================================

  // new function
  InterfaceFunctionM *fun = new InterfaceFunctionM;
  // set the mesh interface
  if(on_nodes) {
    fun->set_mesh_interface_nodeID(_mg_mesh,b,interface_id,order_cmp);
  } else { // on elements (mat)
    fun->set_mesh_interface_elemID(_mg_mesh,b,interface_id);
  }
  // setting the fun in the map _interfaceFunMap at bd_name_id
  _interfaceFunMap[interface_name] = fun;  // added to the map
  fun->set_order(order_cmp);
  // print
#ifdef PRINT_MED
  std::cout << " Added  InterfaceFunction on interface "<< interface_id << "\n";
  fun->printOn(std::cout, interface_id);
#endif
  return;
}

// ============================================================================
// This functon add an interface-function (med-mesh,femus-mesh, field)
//  in the interface-function boundary map (_interfaceFunMap):
//  the field in the interface-function is not assigned here
void EquationSystemsExtendedM::add_interface_fun(       ///<map position(return)
  const int  interface_name,
  InterfaceFunctionM *fun
) { // ========================================================================
  _interfaceFunMap[interface_name] = fun;  // added to the map
#ifdef PRINT_MED
  std::cout << " Added  InterfaceFunction on interface "<< interface_id << "\n";
  fun->printOn(std::cout, interface_id);
#endif
  return;
}

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtendedM::getValuesOnBoundary_nodes(
  int                name,           // int boundary identity   (in)
  const char *system_name,           // system name             (in)
  int               n_cmp,           //  first variable         (in)
  int           first_cmp            // n variables             (in)
)  {// ========================================================================
  // from LibMesh function fct
  InterfaceFunctionM * fct = get_interface_fun(name);
  if(fct == NULL) { return NULL; }
  int nNodes = fct->getSupport()->getNumberOfNodes();
  int n_nodes_mg  = fct->get_n();
  int * map_mg    = fct->get_map_mg();
  int * map_med   = fct->get_map_med();

  // from MGMesh
  int Level=_mg_mesh->_NoLevels-1;
  int offset=_mg_mesh->_NoNodes[Level];

  // from  MGSystem* system;
  MGSolBase * mgsyst=get_eqs(system_name);

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());  f->setName(system_name);  f->setNature(ConservativeVolumic);
  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();  array->alloc(nNodes,n_cmp);  array->fillWithZero();
  int order = fct->get_order();
  int Dim   = fct->getSupport()->getSpaceDimension();
  int npt_elem=fct->getSupport()->getNumberOfNodesInCell(0);  // # of points in one element 
  int Fcc=npt_elem; // # of vertices (linear) in one element
  
  string str(system_name);
  for(int ji=0; ji<n_cmp; ji++) { // --------------------------------------------------------------
    if(((str=="NS0" || str=="FSI0" ||str=="FSIA0"||str=="NSA0")) && ji+first_cmp == Dim) order=1;

    if(order==2) {
      for(int i_mg=0; i_mg < nNodes; i_mg++) {
        int node_mg   = map_mg[i_mg];  // mg  node
        int node_med  = map_med[i_mg];  // med node
        const int kdof_top = mgsyst->_node_dof[Level][node_mg+(ji+first_cmp)*offset];
        double v = (*(mgsyst->x_old[Level]))(kdof_top);
        array->setIJ(node_med,ji,v);
      }
    } 
    else if(order==1) {
      switch(npt_elem) {
      case(3):
        Fcc=2; 
        break;
      case(9):
         Fcc=4; 
        break;
      case(27):
        Fcc=8; 
        break;
      default:
        std::cout<<"\033[0;34m++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\033[0m\n";
        break;
      }
      std::map <int, int> MedMg;
      for(int nMed = 0; nMed < fct->getSupport()->getNumberOfNodes(); nMed++) {
        const int MedNumber = map_med[nMed];
        const int MgNumber  = map_mg[nMed];
        MedMg[MedNumber] = MgNumber;
      }
      std::vector<int> conn;
      for(int i_mg=0; i_mg < fct->getSupport()->getNumberOfCells(); i_mg++) {
        fct -> getSupport() -> getNodeIdsOfCell(i_mg,conn);
        for(int i_node = 0; i_node<Fcc; i_node ++) {
          const int kdof_top  = mgsyst-> _node_dof[Level][MedMg[conn[i_node]] +(ji+first_cmp) * offset];
          double intvar      = (*(mgsyst->x_old[Level]))(kdof_top);
          array->setIJ(conn[i_node],ji,intvar);
        }
        conn.clear();
      }
    }
  }
  for(int kj=0; kj<n_cmp; kj++) { array->setInfoOnComponent(kj,std::to_string(kj+first_cmp)); }  // set info -> array
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//  fct->printOn(std::cout,name);   // print fct

  return f;

}


// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtendedM::getValuesOnBoundary_elem(
  int                  id,           // int boundary identity   (in)el_conn
  const char *system_name,           // system name             (in)
  int               n_cmp,           //  first variable       (in)
  int           first_cmp            // n variables       (in)
)  {// ========================================================================
  // from LibMesh function fct
  InterfaceFunctionM * fct = get_interface_fun(id);
  if(fct == NULL) { return NULL; }
  int nNodes = fct->getSupport()->getNumberOfNodes();
//   std::map<int, int> & NodeID = fct->NodeID();
  int n_nodes_mg = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();

  // from MGMesh
  int Level=_mg_mesh->_NoLevels-1;
  int offset=_mg_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * mgsyst=get_eqs(system_name);
//   int nComp =n_variable;
  std::cerr << "\n EquationSystemsExtendedM nComp = " << n_cmp << std::endl;

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());
  f->setName(system_name);

  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,n_cmp);
  // filling array
  for(int i_mg=0; i_mg < n_nodes_mg; i_mg++) {
    int node_mg   = map_mg[i_mg];  // mg  node
    int node_med  = map_med[i_mg];  // med node

    for(int j=  first_cmp; j< first_cmp+n_cmp; j++) {
      const int kdof_top = mgsyst->_node_dof[Level][node_mg+j*offset];
      double v = (*(mgsyst->x_old[Level]))(kdof_top);
      array->setIJ(node_med,j,v);
    }
  }
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//   fct->printOn(std::cout,id);   // print fct

  return f;

}


// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtendedM::getProcValues(
  ParaMEDMEM::MEDCouplingFieldDouble *NodeMap,
  int InterfaceId,                                                               
  const char *system_name,           // system name             (in)
  int               n_cmp,           //  first variable       (in)
  int           first_cmp,            // n variables       (in)
  int Level2
)  {// ========================================================================


  InterfaceFunctionM * fct = get_interface_fun(InterfaceId);
  if(fct == NULL) { return NULL; }

//   std::map<int, int> & NodeID = fct->NodeID();
  int n_nodes_mg = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();

  // from MGMesh
  
  int Level=_mg_mesh->_NoLevels-1;  
  int offset=_mg_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * mgsyst=get_eqs(system_name);

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(NodeMap->getMesh());
  f->setName(system_name);

  int nNodes = NodeMap->getMesh()->getNumberOfNodes();
  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,n_cmp);
  // filling array
  
 
  
  for(int i_mg=0; i_mg < nNodes; i_mg++) {
    int node_mg   = map_mg[(int)NodeMap->getIJ(i_mg,0)];  // mg  node
    node_mg   =(int)NodeMap->getIJ(i_mg,0);  // mg  node
    
//     int node_med  = map_med[i_mg];  // med node

    

    for(int j=  first_cmp; j< first_cmp+n_cmp; j++) {
      const int kdof_top = mgsyst->_node_dof[Level][node_mg+j*offset];
      double v = (*(mgsyst->x_old[Level]))(kdof_top);
//       if(i_mg== 589 || i_mg == 594 || i_mg == 593 || i_mg == 751 || i_mg == 756 || i_mg == 755 || i_mg == 749 || i_mg == 754 || i_mg == 753)
//         std::cout<<"mg node "<<node_mg<<" for med node "<<i_mg<<" out of "<<nNodes<<" nodes, value "<<v<<std::endl;
      
      array->setIJ(i_mg,j,v);
    }
  }
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//   fct->printOn(std::cout,id);   // print fct

  return f;

}


// =========================================================================
void EquationSystemsExtendedM::write_Boundary_value(
  int id_boundary_name,       /**< identity interface name (in)       */
  std::string mgsystem_name,  /**< system name             (in)       */
  int n_cmp,                  /**<  variable system    (in)           */
  int first_cmp               /**< from variable system      (in)     */
) { // ======================================================================


  InterfaceFunctionM * fct = get_interface_fun(id_boundary_name);
  if(fct == NULL) { return; }


  // interface support (mesh) fct -----------------------------------
  const ParaMEDMEM::MEDCouplingUMesh *support= fct->getSupport();
  int nNodes_med = support->getNumberOfNodes();
  int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();

  // MGSolver -------------------------------------------------------
  MGSolBase * mgsyst=get_eqs(mgsystem_name.c_str());
  int Level=_mg_mesh->_NoLevels-1;
  int offset= _mg_mesh->_NoNodes[Level];
  NumericVectorM *old_sol_top= mgsyst->x_old[Level];
  if(mgsystem_name=="NS2P") old_sol_top= mgsyst->x_nonl[Level];                      // MOD
  // Mesh -----------------------------------------------------------
  MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
  std::vector<int> mat=ext_mesh->_mat_id;
  std::vector<int> bc_id=ext_mesh->_bc_id;

  double bc_value[DIMENSION + 1];
  int order = fct->get_order();
  int NodesPerCell = fct->getSupport()->getNumberOfNodesInCell(0);
  int Dim = fct->getSupport()->getSpaceDimension();

  // setting the boundary values --------------------------------------------------------------------
  for(int kvar=0; kvar<n_cmp; kvar++) {

    if(((mgsystem_name=="NS0" || mgsystem_name=="FSI0" ||mgsystem_name=="FSIA0"||mgsystem_name=="NSA0")) && kvar+first_cmp == Dim) { order=1; }

    if(order==2) { // --------------- order 2 -------------------------------------------------
      for(int i_mg=0; i_mg < n_nodes_mg; i_mg++) { // -----------------------------------
        int gl_node_bd = map_mg[i_mg];  // mg  node
        int node_med  = map_med[i_mg];  // med node

        fct->eval(node_med,n_cmp, bc_value);
        const int kdof_top = mgsyst->_node_dof[Level][gl_node_bd+(kvar+first_cmp)*offset];
        old_sol_top->set(kdof_top , bc_value[kvar]);
      } // ------------------------------------------------------------------------------
    }  // --------------- order 2 --------------------------------------------------------------
    else if(order==1) { // --------------- order 1 --------------------------------------------------------------------
      int Fcc;
      switch(NodesPerCell) {
      case(3):
        Fcc = 2;
        break;
      case(9):
        Fcc = 4;
        break;
      case(27):
        Fcc = 8;
        break;
      default:
        std::cout<<"\033[0;34m+++++++++++++++++++++++++++++++++++++++++++++++++\033[0m\n";
        break;
      }

      std::map <int, int> MedMg;                                    // From med numbering to mg numbering
      for(int nMed = 0; nMed < fct->getSupport()->getNumberOfNodes(); nMed++) {
        const int MedNumber = map_med[nMed];
        const int MgNumber  = map_mg[nMed];
        MedMg[MedNumber] = MgNumber;
      }
      std::vector<int> conn;
      for(int i_mg=0; i_mg < fct->getSupport()->getNumberOfCells(); i_mg++) {                     // Loop over the boundary elements
        fct->getSupport()->getNodeIdsOfCell(i_mg,conn);           // Med numbering of element nodes
        // Loop for the calculation of the field on the gauss point
        for(int i_node = 0; i_node<Fcc; i_node ++) {
          fct->eval(conn[i_node],n_cmp, bc_value);
          const int kdof_top = mgsyst->_node_dof[Level][MedMg[conn[i_node]]+(kvar+first_cmp)*offset];
          old_sol_top->set(kdof_top , bc_value[kvar]);
        }
        conn.clear();
      }
      MedMg.clear();
    } // --------------- order 1 --------------------------------------------------------------------
  }   // end setting boundary values (kvar) -------------------------------------------------------------
  return;
}



void EquationSystemsExtendedM::write_Boundary_value(
  std::string mgsystem_name,  /**< system name             (in)       */
  ParaMEDMEM::MEDCouplingFieldDouble * bcField
) { // ======================================================================
  // MGSolver -------------------------------------------------------
  MGSolBase * mgsyst=get_eqs(mgsystem_name.c_str());
  mgsyst->_ExtField = bcField->deepCpy();
  return;
}

  ParaMEDMEM::MEDCouplingFieldDouble *  EquationSystemsExtendedM::GetField(const std::string &systemName){
    MGSolBase * mgsyst=get_eqs(systemName.c_str());
    ParaMEDMEM::MEDCouplingFieldDouble * ExtField = mgsyst->_ExtField->deepCpy();
    return ExtField;
  }

//get actual support
const MEDCouplingUMesh* EquationSystemsExtendedM::getUMeshCoupling(
  int name
) {// =======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);  // interface-function
//   if(fct == NULL) return;
  const MEDCouplingUMesh *sourceMesh = fct->getSupport();
  return sourceMesh;
}
//get original support
const MEDCouplingUMesh* EquationSystemsExtendedM::getUMeshCoupling_orig(
  int name
) {// =======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);  // interface-function
//   if(fct == NULL) return;
  const MEDCouplingUMesh *sourceMesh = fct->getSupport_orig();
  return sourceMesh;
}

#endif
// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
