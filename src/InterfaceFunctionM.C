#include "InterfaceFunctionM.h"

#ifdef HAVE_MED

#include "MGFE_conf.h"
#include "MeshExtended.h"

#include <iomanip>
#include <limits>


#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"







// using namespace MEDCoupling;
// using namespace libMesh;
// ======================================================
InterfaceFunctionM::~InterfaceFunctionM() {
  delete [] _map_mg; delete [] _map_med;
  if(_field) _field->decrRef(); if(_support_med) _support_med->decrRef();
    if(_support_med_orig) _support_med_orig->decrRef();
  delete [] _mesh_mg;
}

// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionM::set_mesh_BFinterface_nodeID(
  const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id,                      ///< inrface identity  (in)
  const int order_cmp                          ///< order pt (1 or 2) (in)
) { // ========================================================================

  double toll_dist=1.e-12;
//   _nodeID.clear();  // clear map node: FEMus-mesh -> MED-mesh

  // Femus-mesh
  _mesh_mg = mesh;                                  // MG mesh
  int nlevels_mg=_mesh_mg->_NoLevels-1;             // top MG level
  int n_nodes_mg=_mesh_mg->_NoNodes[nlevels_mg];    // top MG n nodes
  int dim_mg= _mesh_mg->_dim;                       // MG dimension

  std::vector<double> xyz_mg;
  std::vector<int> nodeid_mg;

  // Med-mesh
  _support_med = support;                 // MED-mesh
  _support_med_orig = support;                 // MED-mesh
  _n= _support_med->getNumberOfNodes();   // n MED-nodes
  const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples(); // n MED-tuple-cords

//   std::vector<std::pair<double, int > > dist_med;
//   std::vector<std::pair<double, int > >::iterator itdist;
//
//
//   // MED-distance *************************************************************
//   dist_med.resize(n_nodes_med);
//   for(int i_med=0; i_med<n_nodes_med; i_med++) { // i_med= MED-node
//     double dd=0.; // dd=MED-point istance from (0.0)
//     for(int idim=0; idim<dim_mg; idim++) dd  += d->getIJ(i_med, idim)*d->getIJ(i_med, idim);
//     dist_med[i_med] = std::make_pair(dd,i_med);
//   }
//   std::sort(dist_med.begin(), dist_med.end());
//   // **************************************************************************

  // MG-node identity vector **************************************************
  int n_bd_nodes=0; // counter of mg mesh nodes on the boundary interface (interface_id)
  switch(order_cmp) {
  case 2:

    //build xyz_mg and nodeid_mg finding nodes on boundary interface (interface_id) from mg mesh
    for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
      if(fabs(mesh->_bc_id[i]) == interface_id) {
        n_bd_nodes++;
        for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
        nodeid_mg.push_back(i);
      }
    }

    break;
  case 1:
    std::cout<< " \n InterfaceFunctionM::set_mesh_BFinterface_nodeID:"<<
             " order not implemented \n";
//     xyzToSetOfNodes_onBd_linM(dim_mg,n_nodes_mg,mesh->_xyz,mesh->_bc_id,s2,n2);
    break;
  default:
    std::cout<< " \n InterfaceFunctionM::set_mesh_BFinterface_nodeID:"<<
             " order not implemented \n";
    abort();
  } // ************************************************************************


  _map_med = new int [n_bd_nodes];
  _map_mg =  new int [n_bd_nodes];
  _n     =           n_bd_nodes;


  // Building map _nodeID[Mg-node]-> MED-node *********************************
  for(int i_mg=0; i_mg<n_bd_nodes; i_mg++) {// i_mg=MGnode

//     /// MG node distance (dist_mg)
//     double dist_mg=0.;
//     for(int idim=0; idim<dim_mg; idim++)
//       dist_mg  += xyz_mg[i_mg*dim_mg+idim]*xyz_mg[i_mg*dim_mg+idim];
//     // Find the first node (a) in the interval search [a,b]
//     itdist = dist_med.begin();
//     while(fabs(dist_mg - itdist->first) >  toll_dist)  ++itdist;
//     if(itdist == dist_med.end()) {
//       std::cout << "non ci sono punti con toll_dist" << toll_dist << '\n';
//       abort();
//     }
    // finding node (k_min)  with mimimum distance (min)
    int k_min=-1; double min_dist=1.e+15;

    for(int i_med=0; i_med<n_nodes_med; i_med++) { // i_med= MED-node
//     for(; itdist<dist_med.end(); itdist++) { // node med-interface loop
      double dd=0.;// dd=distance(MED-point,MG-point)
      for(int idim=0; idim<dim_mg; idim++)
        dd  += fabs(xyz_mg[i_mg*dim_mg+idim] - d->getIJ(i_med, idim));
      if(dd < min_dist)  {min_dist=dd; k_min=i_med;}
    }
    _map_med[i_mg]=k_min;
    _map_mg[i_mg] =nodeid_mg[i_mg];
//     _nodeID[nodeid_mg[i_mg]] =k_min; // filling the map _nodeID
  } // ************************************************************************

  // clean
  //  std::cout << "==================== mesh node id  =============================== "<<  interface_id <<"\n";
//   std::map<int, int>::const_iterator itN = _nodeID.begin();
//   for(; itN != _nodeID.end(); itN++) {
//     mesh->_node_id[itN->first]=itN->second;
//   }
//  printOn(std::cout,interface_id);
  xyz_mg.clear(); nodeid_mg.clear();
//   dist_med.clear();

  return;
}

  void InterfaceFunctionM::set_maps( ///< Setting the interface-meshes 
  int *map_med,   ///< node/element map index->MGMesh      (size: _n)  
  int *map_mg,    ///< node/element map index->MEDCoupling (size: _n)
  const int nodes    
  ) {
      _map_med = new int[nodes];
      _map_mg = new int[nodes];
      
      for(int i = 0; i < nodes; ++i){
         _map_med[i] = map_med[i];   ///< node/element map index->MGMesh      (size: _n)  
         _map_mg[i]  = map_mg[i];
      }
      return;    
  }

// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionM::set_mesh_interface_nodeID(
  const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id,                      ///< inrface identity  (in)
  const int order_cmp                          ///< order pt (1 or 2) (in)
) { // ========================================================================

  double toll_dist=1.e-10;
//   _nodeID.clear();  // clear map node: FEMus-mesh -> MED-mesh

  // Femus-mesh
  _mesh_mg = mesh;                                  // MG mesh
  int nlevels_mg=_mesh_mg->_NoLevels-1;             // top MG level
  int n_nodes_mg=_mesh_mg->_NoNodes[nlevels_mg];    // top MG n nodes
  int dim_mg= _mesh_mg->_dim;                       // MG dimension

  std::vector<double> xyz_mg;  std::vector<int> nodeid_mg;

  // Med-mesh
  _support_med = support;                 // MED-mesh
    _support_med_orig = support;                 // MED-mesh
//   _n= _support_med->getNumberOfNodes();   // n MED-nodes
  const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples(); // n MED-tuple-cords

  std::vector<std::pair<double, int > > dist_med;
  std::vector<std::pair<double, int > >::iterator itdist;


  // MED-distance *************************************************************
  dist_med.resize(n_nodes_med);
  for(int i_med=0; i_med<n_nodes_med; i_med++) { // i_med= MED-node
    double dd=0.; // dd=MED-point istance from (0.0)
    for(int idim=0; idim<dim_mg; idim++) dd  += d->getIJ(i_med, idim)*d->getIJ(i_med, idim);
    dist_med[i_med] = std::make_pair(dd,i_med);
  }
  std::sort(dist_med.begin(), dist_med.end());
  // **************************************************************************

  // MG-node identity vector **************************************************
  int n_bd_nodes=0; // counter of mg mesh nodes on the boundary interface (interface_id)
  switch(order_cmp) {
  case 2:

    //build xyz_mg and nodeid_mg finding nodes on boundary interface (interface_id) from mg mesh
    for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
      if(fabs(mesh->_bc_id[i]) == interface_id || interface_id < 10) {
        n_bd_nodes++;
        for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
        nodeid_mg.push_back(i);
      }
    }
    break;
    
  case 1:
    
    //build xyz_mg and nodeid_mg finding only linear nodes on boundary interface (interface_id) from mg mesh
    for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
      if(mesh->_bc_id[i] == -interface_id) {
        n_bd_nodes++;
        for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
        nodeid_mg.push_back(i);
      }
    }
    break;
    
   default:
    std::cout<< " \n InterfaceFunctionM::set_mesh_BFinterface_nodeID:"<<
             " order not implemented \n";
    abort();
  } // ************************************************************************


  _map_med = new int [n_bd_nodes];
  _map_mg =  new int [n_bd_nodes];
  _n     =           n_bd_nodes;


  // Building map _nodeID[Mg-node]-> MED-node *********************************
  for(int i_mg=0; i_mg<n_bd_nodes; i_mg++) {// i_mg=MGnode

    /// MG node distance (dist_mg)
    double dist_mg=0.;
    for(int idim=0; idim<dim_mg; idim++)
      dist_mg  += xyz_mg[i_mg*dim_mg+idim]*xyz_mg[i_mg*dim_mg+idim];
    // Find the first node (a) in the interval search [a,b]
    itdist = dist_med.begin();
    while(fabs(dist_mg - itdist->first) >  toll_dist)  ++itdist;
    if(itdist == dist_med.end()) {
      std::cout << "non ci sono punti con toll_dist" << toll_dist << '\n';
      abort();
    }
    // finding node (k_min)  with mimimum distance (min)
    int k_min=-1; double min_dist=1.e+15;
    for(; itdist<dist_med.end(); itdist++) { // node med-interface loop
      double dd=0.;// dd=distance(MED-point,MG-point)
      for(int idim=0; idim<dim_mg; idim++)
        dd  += fabs(xyz_mg[i_mg*dim_mg+idim] - d->getIJ(itdist->second, idim));
      if(dd < min_dist)  {min_dist=dd; k_min=itdist->second;}
      if(fabs(dist_mg - itdist->first) >  toll_dist) break;
    }
    _map_med[i_mg]=k_min;
    _map_mg[i_mg] =nodeid_mg[i_mg];
//     _nodeID[nodeid_mg[i_mg]] =k_min; // filling the map _nodeID
  } // ************************************************************************

  // clean
   
//   std::map<int, int>::const_iterator itN = _nodeID.begin();
//   for(; itN != _nodeID.end(); itN++) {
//     mesh->_node_id[itN->first]=itN->second;
//   }
#ifdef PRINT_MED
 std::cout << "====interface_id= "<<  interface_id << "with nodes" <<  n_bd_nodes << "\n";
//  printOn(std::cout,interface_id);
#endif
  xyz_mg.clear(); nodeid_mg.clear();
  dist_med.clear();

  return;
}





// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionM::set_mesh_interface_elemID(
  const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id                      ///< inrface identity  (in)
) { // ========================================================================

  double toll_dist=1.e-4;
  double x_m[DIMENSION];//tmp buffer

  // Femus-mesh
  _mesh_mg = mesh;                                  // MG mesh
  int Level=_mesh_mg->_NoLevels-1;                  // top MG level
  int n_nodes_mg=_mesh_mg->_NoNodes[Level];         // top MG level  n nodes
  int n_elements_mg=_mesh_mg->_NoElements[0][Level];// n MG elements
  int dim_mg= _mesh_mg->_dim;                       // MG dimension
  int n_nod_el=(dim_mg==DIMENSION)?NDOF_FEM:NDOF_FEMB;

  int ndof_lev=0;
  for(int ilev=0; ilev < Level+1; ilev++) {
    ndof_lev += _mesh_mg->_NoElements[0][ilev];
  }


  double *xyz_mg  = new double [ndof_lev*dim_mg];// MG el centers
  int *elemid_mg = new int [ndof_lev*dim_mg];    // MG el identity

  // Med-mesh
  _support_med = support;                               // MED-mesh
    _support_med_orig = support;                               // MED-mesh
//   _n= _support_med->getNumberOfNodes();                 // n MED-nodes
  int n_elements_med= _support_med->getNumberOfCells(); // n MED-elements
  const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples();             //  MED-nodes
  int dim_med = d->getNumberOfComponents();             //  MED dimension

  double *xyz_med = new double [n_elements_med*dim_mg]; // MED el centers
  std::vector<std::pair<double, int > > dist_med;       // MED el distance
  std::vector<std::pair<double, int > >::iterator itdist;

//   _elemID.clear();  // clear map elements: FEMus-mesh -> MED-mesh
  std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
    for(int idim=0; idim<dim_mg; idim++) x_m[idim]=0.; // zeros
    _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    for(int inode=0; inode< n_nod_el; inode++) {
      for(int idim=0; idim<dim_mg; idim++)
        x_m[idim] += d->getIJ(nodes1[inode], idim);
    } // end inode
    nodes1.clear(); // clear element node vector
    for(int idim=0; idim<dim_mg; idim++)
      xyz_med[dim_mg*ielem+idim]=x_m[idim]/n_nod_el;
  }// *************************************************************************

  // MED-distance from center (dist_med) **************************************
  // build distance between element center and (0,0,0)
  dist_med.resize(n_elements_med);
  for(int i_med=0; i_med<n_elements_med; i_med++) {
    double dd=0.;
    for(int idim=0; idim<dim_mg; idim++)
      dd  += xyz_med[dim_mg*i_med+idim]*xyz_med[dim_mg*i_med+idim];
    dist_med[i_med] = std::make_pair(dd,i_med);
  }
  std::sort(dist_med.begin(), dist_med.end());
  // **************************************************************************

  // Computing the MG center (xyz_mg) *****************************************
  int n_mat_elem=0;
  for(int iproc=0; iproc<_mesh_mg->_n_subdom; iproc++) {
    int nel_b=_mesh_mg->_off_el[0][_mesh_mg->_NoLevels*iproc+Level];
    int nel_e=_mesh_mg->_off_el[0][Level+_mesh_mg->_NoLevels*iproc+1];
    for(int iel=nel_b; iel <nel_e; iel++) {
      if(fabs(mesh->_mat_id[iel]) == interface_id) {
        // average values
        for(int idim=0; idim<dim_mg; idim++) x_m[idim]=0.; // zeros
        for(int  inode=0; inode<n_nod_el; inode++)    {
          int el_conn = _mesh_mg->_el_map[0][iel*n_nod_el+inode];
          for(int  idim=0; idim<dim_mg; idim++)
            x_m[idim] += _mesh_mg->_xyz[el_conn+idim*n_nodes_mg];
        }
        for(int idim=0; idim<dim_mg; idim++)
          xyz_mg[n_mat_elem*dim_mg+idim]=x_m[idim]/n_nod_el;
        elemid_mg[n_mat_elem]=iel;
        n_mat_elem++;
      }
    }// ---- end iel -------------------------
  } // ************************************************************************

  _n=n_mat_elem;
  _map_med = new int [n_mat_elem];
  _map_mg =  new int [n_mat_elem];
//
//
  std::cout  << '\n' << "n mat eleme " << n_mat_elem << '\n';
  std::cout  << '\n' << "n n_elements_med " << n_elements_med << '\n';
  // Building _elemID *********************************************************
  for(int i_mg=0; i_mg<n_mat_elem; i_mg++) {
    double dist_mg=0.; // MG element distance ---------------------------------
    for(int idim=0; idim<dim_mg; idim++)
      dist_mg  += xyz_mg[i_mg*dim_mg+idim]*xyz_mg[i_mg*dim_mg+idim];
    // finding the range [a= dist_mg  <  toll_dist+ itdist->first -------------
    itdist = dist_med.begin();
    while(fabs(dist_mg - itdist->first) >  toll_dist && itdist != dist_med.end())  itdist++;
    if(itdist == dist_med.end()) {
      std::cout << "non ci sono punti con toll_dist" << toll_dist << '\n';
      std::cout << xyz_mg[i_mg*dim_mg] << " " << xyz_mg[i_mg*dim_mg+1] << " " <<  xyz_mg[i_mg*dim_mg+2]<< '\n';
      abort();
    }
    // finding the match for the element center -------------------------------
    int k_min=-1; double min=1.e+15;
    for(; itdist<dist_med.end(); itdist++) {
      double dd=0.;
      for(int idim=0; idim<dim_mg; idim++)
        dd  += fabs(xyz_mg[i_mg*dim_mg+idim] - xyz_med[itdist->second*dim_mg+idim]);//d->getIJ(itdist->second, idim));
      if(dd < min)  {min=dd; k_min=itdist->second;}
      if(fabs(dist_mg - itdist->first) >  toll_dist) break;
    }
    _map_med[i_mg] = k_min;
    _map_mg[i_mg]  = elemid_mg[i_mg];
  } // ************************************************************************


  //  clean *******************************************************************
//   std::cout << "==================== mesh elem id  =============================== "<<  interface_id <<"\n";
//   std::map<int, int>::const_iterator itN = _elemID.begin();
//   for(; itN != _elemID.end(); itN++) {
// //     mesh->_elem_id[itN->first]=itN->second;
//     std::cout << "("<<itN->first << "," <<  itN->second << ") ";
//   }
//   std::cout << '\n';
//   printOn(std::cout, interface_id);
  delete [] xyz_med;   delete [] xyz_mg;
  delete []elemid_mg ; dist_med.clear();
//   d->decrRef();
  return;
}



// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void InterfaceFunctionM::set_analytic_field(
  const char *symbolic_eq,     // symbolic function
  int nComp             // number of componenents
) { // ========================================================================
  if(_field) _field->decrRef();
  //ON CELLS does not work (due to baryc in quad 9 or hex 27)
  MEDCoupling::TypeOfField type = MEDCoupling::ON_NODES;
  int dim=_support_med->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if(dim > 0)    vars[0] = "x";
  if(dim > 1)    vars[1] = "y";
  if(dim > 2)    vars[2] = "z";
  if(symbolic_eq==std::string("")) {
    std::cout<<
             "InterfaceFunctionM::set_analytic_field: NULL function"; abort();
  }
  _field = _support_med->fillFromAnalyticNamedCompo(type, nComp, vars, symbolic_eq);
  _field->setName(symbolic_eq);
  _field->checkConsistencyLight();
//   std::ostream out("prova.field");
  std::cout << "InterfaceFunctionM::set_analytic_field_interface \n";
//   printOn(std::cout);
  return;
}

// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void InterfaceFunctionM::set_analytic_field_elem(
  const char *symbolic_eq,     // symbolic function
  int nComp             // number of componenents
) { // ========================================================================
  if(_field) _field->decrRef();
  //ON CELLS does not work (due to baryc in quad 9 or hex 27)
  MEDCoupling::TypeOfField type = MEDCoupling::ON_NODES;
  int dim=_support_med->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if(dim > 0)    vars[0] = "x";
  if(dim > 1)    vars[1] = "y";
  if(dim > 2)    vars[2] = "z";
  if(symbolic_eq==std::string("")) {
    std::cout<<
             "InterfaceFunctionM::set_analytic_field_elem: NULL function"; abort();
  }



  const MEDCoupling::MEDCouplingFieldDouble * field_nodes =
    _support_med->fillFromAnalyticNamedCompo(type, nComp, vars, symbolic_eq);

  int n_elements_med= _support_med->getNumberOfCells(); // n MED-elements
//   const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
  int n_nodes_med = field_nodes->getNumberOfTuples();             //  MED-nodes
//   int dim = field_nodes->getNumberOfComponents();             //  MED dimension

  int n_nod_el=NDOF_FEM;
  double * sum=new double[nComp];

_field= MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::NO_TIME);
 _field->setMesh(_support_med);

  MEDCoupling::DataArrayDouble *array;
  array=MEDCoupling::DataArrayDouble::New();
  array->alloc(n_elements_med,nComp);

  double *field_av = new double [n_elements_med*nComp]; // MED el centers
//   std::vector<std::pair<double, int > > dist_med;       // MED el distance
//   std::vector<std::pair<double, int > >::iterator itdist;

//   _elemID.clear();  // clear map elements: FEMus-mesh -> MED-mesh
  std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
    for(int icomp=0; icomp<nComp; icomp++) sum[icomp]=0;
    _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    for(int icomp=0; icomp<nComp; icomp++){ 
    for(int inode=0; inode< n_nod_el; inode++)  sum[icomp] += field_nodes->getIJ(nodes1[inode],icomp);
    field_av[ielem*nComp+icomp]=sum[icomp]/n_nod_el;
  }
    nodes1.clear(); // clear element node vector
  }// *************************************************************************
  
//  array->useArray(field_av,true,MEDCoupling::CPP_DEALLOC,n_elements_med,1);
  std::copy(field_av,field_av+n_elements_med*nComp,array->getPointer());

  _field->setArray(array);
  _field->setName(symbolic_eq);
  _field->checkConsistencyLight();
  std::cout << "InterfaceFunctionM::set_analytic_field_elem_interface \n";
//   printOn(std::cout,1);
  return;
}

// ==========================================================================
void InterfaceFunctionM::set_field(
  const MEDCoupling::MEDCouplingFieldDouble *f
) {// =======================================================================
  if(_field!=NULL) _field->decrRef();
  _field = f->deepCopy();
  _field->setName(f->getName());
  return;
}

// ==========================================================================
void InterfaceFunctionM::set_field_source(
  const MEDCoupling::MEDCouplingFieldDouble *f
) {// =======================================================================
  if(_field!=NULL) _field->decrRef();
  _field = f->deepCopy();
  _field->setName(f->getName());
  _field->checkConsistencyLight();

}

// ==========================================================================
MEDCoupling::MEDCouplingFieldDouble * InterfaceFunctionM::getField(
  char const* /*name */) {// =====================================================
//   if (FDEBUG) fDebugPos << "InterfaceFunctionM::getField _field = " << _field << std::endl;
//   if (!_field) return NULL;
//   if (FDEBUG) fDebugPos << "InterfaceFunctionM::getField _field = " << _field << std::endl;
//   _field->incrRef();
//   if (FDEBUG) fDebugPos << "InterfaceFunctionM::getField " << name << std::endl;
//   _field->setName(name);
  return _field;
}
// // ==========================================================================
// void InterfaceFunctionM::eval(
//   int id,
//   std::vector<double> & val
// ) {// =======================================================================
//   if(_field == NULL) {for(int i=0; i<(int)val.size(); i++) val[i]=.0;}
//   else {
//     for(int i=0; i<(int)val.size(); i++) val[i]=_field->getIJ(_elemID[id],i);
//   }
// }
//
// // ========================================================================
// void InterfaceFunctionM::eval_elem(
//   int iel,
//   std::vector<double> & val
// ) { // ====================================================================
//   if(_field == NULL) {for(int i=0; i<(int)val.size(); i++) val[i] =0.0;}
//   else {
// //     int mid = _faceID[std::pair<int, int>(id, side)];
//     for(int i = 0; i<(int)val.size(); i++) val[i] = _field->getIJ(iel, i);
//   }
// }

// ========================================================================
void InterfaceFunctionM::eval(
  int node_med,
  int n_cmp,
  double val[]
) { // ====================================================================
  if(_field == NULL) {for(int i=0; i<n_cmp; i++) val[i] =0.0;}
  else {
//     int mid = _faceID[std::pair<int, int>(id, side)];
    for(int i = 0; i<n_cmp; i++) val[i] = _field->getIJ(node_med, i);
  }
}

// // ========================================================================
// void InterfaceFunctionM::eval(
//   int id,
//   int side,
//   std::vector<double> & val
// ) { // ====================================================================
// //   if(_field == NULL) {for(int i=0; i<(int)val.size(); i++) val[i] =0.0;}
// //   else {
// //     int mid = _faceID[std::pair<int, int>(id, side)];
// //     for(int i = 0; i<(int)val.size(); i++) val[i] = _field->getIJ(mid, i);
// //   }
// }
// =============================================================
// This function prints the interface function
void InterfaceFunctionM::printOn(
  std::ostream & out,                     ///< ostream file
  int id                                  ///< interface name
) const { // ===================================================

  // set up title --------------------------------------------------------
  out << std::endl
      << "================================================== \n"
      << std::endl
      << " InterfaceFunctionM =" << id<<  std::endl << std::endl;
  // ---------------------------------------------------------------------
  // ---------------- Mapping ---------------------------------------------
  //  node map MGMesh <-> MEDCoupling ------------------------------------
  out << "node/element correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<int, int>::const_iterator itN = _nodeID.begin();
  for(int i_mg=0; i_mg < _n; i_mg++) {
    out << "  " << std::setw(4) << _map_mg[i_mg]
        << " <-> " << std::setw(4) << _map_med[i_mg] << std::endl;
  }
  out << std::endl;

//   // element map MGMesh <-> MEDCoupling ----------------------------------
//   out << "element correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<int, int>::const_iterator itE = _elemID.begin();
//   for(; itE != _elemID.end(); itE++)
//     out << "  " << std::setw(4) << itE->first
//         << " <-> " << std::setw(4) << itE->second << std::endl;
//   out << std::endl;
  // face map    MGMesh <-> MEDCoupling ----------------------------------
//   out << "face correspondance (MGMesh <-> MEDCoupling)" << std::endl;
//   std::map<std::pair<int,int>, int>::const_iterator itF = _faceID.begin();
//   for(; itF != _faceID.end(); itF++)
//     out << " (" << std::setw(4) << itF->first.first
//         << ", " << std::setw(4) << itF->first.second << ")"
//         << " <-> " << std::setw(4) << itF->second << std::endl;
//   out << std::endl;

  // ---------------------------------------------------------------------
  // ----------------  field ---------------------------------------------
  if(!_field) {
    out << "field= NULL" << std::endl;
    out << "==================================================" << std::endl;
    return;
  }
  MEDCoupling::TypeOfField type = _field->getTypeOfField();
  const MEDCoupling::DataArrayDouble * d;
  int n;
  const MEDCoupling::MEDCouplingUMesh * m
    = dynamic_cast<const MEDCoupling::MEDCouplingUMesh *>(_field->getMesh());

  out << "name : " << _field->getName() << std::endl;

    int ncf=_field->getNumberOfComponents();
    int nc;
     MEDCoupling::DataArrayDouble * v = _field->getArray();
  
  switch(type) {
  case MEDCoupling::ON_NODES:
    out << "type : ON_NODES" << std::endl;
    d = m->getCoords();
    d->incrRef();
    n = d->getNumberOfTuples();
    out << " " <<n << " nodes";
    
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;

    nc = d->getNumberOfComponents();
  for(int i=0; i<n; i++) {
    for(int j=0; j<nc; j++) out << (j == 0 ? "f( " : ", ")
                                  << std::setw(10) << std::fixed << d->getIJ(i,j);
    out << ") = ";
    for(int j=0; j<ncf; j++)  out << (j == 0 ? "( " : ", ")
                                    << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << ")  ";
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;
    break;
    
    
  case MEDCoupling::ON_CELLS:
    out << "type : ON_CELLS" << std::endl;
//     d = m->getBarycenterAndOwner();
//     n = d->getNumberOfTuples(); 
    nc =m->getNumberOfCells();
    out << " " << nc << " cells";
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;
  
  for(int i=0; i<nc; i++) {
     out << "f( "
                                  << std::setw(4) << std::fixed << i;
    out << ") = ";
    for(int j=0; j<ncf; j++)  out << (j == 0 ? "( " : ", ")
                                    << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << ")  ";
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;
    break;
  default :
    std::cout<< "Error type"; abort();
  }

 

//   v->decrRef();
}

#endif



