#include "InterfaceFunctionDD.h"
#include "DONDRA.h"
#include "mesh_c.h"

#include <iomanip>
#include <limits>
#include <boost/concept_check.hpp>

#ifdef HAVE_MED
// #include "MGFE_conf.h"
// #include "MeshExtended.h"




#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"







// using namespace MEDCoupling;
// using namespace libMesh;
// ======================================================
InterfaceFunctionDD::~InterfaceFunctionDD() {
  delete [] _map_mg; delete [] _map_med;
  if(_field) _field->decrRef(); if(_support_med) _support_med->decrRef();
//   delete [] _mesh_mg;
}

// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionDD::set_mesh_BFinterface_nodeID(
//   const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id,                      ///< inrface identity  (in)
  const int order_cmp                          ///< order pt (1 or 2) (in)
) { // ========================================================================

//   double toll_dist=1.e-12;
// //   _nodeID.clear();  // clear map node: FEMus-mesh -> MED-mesh
//
//   // Femus-mesh
// //   _mesh_mg = mesh;                                  // MG mesh
//   int nlevels_mg=_mesh_mg->_NoLevels-1;             // top MG level
//   int n_nodes_mg=_mesh_mg->_NoNodes[nlevels_mg];    // top MG n nodes
//   int dim_mg= _mesh_mg->_dim;                       // MG dimension
//
//   std::vector<double> xyz_mg;
//   std::vector<int> nodeid_mg;
//
//   // Med-mesh
//   _support_med = support;                 // MED-mesh
//   _n= _support_med->getNumberOfNodes();   // n MED-nodes
//   const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
//   int n_nodes_med = d->getNumberOfTuples(); // n MED-tuple-cords
//
// //   std::vector<std::pair<double, int > > dist_med;
// //   std::vector<std::pair<double, int > >::iterator itdist;
// //
// //
// //   // MED-distance *************************************************************
// //   dist_med.resize(n_nodes_med);
// //   for(int i_med=0; i_med<n_nodes_med; i_med++) { // i_med= MED-node
// //     double dd=0.; // dd=MED-point istance from (0.0)
// //     for(int idim=0; idim<dim_mg; idim++) dd  += d->getIJ(i_med, idim)*d->getIJ(i_med, idim);
// //     dist_med[i_med] = std::make_pair(dd,i_med);
// //   }
// //   std::sort(dist_med.begin(), dist_med.end());
// //   // **************************************************************************
//
//   // MG-node identity vector **************************************************
//   int n_bd_nodes=0; // counter of mg mesh nodes on the boundary interface (interface_id)
//   switch(order_cmp) {
//   case 2:
//
//     //build xyz_mg and nodeid_mg finding nodes on boundary interface (interface_id) from mg mesh
//     for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
//       if(fabs(mesh->_bc_id[i]) == interface_id) {
//         n_bd_nodes++;
//         for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
//         nodeid_mg.push_back(i);
//       }
//     }
//
//     break;
//   case 1:
//     std::cout<< " \n InterfaceFunctionDD::set_mesh_BFinterface_nodeID:"<<
//              " order not implemented \n";
// //     xyzToSetOfNodes_onBd_linM(dim_mg,n_nodes_mg,mesh->_xyz,mesh->_bc_id,s2,n2);
//     break;
//   default:
//     std::cout<< " \n InterfaceFunctionDD::set_mesh_BFinterface_nodeID:"<<
//              " order not implemented \n";
//     abort();
//   } // ************************************************************************
//
//
//   _map_med = new int [n_bd_nodes];
//   _map_mg =  new int [n_bd_nodes];
//   _n     =           n_bd_nodes;
//
//
//   // Building map _nodeID[Mg-node]-> MED-node *********************************
//   for(int i_mg=0; i_mg<n_bd_nodes; i_mg++) {// i_mg=MGnode
//
// //     /// MG node distance (dist_mg)
// //     double dist_mg=0.;
// //     for(int idim=0; idim<dim_mg; idim++)
// //       dist_mg  += xyz_mg[i_mg*dim_mg+idim]*xyz_mg[i_mg*dim_mg+idim];
// //     // Find the first node (a) in the interval search [a,b]
// //     itdist = dist_med.begin();
// //     while(fabs(dist_mg - itdist->first) >  toll_dist)  ++itdist;
// //     if(itdist == dist_med.end()) {
// //       std::cout << "non ci sono punti con toll_dist" << toll_dist << '\n';
// //       abort();
// //     }
//     // finding node (k_min)  with mimimum distance (min)
//     int k_min=-1; double min_dist=1.e+15;
//
//     for(int i_med=0; i_med<n_nodes_med; i_med++) { // i_med= MED-node
// //     for(; itdist<dist_med.end(); itdist++) { // node med-interface loop
//       double dd=0.;// dd=distance(MED-point,MG-point)
//       for(int idim=0; idim<dim_mg; idim++)
//         dd  += fabs(xyz_mg[i_mg*dim_mg+idim] - d->getIJ(i_med, idim));
//       if(dd < min_dist)  {min_dist=dd; k_min=i_med;}
//     }
//     _map_med[i_mg]=k_min;
//     _map_mg[i_mg] =nodeid_mg[i_mg];
// //     _nodeID[nodeid_mg[i_mg]] =k_min; // filling the map _nodeID
//   } // ************************************************************************
//
//   // clean
//   //  std::cout << "==================== mesh node id  =============================== "<<  interface_id <<"\n";
// //   std::map<int, int>::const_iterator itN = _nodeID.begin();
// //   for(; itN != _nodeID.end(); itN++) {
// //     mesh->_node_id[itN->first]=itN->second;
// //   }
// //  printOn(std::cout,interface_id);
//   xyz_mg.clear(); nodeid_mg.clear();
// //   dist_med.clear();

  return;
}



// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionDD::set_mesh_interface_nodeID(
//   const MeshExtended * mesh,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_id,                      ///< inrface identity  (in)
  const int order_cmp                          ///< order pt (1 or 2) (in)
) { // ========================================================================

//   double toll_dist=1.e-12;
// //   _nodeID.clear();  // clear map node: FEMus-mesh -> MED-mesh
//
//   // Femus-mesh
//   _mesh_mg = mesh;                                  // MG mesh
//   int nlevels_mg=_mesh_mg->_NoLevels-1;             // top MG level
//   int n_nodes_mg=_mesh_mg->_NoNodes[nlevels_mg];    // top MG n nodes
//   int dim_mg= _mesh_mg->_dim;                       // MG dimension
//
//   std::vector<double> xyz_mg;  std::vector<int> nodeid_mg;
//
//   // Med-mesh
//   _support_med = support;                 // MED-mesh
// //   _n= _support_med->getNumberOfNodes();   // n MED-nodes
//   const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
//   int n_nodes_med = d->getNumberOfTuples(); // n MED-tuple-cords
//
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
//
//   // MG-node identity vector **************************************************
//   int n_bd_nodes=0; // counter of mg mesh nodes on the boundary interface (interface_id)
//   switch(order_cmp) {
//   case 2:
//
//     //build xyz_mg and nodeid_mg finding nodes on boundary interface (interface_id) from mg mesh
//     for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
//       if(fabs(mesh->_bc_id[i]) == interface_id) {
//         n_bd_nodes++;
//         for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
//         nodeid_mg.push_back(i);
//       }
//     }
//     break;
//
//   case 1:
//
//     //build xyz_mg and nodeid_mg finding only linear nodes on boundary interface (interface_id) from mg mesh
//     for(int i=0; i<n_nodes_mg; i++) { // cycle on mg nodes
//       if(mesh->_bc_id[i] == -interface_id) {
//         n_bd_nodes++;
//         for(int idim=0; idim<dim_mg; idim++)  xyz_mg.push_back(mesh->_xyz[i+idim*n_nodes_mg]); // d->getIJ(i, idim);
//         nodeid_mg.push_back(i);
//       }
//     }
//     break;
//
//    default:
//     std::cout<< " \n InterfaceFunctionDD::set_mesh_BFinterface_nodeID:"<<
//              " order not implemented \n";
//     abort();
//   } // ************************************************************************
//
//
//   _map_med = new int [n_bd_nodes];
//   _map_mg =  new int [n_bd_nodes];
//   _n     =           n_bd_nodes;
//
//
//   // Building map _nodeID[Mg-node]-> MED-node *********************************
//   for(int i_mg=0; i_mg<n_bd_nodes; i_mg++) {// i_mg=MGnode
//
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
//     // finding node (k_min)  with mimimum distance (min)
//     int k_min=-1; double min_dist=1.e+15;
//     for(; itdist<dist_med.end(); itdist++) { // node med-interface loop
//       double dd=0.;// dd=distance(MED-point,MG-point)
//       for(int idim=0; idim<dim_mg; idim++)
//         dd  += fabs(xyz_mg[i_mg*dim_mg+idim] - d->getIJ(itdist->second, idim));
//       if(dd < min_dist)  {min_dist=dd; k_min=itdist->second;}
//       if(fabs(dist_mg - itdist->first) >  toll_dist) break;
//     }
//     _map_med[i_mg]=k_min;
//     _map_mg[i_mg] =nodeid_mg[i_mg];
// //     _nodeID[nodeid_mg[i_mg]] =k_min; // filling the map _nodeID
//   } // ************************************************************************
//
//   // clean
//
// //   std::map<int, int>::const_iterator itN = _nodeID.begin();
// //   for(; itN != _nodeID.end(); itN++) {
// //     mesh->_node_id[itN->first]=itN->second;
// //   }
// #ifdef PRINT_MED
//  std::cout << "====interface_id= "<<  interface_id << "with nodes" <<  n_bd_nodes << "\n";
//  printOn(std::cout,interface_id);
// #endif
//   xyz_mg.clear(); nodeid_mg.clear();
//   dist_med.clear();

  return;
}





// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void InterfaceFunctionDD::set_mesh_interface_elemID(
  const MeshC  &mesh_c,                ///< Femus-mesh        (in)
  const MEDCoupling::MEDCouplingUMesh * support,///< med-mesh          (in)
  const int interface_name                      ///< inrface identity  (in)
) { // ========================================================================



// Med-mesh
  _support_med = support;                  // MED-mesh
  _n= _support_med->getNumberOfCells();    // n MED-elements (_n now is element)
  const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples();//  MED-nodes
  assert(n_nodes_med == _support_med->getNumberOfNodes());  // check  MED-nodes
  int dim_med = d->getNumberOfComponents();//  MED dimension
  int n_nod_el=(dim_med==3)?8:4;           // for Cartesian cell only
  double *xyz_med = new double [ _n*dim_med]; // MED el centers

   // Computing the MED center (xyz_med) -------------------------------------
  {
    std::vector<int> nodes1; // element nodes
    double x_m[3];//tmp buffer

    for(int ielem=0; ielem<  _n; ielem++) {
      for(int idim=0; idim<dim_med; idim++) x_m[idim]=0.; // zeros
      _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
      for(int inode=0; inode< n_nod_el; inode++) {
        for(int im=0; im<dim_med; im++) x_m[im] += d->getIJ(nodes1[inode],im);
      } // end inode
      nodes1.clear(); // clear element node vector
      for(int im=0; im<dim_med; im++) xyz_med[dim_med*ielem+im]=x_m[im]/n_nod_el;
    }
  }
  // Computing the _map_med and _map_mg ---------------------------------------
  _map_med = new int [_n];  _map_mg =  new int [_n];
  int nex=mesh_c.get_nex();  int ney=mesh_c.get_ney();  int nez=mesh_c.get_nez();
  for(int i_m=0; i_m<_n; i_m++) {
    const double xx=xyz_med[dim_med*i_m]; const double yy=xyz_med[dim_med*i_m+1];
    const double zz=xyz_med[dim_med*i_m+2];
    int kjx=0;  while(xx> mesh_c._x[kjx]) kjx++;                       kjx--;
    int kjy=0;  while(yy> mesh_c._y[kjy*(nex+1)+kjx]) kjy++;            kjy--;
    int kjz=0;  while(zz> mesh_c._z[kjz*(ney+1)*(nex+1)+kjx*kjy]) kjz++; kjz--;
    _map_med[i_m] = i_m;                     // med index
    _map_mg[i_m]  = kjx+kjy*nex+kjz*nex*ney;    // mg index
//     std::cout << " kjx " << kjx << " kjy " << kjy << " kjz " << kjz << std::endl;
//         std::cout << " x " << xx << " y " << yy << " z " << zz << std::endl;
  }

  std::cout << "InterfaceFunctionDD::set_mesh_interface_elemID: "<<  interface_name <<"\n";
  std::cout << "Number of elements: "<<  _n << '\n';

  delete [] xyz_med;
//   d->decrRef();
  return;
}



// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void InterfaceFunctionDD::set_analytic_field(
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
             "InterfaceFunctionDD::set_analytic_field: NULL function"; abort();
  }
  _field = _support_med->fillFromAnalyticNamedCompo(type, nComp, vars, symbolic_eq);
  _field->setName(symbolic_eq);
  _field->checkConsistencyLight();
//   std::ostream out("prova.field");
  std::cout << "InterfaceFunctionDD::set_analytic_field_interface \n";
//   printOn(std::cout);
  return;
}

// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void InterfaceFunctionDD::set_analytic_field_elem(
  const char *symbolic_eq,     // symbolic function
  int nComp             // number of componenents
) { // ========================================================================
//   if(_field) _field->decrRef();
//   //ON CELLS does not work (due to baryc in quad 9 or hex 27)
//   MEDCoupling::TypeOfField type = MEDCoupling::ON_NODES;
//   int dim=_support_med->getSpaceDimension();
//   std::vector<std::string> vars(dim);
//   if(dim > 0)    vars[0] = "x";
//   if(dim > 1)    vars[1] = "y";
//   if(dim > 2)    vars[2] = "z";
//   if(symbolic_eq==std::string("")) {
//     std::cout<<
//              "InterfaceFunctionDD::set_analytic_field_elem: NULL function"; abort();
//   }
//
//
//
//   const MEDCoupling::MEDCouplingFieldDouble * field_nodes =
//     _support_med->fillFromAnalytic3(type, nComp, vars, symbolic_eq);
//
//   int n_elements_med= _support_med->getNumberOfCells(); // n MED-elements
// //   const MEDCoupling::DataArrayDouble * d =_support_med->getCoords();// Med-mesh coordinates
//   int n_nodes_med = field_nodes->getNumberOfTuples();             //  MED-nodes
// //   int dim = field_nodes->getNumberOfComponents();             //  MED dimension
//
//   int n_nod_el=NDOF_FEM;
//   double * sum=new double[nComp];
//
// _field= MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::NO_TIME);
//  _field->setMesh(_support_med);
//
//   MEDCoupling::DataArrayDouble *array;
//   array=MEDCoupling::DataArrayDouble::New();
//   array->alloc(n_elements_med,nComp);
//
//   double *field_av = new double [n_elements_med*nComp]; // MED el centers
// //   std::vector<std::pair<double, int > > dist_med;       // MED el distance
// //   std::vector<std::pair<double, int > >::iterator itdist;
//
// //   _elemID.clear();  // clear map elements: FEMus-mesh -> MED-mesh
//   std::vector<int> nodes1; // element nodes
//
//   // Computing the MED center (xyz_med) ***************************************
//   for(int ielem=0; ielem<  n_elements_med; ielem++) {
//     for(int icomp=0; icomp<nComp; icomp++) sum[icomp]=0;
//     _support_med->getNodeIdsOfCell(ielem, nodes1);  // element nodes
//     for(int icomp=0; icomp<nComp; icomp++){
//     for(int inode=0; inode< n_nod_el; inode++)  sum[icomp] += field_nodes->getIJ(nodes1[inode],icomp);
//     field_av[ielem*nComp+icomp]=sum[icomp]/n_nod_el;
//   }
//     nodes1.clear(); // clear element node vector
//   }// *************************************************************************
//
// //  array->useArray(field_av,true,MEDCoupling::CPP_DEALLOC,n_elements_med,1);
//   std::copy(field_av,field_av+n_elements_med*nComp,array->getPointer());
//
//   _field->setArray(array);
//   _field->setName(symbolic_eq);
//   _field->checkCoherency();
//   std::cout << "InterfaceFunctionDD::set_analytic_field_elem_interface \n";
// //   printOn(std::cout,1);
  return;
}

// ==========================================================================
void InterfaceFunctionDD::set_field(
  const MEDCoupling::MEDCouplingFieldDouble *f
) {// =======================================================================
  if(_field) _field->decrRef();
//
  MEDCoupling::TypeOfField type = f->getTypeOfField();
  _field = MEDCoupling::MEDCouplingFieldDouble::New(type);
  *_field=*f;

//    _field->setName(f->getName());
  _field->checkConsistencyLight();

  std::cout << "InterfaceFunctionDD::set_field \n";

  return;
}

// ==========================================================================
void InterfaceFunctionDD::set_field_source(
  const MEDCoupling::MEDCouplingFieldDouble *f
) {// =======================================================================
//   if(_field) _field->decrRef();
// //
//   MEDCoupling::TypeOfField type = f->getTypeOfField();
//   _field = MEDCoupling::MEDCouplingFieldDouble::New(type);
//   *_field=*f;
// //   _field->setMesh(_support_med);
// //   MEDCoupling::DataArrayDouble *array = MEDCoupling::DataArrayDouble::New();
// //   int nTuples;
// //   int nComp;
// //
// //   const MEDCoupling::DataArrayDouble * arrayf = f->getArray(), *coord1, *coord2;
// // //   if (FDEBUG) fDebug.array(arrayf, "arrayf");
// //
// //   nTuples = _support_med->getNumberOfCells();
// //   nComp = f->getNumberOfComponents();
// //   array->alloc(nTuples, nComp);
// //
// //   coord1 = _support_med->getBarycenterAndOwner();
// //   coord2 = f->getMesh()->getBarycenterAndOwner();
// //
// //   array->fillWithValue(std::numeric_limits<double>::max());
// //   std::map<int, int> corresp;
// //   commonPoints(corresp, coord1, coord2);
// //
// //   for (std::map<int, int>::iterator it = corresp.begin();
// //        it != corresp.end(); it++) {
// //     for (int j=0; j<nComp; j++) {
// // //       if (FDEBUG) fDebugPos << " " << it->first << " -> " << it->second
// // //              << " val = " << arrayf->getIJ(it->first, j) << std::endl;
// //       array->setIJ(it->second, j, arrayf->getIJ(it->first, j));
// //     }
// //   }
// //   _field->setArray(array);
//   _field->setName(f->getName());
// //   array->decrRef();
//   _field->checkCoherency();
// //   printOn(std::cout);
}

// ==========================================================================
MEDCoupling::MEDCouplingFieldDouble * InterfaceFunctionDD::getField(
  char const* /*name */) {// =====================================================
//   if (FDEBUG) fDebugPos << "InterfaceFunctionDD::getField _field = " << _field << std::endl;
//   if (!_field) return NULL;
//   if (FDEBUG) fDebugPos << "InterfaceFunctionDD::getField _field = " << _field << std::endl;
//   _field->incrRef();
//   if (FDEBUG) fDebugPos << "InterfaceFunctionDD::getField " << name << std::endl;
//   _field->setName(name);
  return _field;
}
// // ==========================================================================
// void InterfaceFunctionDD::eval(
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
// void InterfaceFunctionDD::eval_elem(
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
void InterfaceFunctionDD::eval(
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
// void InterfaceFunctionDD::eval(
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
void InterfaceFunctionDD::printOn(
  std::ostream & out,                     ///< ostream file
  int id                                  ///< interface name
) const { // ===================================================

  // set up title --------------------------------------------------------
  out << std::endl
      << "================================================== \n"
      << std::endl
      << " InterfaceFunctionDD =" << id<<  std::endl << std::endl;
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



