#include <iostream>
#include <cstdlib>
#include <sstream>
#include <assert.h>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
// #include "Equations_conf.h"
// #include "Equations_tab.h"
#include "MGTimeLoop.h"
// class include
#include "FEMUS.h"


#include "MeshExtended.h"
#include "EquationSystemsExtendedM.h"

#ifdef HAVE_MED

// MED includes
#include "InterfaceFunctionM.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"




// ========================================================================
/// This function gets the Group mesh from the med file (medfile_name)
/// and sets the id and the mesh to the interfaces vector (_interface_mesh_vect[i])
/// through the index of the interface-functions (from EquationSystemsExtendedM)
// ========================================================================
void FEMUS::init_interface (
  const int interface_name,         ///< Interface name
  std::vector<int> IDSvec,          ///< Vector containing group IDS (union)
  int order_cmp,                    ///< Order of the numerical field
  const std::string & medfile_name, ///< Name of the MED file containing the mesh   (in)
  bool on_nodes,                    ///< Flag to indicate if the field is ON_NODES or ON_CELLS
  const int index_medmesh           ///< med-mesh index     (in)
) { // =========================================================================

  //  interface   names (vG[j]) ----------------------------------------------------------
  std::cout << "Interface printed on file:RESU_MED/Interface_mesh.med";
  std::vector<std::string> vG ( IDSvec.size() );
  std::cout << " Creating an interface with id " << interface_name << " containing groups ";

  for ( int j = 0; j < IDSvec.size(); j++ ) {
      vG[j] = to_string ( IDSvec[j] );
      std::cout << IDSvec[j] << " ";
      }

  std::cout << std::endl;
  // med mesh file  info  -----------------------------------------------------------------
  std::string mesh_dir, localFile, filename;
  std::vector<std::string> meshNames, MeshNames, FieldNames;
  GetInfo ( medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );

  // interface mesh file (FieldContainingMap) + map MedToMg (MedToMgMapArray)
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp;




  int FinerLevel = _mg_utils->_geometry["nolevels"] - 1;

  if ( on_nodes ) {
//     FieldContainingMap = MEDCoupling::ReadField(
//                            MEDCoupling::ON_NODES,filename.c_str(),"Mesh_Lev_"+to_string(FinerLevel),
//                            0,"FinerLevelNodeIDS_Lev_"+to_string(FinerLevel), -1,-1
//                          );

      tmp = MEDCoupling::ReadFieldNode (
              filename.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ),
              0, "FinerLevelNodeIDS_Lev_" + to_string ( FinerLevel ), -1, -1
            );
      }
  else {
//     FieldContainingMap = MEDCoupling::ReadField(
//                            MEDCoupling::ON_CELLS, filename.c_str(), "Mesh_Lev_"+to_string(FinerLevel),
//                            0, "MG_cell_id_Lev_"+to_string(FinerLevel), -1,-1
//                          );

      tmp = MEDCoupling::ReadFieldCell (
              filename.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ),
              0, "MG_cell_id_Lev_" + to_string ( FinerLevel ), -1, -1
            );
      }

  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
  FieldContainingMap ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp ) );
  MEDCoupling::DataArrayDouble * MedToMgMapArray =  MedToMgMapArray = FieldContainingMap->getArray();
  // support (mesh interface)
  MEDCoupling::MEDCouplingUMesh * support;
  int id_level = -1;

  if ( IDSvec[0] < 10 ) {
      id_level = 0;
      }

  support = MEDCoupling::ReadUMeshFromGroups (
              localFile.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ), id_level, vG
            );
  // interface_name+order_cmp  -> interface <-support+MedToMgMapArray
  init_interface ( interface_name, order_cmp, support, MedToMgMapArray );

  // print
  if ( get_proc() == 0 )  MEDCoupling::WriteUMesh ( "RESU_MED/Interface_mesh.med", support, true );

  std::cout << "Interface " << interface_name << " on group " << vG[0] << "  on_nodes " << on_nodes;
  std::cout << " NUMBER OF NODES " << support->getNumberOfNodes() << std::endl;

//   FieldContainingMap->decrRef();

  return;
  }

// ========================================================================
/// This function gets the volume mesh from the med file (medfile_name)
/// and sets the id and the mesh to the interfaces vector (_interface_mesh_vect[i])
/// through the index of the interface-functions (from EquationSystemsExtendedM)
// ========================================================================
// ================================================================================================
void FEMUS::init_interface (
  const int interface_name,          ///< Interface name           (in)
  int order_cmp,                     ///< Order of the numerical field    (in)
  const std::string & medfile_name,  ///< Name of the MED file containing the mesh   (in)
  bool on_nodes                      ///< Flag to indicate if the field is ON_NODES or ON_CELLS  (in)
) { // ==================================================================

  // Reading mesh Names from med file  ------------------------------
  std::string mesh_dir, localFile, filename;
  std::vector<std::string> meshNames, MeshNames, FieldNames;
  GetInfo ( medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );

  // interface mesh file (FieldContainingMap) + map MedToMg (MedToMgMapArray)
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField>tmp;
  MEDCoupling::DataArrayDouble * MedToMgMapArray;

  int FinerLevel = _mg_utils->_geometry["nolevels"] - 1;

  if ( on_nodes ) {
      tmp = MEDCoupling::ReadFieldNode (
              /*MEDCoupling::ON_NODES,*/ filename.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ), 0,
              "FinerLevelNodeIDS_Lev_" + to_string ( FinerLevel ), -1, -1
            );
      }
  else {
      tmp = MEDCoupling::ReadFieldCell (
              /*MEDCoupling::ON_CELLS,*/ filename.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ), 0,
              "MG_cell_id_Lev_" + to_string ( FinerLevel ), -1, -1
            );
      }

  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
  FieldContainingMap ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp ) );
  MedToMgMapArray = FieldContainingMap->getArray();

  int id_level = 0;
  MEDCoupling::MEDCouplingUMesh * support = MEDCoupling::ReadUMeshFromFile (
        localFile.c_str(), "Mesh_Lev_" + to_string ( FinerLevel ), id_level
      );

  init_interface ( interface_name, order_cmp, support, MedToMgMapArray );

//   FieldContainingMap->decrRef();

  return;
  }
// ========================================================================
/// This function gets the volume mesh from the med file as default
// ========================================================================
void FEMUS::init_par_interface (
  int order_cmp,    ///< Order (piecewise, linear, quadratic)
  bool on_nodes     ///< Flag to indicate if the field is ON_NODES or ON_CELLS  (in)
) {
  // Reading mesh Names from med file  ------------------------------
  const int levels = _mg_mesh->_NoLevels;
  MEDCoupling::MEDCouplingFieldDouble * NodeMap = NULL;
  MEDCoupling::MEDCouplingUMesh * support = NULL;
  getNodeMapAndProcMeshAtLevel ( levels - 1, support, NodeMap );
  MEDCoupling::DataArrayDouble * MedToMgMapArray = NodeMap->getArray();

  _ParallelInterfaceId = 1235 + get_proc();
  init_interface ( _ParallelInterfaceId, order_cmp, support, MedToMgMapArray );

  return;
  }

// ========================================================================
/// This function gets the NodeMapAndProcMeshAtLeve
// ========================================================================
void FEMUS::getNodeMapAndProcMeshAtLevel (
  int level,
  MEDCoupling::MEDCouplingUMesh *& FemusPar,
  MEDCoupling::MEDCouplingFieldDouble *& NodeMap
) {
  const int proc = get_proc();
  std::string mesh_dir, localFile, filename;
  std::vector<std::string> meshNames, MeshNames, FieldNames;
  GetInfo ( _mg_utils->_interface_mesh, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );

  MEDCoupling::MEDCouplingUMesh * FemusPart  = MEDCoupling::ReadUMeshFromFile (
        filename.c_str(), "Mesh_Lev_" + to_string ( level ), 0
      ); // ORIGINAL FEMUS FROM MED FILE
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp = MEDCoupling::ReadFieldCell (
        /*MEDCoupling::ON_CELLS,*/filename.c_str(), "Mesh_Lev_" + to_string ( level ), 0,
        "Proc_Lev_" + to_string ( level ), -1, -1 );
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
  ProcField ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp ) );

  // PIECEWISE FIELD CONTAINING PROC IDS
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp1 = MEDCoupling::ReadFieldNode (
        /*MEDCoupling::ON_NODES,*/filename.c_str(), "Mesh_Lev_" + to_string ( level ), 0,
        "FinerLevelNodeIDS_Lev_" + to_string ( level ), -1, -1
      );
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
  FinerLevelIDS ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp1 ) );

  double * FinLevId = const_cast<double *> ( FinerLevelIDS->getArray()->getPointer() );
  std::vector<int> CellIdProc;

  for ( int i = 0; i < ProcField->getNumberOfTuples(); i++ ) {
      if ( ( int ) ProcField->getIJ ( i, 0 ) == proc )  CellIdProc.push_back ( i );
      }

  int * CellId = new int[CellIdProc.size()];                                                                     // ARRAY CONTAINING CELL IDS OF SINGLE PROC MESH

  for ( int i = 0; i < CellIdProc.size(); i++ )  CellId[i] = CellIdProc[i];

  FemusPar = FemusPart->buildPartOfMySelf ( CellId, CellId + CellIdProc.size() );  // EXTRACTION OF PROC MESH FROM GLOBAL MESH
  FemusPar->zipCoords();                                                                                         // It changes nodes numbering -> new numbering is from 0 to FemusPar->getNumberOfNodes
  FemusPar->setName ( "FemusProcMesh_Lev_" + to_string ( level ) );

  // CREATING A MAP FROM LOCAL NODE NUMBERING TO GLOBAL NODE NUMBERING - RELATIVE TO MED FILE
  double * NodesId = NULL;
  MEDCoupling::DataArrayDouble * NodeArray = MEDCoupling::DataArrayDouble::New();
  NodeArray->alloc ( FemusPar->getNumberOfNodes(), 1 );
  NodesId = const_cast<double *> ( NodeArray->getPointer() );

  std::vector<int> FeParCellConn, FeTotCellConn;       // VECTOR FOR PARALLEL AND GLOBAL MESH CONNECTIVITIES

  int FemusCells        = FemusPar->getNumberOfCells();
  int FemusNodesPerCell = FemusPar->getNumberOfNodesInCell ( 0 );

  for ( int i_cell = 0; i_cell < FemusCells; i_cell++ ) {
      FemusPar->getNodeIdsOfCell ( i_cell, FeParCellConn );
      FemusPart->getNodeIdsOfCell ( CellId[i_cell], FeTotCellConn );

      for ( int i_cnode = 0; i_cnode < FemusNodesPerCell; i_cnode++ ) {
          NodesId[FeParCellConn[i_cnode]] = FinLevId[FeTotCellConn[i_cnode]];
          }

      FeParCellConn.clear();
      FeTotCellConn.clear();
      }// END MAP ---------------------------------------------------------------------------------

  NodeMap = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES ) ;
  NodeMap->setMesh ( FemusPar );
  NodeMap->setArray ( NodeArray );
  NodeMap->setName ( "NodesID_Lev_" + to_string ( level ) );

  return ;
  }


// =======================================================================================
/// Basic interface construction with support (InterfaceMesh) and  MEDToMg map (MEDToMgMapArray)
void FEMUS::init_interface (
  const int interface_name,                          ///< Interface name
  int order_cmp,                                     ///< Order (piecewise, linear, quadratic)
  MEDCoupling::MEDCouplingUMesh * InterfaceMesh,     ///< MED interface mesh
  MEDCoupling::DataArrayDouble * MEDToMgMapArray,    ///< Array containing the map with MED/MG numeration
  int interface_id                                   ///< Default interface_id
) {
  std::cout << "-------------------------------" << std::endl;
  std::cout << MEDToMgMapArray->getIJ ( 0, 0 ) << std::endl;
  std::cout << "-------------------------------" << std::endl;

  double * MapArray = const_cast<double *> ( MEDToMgMapArray->getPointer() );
  std::map<int, int> MedToMg;

  int celle = InterfaceMesh->getNumberOfCells();
  int NodesPerCell = InterfaceMesh->getNumberOfNodesInCell ( 0 );

  std::vector<int> MgNodesIds, MedNodesIds;

  // Storing cell connectivities inside vector MgNodesIds -> numeration of mesh mg
  for ( int i = 0; i < celle; i++ )  InterfaceMesh->getNodeIdsOfCell ( i, MgNodesIds );

  // Renumbering nodes and cells -> new numeration in accordance with mesh total number of nodes and cells
  InterfaceMesh->zipCoords();

  // Storing cell connectivities inside vector MedNodesIds -> numeration of mesh med
  for ( int i = 0; i < celle; i++ ) {
      InterfaceMesh->getNodeIdsOfCell ( i, MedNodesIds );
      }

  // Creating the map
  for ( int i = 0; i < celle * NodesPerCell; i++ ) {
      int ipp = MgNodesIds[i];
      int iss = ( int ) MapArray[ipp];
      MedToMg[MedNodesIds[i]] = iss ;
      }

  const int NODI = InterfaceMesh->getNumberOfNodes();
  int * map_med = new int[NODI];
  int * map_mg  = new int[NODI];

  for ( int i = 0; i < NODI; i++ ) {
      map_med[i] = i;
      map_mg[i]  = MedToMg[i];
      }

  MgNodesIds.clear();
  MedNodesIds.clear();

  if ( interface_id == 0 ) printf ( "\n \033[1;31m Interface %d support set to mesh without volume group and order %d and nodes %d\033[0m\n", interface_name, order_cmp, celle );
  else std::cout << "\n \033[1;31m Interface " << interface_name
                   << " support set to mesh with interface " << interface_id
                   << " and order " << order_cmp
                   << " and nodes " << NODI
                   << "\033[0m\n";

  InterfaceFunctionM * fun = new InterfaceFunctionM;
  fun->set_maps ( map_med, map_mg, NODI );
  fun->set_order ( order_cmp );
  fun->set_mg_mesh ( _mg_mesh );
  fun->set_support_med ( InterfaceMesh );
  fun->set_NumberOfNodes ( NODI );

  _mg_equations_map->add_interface_fun ( interface_name, fun );
  delete [] map_med;
  delete [] map_mg;
  MapArray = NULL;
  MedToMg.clear();
  return;
  }

// =======================================================================================
void FEMUS::write_Boundary_value (
  int id_boundary_name,       ///< identity interface name (in)
  std::string mgsystem_name, ///< system name          (in)
  int n_cmp,             ///< from variable system (in)
  int first_cmp               ///< to variable system   (in)
) {
  _mg_equations_map->write_Boundary_value ( id_boundary_name, mgsystem_name, n_cmp, first_cmp );
  return;
  }



// =============================================================================
//This function reads and sets the med-mesh
// (also the libmesh calling the other setMesh function)
void FEMUS::setMedMesh() {
  std::string dataFile = _mg_utils->_mesh_dir + _mg_utils->_interface_mesh;
  int l = dataFile.size();
  int proc = get_proc();

  // CHECK THAT MGMESH AND EQUATION MAP ARE INITIALIZED
  MyAssert ( _MgMeshInitialized, "FEMUS::setMedMesh MGMESH not initialized" );
  MyAssert ( dataFile.substr ( l - 4 ) == ".med", "FEMUS::setMesh: " + dataFile + " does not exist!" );

  std::string localFile = dataFile;
  std::vector<std::string> meshNames = MEDCoupling::GetMeshNames ( localFile.c_str() );

  MyAssert ( meshNames.size() > 0.5, " FEMUS::setMesh : no meshes in the file'" );

  _med_mesh = MEDCoupling::ReadUMeshFromFile ( localFile.c_str(), meshNames[0].c_str(), 0 );

  if ( _med_mesh == NULL ) {
      std::cout <<  " FEMUS::setMesh : unable to read the med-mesh'";
      abort();
      }

  // BUILDING LOCAL PROC MESH FROM GLOBAL MED MESH ================================================================
  int FinerLevel = _mg_utils->_geometry["nolevels"] - 1;

  std::cout << " Creating proc mesh for level " << FinerLevel << std::endl;


  for ( int lev = 0; lev < FinerLevel + 1; lev++ ) {
      MEDCoupling::MEDCouplingUMesh * FemusPart       = MEDCoupling::ReadUMeshFromFile ( dataFile, "Mesh_Lev_" + to_string ( lev ), 0 ); // ORIGINAL FEMUS FROM MED FILE
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp = MEDCoupling::ReadFieldCell ( /*MEDCoupling::ON_CELLS,*/
            dataFile,  "Mesh_Lev_" + to_string ( lev ), 0, "Proc_Lev_" + to_string ( lev ), -1, -1 );        // PIECEWISE FIELD CONTAINING PROC IDS

      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
      ProcField ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp ) );

      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp2 =
        MEDCoupling::ReadFieldCell ( /*MEDCoupling::ON_CELLS,*/
          dataFile,
          "Mesh_Lev_" + to_string ( lev ), 0, "MG_cell_id_Lev_" + to_string ( lev ), -1, -1 );
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
      CellField ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp2 ) );



      double * MgCellId = const_cast<double *> ( CellField->getArray()->getPointer() );

      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingField> tmp1 = MEDCoupling::ReadFieldNode ( /*MEDCoupling::ON_NODES,*/
            dataFile,
            "Mesh_Lev_" + to_string ( lev ), 0, "FinerLevelNodeIDS_Lev_" + to_string ( lev ), -1, -1 );
      MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble>
      FinerLevelIDS ( MEDCoupling::DynamicCast<MEDCoupling::MEDCouplingField, MEDCoupling::MEDCouplingFieldDouble> ( tmp1 ) );


      double * FinLevId = const_cast<double *> ( FinerLevelIDS->getArray()->getPointer() );

      std::vector<int> CellIdProc;

      for ( int i = 0; i < ProcField->getNumberOfTuples(); i++ ) {
          if ( ( int ) ProcField->getIJ ( i, 0 ) == proc ) {
              CellIdProc.push_back ( i );
              }
          }

      int * CellId = new int[CellIdProc.size()];                                                                     // ARRAY CONTAINING CELL IDS OF SINGLE PROC MESH

      for ( int i = 0; i < CellIdProc.size(); i++ ) {
          CellId[i] = CellIdProc[i];
          }

      MEDCoupling::MEDCouplingUMesh * FemusPar = FemusPart->buildPartOfMySelf ( CellId, CellId + CellIdProc.size() ); // EXTRACTION OF PROC MESH FROM GLOBAL MESH
      FemusPar->zipCoords();                                                                                         // It changes nodes numbering -> new numbering is from 0 to FemusPar->getNumberOfNodes
      FemusPar->setName ( "FemusProcMesh_Lev_" + to_string ( lev ) );

      // CREATING A MAP FROM LOCAL NODE NUMBERING TO GLOBAL NODE NUMBERING - RELATIVE TO MED FILE
      double * NodesId = NULL;
      MEDCoupling::DataArrayDouble * NodeArray = MEDCoupling::DataArrayDouble::New();
      NodeArray->alloc ( FemusPar->getNumberOfNodes(), 1 );
      NodesId = const_cast<double *> ( NodeArray->getPointer() );

      std::vector<int> FeParCellConn, FeTotCellConn;       // VECTOR FOR PARALLEL AND GLOBAL MESH CONNECTIVITIES

      int FemusCells        = FemusPar->getNumberOfCells();
      int FemusNodesPerCell = FemusPar->getNumberOfNodesInCell ( 0 );

      for ( int i_cell = 0; i_cell < FemusCells; i_cell++ ) {
          FemusPar->getNodeIdsOfCell ( i_cell, FeParCellConn );
          FemusPart->getNodeIdsOfCell ( CellId[i_cell], FeTotCellConn );

          for ( int i_cnode = 0; i_cnode < FemusNodesPerCell; i_cnode++ ) {
              NodesId[FeParCellConn[i_cnode]] = FinLevId[FeTotCellConn[i_cnode]];
              }

          FeParCellConn.clear();
          FeTotCellConn.clear();
          }// END MAP ---------------------------------------------------------------------------------

      _NodeMap.push_back ( MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES ) );
      _NodeMap[lev]->setMesh ( FemusPar );
      _NodeMap[lev]->setArray ( NodeArray );
      _NodeMap[lev]->setName ( "NodesID_Lev_" + to_string ( lev ) );

      MEDCoupling::WriteUMesh ( "RESU_MED/FEMUS_PAR_" + to_string ( proc ) + ".med", FemusPar, false );
      MEDCoupling::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/FEMUS_PAR_" + to_string ( proc ) + ".med", _NodeMap[lev] );

      FemusPar->decrRef();
      FemusPart->decrRef();
      NodeArray->decrRef();
//         CellField->decrRef(); FinerLevelIDS->decrRef();   ProcField->decrRef();
      CellId  = NULL;
      MgCellId = NULL;
      CellIdProc.clear();
      }

  return;
  }
// =============================================================================
void FEMUS::setExtField (
  const std::string & systemName,
  MEDCoupling::MEDCouplingFieldDouble * bcField
) {
  _mg_equations_map->write_Boundary_value ( systemName, bcField );
  return;
  }
// =============================================================================
MEDCoupling::MEDCouplingFieldDouble * FEMUS::GetExtField (
  const std::string & systemName
) {
  MEDCoupling::MEDCouplingFieldDouble * ExtField = _mg_equations_map->GetField ( systemName );
  return ExtField;
  }

// ===================================================================
void FEMUS::setAnalyticSource (
//   const std::string & bcName,           // boundary name
  int interface_name,
  int n_cmp,
  const std::string & bcExpression      // boundary symbolic expr
) {
  _mg_equations_map->setBC ( interface_name, n_cmp, bcExpression.c_str() );
  return;
  }


// =============================================================================
void FEMUS::setFieldSource (
  int interface_name,
  int n_cmp,
  const MEDCoupling::MEDCouplingFieldDouble * srcField ) {

  if ( srcField == NULL ) {
      return;
      }

  _mg_equations_map->setBC ( interface_name, n_cmp, srcField );

  return;
  }



// ============================================================================
/// This function gets all the values on boundary with identity id !!!!!!!!!! old !!!!!!!!!!!!!!!!!!!!!!!!
MEDCoupling::MEDCouplingFieldDouble * FEMUS::getValuesOnInterface ( // field (out)
  int  interface_name,                     // boundary name (char*) (in)
  const std::string & systemName,                  // system name           (in)
  int n_cmp,                                        // component             (in)
  int first_cmp                                        // component             (in)
) {
  return _mg_equations_map->getValuesOnBoundary_nodes ( interface_name, systemName.c_str(), n_cmp, first_cmp );
  }

// ============================================================================
/// This function gets all the values on boundary with identity id
MEDCoupling::MEDCouplingFieldDouble * FEMUS::getValuesOnInterface_from_node ( // field (out)
  int  interface_name,            // boundary name (char*) (in)
  const std::string & systemName, // system name           (in)
  int n_cmp,                      // component             (in)
  int order,                      // order quad=2 lin=1    (in)
  int first_cmp                   // component             (in)
) {
  return _mg_equations_map->getValues_from_node ( interface_name, systemName.c_str(), n_cmp, order, first_cmp );
  }
// ============================================================================
/// This function gets all the values on boundary with identity id
MEDCoupling::MEDCouplingFieldDouble * FEMUS::getValuesOnInterface_from_cell ( // field (out)
  int  interface_name,             // boundary name (char*) (in)
  const std::string & systemName,  // system name           (in)
  int n_cmp,                       // component             (in)
  int order,                       // order quad=2 lin=1 (in)
  int first_cmp                                        // component             (in)
) {
  return _mg_equations_map->getValues_from_cell ( interface_name, systemName.c_str(), n_cmp, order, first_cmp );
  }









// ================================================================================================
MEDCoupling::MEDCouplingFieldDouble * FEMUS::getProcSolution ( // field (out)
  const std::string & systemName, ///< system name  (in)
  int n_cmp,                     ///<  component from  first_cmp  (in)
  int first_cmp,                 ///<  first component   (in)
  int Level                      ///<  level   (in)
) {
  MyAssert ( _NodeMap.size() > Level,
             "FEMUS::getProcSolution _NodeMap not initialized for level " + to_string ( Level )
             + " maximum levels " + to_string ( _NodeMap.size() ) );
  return _mg_equations_map->getProcValues (
           _NodeMap[Level], _GlobInterfaceId, systemName.c_str(), n_cmp, first_cmp, Level
         );
  }

// ============================================================================
/// This function gets the actual mesh of FEMUS problem
const MEDCoupling::MEDCouplingUMesh * FEMUS::getUMesh (
  int name
) {
  return _mg_equations_map->getUMeshCoupling ( name );
  }
// // // ===========================================================================
/// This function gets the original mesh of FEMUS problem
const MEDCoupling::MEDCouplingUMesh * FEMUS::getUMesh_orig (
  int name
) {
  return _mg_equations_map->getUMeshCoupling_orig ( name );
  }




// ==================================================================================
// MESH movement
// ==================================================================================
// ------------------------------------------------
/// This function moves the FEMUS interface according to a given displacement field
void FEMUS::update_interface (
  const int interface_name,
  int n_cmp,
  const   std::vector<MEDCoupling::MEDCouplingFieldDouble *> & srcField
) { // =========================================================================
  MEDCoupling::DataArrayDouble  * mg_disp = MEDCoupling::DataArrayDouble::Meld ( srcField[0]->getArray(),
      srcField[1]->getArray() );

  if ( n_cmp == 3 ) {
      mg_disp = MEDCoupling::DataArrayDouble::Meld ( mg_disp, srcField[2]->getArray() );
      }

  std::cout << "FEMUS::UPDATE: Tuples  src_field  " << mg_disp ->getNumberOfTuples() <<
            "  Comp    " << mg_disp ->getNumberOfComponents()
            << " and x value is " <<  mg_disp->getIJ ( 0, 0 ) << " and y value is " <<  mg_disp->getIJ ( 0, 1 )
            << std::endl;
  std::cout << "FEMUS::UPDATE: support  set to boundary with name " <<
            interface_name  << "\n";


  MEDCoupling::MEDCouplingUMesh * support;
  MEDCoupling::MEDCouplingUMesh * support_up;
//   int id_level=-1;  if(interface_id_1<10) {id_level=0;}
//   support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
  support = ( getUMesh_orig ( interface_name ) )->clone ( 1 );
  support_up = ( getUMesh ( interface_name ) )->clone ( 1 );
// support->zipCoords();

  MEDCoupling::DataArrayDouble * coord;
  MEDCoupling::DataArrayDouble * coord_up;

  coord = support->getCoords();
  coord_up = support_up->getCoords();
  int npt = coord->getNumberOfTuples();
  int ncomp = coord->getNumberOfComponents();

  MEDCoupling::DataArrayDouble * new_cord = MEDCoupling::DataArrayDouble::Add ( coord, mg_disp );

  std::cout << "Tuples    " << npt << "Comp    " << ncomp << std::endl;
  support->setCoords ( new_cord );

  InterfaceFunctionM * fct = _mg_equations_map->get_interface_fun ( interface_name );
  fct->update_support ( support );


  return;
  }


// ============================================================================
/// This function gets all the values on boundary with identity id
MEDCoupling::MEDCouplingFieldDouble * FEMUS::getDisplacement ( // field (out)
  int  interface_name,                     // boundary name (char*) (in)
  const std::string & systemName,                  // system name           (in)
  int n_cmp,                                        // component             (in)
  int first_cmp                                        // component             (in)
) {
  return _mg_equations_map->getDisplacement ( interface_name, systemName.c_str(), n_cmp, first_cmp );
  }



// =====================================================================
void FEMUS::GetInfo (
  std::string medfile_name,  ///<  medfile_name (in)
  std::string & mesh_dir,    ///< mesh_dir (out std::string )
  std::string & localFile,   ///< localFile  (out std::string )
  std::string & filename,    ///< filename   (out std::string  )
  std::vector < std::string > & meshNames, ///< meshNames  (out std::vector<std::string>)
  std::vector < std::string > & MeshNames, ///< MeshNames  (out std::vector<std::string>)
  std::vector < std::string > & FieldNames, ///< FieldNames  (out std::vector<std::string>)
  const int index_medmesh
) { // ---------------------------------------------------------------------------------------
  //  mesh_dir (out) +   localFile  (out)
  mesh_dir = _mg_utils->_mesh_dir;
  localFile = mesh_dir + medfile_name;
  //  mesh_Names (out)
  meshNames = MEDCoupling::GetMeshNames ( localFile.c_str() );

  if ( meshNames.size() < 1 ) std::cout <<  " FEMUS::setMesh : no meshes in the file'";

  //  filename(out)     FieldNames (out)    MeshNames  (out)
  filename = localFile;
  FieldNames = MEDCoupling::GetAllFieldNames ( filename.c_str() );
  MeshNames  = MEDCoupling::GetMeshNames ( filename.c_str() );

  return;
  };


void FEMUS::SetPieceFieldOnYdist (
  MEDCoupling::MEDCouplingFieldDouble * Field,
  MEDCoupling::MEDCouplingFieldDouble * CellMap
) {
  double * FieldVal = const_cast<double *> ( Field->getArray()->getPointer() );
  double * CellArray = const_cast<double *> ( CellMap->getArray()->getPointer() );

  int Cells = Field->getMesh()->getNumberOfCells();

  for ( int i = 0; i < Cells; i++ ) {
      int Cell = ( int ) CellArray[i];
      _mg_mesh->_VolFrac[Cell] = FieldVal[i];
      }

  if ( _mg_mesh->_NoLevels > 1 ) {

      int ChildCells = pow ( 2, DIMENSION );

      for ( int Lev = _mg_mesh->_NoLevels - 2; Lev >= 0; Lev-- ) {

          int ActualLevelOffset = _mg_mesh->_off_el[0][Lev];
          int FinerLevelOffset  = _mg_mesh->_off_el[0][Lev + 1];

          int ActCells = FinerLevelOffset - ActualLevelOffset;

          for ( int cell = 0; cell < ActCells; cell++ ) {
              double val = 0;

              for ( int child = 0; child < ChildCells; child++ )
                val +=  _mg_mesh->_VolFrac[FinerLevelOffset + cell * ChildCells + child];

              _mg_mesh->_VolFrac[ActualLevelOffset + cell] = val / ChildCells;

              }
          }
      }


  return;
  }



#endif



// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 

