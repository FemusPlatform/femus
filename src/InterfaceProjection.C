#include <set>
#include <cstring>
#include <ctime>
#include "InterfaceProjection.h"
// #include<map>

#include "Solverlib_conf.h"
#include "MGFE_conf.h"
#include "MGFE.h"

#ifdef HAVE_MED

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"

// ================================================================================================
InterfaceProjection::InterfaceProjection() : MedUtils() {
  __Filled = 0;
  _AlreadyInitialized = 0;
  }


// ================================================================================================
InterfaceProjection::InterfaceProjection (
  const MEDCoupling::MEDCouplingUMesh * SourceMesh,
  const MEDCoupling::MEDCouplingUMesh * TargetMesh,
  DomainType vol_sur,
  int DiscOrder
) : MedUtils() {
  __Filled = 0;
  _AlreadyInitialized = 0;
  _DiscOrder=DiscOrder;
  FillParameters ( SourceMesh, TargetMesh, vol_sur );
  }

InterfaceProjection::InterfaceProjection (
  const MEDCoupling::MEDCouplingUMesh * SourceMesh,
  const MEDCoupling::MEDCouplingUMesh * TargetMesh,
  int procId,
  DomainType   vol_sur,
  int DiscOrder
) : MedUtils() {

  __Filled = 0;
  _AlreadyInitialized = 0;
  _DiscOrder=DiscOrder;
  setProcId ( procId );
  FillParameters ( SourceMesh, TargetMesh, vol_sur );
  }

InterfaceProjection::InterfaceProjection (
  const MEDCoupling::MEDCouplingUMesh * SourceMesh,
  std::vector<double> coord,
  int procId,
  int DiscOrder
) : MedUtils() {
  __Filled = 0;
  _AlreadyInitialized = 0;
  _DiscOrder=DiscOrder;
  setProcId ( procId );
  FillParametersProbe ( SourceMesh, coord );
  }


// ================================================================================================
void InterfaceProjection::FillParameters (
  const MEDCoupling::MEDCouplingUMesh * SourceMesh,
  const MEDCoupling::MEDCouplingUMesh * TargetMesh,
  int DomainType,
  double XiEtaToll
) {

  // Important!
  // Update the map when interpolation between different element types is available  //
  BuildInterpCoordMap();
//   _DiscOrder = 1;

  if ( __Filled == 1 ) {
      __InMesh->decrRef();
      __OutMesh->decrRef();
      }

  __InMesh = SourceMesh->deepCopy();
  __OutMesh = TargetMesh->deepCopy();
  __Filled = 1;

  std::string name = SourceMesh->getName();
  _SrcCellNodes        = __InMesh->getNumberOfNodesInCell ( 0 );
  _SrcCells = __InMesh -> getNumberOfCells();

  INTERP_KERNEL::NormalizedCellType Type = __InMesh->getTypeOfCell ( 0 );
  _SrcCoordInterpNodes = _InterpCoordsMap[Type][(_DiscOrder+1)%2];
  if ( Type == 6 || Type == 7 || Type == 14 || Type == 20 ) _FamilyType = 0;
  else _FamilyType = 1;

  int TrgCellNodes = __OutMesh->getNumberOfNodesInCell ( 0 );
  _SpaceDim            = __OutMesh->getSpaceDimension();
  _MeshDim             = __OutMesh->getMeshDimension();
  _TrgNodes = __OutMesh-> getNumberOfNodes();               // number of nodes in target mesh
  _TrgCells = __OutMesh-> getNumberOfCells();               // number of cells in target mesh

// Number of nodes of source mesh per mesh element used for the calculation
  // of the target mesh point coordinates in the canonical element
  /////////////////////////////////////////////////////////////////////
  std::clock_t par_time1 = std::clock();

  MEDCoupling::DataArrayDouble * XiEta1 = MEDCoupling::DataArrayDouble::New();
  XiEta1->alloc ( _TrgNodes, _MeshDim );
  XiEta1->fillWithValue ( 0. );

  MEDCoupling::DataArrayInt * BoundingNodes1 = MEDCoupling::DataArrayInt::New();
  BoundingNodes1->alloc ( _TrgNodes, _SrcCellNodes );
  BoundingNodes1->fillWithValue ( -1 );

  MEDCoupling::DataArrayInt * targetArray = MEDCoupling::DataArrayInt::New();
  targetArray -> alloc ( _TrgNodes, 1 );
  targetArray -> fillWithValue ( 1 );

  MEDCoupling::DataArrayDouble * CellBelonging = MEDCoupling::DataArrayDouble::New();
  CellBelonging -> alloc ( _TrgNodes, 1 );
  CellBelonging -> fillWithValue ( -1 );

  //***************************************************************
  // definizioni
  std::vector<double> coord;                              // Vector containing cell node coordinates
  std::vector< int > TrgConn;                             // Vector containing the cell node ids
  double * PointsCoords = new double [TrgCellNodes * DIMENSION]; // Node coord*/inates {xyz xyz xyz}
//     int *nodalConnPerCell = new int[TrgCellNodes];
  double TCbbox[2 * DIMENSION];
  double minmax[2 * DIMENSION];

  if ( _AlreadyInitialized == 0 ) {
      InitFe();
      _AlreadyInitialized = 1;
      }

  const int Quad4[8] = {-1, 1, 1, -1, -1, -1, 1, 1  }; // (xi,eta,zeta)
  const int Edge2[2] = {-1, 1}; // (xi,eta,zeta)

  switch ( DomainType ) {
    case ( Boundary ) :

      // ********* Loop over the target  mesh cells ************** //
      for ( int iCell = 0; iCell < _TrgCells; iCell++ ) {

          __OutMesh -> getNodeIdsOfCell ( iCell, TrgConn );

          for ( int i_node = 0; i_node < TrgCellNodes; i_node++ ) {           // Loop over the element nodes
              __OutMesh -> getCoordinatesOfNode ( TrgConn[i_node], coord );

              for ( int j = 0; j < _SpaceDim; j++ ) {
                  PointsCoords[i_node * _SpaceDim + j] = coord[j];
                  }

              coord.clear();
              }

// ----- Creation of a bounding box for iCell target mesh cell ------

//       int nodalConnPerCell[TrgCellNodes];
//       for(int i = 0; i < TrgCellNodes; i++) { nodalConnPerCell[i]=i; }

          //*********************************************************************//
          // The bounding box can be created in two different ways:
          // -) from the mesh of the single cell: mesh -> getBoundingBox();
          // -) from the DataArrayDouble containing the node coordinates: coordsArr->getMinMaxPerComponent();
          // Some tests should be performed in order to find the fastest method

          /// Method using mesh
//       MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("My2DMesh",2);
//       mesh->allocateCells(1);
//       mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD9,9,nodalConnPerCell);
//       mesh->finishInsertingCells();
//       double TCbbox[2*SpaceDim];
//       DataArrayDouble *coordsArr=DataArrayDouble::New();
//       coordsArr->alloc(SrcCellNodes,SpaceDim);
//       std::copy(PointsCoords,PointsCoords + SrcCellNodes*SpaceDim,coordsArr->getPointer());
//       mesh->setCoords(coordsArr);
//       coordsArr->decrRef();
//       mesh->getBoundingBox(TCbbox);
//       mesh->decrRef();

          /// Method using coordinates array
          // Computing TCbbox
          MEDCoupling::DataArrayDouble * coordsArr = MEDCoupling::DataArrayDouble::New();
          coordsArr->alloc ( _SrcCellNodes, _SpaceDim );
          std::copy ( PointsCoords, PointsCoords + TrgCellNodes * _SpaceDim, coordsArr->getPointer() );
          coordsArr->getMinMaxPerComponent ( TCbbox );
          coordsArr->decrRef();

          MEDCoupling::DataArrayInt * vettore = MEDCoupling::DataArrayInt::New();
          vettore = __InMesh->getCellsInBoundingBox ( TCbbox, 1.e-5 );
          int SrcIntCells = vettore->getNbOfElems();                       // Number of source mesh cells intersecting target mesh cell bounding box
          MEDCoupling::MemArray<int> SrcIntCellsIds = vettore->accessToMemArray();      // Array containing source mesh cell ids


          double ** ExtrCoord = new double*[SrcIntCells];                   // Matrix containing Min and Max coordinates of source mesh cells intersecting target mesh cell bounding box

          for ( int i = 0; i < SrcIntCells; ++i ) {
              ExtrCoord[i] = new double [2 * _SpaceDim] ;
              }

          // cell 0: xmin, xmax, ymin, ymax, zmin, zmax
          /// Cycle over the source mesh cells intersecting the target mesh cell bounding box
          std::vector< int > SrcConn;

          for ( int IntCellNum = 0; IntCellNum < SrcIntCells; IntCellNum++ ) {

              __InMesh ->  getNodeIdsOfCell ( SrcIntCellsIds[IntCellNum], SrcConn );
              int NNodes = SrcConn.size();
              MEDCoupling::DataArrayDouble * XYZ = MEDCoupling::DataArrayDouble::New();
              XYZ->alloc ( NNodes, _SpaceDim );

              /// Min and Max coordinates of source mesh cells intersecting target mesh cell bounding box
              for ( int i_node = 0; i_node < NNodes; i_node++ ) {            // Loop over the element nodes
                  std::vector< double > ElCoord;
                  __InMesh -> getCoordinatesOfNode ( SrcConn[i_node], ElCoord );

                  for ( int dimension = 0; dimension < _SpaceDim; dimension ++ ) {
                      XYZ->setIJ ( i_node, dimension, ElCoord[dimension] );
                      }

                  ElCoord.clear();
                  }

              double minmax[2 * _SpaceDim];
              XYZ->getMinMaxPerComponent ( minmax );

              for ( int i = 0; i < 2 * _SpaceDim; i++ ) {
                  ExtrCoord[IntCellNum][i] = minmax[i];
                  }

              SrcConn.clear();
              XYZ->decrRef();
              }

          // Determination of source mesh cell containing the target mesh node
          int IdSrcCell;
          double NodeCoord[DIMENSION];     // Vector containing target mesh node coordinates

          // Loop over the target mesh cell nodes
          for ( int j = 0; j < TrgCellNodes; j++ ) {            
              if ( targetArray->getIJ ( TrgConn[j], 0 ) > 0.5 ) {
                  std::vector<int> CellIds;                             // Vector where we store the Ids of cells that can contain the target mesh node j

                  for ( int IntCellNum = 0; IntCellNum < SrcIntCells; IntCellNum ++ ) {
                      if ( _MeshDim == 2 ) {
                          /// verifica dell'appartenenza del nodo j alla cella valori[conf]
                          if ( //
                            PointsCoords[j * _SpaceDim]   >= ExtrCoord[IntCellNum][0] - 1.e-9 && PointsCoords[j * _SpaceDim]   <= ExtrCoord[IntCellNum][1] + 1.e-9 &&
                            PointsCoords[j * _SpaceDim + 1] >= ExtrCoord[IntCellNum][2] - 1.e-9 && PointsCoords[j * _SpaceDim + 1] <= ExtrCoord[IntCellNum][3] + 1.e-9 &&
                            PointsCoords[j * _SpaceDim + 2] >= ExtrCoord[IntCellNum][4] - 1.e-9 && PointsCoords[j * _SpaceDim + 2] <= ExtrCoord[IntCellNum][5] + 1.e-9
                          ) {
                              CellIds.push_back ( SrcIntCellsIds[IntCellNum] );
                              }
                          }

                      if ( _MeshDim == 1 ) {
                          /// verifica dell'appartenenza del nodo j alla cella valori[conf]
                          if ( //
                            PointsCoords[j * _SpaceDim]   >= ( ExtrCoord[IntCellNum][0] - 1.e-9 ) && PointsCoords[j * _SpaceDim]   <= ( ExtrCoord[IntCellNum][1] + 1.e-9 ) &&
                            PointsCoords[j * _SpaceDim + 1] >= ( ExtrCoord[IntCellNum][2] - 1.e-9 ) && PointsCoords[j * _SpaceDim + 1] <= ( ExtrCoord[IntCellNum][3] + 1.e-9 )
                          ) {
                              CellIds.push_back ( SrcIntCellsIds[IntCellNum] );
                              }
                          }
                      }

                  int numberOfCell = CellIds.size();

                  std::vector< int > SrcNodesConn;    // connectivity to be filled from getNodeIdsOfCell med function

                  for ( int direction = 0; direction < _SpaceDim; direction++ ) {
                      NodeCoord[direction] = PointsCoords[j * _SpaceDim + direction]; // real coordinates
                      }

                  std::vector<double> XiEtaBound;                               // Vector containing target mesh node reference coordinates

                  // If there are two possible cells containing the target mesh node we search the one having fabs(xi)<=1 and fabs(eta)<=1

                  if ( numberOfCell  > 1 ) { //  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      int FirstTry = 0;
                      bool CellFound = false;

                      while ( !CellFound && FirstTry < numberOfCell ) { // ++++++++ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                          IdSrcCell = CellIds[FirstTry];
                          __InMesh ->  getNodeIdsOfCell ( IdSrcCell, SrcNodesConn );

                          for ( int SNode1 = 0; SNode1 < _SrcCellNodes; SNode1 ++ ) {
                              std::vector< double > ElCoord1;
                              __InMesh -> getCoordinatesOfNode ( SrcNodesConn[SNode1], ElCoord1 );

                              for ( int direction = 0; direction < _SpaceDim; direction++ ) {
                                  _CoordMatrix.push_back ( ElCoord1[direction] );
                                  }

                              ElCoord1.clear();
                              }

                          if ( _MeshDim == 2 ) XiEtaCalc_2D ( NodeCoord, XiEtaBound, Quad4 );

                          if ( _MeshDim == 1 ) XiEtaCalc_1D ( NodeCoord, XiEtaBound );
                                               
//                           if( j == 6 && iCell == 2387){
//                                        std::cout<<"xi0 "<<XiEtaBound[0]<<" xi1 "<<XiEtaBound[1]<<" cell "<<IdSrcCell<<" possible "<<numberOfCell<<"  SrcCellNodes "<<_SrcCellNodes<<std::endl;
//                                        std::cout<<"pos0: "<<NodeCoord[0]<<"  pos1: "<<NodeCoord[1]<< " pos2: "<<NodeCoord[2] <<std::endl;
//                                        for(int i=0; i<_SrcCellNodes; i++){
//                                           std::cout<<"Node "<<i<<" Glob Num: "<<SrcNodesConn[i];
//                                           for(int dim=0; dim<_SpaceDim; dim++)
//                                             std::cout<<"  pos"<<dim<<": "<<_CoordMatrix[i*_SpaceDim + dim];
//                                           std::cout<<"\n";
//                                        }
//                             
//                                    }
                          
                          if ( _MeshDim == 2 ) { // ---------------------------------------------------------------------------------------------
                              if ( fabs ( XiEtaBound[0] ) <= 1. + 1.e-5 &&  fabs ( XiEtaBound[1] ) <= 1. + 1.e-5 ) {
                                  CellFound = true;
                                  CellBelonging->setIJ ( TrgConn[j], 0, IdSrcCell );
//                                   std::cout << " cell found for node *********** " << j << " of cell " << iCell
//                                             << " calculated xi: " << XiEtaBound[0] << " eta: " << XiEtaBound[1] << std::endl;

                                   int val = BoundingNodes1->getIJ ( TrgConn[j], 0 );    
                                   
                                  if(val == -1){
                                     for ( int node = 0; node < _SrcCellNodes; node++ ) {
                                         BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                                         }
                                     }
                                  }
                              else {
//                                   std::cout << " cell not found for node ************ " << j << " of cell " << iCell 
//                                             << " after " << FirstTry + 1 << " out of " << numberOfCell << " attempts "
//                                             << " xi: " << XiEtaBound[0] << " eta: " << XiEtaBound[1] << std::endl;
                                  CellFound = false;
                                  FirstTry ++;

                                  if ( FirstTry != numberOfCell ) {
                                      XiEtaBound.clear();    // clear the coords at the point
                                      }

                                  }
                              }  // ----------------------------------------------------------------------------------------------------------

                          if ( _MeshDim == 1 ) { // ---------------------------------------------------------------------------------------------
                              if ( fabs ( XiEtaBound[0] ) <= 1. + 1.e-5 ) {
                                  CellFound = true;

                                  for ( int node = 0; node < _SrcCellNodes; node++ ) {
                                      BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                                      }
                                  }
                              else {
                                  CellFound = false;
                                  FirstTry ++;
                                  XiEtaBound.clear();
                                  }
                              } // ---------------------------------------------------------------------------------------------------------------

                          SrcNodesConn.clear();   // clear connectivity
                          _CoordMatrix.clear();
                          } // end while ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                      if ( !CellFound ) std::cout << " cell not found for node ************ " << j 
                                                  << " of cell " << iCell << " after " << FirstTry + 1 << " out of " << numberOfCell 
                                                  << " with discretization " << _DiscOrder
                                                  << " attempts " << " xi: " << XiEtaBound[0] << std::endl;

//             else
                      } // end if numeroCelle > 1.5 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                  else {
                      IdSrcCell = CellIds[0];
                      __InMesh ->  getNodeIdsOfCell ( IdSrcCell, SrcNodesConn );
                      CellBelonging->setIJ ( TrgConn[j], 0, IdSrcCell );
                      int srcdim = SrcNodesConn.size();
                      std::map <int, int> Ordering;

                      for ( int nMed = 0; nMed < SrcNodesConn.size(); nMed++ )   
                        if ( _MeshDim == 2 ) {
                            Ordering[nMed] = nMed;
                            }

                      if ( _MeshDim == 1 ) {
                          Ordering[0] = 0;
                          Ordering[1] = 1;
                          Ordering[2] = 2;
                          }

                      _CoordMatrix.clear();

                      for ( int SNode1 = 0; SNode1 < _SrcCoordInterpNodes; SNode1 ++ ) {
                          std::vector< double > ElCoord1;
                          __InMesh -> getCoordinatesOfNode ( SrcNodesConn[Ordering[SNode1]], ElCoord1 );

                          for ( int direction = 0; direction < _SpaceDim; direction++ ) {
                              _CoordMatrix.push_back ( ElCoord1[direction] );
                              }

                          ElCoord1.clear();
                          }

                      if ( _MeshDim == 2 ) {
                          XiEtaCalc_2D ( NodeCoord, XiEtaBound, Quad4 );
                          }

                      if ( _MeshDim == 1 ) {
                          XiEtaCalc_1D ( NodeCoord, XiEtaBound );
                          }
                          
                      int val = BoundingNodes1->getIJ ( TrgConn[j], 0 );    

                            if(val == -1){
                              for ( int node = 0; node < _SrcCellNodes; node++ ) {
                                BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                              }
                            }

                      SrcNodesConn.clear();
                      _CoordMatrix.clear();
                      }


                  for ( int dim = 0; dim < _MeshDim; dim++ ) {
                      XiEta1->setIJ ( TrgConn[j], dim, XiEtaBound[dim] );
                      }

                  targetArray->setIJ ( TrgConn[j], 0, 0 );
                  XiEtaBound.clear();
                  CellIds.clear();
                  }

              }

          TrgConn.clear();

          for ( int i = 0; i < SrcIntCells; ++i )   delete [] ExtrCoord[i];
          delete [] ExtrCoord;
          } // end loop

      break;

    // end case Boundary
    case ( Volume ) :


      std::clock_t par_time12 = std::clock();

      //****** new med search

      std::vector<double> CoordinatesOfNodes;
      std::vector<double> NodeCoord2;

      for ( int TNode = 0; TNode < _TrgNodes; TNode++ ) {
          __OutMesh->getCoordinatesOfNode ( TNode, CoordinatesOfNodes );
          }

//       elts ->     vector returning ids of found cells.
//       eltsIndex-> an array, of length nbOfPoints + 1,
//                   dividing cell ids in elts into groups each referring to one point.
//                   Its every element (except the last one) is an index pointing to the first id of
//                   a group of cells. Cells in contact with the i-th point are described
//                   by following range of indices: [ eltsIndex[ i ], eltsIndex[ i+1 ] )
//                   and the cell ids are elts[ eltsIndex[ i ]], elts[ eltsIndex[ i ] + 1 ].
      MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> elts;
      MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> eltsIndex;
      int dimm = CoordinatesOfNodes.size();
      double Nodi[dimm];

      for ( int dim = 0; dim < dimm; dim++ ) {
          Nodi[dim] = CoordinatesOfNodes[dim];    //
          }

      std::map<int, int> MedLibmesh;
      BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );

//******************************************************************************
// controllare metodo di ricerca ed utilizzo degli array elts e eltsIndex -> guarda classe modificata da daniele
      __InMesh->getCellsContainingPoints ( Nodi, _TrgNodes, 1.e-5, elts, eltsIndex );

      std::vector<double> NodeCoord;
      double * pos = new double[_SpaceDim];
      std::vector<int> SourceConn;
      std::vector< double > ElCoord1;

      int * elts_array = elts->getPointer();
      int * elts_indx_array = eltsIndex->getPointer();

      int count = 0;

      for ( int TNode = 0; TNode < _TrgNodes; TNode++ ) {
          __OutMesh->getCoordinatesOfNode ( TNode, NodeCoord );

          for ( int dir = 0; dir < _SpaceDim; dir++ ) {
              pos[dir] = NodeCoord[dir];
              }

          int NumPossibleCells = elts_indx_array[TNode + 1] - elts_indx_array[TNode];

          bool found = false;

          if ( NumPossibleCells > 0 ) {
//               std::cerr << "\n\n";
              int iCount = 0;

              while ( !found && iCount < NumPossibleCells ) {
                  _XiEtaChi.clear();
                  // cell ids are elts[ eltsIndex[ i ]],..., elts[ eltsIndex[ i ] + NumPossibleCells ].
                  int NewCell = elts_array[elts_indx_array[TNode] + iCount];

                  __InMesh->getNodeIdsOfCell ( NewCell, SourceConn );

                  for ( int SNode = 0; SNode < _SrcCellNodes; SNode++ ) {
                      BoundingNodes1->setIJ ( TNode, SNode, SourceConn[SNode] ); //SourceConn[MedLibmesh[SNode]]=med index

                      if ( SNode < _SrcCoordInterpNodes ) {
                          __InMesh -> getCoordinatesOfNode ( SourceConn[SNode], _CoordMatrix );   // only linear nodes
                          }
                      }

                  _Contained = 0;
//                   std::cerr << "\nTnode " << TNode << " Count " << iCount << " of " << NumPossibleCells << std::endl;
                  IsNodeInsideCell ( pos );

//                   if ( _Contained == 1 ) std::cerr << "Node " << TNode << " _Contained " << _Contained << "  ";

                  XiEtaChiCalc ( pos );

                  int Contained = IsNodeContained();

                  if ( Contained == 0 ) {
                      found = false;
                      _XiEtaChi.clear();
                      _CoordMatrix.clear();
                      SourceConn.clear();
                      }
                  else if ( Contained == 1 ) {
                      CellBelonging->setIJ ( TNode, 0, NewCell );
                      _CoordMatrix.clear();
                      SourceConn.clear();
                      found = true;
                      }

                  iCount++;
                  }// end while

              if ( !found ) {
                  for ( int SNode = 0; SNode < _SrcCellNodes; SNode++ ) BoundingNodes1->setIJ ( TNode, SNode, -1 ); //SourceConn[MedLibmesh[SNode]]=med index
                  }
              else {
                  for ( int Mcomp = 0; Mcomp < _MeshDim; Mcomp++ ) {
                      XiEta1->setIJ ( TNode, Mcomp, _XiEtaChi[Mcomp] );
                      }
                  }

              _CoordMatrix.clear();
              _XiEtaChi.clear();
              }  // end  if cell==-1
          else {
              count ++;
//         std::cout<< " Cell not found for target point "<< TNode << " count "<<count << " x: "<<pos[0]<<" y: "<<pos[1]<<" z: "<<pos[2]<<std::endl;
              }

          NodeCoord.clear();
          SourceConn.clear();
          }

      CoordinatesOfNodes.clear();
      delete [] pos;

      break;
      }

  delete [] PointsCoords;

  targetArray->decrRef();

  MEDCoupling::MEDCouplingFieldDouble * CanonicalPosition = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  CanonicalPosition->setMesh ( __OutMesh );
  CanonicalPosition->setArray ( XiEta1 );
  CanonicalPosition->setName ( "XiEta" );

  MEDCoupling::MEDCouplingFieldDouble * CellS = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  CellS->setMesh ( __OutMesh );
  CellS->setArray ( CellBelonging );
  CellS->setName ( "CellID" );

  std::string FileName = "S" + __InMesh->getName() + "_T" + __OutMesh->getName() + "_MeshCoupling.med";

  if ( _proc == 0 ) {
      MEDCoupling::WriteUMesh ( "RESU_MED/" + FileName, __OutMesh, false );
      MEDCoupling::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/" + FileName, CanonicalPosition );
      MEDCoupling::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/" + FileName, CellS );
      }

  CanonicalPosition->decrRef();
  CellS->decrRef();

  if ( _BoundingNodes != NULL ) _BoundingNodes->decrRef();

  _BoundingNodes = BoundingNodes1;

  if ( _XiEta != NULL ) _XiEta->decrRef();

  _XiEta = XiEta1;

  __BoundNodesPerCell = _SrcCellNodes;
  __TrgNodes = _TrgNodes;
  __Domain = DomainType;

  CellBelonging->decrRef();
  return;
  }

// ================================================================================================
void InterfaceProjection::FillParametersProbe (
  const MEDCoupling::MEDCouplingUMesh * SourceMesh,
  std::vector<double> coord,
  double XiEtaToll
) {
  // Important!
  // Update the map when interpolation between different element types is available  //
  std::map<int, int> InterpCoordNodes;
  InterpCoordNodes[27] = 8;
  InterpCoordNodes[9] = 4;
  InterpCoordNodes[3] = 2;
  InterpCoordNodes[10] = 4;
  InterpCoordNodes[6] = 3;
  InterpCoordNodes[7] = 3;

  if ( __Filled == 1 ) {
      __InMesh->decrRef();
      }

  __InMesh = SourceMesh->deepCopy();
  __Filled = 1;

  std::string name = SourceMesh->getName();
  _SrcCellNodes        = __InMesh->getNumberOfNodesInCell ( 0 );
  _SrcCoordInterpNodes = InterpCoordNodes[_SrcCellNodes];
  _SrcCells = __InMesh -> getNumberOfCells();



  INTERP_KERNEL::NormalizedCellType Type = __InMesh->getTypeOfCell ( 0 );

  if ( Type == 6 || Type == 7 || Type == 14 || Type == 20 ) _FamilyType = 0;
  else _FamilyType = 1;

  _SpaceDim            = __InMesh->getSpaceDimension();
  _MeshDim             = __InMesh->getMeshDimension();

// Number of nodes of source mesh per mesh element used for the calculation
// of the target mesh point coordinates in the canonical element

  std::clock_t par_time1 = std::clock();

  MEDCoupling::DataArrayDouble * XiEta1 = MEDCoupling::DataArrayDouble::New();
  XiEta1->alloc ( 1, 0 );
  XiEta1->fillWithValue ( 0. );

  MEDCoupling::DataArrayInt * BoundingNodes1 = MEDCoupling::DataArrayInt::New();
  BoundingNodes1->alloc ( 1, _SrcCellNodes );
  BoundingNodes1->fillWithValue ( -1 );

  double Bbox[2 * _SpaceDim];
  SourceMesh->getBoundingBox ( Bbox );
  double height = fabs ( Bbox[3] - Bbox[2] );
  double width =  fabs ( Bbox[1] - Bbox[0] );

  double pointA[2];
  pointA[0] = Bbox[0] - 1.e4 * width;
  pointA[1] = Bbox[2] - 1.e4 * height;
  double pointB[2];
  pointB[0] = Bbox[0] - 1.e4 * width;
  pointB[1] = Bbox[3] + 1.e4 * height;

  int SecMatrix[_SrcCells][4];

  MEDCoupling::DataArrayInt * targetArray = MEDCoupling::DataArrayInt::New();
  targetArray -> alloc ( 1, 1 );
  targetArray -> fillWithValue ( 1 );

  MEDCoupling::DataArrayDouble * CellBelonging = MEDCoupling::DataArrayDouble::New();
  CellBelonging -> alloc ( 1, 1 );
  CellBelonging -> fillWithValue ( -1 );

  //***************************************************************
  double TCbbox[2 * DIMENSION];
  double minmax[2 * DIMENSION];

  if ( _AlreadyInitialized == 0 ) {
      InitFe();
      _AlreadyInitialized = 1;
      }

  const int Quad4[8] = {-1, 1, 1, -1, -1, -1, 1, 1  }; // (xi,eta,zeta)
  const int Edge2[2] = {-1, 1}; // (xi,eta,zeta)

  std::clock_t par_time12 = std::clock();

  //****** new med search

//       elts ->     vector returning ids of found cells.
  MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> elts;
  MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> eltsIndex;
  int dimm = coord.size();
  double pos[dimm];

  for ( int dim = 0; dim < dimm; dim++ ) {
      pos[dim] = coord[dim];
      }

  std::map<int, int> MedLibmesh;
  BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );

//******************************************************************************
//         __InMesh->getCellsContainingPoint ( pos, 1.e-5, elts );
  __InMesh->getCellsContainingPoints ( pos, 1, 1.e-5, elts, eltsIndex );

//         Calculation of number of cells that can contain the point
  int NumPossibleCells = elts->getNbOfElems();

  std::vector<int> SourceConn;
  bool found = false;

  if ( NumPossibleCells > 0 ) {
      _OutOfDomain = 0;
      int iCount = 0;

      while ( !found && iCount < NumPossibleCells ) {
          _XiEtaChi.clear();
          int NewCell = elts->getIJ ( iCount, 0 );
          std::cout << "CHECKING CELL No " << NewCell << "\n";
          __InMesh->getNodeIdsOfCell ( NewCell, SourceConn );

          for ( int SNode = 0; SNode < _SrcCellNodes; SNode++ ) {
              BoundingNodes1->setIJ ( 0, SNode, SourceConn[SNode] );
              CellBelonging->setIJ ( 0, 0, NewCell );

              if ( SNode < _SrcCoordInterpNodes ) {
                  __InMesh -> getCoordinatesOfNode ( SourceConn[SNode], _CoordMatrix );   // only linear nodes
                  }
              }

          XiEtaChiCalc ( pos );

//                     int Contained = IsNodeContained();
          if ( _Contained == 0 ) {
              found = false;
              _XiEtaChi.clear();
              _CoordMatrix.clear();
              SourceConn.clear();
              }
          else if ( _Contained == 1 ) found = true;

          iCount++;
          } // end while

      if ( !found ) {
          for ( int SNode = 0; SNode < _SrcCellNodes; SNode++ ) BoundingNodes1->setIJ ( 0, SNode, -1 );
          }
      else {
          for ( int Mcomp = 0; Mcomp < _MeshDim; Mcomp++ ) {
              XiEta1->setIJ ( 0, Mcomp, _XiEtaChi[Mcomp] );
              }
          }

      _CoordMatrix.clear();
      _XiEtaChi.clear();
      }
  else {
      _OutOfDomain = 1;
      }

  SourceConn.clear();

  targetArray->decrRef();

  if ( _BoundingNodes != NULL ) _BoundingNodes->decrRef();

  _BoundingNodes = BoundingNodes1;

  if ( _XiEta != NULL ) _XiEta->decrRef();

  _XiEta = XiEta1;

  __BoundNodesPerCell = _SrcCellNodes;

  CellBelonging->decrRef();
  return;
  }


// ================================================================================================
InterfaceProjection::~InterfaceProjection() {
  __InMesh->decrRef();
  __OutMesh->decrRef();
  _XiEta->decrRef();
  _BoundingNodes->decrRef();
  }


// ================================================================================================
void InterfaceProjection::CheckBelonging ( bool & found ) {
  const double toll = 1.e-20;
  const double xmedium = ( _CoordMatrix[0] + _CoordMatrix[2] + _CoordMatrix[4] + _CoordMatrix[6] ) / 4.;
  const double ymedium = ( _CoordMatrix[1] + _CoordMatrix[3] + _CoordMatrix[5] + _CoordMatrix[7] ) / 4.;
  found = true;
  double S;
  double tildeS;
  double xm, ym, tildexm, tildeym, midX, midY ;

  for ( int dim = 0; dim < _SrcCoordInterpNodes; dim++ ) {
      int A = dim;
      int B = ( 4 - ( dim + 1 ) ) % 4;

      xm = 0.5 * ( xmedium - _pos[0] );
      xm = ( fabs ( xm < 1.e-9 ) ? 1.e-9 : xm );
      ym = 0.5 * ( ymedium - _pos[1] );
      ym = ( fabs ( ym < 1.e-9 ) ? 1.e-9 : ym );
      tildexm = 0.5 * ( _CoordMatrix[B * _SpaceDim] - _CoordMatrix[A * _SpaceDim] );
      tildexm = ( fabs ( tildexm < 1.e-9 ) ? 1.e-9 : tildexm );
      tildeym = 0.5 * ( _CoordMatrix[B * _SpaceDim + 1] - _CoordMatrix[A * _SpaceDim + 1] );
      tildeym = ( fabs ( tildeym < 1.e-9 ) ? 1.e-9 : tildeym );
      midX = 0.5 * ( _pos[0] + xmedium - ( _CoordMatrix[B * _SpaceDim] + _CoordMatrix[A * _SpaceDim] ) );
      midY = 0.5 * ( _pos[1] + ymedium - ( _CoordMatrix[B * _SpaceDim + 1] + _CoordMatrix[A * _SpaceDim + 1] ) );
      S = - ( midX - midY * tildexm / tildeym ) / ( xm - tildexm / tildeym );
      tildeS = S * ( midY + ym ) / tildeym;

      if ( fabs ( fabs ( S ) - 1. ) < toll || fabs ( fabs ( tildeS ) - 1. ) < toll ) {
          found = true;  //coincident nodes
          break;
          }

      if ( ( fabs ( S ) + toll ) < 1. && ( fabs ( tildeS ) + toll ) < 1. ) {
          found = false;    //cross
          break;
          }
      }

  return;
  }



// ===================================================================================
void InterfaceProjection::GetXiEta ( int node, std::vector<double> & XiEtabound ) {
  for ( int dim = 0; dim < _MeshDim; dim++ ) {
      XiEtabound.push_back ( _XiEta->getIJ ( node, dim ) );
      }

  return;
  }
// ===================================================================================
void InterfaceProjection::GetXiEtaChi ( int node, std::vector<double> & XiEtaChi ) {
  for ( int dim = 0; dim < _MeshDim; dim++ ) {
      XiEtaChi.push_back ( _XiEta->getIJ ( node, dim ) );
      }

  return;
  }


// ===================================================================================
void InterfaceProjection::GetInterNodes ( int node, std::vector<int> & VectorContainingNodes ) {
  const int SourceNodes = _BoundingNodes->getNumberOfComponents();

  for ( int i = 0; i < SourceNodes; i++ ) {
      VectorContainingNodes.push_back ( _BoundingNodes->getIJ ( node, i ) );
      }

  return;
  }


// ===================================================================================
int InterfaceProjection::GetBoundNodesPerCell() {
  return __BoundNodesPerCell;
  }


// ===================================================================================
int InterfaceProjection::GetDomainType() {
  return __Domain;
  }

// ===============================================================================================
void InterfaceProjection::XiEtaCalc_2D (
  double NodePos[],                // physical coord (in)
  std::vector<double> & XiEtaBound, // reference coord (out)
  const int  Quad4[],                   // point coordinates for quad4 (x0,y0,x1,y1...)
  int  npt_el                       // number of points for elment

) {// =============================================================================================

  double SquaredError = 1.e10;
  double test_err = SquaredError;
  int test_id = _SrcCellNodes - 1;
  int test = 0;
  int FoundXiEtaChi=0;
  const double Toll = 1.e-10;
  
  XiEtaBound.clear();
  double *XiEtaChi = new double[_MeshDim];
  
  while ( !FoundXiEtaChi && test < _SrcCellNodes ) {
//       if( _MeshDim==3 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordHex27[test + dim * 27];
      if( _MeshDim==2 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordQuad9[test + dim * 9];
//       if( _MeshDim==1 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordEdge3[test + dim * 3];
      
      SquaredError = CalcF_bound ( NodePos, XiEtaChi );

      if ( SquaredError < test_err ) {
          test_id = test;
          test_err = SquaredError;
          }

      if ( SquaredError < Toll ) FoundXiEtaChi = 1;

      test ++;
   }
      
   if(FoundXiEtaChi){
      for ( int dim = 0; dim < _MeshDim; dim++ ) 
          XiEtaBound.push_back ( XiEtaChi[dim] );
   }
   else{
  
        // geometry
        double xm[DIMENSION];
        double normal[DIMENSION];
        double Co[DIMENSION * DIMENSION];
        double Coord[NDOF_FEMB * DIMENSION];
        
        for ( int n_dim = 0; n_dim < DIMENSION; n_dim++ ) {
            xm[n_dim] = 0.;
        
            for ( int nnode = 0; nnode < npt_el; nnode++ ) {
                Coord[ nnode + n_dim * NDOF_FEMB] =  _CoordMatrix[n_dim + nnode * DIMENSION];
                xm[n_dim] += Coord[ nnode + n_dim * NDOF_FEMB] / npt_el;
                }
            }
        
        _fe[1]->normal_g ( Coord, xm, normal ); // normal in the new med surface
        
        // equation surface  coefficients
        //  0=XYZm[0] + Co[1]*psi+Co[2]*eta+Co[0]*psi*eta
        //  0=XYZm[1] + Co[4]*psi+Co[5]*eta+Co[3]*psi*eta
        //  0=XYZm[2] + Co[8]*psi+Co[7]*eta+Co[6]*psi*eta
        double XYZm[DIMENSION];
        double XYZs[DIMENSION];
        double *xieta = new double[_MeshDim]; // reference point coordinates
        
        for ( int kdim = 0; kdim < _SpaceDim; kdim++ ) {
            XYZm[kdim] = 0.;
            Co[kdim * _SpaceDim   ] =  1.e-20;
            Co[kdim * _SpaceDim + 1] =  1.e-20;
            Co[kdim * _SpaceDim + 2] =  1.e-20;
        
            for ( int node = 0; node < npt_el; node++ ) {
                Co[kdim * _SpaceDim   ] += 0.25 * _CoordMatrix[_SpaceDim * node + kdim] * Quad4[node] * Quad4[node + npt_el];
                Co[kdim * _SpaceDim + 1] += 0.25 * _CoordMatrix[_SpaceDim * node + kdim] * Quad4[node];
                Co[kdim * _SpaceDim + 2] += 0.25 * _CoordMatrix[_SpaceDim * node + kdim] * Quad4[node + npt_el];
                XYZm[kdim] += 0.25 * _CoordMatrix[_SpaceDim * node + kdim];
                }
        
            XYZs[kdim] = NodePos[kdim] - XYZm[kdim];
            }
        
        //
        int normal_max_dir = ( ( fabs ( normal[0] ) > fabs ( normal[1] ) ) ?   0 : 1 ) ;
        normal_max_dir = ( fabs ( normal[normal_max_dir] ) > fabs ( normal[2] ) ? normal_max_dir : 2 );
        int id1 = ( normal_max_dir + 1 ) % 3; // non normal direction 1
        int id2 = ( normal_max_dir + 2 ) % 3; // non normal direction 2
        
        // if  Co[1+3*id1] and Co[0+3*id1] are null -> inversion coeff
        if ( fabs ( Co[1 + 3 * id1] ) + fabs ( Co[0 + 3 * id1] ) < 1.e-10 )   {
            int itmp = id1;
            id1 = id2;
            id2 = itmp;
            std::cout << " index inversion " << std::endl;
            }
        
        double alpha_0 = Co[0 + 3 * id1] * Co[2 + 3 * id2] - Co[0 + 3 * id2] * Co[2 + 3 * id1];
        double beta_0 = -XYZs[id2] * Co[0 + 3 * id1] + XYZs[id1] * Co[0 + 3 * id2] + Co[1 + 3 * id1] * Co[2 + 3 * id2] - Co[1 + 3 * id2] * Co[2 + 3 * id1];
        double gamma_0 = -Co[1 + 3 * id1] * XYZs[id2] + Co[1 + 3 * id2] * XYZs[id1];
        
        
        if ( fabs ( alpha_0 ) > 1.e-10 ) { //  alpha_0 != 0 two roots 1.1e-10 ==========================================
            double delta_0 = beta_0 * beta_0 - 4 * alpha_0 * gamma_0;
        
//             if ( delta_0 < 1.e-10 ) {
//                 std::cout << " !!!!!!!!!!!!!!! Error delta " << std::endl;
//                 }
        
            // first root (x_00, y_00)
            xieta[1] = ( -beta_0 + sqrt ( delta_0 ) ) / ( 2.*alpha_0 );
            xieta[0] = - ( -XYZs[id1] + Co[2 + 3 * id1] * xieta[1] ) / ( Co[1 + 3 * id1] + Co[0 + 3 * id1] * xieta[1] );
        
            if ( fabs ( Co[1 + 3 * id1] + Co[0 + 3 * id1]*xieta[1] ) < 1e-10 ) {
                xieta[0] = - ( -XYZs[id2] + Co[2 + 3 * id2] * xieta[1] ) / ( Co[1 + 3 * id2] + Co[0 + 3 * id2] * xieta[1] );
                }
        
            // second root (x_01, y_01)
            if ( fabs ( xieta[0] ) > 1.001 ||  fabs ( xieta[1] ) > 1.001 ) { //  std::cout<<" first not good x1 "<< xieta[0] << " " << xieta[1] <<std::endl;
                xieta[1] = 0.5 * ( -beta_0 - sqrt ( delta_0 ) ) / alpha_0;
                xieta[0] = - ( -XYZs[id1] + Co[2 + 3 * id1] * xieta[1] ) / ( Co[1 + 3 * id1] + Co[0 + 3 * id1] * xieta[1] );
        
                if ( fabs ( Co[1 + 3 * id1] + Co[0 + 3 * id1]*xieta[1] ) < 1e-10 ) {
                    xieta[0] = - ( -XYZs[id2] + Co[2 + 3 * id2] * xieta[1] ) / ( Co[1 + 3 * id2] + Co[0 + 3 * id2] * xieta[1] );
                    }
                }
            }
        else {   // alpha_0 = 0 only oen root =====================================================
            xieta[1] = -gamma_0 / beta_0; //        std::cout << " one root only  " << std::endl;
            xieta[0] = - ( -XYZs[id1] + Co[2 + 3 * id1] * xieta[1] ) / ( Co[1 + 3 * id1] + Co[0 + 3 * id1] * xieta[1] );
            }
        
        for ( int dim = 0; dim < _MeshDim; dim++ ) {
            XiEtaBound.push_back ( xieta[dim] );
            }
            
        delete [] xieta;    
  }
   
  delete [] XiEtaChi;
  return;
  }



// ===============================================================================================
void InterfaceProjection::XiEtaCalc_1D (
  double NodePos[],                // physical coord (in)
  std::vector<double> & XiEtaBound, // reference coord (out)
  int  npt_el
) {// =============================================================================================

  double xm = 0.5 * ( _CoordMatrix[0] + _CoordMatrix[2] );
  double ym = 0.5 * ( _CoordMatrix[1] + _CoordMatrix[3] );
  double a  = 0.5 * ( -_CoordMatrix[0] + _CoordMatrix[2] );
  double b  = 0.5 * ( -_CoordMatrix[1] + _CoordMatrix[3] );
  double xi = ( NodePos[0] + NodePos[1] - ( xm + ym ) ) / ( a + b );
  XiEtaBound.push_back ( xi );
  return;
  };

int InterfaceProjection::IsNodeContained() {
  int Verify = 1;
  double toll = 1.e-3;

  if ( _FamilyType == 1 ) {
      if ( fabs ( _XiEtaChi[0] ) > 1 + toll ||
           fabs ( _XiEtaChi[1] ) > 1 + toll ||
           fabs ( _XiEtaChi[DIMENSION - 1] ) > 1 + toll ) Verify = 0;
      }

  if ( _FamilyType == 0 ) {
      if ( _XiEtaChi[0] < -toll || _XiEtaChi[0] > 1 + toll ||
           _XiEtaChi[1] < -toll || _XiEtaChi[1] > 1 + toll ||
           _XiEtaChi[DIMENSION - 1] < -toll || _XiEtaChi[DIMENSION - 1] > 1 + toll ) Verify = 0;
      }

  return Verify;
  }


void InterfaceProjection::IsNodeInsideCell ( double NodePos[] )
  {
  // Function written using MED nodal connectivity. For more info see
  // see https://docs.salome-platform.org/8/gui/SMESH/connectivity.html?highlight=connectivity 

  if ( _DiscOrder == 2 ) {
      if ( _MeshDim == 3 ) {
          int ins1 = 0, ins2 = 0;
          const double toll = 1.e-15;

//           int MidFaceNodesID[6] = {22, 24, 23, 21, 25, 20};
          const int maxTr = 4;
          const int numCoup = 3;

          int TrianCount = 0;
          int CoupCount = 0;
          int indx1, indx2, indx3;

          int FaceConn[6][4] = {
              {2, 1, 5, 6},
              {0, 3, 7, 4},
              {3, 2, 6, 7},
              {1, 0, 4, 5},
              {4, 7, 6, 5},
              {0, 1, 2, 3}
            };

          double dist = 1.;
          _Contained = 1;

          while ( _Contained && CoupCount < numCoup ) {
              ins1 = 0, ins2 = 0;
              TrianCount = 0;
              dist = 1.;

              while ( dist > 0. && TrianCount < maxTr ) {
                  indx2 = FaceConn[2 * CoupCount][TrianCount % 4];
                  indx1 = FaceConn[2 * CoupCount][ ( TrianCount + 1 ) % 4];
                  indx3 = FaceConn[2 * CoupCount][ ( TrianCount + 2 ) % 4];
                  dist = DistanceFromPlane ( NodePos, indx1, indx2, indx3 );

                  TrianCount ++;
                  }

              if ( dist <= toll )
                ins1 = 1;

              TrianCount = 0;
              dist = 1.;

              while ( dist > 0. && TrianCount < maxTr ) {
                  indx2 = FaceConn[2 * CoupCount + 1][TrianCount % 4];
                  indx1 = FaceConn[2 * CoupCount + 1][ ( TrianCount + 1 ) % 4];
                  indx3 = FaceConn[2 * CoupCount + 1][ ( TrianCount + 2 ) % 4];
                  dist = DistanceFromPlane ( NodePos, indx1, indx2, indx3 );

                  TrianCount ++;
                  }

              if ( dist <= toll )
                ins2 = 1;

              if ( ! ( ins1 * ins2 ) )
                _Contained = 0;

              CoupCount ++;
              }
          }
          if(_MeshDim == 2){
            if(_FamilyType==1){
              const int maxEd = 2;
              const int numCoup = 2;
              int ins1, ins2, indx1, indx2;
              const double toll = -1.e-15;
              
              int EdCount = 0;
              int CoupCount = 0;
              int FaceConn[4][3] = {
                 {0, 4, 1},
                 {1, 5, 2},
                 {2, 6, 3},
                 {3, 7, 0}
              };
              double dist = 1.;
              _Contained = 1;

              while ( _Contained && CoupCount < numCoup ) {
                  ins1 = 0, ins2 = 0;
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2 * CoupCount][EdCount % 3];
                      indx1 = FaceConn[2 * CoupCount][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins1 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2 * CoupCount + 1][EdCount % 3];
                      indx1 = FaceConn[2 * CoupCount + 1][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins2 = 1;
              
                  if ( ! ( ins1 * ins2 ) )
                    _Contained = 0;
              
                  CoupCount ++;
                  }                  
              }
              if(_FamilyType==0){
              const int maxEd = 2;
              int ins1, ins2, ins3, indx1, indx2;
              const double toll = -1.e-15;
              
              int EdCount = 0;
              int CoupCount = 0;
              int FaceConn[3][3] = {
                 {0, 3, 1},
                 {1, 4, 2},
                 {2, 5, 0}
              };
              double dist = 1.;
              _Contained = 1;


                  ins1 = 0, ins2 = 0;
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[0][EdCount % 3];
                      indx1 = FaceConn[0][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins1 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[1][EdCount % 3];
                      indx1 = FaceConn[1][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins2 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2][EdCount % 3];
                      indx1 = FaceConn[2][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins3 = 1;
                  
                  if ( ! ( ins1 * ins2 * ins3 ) )
                    _Contained = 0;
              
                  CoupCount ++;
               
              }   
          }
  }
  else {
    if ( _MeshDim == 3 ) {
          int ins1 = 0, ins2 = 0;
          const double toll = 1.e-15;

//           int MidFaceNodesID[6] = {22, 24, 23, 21, 25, 20};
          const int maxTr = 4;
          const int numCoup = 3;

          int TrianCount = 0;
          int CoupCount = 0;
          int indx1, indx2, indx3;

          int FaceConn[6][4] = {
              {2, 1, 5, 6},
              {0, 3, 7, 4},
              {3, 2, 6, 7},
              {1, 0, 4, 5},
              {4, 7, 6, 5},
              {0, 1, 2, 3}
            };

          double dist = 1.;
          _Contained = 1;

          while ( _Contained && CoupCount < numCoup ) {
              ins1 = 0, ins2 = 0;
              TrianCount = 0;
              dist = 1.;

              while ( dist > 0. && TrianCount < maxTr ) {
                  indx2 = FaceConn[2 * CoupCount][TrianCount % 4];
                  indx1 = FaceConn[2 * CoupCount][ ( TrianCount + 1 ) % 4];
                  indx3 = FaceConn[2 * CoupCount][ ( TrianCount + 2 ) % 4];
                  dist = DistanceFromPlane ( NodePos, indx1, indx2, indx3 );

                  TrianCount ++;
                  }

              if ( dist <= toll )
                ins1 = 1;

              TrianCount = 0;
              dist = 1.;

              while ( dist > 0. && TrianCount < maxTr ) {
                  indx2 = FaceConn[2 * CoupCount + 1][TrianCount % 4];
                  indx1 = FaceConn[2 * CoupCount + 1][ ( TrianCount + 1 ) % 4];
                  indx3 = FaceConn[2 * CoupCount + 1][ ( TrianCount + 2 ) % 4];
                  dist = DistanceFromPlane ( NodePos, indx1, indx2, indx3 );

                  TrianCount ++;
                  }

              if ( dist <= toll )
                ins2 = 1;

              if ( ! ( ins1 * ins2 ) )
                _Contained = 0;

              CoupCount ++;
              }
          }
          if(_MeshDim == 2){
            if(_FamilyType==1){
              const int maxEd = 1;
              const int numCoup = 2;
              int ins1, ins2, indx1, indx2;
              const double toll = -1.e-15;
              
              int EdCount = 0;
              int CoupCount = 0;
              int FaceConn[4][2] = {
                 {0, 1},
                 {1, 2},
                 {2, 3},
                 {3, 0}
              };
              double dist = 1.;
              _Contained = 1;

              while ( _Contained && CoupCount < numCoup ) {
                  ins1 = 0, ins2 = 0;
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2 * CoupCount][EdCount % 3];
                      indx1 = FaceConn[2 * CoupCount][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins1 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2 * CoupCount + 1][EdCount % 3];
                      indx1 = FaceConn[2 * CoupCount + 1][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins2 = 1;
              
                  if ( ! ( ins1 * ins2 ) )
                    _Contained = 0;
              
                  CoupCount ++;
                  }                  
              }
              if(_FamilyType==0){
              const int maxEd = 1;
              int ins1, ins2, ins3, indx1, indx2;
              const double toll = -1.e-15;
              
              int EdCount = 0;
              int CoupCount = 0;
              int FaceConn[3][3] = {
                 {0, 1},
                 {1, 2},
                 {2, 0}
              };
              double dist = 1.;
              _Contained = 1;


                  ins1 = 0, ins2 = 0;
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[0][EdCount % 3];
                      indx1 = FaceConn[0][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins1 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[1][EdCount % 3];
                      indx1 = FaceConn[1][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins2 = 1;
              
                  EdCount = 0;
                  dist = -1.;              
                  while ( dist < 0. && EdCount < maxEd ) {
                      indx2 = FaceConn[2][EdCount % 3];
                      indx1 = FaceConn[2][ ( EdCount + 1 ) % 3];
                      dist = DistanceFromEdge ( NodePos, indx1, indx2 );
              
                      EdCount ++;
                      }              
                  if ( dist >= toll )
                    ins3 = 1;
                  
                  if ( ! ( ins1 * ins2 * ins3 ) )
                    _Contained = 0;
              
                  CoupCount ++;
               
              }   
          }
  }
      

  return;
  }



void InterfaceProjection::XiEtaChiCalc ( double NodePos[] ) {

//     xyz, xyz, ...
//     _CoordMatrix[];


//     _Contained = 0;
//     IsNodeInsideCell(NodePos);

  if ( _Contained ) {
      double * XiEtaChi = new double[_MeshDim];

      for ( int dim = 0; dim < _MeshDim; dim++ ) {
          XiEtaChi[dim] = 0.;
          }

      double * deltaXEC = new double[_MeshDim];
      double SquaredError;
      const double Toll = pow ( 1e-2, 2 * _MeshDim );
      const double Toll2 = pow ( 5e-2, 2 * _MeshDim );

      bool FoundXiEtaChi = false;
      int count = 0;
      SquaredError = CalcF ( NodePos, XiEtaChi );
      double start_err = SquaredError;

      if ( SquaredError < Toll ) {
          FoundXiEtaChi = true;
          }

      double test_err = SquaredError;
      int test_id = _SrcCellNodes - 1;
      int test = 0;

      while ( !FoundXiEtaChi && test < _SrcCellNodes ) {
          if( _MeshDim==3 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordHex27[test + dim * 27];
          if( _MeshDim==2 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordQuad9[test + dim * 9];
          if( _MeshDim==1 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordEdge3[test + dim * 3];
          
          SquaredError = CalcF ( NodePos, XiEtaChi );

          if ( SquaredError < test_err ) {
              test_id = test;
              test_err = SquaredError;
              }

          if ( SquaredError < Toll ) FoundXiEtaChi = true;

          test ++;
          }

      if( _MeshDim==3 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordHex27[test_id + dim * 27];
      if( _MeshDim==2 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordQuad9[test_id + dim * 9];
      if( _MeshDim==1 )for ( int dim = 0; dim < _MeshDim; dim++ ) XiEtaChi[dim] = _CoordEdge3[test_id + dim * 3];

//   std::clock_t par_time1 = std::clock();

      while ( !FoundXiEtaChi ) {

          FirstDerF ( NodePos, XiEtaChi );
          SecondDerF ( NodePos, XiEtaChi );

          // Assembling the inverse of Hessian matrix
          double * InvHex = new double[_MeshDim * _MeshDim];
          double DetHex = 0.;//axx*ayy*azz + 2.*axy*ayz*axz - ayy*axz*axz - azz*axy*axy - axx*ayz*ayz;

          if ( _MeshDim == 3 ) {
              DetHex = _d2F[0] * _d2F[4] * _d2F[8] + 2.*_d2F[1] * _d2F[5] * _d2F[2] - _d2F[4] * _d2F[2] * _d2F[2] - _d2F[8] * _d2F[1] * _d2F[1] - _d2F[0] * _d2F[5] * _d2F[5];

              for ( int row = 0; row < _MeshDim; row++ ) {
                  for ( int col = 0; col < _MeshDim; col++ ) {
                      if ( row <= col ) {
                          int a = ( int ) ( 1 - ( ( row + 1 ) >> 1 ) );
                          int b = ( int ) ( 2 - ( row >> 1 ) );
                          int c = ( int ) ( 1 - ( ( col + 1 ) >> 1 ) );
                          int d = ( int ) ( 2 - ( col >> 1 ) );
                          int sign = ( int ) ( col + row ) % 2;
                          InvHex[row * _MeshDim + col] = ( 1 - 2 * sign ) * ( _d2F[a * _MeshDim + c] * _d2F[b * _MeshDim + d] - _d2F[a * _MeshDim + d] * _d2F[b * _MeshDim + c] ) / DetHex;
                          }
                      else {
                          InvHex[row * _MeshDim + col] = InvHex[col * _MeshDim + row];
                          }
                      }
                  }
              }

          // Per il caso 2D si possono provare anche metodi di soluzione analitica delle coordinate
          if ( _MeshDim == 2 ) {
              DetHex = _d2F[0] * _d2F[3] - _d2F[1] * _d2F[2];

              for ( int row = 0; row < _MeshDim; row++ ) {
                  for ( int col = 0; col < _MeshDim; col++ ) {
                      if ( row <= col ) {
                          int a = ( int ) ( 1 - ( ( row + 1 ) >> 1 ) );
                          int c = ( int ) ( 1 - ( ( col + 1 ) >> 1 ) );
                          int sign = ( int ) ( col + row ) % 2;
                          InvHex[row * _MeshDim + col] = ( 1 - 2 * sign ) * _d2F[a * _MeshDim + c] / DetHex;
                          }
                      else {
                          InvHex[row * _MeshDim + col] = InvHex[col * _MeshDim + row];
                          }
                      }
                  }
              }

          if ( _MeshDim == 1 ) {
              InvHex[0] = 1 / _d2F[0];
              }

          // Calculating the variation of xi, eta, chi
          for ( int dir = 0; dir < _MeshDim; dir ++ ) {
              double SUM = 0.;

              for ( int i = 0; i < _MeshDim; i++ ) {
                  SUM += InvHex[_MeshDim * dir + i] * _dF[i];
                  }

              deltaXEC[dir] = -0.5 * SUM;
              }

          for ( int i = 0; i < _MeshDim; i++ ) {
              XiEtaChi[i] += deltaXEC[i];
              }

          double old_err = SquaredError;
          SquaredError = CalcF ( NodePos, XiEtaChi );
          count ++;
          double rat = SquaredError / old_err;

          if ( SquaredError < Toll ) {
              FoundXiEtaChi = true;
              }

          if ( count > 50 ) {
              FoundXiEtaChi = true;
              }

//         if ( SquaredError < Toll2 || (rat>0.99&&rat<1.01) ){
//             _XiEtaChi.clear();
//             for ( int Mcomp=0; Mcomp<_MeshDim; Mcomp++ ) {
//               _XiEtaChi.push_back ( XiEtaChi[Mcomp] );
//             }
//             int contained = IsNodeContained();
//             if (contained == 0) FoundXiEtaChi = true;
//         }

          delete [] InvHex;
          _dF.clear();
          _d2F.clear();
          }

//       if ( FoundXiEtaChi )
//         std::cerr << " Process stopped, init err " << start_err << " final err " << SquaredError << " count: " << count << "   " << XiEtaChi[0] << " " << XiEtaChi[1] << " " << XiEtaChi[2];

      for ( int Mcomp = 0; Mcomp < _MeshDim; Mcomp++ ) {
          _XiEtaChi.push_back ( XiEtaChi[Mcomp] );
          }

      delete [] XiEtaChi;
      delete [] deltaXEC;

      }

  else {
      for ( int Mcomp = 0; Mcomp < _MeshDim; Mcomp++ ) {
          _XiEtaChi.push_back ( -1.5 );
          }
      }

  return;
  };

// First derivatives of linear 3D test functions
double InterfaceProjection::FirstDer ( double XiEtaChi[], int node, int dir ) {
//   _CoordHex27
  double value ;

  if ( _FamilyType == 1 ) {
//         value = _fe[2]->Rec_Lin_DPhi(_fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim, dir);
      value = _fe[2]->Rec_Quad_DPhi ( _fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim, dir );
      }

  if ( _FamilyType == 0 ) {
    if(_DiscOrder==2){
      if ( _MeshDim == 1 ) value = _fe[2]->Edge_Quad_DPhi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) value = _fe[2]->Tri_2d_QuadraticDerPhi ( node, XiEtaChi, dir );
      if ( _MeshDim == 3 ) value = _fe[2]->Tri_3d_QuadraticDerPhi ( node, XiEtaChi, dir );
    }
    else{
      if ( _MeshDim == 1 ) value = _fe[2]->Edge_Lin_DPhi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) value = _fe[2]->Tri_2d_LinearDerPhi( node, XiEtaChi, dir );
      if ( _MeshDim == 3 ) value = _fe[2]->Tri_3d_LinearDerPhi( node, XiEtaChi, dir );
    }
      }

  return value;
  }

// Second derivatives of linear 3D test functions
double InterfaceProjection::SecondDer ( double XiEtaChi[], int node ) {
  double value;
  // In 2D 0<= node <4 and offset = 4  ->  the eta coordinates are taken from node + 8, so node+_VertCoordOff
  value  = 0.25 * _XiEtaChiVert[node] * _XiEtaChiVert[node + _VertCoordOff]; // direction 0
  return value;
  }
double InterfaceProjection::SecondDer ( double XiEtaChi[], int node, int row, int col ) {

  double value;

  if ( _FamilyType == 1 ) {
//         value = _fe[2]->Rec_Lin_D2Phi(_fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim, row, col);
      value = _fe[2]->Rec_Quad_D2Phi ( _fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim, row, col );
      }

  if ( _FamilyType == 0 ) {
    if(_DiscOrder==2){
      if ( _MeshDim == 1 ) value = _fe[2]->Edge_Quad_D2Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) value = _fe[2]->Tri_2d_QuadraticDer2Phi ( node, XiEtaChi, row, col );
      if ( _MeshDim == 3 ) value = _fe[2]->Tri_3d_QuadraticDer2Phi ( node, XiEtaChi, row, col );
    }
    else{
      value = 0.;
//       if ( _MeshDim == 1 ) value = _fe[2]->Edge_Lin_D2Phi ( node, XiEtaChi[0] );
//       if ( _MeshDim == 2 ) value = _fe[2]->Tri_2d_LinearDer2Phi ( node, XiEtaChi, row, col );
//       if ( _MeshDim == 3 ) value = _fe[2]->Tri_3d_LinearDer2Phi ( node, XiEtaChi, row, col );
    }
      
      }

  return value;
  }


double InterfaceProjection::LinPhi ( int node, double XiEtaChi[] ) {
  double N = 1.;

  if ( _FamilyType == 1 ) {
      N = _fe[2]->Rec_Quad_Phi ( _fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim );
//         N = _fe[2]->Rec_Lin_Phi(_fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim);

      }

  if ( _FamilyType == 0 ) {
    if(_DiscOrder==2){
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Quad_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_QuadraticPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_QuadraticPhi ( node, XiEtaChi );
    }
    else{
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Lin_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_LinearPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_LinearPhi ( node, XiEtaChi );
    }
  }

  return N;
  }

double InterfaceProjection::LinPhi_bound ( int node, double XiEtaChi[] ) {
  double N = 1.;

  if ( _FamilyType == 1 ) {
      N = _fe[2]->Rec_Quad_Phi ( node, XiEtaChi, _MeshDim );
      }

  if ( _FamilyType == 0 ) {
    if(_DiscOrder==2){
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Quad_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_QuadraticPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_QuadraticPhi ( node, XiEtaChi );
    }
    else{
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Lin_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_LinearPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_LinearPhi ( node, XiEtaChi );
    }
  }

  return N;
  }  
  
double InterfaceProjection::InterpPhi ( int node, double XiEtaChi[], int InterpType ) {
  double N = 1.;

  if ( _FamilyType == 1 ) {
    if(InterpType==1){
      if ( _MeshDim == 1 ) N = _fe[2]->Rec_Lin_Phi ( node, XiEtaChi, _MeshDim );
      if ( _MeshDim == 2 ) N = _fe[2]->Rec_Lin_Phi ( node, XiEtaChi, _MeshDim );
      if ( _MeshDim == 3 ) N = _fe[2]->Rec_Lin_Phi ( _fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim );
    }
    else{
      if ( _MeshDim == 1 ) N = _fe[2]->Rec_Quad_Phi ( node, XiEtaChi, _MeshDim );
      if ( _MeshDim == 2 ) N = _fe[2]->Rec_Quad_Phi ( node, XiEtaChi, _MeshDim );
      if ( _MeshDim == 3 ) N = _fe[2]->Rec_Quad_Phi ( _fe[2]->_MedToLib_27[node], XiEtaChi, _MeshDim );
    }
  }

  if ( _FamilyType == 0 ) {
    if(InterpType==2){
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Quad_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_QuadraticPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_QuadraticPhi ( node, XiEtaChi );
    }
    else{
      if ( _MeshDim == 1 ) N = _fe[2]->Edge_Lin_Phi ( node, XiEtaChi[0] );
      if ( _MeshDim == 2 ) N = _fe[2]->Tri_2d_LinearPhi ( node, XiEtaChi );
      if ( _MeshDim == 3 ) N = _fe[2]->Tri_3d_LinearPhi ( node, XiEtaChi );
    }
  }

  return N;
  }  

// Function for the calculation of F, first derivatives of F and second derivatives of F where F is the quantity
//                 F = (xp - xi)*(xp - xi) + (yp - yi)*(yp - yi) + (zp - zi)*(zp - zi)
// (xp,yp,zp) are the coordinates of the point of interest P and (xi,yi,zi) are the one calculated with (xi,eta,chi)
// Refer to "Exact and efficient interpolation using finite elements shape functions" (here the formula is correct)
// ************************ F ******************************
double InterfaceProjection::CalcF ( double NodePos[], double XiEtaChi[] ) {
  double f = 0.;

  for ( int i = 0; i < _SpaceDim; i++ ) {
      double SUM = 0.;

      for ( int j = 0; j < _SrcCoordInterpNodes; j++ ) {
          SUM += _CoordMatrix[j * _SpaceDim + i] * LinPhi ( j, XiEtaChi );
          }

      f += ( NodePos[i] - SUM ) * ( NodePos[i] - SUM );
      }

  return f;
  }

double InterfaceProjection::CalcF_bound ( double NodePos[], double XiEtaChi[] ) {
  double f = 0.;

  for ( int i = 0; i < _SpaceDim; i++ ) {
      double SUM = 0.;

      for ( int j = 0; j < _SrcCellNodes; j++ ) {
          SUM += _CoordMatrix[j * _SpaceDim + i] * LinPhi_bound ( j, XiEtaChi );
          }

      f += ( NodePos[i] - SUM ) * ( NodePos[i] - SUM );
      }

  return f;
  }  
  
  
// *********************** dF ******************************
void InterfaceProjection::FirstDerF ( double NodePos[], double XiEtaChi[] ) {
  // Function for the calculation of first order derivatives of squared error function
  double Error[_SpaceDim];

  for ( int extdir = 0; extdir < _SpaceDim; extdir++ ) {
      Error[extdir] = 0.;
      double xyzInterp = 0.;

      for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
          xyzInterp += _CoordMatrix[nNode * _SpaceDim + extdir] * LinPhi ( nNode, XiEtaChi );
          }

      Error[extdir] += NodePos[extdir] - xyzInterp;
      }

  for ( int comp = 0; comp < _SpaceDim; comp++ ) {
      double derValue = 0.;

      for ( int dir = 0; dir < _SpaceDim; dir++ ) {
          double dxyzInterp = 0.;

          for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
              dxyzInterp += _CoordMatrix[nNode * _SpaceDim + dir] * FirstDer ( XiEtaChi, nNode, comp );
              }

          derValue += -2.*Error[dir] * dxyzInterp;
          }

      _dF.push_back ( derValue );
      }

  return;
  };

// *********************** d2F *****************************
void InterfaceProjection::SecondDerF ( double NodePos[], double XiEtaChi[] ) {

  // Function for the calculation of second order derivatives of squared error function
  // Vector _d2F contains the Hessian matrix of  F
  double Error[_SpaceDim];

  for ( int extdir = 0; extdir < _SpaceDim; extdir++ ) {
      Error[extdir] = 0.;
      double xyzInterp = 0.;

      for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
          xyzInterp += _CoordMatrix[nNode * _SpaceDim + extdir] * LinPhi ( nNode, XiEtaChi );
          }

      Error[extdir] += NodePos[extdir] - xyzInterp;
      }

  for ( int row = 0; row < _MeshDim; row++ ) {
      for ( int col = 0; col < _MeshDim; col++ ) {
          double DerValue = 0.;

          if ( row <= col ) {
              for ( int dir = 0; dir < _MeshDim; dir++ ) {
                  double d2xyzInterp = 0.;
                  double Der1 = 0.;
                  double Der2 = 0.;

                  for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
                      d2xyzInterp += _CoordMatrix[nNode * _SpaceDim + dir] * SecondDer ( XiEtaChi, nNode, row, col );
                      Der1 += _CoordMatrix[nNode * _SpaceDim + dir] * FirstDer ( XiEtaChi, nNode, row );
                      Der2 += _CoordMatrix[nNode * _SpaceDim + dir] * FirstDer ( XiEtaChi, nNode, col );
                      }

                  DerValue += -2.* ( Error[dir] * d2xyzInterp - Der1 * Der2 );
                  }
              }
          else {
              DerValue = _d2F[col * _MeshDim + row];
              }

          _d2F.push_back ( DerValue );
          }
      }

  return;
  }

MEDCoupling::MEDCouplingFieldDouble *
InterfaceProjection::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble * SourceField, const int order ) {

  const int NComp = SourceField->getNumberOfComponents();
  int order1 = order;

  MEDCoupling::DataArrayDouble * targetArray = MEDCoupling::DataArrayDouble::New();
  targetArray -> alloc ( _TrgNodes, NComp );
  std::string EqName = SourceField->getName();

  switch ( __Domain ) {
    case ( Boundary ) :

      // Loop over the target mesh nodes
      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          GetInterNodes ( iNode,  BoundingNodes );
          std::vector<double> XiEtaBound;
          GetXiEta ( iNode, XiEtaBound );
          double *CanPos = new double[_MeshDim];

          for ( int dim = 0; dim < _MeshDim; dim++ ) {
              CanPos[dim] = XiEtaBound[dim];
              }

          //qui
          std::map<int, int> MedLibmesh;
          BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
          double TrgValue = 0.;

          for ( int iComp = 0; iComp < NComp; iComp++ ) {
              TrgValue = 0.;

              if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
//             TrgValue = 0.;
                  for ( int phin = 0; phin < pow ( 2, _MeshDim ); phin ++ ) { /// only edge quad and hex
                      const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                      const double phi = InterpPhi ( phin, CanPos, order );
                      TrgValue +=  val * phi;
                      }
                  }
              else {
//             TrgValue = 0.;
                  for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                      double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                      double phi = InterpPhi ( phin, CanPos, order );
                      TrgValue +=  val * phi;
                      }
                  }

              targetArray->setIJ ( iNode, iComp, TrgValue );
              }

          BoundingNodes.clear();
          XiEtaBound.clear();
          delete [] CanPos;
          }

      break;

    case ( Volume ) :

      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          std::vector<double> XiEtaChi;
          double TrgValue;
          double *CanPos = new double[_MeshDim];
          GetInterNodes ( iNode,  BoundingNodes );
          GetXiEtaChi ( iNode, XiEtaChi );

          if ( BoundingNodes[0] == -1 ) {
              for ( int icomp = 0; icomp < NComp; icomp++ ) {
                  targetArray->setIJ ( iNode, icomp, 0. );  //(node id, comp, value)
                  }

              BoundingNodes.clear(); //       break;
              }
          else {
              for ( int dim = 0; dim < _MeshDim; dim++ ) {
                  CanPos[dim] = XiEtaChi[dim];
                  }

              double val, phi;

              for ( int iComp = 0; iComp < NComp; iComp++ ) {


                  if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < pow ( 2, DIMENSION ); phin ++ ) { /// only edge quad and hex
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = InterpPhi ( phin, CanPos, order );
                          TrgValue +=  val * phi;
                          }
                      }
                  else {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = InterpPhi ( phin, CanPos, order );
                          TrgValue +=  val * phi;
                          }
                      }

                  targetArray->setIJ ( iNode, iComp, TrgValue );
                  BoundingNodes.clear();
                  XiEtaChi.clear();
                  }
              }

              delete [] CanPos;
          }

      break;

      } // end switch

  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  f->setMesh ( __OutMesh );
  f->setArray ( targetArray );
  f->setName ( EqName );
  targetArray->decrRef();
  return f->deepCopy();
  }

//Modifica a InterpolatedField

MEDCoupling::MEDCouplingFieldDouble *
InterfaceProjection::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble * SourceField, const MEDCoupling::MEDCouplingFieldDouble * TargetField, const int order ) {

  const int NComp = SourceField->getNumberOfComponents();

// const DataArrayDouble *TargetFieldArray = TargetField->getArray();

  int order1 = order;
  MEDCoupling::DataArrayDouble * targetArray = MEDCoupling::DataArrayDouble::New();
  targetArray -> alloc ( _TrgNodes, NComp );
  std::string EqName = SourceField->getName();

  switch ( __Domain ) {
    case ( Boundary ) :

      // Loop over the target mesh nodes
      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          GetInterNodes ( iNode,  BoundingNodes );
          std::vector<double> XiEtaBound;
          GetXiEta ( iNode, XiEtaBound );

          if ( BoundingNodes[0] == -1 ) {
              for ( int icomp = 0; icomp < NComp; icomp++ ) {
                  double TargetFieldValue = TargetField->getIJ ( iNode, icomp );
                  targetArray->setIJ ( iNode, icomp, TargetFieldValue );
                  }

              BoundingNodes.clear(); //       break;
              }
          else {

              double CanPos[_MeshDim];

              for ( int dim = 0; dim < _MeshDim; dim++ ) {
                  CanPos[dim] = XiEtaBound[dim];
                  }

              //qui
              std::map<int, int> MedLibmesh;
              BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
              double TrgValue = 0.;

              for ( int iComp = 0; iComp < NComp; iComp++ ) {
                  TrgValue = 0.;

//         TrgValue = 0.;
//         for(int phin = 0; phin< __BoundNodesPerCell; phin ++) {
//           double val = SourceField->getIJ(BoundingNodes[phin], iComp);
//           double phi = QuadPhi(phin, CanPos);
//           TrgValue +=  val*phi;
//         }
                  if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
//             TrgValue = 0.;
                      for ( int phin = 0; phin < pow ( 2, _MeshDim ); phin ++ ) { /// only edge quad and hex
                          const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          const double phi = LinPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }
                  else {
//             TrgValue = 0.;
                      for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = QuadPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }

                  targetArray->setIJ ( iNode, iComp, TrgValue );
                  }

              BoundingNodes.clear();
              XiEtaBound.clear();
              }
          }

      break;

    case ( Volume ) :

      std::map<int, int> InterpCoordNodes;
      InterpCoordNodes[27] = 8;
      InterpCoordNodes[9] = 4;
      InterpCoordNodes[3] = 2;
      InterpCoordNodes[10] = 4;
      InterpCoordNodes[6] = 3;
      InterpCoordNodes[7] = 3;

      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          std::vector<double> XiEtaChi;
          double TrgValue;
          double CanPos[_MeshDim];
          GetInterNodes ( iNode,  BoundingNodes );
          GetXiEtaChi ( iNode, XiEtaChi );

          if ( BoundingNodes[0] == -1 ) {
              for ( int icomp = 0; icomp < NComp; icomp++ ) {
                  double TargetFieldValue = TargetField->getIJ ( iNode, icomp );
                  targetArray->setIJ ( iNode, icomp, TargetFieldValue );
                  }

              BoundingNodes.clear(); //       break;
              }
          else {
              for ( int dim = 0; dim < _MeshDim; dim++ ) {
                  CanPos[dim] = XiEtaChi[dim];
                  }

              double val, phi;

              for ( int iComp = 0; iComp < NComp; iComp++ ) {


                  if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < pow ( 2, DIMENSION ); phin ++ ) { /// only edge quad and hex
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = LinPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }
                  else {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = QuadPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }

                  targetArray->setIJ ( iNode, iComp, TrgValue );
                  BoundingNodes.clear();
                  XiEtaChi.clear();
                  }
              }
          }
      }

  MEDCoupling::MEDCouplingFieldDouble * f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  f->setMesh ( __OutMesh );
  f->setArray ( targetArray );
  f->setName ( EqName );
  targetArray->decrRef();
  return f;
  }

MEDCoupling::MEDCouplingFieldDouble *
InterfaceProjection::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble * SourceField, double DefaultValue, const int order ) {

  const int NComp = SourceField->getNumberOfComponents();

  int order1 = order;
  MEDCoupling::DataArrayDouble * targetArray = MEDCoupling::DataArrayDouble::New();
  targetArray -> alloc ( _TrgNodes, NComp );
  std::string EqName = SourceField->getName();

  switch ( __Domain ) {
    case ( Boundary ) :

      // Loop over the target mesh nodes
      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          GetInterNodes ( iNode,  BoundingNodes );
          std::vector<double> XiEtaBound;
          GetXiEta ( iNode, XiEtaBound );

          if ( BoundingNodes[0] == -1 ) {
              for ( int icomp = 0; icomp < NComp; icomp++ ) {
                  targetArray->setIJ ( iNode, icomp, DefaultValue );
                  }

              BoundingNodes.clear(); //       break;
              }
          else {
              double CanPos[_MeshDim];

              for ( int dim = 0; dim < _MeshDim; dim++ ) {
                  CanPos[dim] = XiEtaBound[dim];
                  }

              //qui
              std::map<int, int> MedLibmesh;
              BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
              double TrgValue = 0.;

              for ( int iComp = 0; iComp < NComp; iComp++ ) {
                  TrgValue = 0.;

//         TrgValue = 0.;
//         for(int phin = 0; phin< __BoundNodesPerCell; phin ++) {
//           double val = SourceField->getIJ(BoundingNodes[phin], iComp);
//           double phi = QuadPhi(phin, CanPos);
//           TrgValue +=  val*phi;
//         }
                  if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
//             TrgValue = 0.;
                      for ( int phin = 0; phin < pow ( 2, _MeshDim ); phin ++ ) { /// only edge quad and hex
                          const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          const double phi = LinPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }
                  else {
//             TrgValue = 0.;
                      for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = QuadPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }



                  targetArray->setIJ ( iNode, iComp, TrgValue );
                  }

              BoundingNodes.clear();
              XiEtaBound.clear();
              }
          }

      break;

    case ( Volume ) :

      std::map<int, int> InterpCoordNodes;
      InterpCoordNodes[27] = 8;
      InterpCoordNodes[9] = 4;
      InterpCoordNodes[3] = 2;
      InterpCoordNodes[10] = 4;
      InterpCoordNodes[6] = 3;
      InterpCoordNodes[7] = 3;

      for ( int iNode = 0; iNode < _TrgNodes; iNode ++ ) {
          // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
          std::vector<int> BoundingNodes;
          std::vector<double> XiEtaChi;
          double TrgValue;
          double CanPos[_MeshDim];
          GetInterNodes ( iNode,  BoundingNodes );
          GetXiEtaChi ( iNode, XiEtaChi );

          if ( BoundingNodes[0] == -1 ) {
              for ( int icomp = 0; icomp < NComp; icomp++ ) {
                  targetArray->setIJ ( iNode, icomp, DefaultValue );
                  }

              BoundingNodes.clear(); //       break;
              }
          else {
              for ( int dim = 0; dim < _MeshDim; dim++ ) {
                  CanPos[dim] = XiEtaChi[dim];
                  }

              double val, phi;

              for ( int iComp = 0; iComp < NComp; iComp++ ) {


                  if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < pow ( 2, DIMENSION ); phin ++ ) { /// only edge quad and hex
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = LinPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }
                  else {
                      TrgValue = 0.;

                      for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                          double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                          double phi = QuadPhi ( phin, CanPos );
                          TrgValue +=  val * phi;
                          }
                      }

                  targetArray->setIJ ( iNode, iComp, TrgValue );
                  BoundingNodes.clear();
                  XiEtaChi.clear();
                  }
              }
          }
      }

  MEDCoupling::MEDCouplingFieldDouble * f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  f->setMesh ( __OutMesh );
  f->setArray ( targetArray );
  f->setName ( EqName );
  targetArray->decrRef();
  return f;
  }

// Function that returns the value of a field (Field) on a desired point (coord) set by FillParametersProbe
double InterfaceProjection::InterpolatedFieldProbe ( const MEDCoupling::MEDCouplingFieldDouble * SourceField,  const int order ) {

  if ( _OutOfDomain == 1 ) {
      std::cout << "\033[038;5;" << 155 << ";1m "
                << "WARNING: Cell not found for target point. Probably the selected point is out of the domain \n \033[0m";
      }

  const int NComp = SourceField->getNumberOfComponents();

  MEDCoupling::DataArrayDouble * targetArray = MEDCoupling::DataArrayDouble::New();
  targetArray -> alloc ( 1, NComp );
  std::string EqName = SourceField->getName();

  std::map<int, int> InterpCoordNodes;
  InterpCoordNodes[27] = 8;
  InterpCoordNodes[9] = 4;
  InterpCoordNodes[3] = 2;
  InterpCoordNodes[10] = 4;
  InterpCoordNodes[6] = 3;
  InterpCoordNodes[7] = 3;
  // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target node
  std::vector<int> BoundingNodes;
  std::vector<double> XiEtaChi;
  double TrgValue;
  double CanPos[_MeshDim];
  GetInterNodes ( 0,  BoundingNodes );
  GetXiEtaChi ( 0, XiEtaChi );

  if ( BoundingNodes[0] != -1 ) {
      for ( int dim = 0; dim < _MeshDim; dim++ ) {
          CanPos[dim] = XiEtaChi[dim];
          }

      double val, phi;

      for ( int iComp = 0; iComp < NComp; iComp++ ) {


          if ( ( ( EqName == "NS0" || EqName == "FSI0" || EqName == "FSIA0" || EqName == "NSA0" ) && iComp == _SpaceDim ) || order == 1 ) {
              TrgValue = 0.;

              for ( int phin = 0; phin < pow ( 2, DIMENSION ); phin ++ ) { /// only edge quad and hex
                  double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                  double phi = LinPhi ( phin, CanPos );
                  TrgValue +=  val * phi;
                  }
              }
          else {
              TrgValue = 0.;

              for ( int phin = 0; phin < __BoundNodesPerCell; phin ++ ) {
                  double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                  double phi = QuadPhi ( phin, CanPos );
                  TrgValue +=  val * phi;
                  }
              }

          targetArray->setIJ ( 0, iComp, TrgValue );
          BoundingNodes.clear();
          XiEtaChi.clear();
          }
      }

  return TrgValue;
  }




double InterfaceProjection::QuadPhi ( int N, double GPoint[] ) {
  double value;

  // General formula for test function: N_i = (1-0.5*|xi_i|)*((2|xi_i|-1)*xi*xi + xi_i*xi + (1-|xi_i|)) for xi, eta, chi

  if ( _FamilyType == 1 ) {
      switch ( _MeshDim ) {
        case ( 1 ) :
          value = ( 1. - 0.5 * fabs ( _CoordEdge3[N] ) ) *
                  ( ( 2.*fabs ( _CoordEdge3[N] ) - 1. ) * GPoint[0] * GPoint[0] + _CoordEdge3[N] * GPoint[0] + ( 1. - fabs ( _CoordEdge3[N] ) ) );
          break;

        case ( 2 ) :
          value = 1.;

          for ( int dim = 0; dim < _MeshDim; dim++ ) value *= ( 1. - 0.5 * fabs ( _CoordQuad9[N + dim * _Quad9Off] ) ) *
                ( ( 2.*fabs ( _CoordQuad9[N + dim * _Quad9Off] ) - 1. ) * GPoint[dim] * GPoint[dim]   + _CoordQuad9[N + dim * _Quad9Off] * GPoint[dim] + ( 1. - fabs ( _CoordQuad9[N + dim * _Quad9Off] ) ) );

          break;

        case ( 3 ) :
          value = 1.;

          for ( int dim = 0; dim < _MeshDim; dim++ ) value *= ( 1. - 0.5 * fabs ( _CoordHex27[N + dim * _Hex27Off] ) ) *
                ( ( 2.*fabs ( _CoordHex27[N + dim * _Hex27Off] ) - 1. ) * GPoint[dim] * GPoint[dim]   + _CoordHex27[N + dim * _Hex27Off] * GPoint[dim] + ( 1. - fabs ( _CoordHex27[N + dim * _Hex27Off] ) ) );

          break;
          }
      }

  if ( _FamilyType == 0 ) {
      double lambda;
      lambda = 1 - ( GPoint[0] + GPoint[1] );
      if ( _MeshDim == 1 ) value = _fe[2]->Edge_Quad_Phi ( N, GPoint[0] );
      if ( _MeshDim == 2 ) value = _fe[2]->Tri_2d_QuadraticPhi ( N, GPoint );
      if ( _MeshDim == 3 ) value = _fe[2]->Tri_3d_QuadraticPhi ( N, GPoint );
  }

  return value;
  };


double InterfaceProjection::DistanceFromPlane ( double NodePos[], int node1, int node2, int node3 )
  {
  double dist;

  double xx[3] = {_CoordMatrix[node1 * 3 ], _CoordMatrix[node2 * 3], _CoordMatrix[node3 * 3]};
  double yy[3] = {_CoordMatrix[node1 * 3 + 1], _CoordMatrix[node2 * 3 + 1], _CoordMatrix[node3 * 3 + 1]};
  double zz[3] = {_CoordMatrix[node1 * 3 + 2], _CoordMatrix[node2 * 3 + 2], _CoordMatrix[node3 * 3 + 2]};

  // plane equation: ax + by + cz - d=0
  double a = 0, b = 0, c = 0, d = 0;

  for ( int i = 0; i < 3; i++ ) {
      a += yy[i] * zz[ ( i + 1 ) % 3] - zz[i] * yy[ ( i + 1 ) % 3] ;
      b += zz[i] * xx[ ( i + 1 ) % 3] - xx[i] * zz[ ( i + 1 ) % 3] ;
      c += xx[i] * yy[ ( i + 1 ) % 3] - yy[i] * xx[ ( i + 1 ) % 3] ;
      d += xx[i] * yy[ ( i + 1 ) % 3] * zz[ ( i + 2 ) % 3] - xx[i] * zz[ ( i + 1 ) % 3] * yy[ ( i + 2 ) % 3];
      }

  double dd = a * xx[1] + b * yy[1] + c * zz[1];
  double mod = sqrt ( a * a + b * b + c * c );

  dist = -1.* ( a * NodePos[0] + b * NodePos[1] + c * NodePos[2] - dd ) / mod;

  if ( fabs ( dist ) < 1.e-15 ) dist = 0.;

  return dist;
  }

double InterfaceProjection::DistanceFromEdge ( double NodePos[], int node1, int node2 )
  {// POSITIVE IF ORIENTED TORWARDS CELL CENTER
  double dist, dist_c;
  std::vector<double> Center;
  if(_FamilyType==1){
    Center.push_back(_CoordMatrix[8 * 2 ]);
    Center.push_back(_CoordMatrix[8 * 2 + 1]);
  }
  else{
    Center.push_back((_CoordMatrix[0] + _CoordMatrix[2] + _CoordMatrix[4])/3.);
    Center.push_back((_CoordMatrix[1] + _CoordMatrix[3] + _CoordMatrix[5])/3.);
  }
  
  double xx[2] = {_CoordMatrix[node1 * 2 ], _CoordMatrix[node2 * 2]};
  double yy[2] = {_CoordMatrix[node1 * 2 + 1], _CoordMatrix[node2 * 2 + 1]};


  // edge equation: ax + by + c=0
  double a = 0, b = 0, c = 0;
  
  a = -(yy[0] - yy[1]);
  b =  (xx[0] - xx[1]);
  c = -b*yy[1] - a*xx[1];
  
  double mod = sqrt ( a * a + b * b );

  dist = ( a * NodePos[0] + b * NodePos[1] + c ) / mod;
  dist_c = ( a * Center[0] + b * Center[1] + c ) / mod;
  
  if ((dist * dist_c) >= 0)
    dist = fabs(dist);
  else 
    dist = -1.*fabs(dist);

  return dist;
  }
  
  
  void InterfaceProjection::BuildInterpCoordMap()
{
  Couple nodes = {2,2};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_SEG2, nodes));
  
  nodes = {2,3};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_SEG3, nodes));
  
  nodes = {3,3};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_TRI3, nodes));
  
  nodes = {4,4};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_QUAD4, nodes));
  
  nodes = {3,6};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_TRI6, nodes));
  
  nodes = {3,6};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_TRI7, nodes));
  
  nodes = {4,8};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_QUAD8, nodes));
  
  nodes = {4,9};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_QUAD9, nodes));
  
  nodes = {4,4};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_TETRA4, nodes));
  
  nodes = {8,8};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_HEXA8, nodes));
  
  nodes = {4,10};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_TETRA10, nodes));
  
  nodes = {8,20};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_HEXA20, nodes));
  
  nodes = {8,27};
  _InterpCoordsMap.insert(std::make_pair(INTERP_KERNEL::NORM_HEXA27, nodes));

  return;  
}

void InterfaceProjection::TestInterpolation(int Order){
  
  const int nSrcNodes = __InMesh->getNumberOfNodes();
  MEDCoupling::DataArrayDouble * Array = MEDCoupling::DataArrayDouble::New();
  Array->alloc(nSrcNodes,1);
  double *Array_pointer = Array->getPointer();
  for(int i=0; i<nSrcNodes; i++)
    Array_pointer[i] = i;
  
  MEDCoupling::MEDCouplingFieldDouble * srcTestField = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
  srcTestField->setArray(Array);
  srcTestField->setMesh(__InMesh);
  srcTestField->setName("Test");
  
  MEDCoupling::MEDCouplingFieldDouble * trgTestField = InterpolatedField(srcTestField, Order);
  if(_proc==0){
    MEDCoupling::WriteField("RESU_MED/srcTestField_"+std::to_string(Order)+".med", srcTestField, true);
    MEDCoupling::WriteField("RESU_MED/trgTestField_"+std::to_string(Order)+".med", trgTestField, true);
  }
  
  srcTestField->decrRef();
  trgTestField->decrRef();
  Array->decrRef();
  
  return;
}

#endif
// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 


