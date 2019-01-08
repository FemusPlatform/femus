#include <set>
#include <cstring>
#include <ctime>
#include "InterfaceProjection.h"
// #include<map>


#ifdef HAVE_MED
namespace MEDCoupling {
class MEDLoader;
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}

#include "InterfaceFunctionM.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDLoader.hxx"
#include "BoundingBox.hxx"
#include "InterpKernelGeo2DNode.hxx"




// ================================================================================================
BoundInterp::BoundInterp() : MMed() {
  __Filled=0;
}


// ================================================================================================
BoundInterp::BoundInterp (
    const MEDCoupling::MEDCouplingUMesh * SourceMesh,
    const MEDCoupling::MEDCouplingUMesh * TargetMesh,
    DomainType vol_sur
) : MMed() {
    __Filled=0;
    _AlreadyInitialized=0;
    FillParameters ( SourceMesh,TargetMesh,vol_sur);
}

BoundInterp::BoundInterp (
    const MEDCoupling::MEDCouplingUMesh * SourceMesh,
    const MEDCoupling::MEDCouplingUMesh * TargetMesh,
    int procId,
    DomainType   vol_sur 
) : MMed() {

   __Filled=0;
    _AlreadyInitialized=0;
    setProcId(procId);
    FillParameters ( SourceMesh,TargetMesh,vol_sur );
}

// ================================================================================================
void BoundInterp::FillParameters (
    const MEDCoupling::MEDCouplingUMesh * SourceMesh,
    const MEDCoupling::MEDCouplingUMesh * TargetMesh,
    int DomainType,
    double XiEtaToll
) {

    // Important!
    // Update the map when interpolation between different element types is available  //
    std::map<int,int> InterpCoordNodes;
    InterpCoordNodes[27]=8;
    InterpCoordNodes[9] =4;
    InterpCoordNodes[3] =2;
    InterpCoordNodes[10]=4;
    InterpCoordNodes[6]=3;
    InterpCoordNodes[7]=3;


    if(__Filled==1){
       __InMesh->decrRef();
       __OutMesh->decrRef();
    }
    __InMesh = SourceMesh->deepCopy();
    __OutMesh = TargetMesh->deepCopy();

    
    __Filled = 1;
    
    string name = SourceMesh->getName();
    _SrcCellNodes        = __InMesh->getNumberOfNodesInCell ( 0 );
    _SrcCoordInterpNodes = InterpCoordNodes[_SrcCellNodes];
    _SrcCells = __InMesh -> getNumberOfCells();
    
    
    
    INTERP_KERNEL::NormalizedCellType Type = __InMesh->getTypeOfCell ( 0 );
    if ( Type==6 || Type==7 || Type==14 || Type==20 ) _FamilyType=0;
    else _FamilyType=1;
    
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
    XiEta1->alloc ( _TrgNodes,_MeshDim );
    XiEta1->fillWithValue ( 0. );

    MEDCoupling::DataArrayInt * BoundingNodes1 = MEDCoupling::DataArrayInt::New();
    BoundingNodes1->alloc ( _TrgNodes,_SrcCellNodes );
    BoundingNodes1->fillWithValue ( -1 );

    double Bbox[2*_SpaceDim];
    SourceMesh->getBoundingBox ( Bbox );
    double height = fabs ( Bbox[3] - Bbox[2] );
    double width =  fabs ( Bbox[1] - Bbox[0] );

    double pointA[2];
    pointA[0] = Bbox[0] - 1.e4*width;
    pointA[1] = Bbox[2] - 1.e4*height;
    double pointB[2];
    pointB[0] = Bbox[0] - 1.e4*width;
    pointB[1] = Bbox[3] + 1.e4*height;

    MEDCoupling::DataArrayDouble *DDA = MEDCoupling::DataArrayDouble::New();
    DDA -> alloc ( _SrcCells,4 ); //  4 ????????????

//   MEDCoupling::DataArrayDouble *DDB = MEDCoupling::DataArrayDouble::New();
//   DDB -> alloc(_SrcCells,2);
//   std::vector<double> DDAM;   std::vector<double> DDAm; std::vector<double> DDBM;  std::vector<double> DDBm;

    const double coeffA = height/width;
    const double coeffB = -height/width;
    double coeffAm = ( Bbox[2] - pointA[1] ) / ( Bbox[1]-pointA[0] );
    double coeffAM = ( Bbox[3] - pointA[1] ) / ( Bbox[0]-pointA[0] );
    double coeffBm = ( Bbox[2] - pointB[1] ) / ( Bbox[0]-pointB[0] );
    double coeffBM = ( Bbox[3] - pointB[1] ) / ( Bbox[1]-pointB[0] );
    double deltaAm = ( coeffA - coeffAm ) /3.;
    double deltaAM = ( coeffAM - coeffA ) /3.;
    double deltaBm = fabs ( coeffA - coeffAm ) /3.;
    double deltaBM = fabs ( coeffAM - coeffA ) /3.;

//   std::cout<< "deltaAm " << deltaAm << " deltaAM " << deltaAM <<std::endl;
//   std::cout<< "deltaBm " << deltaAm << " deltaBM " << deltaAM <<std::endl;

    int SecMatrix[_SrcCells][4];


    MEDCoupling::DataArrayInt *targetArray = MEDCoupling::DataArrayInt::New();
    targetArray -> alloc ( _TrgNodes,1 );
    targetArray -> fillWithValue ( 1 );

    MEDCoupling::DataArrayDouble *CellBelonging = MEDCoupling::DataArrayDouble::New();
    CellBelonging -> alloc ( _TrgNodes,1 );
    CellBelonging -> fillWithValue ( -1 );
//   DataArrayInt * BoundingNodes1 = DataArrayInt::New();
//   DataArrayDouble * XiEta1 = DataArrayDouble::New();
//
//   XiEta1->alloc(_TrgNodes,_MeshDim);
//   BoundingNodes1->alloc(_TrgNodes,_SrcCellNodes);

    //***************************************************************
    // definizioni
    std::vector<double> coord;                              // Vector containing cell node coordinates
    std::vector< int > TrgConn;                             // Vector containing the cell node ids
    double *PointsCoords=new double [TrgCellNodes*DIMENSION];      // Node coord*/inates {xyz xyz xyz}
//     int *nodalConnPerCell = new int[TrgCellNodes];
    double TCbbox[2*DIMENSION];
    double minmax[2*DIMENSION];
    if ( _AlreadyInitialized==0 ) {
        InitFe();
        _AlreadyInitialized=1;
    }
    const int Quad4[8] = {-1, 1, 1,-1, -1,-1, 1, 1  }; // (xi,eta,zeta)
    const int Edge2[2] = {-1, 1}; // (xi,eta,zeta)
    switch ( DomainType ) {
    case ( Boundary ) :
        // ********* Loop over the target  mesh cells ************** //
        for ( int iCell=0; iCell < _TrgCells; iCell++ ) {

            __OutMesh -> getNodeIdsOfCell ( iCell,TrgConn );

            for ( int i_node=0; i_node < TrgCellNodes; i_node++ ) {             // Loop over the element nodes
                __OutMesh -> getCoordinatesOfNode ( TrgConn[i_node],coord );
                for ( int j =0; j<_SpaceDim; j++ ) {
                    PointsCoords[i_node*_SpaceDim + j] = coord[j];
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
            MEDCoupling::DataArrayDouble *coordsArr=MEDCoupling::DataArrayDouble::New();
            coordsArr->alloc ( _SrcCellNodes,_SpaceDim );
            std::copy ( PointsCoords,PointsCoords + TrgCellNodes*_SpaceDim,coordsArr->getPointer() );
            coordsArr->getMinMaxPerComponent ( TCbbox );
            coordsArr->decrRef();

            MEDCoupling::DataArrayInt * vettore = MEDCoupling::DataArrayInt::New();
            vettore = __InMesh->getCellsInBoundingBox ( TCbbox,1.e-5 );
            int SrcIntCells = vettore->getNbOfElems();                       // Number of source mesh cells intersecting target mesh cell bounding box
            MEDCoupling::MemArray<int> SrcIntCellsIds = vettore->accessToMemArray();      // Array containing source mesh cell ids


            double **ExtrCoord=new double*[SrcIntCells];                      // Matrix containing Min and Max coordinates of source mesh cells intersecting target mesh cell bounding box
            for ( int i = 0; i < SrcIntCells; ++i ) {
                ExtrCoord[i] = new double [2*_SpaceDim] ;
            }
            // cell 0: xmin, xmax, ymin, ymax, zmin, zmax
            /// Cycle over the source mesh cells intersecting the target mesh cell bounding box
            std::vector< int > SrcConn;
            for ( int IntCellNum = 0; IntCellNum < SrcIntCells; IntCellNum++ ) {

                __InMesh ->  getNodeIdsOfCell ( SrcIntCellsIds[IntCellNum],SrcConn );
                int NNodes = SrcConn.size();
                MEDCoupling::DataArrayDouble *XYZ=MEDCoupling::DataArrayDouble::New();
                XYZ->alloc ( NNodes,_SpaceDim );
                /// Min and Max coordinates of source mesh cells intersecting target mesh cell bounding box
                for ( int i_node=0; i_node < NNodes; i_node++ ) {              // Loop over the element nodes
                    std::vector< double > ElCoord;
                    __InMesh -> getCoordinatesOfNode ( SrcConn[i_node],ElCoord );
                    for ( int dimension = 0; dimension < _SpaceDim; dimension ++ ) {
                        XYZ->setIJ ( i_node,dimension , ElCoord[dimension] );
                    }
                    ElCoord.clear();
                }

                double minmax[2*_SpaceDim];
                XYZ->getMinMaxPerComponent ( minmax );
                for ( int i = 0; i<2*_SpaceDim; i++ ) {
                    ExtrCoord[IntCellNum][i] = minmax[i];
                }
                SrcConn.clear();
                XYZ->decrRef();
            }



            // Determination of source mesh cell containing the target mesh node
            int IdSrcCell;
            double NodeCoord[DIMENSION];     // Vector containing target mesh node coordinates
            // Loop over the target mesh cell nodes
            for ( int j = 0; j< TrgCellNodes; j++ ) {
                if ( targetArray->getIJ ( TrgConn[j],0 ) >0.5 ) {
                    std::vector<int> CellIds;                             // Vector where we store the Ids of cells that can contain the target mesh node j

                    for ( int IntCellNum = 0; IntCellNum < SrcIntCells; IntCellNum ++ ) {
                        if ( _MeshDim==2 ) {
                            /// verifica dell'appartenenza del nodo j alla cella valori[conf]
                            if ( //
                                PointsCoords[j*_SpaceDim]   >= ExtrCoord[IntCellNum][0]-1.e-9 && PointsCoords[j*_SpaceDim]   <= ExtrCoord[IntCellNum][1]+1.e-9 &&
                                PointsCoords[j*_SpaceDim+1] >= ExtrCoord[IntCellNum][2]-1.e-9 && PointsCoords[j*_SpaceDim+1] <= ExtrCoord[IntCellNum][3]+1.e-9 &&
                                PointsCoords[j*_SpaceDim+2] >= ExtrCoord[IntCellNum][4]-1.e-9 && PointsCoords[j*_SpaceDim+2] <= ExtrCoord[IntCellNum][5]+1.e-9
                            ) {
                                CellIds.push_back ( SrcIntCellsIds[IntCellNum] );
                            }
                        }
                        if ( _MeshDim==1 ) {
                            /// verifica dell'appartenenza del nodo j alla cella valori[conf]
                            if ( //
                                PointsCoords[j*_SpaceDim]   >= ( ExtrCoord[IntCellNum][0]-1.e-9 ) && PointsCoords[j*_SpaceDim]   <= ( ExtrCoord[IntCellNum][1]+1.e-9 ) &&
                                PointsCoords[j*_SpaceDim+1] >= ( ExtrCoord[IntCellNum][2]-1.e-9 ) && PointsCoords[j*_SpaceDim+1] <= ( ExtrCoord[IntCellNum][3]+1.e-9 )
                            ) {
                                CellIds.push_back ( SrcIntCellsIds[IntCellNum] );
                            }
                        }
                    }
                    int numberOfCell = CellIds.size();



                    std::vector< int > SrcNodesConn;    // connectivity to be filled from getNodeIdsOfCell med function
                    for ( int direction = 0; direction<_SpaceDim; direction++ ) {
                        NodeCoord[direction] = PointsCoords[j*_SpaceDim+direction];    // real coordinates
                    }
                    std::vector<double> XiEtaBound;                               // Vector containing target mesh node reference coordinates
                    // If there are two possible cells containing the target mesh node we search the one having fabs(xi)<=1 and fabs(eta)<=1

                    if ( numberOfCell  > 1 ) { //  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        int FirstTry = 0;
                        bool CellFound = false;
                        while ( !CellFound && FirstTry<numberOfCell ) { // ++++++++ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            IdSrcCell = CellIds[FirstTry];
                            __InMesh ->  getNodeIdsOfCell ( IdSrcCell,SrcNodesConn );
                            for ( int SNode1 = 0; SNode1<_SrcCoordInterpNodes; SNode1 ++ ) {
                                std::vector< double > ElCoord1;
                                __InMesh -> getCoordinatesOfNode ( SrcNodesConn[SNode1],ElCoord1 );
                                for ( int direction = 0; direction<_SpaceDim; direction++ ) {
                                    _CoordMatrix.push_back ( ElCoord1[direction] );
                                }
                                ElCoord1.clear();
                            }
                            if ( _MeshDim==2 ) XiEtaCalc_2D ( NodeCoord, XiEtaBound,Quad4 );
                            if ( _MeshDim==1 ) XiEtaCalc_1D ( NodeCoord, XiEtaBound );
                            
                            if ( _MeshDim==2 ) { // ---------------------------------------------------------------------------------------------
                                if ( fabs ( XiEtaBound[0] ) <= 1.+1.e-5 &&  fabs ( XiEtaBound[1] ) <= 1.+1.e-5 ) {
                                    CellFound = true;
                                    CellBelonging->setIJ ( TrgConn[j],0,IdSrcCell );
                                    std::cout<< " cell found for node *********** " <<j<<" of cell "<<iCell<<" calculated xi: "<<XiEtaBound[0]<<" eta: "<<XiEtaBound[1]<<std::endl;
                                    for ( int node = 0; node<_SrcCellNodes; node++ ) {
                                        BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                                    }
                                } else {
                                    std::cout<< " cell not found for node ************ " <<j<<" of cell "<<iCell<<" after "<<FirstTry+1 <<" out of "<<numberOfCell <<" attempts "
                                    <<" xi: "<<XiEtaBound[0]<<" eta: "<<XiEtaBound[1]<<std::endl;
                                    CellFound = false;
                                    FirstTry ++;
                                    if ( FirstTry!=numberOfCell ) {
                                        XiEtaBound.clear();    // clear the coords at the point
                                    }

                                }
                            }  // ----------------------------------------------------------------------------------------------------------
                            if ( _MeshDim==1 ) { // ---------------------------------------------------------------------------------------------
                                if ( fabs ( XiEtaBound[0] ) <= 1.+1.e-5 ) {
                                    CellFound = true;
                                    for ( int node = 0; node<_SrcCellNodes; node++ ) {
                                        BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                                    }
                                } else {
                                    CellFound = false;
                                    FirstTry ++;
                                    XiEtaBound.clear();
                                }
                            } // ---------------------------------------------------------------------------------------------------------------

                            SrcNodesConn.clear();   // clear connectivity
                            _CoordMatrix.clear();
                        } // end while ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       if(!CellFound) std::cout<< " cell not found for node ************ " <<j<<" of cell "<<iCell<<" after "<<FirstTry+1 <<" out of "<<numberOfCell <<" attempts "
                           <<" xi: "<<XiEtaBound[0]<<std::endl;
//             else
                    } // end if numeroCelle > 1.5 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    else {
                        IdSrcCell = CellIds[0];
                        __InMesh ->  getNodeIdsOfCell ( IdSrcCell,SrcNodesConn );
                        CellBelonging->setIJ ( TrgConn[j],0,IdSrcCell );
                        int srcdim = SrcNodesConn.size();
                        std::map <int, int> Ordering;
                        for ( int nMed = 0; nMed < SrcNodesConn.size(); nMed++ )   if ( _MeshDim==2 ) {
                                Ordering[nMed] = nMed;
                            }

                        if ( _MeshDim==1 ) {
                            Ordering[0] = 0;
                            Ordering[1] = 1;
                            Ordering[2] = 2;
                        }
                        _CoordMatrix.clear();
                        for ( int SNode1 = 0; SNode1<_SrcCoordInterpNodes; SNode1 ++ ) {
                            std::vector< double > ElCoord1;
                            __InMesh -> getCoordinatesOfNode ( SrcNodesConn[Ordering[SNode1]],ElCoord1 );
                            for ( int direction = 0; direction<_SpaceDim; direction++ ) {
                                _CoordMatrix.push_back ( ElCoord1[direction] );
                            }
                            ElCoord1.clear();
                        }
                        if ( _MeshDim==2 ) {
                            XiEtaCalc_2D ( NodeCoord, XiEtaBound,Quad4 );
                        }
                        if ( _MeshDim==1 ) {
                            XiEtaCalc_1D ( NodeCoord, XiEtaBound );
                        }
//             std::cout<<" bob calculated xi: "<<XiEtaBound[0]<<" eta: "<<XiEtaBound[1]<<std::endl;
                        for ( int node = 0; node<_SrcCellNodes; node++ ) {
                            BoundingNodes1->setIJ ( TrgConn[j], node, SrcNodesConn[node] );
                        }
                        SrcNodesConn.clear();
                        _CoordMatrix.clear();
                    }


                    for ( int dim=0; dim<_MeshDim; dim++ ) {
                        XiEta1->setIJ ( TrgConn[j],dim,XiEtaBound[dim] );
                    }
                    targetArray->setIJ ( TrgConn[j],0,0 );
                    XiEtaBound.clear();
                    CellIds.clear();
                }

            }
            TrgConn.clear();
//       for (int i = 0; i < SrcIntCells; ++i)   delete [] ExtrCoord[i];
//         delete [] ExtrCoord;
        } // end loop

//     MEDCouplingFieldDouble * CellBelonging = MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
//     CellBelonging->setMesh(__OutMesh);
//     CellBelonging->setArray(CellBelonging);
//     CellBelonging->setName("CellID");


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
        for ( int dim=0; dim<dimm; dim++ ) {
            Nodi[dim] = CoordinatesOfNodes[dim];    //
        }

        std::map<int,int> MedLibmesh;
        BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );

//******************************************************************************
// controllare metodo di ricerca ed utilizzo degli array elts e eltsIndex -> guarda classe modificata da daniele
        __InMesh->getCellsContainingPoints ( Nodi,_TrgNodes,1.e-5,elts,eltsIndex );
//     DataArrayInt * cells = elts
        std::vector<double> NodeCoord;
        double pos[_SpaceDim];
        std::vector<int> SourceConn;
        std::vector< double > ElCoord1;
//     for(int TNode = 0; TNode < _TrgNodes; TNode++) {
//
//
//       __OutMesh->getCoordinatesOfNode(TNode, NodeCoord);
//       for(int dir=0; dir<_SpaceDim; dir++) pos[dir] = NodeCoord[dir];
//       int NewCell = elts->getIJ(eltsIndex->getIJ(TNode,0),0);
// // 	std::cout<<"node "<<TNode<<" cell "<<NewCell<<std::endl;
//
//       __InMesh->getNodeIdsOfCell(NewCell,SourceConn);
//       for(int SNode = 0; SNode<_SrcCellNodes; SNode++) BoundingNodes1->setIJ(TNode, SNode,  SourceConn[SNode]);
//       for(int SNode1 = 0; SNode1<_SrcCoordInterpNodes; SNode1 ++) {
//         __InMesh -> getCoordinatesOfNode(SourceConn[SNode1],_CoordMatrix);
//       }
//
//       XiEtaChiCalc(pos);
//       for(int Mcomp=0; Mcomp<_MeshDim; Mcomp++) XiEta1->setIJ(TNode,Mcomp,_XiEtaChi[Mcomp]);
//       _CoordMatrix.clear();
//       _XiEtaChi.clear();
//       NodeCoord.clear();
//       SourceConn.clear();
//     }
        int count=0;
        for ( int TNode = 0; TNode < _TrgNodes; TNode++ ) {
            __OutMesh->getCoordinatesOfNode ( TNode, NodeCoord );
            for ( int dir=0; dir<_SpaceDim; dir++ ) {
                pos[dir] = NodeCoord[dir];
            }

            int NumPossibleCells = eltsIndex->getIJ ( TNode+1,0 ) - eltsIndex->getIJ ( TNode,0 );
//       int NumPossibleCells =eltsIndex[TNode + 1] - eltsIndex[TNode];

            bool found = false;
            if ( NumPossibleCells>0 ) {
                int iCount =0;
                while ( !found && iCount < NumPossibleCells ) {
                    _XiEtaChi.clear();
                    // cell ids are elts[ eltsIndex[ i ]],..., elts[ eltsIndex[ i ] + NumPossibleCells ].
                    int NewCell = elts->getIJ ( eltsIndex->getIJ ( TNode,0 ) + iCount,0 );
                    __InMesh->getNodeIdsOfCell ( NewCell,SourceConn );

                    for ( int SNode = 0; SNode<_SrcCellNodes; SNode++ ) {
//             BoundingNodes1->setIJ(TNode, SNode, SourceConn[MedLibmesh[SNode]]); //SourceConn[MedLibmesh[SNode]]=med index
//             if(SNode <_SrcCoordInterpNodes)  __InMesh -> getCoordinatesOfNode(SourceConn[MedLibmesh[SNode]],_CoordMatrix); // only linear nodes
                        BoundingNodes1->setIJ ( TNode, SNode, SourceConn[SNode] ); //SourceConn[MedLibmesh[SNode]]=med index
                        CellBelonging->setIJ ( TNode,0,NewCell );
                        if ( SNode <_SrcCoordInterpNodes ) {
                            __InMesh -> getCoordinatesOfNode ( SourceConn[SNode],_CoordMatrix );    // only linear nodes
                        }
                    }
                    XiEtaChiCalc ( pos );
                    int Contained = IsNodeContained();
                    if ( Contained==0 ) {
                        found=false;
                        _XiEtaChi.clear();
                        _CoordMatrix.clear();
                        SourceConn.clear();
                    } else if ( Contained==1 ) found = true;
                    iCount++;
                }// end while

                if ( !found ) {
                    for ( int SNode = 0; SNode<_SrcCellNodes; SNode++ ) BoundingNodes1->setIJ ( TNode, SNode, -1 ); //SourceConn[MedLibmesh[SNode]]=med index
                } else {

                    for ( int Mcomp=0; Mcomp<_MeshDim; Mcomp++ ) {
                        XiEta1->setIJ ( TNode,Mcomp,_XiEtaChi[Mcomp] );
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

        break;
    }

    delete [] PointsCoords;
   
    targetArray->decrRef(); 
    DDA->decrRef(); 

    MEDCoupling::MEDCouplingFieldDouble * CanonicalPosition = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
    CanonicalPosition->setMesh ( __OutMesh );
    CanonicalPosition->setArray ( XiEta1 );
    CanonicalPosition->setName ( "XiEta" );

    MEDCoupling::MEDCouplingFieldDouble * CellS = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
    CellS->setMesh ( __OutMesh );
    CellS->setArray ( CellBelonging );
    CellS->setName ( "CellID" );

    string FileName="S"+__InMesh->getName() +"_T"+__OutMesh->getName() +"_MeshCoupling.med";

    if(_proc==0){
      MEDCoupling::WriteUMesh ( "RESU_MED/"+FileName,__OutMesh,true );
      MEDCoupling::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/"+FileName,CanonicalPosition );
      MEDCoupling::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/"+FileName,CellS );
    }
    CanonicalPosition->decrRef();
    CellS->decrRef();

    if(_BoundingNodes != NULL) _BoundingNodes->decrRef();
    _BoundingNodes = BoundingNodes1;
    if(_XiEta != NULL) _XiEta->decrRef();
    _XiEta = XiEta1;

    __BoundNodesPerCell = _SrcCellNodes;
    __TrgNodes = _TrgNodes;
    __Domain = DomainType;
    
    CellBelonging->decrRef(); 
    return;
}
// ================================================================================================
BoundInterp::~BoundInterp() {
    __InMesh->decrRef();
    __OutMesh->decrRef();
    _XiEta->decrRef();
    _BoundingNodes->decrRef();
}


// ================================================================================================
void BoundInterp::CheckBelonging ( bool &found ) {
    const double toll = 1.e-20;
    const double xmedium = ( _CoordMatrix[0] + _CoordMatrix[2] + _CoordMatrix[4] + _CoordMatrix[6] ) /4.;
    const double ymedium = ( _CoordMatrix[1] + _CoordMatrix[3] + _CoordMatrix[5] + _CoordMatrix[7] ) /4.;
    found = true;
    double S;
    double tildeS;
    double xm, ym, tildexm,tildeym, midX, midY ;
    for ( int dim=0; dim<_SrcCoordInterpNodes; dim++ ) {
        int A = dim;
        int B = ( 4- ( dim+1 ) ) %4;

        xm = 0.5* ( xmedium - _pos[0] );
        xm = ( fabs ( xm<1.e-9 ) ?1.e-9:xm );
        ym = 0.5* ( ymedium - _pos[1] );
        ym = ( fabs ( ym<1.e-9 ) ?1.e-9:ym );
        tildexm = 0.5* ( _CoordMatrix[B*_SpaceDim] - _CoordMatrix[A*_SpaceDim] );
        tildexm = ( fabs ( tildexm<1.e-9 ) ?1.e-9:tildexm );
        tildeym = 0.5* ( _CoordMatrix[B*_SpaceDim+1] - _CoordMatrix[A*_SpaceDim+1] );
        tildeym = ( fabs ( tildeym<1.e-9 ) ?1.e-9:tildeym );
        midX = 0.5* ( _pos[0]+xmedium - ( _CoordMatrix[B*_SpaceDim]+_CoordMatrix[A*_SpaceDim] ) );
        midY = 0.5* ( _pos[1]+ymedium - ( _CoordMatrix[B*_SpaceDim+1]+_CoordMatrix[A*_SpaceDim+1] ) );
        S = - ( midX - midY*tildexm/tildeym ) / ( xm-tildexm/tildeym );
        tildeS = S* ( midY + ym ) /tildeym;

        if ( fabs ( fabs ( S )-1. ) <toll || fabs ( fabs ( tildeS )-1. ) <toll ) {
            found=true;    //coincident nodes
            break;
        }
        if ( ( fabs ( S ) +toll ) <1.&& ( fabs ( tildeS ) +toll ) <1. ) {
            found = false;    //cross
            break;
        }
    }
    return;
}


// ================================================================================================
void BoundInterp::SetSubSector ( int sector, int iCell ) {
    if ( _subsecAM != _subsecAm ) {
        if ( _subsecBM != _subsecBm ) {

            * ( _Scount[sector]+_subsecAm+_subsecBm ) +=  1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBm ) ) *9+_subsecAm+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            }
#endif

            * ( _Scount[sector]+_subsecAM+_subsecBm ) +=1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAM+_subsecBm ) ) *9+_subsecAM+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            }
#endif

            * ( _Scount[sector]+_subsecAm+_subsecBM ) +=1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBM ) ) *9+_subsecAm+_subsecBM ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            }
#endif

            * ( _Scount[sector]+_subsecAM+_subsecBM ) +=1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAM+_subsecBM ) ) *9+_subsecAM+_subsecBM ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAM+_subsecBM,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAM+_subsecBM,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAM+_subsecBM,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAM+_subsecBM,1 );
            }
#endif

        } else {
            * ( _Scount[sector]+_subsecAm+_subsecBm ) +=  1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBm ) ) *9+_subsecAm+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            }
#endif

            * ( _Scount[sector]+_subsecAM+_subsecBm ) +=1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAM+_subsecBm ) ) *9+_subsecAM+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAM+_subsecBm,1 );
            }
#endif
        }
    } else {
        if ( _subsecBM != _subsecBm ) {
            * ( _Scount[sector]+_subsecAm+_subsecBm ) +=  1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBm ) ) *9+_subsecAm+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            }
#endif

            * ( _Scount[sector]+_subsecAm+_subsecBM ) +=1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBM ) ) *9+_subsecAm+_subsecBM ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBM,1 );
            }
#endif

        } else {
            * ( _Scount[sector]+_subsecAm+_subsecBm ) +=  1;
            * ( _SubSec[sector]+ ( * ( _Scount[sector]+_subsecAm+_subsecBm ) ) *9+_subsecAm+_subsecBm ) = iCell;
#ifdef VERIFY_INTERP
            if ( sector==0 ) {
                _sec0->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==1 ) {
                _sec1->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==2 ) {
                _sec2->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            } else if ( sector==3 ) {
                _sec3->setIJ ( iCell,_subsecAm+_subsecBm,1 );
            }
#endif
        }
    }
    return;
}

// ===================================================================================
void BoundInterp::GetXiEta ( int node, std::vector<double> &XiEtabound ) {
    for ( int dim=0; dim<_MeshDim; dim++ ) {
        XiEtabound.push_back ( _XiEta->getIJ ( node,dim ) );
    }
    return;
}
// ===================================================================================
void BoundInterp::GetXiEtaChi ( int node, std::vector<double> &XiEtaChi ) {
    for ( int dim=0; dim<_MeshDim; dim++ ) {
        XiEtaChi.push_back ( _XiEta->getIJ ( node,dim ) );
    }
    return;
}


// ===================================================================================
void BoundInterp::GetInterNodes ( int node, std::vector<int> &VectorContainingNodes ) {
    const int SourceNodes = _BoundingNodes->getNumberOfComponents();
    for ( int i = 0; i<SourceNodes; i++ ) {
        VectorContainingNodes.push_back ( _BoundingNodes->getIJ ( node,i ) );
    }
    return;
}


// ===================================================================================
int BoundInterp::GetBoundNodesPerCell() {
    return __BoundNodesPerCell;
}


// ===================================================================================
int BoundInterp::GetDomainType() {
    return __Domain;
}

// ===============================================================================================
void BoundInterp::XiEtaCalc_2D (
    double NodePos[],                // physical coord (in)
    std::vector<double> &XiEtaBound,  // reference coord (out)
    const int  Quad4[],                   // point coordinates for quad4 (x0,y0,x1,y1...)
    int  npt_el                       // number of points for elment

) {// =============================================================================================

    // geometry
    double  xm[DIMENSION];
    double normal[DIMENSION];
    double Co[DIMENSION*DIMENSION];
    double Coord[NDOF_FEMB*DIMENSION];
    for ( int n_dim =0; n_dim<DIMENSION; n_dim++ ) {
        xm[n_dim]=0.;
        for ( int nnode =0; nnode<npt_el; nnode++ ) {
            Coord[ nnode+n_dim*NDOF_FEMB] =  _CoordMatrix[n_dim+nnode*DIMENSION];
            xm[n_dim] += Coord[ nnode+n_dim*NDOF_FEMB]/npt_el;
        }
    }
    _fe[1]->normal_g ( Coord,xm,normal ); // normal in the new med surface

    // equation surface  coefficients
    //  0=XYZm[0] + Co[1]*psi+Co[2]*eta+Co[0]*psi*eta
    //  0=XYZm[1] + Co[4]*psi+Co[5]*eta+Co[3]*psi*eta
    //  0=XYZm[2] + Co[8]*psi+Co[7]*eta+Co[6]*psi*eta
    double XYZm[DIMENSION];
    double XYZs[DIMENSION];
    double xieta[_MeshDim]; // reference point coordinates
    for ( int kdim=0; kdim<_SpaceDim; kdim++ ) {
        XYZm[kdim]=0.;
        Co[kdim*_SpaceDim   ] =  1.e-20;
        Co[kdim*_SpaceDim +1] =  1.e-20;
        Co[kdim*_SpaceDim +2] =  1.e-20;
        for ( int node=0; node< npt_el; node++ ) {
            Co[kdim*_SpaceDim   ] += 0.25*_CoordMatrix[_SpaceDim*node + kdim]*Quad4[node]*Quad4[node+ npt_el];
            Co[kdim*_SpaceDim +1] += 0.25*_CoordMatrix[_SpaceDim*node + kdim]*Quad4[node];
            Co[kdim*_SpaceDim +2] += 0.25*_CoordMatrix[_SpaceDim*node + kdim]*Quad4[node+ npt_el];
            XYZm[kdim] += 0.25*_CoordMatrix[_SpaceDim*node + kdim];
        }
        XYZs[kdim] = NodePos[kdim] - XYZm[kdim];
    }

    //
    int normal_max_dir= ( ( fabs ( normal[0] ) >fabs ( normal[1] ) ) ?   0:1 ) ;
    normal_max_dir= ( fabs ( normal[normal_max_dir] ) >fabs ( normal[2] ) ? normal_max_dir : 2 );
    int id1= ( normal_max_dir+1 ) %3; // non normal direction 1
    int id2= ( normal_max_dir+2 ) %3; // non normal direction 2

    // if  Co[1+3*id1] and Co[0+3*id1] are null -> inversion coeff
    if ( fabs ( Co[1+3*id1] ) + fabs ( Co[0+3*id1] ) <1.e-10 )   {
        int itmp= id1;
        id1=id2;
        id2=itmp;
        std::cout << " index inversion "<< std::endl;
    }

    double alpha_0= Co[0+3*id1]*Co[2+3*id2]-Co[0+3*id2]*Co[2+3*id1];
    double beta_0= -XYZs[id2]*Co[0+3*id1]+XYZs[id1]*Co[0+3*id2]+Co[1+3*id1]*Co[2+3*id2]-Co[1+3*id2]*Co[2+3*id1];
    double gamma_0= -Co[1+3*id1]*XYZs[id2]+Co[1+3*id2]*XYZs[id1];


    if ( fabs ( alpha_0 ) >1.e-10 ) { //  alpha_0 != 0 two roots 1.1e-10 ==========================================
        double delta_0=beta_0*beta_0-4*alpha_0*gamma_0;
        if ( delta_0<1.e-10 ) {
            std::cout << " !!!!!!!!!!!!!!! Error delta "<< std::endl;
        }
        // first root (x_00, y_00)
        xieta[1]= ( -beta_0+sqrt ( delta_0 ) ) / ( 2.*alpha_0 );
        xieta[0]=- ( -XYZs[id1]+Co[2+3*id1]*xieta[1] ) / ( Co[1+3*id1]+Co[0+3*id1]*xieta[1] );
        if ( fabs ( Co[1+3*id1]+Co[0+3*id1]*xieta[1] ) <1e-10 ) {
            xieta[0]=- ( -XYZs[id2]+Co[2+3*id2]*xieta[1] ) / ( Co[1+3*id2]+Co[0+3*id2]*xieta[1] );
        }
        // second root (x_01, y_01)
        if ( fabs ( xieta[0] ) > 1.001 ||  fabs ( xieta[1] ) > 1.001 ) { //  std::cout<<" first not good x1 "<< xieta[0] << " " << xieta[1] <<std::endl;
            xieta[1]=0.5* ( -beta_0-sqrt ( delta_0 ) ) /alpha_0;
            xieta[0]=- ( -XYZs[id1]+Co[2+3*id1]*xieta[1] ) / ( Co[1+3*id1]+Co[0+3*id1]*xieta[1] );
            if ( fabs ( Co[1+3*id1]+Co[0+3*id1]*xieta[1] ) <1e-10 ) {
                xieta[0]=- ( -XYZs[id2]+Co[2+3*id2]*xieta[1] ) / ( Co[1+3*id2]+Co[0+3*id2]*xieta[1] );
            }
        }
    } else { // alpha_0 = 0 only oen root =====================================================
        xieta[1]=-gamma_0/beta_0;//        std::cout << " one root only  " << std::endl;
        xieta[0]=- ( -XYZs[id1]+Co[2+3*id1]*xieta[1] ) / ( Co[1+3*id1]+Co[0+3*id1]*xieta[1] );
    }
    for ( int dim=0; dim<_MeshDim; dim++ ) {
        XiEtaBound.push_back ( xieta[dim] );
    }

//       if( delta_0<1.e-20){ std::cout << " !!!!!!!!!!!!!!!  error delta "<< std::endl;}
//       xieta[1]=(-beta_0+sqrt(delta_0))/(2*alpha_0);
//      if( fabs( xieta[1])< 1.001){
//          xieta[0]=-(-XYZs[id1]+Co[2+3*id1]*xieta[1])/(Co[1+3*id1]+Co[0+3*id1]*xieta[1]);
//            std::cout<<" Fisrt x1 "<< xieta[0] << " " << xieta[1] <<std::endl;
//         if( fabs( xieta[0])< 1.001){for(int dim=0; dim<_MeshDim; dim++) XiEtaBound.push_back(xieta[dim]);   }
//      }
//        xieta[1]=-(beta_0+sqrt(delta_0))/(2*alpha_0);
//       if( fabs( xieta[1])< 1.001){
//          xieta[0]=-(-XYZs[id1]+Co[2+3*id1]*xieta[1])/(Co[1+3*id1]+Co[0+3*id1]*xieta[1]);
//            std::cout<<" second x1 "<< xieta[0] << " " << xieta[1] <<std::endl;
//         if( fabs( xieta[0])< 1.001){for(int dim=0; dim<_MeshDim; dim++) XiEtaBound.push_back(xieta[dim]);   }
//      }

//     double A = XYZs[0]-XYZs[1]*Co[0]/(Co[3]);
//     double D = XYZs[0]-XYZs[2]*Co[0]/(Co[6]);
//     double B = Co[1]-Co[4]*Co[0]/(Co[3]);
//     double C = Co[2]-Co[5]*Co[0]/(Co[3]);
//     double E = Co[1]-Co[7]*Co[0]/(Co[6]);
//     double F = Co[2]-Co[8]*Co[0]/(Co[6]);
//
// //     xieta[0]  = (A*F-C*D)/(B*F-C*E);    xieta[1] = (B*D-E*A)/(B*F-C*E);
//     double xieta0  = (A*F-C*D)/(B*F-C*E);  double   xieta1 = (B*D-E*A)/(B*F-C*E);
//     if(fabs(B*F-C*E)<1.e-20) { std::cout<<"determinante bob "<<B*F-C*E<<std::endl; }
// //     for(int dim=0; dim<_MeshDim; dim++) { XiEtaBound.push_back(xieta[dim]); }
//     std::cout<<" bob calculated xi: "<<xieta0<<" eta: "<<xieta1<<" against xi: "<<xieta[0]<<" eta: "<<xieta[1]<<std::endl;



    return;
}



// ===============================================================================================
void BoundInterp::XiEtaCalc_1D (
    double NodePos[],                // physical coord (in)
    std::vector<double> &XiEtaBound,  // reference coord (out)
    int  npt_el
) {// =============================================================================================

    double xm = 0.5* ( _CoordMatrix[0]+_CoordMatrix[2] );
    double ym = 0.5* ( _CoordMatrix[1]+_CoordMatrix[3] );
    double a  = 0.5* ( -_CoordMatrix[0]+_CoordMatrix[2] );
    double b  = 0.5* ( -_CoordMatrix[1]+_CoordMatrix[3] );
    double xi = ( NodePos[0]+NodePos[1] - ( xm+ym ) ) / ( a+b );
    XiEtaBound.push_back ( xi );
    return;
};

int BoundInterp::IsNodeContained() {
    int Verify=1;
    double toll = 1.e-3;
    if ( _FamilyType==1 ) {
        if ( fabs ( _XiEtaChi[0] ) > 1+toll ||
                fabs ( _XiEtaChi[1] ) > 1+toll ||
                fabs ( _XiEtaChi[DIMENSION-1] ) > 1+toll ) Verify=0;
    }
    if ( _FamilyType==0 ) {
        if ( _XiEtaChi[0]<-toll || _XiEtaChi[0]>1+toll ||
                _XiEtaChi[1]<-toll || _XiEtaChi[1]>1+toll ||
                _XiEtaChi[DIMENSION-1]<-toll || _XiEtaChi[DIMENSION-1]>1+toll ) Verify=0;
    }

    return Verify;
}


void BoundInterp::XiEtaChiCalc ( double NodePos[] ) {
    double XiEtaChi[_MeshDim];
    for ( int dim=0; dim<_MeshDim; dim++ ) {
        XiEtaChi[dim] = 0.;
    }
    double deltaXEC[_MeshDim];
    double SquaredError;
    const double Toll = 1.e-12;

    bool FoundXiEtaChi = false;
    int count =0;
    SquaredError = CalcF ( NodePos, XiEtaChi );

    if ( SquaredError < Toll ) {
        FoundXiEtaChi = true;
    }

//   std::clock_t par_time1 = std::clock();

    while ( !FoundXiEtaChi ) {

        FirstDerF ( NodePos, XiEtaChi );
        SecondDerF ( NodePos, XiEtaChi );

        // Assembling the inverse of Hessian matrix
        double InvHex[_MeshDim*_MeshDim];
        double DetHex = 0.;//axx*ayy*azz + 2.*axy*ayz*axz - ayy*axz*axz - azz*axy*axy - axx*ayz*ayz;

        if ( _MeshDim==3 ) {
            DetHex = _d2F[0]*_d2F[4]*_d2F[8] + 2.*_d2F[1]*_d2F[5]*_d2F[2] - _d2F[4]*_d2F[2]*_d2F[2] - _d2F[8]*_d2F[1]*_d2F[1] - _d2F[0]*_d2F[5]*_d2F[5];
            for ( int row=0; row<_MeshDim; row++ ) {
                for ( int col=0; col<_MeshDim; col++ ) {
                    if ( row<=col ) {
                        int a = ( int ) ( 1- ( ( row+1 ) >>1 ) );
                        int b = ( int ) ( 2- ( row>>1 ) );
                        int c = ( int ) ( 1- ( ( col+1 ) >>1 ) );
                        int d = ( int ) ( 2- ( col>>1 ) );
                        int sign = ( int ) ( col+row ) %2;
                        InvHex[row*_MeshDim + col] = ( 1-2*sign ) * ( _d2F[a*_MeshDim + c]*_d2F[b*_MeshDim + d] - _d2F[a*_MeshDim + d]*_d2F[b*_MeshDim + c] ) /DetHex;
                    } else {
                        InvHex[row*_MeshDim + col] = InvHex[col*_MeshDim + row];
                    }
                }
            }
        }
        // Per il caso 2D si possono provare anche metodi di soluzione analitica delle coordinate
        if ( _MeshDim==2 ) {
            DetHex = _d2F[0]*_d2F[3] - _d2F[1]*_d2F[2];
            for ( int row=0; row<_MeshDim; row++ ) {
                for ( int col=0; col<_MeshDim; col++ ) {
                    if ( row<=col ) {
                        int a = ( int ) ( 1- ( ( row+1 ) >>1 ) );
                        int c = ( int ) ( 1- ( ( col+1 ) >>1 ) );
                        int sign = ( int ) ( col+row ) %2;
                        InvHex[row*_MeshDim + col] = ( 1-2*sign ) *_d2F[a*_MeshDim + c]/DetHex;
                    } else {
                        InvHex[row*_MeshDim + col] = InvHex[col*_MeshDim + row];
                    }
                }
            }
        }


        // Calculating the variation of xi, eta, chi
        for ( int dir = 0; dir<_MeshDim; dir ++ ) {
            double SUM = 0.;
            for ( int i =0; i<_MeshDim; i++ ) {
                SUM += InvHex[_MeshDim*dir + i]*_dF[i];
            }
            deltaXEC[dir] = -0.5*SUM;
//     if(fabs(deltaXEC[dir])>3){
//       std::cout<<"\033[1;31m Computed delta xi eta chi is greater than 2, now abort \033[0m" <<std::endl;
//       abort();
//     }
        }

        for ( int i =0; i<_MeshDim; i++ ) {
            XiEtaChi[i] += deltaXEC[i];
        }

        SquaredError = CalcF ( NodePos, XiEtaChi );
        count ++;
        if ( SquaredError < Toll ) {
            FoundXiEtaChi = true;
//     std::cout<<" Process stopped, error function: " << SquaredError << " count: "<<count<<std::endl;
        }
        if ( count>1000 ) {
            FoundXiEtaChi = true;

//       std::cout<<" Process stopped, error function: " << SquaredError << " xi: "<< XiEtaChi[0] << " eta " <<XiEtaChi[1];
//       if(_MeshDim==3) { std::cout<<" chi " <<XiEtaChi[2]; }
//        std::cout<< std::endl;
        }
        _dF.clear();
        _d2F.clear();
    }

    for ( int Mcomp=0; Mcomp<_MeshDim; Mcomp++ ) {
        _XiEtaChi.push_back ( XiEtaChi[Mcomp] );
    }

    return;
};

// First derivatives of linear 3D test functions
double BoundInterp::FirstDer ( double XiEtaChi[],int node,int dir ) {
//   _CoordHex27
    double value ;
    if ( _FamilyType==1 ) {
        value = _XiEtaChiVert[node+dir*_VertCoordOff];
        for ( int dim=0; dim<_MeshDim; dim++ ) {
            if ( dim == dir ) {
                value *=0.5;
            } else {
                value *= 0.5* ( 1.+_XiEtaChiVert[node + dim*_VertCoordOff]*XiEtaChi[dim] );
            }
        }
    }
    if ( _FamilyType==0 ) {
        if ( node<2 ) value = ( node==dir ) ? 1.:0.;
        else value=-1.;
    }
    return value;
}

// Second derivatives of linear 3D test functions
double BoundInterp::SecondDer ( double XiEtaChi[],int node ) {
    double value;
    // In 2D 0<= node <4 and offset = 4  ->  the eta coordinates are taken from node + 8, so node+_VertCoordOff
    value  = 0.25*_XiEtaChiVert[node]*_XiEtaChiVert[node+_VertCoordOff]; // direction 0
    return value;
}
double BoundInterp::SecondDer ( double XiEtaChi[],int node, int row, int col ) {

    double value;
    if ( _FamilyType==1 ) {
        if ( row==col ) {
            value=0.;    /// Second derivatives like d2N/dxidxi are equal to zero for linear test functions
        } else {
            value = 0.25*_XiEtaChiVert[node+row*_VertCoordOff]*_XiEtaChiVert[node+col*_VertCoordOff];;
            if ( _MeshDim==3 ) {
                int a = _MeshDim-row-col;
                value *= 0.5* ( 1.+_XiEtaChiVert[node + a*_VertCoordOff]*XiEtaChi[a] );
            }
        }
    }
    if ( _FamilyType==0 ) {
        value=0;
    }
    return value;
}


double BoundInterp::LinPhi ( int node, double XiEtaChi[] ) {
    double N = 1.;
    if ( _FamilyType==1 ) {
        for ( int dir=0; dir<_MeshDim; dir++ ) {
            N *= 0.5* ( 1.+XiEtaChi[dir]*_XiEtaChiVert[node+dir*_VertCoordOff] );
        }
    }
    if ( _FamilyType==0 ) {
        if ( node==0 ) N = XiEtaChi[0];
        if ( node==1 ) N = XiEtaChi[1];
        if ( node==2 ) N = 1-XiEtaChi[0]-XiEtaChi[1];
    }
    return N;
}


// Function for the calculation of F, first derivatives of F and second derivatives of F where F is the quantity
//                 F = (xp - xi)*(xp - xi) + (yp - yi)*(yp - yi) + (zp - zi)*(zp - zi)
// (xp,yp,zp) are the coordinates of the point of interest P and (xi,yi,zi) are the one calculated with (xi,eta,chi)
// Refer to "Exact and efficient interpolation using finite elements shape functions" (here the formula is correct)
// ************************ F ******************************
double BoundInterp::CalcF ( double NodePos[], double XiEtaChi[] ) {
    double f=0.;

    for ( int i = 0; i<_SpaceDim; i++ ) {
        double SUM = 0.;
        for ( int j=0; j<_SrcCoordInterpNodes; j++ ) {
            SUM += _CoordMatrix[j*_SpaceDim + i]*LinPhi ( j, XiEtaChi );
        }
        f += ( NodePos[i] - SUM ) * ( NodePos[i] - SUM );
    }
    return f;
}

// *********************** dF ******************************
void BoundInterp::FirstDerF ( double NodePos[], double XiEtaChi[] ) {
    // Function for the calculation of first order derivatives of squared error function
    double Error[_SpaceDim];

    for ( int extdir=0; extdir<_SpaceDim; extdir++ ) {
        Error[extdir] = 0.;
        double xyzInterp = 0.;
        for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
            xyzInterp += _CoordMatrix[nNode*_SpaceDim + extdir]*LinPhi ( nNode, XiEtaChi );
        }
        Error[extdir] += NodePos[extdir] - xyzInterp;
    }

    for ( int comp = 0; comp<_SpaceDim; comp++ ) {
        double derValue = 0.;
        for ( int dir=0; dir<_SpaceDim; dir++ ) {
            double dxyzInterp = 0.;
            for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
                dxyzInterp += _CoordMatrix[nNode*_SpaceDim + dir]*FirstDer ( XiEtaChi,nNode,comp );
            }
            derValue += -2.*Error[dir]*dxyzInterp;
        }
        _dF.push_back ( derValue );
    }

    return;
};

// *********************** d2F *****************************
void BoundInterp::SecondDerF ( double NodePos[], double XiEtaChi[] ) {

    // Function for the calculation of second order derivatives of squared error function
    // Vector _d2F contains the Hessian matrix of  F
    double Error[_SpaceDim];

    for ( int extdir=0; extdir<_SpaceDim; extdir++ ) {
        Error[extdir] = 0.;
        double xyzInterp = 0.;
        for ( int nNode = 0; nNode < _SrcCoordInterpNodes; nNode++ ) {
            xyzInterp += _CoordMatrix[nNode*_SpaceDim + extdir]*LinPhi ( nNode, XiEtaChi );
        }
        Error[extdir] += NodePos[extdir] - xyzInterp;
    }

    for ( int row=0; row<_MeshDim; row++ ) {
        for ( int col=0; col<_MeshDim; col++ ) {
            double DerValue=0.;
            if ( row<=col ) {
                for ( int dir=0; dir<_MeshDim; dir++ ) {
                    double d2xyzInterp = 0.;
                    double Der1 = 0.;
                    double Der2 = 0.;
                    for ( int nNode=0; nNode<_SrcCoordInterpNodes; nNode++ ) {
                        d2xyzInterp += _CoordMatrix[nNode*_SpaceDim + dir]*SecondDer ( XiEtaChi,nNode,row,col );
                        Der1 += _CoordMatrix[nNode*_SpaceDim + dir]*FirstDer ( XiEtaChi,nNode,row );
                        Der2 += _CoordMatrix[nNode*_SpaceDim + dir]*FirstDer ( XiEtaChi,nNode,col );
                    }
                    if ( col==row ) {
                        DerValue += 2.*Der1*Der2;
                    } else {
                        DerValue += -2.* ( Error[dir]*d2xyzInterp - Der1*Der2 );
                    }
                }
            } else {
                DerValue = _d2F[col*_MeshDim+row];
            }
            _d2F.push_back ( DerValue );
        }
    }

    return;
}

MEDCoupling::MEDCouplingFieldDouble *
BoundInterp::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble* SourceField, const int order ) {

    const int NComp = SourceField->getNumberOfComponents();
    int order1 = order;

    MEDCoupling::DataArrayDouble *targetArray = MEDCoupling::DataArrayDouble::New();
    targetArray -> alloc ( _TrgNodes,NComp );
    std::string EqName = SourceField->getName();
    switch ( __Domain ) {
    case ( Boundary ) :
        // Loop over the target mesh nodes
        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            GetInterNodes ( iNode,  BoundingNodes );
            std::vector<double> XiEtaBound;
            GetXiEta ( iNode, XiEtaBound );
            double CanPos[_MeshDim];
            for ( int dim=0; dim<_MeshDim; dim++ ) {
                CanPos[dim]=XiEtaBound[dim];
            }
            //qui
            std::map<int,int> MedLibmesh;
            BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
            double TrgValue=0.;
            for ( int iComp=0; iComp<NComp; iComp++ ) {
                TrgValue=0.;
//         TrgValue = 0.;
//         for(int phin = 0; phin< __BoundNodesPerCell; phin ++) {
//           double val = SourceField->getIJ(BoundingNodes[phin], iComp);
//           double phi = QuadPhi(phin, CanPos);
//           TrgValue +=  val*phi;
//         }
                if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< pow ( 2,_MeshDim ); phin ++ ) { /// only edge quad and hex
                        const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        const double phi = LinPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                } else {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                        double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        double phi = QuadPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                }
                targetArray->setIJ ( iNode,iComp,TrgValue );
            }
            BoundingNodes.clear();
            XiEtaBound.clear();
        }
        break;

    case ( Volume ) :

        std::map<int,int> InterpCoordNodes;
        InterpCoordNodes[27]=8;
        InterpCoordNodes[9] =4;
        InterpCoordNodes[3] =2;
        InterpCoordNodes[10]=4;
        InterpCoordNodes[6]=3;
        InterpCoordNodes[7]=3;
        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            std::vector<double> XiEtaChi;
            double TrgValue;
            double CanPos[_MeshDim];
            GetInterNodes ( iNode,  BoundingNodes );
            GetXiEtaChi ( iNode, XiEtaChi );
            if ( BoundingNodes[0] == -1 ) {
                for ( int icomp=0; icomp<NComp; icomp++ ) {
                    targetArray->setIJ ( iNode,icomp,0. );    //(node id, comp, value)
                }
                BoundingNodes.clear(); //       break;
            } else {
                for ( int dim=0; dim<_MeshDim; dim++ ) {
                    CanPos[dim] = XiEtaChi[dim];
                }
                double val,phi;
                for ( int iComp=0; iComp<NComp; iComp++ ) {


                    if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< pow ( 2,DIMENSION ); phin ++ ) { /// only edge quad and hex
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = LinPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    } else {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = QuadPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    }
                    targetArray->setIJ ( iNode,iComp,TrgValue );
                    BoundingNodes.clear();
                    XiEtaChi.clear();
                }
            }

        }
        break;

    } // end switch
    MEDCoupling::MEDCouplingFieldDouble * f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
    f->setMesh ( __OutMesh );
    f->setArray ( targetArray );
    f->setName ( EqName );
    targetArray->decrRef();
    return f;
}

//Modifica a InterpolatedField

MEDCoupling::MEDCouplingFieldDouble *
BoundInterp::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble* SourceField, const MEDCoupling::MEDCouplingFieldDouble* TargetField, const int order ) {

    const int NComp = SourceField->getNumberOfComponents();

// const DataArrayDouble *TargetFieldArray = TargetField->getArray();

    int order1 = order;
    MEDCoupling::DataArrayDouble *targetArray = MEDCoupling::DataArrayDouble::New();
    targetArray -> alloc ( _TrgNodes,NComp );
    std::string EqName = SourceField->getName();

    switch ( __Domain ) {
    case ( Boundary ) :
        // Loop over the target mesh nodes
        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            GetInterNodes ( iNode,  BoundingNodes );
            std::vector<double> XiEtaBound;
            GetXiEta ( iNode, XiEtaBound );

            if ( BoundingNodes[0] == -1 ) {
                for ( int icomp=0; icomp<NComp; icomp++ ) {
                    double TargetFieldValue = TargetField->getIJ ( iNode,icomp );
                    targetArray->setIJ ( iNode,icomp,TargetFieldValue );
                }
                BoundingNodes.clear(); //       break;
            } else{
            
            double CanPos[_MeshDim];
            for ( int dim=0; dim<_MeshDim; dim++ ) {
                CanPos[dim]=XiEtaBound[dim];
            }

            //qui
            std::map<int,int> MedLibmesh;
            BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
            double TrgValue=0.;
            for ( int iComp=0; iComp<NComp; iComp++ ) {
                TrgValue=0.;
//         TrgValue = 0.;
//         for(int phin = 0; phin< __BoundNodesPerCell; phin ++) {
//           double val = SourceField->getIJ(BoundingNodes[phin], iComp);
//           double phi = QuadPhi(phin, CanPos);
//           TrgValue +=  val*phi;
//         }
                if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< pow ( 2,_MeshDim ); phin ++ ) { /// only edge quad and hex
                        const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        const double phi = LinPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                } else {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                        double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        double phi = QuadPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                }
                targetArray->setIJ ( iNode,iComp,TrgValue );
            }
            BoundingNodes.clear();
            XiEtaBound.clear();
            }
        }
        break;

    case ( Volume ) :

        std::map<int,int> InterpCoordNodes;
        InterpCoordNodes[27]=8;
        InterpCoordNodes[9] =4;
        InterpCoordNodes[3] =2;
        InterpCoordNodes[10]=4;
        InterpCoordNodes[6]=3;
        InterpCoordNodes[7]=3;

        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            std::vector<double> XiEtaChi;
            double TrgValue;
            double CanPos[_MeshDim];
            GetInterNodes ( iNode,  BoundingNodes );
            GetXiEtaChi ( iNode, XiEtaChi );

            if ( BoundingNodes[0] == -1 ) {
                for ( int icomp=0; icomp<NComp; icomp++ ) {
                    double TargetFieldValue = TargetField->getIJ ( iNode,icomp );
                    targetArray->setIJ ( iNode,icomp,TargetFieldValue );
                }
                BoundingNodes.clear(); //       break;
            } else {
                for ( int dim=0; dim<_MeshDim; dim++ ) {
                    CanPos[dim] = XiEtaChi[dim];
                }
                double val,phi;
                for ( int iComp=0; iComp<NComp; iComp++ ) {


                    if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< pow ( 2,DIMENSION ); phin ++ ) { /// only edge quad and hex
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = LinPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    } else {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = QuadPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    }
                    targetArray->setIJ ( iNode,iComp,TrgValue );
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
BoundInterp::InterpolatedField ( const MEDCoupling::MEDCouplingFieldDouble* SourceField, double DefaultValue, const int order ) {

    const int NComp = SourceField->getNumberOfComponents();

    int order1 = order;
    MEDCoupling::DataArrayDouble *targetArray = MEDCoupling::DataArrayDouble::New();
    targetArray -> alloc ( _TrgNodes,NComp );
    std::string EqName = SourceField->getName();

    switch ( __Domain ) {
    case ( Boundary ) :
        // Loop over the target mesh nodes
        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            GetInterNodes ( iNode,  BoundingNodes );
            std::vector<double> XiEtaBound;
            GetXiEta ( iNode, XiEtaBound );

            if ( BoundingNodes[0] == -1 ) {
                for ( int icomp=0; icomp<NComp; icomp++ ) {
                    targetArray->setIJ ( iNode,icomp,DefaultValue );
                }
                BoundingNodes.clear(); //       break;
            } else{
            double CanPos[_MeshDim];
            for ( int dim=0; dim<_MeshDim; dim++ ) {
                CanPos[dim]=XiEtaBound[dim];
            }

            //qui
            std::map<int,int> MedLibmesh;
            BuildCanonicalElementNodesMap ( _SrcCellNodes, MedLibmesh );
            double TrgValue=0.;
            for ( int iComp=0; iComp<NComp; iComp++ ) {
                TrgValue=0.;
//         TrgValue = 0.;
//         for(int phin = 0; phin< __BoundNodesPerCell; phin ++) {
//           double val = SourceField->getIJ(BoundingNodes[phin], iComp);
//           double phi = QuadPhi(phin, CanPos);
//           TrgValue +=  val*phi;
//         }
                if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< pow ( 2,_MeshDim ); phin ++ ) { /// only edge quad and hex
                        const double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        const double phi = LinPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                } else {
//             TrgValue = 0.;
                    for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                        double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                        double phi = QuadPhi ( phin, CanPos );
                        TrgValue +=  val*phi;
                    }
                }



                targetArray->setIJ ( iNode,iComp,TrgValue );
            }
            BoundingNodes.clear();
            XiEtaBound.clear();
            }
        }
        break;

    case ( Volume ) :

        std::map<int,int> InterpCoordNodes;
        InterpCoordNodes[27]=8;
        InterpCoordNodes[9] =4;
        InterpCoordNodes[3] =2;
        InterpCoordNodes[10]=4;
        InterpCoordNodes[6]=3;
        InterpCoordNodes[7]=3;
        for ( int iNode = 0; iNode<_TrgNodes; iNode ++ ) {
            // Vector where we store the ids of the source mesh nodes that are used for the interpolation of the solution on target mesh node iNode
            std::vector<int> BoundingNodes;
            std::vector<double> XiEtaChi;
            double TrgValue;
            double CanPos[_MeshDim];
            GetInterNodes ( iNode,  BoundingNodes );
            GetXiEtaChi ( iNode, XiEtaChi );

            if ( BoundingNodes[0] == -1 ) {
                for ( int icomp=0; icomp<NComp; icomp++ ) {
                    targetArray->setIJ ( iNode,icomp,DefaultValue );
                }
                BoundingNodes.clear(); //       break;
            } else {
                for ( int dim=0; dim<_MeshDim; dim++ ) {
                    CanPos[dim] = XiEtaChi[dim];
                }
                double val,phi;
                for ( int iComp=0; iComp<NComp; iComp++ ) {


                    if ( ( ( EqName=="NS0" || EqName=="FSI0" ||EqName=="FSIA0"||EqName=="NSA0" ) && iComp == _SpaceDim ) ||order==1 ) {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< pow ( 2,DIMENSION ); phin ++ ) { /// only edge quad and hex
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = LinPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    } else {
                        TrgValue = 0.;
                        for ( int phin = 0; phin< __BoundNodesPerCell; phin ++ ) {
                            double val = SourceField->getIJ ( BoundingNodes[phin], iComp );
                            double phi = QuadPhi ( phin, CanPos );
                            TrgValue +=  val*phi;
                        }
                    }
                    targetArray->setIJ ( iNode,iComp,TrgValue );
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




double BoundInterp::QuadPhi ( int N, double GPoint[] ) {
    double value;

    // General formula for test function: N_i = (1-0.5*|xi_i|)*((2|xi_i|-1)*xi*xi + xi_i*xi + (1-|xi_i|)) for xi, eta, chi

    if ( _FamilyType==1 ) {
        switch ( _MeshDim ) {
        case ( 1 ) :
            value = ( 1.-0.5*fabs ( _CoordEdge3[N] ) ) *
                    ( ( 2.*fabs ( _CoordEdge3[N] )-1. ) *GPoint[0]*GPoint[0] +_CoordEdge3[N]*GPoint[0] + ( 1.-fabs ( _CoordEdge3[N] ) ) );
            break;
        case ( 2 ) :
            value=1.;
            for ( int dim=0; dim<_MeshDim; dim++ ) value *= ( 1.-0.5*fabs ( _CoordQuad9[N + dim*_Quad9Off] ) ) *
                        ( ( 2.*fabs ( _CoordQuad9[N+ dim*_Quad9Off] )-1. ) *GPoint[dim]*GPoint[dim]   +_CoordQuad9[N+ dim*_Quad9Off]*GPoint[dim] + ( 1.-fabs ( _CoordQuad9[N+ dim*_Quad9Off] ) ) );
            break;
        case ( 3 ) :
            value=1.;
            for ( int dim=0; dim<_MeshDim; dim++ ) value *= ( 1.-0.5*fabs ( _CoordHex27[N + dim*_Hex27Off] ) ) *
                        ( ( 2.*fabs ( _CoordHex27[N+ dim*_Hex27Off] )-1. ) *GPoint[dim]*GPoint[dim]   +_CoordHex27[N+ dim*_Hex27Off]*GPoint[dim] + ( 1.-fabs ( _CoordHex27[N+ dim*_Hex27Off] ) ) );
            break;
        }
    }
    if ( _FamilyType==0 ) {
        double lambda;
        lambda=1- ( GPoint[0]+GPoint[1] );
        if ( N==0 ) value= ( 2*GPoint[0]-1 ) *GPoint[0];
        if ( N==1 ) value= ( 2*GPoint[1]-1 ) *GPoint[1];
        if ( N==2 ) value= ( 2*lambda-1 ) *lambda;
        if ( N==3 ) value=4*GPoint[0]*GPoint[1];
        if ( N==4 ) value=4*GPoint[1]*lambda;
        if ( N==5 ) value=4*GPoint[0]*lambda;
    }
    return value;
};
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 


