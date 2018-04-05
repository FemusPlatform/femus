#include "MMed.h"

#ifdef HAVE_MED

#include <set>
#include <cstring>
#include <ctime>


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
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
#include "numeric_vectorM.h"

// // #define USE_FEMUS (1)   uncomment to use

#ifdef USE_FEMUS==1
#include "InterfaceFunctionM.h"
#include "MGFEMap.h"
#include "EquationSystemsExtendedM.h"
#endif


#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDLoader.hxx"
#include "BoundingBox.hxx"
#include "InterpKernelGeo2DNode.hxx"

using namespace ParaMEDMEM;
#endif

#ifdef HAVE_MED
MMed::MMed() {}

MMed::MMed(
  const MEDCouplingUMesh * SourceMesh,
  const MEDCouplingUMesh * TargetMesh,
  int DomainType
) {
}

MMed::~MMed() {
  delete _fe[2];
  delete _fe[1];
  delete _fe[0];
}

//**********************************************************
// Print a vector containing med fields in a unique med file
void MMed::PrintMed(std::vector<ParaMEDMEM::MEDCouplingFieldDouble *>f, std::string FileName) {
  for(int j=0; j<f.size(); j++) { PrintMed(f[j],FileName,j); }
  return;
}
//**********************************************************
/// This function prints a single med field ParaMEDMEM::MEDCouplingFieldDouble *f  in a med file std::string FileName
/// on a new file if  int n !=0
void MMed::PrintMed(ParaMEDMEM::MEDCouplingFieldDouble *f, std::string FileName, int n) {
  bool WriteFromScratch = (n>0.5) ?false:true;
  MEDLoader::WriteField("RESU_MED/"+FileName+".med",f,WriteFromScratch);
  return;
}

#ifdef USE_FEMUS==1
//**********************************************************
// For each FEMus problem stored in the vector it prints
// a med file containing all the fields of the FEMus problem
void MMed::PrintMed(std::vector<FEMUS *>PFem) {
  for(int j=0; j<PFem.size(); j++) { PrintMed(PFem[j],j); }
  return;
}

//**********************************************************
/// This function prints all the fields of a FEMus problem in a unique med file
void MMed::PrintMed(FEMUS *PFem,int nDom) {
  std::vector<std::string> FieldsNames;
  PFem->get_MGExtSystem().getEqsNames(FieldsNames);
  int NumOfComp = 1;
  int FirstComp = 0;
  int Dim = PFem->get_MGMesh()._dim;
  string MeshDir = PFem->GetMeshDir();
  const std::string QuadFieldsFileName = "RESU_MED/SourceField" + std::to_string(nDom) + ".med";

  CreateInterfaces(PFem);

  if(PFem->get_proc() ==0) { MEDLoader::WriteUMeshDep(QuadFieldsFileName,PFem->getUMesh(22),true); }

  for(int i=0; i<FieldsNames.size(); i++) {
    if(FieldsNames[i] != "NSP" && FieldsNames[i] != "NS2P") {
      if(FieldsNames[i]=="NS0" || FieldsNames[i]=="FSI0" || FieldsNames[i]=="FSIA0"||FieldsNames[i]=="NSA0") {NumOfComp = Dim + 1;}

      const MEDCouplingFieldDouble *TrgField=PFem->getValuesOnBoundary(22,FieldsNames[i],NumOfComp,FirstComp);
      TrgField->checkCoherency();

      if(PFem->get_proc() ==0) { MEDLoader::WriteFieldUsingAlreadyWrittenMesh(QuadFieldsFileName,TrgField); }

      std::cout<<"\033[1;36m Field: "<<FieldsNames[i]<<" copied in ";
      std::cout<<"\t in "<<QuadFieldsFileName<<"\033[0m\n";
      std::cout<<"----------------------------------------------------------------\n";
      NumOfComp =1;
    }
    if(FieldsNames[i] == "NSP" || FieldsNames[i] == "NS2P") {
      const MEDCouplingFieldDouble *TrgField= PFem->getValuesOnBoundary(11,FieldsNames[i],NumOfComp ,FirstComp);
      TrgField->checkCoherency();
      if(PFem->get_proc() ==0) { MEDLoader::WriteFieldUsingAlreadyWrittenMesh(QuadFieldsFileName,TrgField); }
      std::cout<<"\033[1;36m Field: "<<FieldsNames[i]<<" copied in ";
      std::cout<<"\t in "<<QuadFieldsFileName<<"\033[0m\n";
      std::cout<<"----------------------------------------------------------------\n";
    }
  }
  FieldsNames.clear();

  return;
}

//**********************************************************
// This function creates two volume interfaces for the PFem FEMUS problem,
// in order to print both quadratic and linear fields
void MMed::CreateInterfaces(
  FEMUS* PFem   // FEMUs problem
) {
  string MeshDir =  PFem->get_MGExtSystem().getMeshMG()->_mgutils._app_dir +  "/MESH/";
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(MeshDir + PFem->getMeshName());
  std::vector<std::string> GroupNames = MEDLoader::GetMeshGroupsNames(MeshDir + PFem->getMeshName(), meshNames[0]);

  int posP            = PFem->getMeshName().find(".");
  std::string MeshFileName_Gen = PFem->getMeshName().substr(0,posP) + "_MedToMg.med" ;
  int nGroups0 = GroupNames.size();

//     if ( nGroups0>0 )
//         {
//         std::vector<int> NumericGroupNames;
//         for ( int i=0; i< nGroups0; i++ ) NumericGroupNames.push_back ( std::stoi ( GroupNames[i] ) );
//         sort ( NumericGroupNames.begin(), NumericGroupNames.end() );
//         bool VolumeGroup = ( (NumericGroupNames[0]<10&&NumericGroupNames[0]!=0) ? true:false );
//
//         if ( !VolumeGroup )
//             {
//             std::cout<<"\033[1;36m Creating interfaces for mesh without volume groups \033[0m\n";
//             PFem->init_interface ( 22,2,MeshFileName_Gen );
//             PFem->init_interface ( 11,1,MeshFileName_Gen );
//             }
//
//         if ( VolumeGroup )
//             {
//             PFem->init_interface ( 22,NumericGroupNames[0],2,MeshFileName_Gen );
//             PFem->init_interface ( 11,NumericGroupNames[0],1,MeshFileName_Gen );
//             }
//         }
//     else
//         {
  std::cout<<"\033[1;36m Creating interfaces for mesh without groups \033[0m\n";
  PFem->init_interface(22,2,MeshFileName_Gen);
  PFem->init_interface(11,1,MeshFileName_Gen);
//         }

  meshNames.clear();
  GroupNames.clear();

  return;
}
#endif

MEDCouplingFieldDouble* MMed::InterpolatedField(const MEDCouplingFieldDouble* SourceField, int order) {
  MEDCouplingFieldDouble* f;
  return f;
}

void MMed::FillParameters(const MEDCouplingUMesh* SourceMesh, const MEDCouplingUMesh* TargetMesh, int DomainType) {
  return;
}







/// This function integrates the field const ParaMEDMEM::MEDCouplingFieldDouble* $Field
/// with method \$method (mean,Aximean,bulk, axibulk,area,l2norm)
double MMed::Integrate(
  const ParaMEDMEM::MEDCouplingFieldDouble* Field, ///< solution field      (in)
  int   order,                                     ///< interpolation order (in)
  int   n_cmp,                                     ///< Number of variables (in)
  int   first_cmp,                                 ///< First component     (in)
  int   method,                                    ///< Method  method  (mean,Aximean,bulk, axibulk,area,norma)     (in)
  const ParaMEDMEM::MEDCouplingFieldDouble* VelField ///< velocity field    (in)
) {
  std::clock_t IntegrationBegin = std::clock();

  int rad = method;

  // Axisym and AxyBulk are not valid with 3D geometries
  if(method == AxiMean && DIMENSION == 3) {
    std::cout<<"\033[1;31m +++++++++ Attention! Axisym not valid with 3D geometry! Method converted to Normal ++++++++++ \033[0m" <<endl; rad = 0;
  }
  if(method == AxiBulk && DIMENSION == 3) {
    std::cout<<"\033[1;31m +++++++++ Attention! AxiBulk not valid with 3D geometry! Method converted to Bulk ++++++++++ \033[0m" <<endl;
    rad = 2;
  }
  bool BulkMedium=_BulkMedium;
  
  IntCoefficients(rad);

  const ParaMEDMEM::MEDCouplingMesh * Mesh = Field->getMesh();
  const int offset= Mesh->getNumberOfNodes();               // offset
  const int dim   = Mesh->getSpaceDimension();              // dimension
  const int nElem_med   = Mesh ->getNumberOfCells();        // # of elements
  const int DimRelToMax = Mesh->getSpaceDimension() - Mesh->getMeshDimension(); // =0 vol ; =1 sup

  // number of nodes in the canonical element
  std::vector<int> nodes;    Mesh->getNodeIdsOfCell(0,nodes);  const int NodesPerCell = nodes.size();  nodes.clear();

  // Number of nodes for the interpolation of the field on the gauss point
  int FieldNodesPerCell;   FieldNodes(NodesPerCell, FieldNodesPerCell, order);
//
  const int FOrd = order;           // Order of the integrated field
  const int VOrd = 2;                // Order of the velocity field
  const int XOrd = 2;                // Order of the coordinates
  const int el_nnodes   = NDOF_FEMB;

  std::vector< double > coord, NodeVar, NodeVel, PointsCoords;
  NodeVar.resize(FieldNodesPerCell);  NodeVel.resize(NodesPerCell*dim);  PointsCoords.resize(NodesPerCell*dim);
  // FEMUS -> Med canonical element connectivity map
  std::map<int,int> MedLibmesh;  BuildCanonicalElementNodesMap(NodesPerCell,MedLibmesh);

  std::vector< int > conn;
  _AREA = _INTEGRAL = _VELOCITY = 0.;

  for(int i_mg=0; i_mg < nElem_med; i_mg++) {                      // Loop over the boundary elements
    Mesh -> getNodeIdsOfCell(i_mg,conn);            // connectivity of cell i_mg
    for(int i_node=0; i_node < NodesPerCell; i_node++) {         // Loop over the element nodes
      Mesh -> getCoordinatesOfNode(conn[MedLibmesh[i_node]],coord);
      for(int j =0; j<dim; j++) {
        PointsCoords[i_node + j*NodesPerCell] = coord[j];
      }
      coord.clear();
    }
    for(int i_node = 0; i_node<FieldNodesPerCell; i_node ++) NodeVar[i_node]=Field->getIJ(conn[MedLibmesh[i_node]],first_cmp); 
    if(BulkMedium) {
      for(int dir=0; dir<dim; dir++)
        for(int i_node = 0; i_node<NodesPerCell; i_node ++) {
          NodeVel[i_node + dir*NodesPerCell]=Field->getIJ(conn[MedLibmesh[i_node]],dir);
        }
    }
           
    GaussLoop(NodeVar,NodeVel,BulkMedium,dim,DimRelToMax,XOrd,FOrd,VOrd,NodesPerCell,PointsCoords,FieldNodesPerCell);
    conn.clear();
  }

  double value_fin = GetIntResult(method);

  std::clock_t IntegrationEnd = std::clock();
  std::cout<<"\033[1;35m  Time for integration  "<<  double(IntegrationEnd-IntegrationBegin) /CLOCKS_PER_SEC <<"\033[0m"<< std::endl;

  return value_fin;
}






void MMed::BuildCanonicalElementNodesMap(int NodesPerCell, std::map<int,int> &MedLibmesh) {

  switch(NodesPerCell) {
  case(3):
    MedLibmesh[0] =  0;
    MedLibmesh[1] =  1;
    MedLibmesh[2] =  2;
    break;
  case(9):
    MedLibmesh[0] =  0;
    MedLibmesh[1] =  1;
    MedLibmesh[2] =  2;
    MedLibmesh[3] =  3;
    MedLibmesh[4] =  4;
    MedLibmesh[5] =  5;
    MedLibmesh[6] =  6;
    MedLibmesh[7] =  7;
    MedLibmesh[8] =  8;
    break;
  case(27):
    // map 1
    /*MedLibmesh[4 ] =  0;  */     /*  MedLibmesh[0] =  5;   */     /*MedLibmesh[0] =  0;       */        MedLibmesh[0] =  4;
    /*MedLibmesh[5 ] =  1;  */     /*  MedLibmesh[1] =  6;   */     /*MedLibmesh[1] =  1;       */        MedLibmesh[1] =  5;
    /*MedLibmesh[6 ] =  2;  */     /*  MedLibmesh[2] =  7;   */     /*MedLibmesh[2] =  2;       */        MedLibmesh[2] =  6;
    /*MedLibmesh[7 ] =  3;  */     /*  MedLibmesh[3] =  4;   */     /*MedLibmesh[3] =  3;       */        MedLibmesh[3] =  7;
    /*MedLibmesh[0 ] =  4;  */     /*  MedLibmesh[4] =  1;   */     /*MedLibmesh[4] =  4;       */        MedLibmesh[4] =  0;
    /*MedLibmesh[1 ] =  5;  */     /*  MedLibmesh[5] =  2;   */     /*MedLibmesh[5] =  5;       */        MedLibmesh[5] =  1;
    /*MedLibmesh[2 ] =  6;  */     /*  MedLibmesh[6] =  3;   */     /*MedLibmesh[6] =  6;       */        MedLibmesh[6] =  2;
    /*MedLibmesh[3 ] =  7;  */     /*  MedLibmesh[7] =  0;   */     /*MedLibmesh[7] =  7;       */        MedLibmesh[7] =  3;
    /*MedLibmesh[16] =  8;  */     /*  MedLibmesh[8] =  13;  */     /*MedLibmesh[8] =  8;       */        MedLibmesh[8] =  12;
    /*MedLibmesh[17] =  9;  */     /*  MedLibmesh[9] =  14;  */     /*MedLibmesh[9] =  9;       */        MedLibmesh[9] =  13;
    /*MedLibmesh[18] =  10; */     /*  MedLibmesh[10] =  15; */     /*MedLibmesh[10] =  10;     */        MedLibmesh[10] =  14;
    /*MedLibmesh[19] =  11; */     /*  MedLibmesh[11] =  12; */     /*MedLibmesh[11] =  11;     */        MedLibmesh[11] =  15;
    /*MedLibmesh[8 ] =  12; */     /*  MedLibmesh[12] =  17; */     /*MedLibmesh[12] =  12;     */        MedLibmesh[12] =  16;
    /*MedLibmesh[9 ] =  13; */     /*  MedLibmesh[13] =  18; */     /*MedLibmesh[13] =  13;     */        MedLibmesh[13] =  17;
    /*MedLibmesh[10] =  14; */     /*  MedLibmesh[14] =  19; */     /*MedLibmesh[14] =  14;     */        MedLibmesh[14] =  18;
    /*MedLibmesh[11] =  15; */     /*  MedLibmesh[15] =  12; */     /*MedLibmesh[15] =  15;     */        MedLibmesh[15] =  19;
    /*MedLibmesh[12] =  16; */     /*  MedLibmesh[16] =  9;  */     /*MedLibmesh[16] =  16;     */        MedLibmesh[16] =  8;
    /*MedLibmesh[13] =  17; */     /*  MedLibmesh[17] =  10; */     /*MedLibmesh[17] =  17;     */        MedLibmesh[17] =  9;
    /*MedLibmesh[14] =  18; */     /*  MedLibmesh[18] =  11; */     /*MedLibmesh[18] =  18;     */        MedLibmesh[18] =  10;
    /*MedLibmesh[15] =  19; */     /*  MedLibmesh[19] =  8;  */     /*MedLibmesh[19] =  19;     */        MedLibmesh[19] =  11;
    /*MedLibmesh[25] =  20; */     /*  MedLibmesh[20] =  25; */     /*MedLibmesh[20] =  20;     */        MedLibmesh[20] =  25;
    /*MedLibmesh[21] =  21; */     /*  MedLibmesh[21] =  22; */     /*MedLibmesh[21] =  21;     */        MedLibmesh[21] =  21;
    /*MedLibmesh[22] =  22; */     /*  MedLibmesh[22] =  23; */     /*MedLibmesh[22] =  22;     */        MedLibmesh[22] =  22;
    /*MedLibmesh[23] =  23; */     /*  MedLibmesh[23] =  24; */     /*MedLibmesh[23] =  23;     */        MedLibmesh[23] =  23;
    /*MedLibmesh[24] =  24; */     /*  MedLibmesh[24] =  21; */     /*MedLibmesh[24] =  24;     */        MedLibmesh[24] =  24;
    /*MedLibmesh[20] =  25; */     /*  MedLibmesh[25] =  20; */     /*MedLibmesh[25] =  25;     */        MedLibmesh[25] =  20;
    /*MedLibmesh[26] =  26; */     /*  MedLibmesh[26] =  26; */     /*MedLibmesh[26] =  26;     */        MedLibmesh[26] =  26;
    break;
  default:
    std::cout<<"\033[0;34m=========================== Map Not Implemented Yet ====================\033[0m\n";
    break;
  }
  return;
}

// ================================================================================================================================
void MMed::IntCoefficients(int rad) {
  _CylCoord   = (1- ((rad&4) >>2)) * (rad%2) + ((rad&8) >>3);
  _BulkMedium = (1- ((rad&4) >>2)) * ((rad&2) >>1);
  _AreaCalc   = ((rad&4) >>2) * (1- (rad%2)) + ((rad&8) >>3);
  _IntegCalc  = ((rad&4) >>2) * (rad%2);
  _L2Norm     = ((rad&16) >>4);
  return;
}
void MMed::FieldNodes(
  const int NodesPerCell,
  int &Fcc,
  int order
) {
  switch(NodesPerCell)  {
  case(3) :
    Fcc = (order==1)? 2:NodesPerCell;
    break;
  case(9) :
    Fcc = (order==1)? 4:NodesPerCell;
    break;
  case(27) :
    Fcc = (order==1)? 8:NodesPerCell;
    break;
  default:
    std::cout<<"\033[0;34m++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\033[0m\n";
    break;
  }
  return;
}
void MMed::InitFe() {
  _fe[0] = new MGFE(0,ELTYPE);
  _fe[0]->init_pie();  ///> Lagrange piecewise constant
  _fe[1] = new MGFE(1,ELTYPE);
  _fe[1]->init_lin();  ///> Lagrange piecewise linear
  _fe[2] = new MGFE(2,ELTYPE);
  _fe[2]->init_qua();  ///> Lagrange piecewise quadratic
  return;
}

double MMed::GetIntResult(int method) {
  double value_fin;
  switch(method) {
  case(Mean) :
    value_fin = (_INTEGRAL/_AREA);
    break;
  case(AxiMean) :
    value_fin = (_INTEGRAL/_AREA);
    break;
  case(Bulk) :
    value_fin = (_INTEGRAL/_VELOCITY);
    break;
  case(AxiBulk) :
    value_fin = (_INTEGRAL/_VELOCITY);
    break;
  case(Integral) :
    value_fin = _INTEGRAL;
    break;
  case(Area) :
    value_fin = _AREA;
    break;
  case(AxiArea) :
    value_fin = _AREA;
    break;
  case(NormL2) :
    value_fin = _INTEGRAL;
    break;
  }
  return value_fin;
}


#ifdef USE_FEMUS==1

/// This function integrates (in $method)  the field with components $n_cmp starting from  $first_cmp
double MMed::Integrate(FEMUS *PFemus,
                       int id,                            /**< Interface name      (in)   */
                       const char *system_name,           /**< Equation name       (in)   */
                       int         n_cmp,                 /**< Number of variables (in)   */
                       int         first_cmp,             /**< First component     (in)   */
                       int         method                 /**< Method  #method (in) */
                      ) {
  std::clock_t IntegrationBegin = std::clock();

  int rad = method;

  // Axisym and AxyBulk are not valid with 3D geometries
  if(method == AxiMean && DIMENSION == 3) {
    std::cout<<"\033[1;31m +++++++++ Attention! Axisym not valid with 3D geometry! Method converted to Normal ++++++++++ \033[0m" <<endl;
    rad = 0;

  }
  if(method == AxiBulk && DIMENSION == 3) {
    std::cout<<"\033[1;31m +++++++++ Attention! AxiBulk not valid with 3D geometry! Method converted to Bulk ++++++++++ \033[0m" <<endl;
    rad = 2;
  }

  IntCoefficients(rad);

  bool BulkMedium = _BulkMedium;

  const double pi = 2.*acos(0);    //3.1415;

//   const int Level = PFemus->get_MGMesh().;
  const int Level= PFemus->get_MGMesh()._NoLevels - 1;
  const int offset= PFemus->get_MGMesh()._NoNodes[Level];
  const int dim   = PFemus->get_MGMesh()._dim;

  MGSolBase *mgsyst=PFemus->get_MGExtSystem().get_eqs(system_name);                               // Field of the variable

  // from LibMesh function fct
  InterfaceFunctionM *fct = PFemus->get_MGExtSystem().get_interface_fun(id);
  if(fct == NULL) { return NULL; }
  const int nElem_med   = fct -> getSupport() ->getNumberOfCells();
  const int nNodes_med  = fct -> getSupport() ->getNumberOfNodes();
  int Dim = fct->getSupport()->getSpaceDimension();
  const double DimRelToMax = fct -> getSupport()->getSpaceDimension()-fct -> getSupport()->getMeshDimension();

  int *map_mg_nodes  = fct -> get_map_mg();
  int *map_med_nodes = fct -> get_map_med();
  const int order = fct->get_order();
// //   if(system_name=="NSP") order =1;

  const int NodesPerCell = fct->getSupport()->getNumberOfNodesInCell(0);
  // Number of nodes for the interpolation of the field on the gauss point
  int FieldNodesPerCell;
  FieldNodes(NodesPerCell, FieldNodesPerCell, order);

  const int FOrd = order; // Order of the integrated field
  const int VOrd = 2;     // Order of the velocity field
  const int XOrd = 2;     // Order of the coordinates

  std::map <int, int> MedMg;                                    // From med numbering to mg numbering
  for(int nMed = 0; nMed < nNodes_med; nMed++) {
    const int MedNumber = map_med_nodes[nMed];
    const int MgNumber  = map_mg_nodes[nMed];
    MedMg[MedNumber] = MgNumber;
  }

  double phi_g[3][NDOF_FEM];
  _AREA = _INTEGRAL = _VELOCITY = 0.;
  double value_fin;
  const int WeightDim = ((DimRelToMax==0) ? dim-1:dim-2);
  const int QPoints = ((DimRelToMax==0) ? _fe[XOrd]->_NoGauss1[WeightDim]:_fe[XOrd]->_NoGauss1[WeightDim]);

  std::vector< int > conn;
  std::vector< double > coord, NodeVar, NodeVel, PointsCoords;

  PointsCoords.resize(NodesPerCell*dim);
  NodeVar.resize(FieldNodesPerCell);
  NodeVel.resize(NodesPerCell*dim);

  std::map<int,int> MedLibmesh;
  BuildCanonicalElementNodesMap(NodesPerCell,MedLibmesh);

  ParaMEDMEM::MEDCouplingFieldDouble *VelField=NULL;
  // The function is called only when a velocity field is available
  if(PFemus->IsVelFieldActive() == 1 && BulkMedium) { VelField=GetVelocityField(PFemus,id,1); }


  for(int i_mg=0; i_mg < nElem_med; i_mg++) {                     // Loop over the boundary elements

    fct -> getSupport() -> getNodeIdsOfCell(i_mg,conn);           // Med numbering of element nodes

    for(int i_node=0; i_node < NodesPerCell; i_node++) {
      fct -> getSupport() -> getCoordinatesOfNode(conn[MedLibmesh[i_node]],coord);
      for(int j =0; j<dim; j++) { PointsCoords[i_node + j*NodesPerCell] = coord[j]; }
      coord.clear();
    }

    // Loop for the calculation of the velocity component normal to the boundary, bulk term
    if(BulkMedium) {
      for(int i_node = 0; i_node<NodesPerCell; i_node ++) {
        VelField->getMesh()->getCoordinatesOfNode(conn[MedLibmesh[i_node]],coord);
        for(int i=0; i<Dim; i++) {
          NodeVel[i_node + NodesPerCell*i] = VelField->getIJ(conn[MedLibmesh[i_node]],i);
        }
        coord.clear();
      }
    }

    for(int i_node = 0; i_node<FieldNodesPerCell; i_node ++) {
      const int kdof_top  = mgsyst-> _node_dof[Level][MedMg[conn[MedLibmesh[i_node]]] + first_cmp*offset];
      NodeVar[i_node]      = (* (mgsyst->x_old[Level]))(kdof_top);
    }

    GaussLoop(NodeVar,NodeVel,BulkMedium,dim,DimRelToMax,XOrd,FOrd,VOrd,NodesPerCell,PointsCoords,FieldNodesPerCell);

    conn.clear();
  }

  value_fin = GetIntResult(method);

  std::clock_t IntegrationEnd = std::clock();
  std::cout<<"\033[1;35m  Time for integration  "<<  double(IntegrationEnd-IntegrationBegin) /CLOCKS_PER_SEC <<"\033[0m"<< std::endl;
  return value_fin;
}
#endif
#endif


//====================================================================================================================
void MMed::GaussLoop(
  vector< double > NodeVar,  ///< field values on cell nodes
  vector< double > Velocity, ///< velocity values on cell nodes
  bool BulkMedium,           ///< flag for bulk integration
  const int dim,             ///< space dimension
  const int DimRelToMax,     ///< volume or surface integration
  const int XOrd,            ///< interpolation order for coordinates
  const int FOrd,            ///< interpolation order for field
  const int VOrd,            ///< interpolation order for velocity field
  const int NodesPerCell,    ///< nodes for cell
  std::vector<double> PointsCoords,   ///< cell node coordinates
  int Fcc                    ///< # of field cell nodes
) {

  const int WeightDim = ((DimRelToMax==0) ? dim-1:dim-2);
  const int QPoints = ((DimRelToMax==0) ? _fe[XOrd]->_NoGauss1[WeightDim ]:_fe[XOrd]->_NoGauss1[WeightDim ]);
  int pte_coord_size= PointsCoords.size();
  double *Coordinates=new double[pte_coord_size]; for(int i=0; i<pte_coord_size; i++) { Coordinates[i]=PointsCoords[i]; }
  double normal[DIMENSION];
  if(BulkMedium) { _fe[XOrd]->normal_g(Coordinates,normal); }
  double phi_g[3][NDOF_FEM];
  const double pi = 2.*acos(0);    //3.1415;
  double InvJac[3][DIMENSION*DIMENSION];
  for(int qp=0; qp<  QPoints; qp++) {
    double det    = (DimRelToMax==0) ? _fe[XOrd] ->Jac(qp,Coordinates,InvJac[XOrd]) :_fe[XOrd] ->JacSur(qp,Coordinates,InvJac[XOrd]);         //
    double JxW_g  = det*_fe[XOrd]->_weight1[WeightDim][qp];                     // weight
    for(int i=0; i<3; i++) { _fe[i]->get_phi_gl_g(dim - DimRelToMax,qp,phi_g[i]); }       // Test functions on Gauss points
    double CylContrib = 0.;    double VarOnGauss = 0.;    double VelOnGauss = 0.;

//     Loop for the calculation of the velocity component normal to the boundary, bulk term
    if(BulkMedium) {  // --------------------------------------------------------------------------
      for(int i_node = 0; i_node<NodesPerCell; i_node ++) {
        double norm_vel = 0.;
        if(Velocity.size() != 0) {
          for(int i=0; i<dim; i++) {
            norm_vel += Velocity[i_node + NodesPerCell*i]*normal[i];
          }
        }
        VelOnGauss += norm_vel * phi_g[VOrd][i_node];
      }
    } // ------------------------------------------------------------------------------------------

    // Loop for the calculation of the correction for axisym geometry
    if(_CylCoord)  for(int i_node = 0; i_node<NodesPerCell; i_node ++) { CylContrib += 2.*pi*(Coordinates[i_node]*phi_g[XOrd][i_node]); }

    // Loop for the calculation of the field on the gauss point
    for(int i_node = 0; i_node<Fcc; i_node ++) { VarOnGauss += NodeVar[i_node] * phi_g[FOrd][i_node] ; }
    VarOnGauss *= (VarOnGauss*_L2Norm + (1.-_L2Norm));

    // --------- Sum over the Gauss points
    _AREA      += JxW_g* (CylContrib * _CylCoord + (1.-_CylCoord)) ;
    _INTEGRAL  += JxW_g* VarOnGauss  * (CylContrib * _CylCoord + (1.-_CylCoord)) * (_BulkMedium * VelOnGauss + (1.-_BulkMedium));
    _VELOCITY  += JxW_g* VelOnGauss  * (CylContrib * _CylCoord + (1.-_CylCoord));
  }
  delete [] Coordinates;
  return;
}

#ifdef USE_FEMUS==1
// ======================================================================================
MEDCouplingFieldDouble *MMed::GetVelocityField(
  FEMUS *PFem,
  int InterfaceId,
  int IsVelCoupled
) {// ===================================================================================
  MEDCouplingFieldDouble *f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  InterfaceFunctionM *Interface = PFem->get_MGExtSystem().get_interface_fun(InterfaceId);
  f->setMesh(Interface->getSupport());

  // IsVelCoupled: 0-> no  1-> yes

  const int Dim = Interface->getSupport()->getSpaceDimension();
  const int NNodes = Interface->getSupport()->getNumberOfNodes();

  DataArrayDouble *VelArray = DataArrayDouble::New();
  DataArrayDouble *newv= DataArrayDouble::New();

  newv->alloc(NNodes,Dim);

  if(!IsVelCoupled) { //Splitted System
    std::vector<string> VelNames;
    VelNames.push_back("NS0X");
    VelNames.push_back("NS0Y");
    if(Dim==3) { VelNames.push_back("NS0Z"); }

    std::vector<const DataArrayDouble *> ArrayVec;
    for(int i=0; i<Dim; i++) {
      MEDCouplingFieldDouble *TempField = PFem->getValuesOnBoundary(InterfaceId,VelNames[i],1,0);
      DataArrayDouble *TempArray = TempField->getArray();
      ArrayVec.push_back(TempArray);
    }
    newv = VelArray->Meld(ArrayVec);
    VelNames.clear();
  }
  if(IsVelCoupled) { //Coupled System
    MEDCouplingFieldDouble *TempField = PFem->getValuesOnBoundary(InterfaceId,"NS0",Dim,0);
    newv = TempField->getArray();
  }

  f->setName("Velocity");
  f->setArray(newv);
  newv->decrRef();
  VelArray->decrRef();
  return f;
}

#endif

MEDCouplingFieldDouble* MMed::GetCellField(const MEDCouplingFieldDouble* SourceField){
  MEDCouplingFieldDouble* f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS);
  const MEDCouplingMesh *SourceMesh = SourceField->getMesh();
  f->setMesh(SourceMesh);
  
  const int MeshDim = SourceMesh->getMeshDimension();
  const int SpaceDim = SourceMesh->getSpaceDimension();
  
  std::vector<int> CellConn;
  SourceMesh->getNodeIdsOfCell(0,CellConn);
  
  const int nCells = SourceMesh->getNumberOfCells();
  const int nNodesPerCell = CellConn.size();
  CellConn.clear();
  
  std::vector<double> CellField;
  CellField.resize(nNodesPerCell);
  std::vector<double> coord;
  std::vector<double> PointsCoords; PointsCoords.resize(nNodesPerCell*SpaceDim);
  std::vector<double> NodeVar; NodeVar.resize(nNodesPerCell);
  std::vector<double> NodeVel;
 
  DataArrayDouble * CellArray = DataArrayDouble::New();
  CellArray->alloc(nCells,SourceField->getNumberOfComponents());
  
  
  for(int i_cell=0; i_cell<nCells; i_cell++){
    SourceMesh->getNodeIdsOfCell(i_cell,CellConn);
    for(int i_cnode=0; i_cnode < nNodesPerCell; i_cnode++) {         // Loop over the element nodes
      SourceMesh -> getCoordinatesOfNode(CellConn[i_cnode],coord);
      for(int j =0; j<SpaceDim; j++) {
        PointsCoords[i_cnode + j*nNodesPerCell] = coord[j];
      }
      coord.clear();
    }
    for(int i_comp=0; i_comp<SourceField->getNumberOfComponents(); i_comp++){
      for(int i_cnode = 0; i_cnode<nNodesPerCell; i_cnode ++) NodeVar[i_cnode]=SourceField->getIJ(CellConn[i_cnode],i_comp); 
//       for(int i_cnode = 0; i_cnode<nNodesPerCell; i_cnode ++) NodeVar[i_cnode]=SourceField->getIJ(CellConn[NDOF_FEM-1],i_comp);
      _INTEGRAL=0.;
      _AREA=0.;
      GaussLoop(NodeVar,NodeVel, false, SpaceDim, SpaceDim-MeshDim, 2, 2, 2,nNodesPerCell,PointsCoords,nNodesPerCell);
      CellArray->setIJ(i_cell,i_comp,_INTEGRAL/_AREA);
//       CellArray->setIJ(i_cell,i_comp,NodeVar[nNodesPerCell-1]);
    }
//    std::cout<<SourceField->getNumberOfComponents()<<std::endl;
//     CellArray->setIJ(i_cell,0,0);
    CellConn.clear();
  }
  f->setArray(CellArray);
  f->setName(SourceField->getName());
  SourceMesh->decrRef();
  CellArray->decrRef();
  return f;
}

#endif

// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
