
// This file was massively overhauled
#include "gambit_io.h"

//C++ includes
#include <fstream>
#include <set>
#include <cstring>
#include <iomanip>

// Local includes
#include "libmesh_config.h"
#include "elem.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "MGFE_conf.h"
// #include "Domain_conf.h"
// lib includes
#include "Solverlib_conf.h"  // lib configure
#include "hdf5.h"            // hdf5 libs


int msh_count=0; // needed in order to use multiple meshes

// anonymous namespace to hold local data
namespace
{
using namespace libMesh;
  
///Defines a structure to hold boundary element information.
struct boundaryElementInfo {
  std::set<unsigned int> nodes;  // nodes
  unsigned int id;               // boundary groups
};

// Defines mapping from libMesh element types to Gmsh element types.
struct elementDefinition {
  std::string label;                // label
  std::vector<unsigned int> nodes;  // nodes
  ElemType type;                    // element type
  unsigned int exptype;
  unsigned int dim;                 // dim   
  unsigned int nnodes;              // number of nodes
  unsigned int nfaces;              // number of faces
};


// maps from a libMesh element type to the proper
// gambit elementDefinition.  Placing the data structure
// here in this anonymous namespace gives us the
// benefits of a global variable without the nasty
// side-effects
//map eletypes_exp: ElemType +  elementDefinition
std::map<ElemType, elementDefinition> eletypes_exp;
//map eletypes_exp:  int +  elementDefinition
std::map<unsigned int, elementDefinition> eletypes_imp;



//  GAMBIT                LIBMESH
//   NAME   FACES         NAME  NUMB
//  ---------- ------ 1 DIM ---------
//  EDGE2                  EDGE2   1
//  EDGE3                  EDGE3   8
//  ------------------ 2 DIM ---------
//  QUAD4                  QUAD4   3
//  QUAD8                  QUAD8   19
//  QUAD9                  QUAD9   10
//  TRI3                   TRI3    2
//  TRI6                   TRI6    9 
// ------------------------------------------------------------
// helper function to initialize the eletypes map
void init_eletypes() {
  if(eletypes_exp.empty() && eletypes_imp.empty()) {
    // This should happen only once.  The first time this method
    // is called, the eletypes data struture will be empty, and
    // we will fill it.  Any subsequent calls will find an initialized
    // eletypes map and will do nothing.

    //==============================
    // setup the element definitions
    elementDefinition eledef;

    // use "swap trick" from Scott Meyer's "Effective STL" to initialize
    // eledef.nodes vector

    // POINT (only Gmsh) -------------------------------------
    {
      eledef.exptype = 15;        // element libmesh number
      eledef.dim     = 0;         // dimension
      eledef.nnodes  = 1;         // number of nodes
      eledef.nodes.clear();       // clean function
      eledef.nfaces  = 0;         // number of faces
      // class definition
      eletypes_imp[15] = eledef;  // element defined by number
    }

    // EDGE2 -> EDGE2 (1) LIBMESH ----------------------------
    {
      eledef.type    = EDGE2;      // element libmesh name
      eledef.dim     = 1;          // dimension
      eledef.nnodes  = 2;          // number of nodes
      eledef.exptype = 1;          // element libmesh number  
      eledef.nodes.clear();        // clean function
      eledef.nfaces  = 2;          // number of faces    
      // class definition (by name and number)
      eletypes_exp[EDGE2] = eledef; // def by name 
      eletypes_imp[1]     = eledef; // def by number
    }

    // EDGE3 -> EDGE3 (8) LIBMESH --------------------------
    {
      eledef.type    = EDGE3;      // element libmesh name
      eledef.dim     = 1;          // dimension
      eledef.nnodes  = 3;          // number of nodes
      eledef.exptype = 8;          // element libmesh number
      eledef.nodes.clear();        // clean function
      eledef.nfaces  = 2;          // number of faces
      // class definition (by name and number)
      eletypes_exp[EDGE3] = eledef; // def by name
      eletypes_imp[8]     = eledef; // def by number
    }

    // TRI3 Gambit  -> TRI3 (2) LIBMESH ------------------------
    {
      eledef.type    = TRI3;      // element libmesh name
      eledef.dim     = 2;         // dimension 
      eledef.nnodes  = 3;         // number of nodes
      eledef.exptype = 2;         // element libmesh number
      eledef.nfaces  = 3;         // number of faces
      // map
      const unsigned int nodes[] = {0,1,2};  
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // class definition (by name and number)
      eletypes_exp[TRI3] = eledef; // def by name
      eletypes_imp[2] = eledef;    // def by number
     
    }

    // TRI6 Gambit -> TRI6 (9) LIBMESH -----------------------
    {
      eledef.type    = TRI6;      // element libmesh name
      eledef.dim     = 2;         // dimension 
      eledef.nnodes  = 6;         // number of nodes
      eledef.exptype = 9;         // element libmesh number
      eledef.nfaces  = 3;         // number of faces
      //map
      const unsigned int nodes[] = {0,3,1,4,2,5}; 
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // class definition (by name and number)
      eletypes_exp[TRI6] = eledef; // def by name
      eletypes_imp[9]    = eledef; // def by number 
    }

    // QUAD4 Gambit  -> QUAD4 (3) LIBMESH -----------------------
    {
      eledef.type    = QUAD4;    // element libmesh name
      eledef.dim     = 2;        // dimension 
      eledef.nnodes  = 4;        // number of nodes
      eledef.exptype = 3;        // element libmesh number
      eledef.nfaces  = 4;        // number of faces
      //map
      const unsigned int nodes[] = {0,1,2,3};
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);     
      eletypes_exp[QUAD4] = eledef; // def by name
      eletypes_imp[3]     = eledef; // def by number
    }

    // QUAD8 Gambit -> QUAD8 (100) LIBMESH -----------------------
    {
      eledef.type    = QUAD8;     // element libmesh name
      eledef.dim     = 2;         // dimension 
      eledef.nnodes  = 8;         // number of nodes
      eledef.exptype = 19;       // element libmesh number
      eledef.nfaces  = 4;         // number of faces 
      // map
      const unsigned int nodes[] = {0,4,1,5,2,6,3,7};
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[QUAD8] = eledef; // def by name
      eletypes_imp[19]    = eledef; // def by number
    }

    // QUAD9 Gambit -> QUAD4 (10) LIBMESH -----------------------
    {
      eledef.type    = QUAD9;     // element libmesh name
      eledef.dim     = 2;         // dimension 
      eledef.nnodes  = 9;         // number of nodes
      eledef.exptype = 10;        // element libmesh number
      eledef.nfaces  = 4;         // number of faces       
      // map
      const unsigned int nodes[] =  {0,4,1,5,2,6,3,7,8};
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[QUAD9] = eledef; // def by name
      eletypes_imp[10]    = eledef; // def by number
    }

    // HEX8 Gambit  -> HEX8 (5) LIBMESH -----------------------
    {
      eledef.type    = HEX8;     // element libmesh name
      eledef.dim     = 3;        // dimension 
      eledef.nnodes  = 8;        // number of nodes
      eledef.exptype = 5;        // element libmesh number
      eledef.nfaces  = 6;        // number of faces       
      // eledef.nodes.clear();
      //  map
      const unsigned int nodes[] = {4,5,0,1,7,6,3,2};
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[HEX8] = eledef;// def by name
      eletypes_imp[5]    = eledef;// def by number
      
    }

    // HEX20 Gambit-> HEX20 (101) LIBMESH -----------------------
    {
      eledef.type    = HEX20;             // element libmesh name
      eledef.dim     = 3;                 // dimension 
      eledef.nnodes  = 20;                // number of nodes
      eledef.exptype = 101;               // element libmesh number
      eledef.nfaces  = 6;                 // number of faces       
      //	const unsigned int nodes[] = {0,8,1,11,9,3,10,2,12,13,15,14,4,16,5,19,17,7,18,6};
      //  map
      const unsigned int nodes[] = {4,16,5,12,13,0,8,1,19,17,11,9,7,18,6,15,14,3,10,2};
      // !! negative jacobian ???????????????????
      const unsigned int nnodes = sizeof(nodes) /sizeof(nodes[0]);
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[HEX20] = eledef;// def by name
      eletypes_imp[20]    = eledef;// def by number
 
    }

    // HEX27 Gambit -> HEX27 (12) LIBMESH -----------------------
    {
      eledef.type    = HEX27;      // element libmesh name
      eledef.dim     = 3;          // dimension 
      eledef.nnodes  = 27;         // number of nodes
      eledef.exptype = 12;         // element libmesh number
      eledef.nfaces  = 6;          // number of faces       
      //  map
      const unsigned int nodes[] = {0,8,1,11,20,9,3,10,2,12,21,13,24,26,22,15,
                                    23,14,4,16,5,19,25,17,7,18,6};
      // inverse map           
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[HEX27] = eledef;// def by name
      eletypes_imp[12]    = eledef;// def by number
    }

    // TET4 Gambit
    {
      eledef.type    = TET4;
      eledef.dim     = 3;
      eledef.nnodes  = 4;
      eledef.exptype = 4;
      //eledef.nodes.clear();
      const unsigned int nodes[] = {0,1,2,3};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[TET4] = eledef;
      eletypes_imp[4]    = eledef;
      eledef.nfaces  = 4;
    }

    // TET10 Gambit
    {
      eledef.type    = TET10;
      eledef.dim     = 3;
      eledef.nnodes  = 10;
      eledef.exptype = 11;
      const unsigned int nodes[] = {0,6,2,7,9,3,4,5,8,1};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[TET10] = eledef;
      eletypes_imp[11]    = eledef;
      eledef.nfaces  = 4;
    }

    // PRISM6 Gambit
    {
      eledef.type    = PRISM6;
      eledef.dim     = 3;
      eledef.nnodes  = 6;
      eledef.exptype = 6;
      //eledef.nodes.clear();
      const unsigned int nodes[] = {2,0,1,5,3,4};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[PRISM6] = eledef;
      eletypes_imp[6]      = eledef;
      eledef.nfaces  = 5;
    }

    // PRISM15 Gambit
    // TODO: what should be done with this on writing?
    {
      eledef.type    = PRISM15;
      eledef.dim     = 3;
      eledef.nnodes  = 15;
      eledef.exptype = 103;
      // eledef.nodes.clear();
      const unsigned int nodes[] = {4,12,3,13,14,5,10,9,11,1,6,0,7,8,2};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[PRISM15] = eledef;
      eletypes_imp[16] = eledef;
      eledef.nfaces  = 5;
    }

    // PRISM18 Gambit
    {
      eledef.type    = PRISM18;
      eledef.dim     = 3;
      eledef.nnodes  = 18;
      eledef.exptype = 13;
      //   const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,
      //                            12,14,13,15,17,16};
      const unsigned int nodes[] = {5,13,4,14,12,3,11,16,10,17,15,9,2,7,1,8,6,0};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[PRISM18] = eledef;
      eletypes_imp[13]      = eledef;
      eledef.nfaces  = 5;
    }

    // PYRAMID5
    {
      eledef.type    = PYRAMID5;
      eledef.dim     = 3;
      eledef.nnodes  = 5;
      eledef.exptype = 7;
      // eledef.nodes.clear();

      const unsigned int nodes[] = {0,1,2,3,4};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[PYRAMID5] = eledef;
      eletypes_imp[7]        = eledef;
      eledef.nfaces  = 5;
    }
    //==============================
  }
}
} // end anonymous namespace
// ====================================================
/// Print int data into dhdf5 file
hid_t GambitIO::print_Ihdf5(
  hid_t file,                // file to print <-
  const std::string & name,  // dataset name <- 
  hsize_t dimsf[],           // dataset dim  <-
  int data[]                 // vector data  <-
) { // ================================================

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                           );			   
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  // clean
  H5Sclose(dataspace);  H5Dclose(dataset);
  return status;
}

// =========================================================
// GambitIO  members
void GambitIO::read_mesh(
  std::istream& /*in*/
) { // =====================================================
  std::cout<< "Use read"; abort();
  return;
}


// =======================================================
void GambitIO::read(
   const std::string& name  // namefile <-
) {  // ==================================================

  std::set<boundary_id_type> bc_type;
 

  // initialize the map with element types
  init_eletypes();
  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh(); 
  mesh.clear();

  // Open the mesh file ---------------
  std::ifstream in(name.c_str());
  libmesh_assert(in.good());
    
  //  variables to read
  const int  bufLen = 256;  char  buf[bufLen+1];
  int format=0;
  int bcmat_flag=0;   // bcmat_flag 0(no) 1(YES) from file
  int size=0;
  int numElem = 0;   int numNodes = 0;
  int numGroups = 0; int numBCsets = 0;
  int numCoords = 0; int numfaces = 0;
  int numElem_nodes = 0;
  
  /// A) Setup
  while(strncmp(buf,"NDFVL",5) != 0) in >> buf;

  std::cout << " Reading from Gambit neutral mesh file" << std::endl;
  in >> numNodes >>  numElem >> numGroups >> numBCsets >> numCoords;
  std::cout << std::setw(12) << numNodes << std::setw(12)<< 
               numElem << std::setw(12)<< numCoords << std::endl;
  mesh.reserve_nodes(numNodes);  // mesh dim -> nodes
  mesh.reserve_elem(numElem);    // mesh dim -> elements
  int * mat_flag;  mat_flag = new int[numElem]; //material flag
  int * bc_flag;   bc_flag = new int[numNodes]; //boundary condition flag:

  /// B) coordinates *******************************************************
  Real x, y, z;        // cordinates
  unsigned int id;     // node id
  // map nodetrans to hold the node numbers for translation
  // note that the nodes can be non-consecutive or double
  std::map<unsigned int, unsigned int> nodetrans;

  // finding the location in the file
  while(strncmp(buf,"COORDINATES",11) != 0) in >> buf;
  in >> buf; //read version

  // read in the nodal coordinates and form points
  // and add the nodal coordinates to the mesh
  for(unsigned int i=0; i<numNodes; ++i) {
    if(numCoords==2) {        //  2D 
      in >>  id >> x >> y;
      mesh.add_point(Point(x, y),i);     
    } else if(numCoords==3) { //  3D
      in >> id >> x >> y >> z;
      mesh.add_point(Point(x, y, z),i);    
    }
    nodetrans[id-1] = i; // if not in order
    bc_flag[i]=0;        // bc_flag init
  }
  /// c) Elements *******************************************************
  while(strncmp(buf,"ELEMENTS/CELLS",14) != 0) in >> buf;
  in >> buf;

  // read the elements
  for(unsigned int iel=0; iel<numElem; ++iel) {
    unsigned int id, type, face, elementary=0, partition = 1, numTags;
    in >> id >> face >>  type;
    // number of nodes
    int numElem_nodes=type;
    int itype=0;
    // consult the import element table which element to build
    switch(face) {
    case 1: //"Edge"
      switch(type) {
      case 2: itype=1;  break; // edge 2
      case 3: itype=8;  break; // edge 4
      } break;

    case 2: //"Quadrilateral"
      switch(type) {
      case 4: itype=3;  break; // quad 4
      case 8: itype=19; break; // quad 8
      case 9: itype=10; break; // quad 9
      } break;

    case 3: //"Triangle"
      switch(type) {
      case 3: itype=2; break;  // tri 3
      case 6: itype=9; break;  // tri 6
//           case 7: itype=; break; // tri 7
      } break;

    case 4: //"Brick"
      switch(type) {
      case 8:  itype=5;  break;  // brick 8
      case 20: itype=20; break; // brick 20
      case 27: itype=12; break; // brick 27
      } break;
    case 5: //"Wedge (prism)"
      switch(type) {
      case 6:  itype=6;  break;  //
      case 9:  itype=103; break; //
      case 18: itype=13; break; //
      } break;

    case 6: //"Tetrahedron"
      switch(type) {
      case 4:  itype=4;  break; // tet 4
      case 10: itype=11; break; // tet 10
      } break;

    case 7: //"Pyramid"
      std::cout << "Pyramid geom elem not yet implemented" << std::endl; abort(); break;

    default: std::cout << "error element" << std::endl; abort();  break;
    }//end switch

    const elementDefinition& eletype = eletypes_imp[itype]; //alias of element type
    numfaces = eletype.nfaces;                              // face number
    // add the elements to the mesh
    Elem* elem = mesh.add_elem(Elem::build(eletype.type).release());

   // Connectivity
      int in_node=0; //  int* vecnode=new int [numElem_nodes];
      for(unsigned int i=0; i<numElem_nodes; i++)  {
        in>> in_node;
        elem->set_node(eletype.nodes[i])= mesh.node_ptr(nodetrans[in_node-1]);
      }
   mat_flag[iel]=1;   // element material flag (init ->1)
  } //end loop elem

  // bc and material *****************************************
  // material flag reading from file   (numGroups)
  int c_fl=0; int c_tmp=1;  int numb_c_fl=0;
  for(int i=0; i<(int) numGroups; ++i) {
    while(strncmp(buf,"GROUP:",6) != 0) in >> buf;
    in >> buf; in >> buf; in >> numb_c_fl;
    in >> buf;  in >> c_fl;
    in >> buf; in >> buf; in >> buf; in >> buf;
    for(int k=0; k<numb_c_fl; ++k) {
      in >> c_tmp;
      mat_flag[c_tmp-1]=c_fl;  //  -1 because C start from 0 and not from 1
    }
  }
  //boundary condition flag reading from file   (numBCsets)
  for(int i=0; i<(int) numBCsets; ++i) {
    while(strncmp(buf,"CONDITIONS",10) != 0) in >> buf;
    in >> buf; in >> c_fl; in >> buf; in >> numb_c_fl; in >> buf; in >> buf;
    for(int k=0; k<numb_c_fl; ++k) {
      in >> c_tmp;
      bc_flag[nodetrans[c_tmp-1]]=c_fl;
    }
  }

  
  // *****************************************************************************
  // Print material and boundary cond flag in ta file bc_case.mesh1.h5
  std::string femus_dir = getenv("FEMUS_DIR");
  std::string myapp_name = getenv("FM_MYAPP");

  std::ostringstream bc_casedir;
  if(msh_count==0)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh1.h5";
  else if(msh_count==1)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh2.h5";
  else if(msh_count==2)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh3.h5";
  else if(msh_count==3)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh4.h5";
  msh_count++;

  //    create a new file for storing material and boundary conditions flag
  hid_t file_id = H5Fcreate(bc_casedir.str().c_str(), 
			    H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   
  hsize_t     dimsf[2]; 
  dimsf[0]=numNodes; dimsf[1]=1;
  print_Ihdf5(file_id,"bc_flag",dimsf,bc_flag);
  dimsf[0]=numElem; dimsf[1]=1;
  print_Ihdf5(file_id,"material",dimsf,mat_flag);
  H5Fclose(file_id);
  // ******************************************************************
 
  // clean
  delete [] mat_flag;
  delete [] bc_flag;
  nodetrans.clear();
return;
}