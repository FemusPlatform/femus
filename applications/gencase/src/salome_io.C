// C++ includes ---------------------------------
#include <fstream>
#include <set>
#include <cstring>
#include <iomanip>
#include <cassert>

// config includes -------------------------------
#include "Solverlib_conf.h"    // lib configure 
#include "libmesh_config.h"
#include "MGFE_conf.h"

// Local includes -------------------------------
#include "elem.h"
#include "mesh_base.h"
#include "boundary_info.h"

// library includes -----------------------------
#include <hdf5.h>              // hdf5 libs

// local includes -------------------------------
#include "salome_io.h"

int msh_count2=0; // needed in order to use multiple meshes
// anonymous namespace to hold local data

namespace {
using namespace libMesh;

// Defines a structure to hold boundary element information.
struct boundaryElementInfo {
  std::set<unsigned int> nodes; // nodes
  unsigned int id;              // boundary groups
};

// Defines mapping from libMesh element types to Gmsh element types.
struct elementDefinition {
  std::string label;                // label
  std::vector<unsigned int> nodes;  // nodes
  ElemType type;                    // element by name
  unsigned int exptype;             // element by number
  unsigned int dim;                 // dim
  unsigned int nnodes;              // number of nodes
  unsigned int nfaces;              // number of faces
};


// maps from a libMesh element type to the proper
// salome elementDefinition.  Placing the data structure
// here in this anonymous namespace gives us the
// benefits of a global variable without the nasty
// side-effects
//map eletypes_exp: ElemType +  elementDefinition
std::map<ElemType, elementDefinition> eletypes_exp;
//map eletypes_exp:  int +  elementDefinition
std::map<unsigned int, elementDefinition> eletypes_imp;

//  Salome      LIBMESH
//             NAME  NUMB
//  ----------  1 DIM ----------------------------
//  SE2        EDGE2   1
//  SE3        EDGE3   8
//  ----------- 2 DIM ----------------------------
//  QU4        QUAD4   3
//  QU8        QUAD8   19
//  QU9        QUAD9   10
//  TR3        TRI3    2
//  TR6        TRI6    9


// ------------------------------------------------------------
// helper function to initialize the eletypes map
void init_eletypes() {
  if(eletypes_exp.empty() && eletypes_imp.empty()) {
    // This should happen only once.  The first time this method
    // is called the eletypes data struture will be empty, and
    // we will fill it.  Any subsequent calls will find an initialized
    // eletypes map and will do nothing.

    //==============================
    // setup the element definitions
    elementDefinition eledef;

    // use "swap trick" from Scott Meyer's "Effective STL" to initialize
    // eledef.nodes vector
    // 0 DIMENSIONAL GEOM ***************************************************************
    {
      eledef.dim     = 0;// POINT (only Gmsh) -------------------------------------------
      eledef.exptype = 15;        // element libmesh number
      eledef.nnodes  = 1;    eledef.nfaces  = 0;       // number of nodes/faces
      eledef.nodes.clear();       // clean function
      // class definition
      eletypes_imp[15] = eledef;  // element defined by number
    }
    // 1 DIMENSIONAL GEOM ***************************************************************
    {
      eledef.dim     = 1;   // SE2 -> EDGE2 (1) LIBMESH ---------------------------------
      eledef.type    = EDGE2;  eledef.exptype = 1; // element libmesh type/name
      eledef.nnodes  = 2;      eledef.nfaces  = 2; // number of nodes/faces
      eledef.nodes.clear();       // clean function
      // class definition (by name and number)
      eletypes_exp[EDGE2] = eledef; // def by name
      eletypes_imp[1]     = eledef; // def by number
    }
    {
      eledef.dim     = 1;    // SE3 Salome -> EDGE3 (8) LIBMESH ------------------------
      eledef.type    = EDGE3;  eledef.exptype = 8;     // element libmesh name/number
      eledef.nnodes  = 3;   eledef.nfaces  = 2;        // number of nodes/faces
      eledef.nodes.clear();      // clean function
      const unsigned int nodes[] = {0,2,1}; // map
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // class definition (by name and number)
      eletypes_exp[EDGE3] = eledef; // def by name
      eletypes_imp[8]     = eledef; // def by number
    }
    // 2 DIMENSIONAL GEOM ***************************************************************
    {
      eledef.dim     = 2;  // TR3 Salome -> TRI3 (2)  LIBMESH ---------------------------
      eledef.type    = TRI3;   eledef.exptype = 2;//element libmesh name/number
      eledef.nnodes  = 3;      eledef.nfaces  = 3;// number of nodes/faces
      // map
      const unsigned int nodes[] = {1,2,0};
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // class definition (by name and number)
      eletypes_exp[TRI3] = eledef; // def by name
      eletypes_imp[2] = eledef;    // def by number
    }
    {
      eledef.dim     = 2;  // TR6 Salome  -> TRI6 (9) LIBMESH ---------------------------
      eledef.type    = TRI6;   eledef.exptype = 9;  // element libmesh name/number
      eledef.nnodes  = 6;    eledef.nfaces  = 3;     // number of nodes/faces
      const unsigned int nodes[] = {1,2,0,4,5,3};//map
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // class definition (by name and number)
      eletypes_exp[TRI6] = eledef;
      eletypes_imp[9]    = eledef;
    }
    {
      eledef.dim     = 2; // QU4 Gambit  -> QUAD4 (3) LIBMESH ---------------------------
      eledef.type    = QUAD4;    eledef.exptype = 3;   // element libmesh name/number
      eledef.nnodes  = 4;    eledef.nfaces  = 4;      // number of nodes/faces
      const unsigned int nodes[] = {3,0,1,2}; //map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      // inverse map
      eletypes_exp[QUAD4] = eledef; // def by name
      eletypes_imp[3]     = eledef; // def by number
    }
    {
      eledef.dim     = 2;  // QU8 Salome -> QUAD4 (19) LIBMESH -------------------------
      eledef.type    = QUAD8;  eledef.exptype = 100;  // element libmesh name/number
      eledef.nnodes  = 8;     eledef.nfaces  = 4;   // number of nodes/faces
      const unsigned int nodes[] = {3,0,1,2,7,4,5,6};  // map
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[QUAD8] = eledef; // def by name
      eletypes_imp[19]    = eledef; // def by number
    }
    {
      eledef.dim     = 2;// QU9 Salome -> QUAD9 (10) LIBMESH ----------------------------
      eledef.type    = QUAD9;   eledef.exptype = 10; // element libmesh name
      eledef.nnodes  = 9;   eledef.nfaces  = 4;     // number of nodes
      const unsigned int nodes[] =  {3,0,1,2,7,4,5,6,8};// map
      // inverse map
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[QUAD9] = eledef; // def by name
      eletypes_imp[10]    = eledef; // def by number
    }
    // 3 DIMENSIONAL GEOM ***************************************************************
    {
      eledef.dim     = 3;  // HEX8 Gambit -----------------------------------------------
      eledef.type    = HEX8; eledef.exptype = 5;
      eledef.nnodes  = 8;
      const unsigned int nodes[] = {4,7,3,0,5,6,2,1};
//        const unsigned int nodes[] = {5,2,3,4,1,6,7,0};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[HEX8] = eledef;
      eletypes_imp[5]    = eledef;
    }
    {
      eledef.dim     = 3; // HEX20 Gambit -----------------------------------------------
      eledef.type    = HEX20; eledef.exptype = 101;
      eledef.nnodes  = 20;
//       const unsigned int nodes[] = {4,7,3,0,5,6,2,1,19,15,11,12,17,14,9,13,16,18,10,8};
      const unsigned int nodes[] = {4,7,3,0,5,6,2,1,19,15,11,12,17,14,9,13,16,18,10,8};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[HEX20] = eledef;
      eletypes_imp[20]    = eledef;
    }
    {
      eledef.dim     = 3;// HEX27 Salome ------------------------------------------------
      eledef.type    = HEX27; eledef.exptype = 12;
      eledef.nnodes  = 27;
      const unsigned int nodes[] = {4,7,3,0,5,6,2,1,19,15,11,12,17,14,9,13,16,18,
                                    10,8,24,25,23,20,21,22,26
                                   };
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[HEX27] = eledef;
      eletypes_imp[12]    = eledef;
    }
    {
      eledef.dim     = 3; // TET4 Gambit -------------------------------------------------
      eledef.type    = TET4; eledef.exptype = 4;
      eledef.nnodes  = 4;
      eledef.nodes.clear();
      const unsigned int nodes[] = {2,3,1,0};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[TET4] = eledef;
      eletypes_imp[4]    = eledef;
    }
    {
      eledef.dim     = 3;  // TET10 Gambit ----------------------------------------------
      eledef.type    = TET10; eledef.exptype = 11;
      eledef.nnodes  = 10;
      const unsigned int nodes[] = {2,3,1,0,9,8,5,6,7,4};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[TET10] = eledef;
      eletypes_imp[11]    = eledef;
    }
    {
      eledef.dim     = 3;// PRISM6 Gambit ----------------------------------------------
      eledef.type    = PRISM6; eledef.exptype = 6;
      eledef.nnodes  = 6;
      eledef.nodes.clear();
      const unsigned int nodes[] = {2,0,1,5,3,4};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[PRISM6] = eledef;
      eletypes_imp[6]      = eledef;
    }
    {
      eledef.dim     = 3; // PRISM15 Gambit -------------------------------------------------------
      eledef.type    = PRISM15; eledef.exptype = 103;
      eledef.nnodes  = 15;
      eledef.nodes.clear();
      const unsigned int nodes[] = {4,12,3,13,14,5,10,9,11,1,6,0,7,8,2};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);
      eletypes_exp[PRISM15] = eledef;
      eletypes_imp[16] = eledef;
    }
    {
      eledef.dim     = 3;// PRISM18 Gambit ----------------------------------------------
      eledef.type    = PRISM18; eledef.exptype = 13;
      eledef.nnodes  = 18;
      //   const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,
      //                            12,14,13,15,17,16};
      const unsigned int nodes[] = {5,13,4,14,12,3,11,16,10,17,15,9,2,7,1,8,6,0};
      std::vector<unsigned int> (nodes, nodes+eledef.nnodes).swap(eledef.nodes);

      eletypes_exp[PRISM18] = eledef;
      eletypes_imp[13]      = eledef;
    }
    {
      eledef.dim     = 3;// PYRAMID5 ----------------------------------------------------
      eledef.type    = PYRAMID5;  eledef.exptype = 7;
      eledef.nnodes  = 5;
      eledef.nodes.clear();
      eletypes_exp[PYRAMID5] = eledef;  eletypes_imp[7] = eledef;

    }

    //==============================
  }
}
} // end anonymous namespace
// =========================================================


// =========================================================
void SalomeIO::read_mesh(
  std::istream&    // istream file to read  <-
) { // =====================================================
  std::cout<< "Use read"; abort();
  return;
}

// ============================================================================
// SalomeIO  members
void SalomeIO::read(
  const std::string& namefile  // namefile <-
) {// =========================================================================

  // *************  setup *****************************************************
  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();
  // Open the mesh file. ---------------
  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dims[2]; double xyz3[DIMENSION];
  std::string el_fem_type_vol(""); std::string el_fem_type_bd("");



  // LIBMESH volume fem type and initialization***********************************
  // reading  fem type from hdf5-med file
  int itype_vol=  read_fem_type(file_id,el_fem_type_vol,el_fem_type_bd);
  // libmesh element -> eletype+ Node_el from above elementDefinition structure
  init_eletypes(); // initalization
  const elementDefinition& eletype = eletypes_imp[itype_vol];  // libmesh element type
  int Node_el=eletype.nnodes;                              // libmesh nodes for element
  std::cout << " Family element med ibmesh = " << itype_vol << std::endl;

  // Coordinates -----------------------------------------------
  // Getting dataset xyz
  hid_t dtset = H5Dopen(file_id,COORD_NAME_DIR
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) std::cerr << "SalomeIO::read dims not found";
  // reading xyz_med
  unsigned int n_nodes =dims[0]/DIMENSION;
  double   *xyz_med=new double[dims[0]];
  std::cout << " Number of points in med file =" <<  n_nodes << " " <<  std::endl;
  status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xyz_med);
  H5Dclose(dtset);
  // libmesh storage
  mesh.reserve_nodes(n_nodes);
  for(int j=0; j<(int)n_nodes; j++) {
    for(int idim=0; idim<DIMENSION; idim++) xyz3[idim]= xyz_med[j+n_nodes*idim];  // coordinate table conversion
    if(DIMENSION==3) mesh.add_point(Point(xyz3[0],xyz3[1],xyz3[2]),j);     //3D
    else mesh.add_point(Point(xyz3[0],xyz3[1]),j); //2D
  }
  delete[] xyz_med;


  // Connectivity ---------------------------------------------------
  // Getting  connectivity structure from file *.med
  std::string node_name_dir=MESH_NAME_DIR+el_fem_type_vol/*(std::string) el_fem_type[index_vol]*/
                            +"/NOD";
  dtset = H5Dopen(file_id,node_name_dir.c_str()
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT
#endif
                 );
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) {std::cerr << "SalomeIO::read dims not found"; abort();}
  const int dim_conn=dims[0];

  // reading connectivity
  unsigned int n_elements =dim_conn/Node_el; int   *conn_map5=new  int[dim_conn];
  std::cout << " Number of elements " <<  n_elements << " " <<  std::endl;
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,conn_map5);
  if(status !=0) {std::cout << "SalomeIO::read: connectivity not found"; abort();}
  H5Dclose(dtset);

  // Adding Connectivity to mesh structure (Libmesh)
  mesh.reserve_elem(n_elements);
  // read the elements
  for(int iel=0; iel<(int)n_elements; ++iel) {
    // add the elements to the mesh
    Elem* elem = mesh.add_elem(Elem::build(eletype.type).release());
    // add node pointers to the elements
    for(int i=0; i<Node_el; i++) {
//       elem->set_node(eletype.nodes[i])= mesh.node_ptr(conn_map[Node_el*iel+i]);
      elem->set_node(eletype.nodes[i])= mesh.node_ptr(conn_map5[iel+i*n_elements]-1);
    }
  }
  // clean
  delete[] conn_map5;  /*delete []conn_map;*/


 
  // reading bc mad mat groups from salome ************************************
  read_bc_mat(file_id,el_fem_type_vol /*(std::string) el_fem_type[index_vol]*/,
              el_fem_type_bd /*(std::string) el_fem_type[index_bd]*/,
//               bc_flag,
//               mat_flag, 
n_elements,
 n_nodes
             );
  //   for (int i=0; i<n_elements; i++) std::cout<<i+1<<"  " <<bc_flag[i]<<"\n";
  //   for (int i=0; i<n_nodes; i++)    std::cout<<i+1<<"  "<<mat_flag[i]<<"\n";
  
 
  H5Fclose(file_id);

  // mesh count
//   std::string femus_dir = getenv("FEMUS_DIR");
//   std::string myapp_name = getenv("FM_MYAPP");
//   std::ostringstream bc_casedir;
//   if(msh_count2%2==0)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh1.h5";
//   else if(msh_count2%2==1)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh2.h5";
//   else if(msh_count2%2==2)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh3.h5";
//   else if(msh_count2%2==3)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh4.h5";
//   msh_count2++;



  return;
}
// // ============================================================================
// // This function reads the boundary condition and mat groups
// void SalomeIO::read_bc_mat(
//   hid_t file_id,             // file id (hdf5)
//   std::string  el_fem_type,  // vol fem salome
//   std::string elbd_fem_type,  // boundary fem salome
//   int n_elements,            // n of elements
//   int n_nodes,                // n of nodes
//   int & bc_flag,
// 	 int &     mat_flag
// ) { // =========================================================================
//
//   // material and bc flags
//   int *mat_flag = new int[n_elements];      // material flag
//   int *bc_flag = new int[n_nodes];       // boundary condition flag
//   for(int k=0; k<(int) n_nodes; k++) bc_flag[k]=0;  // initializing
//   for(int k=0; k<(int) n_elements; k++) mat_flag[k]=1;
//
//   // Boundary groups (boundary +material) ****************************************
//   // test if boundary groups are in the file ------------------------------------
// //   int bcmat_flag=0;                          // bcmat_flag 0(no) 1(YES) from file
//   std::string top_name="/FAS/Mesh_1/";
//   hsize_t     n_bdgroup=0;
//   std::map<int, int > bdgroup_name; // map of boundary names --> id
//   char name_dir[20], name_sub[20];  // link name (Group**_**_***)
//   char s2[81];           // char
//   hsize_t    dim[] = {1};   /* Dataspace dimensions */
//   hsize_t    array_dim[] = {80};   /* Array dimensions */
//
//   hid_t gid=H5Gopen(file_id, top_name.c_str()); // group identity
//   hid_t status= H5Gget_num_objs(gid, &n_bdgroup);  // number of links
//   if(status !=0) std::cerr <<  "SalomeIO::read: number of bdgroup not found \n";
//   hsize_t     n_bdgroup_sub=0;
//   for(int iname_dir=0; iname_dir<(int)n_bdgroup; iname_dir++) {
//     H5Gget_objname_by_idx(gid, iname_dir ,&name_dir[0], 20);
//     std::string a1=top_name+name_dir;
//     int gid_sub=H5Gopen(file_id,a1.c_str()); // group identity
//     status= H5Gget_num_objs(gid_sub, &n_bdgroup_sub);  // number of links
//     if(status !=0) std::cerr <<  "SalomeIO::read: number of bdgroup not found \n";
//     if(n_bdgroup_sub>0) {
//       for(int iname_sub=0; iname_sub<(int)n_bdgroup_sub; iname_sub++) {
//         H5Gget_objname_by_idx(gid_sub, iname_sub ,&name_sub[0], 20);
//         std::string a2=top_name+name_dir+"/"+name_sub+"/GRO/";
//         int gid3=H5Gopen(file_id,a2.c_str());
//
//         //  char array format
//         hid_t array_tid = H5Tarray_create(H5T_NATIVE_CHAR, dim[0],array_dim, NULL);
//         hid_t  dataset = H5Dopen(gid3,"NOM"
// #if HDF5_VERSIONM != 1808
//                                  ,H5P_DEFAULT
// #endif
//                                 );
// //     status = H5Dwrite(dataset, array_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s2);
//         status = H5Dread(dataset,array_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &s2);
//         std::string str2(s2);
//         std::string as(name_sub);                 // name format  ***_**_***
//         unsigned pos1 = as.find("_");         // position of "_" in link name
//         std::string str1 = as.substr(pos1+1); // string to the end  **_***
//         unsigned pos2 = str1.find("_");       // position of "_" in link name
//         std::string str3=as.substr(pos1+1,pos2);
//
//         bdgroup_name[atoi(str3.c_str())]=atoi(str2.c_str()) ;
//         std::cout << "SalomeIO::read Group name " << as << "= group " << str3 << " id " <<  bdgroup_name[atoi(str3.c_str())] << '\n';
//
//       }
//     }
//   }
//   // Reading surface structure (bc) ******************************************
//   //   FEMEL (SURF)
//   //   |- FAM
//   //   |- NOD
//   //   |- NUM
//   //  getting  FAM dataset ------------------------------------
//   std::string path_fam= MESH_NAME_DIR + (std::string) elbd_fem_type+(std::string) "/FAM";
//   hid_t dtset = H5Dopen(file_id,path_fam.c_str()
// #if HDF5_VERSIONM != 1808
//                         ,H5P_DEFAULT
// #endif
//                        );
//   hsize_t dims1[2];
//   hid_t  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
//   status  = H5Sget_simple_extent_dims(filespace, dims1, NULL);
//   int *fam_vector=new int[dims1[0]];
//   status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,fam_vector);
//   if(status !=0) std::cout << "SalomeIO::read: FAM not found";
//   H5Dclose(dtset);
//
//   //  getting  NOD dataset -----------------------------------
//   std::string path_nod= MESH_NAME_DIR + elbd_fem_type+(std::string) "/NOD";
//   dtset = H5Dopen(file_id,path_nod.c_str()
// #if HDF5_VERSIONM != 1808
//                   ,H5P_DEFAULT
// #endif
//                  );
//   hsize_t dims2[2];
//   filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
//   status    = H5Sget_simple_extent_dims(filespace, dims2, NULL);
//   int *nod_vector=new int[dims2[0]];
//   status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,nod_vector);
//   if(status !=0) std::cout << "SalomeIO::read: NOD not found";
//   H5Dclose(dtset);
//   int n_node_bdel=dims2[0]/dims1[0];
//
//   // boundary flag -------------------------------
//   for(int j=0; j<(int) n_node_bdel; j++) {
//     for(int i=0; i<(int) dims1[0]; i++) {
//       int knode=nod_vector[i+j*dims1[0]]-1;//order[nod_vector[i+j*dims1[0]]-1];
//       if(bc_flag[knode]<bdgroup_name[fam_vector[i]] &&bdgroup_name[fam_vector[i]] >9)
//         bc_flag[knode]= bdgroup_name[fam_vector[i]];
//     }
//   }
//   // clean --------------------------------------
//   delete []fam_vector;  delete []nod_vector;
//
// // Reading volume structure (mat)***********************************************
//   //   FEMEL (VOL)
//   //   |- FAM
//   //   |- NOD
//   //   |- NUM
//   //  getting  FAM dataset for material ------------------------
//   path_fam= MESH_NAME_DIR + el_fem_type+(std::string)"/FAM";
//   dtset = H5Dopen(file_id,path_fam.c_str()
// #if HDF5_VERSIONM != 1808
//                   ,H5P_DEFAULT
// #endif
//                  );
//   filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
//   status  =   H5Sget_simple_extent_dims(filespace, dims1, NULL);
//   fam_vector=new int[dims1[0]];
//   status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,fam_vector);
//   if(status !=0) std::cout << "SalomeIO::read: VOLEL/FAM not found";
//   H5Dclose(dtset);
//
//   // material flag ----------------------------
//   for(int i=0; i<(int) dims1[0]; i++) {
//     int knode =i;
//     if(mat_flag[knode]<bdgroup_name[fam_vector[i]] &&bdgroup_name[fam_vector[i]] <9)
//       mat_flag[knode]= bdgroup_name[fam_vector[i]];
//   }
//
//   // clean
//   delete []fam_vector; bdgroup_name.clear();
// //   for (int i=0; i<n_elements; i++) std::cout<<i+1<<"  " <<bc_flag[i]<<"\n";
// //   for (int i=0; i<n_nodes; i++)    std::cout<<i+1<<"  "<<mat_flag[i]<<"\n";
//
//   std::string femus_dir = getenv("FEMUS_DIR");
//   std::string myapp_name = getenv("FM_MYAPP");
//   std::ostringstream bc_casedir;
//   if(msh_count2%2==0)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh1.h5";
//   else if(msh_count2%2==1)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh2.h5";
//   else if(msh_count2%2==2)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh3.h5";
//   else if(msh_count2%2==3)  bc_casedir << femus_dir <<  "/USER_APPL/"  <<  myapp_name  <<"/RESU/bc_case.msh4.h5";
//   msh_count2++;
//
//   //  create a new file for storing material and boundary conditions flag
//   file_id = H5Fcreate(bc_casedir.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   hsize_t     dimsf[2];
//   dimsf[0]=n_nodes; dimsf[1]=1;
//   print_Ihdf5(file_id,"bc_flag",dimsf,bc_flag);
//   dimsf[0]=n_elements; dimsf[1]=1;
//   print_Ihdf5(file_id,"material",dimsf,mat_flag);
//
//   //clean
//   delete [] mat_flag;  delete [] bc_flag;
//   H5Fclose(file_id);
//
//   return;
// }












// ============================================================================
// This function reads the boundary condition and mat groups
void SalomeIO::read_bc_mat(
  hid_t file_id,             // file id (hdf5)
  std::string  el_fem_type,  // vol fem salome
  std::string elbd_fem_type,  // boundary fem salome
  int n_elements,
  int n_nodes
//   int *bc_flag,
//   int *mat_flag,
//   std::map<int, int > bdgroup_name // map of boundary names --> id
) { // =========================================================================


 // material and bc flags
  int *mat_flag = new int[n_elements];      // material flag
  int *bc_flag = new int[n_nodes];          // boundary condition flag
  for(int k=0; k<(int) n_nodes; k++) bc_flag[k]=0;  // initializing
  for(int k=0; k<(int) n_elements; k++) mat_flag[k]=1;
  std::map<int,int> bdgroup_name; // map of boundary names --> id
  // Boundary groups (boundary +material) ****************************************
  // test if boundary groups are in the file ------------------------------------
//   int bcmat_flag=0;                          // bcmat_flag 0(no) 1(YES) from file
  std::string top_name="/FAS/Mesh_1/";
  hsize_t     n_bdgroup=0;
//   std::map<int, int > bdgroup_name; // map of boundary names --> id
  char name_dir[20], name_sub[20];  // link name (Group**_**_***)
  char s2[81];           // char
  hsize_t    dim[] = {1};   /* Dataspace dimensions */
  hsize_t    array_dim[] = {80};   /* Array dimensions */

  hid_t gid=H5Gopen(file_id, top_name.c_str(),
#if           HDF5_VERSIONM!=188          
                    H5P_DEFAULT
#endif                    
                   ); // group identity
  hid_t status= H5Gget_num_objs(gid, &n_bdgroup);  // number of links
  if(status !=0) std::cerr <<  "SalomeIO::read: number of bdgroup not found \n";
  hsize_t     n_bdgroup_sub=0;
  for(int iname_dir=0; iname_dir<(int)n_bdgroup; iname_dir++) {
    H5Gget_objname_by_idx(gid, iname_dir ,&name_dir[0], 20);
    std::string a1=top_name+name_dir;
    int gid_sub=H5Gopen(file_id,a1.c_str(),
#if HDF5_VERSIONM!=188          
                    H5P_DEFAULT
#endif    
    ); // group identity
    status= H5Gget_num_objs(gid_sub, &n_bdgroup_sub);  // number of links
    if(status !=0) std::cerr <<  "SalomeIO::read: number of bdgroup not found \n";
    if(n_bdgroup_sub>0) {
      for(int iname_sub=0; iname_sub<(int)n_bdgroup_sub; iname_sub++) {
        H5Gget_objname_by_idx(gid_sub, iname_sub ,&name_sub[0], 20);
        std::string a2=top_name+name_dir+"/"+name_sub+"/GRO/";
        int gid3=H5Gopen(file_id,a2.c_str(),
#if HDF5_VERSIONM!=188          
                    H5P_DEFAULT
#endif
        );

        //  char array format
        hid_t array_tid = H5Tarray_create(H5T_NATIVE_CHAR, dim[0],array_dim
#if HDF5_VERSIONM==188          
                    ,NULL
#endif
         );
        hid_t  dataset = H5Dopen(gid3,"NOM"
#if HDF5_VERSIONM != 188
                                 ,H5P_DEFAULT
#endif
                                );
//     status = H5Dwrite(dataset, array_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s2);
        status = H5Dread(dataset,array_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &s2);
        std::string str2(s2);
        std::string as(name_sub);                 // name format  ***_**_***
        unsigned pos1 = as.find("_");         // position of "_" in link name
        std::string str1 = as.substr(pos1+1); // string to the end  **_***
        unsigned pos2 = str1.find("_");       // position of "_" in link name
        std::string str3=as.substr(pos1+1,pos2);

        bdgroup_name[atoi(str3.c_str())]=atoi(str2.c_str()) ;
        std::cout << "SalomeIO::read Group name " << as << "= group " << str3 << " id " <<  bdgroup_name[atoi(str3.c_str())] << '\n';

      }
    }
  }
  // Reading surface structure (bc) ******************************************
  //   FEMEL (SURF)
  //   |- FAM
  //   |- NOD
  //   |- NUM
  //  getting  FAM dataset ------------------------------------
  std::string path_fam= MESH_NAME_DIR + (std::string) elbd_fem_type+(std::string) "/FAM";
  hid_t dtset = H5Dopen(file_id,path_fam.c_str()
#if HDF5_VERSIONM != 1808
                        ,H5P_DEFAULT
#endif
                       );
  hsize_t dims1[2];
  hid_t  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  status  = H5Sget_simple_extent_dims(filespace, dims1, NULL);
  int *fam_vector=new int[dims1[0]];
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,fam_vector);
  if(status !=0) std::cout << "SalomeIO::read: FAM not found";
  H5Dclose(dtset);

  //  getting  NOD dataset -----------------------------------
  std::string path_nod= MESH_NAME_DIR + elbd_fem_type+(std::string) "/NOD";
  dtset = H5Dopen(file_id,path_nod.c_str()
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT
#endif
                 );
  hsize_t dims2[2];
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  status    = H5Sget_simple_extent_dims(filespace, dims2, NULL);
  int *nod_vector=new int[dims2[0]];
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,nod_vector);
  if(status !=0) std::cout << "SalomeIO::read: NOD not found";
  H5Dclose(dtset);
  int n_node_bdel=dims2[0]/dims1[0];

  // boundary flag -------------------------------
  for(int j=0; j<(int) n_node_bdel; j++) {
    for(int i=0; i<(int) dims1[0]; i++) {
      int knode=nod_vector[i+j*dims1[0]]-1;//order[nod_vector[i+j*dims1[0]]-1];
      if(bc_flag[knode]<bdgroup_name[fam_vector[i]] &&bdgroup_name[fam_vector[i]] >9)
        bc_flag[knode]= bdgroup_name[fam_vector[i]];
    }
  }
  // clean --------------------------------------
  delete []fam_vector;  delete []nod_vector;

// Reading volume structure (mat)***********************************************
  //   FEMEL (VOL)
  //   |- FAM
  //   |- NOD
  //   |- NUM
  //  getting  FAM dataset for material ------------------------
  path_fam= MESH_NAME_DIR + el_fem_type+(std::string)"/FAM";
  dtset = H5Dopen(file_id,path_fam.c_str()
#if HDF5_VERSIONM != 1808
                  ,H5P_DEFAULT
#endif
                 );
  filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  status  =   H5Sget_simple_extent_dims(filespace, dims1, NULL);
  fam_vector=new int[dims1[0]];
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,fam_vector);
  if(status !=0) std::cout << "SalomeIO::read: VOLEL/FAM not found";
  H5Dclose(dtset);

  // material flag ----------------------------
  for(int i=0; i<(int) dims1[0]; i++) {
    int knode =i;
    if(mat_flag[knode]<bdgroup_name[fam_vector[i]] &&bdgroup_name[fam_vector[i]] <9)
      mat_flag[knode]= bdgroup_name[fam_vector[i]];
  }

   // print bc and mat conditions into the file "bc_case*.*"-----------------------------
  // the file name
  std::string app_path = getenv("APP_PATH");
  std::ostringstream bc_casedir;


  int num_bc_case = msh_count2/2 + 1;
  std::string num_bc = std::to_string(num_bc_case);
  std::string filename = "/RESU/bc_case.msh" + num_bc + ".h5";
  
  bc_casedir << app_path << filename;
  msh_count2++;
  
  
  //  create a new file for storing material and boundary conditions flag
  hid_t file_id2 = H5Fcreate(bc_casedir.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t     dimsf[2];
  // print bc_flag
  dimsf[0]=n_nodes; dimsf[1]=1;  print_Ihdf5(file_id2,"bc_flag",dimsf,bc_flag);
   // print bc_flag
  dimsf[0]=n_elements; dimsf[1]=1;  print_Ihdf5(file_id2,"material",dimsf,mat_flag);
  
  // family identities and group names
  dimsf[0]=bdgroup_name.size(); dimsf[1]=1; int icount=0;
  // family vectors
  int *family_id=new int[dimsf[0]];  int *group_names=new int[dimsf[0]];
  std::map<int,int>::iterator itr;
  for(itr=bdgroup_name.begin();itr!=bdgroup_name.end();++itr){
    family_id[icount]=itr->first;   group_names[icount]=itr->second;icount++;
  }
  // print family identities and group names
  print_Ihdf5(file_id2,"family_id",dimsf,family_id);
  print_Ihdf5(file_id2,"group_names",dimsf,group_names);
  
  //clean
  delete[] mat_flag; delete[] bc_flag;
  delete family_id; delete[]group_names;
   H5Fclose(file_id2);
  
  // clean
  delete []fam_vector; bdgroup_name.clear();
//   for (int i=0; i<n_elements; i++) std::cout<<i+1<<"  " <<bc_flag[i]<<"\n";
//   for (int i=0; i<n_nodes; i++)    std::cout<<i+1<<"  "<<mat_flag[i]<<"\n";


  return;
}









// ============================================================================
// This function reads the boundary condition and mat groups
int  SalomeIO::read_fem_type(
  hid_t file_id,
  std::string & el_fem_type_vol,
  std::string & el_fem_type_bd
) {

// ************  reading  hdf5-med file  ************************************

  // ---------  Reading Element type (From salome to LIBMESH (itype_vol))------
  std::map< std::string, int > fem_type_vol;// Salome fem table name (vol)
  fem_type_vol["HE8"] = 5;  fem_type_vol["H20"] = 20; fem_type_vol["H27"] = 12;
  fem_type_vol["TE4"] = 4;  fem_type_vol["T10"] = 11;
  fem_type_vol["QU4"] = 3;  fem_type_vol["QU8"] = 19; fem_type_vol["QU9"] = 10;
  fem_type_vol["TR3"] = 2;  fem_type_vol["TR6"] = 9;
  fem_type_vol["SE2"] = 0;  fem_type_vol["SE3"] = 0;  // not valid in 3D
  if(DIMENSION==3) {    // not valid in 3D as boundary
    fem_type_vol["QU4"] = 0; fem_type_vol["QU8"] = 0; fem_type_vol["QU9"] = 0;
    fem_type_vol["TR3"] = 0; fem_type_vol["TR6"] = 0;
  }
  std::map< std::string, int > fem_type_bd; // Salome fem table name (surface)
  fem_type_bd["HE8"] = 0;  fem_type_bd["H20"] = 0;  fem_type_bd["H27"] = 0;
  fem_type_bd["TE4"] = 0;  fem_type_bd["T10"] = 0;  // not valid in 2D
  fem_type_bd["QU4"] = 3;  fem_type_bd["QU8"] = 19;  fem_type_bd["QU9"] = 10;
  fem_type_bd["TR3"] = 2;  fem_type_bd["TR6"] = 9;
  fem_type_bd["SE2"] = 1;  fem_type_bd["SE3"] = 8;
  if(DIMENSION==3) {    // not valid in 3D as boundary
    fem_type_bd["SE2"] = 0;  fem_type_bd["SE3"] = 0;
  }
    if(DIMENSION==2) {    // not valid in 3D as boundary
    fem_type_bd["TR6"] = 0;fem_type_bd["QU8"] = 0;fem_type_bd["QU9"] = 0;
    fem_type_bd["QU4"] = 0; fem_type_bd["TR3"] = 0;
  }

  hid_t       gid=H5Gopen(file_id,MESH_NAME_DIR,
#if HDF5_VERSIONM!=188          
                    H5P_DEFAULT
#endif
  ); // group identity
  hsize_t     n_fem_type;
  hid_t status= H5Gget_num_objs(gid, &n_fem_type);  // number of links
  if(status !=0) {std::cout << "SalomeIO::read_fem_type:   H5Gget_num_objs not found"; abort();}

  // Get the element name from MESH_NAME_DIR in the file med (el_name)
  char **el_fem_type=new char*[n_fem_type];
  int index_vol=0;  int index_bd=0;
  for(int i=0; i<(int)n_fem_type; i++) {
    el_fem_type[i]=new char[4];
    H5Lget_name_by_idx(file_id,MESH_NAME_DIR, H5_INDEX_NAME, H5_ITER_INC,i,el_fem_type[i],4, H5P_DEFAULT);
    if(fem_type_vol[el_fem_type[i]]!=0) index_vol=i;
    if(fem_type_bd[el_fem_type[i]]!=0) index_bd=i;
  }

  // LIBMESH volume fem type (itype_vol) and MEDfem (el_fem_type_vol-el_fem_type_bd)
  int itype_vol= fem_type_vol[el_fem_type[index_vol]]; assert(itype_vol!=0);
  el_fem_type_vol=el_fem_type[index_vol];
  el_fem_type_bd=el_fem_type[index_bd];
  // clean
  for(int i=0; i<(int)n_fem_type; i++) delete[] el_fem_type[i];
  delete[] el_fem_type; fem_type_vol.clear(); fem_type_bd.clear();
  return itype_vol;
}


// ============================================================================
/// Print int data into dhdf5 file
hid_t SalomeIO::print_Ihdf5(
  hid_t file,               // file to print <-
  const std::string & name, // dataset name <-
  hsize_t dimsf[],          // dataset dim  <-
  int data[]                // vector data  <-
) {// =========================================================================

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                           );
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  /* clean */  H5Sclose(dataspace);  H5Dclose(dataset);
  return status;
}


