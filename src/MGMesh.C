// std lib
#include <iomanip>
#include <sstream>
/*#include <iomanip>
#include <sstream>*/
#include <cmath>

// conf includes ------------------
#include "Solverlib_conf.h"  // petsc conf
#include "MGFE_conf.h"        // FEM conf
#include "Printinfo_conf.h"  // print config file

// class includes ---------------------
#include "MGMesh.h"
#include "MGUtils.h"
#include "MGGeomEl.h"
#include "MeshExtended.h"

// // #include "parallelM.h"
#ifdef HAVE_PETSCM
# include <mpi.h>  // This is needed in the constructor 
#endif


// ====================================================
/// This function is the mesh constructor
MGMesh::MGMesh (
  const ParallelObjectM & comm,
  MGUtils & mgutils_in, // input file name class
  MGGeomEl & geomel_in // element geometry class input
//   const double Lref    // reference length input
) :
  _n_GeomEl ( NUMGEOMELS ), // element geometry class
  _GeomEl ( geomel_in ), //  element geometry class
  _mgutils ( mgutils_in ), // file name class
  _comm ( comm ),       // communicator
  _Lref ( LREF )
  { // reference length
  // ==================================

  _ord_FEM = new int [27];
  idFEMord ( 27 );

  // set subdomain number
  int i;
#ifdef HAVE_PETSCM
  MPI_Comm_rank ( MPI_COMM_WORLD, &i );
  _iproc = static_cast< int> ( i );
#ifdef PRINT_INFO
  std::cout <<  " ===== Mesh: ===== \n communicator: " << MPI_COMM_WORLD << std::endl;
  printf ( " \n %d proc=: ", _iproc );
#endif
#else
  _iproc = 0;
#endif
  this->read_c(); // reading
  return;
  }

// ====================================================
/// This function is the mesh destructor
MGMesh::~MGMesh (
)  {// ================================================

  clear();
  delete[] _NoElements;
  delete[] _xyz;
  delete[] _xyzo;
  delete[] _dxdydz;
  delete[] _NoNodes;

  delete[] _type_FEM;
  delete[] _el_map;  //  delete[]_conn_map;
  delete[] _off_el;
//   delete[] _off_nd;
//   delete[] _off_nd1;
  delete [] _ord_FEM;
  delete [] _dist;
  delete [] _VolFrac;
  delete [] _ctrl_dom;
  delete [] _el_neighbor;
  delete [] _node_map;
  delete [] _NodeBDgroup;
  }


// ====================================================
/// This function destructs the mesh substructures
void MGMesh::clear (
)  {// ============================================

  for ( int  imesh = 0; imesh < _NoLevels + 1; imesh++ )    delete []_node_map[imesh];

//   for(int  imesh =0; imesh<_NoLevels; imesh++)   delete []_dist; //[imesh];

  delete[] _off_nd[0];
  delete[] _off_nd[1];

  for ( int  imesh = 0; imesh < _NoFamFEM; imesh++ ) {
      delete [] _el_map[imesh];
      delete [] _off_el[imesh];
      delete []_NoElements[imesh];
      }

  for ( int  imesh = 0; imesh < _NoFamFEM - 1; imesh++ ) {
      delete [] _el_neighbor[imesh]; //now _el_neighbor is allocated only for the volume family
      }

  return;
  }
// ====================================================
/// This function computes the element center
void MGMesh::get_el_ctr (
  const int  el_nnodes,   // element nodes  (<- input)=_GeomEl.n_q[bdry]
  const int  bdry,        // zone vol/bd    (<- input)=(1)/(0) flag
  const double * xx_nds,  // element coords (-> output)
  double * el_xm          // element center (-> output)
) const {// =======================================

  for ( int  idim = 0; idim < _dim; idim++ ) {
      el_xm[idim] = 0.; // zero

      for ( int  eln = 0; eln < el_nnodes; eln++ )   el_xm[idim]   += xx_nds[eln + idim * el_nnodes];

      el_xm[idim] = el_xm[idim] / el_nnodes; // normalization
      }

  return;
  }

// =======================================
/// This function gets the global node numbers for that element
/// and their coordinates
void MGMesh::get_el_nodes (
  const int  el_nnodes,  // element nodes(<- input)= _GeomEl.n_q[bdry];
  const int  bdry,      // zone vol/bd   (<- input)=(1)/(0) flag
  const int  Level,     // level         (<- input)
  const int  iel,       // element       (<- input)
  double xx[]           // element coords(-> output)
) const {// =======================================


  // element  definition
  for ( int  n = 0; n < el_nnodes; n++ )    {
      //get the global node number
      int  el_conn = _el_map[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * _iproc] ) * el_nnodes + n];

      // get the element coordinades
      for ( int  idim = 0; idim < _dim; idim++ ) {
          const int  indxn = n + idim * el_nnodes;
          xx[indxn] = _xyz[el_conn + idim * _NoNodes[_NoLevels - 1]];
          }
      }

  return;
  }
// ===============================================================
/// This function return the element connectivity and  coordinates
// ==============================================================
void MGMesh::get_el_nod_disp (
  const int  bdry,  // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
  const int  Level, // level                (<- input)
  const int  iel,   // element              (<- input)
//   int  el_conn[],   // element connectivity (<- intput)
  double xx[]       // element coordinates  (-> output)
) const {// =================================

  const int  el_nnodes = _GeomEl.n_q[bdry]; // element nodes
//   const int  offset   =_NoNodes[_NoLevels-1];
//   // element  definition
//   for(int  inode=0; inode<el_nnodes; inode++)    {
//     //get the global node number
// //     el_conn[inode] = _el_map[bdry][(iel+_off_el[bdry][Level+_NoLevels*_iproc])*el_nnodes+inode];
//     // get the element coordinades
//     for(int  idim=0; idim<_dim; idim++) {
//       xx[inode+idim*el_nnodes] = _dxdydz[el_conn[inode]+idim*offset];
// //       xx[inode+idim*el_nnodes] = _xyz[el_conn[inode]+idim*offset]+_dxdydz[el_conn[inode]+idim*offset];
//
// //        xx[inode+idim*el_nnodes] = _dxdydz[el_conn[inode]+idim*offset];
//     }
//   }

  // element  definition
  for ( int  n = 0; n < el_nnodes; n++ )    {
      //get the global node number
      int  el_conn = _el_map[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * _iproc] ) * el_nnodes + n];

      // get the element coordinades
      for ( int  idim = 0; idim < _dim; idim++ ) {
          const int  indxn = n + idim * el_nnodes;
          xx[indxn] = _dxdydz[el_conn + idim * _NoNodes[_NoLevels - 1]];
          }
      }

  return;
  }
// ===============================================================
/// This function return the element connectivity and  coordinates
// ==============================================================
void MGMesh::get_el_nod_conn (
  const int  bdry,  // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
  const int  Level, // level                (<- input)
  const int  iel,   // element              (<- input)
  int  el_conn[],   // element connectivity (-> output)
  double xx[],       // element coordinates (-> output)
  int  iproc
) const {// =================================

  const int  el_nnodes = _GeomEl.n_q[bdry]; // element nodes
  const int  offset   = _NoNodes[_NoLevels - 1];

  // element  definition
  for ( int  inode = 0; inode < el_nnodes; inode++ )    {
      //get the global node number
      el_conn[inode] = _el_map[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * iproc] ) * el_nnodes + inode];

      // get the element coordinades
      for ( int  idim = 0; idim < _dim; idim++ ) {
          xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
          }
      }

  return;
  }

// ===============================================================
/// This function returns the element connectivity and coordinates
// ==============================================================
void MGMesh::get_el_nod_conn (
  const int  bdry,  // bd(1)/vol(0) flag    (<- input)=(1)/(0) flag
  const int  Level, // level                (<- input)
  const int  iel,   // element              (<- input)
  int  el_conn[],   // element connectivity (-> output)
  double xx[]       // element coordinates  (-> output)
) const {// =================================

  const int  el_nnodes = _GeomEl.n_q[bdry]; // element nodes
  const int  offset   = _NoNodes[_NoLevels - 1];

  // element  definition
  for ( int  inode = 0; inode < el_nnodes; inode++ )    {
      //get the global node number
      el_conn[inode] = _el_map[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * _iproc] ) * el_nnodes + inode];

      // get the element coordinates
      for ( int  idim = 0; idim < _dim; idim++ ) {
          xx[inode + idim * el_nnodes] = _xyz[el_conn[inode] + idim * offset];
          }
      }

  return;
  }

// ===============================================================
/// This function return the element neigbours
// ==============================================================
void MGMesh::get_el_neighbor (
  const int  el_sides, // element sides        (<- input) _GeomEl._n_sides[bdry];
  const int  bdry,     // bd(1)/vol(0) flag    (<- input) (1)/(0) flag
  const int  Level,    // level                (<- input)
  const int  iel,      // element              (<- input)
  int el_neigh[],      // element connectivity (-> output)
  int iproc
) const {// =================================

//   const int  el_sides= _GeomEl._n_sides[bdry]; // element nodes
  for ( int  iside = 0; iside < el_sides; iside++ )    {
      //get the global node number
      el_neigh[iside] = _el_neighbor[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * iproc] ) * el_sides + iside];
      }

  return;
  }

// ===============================================================
/// This function return the element neigbours
// ==============================================================
void MGMesh::get_el_neighbor (
  const int  el_sides, // element sides        (<- input) _GeomEl._n_sides[bdry];
  const int  bdry,     // bd(1)/vol(0) flag    (<- input) (1)/(0) flag
  const int  Level,    // level                (<- input)
  const int  iel,      // element              (<- input)
  int el_neigh[]       // element connectivity (-> output)
) const {// =================================

//   const int  el_sides= _GeomEl._n_sides[bdry]; // element nodes
  for ( int  iside = 0; iside < el_sides; iside++ )    {
      //get the global node number
      el_neigh[iside] = _el_neighbor[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * _iproc] ) * el_sides + iside];
      }

  return;
  }

// ===============================================
/// This function gets the element connectivity
void MGMesh::get_el_conn (
  const int  el_nnodes, // element nodes    (<- input) = _GeomEl.n_q[bdry];
  const int  bdry,      // vol/surf         (<- input)=(1)/(0) flag
  const int  Level,     // Level            (<- input)
  const int  iel,       // element          (<- input)
  int  el_conn[]        // connectivity map (-> output)
) const {// ==========

  //get the global node number
  for ( int  n = 0; n < el_nnodes; n++ )    {
      el_conn[n] = _el_map[bdry][ ( iel + _off_el[bdry][Level + _NoLevels * _iproc] ) * el_nnodes + n];
      }

  return;
  }

// ===============================================
/// identity
void MGMesh::idFEMord (
  const  int n
) {
  for ( int  ik = 0; ik < n; ik++ ) _ord_FEM[ik] = ik;
  }

// =========================================
/// This function returns the coordinate vector
void MGMesh::xcoord (
  double cvect[],
  const int  n_nodes,
  const int  offset
) const { // =========================================
  // computation xyz vector in variable order
  for ( int  i = 0; i < n_nodes; i++ )
    cvect[i] = _xyz[i + offset];
  }

// =========================================
void MGMesh::nodesxyz (
  double xyzvect[],  // coord vector ->
  const int  n_nodes // node number  <-
) const {// ====================================
  // computation xyz vector in node order
  for ( int  i = 0; i < n_nodes; i++ )    {
      xyzvect[i * 3] = _xyz[i];
      xyzvect[i * 3 + 1] = 0.;
      xyzvect[i * 3 + 2] = 0.;
#if DIMENSION==2
      xyzvect[i * 3 + 1] = _xyz[i + n_nodes];
#endif
#if DIMENSION==3
      xyzvect[i * 3 + 1] = _xyz[i + n_nodes];
      xyzvect[i * 3 + 2] = _xyz[i + 2 * n_nodes];
#endif

      }
  }

// =============================================
/// This function computes the subconnectivity
int  MGMesh::sub_conn (
  int  gl_conn[],        // global sub-connectivity map
  const int  ifem,       // fem type
  const int  ilev,       // Level
  const int  n_points_el // point number
) const {// =============================================

  int  n_elements = _NoElements[ifem][ilev];
  int  type_family_in = _type_FEM[ilev + ifem * _NoLevels];

  int  icount = 0;

  for ( int  ik = 0; ik < n_elements; ik++ )  {
      for ( int  inode = 0; inode < n_points_el; inode++ )  {
          gl_conn[icount] = _el_map[ifem][inode + ik * type_family_in];
          icount++;
          }
      }

#ifdef PRINT_INFO
  std::cout << icount << std::endl;
#endif
  return icount;

  }

// =============================================
/// This function computes the connectivity
void MGMesh::conn (
  int  gl_conn[],       // global sub-connectivity map
  const int  ifem,      // fem type
  const int  indx_mesh  // mesh
) const {// =============================================
  for ( int  ik = 0; ik < _NoElements[ifem][indx_mesh]*_type_FEM[indx_mesh]; ik++ )  {
      gl_conn[ik] = _el_map[ifem][ik];
      }

  return;
  }








// ===============================================
/// This function manages the printing in Xdmf format
void MGMesh::print (
  const int  Level, // level
  const int  t_step
)  { // ============================================

  const int  iproc = _iproc;

  if ( iproc == 0 ) {
      print_conn_lin_hf5 ( Level ); // print new connection mesh in hdf5 format
      print_vol ( Level, t_step ); // print new connection mesh in xdmf format (vol+bd)
      }

  return;
  }

// ==============================================
/// This funcTion prints the volume/boundary
/// Mesh (connectivity) in Xdmf format
void MGMesh::print_vol (
  const int  Level, // level
  const int  /*t_step */ //time
)  {// =========================================

  std::string    inout_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file ( "BASEMESH" );
  std::string    connlin = _mgutils.get_file ( "CONNLIN" );

  std::ostringstream namefile;
  namefile << inout_dir << basemesh <<  ".xmf";
  std::ostringstream conn_file;
  conn_file/* << inout_dir << "/"*/ << basemesh;
  std::ostringstream topol_file;
  topol_file << /*inout_dir << */ basemesh << connlin << ".h5";
  conn_file << ".h5";

  std::ofstream out ( namefile.str().c_str() );
//   int nvrt[2];  std::string mtype[2];
  std::string grid_mesh[2];
  grid_mesh[0] = "Mesh";
  grid_mesh[1] = "Boundary";

  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
  out << " [ <!ENTITY HeavyData \"mesh.h5 \"> ] ";
  out << ">\n";
  out << " \n";
  out << "<Xdmf> \n" << "<Domain> \n";

  for ( int  ik = 0; ik < _NoFamFEM; ik++ ) {
      out << "<Grid Name=\"" << grid_mesh[ik].c_str() << "\"> \n";
      // Topology
      print_xmf_topology ( out, topol_file.str(), Level, ik );
      // Geometry
      print_xmf_geometry ( out, conn_file.str(), Level, ik );
      out << "</Grid> \n";
      }

  out << "</Domain> \n" << "</Xdmf> \n";
  out.close();
  return;
  }
// ==============================================
/// This funcTion prints the Xdmf format mesh topology
void MGMesh::print_xmf_topology (
  std::ofstream & out, //  file xdmf
  std::string store_file, // file where the topology is
  int   Level,           // Level
  int  ik              // Fem family (boundary (1) /volume (0) mesh)
)  {// =========================================
  // Topology
  int tot_el = _NoElements[ik][Level] * ( _GeomEl.n_se[ik] );
  out << "<Topology Type=\"" << _GeomEl.pname[ik]
      << "\"  Dimensions=\"" << tot_el
      << "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\" "
      <<  tot_el << "  " << _GeomEl.n_l[ik]
      << "\" Format=\"HDF\">  \n";
  out << store_file    << ":MSH" << ik  << "CONN \n";
  out << "</DataStructure> \n" << "</Topology> \n";
  return;
  }

// ==============================================
/// This function prints mesh Xdmf format geometry
void MGMesh::print_xmf_geometry (
  std::ofstream & out, //  file xdmf
  std::string store_file, // file where the geometry is
  int   Level,           // Level
  int  /*ik*/              // Fem family (boundary (1) /volume (0) mesh)
)  {// =========================================
  // Geometry
  out << "<Geometry Type=\"X_Y_Z\"> \n";

  for ( int  ix = 1; ix < 4; ix++ ) {
      out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
          <<  _NoNodes[Level] << "  " << 1
          <<  "\" Format=\"HDF\">  \n";
      out << store_file  << ":/NODES/COORD/X" << ix << "\n";
      out << "</DataStructure> \n";
      }

  out << " </Geometry>\n";
  return;
  }




// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGMesh::Translate ( const int dir, NumericVectorM & x_old ) { // ------------------------------------


  const int flag_moving_mesh = ( int ) _mgutils._geometry["moving_mesh"];

  if ( flag_moving_mesh ) {
// #ifdef FINE_MOVEMENT
      /// E) mesh update
      double ww = 1.;
      const int n_nodes = _NoNodes[_NoLevels - 1];
      const int offsetp = dir * n_nodes;


      for ( int inode = 0; inode < n_nodes; inode++ ) {
          _xyz[inode + offsetp] = _xyzo[inode + offsetp] + ww * x_old ( inode );
          }

// #endif
// #ifdef COARSE_MOVEMENT
//     MoveMesh(_NoLevels-1);
// // disp[_NoLevels-1]->zero();
//     for(int inode=0; inode<n_nodes; inode++) {
//       _mgmesh._xyz[inode+offsetp] += (*disp[_NoLevels-1])(inode);
//       _mgmesh._dxdydz[inode+offsetp]= (*disp[_NoLevels-1])(inode);
//     }
//     MoveMesh(_NoLevels-1);
// #endif
      }

  // ==============================================================
  return;
  }



// ////NEW MOVEMESH
// #include "dense_vectorM.h"   // algebra dense vectors
// //// =========================================================================================
// //// This function controls the assembly and the solution of the T_equation system:
// void MGMesh::MoveMesh(const int dir, NumericVectorM &x_old, int nvars[], int n_fe[]) { // ------------------------------------
//
//
//     const int flag_moving_mesh = (int)_mgutils._geometry["moving_mesh"];
//
//     const int n_nodes=_NoNodes[_NoLevels-1];
//     const int n_elem=_NoElements[0][_NoLevels-1];
//     const int offsetp=dir*n_nodes;
//     const int  offset = _NoNodes[_NoLevels-1];                     // mesh nodes
//     double _xx_qnds[NDOF_FEM *DIMENSION];
//     int   bc_vol[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags (vol int)
//     int   bc_bd[NDOF_FEM* ( DIMENSION+1 )]; ///<  element  b.cond flags  (bd int)
//     int        el_conn[NDOF_FEM];                     // element connectivity
//     int el_mat_nrows =0;
//     int el_ndof[3];
//     el_ndof[0]=0;
//     const int  el_sides= _GeomEl._n_sides[0];                      // element nodes
//     const int Level=_NoLevels-1;
//     int        sur_toply[NDOF_FEMB]; // boundary topology
//     const ParallelM::Communicator &comm1=_comm.comm();
// //     int n_glob=*_NoNodes;
//     std::vector<NumericVectorM *> disp; ///< displacement solution x
//     disp.resize(_NoLevels);// displacement for mesh
//     disp[Level] = NumericVectorM::build(comm1).release();
//     disp[Level]->init(offset,false, SERIALM);
//     disp[Level]->zero();
// //     disp[0] = NumericVectorM::build(comm1).release();
// //     disp[0]->init(offset,false, SERIALM);
// //     disp[0]->zero();
//     for(int ideg=1; ideg<3; ideg++) {                               //     ...
//         el_ndof[ideg]=((nvars[ideg]>0)?    n_fe[ideg]:0);                    //   computing
//         el_mat_nrows +=nvars[ideg]*el_ndof[ideg];
//     };
//
//     int el_mat_ncols = el_mat_nrows;                     //square matrixbc_Neum[0]=0;
//     std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector
//     MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(this); //??????
//     int Level_c=0;//_NoLevels-1;
//     DenseVectorM LoD;//Local Displacement
//     LoD.resize(el_ndof[2]);
//
//
//
//
//
//
//
//
//
//      // print linear -----------------------------------
//   double *sol_c=new double[_GeomEl.n_l[0]];
// //   double *Prol=_mgmesh._GeomEl.Prol;
// //   const int  lev_l=(Level >0) ?Level-1:_NoLevels;
//
// //   for(int ivar=_nvars[2]; ivar<_nvars[2]+_nvars[1]; ivar++)   {
// //     std::string var_name = _var_names[ivar];
//     // class writing routine
//     // 2bA projection of the fine nodes -------------------------
// //     for (int  isubdom=0; isubdom<_mgmesh._n_subdom; isubdom++) {
// //       for ( int i= _mgmesh._off_nd[isubdom*_NoLevels];
// //            i< _mgmesh._off_nd[isubdom*_NoLevels]+
// //            _mgmesh._off_nd1[Level+1+isubdom*_NoLevels]-
// //            _mgmesh._off_nd1[isubdom*_NoLevels];i++) {
// //         sol[i]=(*x_old[Level])(_node_dof[Level][i+ivar*offset]);
// //       }
// //     } // 2bA end proj fine grid ----------------------------------
//   int ndof_lev=0;
// //         for(int pr=0; pr <proc; pr++) {
// //             int delta =_off_el[0][pr*_NoLevels+Level_c+1]-_off_el[0][pr*_NoLevels+Level_c];
// //             ndof_lev +=delta;
// //         }
//     //  2bB element interpolation over the fine mesh -----------------------
//     for(int iproc=0; iproc<_n_subdom; iproc++) {
//
//        int n_iel0=_off_el[0][_NoLevels*iproc+0];
//       int n_ielf=_off_el[0][0+_NoLevels*iproc+1];
//       for(int iel=0; iel <n_ielf-n_iel0; iel++) {
//          int e_indx=(iel+n_iel0)*NDOF_FEM;
//
//        const int ivar =0;
//         for(int in=0; in<NDOF_P; in++) {
//           int gl_i=_el_map[0][ e_indx+in];
//           double val= x_old(/*_node_dof[Level][*/gl_i+ivar*offset/*]*/);
//           sol_c[in]= val;
//         }
//
//         for(int in=0; in<NDOF_FEM; in++) {    // mid-points
//           double sum=0.;
//           for(int jn=0; jn<NDOF_P; jn++) { sum += _GeomEl.Prol[in*NDOF_P+jn]*sol_c[jn]; }
//           disp[Level]->set(/*_node_dof[Level][*/ _el_map[0][ e_indx+in]/*]*/,sum);
//         }
//       }
//
//       for (int ilev=1;ilev<=Level;ilev++){
//       int n_iel0=_off_el[0][_NoLevels*iproc+ilev];
//       int n_ielf=_off_el[0][ilev+_NoLevels*iproc+1];
//       for(int iel=0; iel <n_ielf-n_iel0; iel++) {
//          int e_indx=(iel+n_iel0)*NDOF_FEM;
//         // element iel ----------------------------
//         // vertices
// //         for (int  in=0; in<NDOF_P; in++) sol_c[in]= sol[_mgmesh._el_map[0][indx+in]]*_refvalue[ivar];
//         const int ivar =0;
//         for(int in=0; in<NDOF_P; in++) {
//           int gl_i=_el_map[0][ e_indx+in];
//           double val= (*disp[Level])(/*_node_dof[Level][*/gl_i+ivar*offset/*]*/);
//           sol_c[in]= val;
//         }
//
//         for(int in=0; in<NDOF_FEM; in++) {    // mid-points
//           double sum=0.;
//           for(int jn=0; jn<NDOF_P; jn++) { sum += _GeomEl.Prol[in*NDOF_P+jn]*sol_c[jn]; }
//           disp[Level]->set(/*_node_dof[Level][*/ _el_map[0][ e_indx+in]/*]*/,sum);
//         }
//       } // ---- end iel -------------------------
//     } // 2bB end interpolation over the fine mesh ------------------------
// //     _mgutils.print_Dhdf5(file_id,var_name,dimsf,sol);
// //   } // ivar
//     }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// //     for(int proc=0; proc<_n_subdom; proc++) {
// //         int ndof_lev=0;
// //         for(int pr=0; pr <proc; pr++) {
// //             int delta =_off_el[0][pr*_NoLevels+Level_c+1]-_off_el[0][pr*_NoLevels+Level_c];
// //             ndof_lev +=delta;
// //         }
// //
// //
// //
// //         //find the coarse nodes
// //         const int coarse_nel_e =_off_el[0][Level_c+_NoLevels*proc+1]; // start element
// //         const int coarse_nel_b =_off_el[0][Level_c+_NoLevels*proc];   // stop element
// //         int       coarse_node[NDOF_FEM];//, elb_conn[NDOF_FEMB];                     // element connectivity
// //         std::cout;
// //
// //         for(int iel=0; iel < (coarse_nel_e - coarse_nel_b); iel++) {
// //             LoD.zero();
// //             //getting connectivity and coordinates
// //             get_el_nod_conn(0,Level_c,iel,el_conn,_xx_qnds,proc);
// //
// //             for(int id=0; id<NDOF_FEM; id++)  {
// //                 // quadratic -------------------------------------------------
// //                 if(id <el_ndof[2])  for(int  ivar=0; ivar<nvars[2]; ivar++) {     //ivarq is like idim
// //                         const int indx_loc_ql = id +ivar*el_ndof[2];
// // //         const int indx_glob= el_conn[id]+ivar*offset;
// //
// //                         el_dof_indices[indx_loc_ql]=el_conn[id]+ivar*offset;     //from mesh to dof
// //                     } // end quadratic ------------------------------------------------
// //             }
// //
// //
// //             int count=0;
// //             int written[NDOF_FEM]; //needed to have fine displacement on the solid
// //             for(int inode=0; inode<el_ndof[2]; inode++)     {
// //                 coarse_node[inode]=  _node_map[Level_c][el_dof_indices[inode]];
// //                 int glob_indx=_node_map[Level_c][el_dof_indices[inode]];
// //                 if(ext_mesh->_mat_id[iel]==4) {
// //                     LoD(inode) =x_old(glob_indx);
// //                     written[count]=inode;
// //                     count++;
// //                 }//end solid
// //                 else { //liquid
// //                     if(inode < NDOF_P)
// //                     {
// //                         LoD(inode)= x_old(glob_indx);
// //                     }
// //                 }
// //             } //end inode
// //
// // #if DIMENSION==2
// //             for(int iside=0; iside< el_sides; iside++)  {
// //                 int check=0;
// //                 for(int idof=0; idof<NDOF_FEMB; idof++) {
// //                     sur_toply[idof]=_GeomEl._surf_top[idof+NDOF_FEMB*iside];    // local nodes
// //                 }
// //                 check=0;
// //                 for(int i=0; i<count; i++) if(written[i]==sur_toply[NDOF_FEMB-1]) {
// //                         check =1;
// //                     }
// //                 if(check==0) {
// //                     LoD(sur_toply[NDOF_FEMB-1]) = 0.5*(LoD(sur_toply[NDOF_FEMB-3]))
// //                                                   + 0.5*(LoD(sur_toply[NDOF_FEMB-2]));
// //                 }//end if-check
// //             }//end iside
// // #endif
// // #if DIMENSION==3
// //             for(int iside=0; iside< el_sides; iside++)  {
// // //                 int check=0;
// //                 for(int idof=0; idof<NDOF_FEMB; idof++) {
// //                     sur_toply[idof]=_GeomEl._surf_top[idof+NDOF_FEMB*iside];    // local nodes
// //                 }
// //                if (ext_mesh->_mat_id[iel]==2){
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
// // //                 if (check==0) {
// // //                 if(_bc_vol[sur_toply[4]]<1.5)
// //       LoD(sur_toply[4]) = 0.5*(LoD(sur_toply[0]))
// //                             + 0.5*(LoD(sur_toply[1]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
// // //                 if (check==0) {
// // //                 if(_bc_vol[sur_toply[5]]<1.5)
// //       LoD(sur_toply[5]) = 0.5*(LoD(sur_toply[1]))
// //                             + 0.5*(LoD(sur_toply[2]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
// // //                 if (check==0) {
// // //                 if(_bc_vol[sur_toply[6]]<1.5)
// //       LoD(sur_toply[6]) = 0.5*(LoD(sur_toply[2]))
// //                             + 0.5*(LoD(sur_toply[3]));
// // //                 if(_bc_vol[sur_toply[7]]<1.5)
// //       LoD(sur_toply[7]) = 0.5*(LoD(sur_toply[3]))
// //                             + 0.5*(LoD(sur_toply[0]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
// // //                 if (check==0) {
// // //                 if(_bc_vol[sur_toply[8]]<1.5)
// //       LoD(sur_toply[8]) = 0.25*(LoD(sur_toply[0]))
// //                             +0.25*(LoD(sur_toply[1]))
// //                             + 0.25*(LoD(sur_toply[2]))
// //                             + 0.25*(LoD(sur_toply[3]));
// // //                 }
// //          }
// //             }//end iside
// // #endif
// //             int check=0;
// //             for(int i=0; i<count; i++) if(written[i]==NDOF_FEM-1) {
// //                     check =1;    //last node
// //                 }
// //             if(check==0) {
// //                 LoD(NDOF_FEM-1)=0;
// //                 for(int jnode=0; jnode<NDOF_P; jnode++) {
// //                     LoD(NDOF_FEM-1)+=LoD(jnode)/NDOF_P;
// //                 }
// //             }
// //             for(int inode=0; inode < el_ndof[2] ; inode++) {
// //                 disp[Level]->set(_node_map[Level_c][el_dof_indices[inode]],LoD(inode));    //setting the global displacement
// //             }
// //         } //end iel
// //         el_dof_indices.clear();
// //     } //end proc
// //     //      end coarse element loop
// // //        end find the coarse nodes
// //     //finer level
// //     for(int ilev=Level_c+1; ilev<_NoLevels; ilev++) {
// //         for(int proc=0; proc<_n_subdom; proc++) {
// //             int ndof_lev=0;
// //             for(int pr=0; pr <proc; pr++) {
// //                 int delta =_off_el[0][pr*_NoLevels+ilev+1]-_off_el[0][pr*_NoLevels+ilev];
// //                 ndof_lev +=delta;
// //             } //end pr
// //             const int nel_e =_off_el[0][ilev+_NoLevels*proc+1]; // start element
// //             const int nel_b =_off_el[0][ilev+_NoLevels*proc];   // stop element
// //             for(int iel=0; iel < (nel_e - nel_b); iel++) {
// //                 LoD.zero();
// //                 get_el_nod_conn(0,ilev,iel,el_conn,_xx_qnds,proc);
// //
// //                 for(int id=0; id<NDOF_FEM; id++)  {
// //                     // quadratic -------------------------------------------------
// //                     if(id <el_ndof[2])  for(int  ivar=0; ivar<nvars[2]; ivar++) {     //ivarq is like idim
// //                             const int indx_loc_ql = id +ivar*el_ndof[2];
// //                             el_dof_indices[indx_loc_ql]=_node_dof[ilev][el_conn[id]+ivar*offset];     //from mesh to dof
// //                         } // end quadratic ------------------------------------------------
// //                 }
// //                 int count=0;
// //                 int written[NDOF_FEM]; //needed to have fine displacement on the solid
// //                 int countg=0;
// //                 for(int inode=0; inode<el_ndof[2]; inode++)     {
// //                     int dof_idx=_node_map[ilev][el_dof_indices[inode]];
// //                     if(ext_mesh->_mat_id[iel+nel_b]==4) {
// //                         written[count]=inode;
// //                         count++;
// //                         LoD(inode) =x_old(dof_idx);
// //                     }//end solid
// //                     else { //liquid
// //                         if(inode < NDOF_P) {
// //                             LoD(inode)=(*disp[Level])(dof_idx);
// //                         }
// //                     }//end liquid
// //                 }//end loop inode
// // #if DIMENSION==2
// //                 for(int iside=0; iside< el_sides; iside++)  {
// //                     int check=0;
// //                     for(int idof=0; idof<NDOF_FEMB; idof++) {
// //                         sur_toply[idof]=_GeomEl._surf_top[idof+NDOF_FEMB*iside];    // local nodes
// //                     }
// //                     check=0;
// //                     for(int i=0; i<count; i++) if(written[i]==sur_toply[NDOF_FEMB-1]) {
// //                             check =1;
// //                         }
// //                     if(check==0) {
// //                         LoD(sur_toply[NDOF_FEMB-1]) = 0.5*(LoD(sur_toply[NDOF_FEMB-3]))
// //                                                       + 0.5*(LoD(sur_toply[NDOF_FEMB-2]));
// //                     }
// //
// //                 }
// // #endif
// // #if DIMENSION==3
// //                 for(int iside=0; iside< el_sides; iside++)  {
// // //                 int check=0;
// //                     for(int idof=0; idof<NDOF_FEMB; idof++) {
// //                         sur_toply[idof]=_GeomEl._surf_top[idof+NDOF_FEMB*iside];    // local nodesdet
// //                     }
// //                     if(ext_mesh->_mat_id[iel+nel_b]==4) {
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
// // //                 if (check==0) {
// // //                     if(_bc_vol[sur_toply[4]]<1.5)
// //           LoD(sur_toply[4]) = 0.5*(LoD(sur_toply[0]))
// //                                 + 0.5*(LoD(sur_toply[1]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
// // //                 if (check==0) {
// // //                     if(_bc_vol[sur_toply[5]]<1.5)
// //           LoD(sur_toply[5]) = 0.5*(LoD(sur_toply[1]))
// //                                 + 0.5*(LoD(sur_toply[2]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
// // //                 if (check==0) {
// // //                     if(_bc_vol[sur_toply[6]]<1.5)
// //           LoD(sur_toply[6]) = 0.5*(LoD(sur_toply[2]))
// //                                 + 0.5*(LoD(sur_toply[3]));
// // //                     if(_bc_vol[sur_toply[7]]<1.5)
// //           LoD(sur_toply[7]) = 0.5*(LoD(sur_toply[3]))
// //                                 + 0.5*(LoD(sur_toply[0]));
// // //                 }
// // //                 check=0;
// // //                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
// // //                 if (check==0) {
// // //                     if(_bc_vol[sur_toply[8]]<1.5)
// //           LoD(sur_toply[8]) = 0.25*(LoD(sur_toply[0]))
// //                                 +0.25*(LoD(sur_toply[1]))
// //                                 + 0.25*(LoD(sur_toply[2]))
// //                                 + 0.25*(LoD(sur_toply[3]));
// //                 }
// // det
// //                 }//end iside
// // #endif
// //                 int check=0;
// //                 for(int i=0; i<count; i++) if(written[i]==NDOF_FEM-1) {
// //                         check =1;    //last node
// //                     }
// //                 if(check==0) {
// //                     for(int jnode=0; jnode<NDOF_P; jnode++) {
// //                         LoD(NDOF_FEM-1)+=LoD(jnode)/NDOF_P;
// //                     }
// //                 }
// //                 int check2=0;
// //                 for(int inode=0; inode < el_ndof[2] ; inode++) {
// //                     disp[Level]->set(_node_map[ilev][el_dof_indices[inode]],LoD(inode));
// //                 }
// //             }//end iel
// //
// //             el_dof_indices.clear();
// //         } //end iproc
// //
// //     }//end ilev
//     disp[Level]->close();
//
//     if(flag_moving_mesh) {
//         double normttt=x_old.l2_norm();
//         std::cout << "  ---------------------------- kkkkk  " << normttt << " \n" ;
//         double ww=1.;    //under relaxation
//         double sum1 =0.;
//         double sum2 =0.;double sum3 =0.;double sum4 =0.;
//         const int n_nodes=_NoNodes[Level];  //nodes at fine level
//         const int offsetp=dir*n_nodes;
//         for(int inode=0; inode<n_nodes; inode++) {
//
//             _xyz[inode+offsetp] =_xyzo[inode+offsetp]+ww*(*disp[Level])(inode);   //mesh update  SE METTO X_OLD INVECE CHE DISP
//
//                      //ALLORA FUNZIONA DOPO CHE SI SPACCA
//             sum1 += _xyz[inode+offsetp]*_xyz[inode+offsetp];
//             sum2 +=_xyzo[inode+offsetp]*_xyzo[inode+offsetp];
//              sum3 +=x_old(inode)*x_old(inode);
//              sum4 +=(*disp[Level])(inode)*(*disp[Level])(inode);
// //       std::cout<<(*disp[_NoLevels-1])(inode)<<"   <-disp x_old->   "<<x_old(inode)<<"\n";
// //       if (inode%5==0) std::cout<<std::endl;
//         }
//             std::cout << "  ---------------------------- kkkkk xxxx  " << sqrt(sum1) << " " << sqrt(sum2) << " " << sqrt(sum3)<<" " <<sqrt(sum4)  <<   " \n" ;
//     }
// //     disp[Level]->zero();
//     // ==============================================================
//     return;
// }




// ===============================================
//           HDF5 FUNCTION
// ===============================================

/// Read mesh from hdf5 file (namefile)
///           as MGMesh class (MGMesh.h):
/// data parameters : DFL,_type_FEM
/// Verteces        : _NoNodes,_xyz
/// Connectivity    : _NoElements,_conn_map

void MGMesh::read_c() {  // ======================

//   const int  restart = _mgutils.get_par("restart");
  const int  restart = stoi ( _mgutils._sim_config["restart"] );

  if ( restart == 0 )  { //=================================
      std::string    input_dir = _mgutils._inout_dir;
      std::string    basemesh = _mgutils.get_file ( "BASEMESH" );

      // Open an existing file. ---------------
      std::ostringstream meshname;
      meshname << input_dir << basemesh << ".h5";
      std::cout << " Reading mesh from= " <<  meshname.str() <<  std::endl;
      hid_t  file_id = H5Fopen ( meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      hid_t  status = 0;
      // Reading DFL -------------------------
      int  topdata[4];
      status = H5Dread ( H5Dopen ( file_id, "/DFLS"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata );
      _dim = topdata[0];        //dimension
      _NoFamFEM = topdata[1];   //Volume/Boundary
      _NoLevels = topdata[2];   //number of levels
      _n_subdom = topdata[3];   //number of subdomains

      // Reading _type_FEM ------------
      int  n_meshes = _NoFamFEM * _NoLevels;
      _type_FEM = new int [n_meshes];
      status = H5Dread ( H5Dopen ( file_id, "/FEM"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM );
#ifdef PRINT_INFO
      std::cout << " NoFamFEM= " << _NoFamFEM << " NoLevels= " << _NoLevels << std::endl;
#endif
      // Reading _NoNodes ------------------
      _NoNodes = new int [n_meshes + 1];
      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/NDxLEV"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoNodes );
      // Reading  _xyz ----------------
      double ILref = 1. / _Lref;
      int  n_nodes = _NoNodes[_NoLevels - 1];
      _xyz = new double[_dim * n_nodes];        //double scale[3]
      _xyzo = new double[_dim * n_nodes];        //double scale[3]
      _dxdydz = new double[_dim * n_nodes];
      double * coord;
      coord = new double[n_nodes];

      for ( int  inode = 0; inode < DIMENSION * n_nodes; inode++ ) {
          _xyz[inode] = 0.; //nondimensionalization
          _xyzo[inode] = 0.; //nondimensionalization
          _dxdydz[inode] = 0.; //nondimensionalization
          }

      std::cout << "-----------------------------------------------\n";
      std::cout << "-----------------------------------------------\n";
      std::cout << "        DIMENSION   " << DIMENSION << std::endl;
      std::cout << "-----------------------------------------------\n";
      std::cout << "-----------------------------------------------\n";
      double para_rot = -6.;

      for ( int  kc = 0; kc < _dim; kc++ ) {
          std::ostringstream Name;
          Name << "NODES/COORD/X" << kc + 1;
          status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                       , H5P_DEFAULT
#endif
                                     ), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );



          for ( int  inode = 0; inode < n_nodes; inode++ ) {

//         _xyz[inode+0*n_nodes] +=(1-kc)*cos(-3.141592/para_rot)*coord[inode]*ILref
//                                   -(kc)*sin(-3.141592/para_rot)*coord[inode]*ILref; //nondimensionalization
//          _xyzo[inode+0*n_nodes] +=(1-kc)*cos(-3.141592/para_rot)*coord[inode]*ILref
//                                   -(kc)*sin(-3.141592/para_rot)*coord[inode]*ILref; //nondimensionalization
//            _xyz[inode+1*n_nodes] +=(1-kc)*sin(-3.141592/para_rot)*coord[inode]*ILref
//                                   +(kc)*cos(-3.141592/para_rot)*coord[inode]*ILref; //nondimensionalization
//          _xyzo[inode+1*n_nodes] +=(1-kc)*sin(-3.141592/para_rot)*coord[inode]*ILref
//                                   +(kc)*cos(-3.141592/para_rot)*coord[inode]*ILref; //nondimensionalization
//
//           _dxdydz[inode+kc*n_nodes]=0.; //nondimensionalization

              _xyz[inode + kc * n_nodes] = coord[inode] * ILref; //nondimensionalization
              _xyzo[inode + kc * n_nodes] = coord[inode] * ILref; //nondimensionalization
              _dxdydz[inode + kc * n_nodes] = 0.; //nondimensionalization
              }
          }

      delete []coord;
      std::cout << " Reading Multimesh with  " << n_meshes << " meshes, " << n_nodes <<
                " nodes and " << DIMENSION * n_nodes << " coordinates" << std::endl;
      _off_nd[0] = new int[_n_subdom * _NoLevels + 1];

      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/OFF_ND"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_nd[0] );
      _off_nd[1] = new int[_n_subdom * _NoLevels + 1];
      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/OFF_ND1"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_nd[1] );
      // Reading  _NoElements --------------
      _off_el = new int * [_NoFamFEM];             //volume and boundary
      _off_el[0] = new int [_n_subdom * _NoLevels + 1]; //offset for VOLUME
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/OFF_EL"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[0] );
      _off_el[1] = new int [_n_subdom * _NoLevels + 1]; //offset for boundary
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM1/OFF_EL"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[1] );
      // Number of elements (for each level)
      _NoElements = new int * [_NoFamFEM]; //quadratic and linear

      for ( int  ifem = 0; ifem < _NoFamFEM; ifem++ ) {
          _NoElements[ifem] = new int [_NoLevels];
          std::ostringstream Name;
          Name << "/ELEMS/FEM" << ifem  << "/NExLEV";
          status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()

#if HDF5_VERSIONM != 1808
                                       , H5P_DEFAULT
#endif
                                     ),
                             H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoElements[ifem] );
          }

      // node mapping ------------------------------------------------------------

      _node_map = new int * [_NoLevels + 1];      // _NoLevels+1 inverse node maps
      int  n_nodes_top = _NoNodes[_NoLevels - 1]; //  same length in file (n_nodes_top)
      int * temp = new int [n_nodes_top];         // temp to read

      for ( int  ilev = 0; ilev < _NoLevels; ilev++ ) { // maps from o to 2*_NoLevels-1
          _node_map[ilev] = new int  [_NoNodes[ilev]]; // inverse node map

          // hdf5 reading
          std::ostringstream Name;
          Name << "/NODES/MAP/MAP" << ilev;          // filename
          status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                       , H5P_DEFAULT
#endif

                                     ),  // read from file
                             H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp );

          // inverse map
          for ( int  inode = 0; inode < n_nodes_top; inode++ )
            if ( temp[inode] != -1 ) _node_map[ilev][temp[inode]] = inode;
          }

      // linear coarse map (storage in [2*_NoLevels])
      _node_map[_NoLevels] = new int  [_NoNodes[2 * _NoLevels]];
      // hdf5 reading
      std::ostringstream Name;
      Name << "/NODES/MAP/MAP" << _NoLevels;     // filename
      status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp ); // read from file

      // inverse map
      for ( int  inode = 0; inode < n_nodes_top; inode++ )
        if ( temp[inode] != -1 ) _node_map[_NoLevels][temp[inode]] = inode;

      //connectivity map -----------------------------------------------------
      _el_map = new int * [_NoFamFEM];
      _el_map[0] = new int  [_off_el[0][_NoLevels * _n_subdom]*NDOF_FEM]; // Volume
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/MSH"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),// Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[0] );

      _el_map[1] = new int  [_off_el[1][_NoLevels * _n_subdom]*NDOF_FEMB]; // Boundary
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM1/MSH"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), // Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[1] );


      // neighbour element map -----------------------------------------------------
//     const int n_sides=_GeomEl._n_sides[0];
      _el_neighbor = new int * [_NoFamFEM - 1]; //now only for the volume family
      _el_neighbor[0] = new int [_off_el[0][_NoLevels * _n_subdom]*NDOF_FEM]; // Volume
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/EL_NEIG"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),// Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_neighbor[0] );

      // Close the file. ---------------------
      H5Fclose ( file_id );
      delete []temp;

      if ( status < 0 ) {
          std::cout << " File Mesh input in data_in is missing";
          abort();
          }

      }


  else {  //======restart != 0 ======================

      std::string    inout_dir = _mgutils._inout_dir;
      std::string    basemesh = _mgutils.get_file ( "BASEMESH" );
      const int ndigits = stoi ( _mgutils._sim_config["ndigits"] );
      std::ostringstream meshname_xyz;
      meshname_xyz <<  inout_dir << basemesh << "." << std::setw ( ndigits ) << std::setfill ( '0' ) << 0 << ".h5";
      std::ostringstream meshname;
      meshname << inout_dir << basemesh << ".h5";
      std::cout << " Reading mesh from= " <<  meshname.str() <<  std::endl;
      std::cout << " Reading mesh coordinates from = " <<  meshname_xyz.str() <<  std::endl;
      hid_t  file_id     = H5Fopen ( meshname.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      hid_t  file_id_xyz = H5Fopen ( meshname_xyz.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      hid_t  status = 0;
      // Reading DFL -------------------------
      int  topdata[4];
      status = H5Dread ( H5Dopen ( file_id, "/DFLS"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif

                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata );
      _dim = topdata[0];        //dimension
      _NoFamFEM = topdata[1];   //Volume/Boundary
      _NoLevels = topdata[2];   //number of levels
      _n_subdom = topdata[3];   //number of subdomains

      // Reading _type_FEM ------------
      int  n_meshes = _NoFamFEM * _NoLevels;
      _type_FEM = new int [n_meshes];
      status = H5Dread ( H5Dopen ( file_id, "/FEM"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM );
#ifdef PRINT_INFO
      std::cout << " NoFamFEM= " << _NoFamFEM << " NoLevels= " << _NoLevels << std::endl;
#endif
      // Reading _NoNodes ------------------
      _NoNodes = new int [n_meshes + 1];
      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/NDxLEV"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoNodes );

      // Reading  _xyz ----------------
      double ILref = 1. / _Lref;
      int  n_nodes = _NoNodes[_NoLevels - 1];
      _xyz = new double[_dim * n_nodes];        //double scale[3]
      _xyzo = new double[_dim * n_nodes];        //double scale[3]
      _dxdydz = new double[_dim * n_nodes];
      double * coord;
      coord = new double[n_nodes];

      for ( int  kc = 0; kc < _dim; kc++ ) {
          std::ostringstream Name;
          Name << "NODES/COORD/X" << kc + 1;
          status = H5Dread ( H5Dopen ( /*file_id*/file_id_xyz, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                                  , H5P_DEFAULT
#endif

                                     ),
                             H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );

          for ( int  inode = 0; inode < n_nodes; inode++ ) {
              _xyz[inode + kc * n_nodes] = coord[inode] * ILref; //nondimensionalization
              _xyzo[inode + kc * n_nodes] = coord[inode] * ILref; //nondimensionalization
              _dxdydz[inode + kc * n_nodes] = 0.; //nondimensionalization
              }
          }

      delete []coord;
#ifdef PRINT_INFO
      std::cout << " Reading Multimesh with  " << n_meshes << " meshes, " << n_nodes <<
                " nodes and " << DIMENSION * n_nodes << " coordinates" << std::endl;
#endif
      _off_nd[0] = new int[_n_subdom * _NoLevels + 1];

      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/OFF_ND"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_nd[0] );
      _off_nd[1] = new int[_n_subdom * _NoLevels + 1];
      status = H5Dread ( H5Dopen ( file_id, "/NODES/MAP/OFF_ND1"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_nd[1] );


      // Reading  _NoElements --------------

      _off_el = new int * [_NoFamFEM];             //volume and boundary
      _off_el[0] = new int [_n_subdom * _NoLevels + 1]; //offset for VOLUME
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/OFF_EL"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[0] );
      _off_el[1] = new int [_n_subdom * _NoLevels + 1]; //offset for boundary
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM1/OFF_EL"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _off_el[1] );

      // Number of elements (for each level)
      _NoElements = new int * [_NoFamFEM]; //quadratic and linear

      for ( int  ifem = 0; ifem < _NoFamFEM; ifem++ ) {
          _NoElements[ifem] = new int [_NoLevels];
          std::ostringstream Name;
          Name << "/ELEMS/FEM" << ifem  << "/NExLEV";
          status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                       , H5P_DEFAULT
#endif
                                     ),
                             H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoElements[ifem] );
          }

      // node mapping ------------------------------------------------------------

      _node_map = new int * [_NoLevels + 1];      // _NoLevels+1 inverse node maps
      int  n_nodes_top = _NoNodes[_NoLevels - 1]; //  same length in file (n_nodes_top)
      int * temp = new int [n_nodes_top];         // temp to read

      for ( int  ilev = 0; ilev < _NoLevels; ilev++ ) { // maps from o to 2*_NoLevels-1
          _node_map[ilev] = new int  [_NoNodes[ilev]]; // inverse node map

          // hdf5 reading
          std::ostringstream Name;
          Name << "/NODES/MAP/MAP" << ilev;          // filename
          status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                       , H5P_DEFAULT
#endif
                                     ),  // read from file
                             H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp );

          // inverse map
          for ( int  inode = 0; inode < n_nodes_top; inode++ )
            if ( temp[inode] != -1 ) _node_map[ilev][temp[inode]] = inode;
          }

      // linear coarse map (storage in [2*_NoLevels])
      _node_map[_NoLevels] = new int  [_NoNodes[2 * _NoLevels]];
      // hdf5 reading
      std::ostringstream Name;
      Name << "/NODES/MAP/MAP" << _NoLevels;     // filename
      status = H5Dread ( H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp ); // read from file

      // inverse map
      for ( int  inode = 0; inode < n_nodes_top; inode++ )
        if ( temp[inode] != -1 ) _node_map[_NoLevels][temp[inode]] = inode;

      //connectivity map -----------------------------------------------------
      _el_map = new int * [_NoFamFEM];
      _el_map[0] = new int  [_off_el[0][_NoLevels * _n_subdom]*NDOF_FEM]; // Volume
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/MSH"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),// Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[0] );

      _el_map[1] = new int  [_off_el[1][_NoLevels * _n_subdom]*NDOF_FEMB]; // Boundary
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM1/MSH"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ), // Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_map[1] );


      // neighbour element map -----------------------------------------------------
//     const int n_sides=_GeomEl._n_sides[0];
      _el_neighbor = new int * [_NoFamFEM - 1]; //now only for the volume family
      _el_neighbor[0] = new int [_off_el[0][_NoLevels * _n_subdom]*NDOF_FEM]; // Volume
      status = H5Dread ( H5Dopen ( file_id, "/ELEMS/FEM0/EL_NEIG"
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT
#endif
                                 ),// Read from file
                         H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _el_neighbor[0] );



      // Close the files. ---------------------
      H5Fclose ( file_id );
      H5Fclose ( file_id_xyz );
      delete []temp;

      if ( status < 0 ) {
          std::cout << " File Mesh input in data_in is missing";
          abort();
          }




      } //end restart != 0


  // READ AND STORE BOUNDARY GROUP IDS

  const int   offset     = _NoNodes[_NoLevels - 1];
  std::ostringstream file_bc;
  file_bc  << _mgutils._inout_dir << _mgutils.get_file ( "INMESH" ); //"/mesh.h5";
#ifdef PRINT_INFO
  std::cout << " Reading bc_id from= " <<  file_bc.str() <<  std::endl;
#endif
  hid_t  file_id = H5Fopen ( file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dims[2];

  // face id vector ---------------------------------------------------------
  _NodeBDgroup = new int [offset];

  for ( int i = 0; i < offset; i++ ) {
      _NodeBDgroup[i] = 0;
      }

  // Getting dataset
  std::ostringstream Name;
  Name << "NODES/COORD/BC";
  hid_t dtset = H5Dopen ( file_id, Name.str().c_str()
#if HDF5_VERSIONM != 1808
                          , H5P_DEFAULT
#endif
                        );
  hid_t filespace = H5Dget_space ( dtset ); /* Get filespace handle first. */
  hid_t status  = H5Sget_simple_extent_dims ( filespace, dims, NULL );

  if ( status < 0 ) {
      std::cerr << "GenIc::read dims not found";
      }
  else {   // reading
      assert ( ( int ) dims[0] == offset );
      status = H5Dread ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_NodeBDgroup[0] );
      }

  H5Dclose ( dtset );
  H5Sclose ( filespace );

  return;

  }



/// Write mesh to hdf5 file (namefile)
///              as MGMesh class (MGMesh.h):
/// data parameters : DFL,_type_FEM
/// Vertices        : _NoNodes,_xyz
/// Connectivity    : _NoElements,_conn_map

void MGMesh::write_c ( const int  t_step ) {
  //this is my version commented. His version commented is there.
  //why is that commented?

  std::string    inout_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file ( "BASEMESH" );

  const int  ndigits  = stoi ( _mgutils._sim_config["ndigits"] );

  // Open an existing file. ---------------
  std::ostringstream namefile;
  namefile << inout_dir << basemesh << "." << std::setw ( ndigits ) << std::setfill ( '0' ) << t_step << ".h5";

  // Open file to write multilevel mesh
  std::cout << " Mesh written to=     " << namefile.str() <<  std::endl;
  hid_t file = H5Fcreate ( namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  int  n_meshes = _NoFamFEM * _NoLevels;

  // Writing DFL ---------------------
  int  topdata[3];
  topdata[0] = _dim;
  topdata[1] =  _NoFamFEM;
  topdata[2] = _NoLevels;
  hsize_t dimsf[2];
  dimsf[0] = 3;
  dimsf[1] = 1;
  hid_t dtsp = H5Screate_simple ( 2, dimsf, NULL );
  hid_t dtset = H5Dcreate ( file, "DFL", H5T_NATIVE_INT, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                          );
  H5Dwrite ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, topdata );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );

  // Writing _type_FEM --------------
  dimsf[0] = n_meshes;
  dimsf[1] = 1;
  dtsp = H5Screate_simple ( 2, dimsf, NULL );
  dtset = H5Dcreate ( file, "FEM", H5T_NATIVE_INT, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                      , H5P_DEFAULT, H5P_DEFAULT
#endif
                    );
  H5Dwrite ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _type_FEM );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );



  // Create a group named "/MyGroup" in the file.
  hid_t   group_id = H5Gcreate ( file, "NODES", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                 , H5P_DEFAULT, H5P_DEFAULT
#endif
                               );
  hid_t   group_id2 = H5Gcreate ( file, "NODES/COORD", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                  , H5P_DEFAULT, H5P_DEFAULT
#endif
                                );

  // writing _NoNodes ---------------------
  dimsf[0] = n_meshes;
  dimsf[1] = 1;
  dtsp = H5Screate_simple ( 2, dimsf, NULL );
  dtset = H5Dcreate ( file, "/NNODES", H5T_NATIVE_INT, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                      , H5P_DEFAULT, H5P_DEFAULT
#endif
                    );
  H5Dwrite ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _NoNodes );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );

  //  Writing  _xyz -------------------
  int  n_nodes = _NoNodes[_NoLevels - 1];
  dimsf[0] = n_nodes ;
  dimsf[1] = 1;
  double * coord;
  coord = new double[n_nodes];

  for ( int  kc = 0; kc < _dim; kc++ ) {
      std::ostringstream Name;
      Name << "NODES/COORD/X" << kc + 1;
      dtsp = H5Screate_simple ( 2, dimsf, NULL );

      for ( int  v = 0; v < n_nodes; v++ ) coord[v] = _xyz[v + kc * n_nodes] /*+ _dxdydz[v+kc*n_nodes]*/;

      dtset = H5Dcreate ( file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                          , H5P_DEFAULT, H5P_DEFAULT
#endif
                        );
      H5Dwrite ( dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );
      H5Sclose ( dtsp );
      H5Dclose ( dtset );
      }


#if DIMENSION==1
  std::ostringstream Name;
  Name << "NODES/COORD/X" << 2;
  dtsp = H5Screate_simple ( 2, dimsf, NULL );

  for ( int  v = 0; v < n_nodes; v++ ) coord[v] = 0.;

  dtset = H5Dcreate ( file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite ( dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );


  std::ostringstream Nameb;
  Nameb << "NODES/COORD/X" << 3;
  dtsp = H5Screate_simple ( 2, dimsf, NULL );

  for ( int  v = 0; v < n_nodes; v++ ) coord[v] = 0.;

  dtset = H5Dcreate ( file, Nameb.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite ( dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );
#endif


#if DIMENSION==2
  std::ostringstream Name;
  Name << "NODES/COORD/X" << 3;
  dtsp = H5Screate_simple ( 2, dimsf, NULL );

  for ( int  v = 0; v < n_nodes; v++ ) coord[v] = 0.;

  dtset = H5Dcreate ( file, Name.str().c_str(), H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                      , H5P_DEFAULT, H5P_DEFAULT
#endif
                    );
  H5Dwrite ( dtset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord );
  H5Sclose ( dtsp );
  H5Dclose ( dtset );
#endif


// Close the group
  H5Gclose ( group_id2 );
  H5Gclose ( group_id );
  delete []coord;




// // //
// // //
// // //   // writing _NoElements --------------
// // //   dimsf[0] = n_meshes;  dimsf[1] = 1;
// // //   dtsp =H5Screate_simple (2, dimsf, NULL);
// // //   dtset=H5Dcreate(file,"/CONN/NELEMENTS",H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
// // //   H5Dwrite (dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, _NoElements);
// // //   H5Sclose (dtsp);  H5Dclose (dtset);
// // //
// // //   // Writing  Connectivity ---------------
// // //   for (int  ifem=0; ifem<_NoFamFEM;ifem++){
// // //     for (int  ilev=0; ilev<_NoLevels;ilev++) {
// // //       // int   *dconn;dconn=new int [_NoElements[ilev+ifem*_NoLevels]*8];
// // //       // int  n_len=sub_conn(dconn,ifem,ilev,8);
// // //       //std::cout <<  n_len << " n_len " << _NoElements[ilev+ifem*_NoLevels]*8<< std::endl;
// // //       std::ostringstream Name; Name << "/CONN/FEM"<< ifem+1 << "/MSH"<< ifem*_NoLevels+ilev;
// // //       dimsf[0] = _NoElements[ifem][ilev]*_type_FEM[ifem];
// // //       dimsf[1] = 1;
// // //       //dimsf[0] = n_len;  dimsf[1] = 1;
// // //       dtsp = H5Screate_simple (2, dimsf, NULL);
// // //       dtset = H5Dcreate (file,Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
// // //       H5Dwrite(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,_el_map[ifem]);
// // //       // H5Dwrite(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,dconn);
// // //       H5Sclose (dtsp);  H5Dclose (dtset);
// // //     }
// // //   }



  // Close the file -------------
  H5Fclose ( file );



  return;
  }


// ========================================================
/// This function prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM
void MGMesh::print_conn_lin_hf5 (
  const int  Level      // Level <-
)  { // ================================================

  std::string    basemesh = _mgutils.get_file ( "BASEMESH" );
  std::string    input_dir = _mgutils._inout_dir;
  std::string    connlin = _mgutils.get_file ( "CONNLIN" );

  int conn[8][8];
  int * gl_conn;
  // storage in hf5 (Xdmf)
  std::ostringstream namefile;
  namefile << input_dir << basemesh << connlin << ".h5";
#ifdef PRINT_INFO
  std::cout << "Print mesh connectivity to =" << namefile.str() << std::endl;
#endif
  hid_t file = H5Fcreate ( namefile.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );


  for ( int  ik = 0; ik < _NoFamFEM; ik++ ) {
      int  icount = 0;
      int  indx_mesh = Level + _NoLevels * ik;
      int  mode = _type_FEM[indx_mesh];
      int  n_elements = _NoElements[ik][Level];//  int  n_elements = _NoElements[indx_mesh];
      int  nsubel, nnodes;

      switch ( mode ) {
// #if DIMENSION==2
        // -----------------------------------
        case 9: // Quad 9  0-4-1-5-2-6-3-7-8
          nsubel = 4;
          nnodes = 4; // Quad9 -> 4 Quad4
          gl_conn = new int [n_elements * 4 * 4];
          conn[0][0] = 0;
          conn[0][1] = 4;
          conn[0][2] = 8;
          conn[0][3] = 7;// quad4  0-4-8-7
          conn[1][0] = 4;
          conn[1][1] = 1 ;
          conn[1][2] = 5;
          conn[1][3] = 8;// quad4  4-1-5-8
          conn[2][0] = 8;
          conn[2][1] = 5;
          conn[2][2] = 2;
          conn[2][3] = 6;// quad4  8-5-2-6
          conn[3][0] = 7;
          conn[3][1] = 8;
          conn[3][2] = 6;
          conn[3][3] = 3;// quad4  7-8-6-3
          break;

//=====================
        case 6: //  Tri6  0-3-1-4-2-5
          nsubel = 4;
          nnodes = 3; //  Tri6 -> 4 Tri3
          gl_conn = new int [n_elements * 4 * 3];
          conn[0][0] = 0;
          conn[0][1] = 3;
          conn[0][2] = 5; // quad4  0-4-8-7
          conn[1][0] = 3;
          conn[1][1] = 4;
          conn[1][2] = 5; // quad4  4-1-5-8
          conn[2][0] = 3;
          conn[2][1] = 1;
          conn[2][2] = 4; // quad4  8-5-2-6
          conn[3][0] = 4;
          conn[3][1] = 2;
          conn[3][2] = 5; // quad4  7-8-6-3
          break;

//===================
        case 3: // Edge 3  0-2-1
          nsubel = 2;
          nnodes = 2; // Edge 3 -> 2 Edge 2
          gl_conn = new int [n_elements * 2 * 2];
          conn[0][0] = 0;
          conn[0][1] = 2;   // element 0-2
          conn[1][0] = 2;
          conn[1][1] = 1;   // element 1-2
          break;

// #else
        // ----------------------
        case  27: //  Hex 27 (8 Hex8)
          nsubel = 8;
          nnodes = 8; //  Hex 27 -> 8 Hex8
          gl_conn = new int [n_elements * 8 * 8];
          conn[0][0] = 0;
          conn[0][1] = 8;
          conn[0][2] = 20;
          conn[0][3] = 11;
          conn[0][4] = 12;
          conn[0][5] = 21;
          conn[0][6] = 26;
          conn[0][7] = 24;
          conn[1][0] = 8;
          conn[1][1] = 1;
          conn[1][2] = 9;
          conn[1][3] = 20;
          conn[1][4] = 21;
          conn[1][5] = 13;
          conn[1][6] = 22;
          conn[1][7] = 26;
          conn[2][0] = 11;
          conn[2][1] = 20;
          conn[2][2] = 10;
          conn[2][3] = 3;
          conn[2][4] = 24;
          conn[2][5] = 26;
          conn[2][6] = 23;
          conn[2][7] = 15;
          conn[3][0] = 20;
          conn[3][1] = 9;
          conn[3][2] = 2;
          conn[3][3] = 10;
          conn[3][4] = 26;
          conn[3][5] = 22;
          conn[3][6] = 14;
          conn[3][7] = 23;
          conn[4][0] = 12;
          conn[4][1] = 21;
          conn[4][2] = 26;
          conn[4][3] = 24;
          conn[4][4] = 4;
          conn[4][5] = 16;
          conn[4][6] = 25;
          conn[4][7] = 19;
          conn[5][0] = 21;
          conn[5][1] = 13;
          conn[5][2] = 22;
          conn[5][3] = 26;
          conn[5][4] = 16;
          conn[5][5] = 5;
          conn[5][6] = 17;
          conn[5][7] = 25;
          conn[6][0] = 24;
          conn[6][1] = 26;
          conn[6][2] = 23;
          conn[6][3] = 15;
          conn[6][4] = 19;
          conn[6][5] = 25;
          conn[6][6] = 18;
          conn[6][7] = 7;
          conn[7][0] = 26;
          conn[7][1] = 22;
          conn[7][2] = 14;
          conn[7][3] = 23;
          conn[7][4] = 25;
          conn[7][5] = 17;
          conn[7][6] = 6;
          conn[7][7] = 18;
          break;

        // ---------------------------------------
        case  10: // Tet10 -> 8 Tet4
          nsubel = 8;
          nnodes = 4; // Tet10 -> 8 Tet4
          gl_conn = new int [n_elements * 8 * 4];
          conn[0][0] = 0;
          conn[0][1] = 4;
          conn[0][2] = 6;
          conn[0][3] = 7;
          conn[1][0] = 4;
          conn[1][1] = 1;
          conn[1][2] = 5;
          conn[1][3] = 8;
          conn[2][0] = 5;
          conn[2][1] = 2;
          conn[2][2] = 6;
          conn[2][3] = 9;
          conn[3][0] = 7;
          conn[3][1] = 8;
          conn[3][2] = 9;
          conn[3][3] = 3;
          conn[4][0] = 4;
          conn[4][1] = 8;
          conn[4][2] = 6;
          conn[4][3] = 7;
          conn[5][0] = 4;
          conn[5][1] = 5;
          conn[5][2] = 6;
          conn[5][3] = 8;
          conn[6][0] = 5;
          conn[6][1] = 9;
          conn[6][2] = 6;
          conn[6][3] = 8;
          conn[7][0] = 7;
          conn[7][1] = 6;
          conn[7][2] = 9;
          conn[7][3] = 8;
          break;

        // ---------------------------------------
        // -----------------------------------------
        default:   // interior 3D
          nsubel = 1;
          nnodes = mode;
          gl_conn = new int [n_elements * nsubel * nnodes];

          for ( int  n = 0; n < mode; n++ ) conn[0][n] = n;

          break;
          }

      // mapping
      for ( int  iproc = 0; iproc < _n_subdom; iproc++ ) {
          for ( int el = _off_el[ik][iproc * _NoLevels + Level];
                el < _off_el[ik][iproc * _NoLevels + Level + 1]; el++ ) {
              for ( int  se = 0; se < nsubel; se++ ) {
                  for ( int  i = 0; i < nnodes; i++ ) {
                      gl_conn[icount] = _el_map[ik][el * mode + conn[se][i]];
                      icount++;
                      }
                  }
              }
          }

      // Print mesh in hdf files
      std::ostringstream Name;
      Name << "MSH" << ik << "CONN";
      hsize_t dimsf[2];
      dimsf[0] = n_elements * nsubel * nnodes;
      dimsf[1] = 1;
      hid_t dtsp = H5Screate_simple ( 2, dimsf, NULL );
      hid_t dtset = H5Dcreate ( file, Name.str().c_str(), H5T_NATIVE_INT, dtsp, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                , H5P_DEFAULT, H5P_DEFAULT
#endif
                              );
      H5Dwrite ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gl_conn );
      H5Sclose ( dtsp );
      H5Dclose ( dtset );

      delete[]gl_conn;
      }

  H5Fclose ( file );
  return;
  }

// ====================================================
/// This function prints the subdomain subdivision (proc)
void MGMesh::print_subdom_hf5 (
  std::string filename // filename
) const { // ==========================================

  // setup ucoord
  int       n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = ( DIMENSION > 1 ) ? 4 * ( DIMENSION - 1 ) : 2;
  double * ucoord;
  ucoord = new double[n_subdivision * n_elements];
  // storage PID -> ucoord -----------------------------
  int cel = 0;

  for ( int  iproc = 0; iproc < _n_subdom; iproc++ ) {
      for ( int iel = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];
            iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1]; iel++ ) {
          for ( int  is = 0; is < n_subdivision; is++ )
            ucoord[cel * n_subdivision + is] = iproc;

          cel++;
          }
      }

  // print to hdf5 format ------------------------------
  hid_t file_id = H5Fopen ( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dimsf[2];
  dimsf[0] = n_subdivision * n_elements;
  dimsf[1] = 1;
  hid_t dataspace = H5Screate_simple ( 2, dimsf, NULL );
  hid_t dataset = H5Dcreate ( file_id, "PID", H5T_NATIVE_DOUBLE,
                              dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                              , H5P_DEFAULT, H5P_DEFAULT
#endif
                            );
  hid_t  status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord );
  assert ( status == 0 );
  H5Sclose ( dataspace );
  H5Dclose ( dataset );
  H5Fclose ( file_id );

  // clean ----------------------------------------------
  delete []ucoord;
//   }

  return;
  }

// ====================================================
/// This function prints the distance from the wall
void MGMesh::print_dist_hf5 (
  std::string filename, // filename
  double dist[],
  std::string dir_name
) const { // ==========================================

  // setup ucoord ---------------------------------------
  int    n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = ( DIMENSION > 1 ) ? 4 * ( DIMENSION - 1 ) : 2;
  double * ucoord;
  ucoord = new double[n_elements * n_subdivision];
  // storage DIST -> ucoord -----------------------------
  int cel = 0;

  for ( int  iproc = 0; iproc < _n_subdom; iproc++ ) {
      const int  nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];

      for ( int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++ ) {
          for ( int is = 0; is < n_subdivision; is++ )
            ucoord[cel * n_subdivision + is] = dist[iel + nel_b];

          cel++;
          }
      }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen ( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dimsf[2];
  dimsf[0] = n_elements * n_subdivision;
  dimsf[1] = 1;

  hid_t dataspace = H5Screate_simple ( 2, dimsf, NULL );
  hid_t dataset = H5Dcreate ( file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE,
                              dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                              , H5P_DEFAULT, H5P_DEFAULT
#endif
                            );
  hid_t  status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord );
  assert ( status == 0 );
  H5Sclose ( dataspace );
  H5Dclose ( dataset );
  H5Fclose ( file_id );

  // clean --------------------------------------------------------
  delete []ucoord;

  return;
  }


void MGMesh::print_VolFrac_hf5 (
  std::string filename, // filename
  std::string dir_name
) const { // ==========================================

  // setup ucoord ---------------------------------------
  int    n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = ( DIMENSION > 1 ) ? 4 * ( DIMENSION - 1 ) : 2;
  double * ucoord;
  ucoord = new double[n_elements * n_subdivision];
  // storage DIST -> ucoord -----------------------------
  int cel = 0;

  for ( int  iproc = 0; iproc < _n_subdom; iproc++ ) {
      const int  nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];

      for ( int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++ ) {
          for ( int is = 0; is < n_subdivision; is++ )
            ucoord[cel * n_subdivision + is] = _VolFrac[iel + nel_b];

          cel++;
          }
      }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen ( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dimsf[2];
  dimsf[0] = n_elements * n_subdivision;
  dimsf[1] = 1;

  hid_t dataspace = H5Screate_simple ( 2, dimsf, NULL );
  hid_t dataset = H5Dcreate ( file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE,
                              dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                              , H5P_DEFAULT, H5P_DEFAULT
#endif
                            );
  hid_t  status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord );
  assert ( status == 0 );
  H5Sclose ( dataspace );
  H5Dclose ( dataset );
  H5Fclose ( file_id );

  // clean --------------------------------------------------------
  delete []ucoord;

  return;
  }

// ====================================================
/// This function prints the control domain
void MGMesh::print_ctrl_dom_hf5 (
  std::string filename, // filename
  double ctrl_dom[],
  std::string dir_name
) const { // ==========================================

  // setup ucoord ---------------------------------------
  int    n_elements = _NoElements[0][_NoLevels - 1];
  int n_subdivision = ( DIMENSION > 1 ) ? 4 * ( DIMENSION - 1 ) : 2;
  double * ucoord;
  ucoord = new double[n_elements * n_subdivision];
  // storage ctrl_dom -> ucoord -----------------------------
  int cel = 0;

  for ( int  iproc = 0; iproc < _n_subdom; iproc++ ) {
      const int  nel_b = _off_el[0][_NoLevels - 1 + iproc * _NoLevels];

      for ( int iel = 0; iel < _off_el[0][_NoLevels - 1 + iproc * _NoLevels + 1] - nel_b; iel++ ) {
          for ( int is = 0; is < n_subdivision; is++ )
            ucoord[cel * n_subdivision + is] = ctrl_dom[iel + nel_b];

          cel++;
          }
      }

  // print to hdf5 format -----------------------------------------
  hid_t file_id = H5Fopen ( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dimsf[2];
  dimsf[0] = n_elements * n_subdivision;
  dimsf[1] = 1;

  hid_t dataspace = H5Screate_simple ( 2, dimsf, NULL );
  hid_t dataset = H5Dcreate ( file_id, dir_name.c_str(), H5T_NATIVE_DOUBLE,
                              dataspace, H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                              , H5P_DEFAULT, H5P_DEFAULT
#endif
                            );
  hid_t  status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ucoord );
  assert ( status == 0 );
  H5Sclose ( dataspace );
  H5Dclose ( dataset );
  H5Fclose ( file_id );

  // clean --------------------------------------------------------
  delete []ucoord;

  return;
  }

