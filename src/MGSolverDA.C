// #include "Equations_tab.h"
// #include "Equations_conf.h"
// ===============================
// ===============================
//std lib
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>

// configuration files ----------
#include "Printinfo_conf.h"

// class files ----------------
#include "MGSolverDA.h"
// #include "MGSclass_conf.h"

// local includes --------------
#include "MeshExtended.h"
#include "MGGeomEl.h"
#include "MGFE_conf.h"
#include "dense_set.h"
#include "MGSystem.h"
#include "EquationSystemsExtendedM.h"

#include "MGUtils.h"
#include "MGEquationsSystem.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGGraph.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif

// algebric includes -----------
#include "numeric_vectorM.h"
#include "sparse_matrixM.h"
#include "sparse_MmatrixM.h"
#include "dense_vectorM.h"
#include "dense_matrixM.h"



// ==========================================================
//               MGSolDA functions
// ==========================================================

/// This function is the MGSolDA constructor I
MGSolDA::MGSolDA (
  MGEquationsSystem & mg_equations_map_in, //
  const int nvars_in[],                // # of quad variables
//   const int nvarsl_in,                // # of linear variables
  std::string eqname_in,               // equation name
  std::string /*varname_in*/ )             // basic variable name
  : MGSolBase ( mg_equations_map_in, nvars_in, eqname_in ) { // ===============================================
  //piecewise constant linear, quadratic and (...cubic)
  _fe[0] =   _mgfemap.get_FE ( 0 ); ///> Lagrange piecewise constant
  _fe[1] =   _mgfemap.get_FE ( 1 ); ///> Lagrange piecewise linear
  _fe[2] =   _mgfemap.get_FE ( 2 ); ///> Lagrange piecewise quadratic

  _el_dof[0] = ( ( nvars_in[0] > 0 ) ? NDOF_K : 0 ); //Lagrange piecewise constant variables
  _el_dof[1] = ( ( nvars_in[1] > 0 ) ? NDOF_P : 0 ); //Lagrange piecewise linear variables
  _el_dof[2] = ( ( nvars_in[2] > 0 ) ? NDOF_FEM : 0 ); //Lagrange piecewise linear variables

  // System names and units

  _var_names = new std::string[_n_vars];           // names
  _refvalue = new double[_n_vars];
  //output string stream
  int iname = 0;

  for ( ; iname < nvars_in[2]; iname++ ) { //quadratic
      std::ostringstream ostr;
      ostr << "qua" << iname + 1; //use the string stream just like cout,
      _var_names[iname] = ostr.str();
      _refvalue[iname] = 1;
      }

  for ( ; iname < ( nvars_in[1] + nvars_in[2] ); iname++ ) { //linear
      std::ostringstream ostr;
      ostr << "lin" << iname - nvars_in[2] + 1; //use the string stream just like cout,
      _var_names[iname] = ostr.str();
      _refvalue[iname] = 1;
      }

  for ( ; iname < ( nvars_in[0] + nvars_in[1] + nvars_in[2] ); iname++ ) { //piecewise
      std::ostringstream ostr;
      ostr << "pie" << iname - nvars_in[2] - nvars_in[1] + 1; //use the string stream just like cout,
      _var_names[iname] = ostr.str();
      _refvalue[iname] = 1;
      }

  _NumRestartSol = 1;
  return;
  }


// =========================================
/// This function sets up data structures for each problem class
// =========================================
void MGSolDA::setUpExtFieldData() {
  /// A) set up _mg_eqs
  // external system and index vectors
  for ( int deg = 0; deg < 3; deg++ ) {
      for ( int kl = 0; kl < _data_eq[deg].max_neqs; kl++ ) {
          _data_eq[deg].n_eqs = 0;
          _data_eq[deg].mg_eqs[kl] = NULL;
          _data_eq[deg].indx_ub[kl] = -1;
          _data_eq[deg].tab_eqs[kl] = -1;

          for ( int kk = 0; kk < NDOF_FEM; ++kk )  {
              _data_eq[deg].ub[kk + kl * NDOF_FEM] = 0.; // data
              }
          }
      }

// start index K from 0, L from 0, Q from DIMENSION (coordinates+ q variable)
  _data_eq[0].indx_ub[0] = 0;       //_data_eq[0].n_eqs=0; // piecewice constant  (0)
  _data_eq[1].indx_ub[0] = 0;       //_data_eq[1].n_eqs=0; // piecewice linear    (1)
  _data_eq[2].indx_ub[0] = 0;       //_data_eq[2].n_eqs=0; // piecewice quadratic (2)
  return;
  }

// =========================================
/// This function is the destructor
// =========================================
MGSolDA::~MGSolDA() {
//     delete []_data_eq;
  }



// =========================================
/// This function is the main initialization function:
// =========================================
void MGSolDA::init ( const int Level ) {

  std::string    f_matrix = _mgutils.get_file ( "F_MATRIX" );
  std::string    f_rest = _mgutils.get_file ( "F_REST" );
  std::string    f_prol = _mgutils.get_file ( "F_PROL" );

  const ParallelM::Communicator & comm1 = _mgmesh._comm.comm();
  // number of partial dofs
  int off_proc = _NoLevels * _iproc;
//   const int  n_glob_q=_mgmesh._NoNodes[Level];
//     int n_glob_l=_mgmesh._NoNodes[_mgmesh._NoFamFEM*_NoLevels];
//     if (Level>0) n_glob_l=_mgmesh._NoNodes[Level-1];
  int n_local_q = ( _el_dof[2] > 0 ) ? _mgmesh._off_nd[0][Level + 1 + off_proc] - _mgmesh._off_nd[0][off_proc] : 0; //0=QUADRATIC
  int n_local_l = ( _el_dof[1] > 0 ) ? _mgmesh._off_nd[1][Level + 1 + off_proc] - _mgmesh._off_nd[1][off_proc] : 0; //1=LINEAR
  int n_local_k = ( _el_dof[0] > 0 ) ? ( _mgmesh._off_el[0][Level + 1 + off_proc] - _mgmesh._off_el[0][Level + off_proc] ) * _el_dof[0] : 0; //0=Volume

  // global dofs
  int n_glob = _Dim[Level];
  int n_local = _nvars[1] * n_local_l + _nvars[2] * n_local_q + _nvars[0] * n_local_k;

  // quad-linear matrices


  //   matrix -----------------------------
  std::ostringstream filename ( "" );
//   filename.str ( "" );
  filename << _mgutils._inout_dir << f_matrix << Level << ".h5";
  A[Level] = SparseMatrixM::build ( _mgmesh._comm.comm() ).release();
  A[Level]->init ( n_glob, n_glob, n_local, n_local );

  ReadMatrix ( Level, filename.str(), *A[Level], _nvars );

//   int  nel_glob=_mgmesh._NoElements[0][Level];
//   int nel_local=_mgmesh._off_el[0][Level+off_proc+1]-_mgmesh._off_el[0][Level+off_proc];

  // vectors --------------------------------
  b[Level] = NumericVectorM::build ( comm1 ).release();
  b[Level]->init ( n_glob, n_local, false, AUTOMATICM );

  res[Level] = NumericVectorM::build ( comm1 ).release();
  res[Level]->init ( n_glob, n_local, false, AUTOMATICM );
  x[Level] = NumericVectorM::build ( comm1 ).release();
  x[Level]->init ( n_glob, n_local, false, AUTOMATICM );
  x_old[Level] = NumericVectorM::build ( comm1 ).release();
  x_old[Level]->init ( n_glob, false, SERIALM );
  x_oold[Level] = NumericVectorM::build ( comm1 ).release();
  x_oold[Level]->init ( n_glob, false, SERIALM );
  x_nonl[Level] = NumericVectorM::build ( comm1 ).release();
  x_nonl[Level]->init ( n_glob, false, SERIALM );
  disp[Level] = NumericVectorM::build ( comm1 ).release();
  disp[Level]->init ( n_glob, false, SERIALM );
  disp_old[Level] = NumericVectorM::build ( comm1 ).release();
  disp_old[Level]->init ( n_glob, false, SERIALM );
  disp_oold[Level] = NumericVectorM::build ( comm1 ).release();
  disp_oold[Level]->init ( n_glob, false, SERIALM );
  disp_ooold[Level] = NumericVectorM::build ( comm1 ).release();
  disp_ooold[Level]->init ( n_glob, false, SERIALM );
  x_ooold[Level] = NumericVectorM::build ( comm1 ).release();
  x_ooold[Level]->init ( n_glob, false, SERIALM );
  x_oooold[Level] = NumericVectorM::build ( comm1 ).release();
  x_oooold[Level]->init ( n_glob, false, SERIALM );


  if ( Level < _NoLevels - 1 ) { // Restrictor
      Rst[Level] = SparseMMatrixM::build ( comm1 ).release();
      Rst[Level]->init ( n_glob, _Dim[Level + 1], n_glob, _Dim[Level + 1] );
//     int n_nodes_qp1= _mgmesh._NoNodes[Level+1];
      filename.str ( "" );
      filename << _mgutils._inout_dir  << f_rest << Level + 1  << "_" << Level << ".h5";
      ReadRest ( Level, filename.str(), *Rst[Level], _nvars, _node_dof[Level + 1], _node_dof[Level], _node_dof[_NoLevels - 1] );
      Rst[Level]->close();
      }

  if ( Level > 0 ) {  // Prolongation
      Prl[Level] = SparseMMatrixM::build ( comm1 ).release();
      Prl[Level]->init ( n_glob, _Dim[Level - 1], _Dim[Level], _Dim[Level - 1] );
      filename.str ( "" );
      filename << _mgutils._inout_dir  << f_prol << Level - 1  << "_" << Level << ".h5";
      ReadProl ( Level, filename.str(), *Prl[Level], _nvars, _node_dof[Level], _node_dof[Level - 1] );
      Prl[Level]->close();
//  Prl[Level]->print_personal(filename);//print on screen Prol
      }

  return;
  }






// ============================================================================
/// This function generates the initial conditions for a Base system:
void MGSolDA::GenIc() {

  std::string input_dir = _mgutils._inout_dir;
  std::string ibc = _mgutils.get_file ( "IBC" );

  std::ostringstream ibc_file;

  ibc_file << input_dir << ibc << ".h5";

  std::ifstream in ( ( ibc_file.str() ).c_str() );

  if ( !in ) {  // ------- function initial condition -> user fuinction -> ic_read

      // Get a constant reference to the mesh objects.
      const double * xyz_glob = _mgmesh._xyz;               // global coordinates
      const int * map_nodes  = _mgmesh._el_map[0];          // node map (elem,loc_nodes)->gl_nodes
      const int   offset     = _mgmesh._NoNodes[_NoLevels - 1]; // offset variables (on level _NoLevels-1)
      const int  * off_el    = _mgmesh._off_el[0];          // offset element
      int face_id_node = 0;
      int mat_id_elem = 0;
      NumericVectorM & sol_top = *x[_NoLevels - 1]; // solution (top level)
      NumericVectorM & old_sol_top = *x_old[_NoLevels - 1];
      NumericVectorM & nl_sol_top = *x_nonl[_NoLevels - 1];
      NumericVectorM & oold_sol_top = *x_oold[_NoLevels - 1];

      const int * node_dof_top = &_node_dof[_NoLevels - 1][0];
      int ntot_elements = 0;

      for ( int ilev = 0; ilev < _NoLevels; ilev++ ) {
          ntot_elements += _mgmesh._NoElements[0][ilev];
          }

      // temp vect
      double * u_value = new double[_n_vars];
      double xp[DIMENSION];

      // loop reading
      int off_proc = _iproc * _NoLevels;
      int ndof_lev = 0;

      // Reading  bc_id *********************************************************
      // Open an existing file -------------------------------------------------
      std::ostringstream file_bc;
      file_bc  << _mgutils._inout_dir << _mgutils.get_file ( "INMESH" ); //"/mesh.h5";
#ifdef PRINT_INFO
      std::cout << " Reading bc_id from= " <<  file_bc.str() <<  std::endl;
#endif
      hid_t  file_id = H5Fopen ( file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
      hsize_t dims[2];

      // face id vector ---------------------------------------------------------
      int * face_id_vect;
      face_id_vect = new int [offset];

      for ( int i = 0; i < offset; i++ ) {
          face_id_vect[i] = 0;
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
          status = H5Dread ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &face_id_vect[0] );
          }

      H5Dclose ( dtset );
      H5Sclose ( filespace );
      // Reading  mat_id ********************************************************
      // mat id vector
      int * mat_id_vect = new int [ntot_elements];

      for ( int i = 0; i < ntot_elements; i++ ) {
          mat_id_vect[i] = 1;
          }

      // Getting dataset
      std::ostringstream Name1;
      Name1 << "ELEMS/SUB/MAT";
      dtset = H5Dopen ( file_id, Name1.str().c_str()
#if HDF5_VERSIONM != 1808
                        , H5P_DEFAULT
#endif
                      );
      filespace = H5Dget_space ( dtset ); /* Get filespace handle first. */
      status  = H5Sget_simple_extent_dims ( filespace, dims, NULL );

      if ( status < 0 ) {
          std::cerr << "GenIc::read mat dims not found";
          }
      else {   // reading
          status = H5Dread ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, &mat_id_vect[0] );
          }  // end else

      H5Sclose ( filespace );
      H5Dclose ( dtset );
      H5Fclose ( file_id );
//  end Reading  mat_id ******************************************************

      for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
          int delta = off_el[pr * _NoLevels + _NoLevels - 1 + 1] - off_el[pr * _NoLevels + _NoLevels - 1];
          ndof_lev += delta;
          }

      for ( int iel = 0; iel < off_el[off_proc + _NoLevels] - off_el[off_proc + _NoLevels - 1]; iel++ ) {
          int elem_gidx = ( iel + off_el[_iproc * _NoLevels + _NoLevels - 1] ) * NDOF_FEM;
          int elem_indx = ( iel + ndof_lev ) * _el_dof[0];
          mat_id_elem = mat_id_vect[iel + ndof_lev];

          // the local nodes
          for ( int i = 0; i < NDOF_FEM; i++ ) {
              // coordinates
              int k = map_nodes[elem_gidx + i];

              for ( int idim = 0; idim < DIMENSION; idim++ ) {
                  xp[idim] = xyz_glob[k + idim * offset];
                  }

              face_id_node = face_id_vect[k];

              // ===================================================
              // user definition reading function ----------------==
              ic_read ( face_id_node, mat_id_elem, xp, iel + off_el[off_proc + _NoLevels - 1], u_value ); //-----------------------==
              // -------------------------------------------------==
              // ===================================================

              // Set discontinuous fields
              if ( i == NDOF_FEM - 1 )     for ( int ivar = 0; ivar < _nvars[0]; ivar++ ) {
                    sol_top.set ( node_dof_top[elem_indx + ( ivar + _nvars[2] + _nvars[1] ) *offset], u_value[_nvars[2] + _nvars[1] + ivar] );

                    for ( int kdof0 = 1; kdof0 < _el_dof[0]; kdof0++ ) {
                        sol_top.set ( node_dof_top[ kdof0 + elem_indx + ( ivar + _nvars[2] + _nvars[1] ) *offset], 0. );
                        }
                    }

              // Set the quadratic and linear fields
              if ( i < _el_dof[1] ) for ( int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++ ) {
                    sol_top.set ( node_dof_top[k + ivar * offset], u_value[ivar] );
                    int a = 1;
                    }
              else  for ( int ivar = 0; ivar < _nvars[2]; ivar++ ) {
                    int irrr = node_dof_top[k + ivar * offset];
                    sol_top.set ( irrr, u_value[ivar] );
                    }
              }
          } // end of element loop

//     // delocalize
      sol_top.localize ( old_sol_top );
      sol_top.localize ( nl_sol_top );
      sol_top.localize ( oold_sol_top );

      //clean
      delete []u_value;
#ifdef PRINT_INFO
      std::cout << "\n GenIc(DA): Initial solution defined by    ic_read(xp,u_value)" << "\n \n";
#endif
      delete [] face_id_vect;
      delete [] mat_id_vect;
      }
  else {  // -------------------- file reading --> data_in/case.h5
      const int restart_lev = stoi ( _mgutils._sim_config["restart_lev"] ); // restart label
      read_u ( ibc_file.str(), restart_lev );
#ifdef PRINT_INFO
      std::cout << "\n GenIc(DA): Initial solution defined by " << ibc_file.str() << "\n \n";
#endif
      } // +++++++++++++++++++++++++++++++++++++++++++++


  in.close();
  return;
  }

// ============================================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::GenBc (
) { // ========================================================================
  /// A) Set up: mesh,dof, bc
  //mesh ----------------------------------------------------------------------
  const int offset    = _mgmesh._NoNodes[_NoLevels - 1];
  int ntot_elements = 0;

  for ( int ilev = 0; ilev < _NoLevels; ilev++ ) {
      ntot_elements += _mgmesh._NoElements[0][ilev];
      }

  // Dof ----------------------------------------------------------------------
  const int n_kb_dofs = ( ( _nvars[0] > 0 ) ? DIMENSION + 1 : 0 ); // surface dofs
  const int n_pb_dofs = ( ( _nvars[1] > 0 ) ? _fe[1]->_NoShape[DIMENSION - 2] : 0 ); //get_n_shapes(DIMENSION-2);
  const int n_ub_dofs = ( ( _nvars[2] > 0 ) ? _fe[2]->_NoShape[DIMENSION - 2] : 0 ); //get_n_shapes(DIMENSION-2);
//   const int  n_dofs =  n_pb_dofs*_nvars[1] + n_ub_dofs*_nvars[2];
  const int n_k_dofs = ( ( _nvars[0] > 0 ) ?  DIMENSION + 1 : 0 ); // volume dofs
  const int n_l_dofs = ( ( _nvars[1] > 0 ) ? _fe[1]->_NoShape[DIMENSION - 1] : 0 ); //get_n_shapes(DIMENSION-1);
  const int n_u_dofs = ( ( _nvars[2] > 0 ) ? _fe[2]->_NoShape[DIMENSION - 1] : 0 ); //get_n_shapes(DIMENSION-1);

  // set 1 all the points for  bc (boundary condition) ------------------------
  for ( int i1 = 0; i1 < _Dim[_NoLevels - 1]; i1++ ) {
      _bc[0][i1] = 1;
      _bc[1][i1] = 1;
      }

  // **************************************************************************
  // B) Reading  face_id vector (boundary zones) if the dataset exists
  // Open an existing file ----------------------------------------------------
  std::ostringstream file_bc;
  file_bc  << _mgutils._inout_dir << _mgutils.get_file ( "INMESH" ); //"/mesh.h5";
#ifdef PRINT_INFO
  std::cout << " Reading bc_id from= " <<  file_bc.str() <<  std::endl;
#endif
  hid_t  file_id = H5Fopen ( file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t dims[2];
  // face id vector initial setting (set 0) -----------------------------------
  int * face_id_vect = new int [offset];

  for ( int i = 0; i < offset; i++ ) {
      face_id_vect[i] = 0;
      }

  // Getting dataset ----------------------------------------------------------
  std::ostringstream Name;
  Name << "NODES/COORD/BC";
  hid_t dtset = H5Dopen ( file_id, Name.str().c_str(), H5P_DEFAULT );
  hid_t filespace = H5Dget_space ( dtset ); /* Get filespace handle first. */
  hid_t status  = H5Sget_simple_extent_dims ( filespace, dims, NULL );

  if ( status < 0 ) {
      std::cerr << "GenIc::read dims not found";
      }
  else {   // reading (otherwise it stays 0)
      assert ( ( int ) dims[0] == offset );
      status = H5Dread ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &face_id_vect[0] );
      }

  H5Dclose ( dtset );
  H5Sclose ( filespace );
  // **************************************************************************
  // C) Reading  mat_id (volume zones) if the dataset exists
  // mat id vector initialization ---------------------------------------------
  int * mat_id_vect = new int [ntot_elements];

//   int *mat_id_vect_lev; int icount=0;
  for ( int i = 0; i < ntot_elements; i++ ) {
      mat_id_vect[i] = 1;
      }

  // level loop (in the file are written for each level)
  // Getting dataset --------------------------------------------------------
  std::ostringstream Name1;
  Name1 << "ELEMS/SUB/MAT";
  dtset = H5Dopen ( file_id, Name1.str().c_str(), H5P_DEFAULT );


  filespace = H5Dget_space ( dtset ); /* Get filespace handle first. */
  status  = H5Sget_simple_extent_dims ( filespace, dims, NULL );

  if ( status < 0 ) {
      std::cerr << "GenIc::read mat dims not found";
      }
  else {   // reading if the dataset exists (otherwise it stays 1) ------------
      status = H5Dread ( dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat_id_vect[0] );
      }  // end reading ---------------------------------------------------------

  // clean --------------------------------------------------------------------
  H5Sclose ( filespace );
  H5Dclose ( dtset );
  H5Fclose ( file_id );
  // *************************************************************************
  /// C  reading bc from function --> bc_read

  GenBc_loop ( 0, NDOF_FEM, n_u_dofs, n_l_dofs, n_k_dofs,
               face_id_vect, mat_id_vect );


//   GenBc_loop(1,NDOF_FEMB,n_ub_dofs,n_pb_dofs,n_kb_dofs,
//              face_id_vect,mat_id_vect);
  // clean
  delete []face_id_vect;
  delete []mat_id_vect;
#ifdef PRINT_INFO
  std::cout << "\n GenBc(DA): boundary conditions defined by  bc_read(xp,normal,bc_value) \n" ;
#endif


  return;
  }


// ========================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::GenBc_loop (
  const int vvvb,         // 0-> volume   1-> boundary
  const int ndof_femv,  // number of elment nodes
  const int n_u_dofs,   // element quad dofs
  const int n_l_dofs,   // element linear dofs
  const int /*n_k_dofs*/,    // element konst dofs
  int face_id_vect[],   //bc from gambit
  int mat_id_vect[]    //bc from gambit

) { // ======================================================
  const int  n_dofs = _nvars[0] + _nvars[1] + _nvars[2];
  /// A) Set up: mesh,dof, bc
  //mesh
  const double * xyzgl  = _mgmesh._xyz;
  const int  offset    = _mgmesh._NoNodes[_NoLevels - 1];
//   double normal[DIMENSION];
  double xp[DIMENSION];
  double xxb_qnds[NDOF_FEM * DIMENSION];
  double xx_qnds[NDOF_FEM * DIMENSION];
  double normal[DIMENSION];
  int   el_neigh[NDOF_FEM];
  int   el_flag[NDOF_FEM];
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
//   int  *bc_Neu  =new int[n_dofs]; // element bc
//   int  *bc_value=new int[n_dofs]; // element bc
  int  bc_Neu[DIMENSION]; //new int[n_dofs]; // element bc
  int bc_value[DIMENSION];  //new int[n_dofs]; // element bc
  int face_id_node = 0;
  int mat_id_elem = 0;
  int ndof_lev = 0;
  const int  el_sides = _mgmesh._GeomEl._n_sides[0];

  for ( int isub = 0; isub < _mgmesh._n_subdom; ++isub ) {
      int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
      int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
      int delta = ielf - iel0;

      for ( int iel = 0; iel < delta; iel++ ) { // element loop
          _mgmesh.get_el_nod_conn ( 0, _NoLevels - 1, iel, el_conn, xx_qnds, isub );
          _mgmesh.get_el_neighbor ( el_sides, 0, _NoLevels - 1, iel, el_neigh, isub );

          for ( int i = 0; i < ndof_femv; i++ ) { // node lement loop
              const int k = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + i]; // global node

              // coordinates
              for ( int idim = 0; idim < DIMENSION; idim++ ) {
                  xp[idim] = xyzgl[k + idim * offset];
                  }

              if ( _node_dof[_NoLevels - 1][k] > -1 ) {
                  _bc[0][_node_dof[_NoLevels - 1][k]]  = -1000;
                  }

              }
          }
      }

  /// B) Element Loop to set  bc[] (which is a node vector)
  for ( int isub = 0; isub < _mgmesh._n_subdom; ++isub ) {
      int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
      int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
      int delta = ielf - iel0;

      for ( int iel = 0; iel < delta; iel++ ) { // element loop
          _mgmesh.get_el_nod_conn ( 0, _NoLevels - 1, iel, el_conn, xx_qnds, isub );
          _mgmesh.get_el_neighbor ( el_sides, 0, _NoLevels - 1, iel, el_neigh, isub );

          mat_id_elem = 0; // v==1 (boundari) no mat_id

          if ( vvvb == 0 ) {
              mat_id_elem = mat_id_vect[iel + iel0];
              }

          for ( int i = 0; i < ndof_femv; i++ ) { // node lement loop
              const int k = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + i]; // global node

//                 std::cout<<"i  "<<i<<"  k "<<k<<"   _node_dof[_NoLevels-1][k]  "<<_node_dof[_NoLevels-1][k]<<std::endl;
              // coordinates
              if ( _node_dof[_NoLevels - 1][k] > -1 ) {
                  for ( int idim = 0; idim < DIMENSION; idim++ ) {
                      xp[idim] = xyzgl[k + idim * offset];
                      }

                  face_id_node = face_id_vect[k];
//          for(int ivar=0; ivar<n_dofs; ivar++)  {
                  bc_value[0] = 1;
                  bc_Neu[0] = 11; // variable loop
                  // boundary (face_id_vect) and volume zones (mat_id_elem)

                  /// Calling the local point functions
                  bc_intern_read ( face_id_node, mat_id_elem, xp, bc_Neu, bc_value );

                  int old_val =  _bc[0][_node_dof[_NoLevels - 1][k]] ;
                  _bc[0][_node_dof[_NoLevels - 1][k]]  = ( _bc[0][_node_dof[_NoLevels - 1][k]] == -1000 ) ? bc_Neu[0] : ( _bc[0][_node_dof[_NoLevels - 1][k]] );

                  // sharing boundary nodes on the same element
                  for ( int ivar = 0; ivar < n_dofs; ivar++ )  {
                      int kdof = _node_dof[_NoLevels - 1][k + ivar * offset]; // kdof <-k
                      int number = _bc[1][kdof] / 10000;   // number of old count for double pt in an element

                      if ( abs ( number ) == 1 ) {
                          _bc[1][kdof] = _bc[1][kdof] - _bc[1][kdof] * 10000 / abs ( _bc[1][kdof] ); // if the count is 1 set 0 if >1 leave
                          }
                      }
                  }
              }  // i-loop


          // -----------------------------  Boundary ------------------------
          for ( int  iside = 0; iside < el_sides; iside++ ) {
              if ( el_neigh[iside] == -1 ) {
                  // setup boundary element -> connectivity+coordinates

//    int face_mid=;
                  for ( int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++ ) {
                      int lnode = _mgmesh._GeomEl._surf_top[lbnode + NDOF_FEMB * iside]; // local nodes
                      sur_toply[lbnode] = lnode;        // lbnode -> lnode
                      const int k = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + sur_toply[lbnode]];

                      if ( _node_dof[_NoLevels - 1][k] > -1 ) {
                          face_id_node = face_id_vect[k]; // boundary (face_id_vect)
//           elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn

                          for ( int idim = 0; idim < DIMENSION; idim++ ) {
                              xp[idim] = xxb_qnds[idim * NDOF_FEMB + lbnode] = xx_qnds[idim * NDOF_FEM + lnode];
                              } // coordinates

//           }
//           for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
//             const int k=_mgmesh._el_map[0][(iel+iel0)*ndof_femv+ sur_toply[lbnode]];
//             face_id_node=face_id_vect[k]; // boundary (face_id_vect)
                          bc_value[0] = 1;
                          bc_Neu[0] = 11;
                          bc_read ( face_id_node, mat_id_elem, xp, bc_Neu, bc_value );
                          _bc[0][_node_dof[_NoLevels - 1][k]] = bc_Neu[0];
                          }
                      }

                  // normal
                  _fe[2]->normal_g ( xxb_qnds, normal );
                  int dir_maxnormal = ( fabs ( normal[0] ) > fabs ( normal[1] ) ) ? 0 : 1 ;
                  dir_maxnormal = ( fabs ( normal[dir_maxnormal] ) > fabs ( normal[DIMENSION - 1] ) ) ? dir_maxnormal : DIMENSION - 1;
                  // computation of k_el and face_id_el
                  int k_el = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + sur_toply[NDOF_FEMB - 1]]; // global node
                  int k_el_dof = _node_dof[_NoLevels - 1][k_el];
                  int bc_el = ( int ) _bc[0][k_el_dof] % 100;
                  int score_old = 1000;
//
//
// //         if(ndof_femv-NDOF_P==3) { // tet 10 only
// //           for(int ik=0; ik<3; ik++) {
// //             int score=0;
// //             int  bc_id=bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+ndof_femv-1-ik]]];
// //             score +=(bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0)*ndof_femv+2-ik]]]    == bc_id)?1:0;
// //             score +=(bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+(3-ik)%3]]]== bc_id)?1:0;
// //             if(score<score_old) {
// //               score_old=score; k_el=_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+ndof_femv-1-ik];
// //               bc_el=bc[0][ _node_dof[_NoLevels-1][k_el]];
// //             }
// //           }
// //         }
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++
                  const int k_face = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + sur_toply[NDOF_FEMB - 1]]; // global node
                  int k0dof =  _node_dof[_NoLevels - 1][k_face];
                  int face_id_mid_face = face_id_vect[k0dof];
                  int bc_face = ( int ) _bc[0][k0dof];
                  bc_face = bc_face % 100;

                  for ( int i = 0; i < NDOF_FEMB; i++ ) { // node lement loop
                      const int k = _mgmesh._el_map[0][ ( iel + iel0 ) * ndof_femv + sur_toply[i]]; // global node
                      int k0dof =  _node_dof[_NoLevels - 1][k];

                      if ( k0dof > -1 ) {
                          face_id_node = face_id_vect[k0dof];
                          int bc_id = ( int ) _bc[0][k0dof]; // label surface pt 00
                          int mynormal = dir_maxnormal; // dir_normal=geometrical normal

                          int sign = ( bc_id == 0 ) ? 1 : ( bc_id / ( abs ( bc_id ) ) );

                          if ( ( bc_id % 1000 ) / 100 > 0 ) {
                              mynormal = ( bc_id % 1000 ) / 100 - 1;
                              }

                          if ( ( bc_id % 1000 ) / 100 > 3 ) { //   pts label x00 ----------------------------------------------------------
                              std::cout << "\n ... single pressure pts .. \n";

                              if ( i < NDOF_PB ) {
                                  _bc[1][_node_dof[_NoLevels - 1][k + DIMENSION * offset]] = 4;
                                  }

                              mynormal = abs ( ( bc_id % 1000 ) / 100 ) % 4; // force the normal
                              }

                          bc_id = bc_id % 100; // pts label x00 forced normal

                          //   pts label 00 ----------------------------------------------------------
                          // set the local boundary conditions into global vector bc[]
                          for ( int ivar = 0; ivar < _nvars[2]; ivar++ ) { // quad el --
                              int kdof = _node_dof[_NoLevels - 1][k + ( ivar ) * offset];
                              int number = abs ( _bc[1][kdof] / 10000 ) + 1;
                              _bc[1][kdof] += sign * 10000; // updating number of common nodes

                              if ( abs ( _bc[1][kdof] ) < 10000 || bc_id == bc_face ) {
                                  _bc[1][kdof] = sign * ( abs ( bc_id ) + ( mynormal + 1 ) * 1000 + number * 10000 );
                                  }
                              } // ----------------------------------------------------------------
                          }
                      } // loop i +++++++++++++++++++++++++++++++++++++++++++++
                  } // iside -1
              }  // -----------------------------  End Boundary -------------------------------------
          } // end of element loop

      ndof_lev += delta;
      } // i-sub


//   if(vb==1) {
//
//     ndof_lev =0;
//     /// B) Element Loop to set  bc[] (which is a node vector)
//     for(int isub=0; isub<_mgmesh._n_subdom; ++isub) {
//       int iel0=_mgmesh._off_el[vb][_NoLevels-1+_NoLevels*isub];
//       int ielf=_mgmesh._off_el[vb][_NoLevels-1+_NoLevels*isub+1];
//       int delta =ielf-iel0;
//
//       for(int iel=0; iel <delta; iel++) { // element loop
//
//         for(int i=0; i< ndof_femv; i++)  // node lement loop
//           for(int idim=0; idim< DIMENSION; idim++) { xxb_qnds[i+ndof_femv*idim]=  xyzgl[_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+i]+idim*offset]; }
//
//
//         // normal
//         _fe[2]->normal_g(xxb_qnds,normal);
//         int dir_maxnormal = (fabs(normal[0])>fabs(normal[1]))?0:1 ;
//         dir_maxnormal= (fabs(normal[dir_maxnormal])>fabs(normal[DIMENSION-1]))? dir_maxnormal:DIMENSION-1;
//
//         // computation of k_el and face_id_el
//         int k_el=_mgmesh._el_map[1][(iel+iel0)*ndof_femv+ndof_femv-1]; // global node
//         int k_el_dof=_node_dof[_NoLevels-1][k_el];
//         int bc_el=(int)bc[0][k_el_dof];
//         int score_old=1000;
//
//
//         if(ndof_femv-NDOF_P==3) { // tet 10 only
//           for(int ik=0; ik<3; ik++) {
//             int score=0;
//             int  bc_id=bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+ndof_femv-1-ik]]];
//             score +=(bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0)*ndof_femv+2-ik]]]    == bc_id)?1:0;
//             score +=(bc[0][ _node_dof[_NoLevels-1][_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+(3-ik)%3]]]== bc_id)?1:0;
//             if(score<score_old) {
//               score_old=score; k_el=_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+ndof_femv-1-ik];
//               bc_el=bc[0][ _node_dof[_NoLevels-1][k_el]];
//             }
//           }
//         }
//
//
//         for(int i=0; i< ndof_femv; i++) { // node lement loop
//           const int k=_mgmesh._el_map[vb][(iel+iel0) *ndof_femv+i]; // global node
//           int k0dof=  _node_dof[_NoLevels-1][k];
//           face_id_node=face_id_vect[k0dof];
//           int bc_id=(int)bc[0][k0dof]; // label surface pt 00
//           int bc_id_pts= bc_id/100;  // pts label x00
//
//           if(bc_id > 99) { //   pts label x00 ----------------------------------------------------------
//             std::cout << "\n ... single pts  .. \n";
//             bc_id= bc_id%100;
//             if(bc_id_pts ==4) if(i<NDOF_P) { bc[1][_node_dof[_NoLevels-1][k+DIMENSION* offset]] =4; }
//             for(int ivar=0; ivar<_nvars[2]; ivar++) {// quad el --
//               int kdof=_node_dof[_NoLevels-1][k+(ivar)* offset];
//               if(ivar+_dir == (bc_id_pts-1)%DIMENSION) { bc[1][kdof] =(bc_id/10)+(ivar+_dir+1)*10; }  // normal= bc_id_pts-1
//               else {   bc[1][kdof] =bc_id%10; }  // tang
//             }
//
//           }  // end ------(bc_id > 99)
//           else { //   pts label 00 ----------------------------------------------------------
//
//
// //           bc_read(face_id_node,mat_id_elem,xp,bc_Neu,bc_value);
//             // set the local boundary conditions into global vector bc[]
//             if(bc_id==bc_el) {
//               if(i<n_u_dofs && bc_id==bc_el)   {
//                 for(int ivar=0; ivar<_nvars[2]; ivar++) {// quad el --
//                   int kdof=_node_dof[_NoLevels-1][k+(ivar)* offset]; int tdof=  _node_dof[_NoLevels-1][k];
//                   if(ivar+_dir == dir_maxnormal)  {  bc[1][kdof] =(bc_el/10)+(dir_maxnormal+1)*10;}
//                   else {   bc[1][kdof] =bc_el%10; }
//                 } // ----------------------------------------------------------------
//               }
//             }
// //             _node_dof[_NoLevels-1][i+(ivar+_dir)*offset]]
// //             if(i<n_l_dofs) { // Set the linear fields -----------------------------
// //               for(int ivarp=_nvars[2]; ivarp<(_nvars[1]+_nvars[2]); ivarp++) {
// //                 int kdof= _node_dof[_NoLevels-1][k+ (ivarp) * offset]; int tdof=  _node_dof[_NoLevels-1][k];
// //                 bc[1][kdof] =(int)(bc[0][tdof]/10)+10;
// //               }
// //             }// end if(i<n_l ------------------------------------------------------
// //           }// bc_id==bc_el
// //
// //             if(i == ndof_femv-1) { // Set the constant piecewise fields---
// //           for(int ivarp=_nvars[2]+_nvars[1]; ivarp<(_nvars[0]+_nvars[2]+_nvars[1]); ivarp++) {
// //             for(int idof=0; idof<_el_dof[0]; idof++) {
// //               int kdof= _node_dof[_NoLevels-1][idof+(iel+ndof_lev)*_el_dof[0]+ivarp*offset];
// //               int tdof=  _node_dof[_NoLevels-1][k];
// //               bc[1][kdof] =(int)(bc[0][tdof]/10)+10;
// // //               bc[0][kdof]=bc_Neu[ivarp];
// // //               bc[1][kdof]=bc_value[ivarp];
// //             }
// //           }
//           } // end  if(i
//
//         } // loop i
//
//
//
//
//
//
//       } // end of element loop
//       ndof_lev +=delta;
//     }
//   }


// clean
// //   delete []bc_value;  delete []bc_Neu;
  return;
  }


// ========================================
/// This function  defines the boundary conditions for the system:
void MGSolDA::bc_intern_read (
  int /*face_id_node*/,  ///<  face identity           (in)
  int  /*mat_flag*/,     ///<  volume identity         (in)
  double /*xp*/[],       ///< xp[] node coordinates    (in)
  int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
  int bc_flag[]          ///< boundary condition flag  (out)
) { // ===================================
  /// Default: all Neumann
//   for(int ivar=0; ivar<_n_vars; ivar++) {
  bc_flag[0] = 1;
  bc_Neum[0] = 11;

//   }
  return;
  }
// ===========================================================
/// This function initializes the system degrees of freedom (dof)
void MGSolDA::init_dof (
  const int Level  // Level
) { // ========================================================

  // Set up from mesh -----------------------------
  const int n_subdom = _mgmesh._n_subdom;
  const int n_nodes  = _mgmesh._NoNodes[Level];
  const int n_elem   = _mgmesh._NoElements[0][Level]; //0=volume

  const  int offset = _mgmesh._NoNodes[_NoLevels - 1];
  const int * off_nd_q = _mgmesh._off_nd[0];
  const int * off_nd_l = _mgmesh._off_nd[1];
  const int * off_el  = _mgmesh._off_el[0]; //volume elements

// number of total dofs
  int n_nodes_l = _mgmesh._NoNodes[_mgmesh._NoFamFEM * _NoLevels];

  if ( Level > 0 ) {
      n_nodes_l = _mgmesh._NoNodes[Level - 1];
      }

  _Dim[Level] = _nvars[0] * n_elem * _el_dof[0] + _nvars[1] * n_nodes_l + _nvars[2] * n_nodes;

#ifdef PRINT_INFO
  std::cout << Level << " node " << n_nodes << " dof " << _Dim[Level] << " pres " << n_nodes_l << std::endl;
#endif

// Set up boundary conditions ++++++++++++++++++++++++++++
  if ( Level == _NoLevels - 1 ) {
      _bc[0] = new int[_Dim[Level]];
      _bc[1] = new int[_Dim[Level]];

      for ( int k1 = 0; k1 < _Dim[Level]; k1++ ) {
          _bc[0][k1] = 1;
          _bc[1][k1] = 1;
          }
      }

  // construction dof node vector(_node_dof) +++++++++++++++++++++
  _node_dof[Level] = new int[_n_vars * offset];

  for ( int k1 = 0; k1 < offset * _n_vars; k1++ ) {
      _node_dof[Level][k1] = -1;
      }

  int count = 0;
  int ndof_lev = 0;

  for ( int isubdom = 0; isubdom < n_subdom; isubdom++ ) {
      // quadratic -----------------------------------
      int off_proc = isubdom * _NoLevels;

      for ( int ivar = 0; ivar < _nvars[2]; ivar++ ) {
          for ( int k1 = off_nd_q[off_proc]; k1 < off_nd_q[off_proc + Level + 1]; k1++ ) {
              _node_dof[Level][k1 + ivar * offset] = count;
              count++;
              }
          }

      // linear -----------------------------------
      for ( int ivar = 0; ivar < _nvars[1]; ivar++ ) {
          for ( int k1 = off_nd_q[off_proc];
                k1 < off_nd_q[off_proc] + off_nd_l[Level + 1 + off_proc] - off_nd_l[off_proc]; k1++ ) {
              _node_dof[Level][k1 + ( _nvars[2] + ivar ) *offset] = count;
              count++;
              }
          }

      // konstant polynomial of order _el_dof[0] -----------------------------------
      int delta_el = off_el[off_proc + Level + 1] - off_el[off_proc + Level];

      for ( int ivar = 0; ivar < _nvars[0]; ivar++ ) {
          for ( int iel = 0; iel < delta_el; iel++ ) {
              for ( int idof = 0; idof < _el_dof[0]; idof++ ) {
                  _node_dof[Level][idof + ( iel + ndof_lev ) *_el_dof[0] + ( _nvars[2] + _nvars[1] + ivar ) *offset] = count;
                  count++;
                  }
              }
          }

      ndof_lev += delta_el;
      }

#ifdef PRINT_INFO
  std::cout << "MGSol::init_dof(D)   Level= " << Level  << std::endl;
#endif
  return;
  }


// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolDA::get_el (
  const int Level,       // level
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  std::vector<int>  &  el_dof_indices, // element DOFs->
  int   bc_dofs[][NDOF_FEM],        // element boundary cond flags ->
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for ( int
        id = 0; id < el_nds; id++ )    {
      // quadratic -------------------------------------------------
      for ( int
            ivar = ivar0; ivar < ivar0 + nvars; ivar++ ) { //ivarq is like idim
          const int
          indx_loc = id + ivar * NDOF_FEM;     // local (element) index
          const int indx_glob = el_conn[id] + ivar * offset; // global (mesh) index
          const int  kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level

          el_dof_indices[indx_loc] = _node_dof[Level][indx_glob];    //from mesh to dof
          bc_dofs[0][indx_loc]       = _bc[0][kdof_top];                    // element bc
          bc_dofs[1][indx_loc]       = _bc[1][kdof_top];                    // element bc
          uold[indx_loc]          = ( *x_old[_NoLevels - 1] ) ( kdof_top ); // element sol
          } // end quadratic ------------------------------------------------
      }

  return;
  }

// ==============================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void  MGSolDA::get_el_dof_bc (
  const int Level,       // level
  const int iel,         // element subdomain
//   const int nvars[],       // # of variables to get  <-
  const int el_nds[],      // # of element nodes for this variable  <-
  const int el_conn[],   // connectivity <-
  const int  offset,      // offset for connectivity <-
  std::vector<int>  &  el_dof_indices, // element connectivity ->
  int  bc_vol[],        // element boundary cond flags ->
  int  bc_bd[]        // element boundary cond flags ->
)  const { // ==============================================================

  for ( int id = 0; id < NDOF_FEM; id++ )  {
      // quadratic -------------------------------------------------
      if ( id < el_nds[2] )  for ( int  ivar = 0; ivar < _nvars[2]; ivar++ ) { //ivarq is like idim
            const int indx_loc = id + ivar * NDOF_FEM;
            const int indx_loc_ql = id + ivar * el_nds[2];
            const int indx_glob = el_conn[id] + ivar * offset;
            const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level

            el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
            bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
            bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
            } // end quadratic ------------------------------------------------

//     // linear -----------------------------
      if ( id < el_nds[1] )    for ( int ivar = 0; ivar < _nvars[1]; ivar++ ) { //ivarq is like idim
//        const int  indx_loc_l = id +ivar*el_nds[1];
            const int indx_loc = id + ( ivar + _nvars[2] ) * NDOF_FEM;
            const int indx_loc_ql = id + ivar * el_nds[1] + _nvars[2] * el_nds[2];
            const int indx_glob = el_conn[id] + ( ivar + _nvars[2] ) * offset;
            const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level

            el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
            bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
            bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
            } // end quadratic ------------------------------------------------


      //     // piecewise -----------------------------
      if ( id < el_nds[0] )    for ( int ivar = 0; ivar < _nvars[0]; ivar++ ) { //ivarq is like idim
//        const int  indx_loc_l = id +ivar*el_nds[1];
            const int indx_loc = id + ( ivar + _nvars[2] + _nvars[1] ) * NDOF_FEM;
            const int indx_loc_ql = id + ivar * el_nds[0] + _nvars[2] * el_nds[2] + _nvars[1] * el_nds[1];
            const int indx_glob = id + iel * el_nds[0] + ( ivar + _nvars[2] + _nvars[1] ) * offset;
            const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level

            el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
            bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
            bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
            } // end piecewise ------------------------------------------------
      }

  return;
  }

// =================================================================================
/// This function interpolates a vector field and its derivative over the fem element
void  MGSolDA::interp_el_gdx (
  double uold_b[],       // node values <-
  const int ivar0,       // init variable  <-
  const int nvars,       // # of variables  <-
  const double dphi[],   // derivatives of the shape functions  <-
  const int n_shape,     // # of shape functions  <-
  double uold_dx[]       // interpolated derivatives ->
)  const {  // =================================================================
  // All data vectors are NDOF_FEM long

  // variable loop

  // dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz

  for ( int ivar = 0; ivar < nvars; ivar++ ) {
      for ( int jdim = 0; jdim < DIMENSION; jdim++ ) {
          uold_dx[ivar * DIMENSION + jdim] = 0.; // set zero
          }

      // interpolation with shape functions
      for ( int eln = 0; eln < n_shape; eln++ ) {
          for ( int jdim = 0; jdim < DIMENSION; jdim++ ) {
              uold_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * uold_b[eln + ( ivar + ivar0 ) * NDOF_FEM];
              }
          }
      }

  return;
  }

// =================================================================================
/// This function interpolates a vector field and its second derivatives over the fem element
void  MGSolDA::interp_el_gddx (
  double uold_b[], // node values <-
  const int ivar0,      // init variable  <-
  const int nvars,      // # of variables  <-
  const double dphi[],   // derivatives of the shape functions  <-
  const int n_shape,    // # of shape functions  <-
  double uold_dx[]       // interpolated derivatives ->
)  const {  // =================================================================
  // All data vectors are NDOF_FEM long

  // variable loop
  for ( int ivar = 0; ivar < nvars; ivar++ ) {
      for ( int jdim = 0; jdim < DIMENSION * DIMENSION; jdim++ ) {
          uold_dx[ivar * DIMENSION * DIMENSION + jdim] = 0.; // set zero
          }

      // interpolation with shape functions
      for ( int eln = 0; eln < n_shape; eln++ ) {
          for ( int jdim = 0; jdim < DIMENSION * DIMENSION; jdim++ ) {
              uold_dx[ivar * DIMENSION * DIMENSION + jdim] += dphi[eln * DIMENSION * DIMENSION + jdim] * uold_b[eln + ( ivar + ivar0 ) * NDOF_FEM];
              }
          }
      }

  return;
  }


///
void  MGSolDA::compute_jac (
  const int
  j,
  const int
  idim,
  double uold_b[], // node values <-
  const int
  nvars,      // # of variables  <-
  const double phi[],    // shape functions  <-
  const double dphi[],   // derivatives of the shape functions  <-
  const int
  n_shape,    // # of shape functions  <-
  double u_forw[],         // interpolated function ->
  double u_back[],         // interpolated function ->
  double u_forw_dx[],       // interpolated derivatives ->
  double u_back_dx[]       // interpolated derivatives ->
)  const {  // =================================================================
  // All data vectors are NDOF_FEM long

  const double alfa = 1.e-08;

  // variable loop
  for ( int
        ivar = 0; ivar < nvars; ivar++ ) {
      u_forw[ivar] = 0.; // set zero
      u_back[ivar] = 0.; // set zero

      for ( int
            jdim = 0; jdim < DIMENSION; jdim++ ) {
          u_forw_dx[ivar * DIMENSION + jdim] = 0.; // set zero
          u_back_dx[ivar * DIMENSION + jdim] = 0.; // set zero
          }

      // interpolation with shape functions
      for ( int
            eln = 0; eln < n_shape; eln++ )   {
          const int
          indx  = eln + ivar * NDOF_FEM;

          if ( indx == j + idim * NDOF_FEM ) {
              u_forw[ivar] += phi[eln] * ( uold_b[indx] + alfa );
              u_back[ivar] += phi[eln] * ( uold_b[indx] - alfa );
              }
          else {
              u_forw[ivar] += phi[eln] * uold_b[indx];
              u_back[ivar] += phi[eln] * uold_b[indx];
              }

          for ( int
                jdim = 0; jdim < DIMENSION; jdim++ ) {

              if ( indx == j + idim * NDOF_FEM ) {
                  u_forw_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * ( uold_b[indx] + alfa );
                  u_back_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * ( uold_b[indx] - alfa );
                  }
              else {
                  u_forw_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
                  u_back_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
                  }
              }
          }
      }

  return;
  }

// =================================================================================
/// This function interpolates a vector field over the fem element
void  MGSolDA::interp_el_sol (
  const double uold_b[],     // node values <-
  const int  ivar0,          // init variable  <-
  const int  nvars,          // # of variables  <-
  const double phi[],        // shape functions  <-
  const int  n_shape,        // # of shape functions  <-
  double uold[]              // interpolated function ->
)  const { // =======================================
  for ( int ivar = ivar0; ivar < ivar0 + nvars; ivar++ ) {
      uold[ivar] = 0.;

      for ( int eln = 0; eln < n_shape; eln++ ) {
          uold[ivar] += phi[eln] * uold_b[eln + ivar * n_shape];
          }
      }

  return;
  }



// =================================================================================
/// This function interpolates a vector field over the fem element
void  MGSolDA::interp_el_bd_sol (
  const double uold_b[],     // node values <-
  const int sur_tpgly[],     //surface nodes topology <-
  const int el_ndof,     //surface nodes topology <-
  const int ivar0,          // init variable  <-
  const int nvars,          // # of variables  <-
  const double phi[],        // shape functions  <-
  const int n_shape,        // # of shape functions  <-
  double uold[]             // interpolated function ->
)  const { // =======================================
  for ( int ivar = ivar0; ivar < ivar0 + nvars; ivar++ ) {
      uold[ivar] = 0.;

      for ( int eln = 0; eln < n_shape; eln++ )  {
          uold[ivar] += phi[eln] * uold_b[sur_tpgly[eln] + ivar * el_ndof];
          }
      }

  return;
  }



// =================================================================================
/// This function interpolates a vector field and its derivative over a boundary
void  MGSolDA::interp_el_bd_gdx (
  const double uold_b[],       // node values <-
  const int sur_tpgly[],     //surface nodes topology <-
  const int el_ndof,     //surface nodes topology <-
  const int ivar0,       // init variable  <-
  const int nvars,       // # of variables  <-
  const double dphi[],   // derivatives of the shape functions  <-
  const int n_shape,     // # of shape functions  <-
  double uold_dx[]       // interpolated derivatives ->
)  const {  // =================================================================
  // variable loop

  int dim_b = DIMENSION ;

  for ( int ivar = 0; ivar < nvars; ivar++ ) {      // loop over variables (u,v,w,...)
      for ( int jdim = 0; jdim < dim_b; jdim++ ) {
          uold_dx[ivar * dim_b + jdim] = 0.; // set zero
          }

      // interpolation with shape functions
      for ( int eln = 0; eln < n_shape; eln++ ) {
          for ( int jdim = 0; jdim < dim_b; jdim++ ) { // loop over directions (x,y,z,...)
              uold_dx[ivar * dim_b + jdim] += dphi[eln + jdim * n_shape] * uold_b[sur_tpgly[eln] + ( ivar + ivar0 ) * el_ndof];
              }
          }
      }

  return;
  }



// =============================================
/// Print xml attrib
void MGSolDA::print_xml_attrib (
  std::ofstream & out, //  file xdmf
  int nodes,
  int nelems,
  std::string file_name
) const { // ================================

  for ( int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++ )   {
      std::string var_name = _var_names[ivar];
      out << "<Attribute Name=\"" << var_name << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
      out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
          << nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
      out << file_name
//       femus_dir << "/" << output_dir << basesol << "."
//       << setw(ndigits) << setfill('0') << t_step << ".h5"
          << ":" << var_name << "\n";
      out << "</DataItem>\n" << "</Attribute>";
      }

//   if(DIMENSION==3) nelems*=2;
  for ( int ivar = 0; ivar < _nvars[0]; ivar++ )   {
      std::string var_name = _var_names[ivar + _nvars[2] + _nvars[1]];
      out << "<Attribute Name=\"" << var_name << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
      out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
          << nelems << "  " << 1 << "\" Format=\"HDF\">  \n";
      out << file_name
//       femus_dir << "/" << output_dir << basesol << "."
//       << setw(ndigits) << setfill('0') << t_step << ".h5"
          << ":" << var_name << "\n";
      out << "</DataItem>\n" << "</Attribute>";
      }

  return;
  }
// =============================================
// ============================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_u (
  std::string namefile, const int Level
) { // ===============================================
  //  set up
  const int n_nodes = _mgmesh._NoNodes[Level];
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  double * sol = new double[n_nodes + 1];

  hid_t file_id = H5Fopen ( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t  dimsf[2];
  dimsf[0] = n_nodes;
  dimsf[1] = 1;
  // loop over all systems ---------------------------

  int SolToPrint = _NumRestartSol;
  std::string SolSuffix[3];
  SolSuffix[0] = "";
  SolSuffix[1] = "_old";
  SolSuffix[2] = "_oold";
  
  // print quad -------------------------------------
  for ( int ivar = 0; ivar < _nvars[2]; ivar++ )        {
      for( int timeStep = 0; timeStep < SolToPrint; timeStep ++ ){
      std::string var_name = _var_names[ivar] + SolSuffix[timeStep];

      if(timeStep == 0)
        for ( int i = 0; i < n_nodes; i++ ) 
          sol[i]  = ( *x_old[Level] ) ( _node_dof[Level][i + ivar * offset] ) * _refvalue[ivar];
      else if(timeStep == 1)
        for ( int i = 0; i < n_nodes; i++ ) 
          sol[i]  = ( *x_oold[Level] ) ( _node_dof[Level][i + ivar * offset] ) * _refvalue[ivar];    
      else if(timeStep == 2)
        for ( int i = 0; i < n_nodes; i++ ) 
          sol[i]  = ( *x_ooold[Level] ) ( _node_dof[Level][i + ivar * offset] ) * _refvalue[ivar];      

      _mgutils.print_Dhdf5 ( file_id, var_name, dimsf, sol );
      }
      }

  // print linear -----------------------------------
  double * sol_c = new double[_mgmesh._GeomEl.n_l[0]];
//   double *Prol=_mgmesh._GeomEl.Prol;
//   const int  lev_l=(Level >0) ?Level-1:_NoLevels;

  for ( int ivar = _nvars[2]; ivar < _nvars[2] + _nvars[1]; ivar++ )   {

      for ( int timeStep = 0; timeStep < SolToPrint; timeStep ++ ) {
          std::string var_name = _var_names[ivar] + SolSuffix[timeStep];

          //  2bB element interpolation over the fine mesh -----------------------
          for ( int iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
              int start = _mgmesh._off_el[0][_NoLevels * iproc + Level];
              int end   = _mgmesh._off_el[0][Level + _NoLevels * iproc + 1];
              for ( int iel = start; iel < end; iel++ ) {
                  int indx = iel * NDOF_FEM;
                  // element iel ----------------------------
                  // vertices
                  for ( int in = 0; in < NDOF_P; in++ ) {
                      int gl_i = _mgmesh._el_map[0][indx + in];
                      double val = 0.;
                      if(timeStep == 0) val = ( *x_old[Level] ) ( _node_dof[Level][gl_i + ivar * offset] );
                      if(timeStep == 1) val = ( *x_oold[Level] ) ( _node_dof[Level][gl_i + ivar * offset] );
                      if(timeStep == 2) val = ( *x_ooold[Level] ) ( _node_dof[Level][gl_i + ivar * offset] );
                      
                      sol_c[in] = val * _refvalue[ivar];
                      }

                  for ( int in = 0; in < NDOF_FEM; in++ ) { // mid-points
                      double sum = 0;

                      for ( int jn = 0; jn < NDOF_P; jn++ ) {
                          sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn];
                          }

                      sol[_mgmesh._el_map[0][indx + in]] = sum;
                      }
                  } // ---- end iel -------------------------
              } // 2bB end interpolation over the fine mesh ------------------------

          _mgutils.print_Dhdf5 ( file_id, var_name, dimsf, sol );
          }
      } // ivar


  // print constant -----------------------------------
  int size = 0;

  for ( int iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
      int delta_c = _mgmesh._off_el[0][Level + 1 + _NoLevels * iproc] - _mgmesh._off_el[0][_NoLevels * iproc + Level];
      size += delta_c;
      }

  size *= NSUBDOM;

  double * sol_p = new double[size];

  for ( int ivar = 0; ivar < _nvars[0]; ivar++ )   {
      std::string var_name = _var_names[ivar + _nvars[2] + _nvars[1]];

      //  2bB element interpolation over the fine mesh -----------------------
      int eldof_lev = 0;

      for ( int iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
          int delta = _mgmesh._off_el[0][Level + 1 + _NoLevels * iproc] - _mgmesh._off_el[0][_NoLevels * iproc + Level];

          for ( int iel = 0; iel < delta; iel++ ) {
              int indx = iel + eldof_lev; // print only the first dof

              double val = ( *x_old[Level] ) ( _node_dof[Level][indx * _el_dof[0] + ( ivar + _nvars[2] + _nvars[1] ) * offset] ) * _refvalue[_nvars[0]];

              for ( int isubcell = 0; isubcell < NSUBDOM; isubcell++ ) {
                  sol_p[indx * NSUBDOM + isubcell] = val * _refvalue[DIMENSION];
                  }

              } // ---- end iel -------------------------

          eldof_lev += delta;
          } // 2bB end interpolation over the fine mesh ------------------------

      dimsf[0] = eldof_lev * NSUBDOM;
      dimsf[1] = 1;
      _mgutils.print_Dhdf5 ( file_id, var_name, dimsf, sol_p );
      } // ivar




  // clean and close -----------------

  H5Fclose ( file_id );
  delete []sol_c;
  delete []sol;
  delete []sol_p;

  return;
  }

// void MGSolDA::set_ext_fields(const std::vector<FIELDS> &pbName) {
//    int a =1;
//    return;
// }
// =========================================
void  MGSolDA::set_xooold2x() {
  for ( int Level = 1; Level <= _NoLevels; Level++ ) {
/// A. Setup
      const int offset = _mgmesh._NoNodes[Level - 1]; // fine level # of nodes
      const int pie_offset = ( NSUBDOM ) * ( _mgmesh._NoElements[0][Level - 1] ); // fine level # of nodes

      // field
      // file to read
      // reading loop over system varables
      for ( int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++ ) {
          int  el_nds = NDOF_FEM;

          if ( ivar >= _nvars[2] ) {
              el_nds = NDOF_P;  // quad and linear
              }

          // reading ivar param
          double Irefval = 1. / _refvalue[ivar]; // units

          // storing  ivar variables (in parallell)
          for ( int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + Level] -
                _mgmesh._off_el[0][_iproc * _NoLevels + Level - 1]; iel++ ) {
              int  elem_gidx = ( iel + _mgmesh._off_el[0][_iproc * _NoLevels + Level - 1] ) * NDOF_FEM;

              for ( int  i = 0; i < el_nds; i++ ) { // linear and quad
                  int k = _mgmesh._el_map[0][elem_gidx + i]; // the global node
                  const double value = ( *x_ooold[Level - 1] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                  x[Level - 1]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                  }
              }
          }

      int ndof_lev = 0;

      for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
          int delta = _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
          ndof_lev += delta;
          }

      /// D. delocalization and clean
      x[Level - 1]->localize ( *x_old[Level - 1] );
      x[Level - 1]->localize ( *x_oold[Level - 1] );
      }

  return;
  }
// ==========================================================================
/// This function copies values from "vec_from" to "vec_to" vectors.
/// Available vectors are ordered as: x->0, x_old->1, x_oold->2, x_ooold->3,
/// x_nonl->4, disp->5, disp_old->6,  disp_oold->7
// ==========================================================================
void  MGSolDA::set_vector ( const int & vec_from, const int & vec_to ) {

  for ( int Level = 0; Level <= _NoLevels - 1; Level++ ) {
      const int offset = _mgmesh._NoNodes[Level]; // fine level # of nodes

      for ( int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++ ) {
          int  el_nds = NDOF_FEM;

          if ( ivar >= _nvars[2] ) {
              el_nds = NDOF_P;  // quad and linear
              }

          for ( int  iproc = 0; iproc < _mgmesh._n_subdom; iproc++ )
            for ( int iel = 0; iel < _mgmesh._off_el[0][iproc * _NoLevels + Level + 1] -
                  _mgmesh._off_el[0][iproc * _NoLevels + Level]; iel++ ) {
                int  elem_gidx = ( iel + _mgmesh._off_el[0][iproc * _NoLevels + Level] ) * NDOF_FEM;

                for ( int  i = 0; i < el_nds; i++ ) { // linear and quad
                    int k = _mgmesh._el_map[0][elem_gidx + i]; // the global node
                    double value = 0.;

                    switch ( vec_from ) {
                      case 1:
                        value = ( *x_old[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 2:
                        value = ( *x_oold[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 3:
                        value = ( *x_ooold[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 4:
                        value = ( *x_nonl[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 5:
                        value = ( *disp[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 6:
                        value = ( *disp_old[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 7:
                        value = ( *disp_oold[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 8:
                        value = ( *disp_ooold[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      case 9:
                        value = ( *x_oooold[Level] ) ( _node_dof[_NoLevels - 1][k + ivar * offset] );
                        break;

                      default:
                        cout << "Incorrect vec_from number in set_uoold function" << endl;
                        break;
                        }

                    switch ( vec_to ) {
                      case 1:
                        x_old[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 2:
                        x_oold[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 3:
                        x_ooold[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 4:
                        x_nonl[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 5:
                        disp[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 6:
                        disp_old[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 7:
                        disp_oold[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 8:
                        disp_ooold[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      case 9:
                        x_oooold[Level]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], value ); // set the field
                        break;

                      default:
                        cout << "Incorrect vec_to number in set_uoold function" << endl;
                        break;
                        }
                    }
                }
          }

      if ( vec_from == 1 && vec_to == 3 ) {
          x[Level]->localize ( *x_oold[Level] );    //for backward compatibility
          }
      }

  int ndof_lev = 0;

  for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
      int delta = _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
      ndof_lev += delta;
      }

  return;
  }

#ifdef HAVE_MED
// ============================================================
/// This function prints the controlled domain in med format
void MGSolDA::print_weight_med (
  std::string namefile, const int /* Level1*/
) { // ===============================================

  int n_nodes = _mgmesh._NoNodes[_NoLevels - 1];
  int n_elements = _mgmesh._NoElements[0][_NoLevels - 1];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * ( _mgmesh._n_subdom - 1 )];
  int nodes_el = _mgmesh._type_FEM[0];

#if ELTYPE==27
#if DIMENSION==3
  const unsigned int nodesinv[]   = {7, 4, 5, 6, 3, 0, 1, 2, 19, 16, 17, 18, 11, 8, 9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};
  const unsigned int nodesinvbd[] = {3, 0, 1, 2, 7, 4, 5, 6, 8};
#else
  const unsigned int nodesinv[]   = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const unsigned int nodesinvbd[] = {0, 1, 2};
#endif
#endif
#if ELTYPE==10
#if DIMENSION==3
  const unsigned int nodesinv[]   = {2, 3, 1, 0, 9, 8, 5, 6, 7, 4};
  const unsigned int nodesinvbd[] = {1, 2, 0, 4, 5, 3};
#else
  const unsigned int nodesinv[]   = {1, 2, 0, 4, 5, 3};
  const unsigned int nodesinvbd[] = {0, 2, 1};
#endif
#endif

  double * coord;
  coord = new double[n_nodes * 3];

  for ( int i = 0; i < n_nodes; i++ ) {
      coord[i * 3 + 1] = 0.;
      coord[i * 3 + 2] = 0.;

      for ( int idim = 0; idim < _mgmesh._dim; idim++ ) {
          coord[i * 3 + idim] = _mgmesh._xyz[i + idim * n_nodes];
          }
      }

  int  Level = _mgmesh._NoLevels - 1;
  int icount = 0;
  int * conn;
  conn = new int [n_elements * nodes_el];

  for ( int  iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
      for ( int el = _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level];
            el < _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level + 1]; el++ ) {
          for ( int  i = 0; i < nodes_el; i++ ) {
              conn[icount] = _mgmesh._el_map[0][el * nodes_el + nodesinv[i]];
              icount++;
              }
          }
      }

  MEDCoupling::MEDCouplingUMesh * mesh = MEDCoupling::MEDCouplingUMesh::New ( "Mesh_1", _mgmesh._dim );
  mesh->allocateCells ( n_elements );

  for ( int  i = 0; i < n_elements; i++ ) {
      mesh->insertNextCell ( MED_EL_TYPE, nodes_el, conn + i * nodes_el );
      }

  mesh->finishInsertingCells();

  MEDCoupling::DataArrayDouble * coordarr = MEDCoupling::DataArrayDouble::New();
  coordarr->alloc ( n_nodes, 3 );
  std::copy ( coord, coord + n_nodes * 3, coordarr->getPointer() );
  mesh->setCoords ( coordarr );
//     MEDLoader::WriteUMesh(namefile.c_str(), mesh, true);

  delete[] conn;
  delete[] coord;
  //  set up
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  double * sol = new double[n_nodes + 1];

  MEDCoupling::MEDCouplingFieldDouble * f  = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_CELLS );
  f->setMesh ( mesh );
  f->setName ( "weight" );
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array -> alloc ( n_elements, 1 );
  const int nel_e = _mgmesh._off_el[0][_NoLevels]; // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1]; // stop element
  int i = 0;

  // print quad -------------------------------------
  for ( int ivar = 0; ivar < 1; ivar++ )        {
      std::string var_name = "weight";

      for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {
          sol[i]  = _weight_ctrl[iel + nel_b];
          array ->setIJ ( i, 0, sol[i] );
          i++;
          }
      }

  f->setArray ( array );
  std::string s = "weight_control";
  s += ".med";
  MEDCoupling::WriteField ( s.c_str(), f, true );
  coordarr->decrRef();
  mesh->decrRef();
  f->decrRef();
  array->decrRef();
  delete[] sol;

  return;
  }
#endif
#ifdef HAVE_MED
// ============================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_u_med (
  std::string namefile, const int /* Level1*/
) { // ===============================================

  int n_nodes = _mgmesh._NoNodes[_NoLevels - 1];
  int n_elements = _mgmesh._NoElements[0][_NoLevels - 1];
  int nodes_el = _mgmesh._type_FEM[0];

  double * coord;
  coord = new double[n_nodes * 3];

  for ( int i = 0; i < n_nodes; i++ ) {
      coord[i * 3 + 1] = 0.;
      coord[i * 3 + 2] = 0.;

      for ( int idim = 0; idim < _mgmesh._dim; idim++ ) {
          coord[i * 3 + idim] = _mgmesh._xyz[i + idim * n_nodes];
//        std::cout << " "<< coord[i*3+idim];
          }

//      std::cout<< " "<< coord[i*3+2] << std::endl;
      }

  std::cout << " " << n_nodes << " " << n_elements << " " << nodes_el << std::endl;
  int  Level = _mgmesh._NoLevels - 1;
  int icount = 0;

  int * conn;
  conn = new int [n_elements * nodes_el];

  for ( int  iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
      for ( int el = _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level];
            el < _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level + 1]; el++ ) {
          for ( int  i = 0; i < nodes_el; i++ ) {
              conn[icount] = _mgmesh._el_map[0][el * nodes_el + i];
              icount++;
              }
          }
      }



  MEDCoupling::MEDCouplingUMesh * mesh = MEDCoupling::MEDCouplingUMesh::New ( "Mesh_1", _mgmesh._dim );
  mesh->allocateCells ( n_elements );

  for ( int  i = 0; i < n_elements; i++ ) {
      mesh->insertNextCell ( MED_EL_TYPE, nodes_el, conn + i * nodes_el );
      }

  mesh->finishInsertingCells();

  MEDCoupling::DataArrayDouble * coordarr = MEDCoupling::DataArrayDouble::New();
  coordarr->alloc ( n_nodes, 3 );
  std::copy ( coord, coord + n_nodes * 3, coordarr->getPointer() );
  mesh->setCoords ( coordarr );
  MEDCoupling::WriteUMesh ( namefile.c_str(), mesh, true );
//   coordarr->decRef();
//   mesh->decRef();

  delete[] coord;

  //  set up
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  double * sol = new double[n_nodes + 1];

  MEDCoupling::MEDCouplingFieldDouble * f  = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  f->setMesh ( mesh );
  f->setName ( _var_names[0].c_str() );
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array -> alloc ( n_nodes, 1 );

  // print quad -------------------------------------
  for ( int ivar = 0; ivar < 1; ivar++ )        {
      std::string var_name = _var_names[ivar];

      for ( int i = 0; i < n_nodes; i++ ) {
          sol[i]  = ( *x_old[Level] ) ( _node_dof[Level][i + ivar * offset] ) * _refvalue[ivar];
          array ->setIJ ( i, 0, sol[i] );

          }

      }

  f->setArray ( array );
  std::string s = "my_test3_sol_1";
  s += ".med";
  MEDCoupling::WriteField ( s.c_str(), f, true );

  delete [] conn;
  delete [] sol;
  return;
  }
#endif

// ================================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_bc ( std::string namefile, const int Level ) {

  //  setup
  const int n_nodes = _mgmesh._NoNodes[Level];
//   const int  n_nodes_l=(Level >0) ?_mgmesh._NoNodes[Level-1]:_mgmesh._NoNodes[2*_NoLevels];
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
//   int nt= _nvars[0]*n_nodes+_nvars[1]*n_nodes_l;
  int * sol = new int[n_nodes];

  // file hdf5
  hid_t file_id = H5Fopen ( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  hsize_t     dimsf[2];
  dimsf[0] = n_nodes;
  dimsf[1] = 1;

  // loop over all systems ---------------------------
  // print quad -------------------------------------
  for ( int ivar = 0; ivar < _nvars[2]; ivar++ ) {
      std::string var_name = _var_names[ivar] + "bd";

      for ( int i = 0; i < n_nodes; i++ ) {
          sol[i] = _bc[1][_node_dof[_NoLevels - 1][i + ( ivar ) * offset]];
          }

      _mgutils.print_Ihdf5 ( file_id, var_name, dimsf, sol );
      std::string var_name2 = _var_names[ivar] + "vl";

      for ( int i = 0; i < n_nodes; i++ ) {
          sol[i] = _bc[0][_node_dof[_NoLevels - 1][i + ivar * offset]];
          }

      _mgutils.print_Ihdf5 ( file_id, var_name2, dimsf, sol );
      }

  // print linear -----------------------------------
  double * sol_c = new double[_mgmesh._GeomEl.n_l[0]];
//   const int  lev_l=(Level >0) ?Level-1:_NoLevels;


  //put the attention to the end of the for cycle (_nvars[1]+_nvars[2])
  for ( int ivar = _nvars[2]; ivar < _nvars[1] + _nvars[2]; ivar++ ) {
      std::string var_name = _var_names[ivar] + "bd";

      // class writing routine
      // 2bA projection of the fine nodes -------------------------
      for ( int isubdom = 0; isubdom < _mgmesh._n_subdom; isubdom++ ) {
          for ( int i = _mgmesh._off_nd[0][isubdom * _NoLevels];
                i < _mgmesh._off_nd[0][isubdom * _NoLevels] +
                _mgmesh._off_nd[1][Level + 1 + isubdom * _NoLevels] -
                _mgmesh._off_nd[1][isubdom * _NoLevels]; i++ ) {

              sol[i] =  _bc[1][_node_dof[_NoLevels - 1][i + ivar * offset]];

              }
          } // 2bA end proj fine grid ----------------------------------

      //  2bB element interpolation over the fine mesh -----------------------
      for ( int iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
          for ( int iel = _mgmesh._off_el[0][_NoLevels * iproc + Level];
                iel < _mgmesh._off_el[0][Level + _NoLevels * iproc + 1]; iel++ ) {
              int indx = iel * NDOF_FEM;

              // element iel ----------------------------
              // vertices
              for ( int in = 0; in < NDOF_P; in++ ) {
                  sol_c[in] = sol[_mgmesh._el_map[0][indx + in]];
                  }

              for ( int in = 0; in < NDOF_FEM; in++ ) { // mid-points
                  double sum = 0;

                  for ( int jn = 0; jn < NDOF_P; jn++ ) {
                      sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn];
                      }

                  sol[_mgmesh._el_map[0][indx + in]] = sum;
                  }
              } // ---- end iel -------------------------
          } // 2bB end interpolation over the fine mesh ------------------------

      _mgutils.print_Ihdf5 ( file_id, var_name, dimsf, sol );



      std::string var_name2 = _var_names[ivar] + "vl";

      // class writing routine
      // 2bA projection of the fine nodes -------------------------
      for ( int isubdom = 0; isubdom < _mgmesh._n_subdom; isubdom++ ) {
          for ( int i = _mgmesh._off_nd[0][isubdom * _NoLevels];
                i < _mgmesh._off_nd[0][isubdom * _NoLevels] +
                _mgmesh._off_nd[1][Level + 1 + isubdom * _NoLevels] -
                _mgmesh._off_nd[1][isubdom * _NoLevels]; i++ ) {

              sol[i] =  _bc[0][_node_dof[_NoLevels - 1][i + ivar * offset]];

              }
          } // 2bA end proj fine grid ----------------------------------

      //  2bB element interpolation over the fine mesh -----------------------
      for ( int iproc = 0; iproc < _mgmesh._n_subdom; iproc++ ) {
          for ( int iel = _mgmesh._off_el[0][_NoLevels * iproc + Level];
                iel < _mgmesh._off_el[0][Level + _NoLevels * iproc + 1]; iel++ ) {
              int indx = iel * NDOF_FEM;

              // element iel ----------------------------
              // vertices
              for ( int in = 0; in < NDOF_P; in++ ) {
                  sol_c[in] = sol[_mgmesh._el_map[0][indx + in]];
                  }

              for ( int in = 0; in < NDOF_FEM; in++ ) { // mid-points
                  double sum = 0;

                  for ( int jn = 0; jn < NDOF_P; jn++ ) {
                      sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn];
                      }

                  sol[_mgmesh._el_map[0][indx + in]] = sum;
                  }
              } // ---- end iel -------------------------
          } // 2bB end interpolation over the fine mesh ------------------------

      _mgutils.print_Ihdf5 ( file_id, var_name2, dimsf, sol );

      } // ivar

  // clean -------------------------------
  delete []sol_c;
  delete []sol;
  H5Fclose ( file_id );

  return;
  }

// ===========================================================================================
/// This function reads the MSolDA system solution from namefile.h5
void MGSolDA::read_u (
  std::string namefile,  // filename (with path)
  int Level_restart      // restart Level-1 flag
) { //========================================================================================


/// A. Setup
  // mesh
  const int offset = _mgmesh._NoNodes[_NoLevels - 1]; // fine level # of nodes
  const int pie_offset = ( NSUBDOM ) * ( _mgmesh._NoElements[0][_NoLevels - 1] ) ; // fine level # of nodes
  // field
  double * sol;
  double * sol_pie;
  hid_t  file_sol;

  /// B. Read from  Level-1  (NoLevels-2)
  if ( Level_restart != 0 ) {

      const int ndigits = stoi ( _mgutils._sim_config["ndigits"] ); // digit for namefile
      const int restart = stoi ( _mgutils._sim_config["restart"] ); // restart label

      std::string resu_dir = _mgutils._inout_dir;
//     std::string i_aux_dir = _mgutils.get_file("MESH_AUX_DIR");    // mesh aux dir
      std::ostringstream namefile_sol;
      namefile_sol << resu_dir <<  "/sol.msh1." << setw ( ndigits ) << setfill ( '0' ) <<  restart << ".h5";    //  namefile_sol
      file_sol = H5Fopen ( namefile_sol.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

      std::ostringstream namefile_mesh;                               //   namefile_mesh
      namefile_mesh << "/home/filippo/software/new_femus_inst/femus/USER_APPL/turb_1/RESU_AUX/mesh.msh1.h5";
      hid_t  file_mesh = H5Fopen ( namefile_mesh.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

      int Level = _NoLevels - 2;                                    // Reading from Level

      const int n_elem = _mgmesh._NoElements[0][Level];             // coarse level # of elements
      const int offset_lev = _mgmesh._NoNodes[Level];               // coarse level # of nodes
      int * map_f2c = new int[offset];                              // node map fine to coarse

      std::cout << " Reading Level-1 sol from file " <<  namefile_sol.str() << "\n";
      std::cout << " Reading Level-1 connectivity from file " <<  namefile_mesh.str() << "\n";

      // file to read
      double sol_c[NDOF_FEM];
      double sol_f[NDOF_FEM];   // element connectivity
      sol = new double[offset_lev];             // coarse solution     (sol.xxx.h5)
      sol_pie = new double[1];

      std::ostringstream namedir_conn;                  // coarse connectivity (-> mesh.h5)
      namedir_conn << "/ELEMS/FEM0/MSH" << Level;
      int * map_coarse = new int[NDOF_FEM * n_elem];
      _mgutils.read_Ihdf5 ( file_mesh, namedir_conn.str().c_str(), map_coarse );

      int offel_coarse = 0; // # of element at Level for proc < _iproc

      for ( int jproc = 0; jproc < _iproc; jproc++ ) {
          offel_coarse += _mgmesh._off_el[0][jproc * _NoLevels + Level + 1] - _mgmesh._off_el[0][jproc * _NoLevels + Level];
          }


      for ( int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + Level + 1] -
            _mgmesh._off_el[0][_iproc * _NoLevels + Level]; iel++ ) {

          int    elem_gidx = ( iel + _mgmesh._off_el[0][_iproc * _NoLevels + Level] ) * NDOF_FEM;
          int    elem_coarse_gidx = ( iel + offel_coarse ) * NDOF_FEM;

          for ( int in = 0; in < NDOF_FEM; in++ ) {
              map_f2c[_mgmesh._el_map[0][elem_gidx + in]] = map_coarse[elem_coarse_gidx + in];
              }
          }


      // reading loop over system varables
      for ( int  ivar = 0; ivar < _n_vars; ivar++ ) {
          int  el_nds = NDOF_FEM;

          if ( ivar >= _nvars[2] ) {
              el_nds = NDOF_P;  // quad and linear
              }

          // reading ivar param
          _mgutils.read_Dhdf5 ( file_sol, "/" + _var_names[ivar], sol ); // reading coarse quad solution (sol.xxx.h5)
          double Irefval = 1. / _refvalue[ivar]; // units

          // storing  ivar variables (in parallell)
          for ( int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels] -
                _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]; iel++ ) {
              int    elem_gidx = ( iel + _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1] ) * NDOF_FEM;

              // vertex storage
              for ( int in = 0; in < NDOF_P; in++ ) { //  only vertices
                  int k_fine = _mgmesh._el_map[0][elem_gidx + in];
                  sol_c[in] = sol[ map_f2c[k_fine]];
                  }

              // interpolation
              for ( int in = 0; in < el_nds; in++ ) { // all element points
                  double sum = 0;

                  for ( int jn = 0; jn < NDOF_P; jn++ ) {
                      sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn];
                      }

                  sol_f[in] = sum;
                  }

              for ( int    i = 0; i < el_nds; i++ ) { // linear and quad
                  int k = _mgmesh._el_map[0][elem_gidx + i]; // the global node
                  x[_NoLevels - 1]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], sol_f[i]*Irefval ); //sol[k]*Irefval);    // set the field
                  }
              }
          }

      x[_NoLevels - 1]->localize ( *x_old[_NoLevels - 1] );

      // clean
      H5Fclose ( file_mesh );
      delete []map_coarse;
      }
  /// C. Read from  the same Level (NoLevels-1)
  else {

      int SolToRead = _NumRestartSol;
      std::string SolSuffix[3];
      SolSuffix[0] = "";
      SolSuffix[1] = "_old";
      SolSuffix[2] = "_oold";

      // file to read
      sol = new double[offset]; // temporary vector
      sol_pie = new double[pie_offset];
      file_sol = H5Fopen ( namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );

      for ( int oldTstep = 0; oldTstep < SolToRead ; oldTstep++ ) {
          // reading loop over system varables
          for ( int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++ ) {
              int  el_nds = NDOF_FEM;

              if ( ivar >= _nvars[2] ) {
                  el_nds = NDOF_P;  // quad and linear
                  }

              // reading ivar param
              _mgutils.read_Dhdf5 ( file_sol, "/" + _var_names[ivar] + SolSuffix[oldTstep], sol );
              double Irefval = 1. / _refvalue[ivar]; // units

              // storing  ivar variables (in parallell)
              for ( int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels] -
                    _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]; iel++ ) {
                  int  elem_gidx = ( iel + _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1] ) * NDOF_FEM;

                  for ( int  i = 0; i < el_nds; i++ ) { // linear and quad
                      int k = _mgmesh._el_map[0][elem_gidx + i]; // the global node
                      x[_NoLevels - 1]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], sol[k]*Irefval ); // set the field
                      }
                  }
              }

          int ndof_lev = 0;

          for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
              int delta = _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
              ndof_lev += delta;
              }

          // PIECEWISE VARIABLES
          for ( int ivar = _nvars[2] + _nvars[1]; ivar < _n_vars; ivar++ ) {
              int  el_nds = NDOF_FEM;
              el_nds = NSUBDOM; // quad and linear
              // reading ivar param
              _mgutils.read_Dhdf5 ( file_sol, "/" + _var_names[ivar], sol_pie );
              double Irefval = 1. / _refvalue[ivar]; // units

              // storing  ivar variables (in parallell)
              for ( int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels] -
                    _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]; iel++ ) {
                  int  k = ( iel + ndof_lev );
//                     int k=_mgmesh._el_map[0][elem_gidx*el_nds];   // the global node
                  x[_NoLevels - 1]->set ( _node_dof[_NoLevels - 1][k + ivar * offset], sol_pie[k * el_nds]*Irefval ); // set the field
                  }
              }

          if ( oldTstep == 0 ) x[_NoLevels - 1]->localize ( *x_old[_NoLevels - 1] );
          if ( oldTstep == 1 ) x[_NoLevels - 1]->localize ( *x_oold[_NoLevels - 1] );
          if ( oldTstep == 2 ) x[_NoLevels - 1]->localize ( *x_ooold[_NoLevels - 1] );
          }
      } //end if

  /// D. delocalization and clean

  H5Fclose ( file_sol );
  delete []sol;
  delete []sol_pie;
  return;
  }




/// ==================================================
///      ----------OPERATORS  READING-------------
///  -------------------------------------------------



/// This function read the Matrix linear-quad operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadMatrix ( const int Level,          //Level
                           const  std::string & namefile, // file name
                           SparseMatrixM & Mat,         // Matrix to read
                           const int nvars_in[]          // # quad variables
//                            const int nvars_l           // # linear variables
                         ) { // -p----------------------------------------------------

#ifdef PRINT_INFO // ------------------------------------------
  std::cout << " ReadMatrix(DA): start matrix reading    "  << std::endl;
#endif
  int NoLevels = _mgmesh._NoLevels; // max number of levels

  ///-------------------------------------------------------------------------------
  /// [1] Reading matrix from file hdf5---------------------------------------------
  ///-------------------------------------------------------------------------------

  // reading matrix dimensions
  //initialize dimension
  int dim_qlk[3][3], ldim[2];// [*][*] row variable type column variable type

  for ( int aa = 0; aa < 3; aa++ ) for ( int bb = 0; bb < 3; bb++ ) {
        dim_qlk[aa][bb] = 0;
        }

  for ( int iql = 0; iql < 3; iql++ ) {
      for ( int jql = 0; jql < 3; jql++ ) {
          std::ostringstream mode;
          mode << ( iql ) * 10 + ( jql ) * 1;
          int status = Mat.read_dim_hdf5 ( namefile.c_str(), mode.str().c_str(), ldim ); // reading dimensions

          if ( status == 0 ) {
              dim_qlk[iql][jql] = ldim[0];  //if "mode" found -> put the row number in dim_qlk[][]
              }
          }
      }

  // setup len and len_off vectors
  int ** length_row   = new int * [9]; // total length entry sparse struct
  int ** length_offrow = new int * [9]; // offdiag entry sparse matrix struct

  //length and offlength structure
  //    [0]kk [1]kl [2]kq
  //    [3]lk [4]ll [5]lq
  //    [6]qk [7]ql [8]qq

#ifdef HAVE_LASPACKM
  int ** pos_row      = new int * [9]; // pos matrix  struct
#endif

  // reading diagonal row length and off row length (0=k,1=l,2=q)
  for ( int iql = 0; iql < 3; iql++ ) { //cycle on the row blocks
      if ( nvars_in[iql] != 0 ) { // check if the matrix is needed for this row block
          for ( int jql = 0; jql < 3; jql++ ) { //cycle on the column blocks
              if ( nvars_in[jql] != 0 ) { // check if the matrix is needed for this column block
                  length_row[iql * 3 + jql] = new int[dim_qlk[iql][jql] + 1]; //prepare vectors of correct dimension
                  length_offrow[iql * 3 + jql] = new int[dim_qlk[iql][jql] + 1];
                  std::ostringstream mode;
                  mode << ( iql ) * 10 + ( jql ) * 1;
                  Mat.read_len_hdf5 ( namefile.c_str(), mode.str().c_str(), length_row[iql * 3 + jql], length_offrow[iql * 3 + jql] ); // reading sparse struct
#ifdef HAVE_LASPACKM
                  // reading compressed column index
                  pos_row[ql] = new int[length_row[ql][dim_ql[ql]]];      //  matrix pos
                  Mat.read_pos_hdf5 ( namefile.c_str(), ql, pos_row[ql] ); // reading sparse pos
#endif
                  }
              else {   //if there is no variable of jlq type assign -1 in the length vector
                  length_row[iql * 3 + jql] = new int[1];
                  length_row[iql * 3 + jql][0] = -1;
                  length_offrow[iql * 3 + jql] = new int[1];
                  length_offrow[iql * 3 + jql][0] = -1;
                  }
              }//end cycle on the column blocks
          }//if there is no variable of ilq type assign -1 in the length vector for all the column blocks
      else {
          for ( int kql = 0; kql < 3; kql++ ) {
              length_row[iql * 3 + kql] = new int[1];
              length_row[iql * 3 + kql][0] = -1;
              length_offrow[iql * 3 + kql] = new int[1];
              length_offrow[iql * 3 + kql][0] = -1;
              }
          }
      }//end cycle on the row blocks


  ///-------------------------------------------------------------------------------------
  /// [2] Set up initial matrix positions for levels, processors and type of variables----
  ///-------------------------------------------------------------------------------------
  // mesh index vector
  int * off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0]; // quad offset index(subdomains and levels)
  off_nd[1] = _mgmesh._off_nd[1]; // linear offset index (subdomains and levels)

  //prepare row indeces
  int ml[3];
  int m = 0;

  for ( int l = 0; l < 3; l++ ) {
      m += nvars_in[l] * dim_qlk[l][l];  // global last row
      }

  int n = m; // global last column
  ml[0] = ( _mgmesh._off_el[0][_iproc * NoLevels + Level + 1] -
            _mgmesh._off_el[0][_iproc * NoLevels + Level] ) * _el_dof[0]; // const local row
  ml[1] = off_nd[1][_iproc * NoLevels + Level + 1] - off_nd[1][_iproc * NoLevels]; // linear local row
  ml[2] = off_nd[0][_iproc * NoLevels + Level + 1] - off_nd[0][_iproc * NoLevels]; // quad local row
  int m_l = ml[0] * nvars_in[0] + ml[1] * nvars_in[1] + ml[2] * nvars_in[2]; // quad+lin+const last local row
  int n_l = m_l;   // quad+lin+const last local column

  // setting levels
  const  int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev_c = ( Level - 1 + _NoLevels + 1 ) % ( _NoLevels + 1 );
  int lev[3];
  lev[0] = Level;
  lev[1] = lev_c;
  lev[2] = Level;

  // initial local matrix positions
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0; // initialize ml_init

  for ( int isubd = 0; isubd < _iproc; isubd++ ) {
      ml_init[0] += ( _mgmesh._off_el[0][isubd * NoLevels + Level + 1] -
                      _mgmesh._off_el[0][isubd * NoLevels + Level] ) * _el_dof[0];   // const
      ml_init[1] += off_nd[1][Level + 1 + isubd * _NoLevels] - off_nd[1][isubd * _NoLevels]; // linear
      ml_init[2] += off_nd[0][Level + 1 + isubd * _NoLevels] - off_nd[0][isubd * _NoLevels]; // quad
      }

  int top_node_iproc = off_nd[0][_iproc * NoLevels];
  int ml_start = _node_dof[Level][top_node_iproc];

  ///-------------------------------------------------------------------------------------
  /// [3] Write matrix sparse structure---------------------------------------------------
  ///-------------------------------------------------------------------------------------

#ifdef HAVE_LASPACKM  // ------------------- only laspack ----------------

  // Graph dimension
  Graph graph;
  graph.resize ( m );
  graph._m = m;
  graph._n = n;
  graph._ml = m_l;
  graph._nl = n_l;

  graph._ml_start = ml_start;

  // quadratic
  for ( int  ivar = 0; ivar < nvar[0]; ivar++ ) {
      for ( int i = ml_init[0]; i < ml_init[0] + ml[2]; i++ ) {
          // row index
          len[0] = ( length_row[0][i + 1] - length_row[0][i] );
          int irowl = _node_dof[Level][_mgmesh._node_map[lev[2]][i] + ( ivar + 0 * nvars_in[2] ) * offset] - ml_start;
          len[1] = ( length_row[1][i + 1] - length_row[1][i] );

          graph[irowl].resize ( nvars_in[2]*len[0] + nvars_in[1]*len[1] + 1 );

//       int  jlq=0;
          for ( int  jlq = 0; jlq < 2; jlq++ ) { // ilq=0 quadratic-quadratic ilq=1 quadratic-linear
              for ( int  jvar = 0; jvar < nvars_in[jlq]; jvar++ ) {
                  for ( int j = 0; j < len[jlq]; j++ ) {
                      graph[irowl][j + jvar * len[jlq] + jlq * nvar[0]*len[0]] = _node_dof[Level][
                            _mgmesh._node_map[lev[jlq]][pos_row[jlq][j + length_row[jlq][i]]] + ( jvar + jlq * nvar[0] ) * offset];
                      }
                  }

              // last stored value is the number of in-matrix nonzero values
              graph[irowl][nvar[0]*len[0] + nvar[1]*len[1]] = nvar[0] * ( length_offrow[0][i + 1] - length_offrow[0][i] )
                  + nvar[1] * ( length_offrow[1][i + 1] - length_offrow[1][i] );

              }
          } // end quadratic ++++++++++++++++++++++++++++++++++++
      }

  // linear part ++++++++++++++++++++++++++++++++++++++++++++++++++
  for ( int  ivar = 0; ivar < nvar[1]; ivar++ ) {
      for ( int  i = ml_init[1]; i < ml_init[1] + ml[1]; i++ ) {
          // row index
          len[3] = ( length_row[3][i + 1] - length_row[3][i] );
          int irowl = _node_dof[Level][_mgmesh._node_map[lev_c][i] + ( ivar + nvar[0] ) * offset];
          len[2] = ( length_row[2][i + 1] - length_row[2][i] );

          graph[irowl].resize ( nvar[1]*len[3] + nvar[0]*len[2] + 1 );

          for ( int  ilq = 2; ilq < 4; ilq++ ) { // ilq=3 linear-linear ilq=2 linear-quadratic
              for ( int  jvar = 0; jvar < nvar[ilq - 2]; jvar++ ) {
                  for ( int j = 0; j < len[ilq]; j++ ) {
                      graph[irowl][j + jvar * len[ilq] + ( ilq - 2 ) *nvar[0]*len[2]] =
                        _node_dof[Level][_mgmesh._node_map[lev[ilq - 2]][pos_row[ilq][j + length_row[ilq][i]]] + ( jvar + ( ilq - 2 ) * nvar[0] ) * offset];
                      }
                  }
              }

          graph[irowl][nvar[1]*len[3] + nvar[0]*len[2]] = nvar[1] * ( length_offrow[3][i + 1] - length_offrow[3][i] )
              + nvar[0] * ( length_offrow[2][i + 1] - length_offrow[2][i] );

          }
      } // end linear

  // Update sparsity pattern -----------------------------
  Mat.update_sparsity_pattern ( graph );
  graph.clear();
#endif   // ------------------ ----------------------------------------
#ifdef HAVE_PETSCM  // -------------------  only petsc ----------------

  int cont_var = 0; //shifter for variables of all type
  std::vector< int> n_dz ( m_l ); // # diagonal entries
  std::vector< int> n_oz ( m_l ); // # offset entries

  for ( int  i = 0; i < m_l; i++ ) {
      n_dz[i] = 0;  //initialize to zero
      n_oz[i] = 0;
      }

  for ( int  ilq = 2; ilq >= 0; ilq-- ) { //cycle on type of variable (2=q,1=l,0=k) - row block
      for ( int ivar = 0; ivar < nvars_in[ilq]; ivar++ ) { //cycle on variables of ilq type
          for ( int  i1 = 0; i1 < ml[ilq]; i1++ ) { //cycle on the rows
              int i = ml_init[ilq] + i1;
              int i_top;

              if ( ilq != 0 ) {
                  i_top = _mgmesh._node_map[lev[ilq]][i];
                  }
              else {
                  i_top = i;  //for piecewise only
                  }

              int ind = _node_dof[Level][i_top + ( ivar + cont_var ) * offset] - ml_start; // ind is local (per proc) row index

              for ( int jlq = 2; jlq >= 0; jlq-- ) { //cycle on type of variable (2=q,1=l,0=k) - column block
                  int jlq1 = ( ilq + jlq ) % 3;

                  for ( int jvar = 0; jvar < nvars_in[jlq1]; jvar++ ) { //cycle on variables of jlq1 type
                      n_dz[ind] += ( length_row[ilq * 3 + jlq1][i + 1]   - length_row[ilq * 3 + jlq1][i] );
                      n_oz[ind] += ( length_offrow[ilq * 3 + jlq1][i + 1] - length_offrow[ilq * 3 + jlq1][i] );
                      }//jlq1 variables loop
                  }//column block loop
              } // row loop
          }//ilq variables loop

      cont_var += nvars_in[ilq]; // shift in the _node_dof for all the variables written till now (q, then l, then k)
      }

  for ( int  i = 0; i < m_l; i++ ) {
      n_dz[i] -= n_oz[i];   //count correctly diagonal entries
      }

  // Update sparsity pattern -----------------------------
  Mat.update_sparsity_pattern ( m, n, m_l, n_l, n_dz, n_oz );

#endif

  //  clean ----------------------------------------------------
  n_dz.clear();
  n_oz.clear();

  for ( int ql = 0; ql < 9; ql++ ) {
#ifdef HAVE_LASPACKM
      delete []pos_row[ql];  // sparse compressed postion
#endif
      delete []length_row[ql];
      delete []length_offrow[ql];
      }

#ifdef HAVE_LASPACKM
  delete []pos_row;  // sparse compressed postion
#endif
  delete []length_row;
  delete []length_offrow;

#ifdef PRINT_INFO // ------------------------------------------
  std::cout << " ReadMatrix(DA): end matrix reading    "  << std::endl;
#endif
  }


// =================================================================
/// This function read the Prolongation operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadProl (
  const int          Level,       // Level
  const std::string & name,       // file name (reading from)
  SparseMMatrixM  &  Prol,        // Prolongation matrix
  const int         nvars_in[],   // # quad, linear and const variables
  int               node_dof_f[],
  int               node_dof_c[]
) {

#ifdef PRINT_INFO // ------------------------------------------
  std::cout << " ReadProl(DA): start prolongator reading    "  << std::endl;
#endif
  ///-------------------------------------------------------------------------------
  /// [1] Reading prolongator from file hdf5----------------------------------------
  ///-------------------------------------------------------------------------------

  // reading dimensions
  int ldim[2];
  int dim_qlk[3][2]; // [0][0] piecewise fine, [2][1] quad coarse

  for ( int aa = 0; aa < 3; aa++ ) for ( int bb = 0; bb < 2; bb++ ) {
        dim_qlk[aa][bb] = 0;
        }

  for ( int jql = 0; jql < 3; jql++ ) {
      std::ostringstream mode;
      mode << jql;
      int status = Prol.read_dim_hdf5 ( name.c_str(), mode.str().c_str(), ldim );

      if ( status == 0 ) {
          dim_qlk[jql][0] = ldim[0];
          dim_qlk[jql][1] = ldim[1];
          }
      }

  // setup length (0=k,1=l,2=q)
  int ** length_row    = new int * [3];        // total length entry sparse struct
  int ** length_off_row = new int * [3];       // offdiag entry sparse matrix struct
  int ** pos_row       = new int * [3];        // compressed matrix  positions
  double ** val_row    = new double*[3];       // compressed matrix  values


  for ( int ql1 = 0; ql1 < 3; ql1++ ) //type of variable loop (0=k,1=l,2=q)
    if ( nvars_in[ql1] != 0 ) { //check if the prolongator is needed for this variable
        std::ostringstream mode;
        mode << ql1;
        // reading diagonal row length and off row length
        length_row[ql1] = new int[dim_qlk[ql1][0] + 1]; //prepare vectors of correct dimension
        length_off_row[ql1] = new int[dim_qlk[ql1][0] + 1];
        Prol.read_len_hdf5 ( name.c_str(), mode.str().c_str(), length_row[ql1], length_off_row[ql1] );
        // reading position and values of the prologator
        pos_row[ql1] = new int[length_row[ql1][dim_qlk[ql1][0]]]; //prepare vectors of correct dimension
        val_row[ql1] = new double[length_row[ql1][dim_qlk[ql1][0]]];
        Prol.read_pos_hdf5 ( name.c_str(), mode.str().c_str(), pos_row[ql1], val_row[ql1] );
        }
    else {   //if there is no variable of ql1 type assign -1 in the length vector
        length_row[ql1] = new int[1];
        length_row[ql1][0] = -1;
        length_off_row[ql1] = new int[1];
        length_off_row[ql1][0] = -1;
        pos_row[ql1] = new int[1];
        pos_row[ql1][0] = -1;
        val_row[ql1] = new double[1];
        val_row[ql1][0] = -1.;
        }

  ///------------------------------------------------------------------------------------------
  /// [2] Set up initial prolongator positions for levels, processors and type of variables----
  ///------------------------------------------------------------------------------------------
  // mesh index vector
  int * off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0]; // quad offset index   (subdomains and levels)
  off_nd[1] = _mgmesh._off_nd[1]; // linear offset index (subdomains and levels)

  // Set up Projector pattern -------------------------------------------
  int off_proc = _NoLevels * _iproc;
  int ml[3];   // row processor dofs
  ml[0] = ( _mgmesh._off_el[0][off_proc + Level + 1] -
            _mgmesh._off_el[0][off_proc + Level] ) * _el_dof[0]; //row constant processor dofs
  ml[1] = off_nd[1][off_proc + Level + 1] - off_nd[1][off_proc]; //row linear processor dofs
  ml[2] = off_nd[0][off_proc + Level + 1] - off_nd[0][off_proc]; //row quadratic processor dofs

  int nl[3]; // column processor dofs
  nl[0] = ( _mgmesh._off_el[0][off_proc + Level] -
            _mgmesh._off_el[0][off_proc + Level - 1] ) * _el_dof[0]; // column constant processor dofs
  nl[1] = off_nd[1][off_proc + Level] - off_nd[1][off_proc]; // column linear processor dofs
  nl[2] = off_nd[0][off_proc + Level] - off_nd[0][off_proc]; // column quadratic processor dofs

  // initial prolongator entry (parallel)
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0; //proc initial dofs [0]=const [1]=linear [2] quad

  for ( int isubd = 0; isubd < _iproc; isubd++ ) {
      ml_init[0] += ( _mgmesh._off_el[0][Level + 1 + isubd * _NoLevels] -
                      _mgmesh._off_el[0][Level + isubd * _NoLevels] ) * _el_dof[0];
      ml_init[1] += off_nd[1][Level + 1 + isubd * _NoLevels] - off_nd[1][isubd * _NoLevels];
      ml_init[2] += off_nd[0][Level + 1 + isubd * _NoLevels] - off_nd[0][isubd * _NoLevels];
      }

  // Local prolongator dimension (parallel)
  int m_l = ml[0] * nvars_in[0] + ml[1] * nvars_in[1] + ml[2] * nvars_in[2];
  int n_l = nvars_in[0] * nl[0] + nvars_in[1] * nl[1] + nvars_in[2] * nl[2];

  // global indices
  int m = nvars_in[0] * dim_qlk[0][0] + nvars_in[1] * dim_qlk[1][0] + nvars_in[2] * dim_qlk[2][0];
  int n = nvars_in[0] * dim_qlk[0][1] + nvars_in[1] * dim_qlk[1][1] + nvars_in[2] * dim_qlk[2][1];

  // pattern operator storage -----------------------------
  const  int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev[3];
  lev[0] = Level - 1;
  lev[1] = _NoLevels;
  lev[2] = Level - 1;

  if ( Level > 1 ) {
      lev[1] = Level - 2;
      }

  ///-------------------------------------------------------------------------------------
  /// [3] Write prolongator sparse structure----------------------------------------------
  ///-------------------------------------------------------------------------------------

#ifdef HAVE_LASPACKM // ------------ only laspack -----------------
  // pattern dimension
  Graph pattern;
  pattern.resize ( m );
  pattern._n = n;
  pattern._m = m; // global
  pattern._ml = m_l;
  pattern._nl = n_l; // local
  int  ml_start = _node_dof[Level][off_nd[0][off_proc]];
  pattern._ml_start = ml_start;

  // pattern structure
  for ( int  iql = 0; iql < 2; iql++ ) { // [0]=quad [1]=linear
      for ( int  ivar = 0; ivar < nvar[iql]; ivar++ ) {
          int  off_varl = ( ivar + iql * nvar[0] ) * offset;

          for ( int  i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++ ) { // element loop
              int irow = _node_dof[Level][_mgmesh._node_map[Level - iql][i] + off_varl];
              int  ncol = length_row[iql][i + 1] - length_row[iql][i];
              // pattern dimensions
              pattern[irow].resize ( ncol + 1 );
              pattern[irow][ncol] = length_off_row[iql][i + 1] - length_off_row[iql][i];

              // fill the pattern
              for ( int  j = 0; j < ncol; j++ ) {
                  pattern[irow][j] = _node_dof[Level - 1][
                                       _mgmesh._node_map[lev[iql] - 1 + iql][pos_row[iql][j + length_row[iql][i]]] + off_varl];
                  }
              }
          }// end var loop -------------------------------------
      }

  // Update sparsity pattern
  Prol.update_sparsity_pattern ( pattern );
  // clean
  pattern.clear();
#endif  // ------------ end only laspack -----------------
#ifdef HAVE_PETSCM
  std::vector< int> n_nz ( m_l ); // # diagonal entries
  std::vector< int> n_oz ( m_l ); // # offset entries

  for ( int  ilq = 0; ilq < 3; ilq++ ) { //type of variable loop (0=k,1=l,2=q)
      int mll_var = 0; //counter for preceding variables

      for ( int count = ilq + 1; count < 3; count++ ) {
          mll_var += nvars_in[count] * ml[count];
          }

      for ( int ivar = 0; ivar < nvars_in[ilq]; ivar++ ) { //loop on variables of ilq type
          for ( int  i = ml_init[ilq]; i < ml_init[ilq] + ml[ilq]; i++ ) { //row loop
              int len = ( length_row[ilq][i + 1] - length_row[ilq][i] );       // diag row length
              int len_off = ( length_off_row[ilq][i + 1] - length_off_row[ilq][i] ); // off row length
              n_oz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len_off; // set petsc off-row length
              n_nz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len - len_off ; // set petsc diag-row length
              } // row loop
          } // variables loop of ilq type
      }//type of variable loop

  Prol.update_sparsity_pattern ( m, n, m_l, n_l, n_oz, n_nz );
#endif

  ///-------------------------------------------------------------------------------------
  /// [4] Assemblying the Prolongator-----------------------------------------------------
  ///-------------------------------------------------------------------------------------

  DenseMatrixM * valmat;
  std::vector<int> tmp ( 1 );

  for ( int iql = 0; iql < 3; iql++ ) { //type of variable loop (0=k,1=l,2=q)
      int iql1 = 2 - iql; //to set the correct level (0 for quad, 1 for linear)
      int mll_var = 0; //counter for preceding variables

      for ( int  count = iql + 1; count < 3; count++ ) {
          mll_var += nvars_in[count];
          }

      for ( int ivar = 0; ivar < nvars_in[iql]; ivar++ ) { //loop on variables of iql type
          int off_varl = ( ivar + mll_var ) * offset;

          for ( int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++ ) { //row loop
              int ind_irow = i;

              if ( iql > 0 ) {
                  ind_irow = _mgmesh._node_map[Level - iql1][i];
                  }

              int irow = node_dof_f[ind_irow + off_varl];         // dof of the row
              int ncol = length_row[iql][i + 1] - length_row[iql][i]; // row length

              valmat = new DenseMatrixM ( 1, ncol );              // row storage
              tmp[0] = irow;
              std::vector<int> ind ( ncol );         // row-col indeces




              for ( int j = 0; j < ncol; j++ ) { // fill the row with values and the indices - column loop
                  int jcol = j + length_row[iql][i];
                  int jpos = pos_row[iql][jcol]; //set the position
                  int ind_jpos = jpos;

                  if ( iql > 0 ) {
                      ind_jpos = _mgmesh._node_map[lev[iql]][jpos];
                      }

                  ind[j] = node_dof_c[ind_jpos + off_varl];
//      if (iql==1) std::cout<<"  bc[0]  "<< bc[0][ind[j]]<<"  pos  "<< ind[j]<<"\n";
                  }//column loop




              for ( int j = 0; j < ncol; j++ ) { // fill the row with values and the indices - column loop
                  int jcol = j + length_row[iql][i];
                  double tmp = 1.;

                  if ( iql == 1 && _bc[0][irow] > 1.5 && _bc[0][irow] < 2.5 ) {
                      tmp = 1.;

                      }

                  ( *valmat ) ( 0, j ) = tmp * val_row[iql][jcol]; //assign the value in the matrix
                  int jpos = pos_row[iql][jcol]; //set the position
                  int ind_jpos = jpos;

                  if ( iql > 0 ) {
                      ind_jpos = _mgmesh._node_map[lev[iql]][jpos];
                      }

                  ind[j] = node_dof_c[ind_jpos + off_varl];
                  }//column loop

              Prol.add_matrix ( *valmat, tmp, ind ); //insert row in the Prolongator
              delete  valmat;                     //clean the row
              }//row loop
          }//loop on variables of iql type
      }//type of variable loop

  //  clean ----------------------------------------------
  for ( int ql = 0; ql < 3; ql++ ) {
      delete []pos_row[ql];
      delete []val_row[ql];
      delete []length_off_row[ql];
      delete []length_row[ql];
      }

  delete []pos_row;
  delete []val_row;
  delete []length_off_row;
  delete []length_row;

#ifdef PRINT_INFO
  std::cout << " ReadProl(DA): end reading Prol " << name.c_str() << std::endl;
#endif

  }

// =================================================================
/// This function read the Restriction linear-quad operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadRest (
  const int Level,      // Level
  const std::string & name, // file name (reading from)
  SparseMMatrixM & Rest,   // Restriction Matrix
  const int nvars_in[],      // # const,linear quad variables
  int  node_dof_f[],         // dof map fine mesh
  int  node_dof_c[],        // dof map coarse
  int  /*_node_dof_top*/[]      // dof map top level
) {  // -----------------------------------------------------------------

#ifdef PRINT_INFO // ------------------------------------------
  std::cout << " ReadRest(DA): start restrictor reading    "  << std::endl;
#endif
  ///-------------------------------------------------------------------------------
  /// [1] Reading restrictor from file hdf5-----------------------------------------
  ///-------------------------------------------------------------------------------

  // reading dimensions
  int ldim[2];
  int dim_qlk[3][2]; // [0][0] piecewise fine, [2][1] quad coarse

  for ( int aa = 0; aa < 3; aa++ ) for ( int bb = 0; bb < 2; bb++ ) {
        dim_qlk[aa][bb] = 0;
        }

  for ( int jql = 0; jql < 3; jql++ ) {
      std::ostringstream mode;
      mode << jql;
      int status = Rest.read_dim_hdf5 ( name.c_str(), mode.str().c_str(), ldim );

      if ( status == 0 ) {
          dim_qlk[jql][0] = ldim[0];
          dim_qlk[jql][1] = ldim[1];
          }
      else {
          std::cout << "  MGSolDA::ReadRest: error reading dimension ";
          abort();
          }
      }

  // setup length (0=k,1=l,2=q)
  int ** length_row    = new int * [3];        // total length entry sparse struct
  int ** length_off_row = new int * [3];       // offdiag entry sparse matrix struct
  int ** pos_row       = new int * [3];        // compressed matrix  positions
  double ** val_row    = new double*[3];       // compressed matrix  values

  // reading diagonal row length and off row length (0=qq.1=ql,3=lq,3=ll)
  for ( int ql1 = 0; ql1 < 3; ql1++ ) { //type of variable loop (0=k,1=l,2=q)
      if ( nvars_in[ql1] != 0 ) { //check if the restrictor is needed for this variable
          std::ostringstream mode;
          mode << ql1;
          // reading diagonal row length and off row length
          length_row[ql1] = new int[dim_qlk[ql1][0] + 1]; //prepare vectors of correct dimension
          length_off_row[ql1] = new int[dim_qlk[ql1][0] + 1];
          Rest.read_len_hdf5 ( name.c_str(), mode.str().c_str(), length_row[ql1], length_off_row[ql1] );
          // reading position and values of the restrictor
          pos_row[ql1] = new int[length_row[ql1][dim_qlk[ql1][0]]];  //prepare vectors of correct dimension
          val_row[ql1] = new double[length_row[ql1][dim_qlk[ql1][0]]];
          Rest.read_pos_hdf5 ( name.c_str(), mode.str().c_str(), pos_row[ql1], val_row[ql1] );
          }
      else {   //if there is no variable of ql1 type assign -1 in the length vector
          length_row[ql1] = new int[1];
          length_row[ql1][0] = -1;
          length_off_row[ql1] = new int[1];
          length_off_row[ql1][0] = -1;
          pos_row[ql1] = new int[1];
          pos_row[ql1][0] = -1;
          val_row[ql1] = new double[1];
          val_row[ql1][0] = -1.;
          }
      }

  // end reading file -------------------------------------------------------

  ///------------------------------------------------------------------------------------------
  /// [2] Set up initial restrictor positions for levels, processors and type of variables-----
  ///------------------------------------------------------------------------------------------
  // mesh index vector
  int * off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0]; // quad offset index(subdomainsand levels)
  off_nd[1] = _mgmesh._off_nd[1]; // linear offset index (subdomains and levels)

  // Set up Projector pattern -------------------------------------------
  // Matrix dimension
  // local indeces (parallel matrix)
  int off_proc = _NoLevels * _iproc;
  int ml[3];
  int nl[3]; //([2]=quad [1]=linear [0]=const)

  ml[0] = ( _mgmesh._off_el[0][off_proc + Level + 1] - _mgmesh._off_el[0][off_proc + Level] )
          * _el_dof[0];                   // local const m
  ml[1] = off_nd[1][off_proc + Level + 1] - off_nd[1][off_proc]; // local linear m
  ml[2] = off_nd[0][off_proc + Level + 1] - off_nd[0][off_proc]; // local quadratic m

  nl[0] = ( _mgmesh._off_el[0][off_proc + Level + 2] -
            _mgmesh._off_el[0][off_proc + Level + 1] ) * _el_dof[0]; // local const n
  nl[1] = off_nd[1][off_proc + Level + 2] - off_nd[1][off_proc]; // local linear n
  nl[2] = off_nd[0][off_proc + Level + 2] - off_nd[0][off_proc]; // local quadratic n

  // initial matrix entry (parallel)
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0; //proc initial dofs [0]=const [1]=linear [2] quad

  for ( int isubd = 0; isubd < _iproc; isubd++ ) {
      ml_init[0] += ( _mgmesh._off_el[0][Level + 1 + isubd * _NoLevels] -
                      _mgmesh._off_el[0][Level + isubd * _NoLevels] ) * _el_dof[0];
      ml_init[1] += off_nd[1][Level + 1 + isubd * _NoLevels] - off_nd[1][isubd * _NoLevels];
      ml_init[2] += off_nd[0][Level + 1 + isubd * _NoLevels] - off_nd[0][isubd * _NoLevels];
      }

  // Local restrictor dimension (parallel)
  int m_l = ml[0] * nvars_in[0] + ml[1] * nvars_in[1] + ml[2] * nvars_in[2];     //  local m
  int n_l = nvars_in[0] * nl[0] + nvars_in[1] * nl[1] + nvars_in[2] * nl[2];     //  local n

  // global indices
  int m = nvars_in[0] * dim_qlk[0][0] + nvars_in[1] * dim_qlk[1][0] + nvars_in[2] * dim_qlk[2][0];
  int n = nvars_in[0] * dim_qlk[0][1] + nvars_in[1] * dim_qlk[1][1] + nvars_in[2] * dim_qlk[2][1];

  // pattern operator storage -----------------------------
  const  int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev_c = ( Level - 1 + _NoLevels + 1 ) % ( _NoLevels + 1 ); // coarse linear level
  int lev[3];
  lev[0] = Level;
  lev[1] = lev_c;
  lev[2] = Level;

  ///-------------------------------------------------------------------------------------
  /// [3] Write restrictor sparse structure-----------------------------------------------
  ///-------------------------------------------------------------------------------------

#ifdef HAVE_PETSCM  // ----------- only petsc -----------------------------
  std::vector< int> n_nz ( m_l ); // # diagonal entries
  std::vector< int> n_oz ( m_l ); // # offset entries

//   if (nvars_in[0]==0) {
  // linear + quadratic  operator
  for ( int  ilq = 0; ilq < 3; ilq++ ) { //type of variable loop (0=k,1=l,2=q)
      int mll_var = 0; //counter for preceding variables

      for ( int  count = ilq + 1; count < 3; count++ ) {
          mll_var += nvars_in[count] * ml[count];
          }

      for ( int ivar = 0; ivar < nvars_in[ilq]; ivar++ ) { //loop on variables of ilq type
          for ( int  i = ml_init[ilq]; i < ml_init[ilq] + ml[ilq]; i++ ) { //row loop
              int len = ( length_row[ilq][i + 1] - length_row[ilq][i] );       // diag row length
              int len_off = ( length_off_row[ilq][i + 1] - length_off_row[ilq][i] ); // off row length
              n_oz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len_off ;
              n_nz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len - len_off ;
              } //row loop
          } // variables loop of ilq type
      }//type of variable loop

  Rest.update_sparsity_pattern ( m, n, m_l, n_l, n_oz, n_nz );
#endif // ----------- end only petsc --------------------------------------

#ifdef HAVE_LASPACKM    // ----------- only laspack  -----------------------------
  // pattern structure dimension
  Graph pattern;
  pattern.resize ( m );   // values.resize(nrowt);
  pattern._m = m;
  pattern._n = n;      // global dim _m x _n
  pattern._nl = n_l; //  local _n
  pattern._ml = m_l; //  local _m
  int  ml_start = _node_dof[Level][off_nd[0][off_proc]]; //  offset proc nodes
  pattern._ml_start = ml_start; // starting indices for local matrix

  // quadratic matrix dof (mlq_init-mlq_init+ml_q)
  for ( int  iql = 0; iql < 2; iql++ ) {      // [0]=quad [1]=linear
      for ( int  ivar = 0; ivar < nvar[iql]; ivar++ ) { // variable loop
          int  off_var_q = ( ivar + iql * nvar[0] ) * offset;

          for ( int i = ml_init[ iql]; i < ml_init[ iql] + ml[ iql]; i++ ) { // element loop
              int  ncol = length_row[ iql][i + 1] - length_row[ iql][i];          // # columns
              int irow = _node_dof[Level][_mgmesh._node_map[lev[iql]][i] + off_var_q]; // row
              // pattern data
              pattern[irow].resize ( ncol + 1 );                                     // dimension
              pattern[irow][ncol] = length_off_row[ iql][i + 1] - length_off_row[iql][i]; // # entry off diag

              for ( int  j = 0; j < ncol; j++ ) { // inserting  row data
                  pattern[irow][j] = _node_dof[Level + 1][
                                       _mgmesh._node_map[Level + 1][pos_row[ iql][j + length_row[ iql][i]]] + off_var_q];
                  } //  row data
              } // element loop
          } // variable loop
      } // [0]=quad [1]=linear

  //  update sparsity pattern for Rst ------------------
  Rst[Level]->update_sparsity_pattern ( pattern );
  pattern.clear();
#endif  // ----------- end only laspack  -----------------------------

  ///-------------------------------------------------------------------------------------
  /// [4] Assemblying the Restrictor------------------------------------------------------
  ///-------------------------------------------------------------------------------------

  DenseMatrixM * valmat;
  std::vector<int> tmp ( 1 );

  for ( int iql = 0; iql < 3; iql++ ) { //type of variable loop (0=k,1=l,2=q)
      int iql1 = 2 - iql; //to set the correct level (0 for quad, 1 for linear)
      int mll_var = 0; //counter for preceding variables

      for ( int count = iql + 1; count < 3; count++ ) {
          mll_var += nvars_in[count];
          }

      for ( int ivar = 0; ivar < nvars_in[iql]; ivar++ ) { //loop on variables of iql type
          int off_val_l = ( ivar + mll_var ) * offset;

          for ( int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++ ) { //row loop
              int top_node = i;

              if ( iql > 0 ) {
                  top_node = _mgmesh._node_map[lev[iql]][i];
                  }

              int irow = node_dof_c[top_node + off_val_l];      // dof of the row
              int ncol = length_row[iql][i + 1] - length_row[iql][i]; // row length - # colunms
              int irow_top = _node_dof[_NoLevels - 1][top_node + off_val_l];
              valmat = new DenseMatrixM ( 1, ncol );      // row storage
              tmp[0] = irow;
              std::vector< int> ind ( ncol ); // row-col indeces

              for ( int j = 0; j < ncol; j++ ) { // fill the row with values and the indices - column loop
                  int jcol = j + length_row[iql][i];
                  int jpos = pos_row[iql][jcol]; //set the position
                  int ind_jpos = jpos;

                  if ( iql > 0 ) {
                      ind_jpos = _mgmesh._node_map[Level + 1 - iql1][jpos];
                      }

                  ind[j] = node_dof_f[ind_jpos + off_val_l];
                  ( *valmat ) ( 0, j ) =
                    ( _bc[0][irow_top] % 2 ) *
                    val_row[iql][j + length_row[iql][i]]; //assign the value in the matrix
                  }//column loop

              Rest.add_matrix ( *valmat, tmp, ind ); //insert row in the Restrictor
              delete  valmat;                    //clean the row
              }//row loop
          }//loop on variables of iql type
      }//type of variable loop

  //  clean   -----------------------------------
  for ( int ql = 0; ql < 3; ql++ ) {
      delete []pos_row[ql];
      delete []val_row[ql];
      delete []length_off_row[ql];
      delete []length_row[ql];
      }

  delete []pos_row;
  delete []val_row;
  delete []length_off_row;
  delete []length_row;

#ifdef PRINT_INFO
  std::cout << " ReadRest(DA): end reading Rest " << name.c_str() << std::endl;
#endif
  return;
  }









/// ==================================================
///   ==============SOLVERS  SOURCE================
/// ==================================================



/// ======================================================
/// This function controls the time step operations:
/// ======================================================
void MGSolDA::MGTimeStep ( const double time, const int ) {

  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if  PRINT_TIME==1
  std::clock_t start_time = std::clock();
#endif
  GenMatRhs ( time, _NoLevels - 1, 1 );

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
      GenMatRhs ( time, Level, 0 );
      }

#if    PRINT_TIME==1
  std::clock_t end_time = std::clock();
  std::cout << " Assembly time =" << double ( end_time - start_time ) / CLOCKS_PER_SEC
            << " s " << std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve ( 1.e-6, 40 );
#if    PRINT_TIME==1
  std::clock_t end_timef = std::clock();
  std::cout << " Assembly+solution time =" << double ( end_timef - start_time ) / CLOCKS_PER_SEC
            << "s " << std::endl;
#endif
  /// [d] Update of the old solution at the top Level
//   x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
  return;
  }



///  ====================================================
/// This function assembles the matrix and the rhs:
///  ====================================================
void  MGSolDA::GenMatRhs (
  const double /*time*/,   // time  <-
  const int    Level,  // Level <-
  const int    mode    // mode  <- (1=rhs+matrix) (0=only matrix)
) {  // ===============================================

  /// Set up
  // geometry and bc---------------------------------------------------------------------------------
  const int ndim = DIMENSION;                         // dimension
  const int offset = _mgmesh._NoNodes[_NoLevels - 1]; // mesh nodes
  const int el_sides = _mgmesh._GeomEl._n_sides[0];   // element sides
  int       el_conn[NDOF_FEM];// elb_conn[NDOF_FEMB];   // element connectivity
  int       el_neigh[NDOF_FEM];                       // bd element connectivity
  int       elb_conn[NDOF_FEMB];// elb_conn[NDOF_FEMB];   // element connectivity
  int       sur_toply[NDOF_FEMB];                     // boundary topology
  double    xx_qnds[DIMENSION * NDOF_FEM];            // element node coords
  double    xxb_qnds[DIMENSION * NDOF_FEMB];          // boundary of the element node coords
  int       _bc_vol[NDOF_FEM * 10];              // element  b.cond flags (Neu or Dir)
  int       _bc_bd[NDOF_FEM * 10];               // element  b.cond flags (different possibilities)

  // gauss integration  -----------------------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim - 1];       // quadratic element gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION - 2]; // quadratic bd elem gauss points
  double det[3], JxW_g[3], InvJac[3][DIMENSION * DIMENSION]; // determinant of Jacobean and det*gauss_weigth
  double dphijdx_g[3][DIMENSION];                          // global derivatives at gauss point
  double dphiidx_g[3][DIMENSION];                          // global derivatives at gauss point
//   double _dt =0.1;

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int elb_ndof[3];
  elb_ndof[0] = NDOF_BK;
  elb_ndof[1] = NDOF_PB;
  elb_ndof[2] = NDOF_FEMB; // number of element boundary dofs
  int el_mat_nrows = 0;                                             // number of matrix rows (dofs)

  for ( int ideg = 0; ideg < 3; ideg++ ) {
      el_mat_nrows += _nvars[ideg] * _el_dof[ideg];
      }

  int el_mat_ncols = el_mat_nrows;                    // square matrix
  std::vector<int> el_dof_indices ( el_mat_ncols );   // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  int DA_f = _data_eq[2].indx_ub[_data_eq[2].tab_eqs[DA_F]]; // Quadratic
  double _ub_g[3][14];                                     // values of external fields

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix+rhs) ---------------------------
  A[Level]->zero();

  if ( mode == 1 ) {
      b[Level]->zero();    // global matrix A and rhs b
      }

  DenseMatrixM KeM;
  DenseVectorM FeM;                              // local  matrix KeM and rhs FeM
  KeM.resize ( el_mat_nrows, el_mat_ncols );
  FeM.resize ( el_mat_nrows ); // resize local matrix and rhs

  // number of total elements for level
  int ndof_lev = 0;

  for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
      int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
      ndof_lev += delta;
      }

  // test pass
// #ifdef HAVE_COUPLING
//     MGSystemExtended * ext_es= static_cast<MGSystemExtended *>(&_mgphys);
//      double test_pass= ext_es->_nComp;
// #endif


  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

  for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {

      KeM.zero();
      FeM.zero();      // set to zero matrix and rhs

      // geometry and element  fields ------------------------------------

      _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, xx_qnds ); // Element Connectivity (el_conn) and coordinates (xx_qnds)
      _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh ); // Neighbors of the element

      // set element-nodes variables  bc (bc_q_dofs)
      get_el_dof_bc ( Level, iel + ndof_lev, _el_dof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

      // grid and gaussian points
      for ( int idim = 0; idim < DIMENSION; idim++ ) {
          for ( int d = 0; d < NDOF_FEM; d++ ) {
              _data_eq[2].ub[idim * NDOF_FEM + d] = xx_qnds[idim * NDOF_FEM + d]; // element nodes xxg (DIM)
              }

          // element grid distance
          }

      // element field values
//         for (int deg=0; deg<3; deg++) {
      for ( int eq = 0; eq < _data_eq[2].n_eqs; eq++ ) {
          _data_eq[2].mg_eqs[eq]->get_el_sol ( 0, _data_eq[2].indx_ub[eq + 1] - _data_eq[2].indx_ub[eq],
                                               _el_dof[2], el_conn, offset, _data_eq[2].indx_ub[eq], _data_eq[2].ub );
          }

//         }
      //linear field
      _data_eq[1].mg_eqs[0]->get_el_sol ( _nvars[2], _data_eq[1].indx_ub[0], _el_dof[1], el_conn, offset, 0, _data_eq[1].ub );
      //  external node quantities -------------------------------------

      _data_eq[0].mg_eqs[0]->get_el_sol_piece ( _nvars[2] + _nvars[1], _data_eq[0].indx_ub[0], _el_dof[0], iel + ndof_lev, offset, 0, _data_eq[0].ub );


      // ======================================================================
      // Volume =============================================================
      // ======================================================================
      // ---------------------------------------------
      /// c) gaussian integration loop (n_gauss)
      // --------------------------------------------
      for ( int qp = 0; qp < el_ngauss; qp++ ) {

          // shape functions at gaussian points -----------------------------------
          for ( int ideg = 1; ideg < 3; ideg++ ) { // linear-quadratic  [1]=linear [2]=quad
              det[ideg]      = _fe[ideg]->Jac ( qp, xx_qnds, InvJac[ideg] ); // Jacobian
              JxW_g[ideg] = det[ideg] * _fe[ideg]->_weight1[ndim - 1][qp];  // weight
              _fe[ideg]->get_phi_gl_g ( ndim, qp, _phi_g[ideg] );          // shape funct
              _fe[ideg]->get_dphi_gl_g ( ndim, qp, InvJac[ideg], _dphi_g[ideg] ); // global coord deriv
              }

          JxW_g[0] = JxW_g[2];
          _fe[0]->get_phi_gl_g ( ndim, qp, _phi_g[0] );          // shape function piecewise


          //  fields -----------------------------------------------------------
          interp_el_sol ( _data_eq[2].ub, 0, _data_eq[2].indx_ub[_data_eq[2].n_eqs],
                          _phi_g[2], _el_dof[2], _ub_g[2] ); // quadratic
          interp_el_sol ( _data_eq[1].ub, 0, _data_eq[1].indx_ub[0], _phi_g[1], _el_dof[1], _ub_g[1] ); // linear
          interp_el_sol ( _data_eq[0].ub, 0, _data_eq[0].indx_ub[0], _phi_g[0], _el_dof[0], _ub_g[0] ); //constant


          /// d) Local (element) assemblying DA equation
          // *********************** *******************************

          for ( int i = 0; i < _el_dof[2]; i++ )     { //  --- QUADRATIC ---
              // set up row i
              const double phii_g = _phi_g[2][i];

              for ( int idim = 0; idim < ndim; idim++ ) {
                  dphiidx_g[2][idim] = _dphi_g[2][i + idim * _el_dof[2]];
                  }

//         int gl_node=el_conn[i];

              for ( int ivar = 0; ivar < _nvars[2]; ivar++ ) {

                  double dtxJxW_g = JxW_g[2] * _bc_vol[i + ivar * _el_dof[2]];
                  int index = i + ivar * _el_dof[2];

                  // Rhs Assemblying  ----------------------------------------------------------
                  if ( mode == 1 ) {

                      // rhs

                      FeM ( index ) += dtxJxW_g * (
                                         _ub_g[2][DA_f + ivar] * phii_g / _dt // time
                                         + 2.*phii_g                     // heat source
                                       );
                      }

                  // Matrix Assemblying ---------------------------
                  for ( int j = 0; j < _el_dof[2]; j++ ) {
                      double phij_g = _phi_g[2][j];
                      double Lap = 0.;

                      for ( int idim = 0; idim < ndim; idim++ ) {
                          dphijdx_g[2][idim] = _dphi_g[2][j + idim * _el_dof[2]];
                          Lap += dphijdx_g[2][idim] * dphiidx_g[2][idim];         // Laplacian
                          }

                      // energy-equation
                      KeM ( index, j + ivar * _el_dof[2] ) += dtxJxW_g * (
                          phii_g * phij_g / _dt // time term
                          + Lap                   //diff
                                                              );
                      }
                  } // quadratic variable cycle
              } // ----------------------------------------END QUADRATIC


          for ( int i = 0; i < _el_dof[1]; i++ )     { //    --- LINEAR ---
              // set up row i
              const double phii_g = _phi_g[1][i];

              for ( int idim = 0; idim < ndim; idim++ ) {
                  dphiidx_g[1][idim] = _dphi_g[1][i + idim * _el_dof[1]];
                  }

              for ( int ivar = 0; ivar < _nvars[1]; ivar++ ) {

                  double dtxJxW_g = JxW_g[1] * _bc_vol[i + ( ivar + _nvars[2] ) * NDOF_FEM];
                  int index = i + ivar * _el_dof[1] + _el_dof[2] * _nvars[2];

                  // Rhs Assemblying  ----------------------------------------------------------
                  if ( mode == 1 ) {

                      // rhs

                      FeM ( index ) += dtxJxW_g * (
                                         _ub_g[1][ivar] * phii_g / _dt // time
                                         + 1.*phii_g                   // heat source
                                       );
                      }


                  // Matrix Assemblying ---------------------------
                  for ( int j = 0; j < _el_dof[1]; j++ ) {
                      double phij_g = _phi_g[1][j];
                      int jndex = j + _el_dof[2] * _nvars[2] + ivar * _el_dof[1];
                      double Lap = 0.;

                      for ( int idim = 0; idim < ndim; idim++ ) {
                          dphijdx_g[1][idim] = _dphi_g[1][j + idim * _el_dof[1]];
                          Lap += dphijdx_g[1][idim] * dphiidx_g[1][idim];         // Laplacian
                          }

                      // energy-equation
                      KeM ( index, jndex ) += dtxJxW_g * (
                                                phii_g * phij_g / _dt // time term
                                                + Lap                   //diff
                                              );
                      }
                  } //// linear variable cycle
              } // ----------------------------------END LINEAR


          for ( int i = 0; i < _el_dof[0]; i++ )     { //    --- Piecewise ---
              // set up row i
              const double phii_g = _phi_g[0][i];

              for ( int ivar = 0; ivar < _nvars[0]; ivar++ ) {

                  double dtxJxW_g = JxW_g[0] * _bc_vol[i + ( ivar + _nvars[1] + _nvars[2] ) * NDOF_FEM];
                  int index = i + ivar * _el_dof[0] + _el_dof[2] * _nvars[2] + _el_dof[1] * _nvars[1];

                  // Rhs Assemblying  ----------------------------------------------------------
                  if ( mode == 1 ) {
                      // rhs
                      FeM ( index ) += dtxJxW_g * (
//                        _ub_g[2][0]*phii_g
                                         _ub_g[0][ivar] * phii_g / _dt // time
                                         + _data_eq[2].ub[0] * phii_g * ( 1 - ivar ) //source x
                                         + _data_eq[2].ub[0] * _data_eq[2].ub[1] * _data_eq[2].ub[2] * phii_g * ivar //source x*y*z
                                       );
                      }

                  // Matrix Assemblying ---------------------------
                  for ( int j = 0; j < _el_dof[0]; j++ ) {
                      double phij_g = _phi_g[0][j];
                      int jndex = j + ivar * _el_dof[0] + _el_dof[2] * _nvars[2] + _el_dof[1] * _nvars[1];

                      KeM ( index, jndex ) += dtxJxW_g * ( phii_g * phij_g / _dt ); // time term
                      }
                  } ////  variable cycle
              } // ----------------------------------END piecewise





          } // end of the quadrature point qp-loop ***********************


      // ======================================================================
      // ====================== boundary ======================================
      // ======================================================================

      for ( int iside = 0; iside < el_sides; iside++ )  {
          if ( el_neigh[iside] == -1 ) {

              for ( int idof = 0; idof < NDOF_FEMB; idof++ ) {                   //idof -> boundary node
                  sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside]; //use map to find global node
                  int idofb = sur_toply[idof];                                     //idofb -> element node
                  elb_conn[idof] = el_conn[idofb];                                 //connectivity vector

                  for ( int idim = 0; idim < DIMENSION; idim++ ) {
                      xxb_qnds[idim * NDOF_FEMB + idof] = xx_qnds[idim * NDOF_FEM + idofb]; //get boundary coordinates
                      _data_eq[2].ub[idim * NDOF_FEM + idof] = xxb_qnds[idim * NDOF_FEMB + idof];
                      }
                  }

              for ( int iql = 2; iql < 3; iql++ )    // --- QUADRATIC ---
                for ( int ivar = 0; ivar < _nvars[iql]; ivar++ )    {
                    // Dirichlet boundary conditions  ***********************************
                    if ( _bc_vol[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM] == 0 ) {

                        //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))
                        int bc_s = ( int ) _bc_bd[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM]; // b cond
                        double det = _fe[iql]->JacSur ( elb_ngauss - 1, xxb_qnds, InvJac[2] ); // jacobian
                        double Ipenalty = det / _dt;                           // Dirichlet bc flag

                        // local boundary loop   ---------------------------------------
                        for ( int lb_node = 0; lb_node < elb_ndof[iql]; lb_node++ ) {
                            int index = sur_toply[lb_node] + ivar * NDOF_FEM; // local vol index
                            int gl_node_bd = elb_conn[lb_node];
                            // flag setup (\int bc_var*T+bc_var*val)
                            int  bc_val = ( int ) ( ( bc_s & 2 ) >> 1 ); // (1?)non homogeneous
                            int  bc_var = ( int ) ( bc_s % 2 ); // (?1) variable

                            // Assemblying  Matrix & rhs
                            if ( mode == 1 ) {
                                FeM ( index ) += bc_val * Ipenalty * _data_eq[iql].ub[DA_f * NDOF_FEM + index];
                                FeM ( index ) += ( 1 - bc_val ) * bc_var * Ipenalty * ( 0. );
                                FeM ( index ) += Ipenalty * ( 1. );
                                }

                            KeM ( index, index ) += Ipenalty; //  Dirichlet bc
                            }// lb_node -end  local boundary loop -------------------------
                        } // end if Dirichlet  boundary conditions

                    // **********************************************************************

                    else if ( _bc_vol[sur_toply[NDOF_FEMB - 1] + ivar * NDOF_FEM] != 0 ) {

                        // Non homogenous Neumann boundary conditions  ***********************************

                        // gaussian integration loop (n_gauss)
                        // -----------------------------------------------
                        for ( int qp = 0; qp <  elb_ngauss; qp++ ) {

                            // quad/linear  [2]=quad [1]=linear------------------------------------
                            det[iql]  = _fe[iql]->JacSur ( qp, xxb_qnds, InvJac[2] ); // local coord _phi_g and jac
                            JxW_g[iql] = det[iql] * _fe[iql]->_weight1[ndim - 2][qp]; // weight
                            _fe[iql]->get_phi_gl_g ( ndim - 1, qp, _phi_g[iql] ); // global coord _phi_g


#ifdef AXISYM   // axisymmetric  (index ->0)
                            interp_el_sol ( _data_eq[iql].ub, 0, DIMENSION, _phi_g[iql], elb_ndof[iql], _ub_g[iql] );
                            JxW_g[iql]  *= _ub_g[iql][0];

#endif

                            // local side loop (over the node face)
                            for ( int lsi_node = 0; lsi_node < elb_ndof[iql]; lsi_node++ ) {

                                // set up row i
                                const double phii_g = _phi_g[iql][lsi_node]; // boundary test function
                                int index = sur_toply[lsi_node] + ivar * NDOF_FEM; // local element index
                                int bc_s = ( int ) _bc_bd[index];
                                int bc_v = ( int ) _bc_vol[index];
                                double dtxJxW_g = JxW_g[iql] * bc_v; // Neumann bc flag and Jac

                                // flag setup: +++++++++++++++++++++++++++++++++++++++
                                //  \int (bc_var*T+bc_val*val)ds
                                int  bc_val = ( int ) ( ( bc_s & 2 ) >> 1 ); // (1?) non-homogeneous
                                int  bc_var = ( int ) ( bc_s % 2 ); // (?1) tau=A*variable
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                                // Assemblying rhs ----------------------------
                                if ( mode == 1 ) {
                                    FeM ( index ) += bc_val * ( 1 - bc_var ) * dtxJxW_g * phii_g * 1.;
                                    }

                                // Assemblying Matrix ---------------------------------
                                for ( int lsj_node = 0; lsj_node < elb_ndof[iql];  lsj_node++ ) {
                                    int jndex = sur_toply[lsj_node] + ivar * NDOF_FEM;
                                    KeM ( index, jndex ) += dtxJxW_g * bc_var * phii_g * _phi_g[iql][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
                                    }// end j  ---------------------------------------

                                }// i   +++++++++++++++++++++++++++++++
                            } // end of the quadrature point qp-loop **********************


                        } // Neumann non homog

                    }    // ---END QUADRATIC ---


//        } //end if side


              for ( int iql = 1; iql < 2; iql++ )  { // --- LINEAR ---
                  for ( int ivar = 0; ivar < _nvars[iql]; ivar++ )    {
                      // Dirichlet boundary conditions  ***********************************
//         int dn_flag=0;
//            for (int lb_node=0; lb_node< elb_ndof[iql]; lb_node++)
//       if (_bc_vol[sur_toply[lb_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM] ==1) dn_flag++;


//       if (dn_flag< 2 ){ // Dirichlet  boundary conditions <-- ONLY THIS
                      //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))

                      double det = _fe[iql]->JacSur ( 0, xxb_qnds, InvJac[2] ); // jacobian
                      double Ipenalty = det / _dt;                           // Dirichlet bc flag

                      // local boundary loop   ---------------------------------------
                      for ( int lb_node = 0; lb_node < elb_ndof[iql]; lb_node++ ) {
                          if ( _bc_vol[sur_toply[lb_node] + ( ivar + _nvars[2] ) *NDOF_FEM] == 0 ) {

                              int index_bd = sur_toply[lb_node] + ( ivar + _nvars[2] ) * NDOF_FEM; // local vol index
                              int index_fem = sur_toply[lb_node] + ivar * _el_dof[iql] + _nvars[2] * NDOF_FEM;
                              // flag setup (\int bc_var*T+bc_var*val)
                              int bc_s = ( int ) _bc_bd[index_bd]; // b cond
//                             int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                              int  bc_var = ( int ) ( bc_s % 2 ); // (?1) variable

                              // Assemblying  Matrix & rhs
                              if ( mode == 1 ) {
                                  FeM ( index_fem ) += ( 1 - bc_var ) * Ipenalty * _data_eq[iql].ub[ivar * NDOF_FEM + sur_toply[lb_node]];
                                  FeM ( index_fem ) += bc_var * Ipenalty * ( 5. );
                                  }

                              KeM ( index_fem, index_fem ) += Ipenalty; //  Dirichlet bc
                              }//end if on bc_vol
                          }// lb_node -end  local boundary loop -------------------------

//          } // end if Dirichlet  boundary conditions
                      // **********************************************************************

//         else if (_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {
//           if (dn_flag>0){
                      // Non homogenous Neumann boundary conditions  ***********************************

                      // gaussian integration loop (n_gauss)
                      // -----------------------------------------------
//           for (int qp=0; qp<  elb_ngauss; qp++) {
//
//             // quad/linear  [2]=quad [1]=linear------------------------------------
//             det[iql]  = _fe[iql]->JacSur(qp,xxb_qnds);    // local coord _phi_g and jac
//             JxW_g[iql]=det[iql]*_fe[iql]->_weight1[ndim-2][qp];// weight
//             _fe[iql]->get_phi_gl_g(ndim-1,qp,_phi_g[iql]);   // global coord _phi_g
//
//
// #ifdef AXISYM   // axisymmetric  (index ->0)
//             interp_el_sol(_data_eq[iql].ub,0,DIMENSION,_phi_g[iql],elb_ndof[iql],_ub_g[iql]);
//             JxW_g[iql]  *=_ub_g[2][0];
//
// #endif
//             // local side loop (over the node face)
//             for (int lsi_node=0; lsi_node< elb_ndof[iql]; lsi_node++) {
//
//               // set up row i
//               const double phii_g=_phi_g[iql][lsi_node]; // boundary test function
//               int index= sur_toply[lsi_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM; // local element index
//               int bc_s=(int)_bc_bd[index];
//               int bc_v=(int)_bc_vol[index];
//               double dtxJxW_g=JxW_g[iql]*bc_v; // Neumann bc flag and Jac
//
//               // flag setup: +++++++++++++++++++++++++++++++++++++++
//               //  \int (bc_var*T+bc_val*val)ds
//               int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
//               int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
//               // ++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//               // Assemblying rhs ----------------------------
//               if (mode == 1){
//             FeM(index) += bc_val*(1-bc_var)*dtxJxW_g*phii_g*1.;
//        }
//
//               // Assemblying Matrix ---------------------------------
//               for (int lsj_node=0; lsj_node< elb_ndof[iql];  lsj_node++) {
//    int jndex=sur_toply[lsj_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM;
//                  KeM(index,jndex) += dtxJxW_g*bc_var*phii_g*_phi_g[iql][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
//               }// end j  ---------------------------------------
//
//             }// i   +++++++++++++++++++++++++++++++
//           } // end of the quadrature point qp-loop **********************
//         } // Neumann non homog
                      } //end for variables
                  }        // ---END LINEAR ---


              for ( int iql = 0; iql < 1; iql++ )  { // --- PIECEWISE ---
                  for ( int ivar = 0; ivar < _nvars[iql]; ivar++ )    {
                      // Dirichlet boundary conditions  ***********************************
                      double det = _fe[2]->JacSur ( elb_ngauss - 1, xxb_qnds, InvJac[2] ); // jacobian
                      double Ipenalty = det / _dt;                           // Dirichlet bc flag

                      if ( _bc_vol[ ( ivar + _nvars[1] + _nvars[2] ) *NDOF_FEM] == 0 ) {
                          int index = ivar * _el_dof[iql] + _nvars[1] * _el_dof[1] + _nvars[2] * _el_dof[2]; // local vol index

                          // Assemblying  Matrix & rhs
                          if ( mode == 1 ) {
                              FeM ( index ) += Ipenalty * 0.; //_data_eq[iql].ub[ivar*NDOF_FEM];
                              }

                          KeM ( index, index ) += Ipenalty; //  Dirichlet bc
                          }//end if on bc_vol
                      } //end loop variables
                  } // -- END PIECEWISE --


              } //end if side

          } // ======================  end for boundary =======================================

//     std::cout << KeM << "\n";
//     std::cout << FeM << "\n";
      /// e) Global assemblying energy equation
      A[Level]->add_matrix ( KeM, el_dof_indices );              // global matrix

      if ( mode == 1 ) {
          b[Level]->add_vector ( FeM, el_dof_indices ); // global rhs
          }

      } // end of element loop

  // clean
  el_dof_indices.clear();
  A[Level]->close();

//      A[Level]->print();
  if ( mode == 1 ) {
      b[Level]->close();
      }

#ifdef PRINT_INFO
  std::cout << " Matrix Assembled(DA)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

  return;
  }


///---------USER SOURCE---------------
// ============================================================================
/// This function reads Boundary conditions  from function
void MGSolDA::bc_read (
  int /*face_id_node*/,  ///<  face identity          (in)
  int  /*mat_flag*/,     ///<  volume identity         (in)
  double /*xp*/[],       ///< xp[] node coordinates    (in)
  int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
  int bc_flag[]          ///< boundary condition flag  (out)
) { // =========================================================================

//   for(int ivar=0; ivar<_nvars[1]+_nvars[2]; ivar++) {// lin+quad
  bc_Neum = 0;
  bc_flag = 0;
//  if(xp[0]>0.9 || xp[1]>0.9) {bc_Neum[ivar]=1;
//         bc_flag[ivar]=0;}
//   }
//   for(int ivar=_nvars[1]+_nvars[2]; ivar<_n_vars; ivar++) {// pw constant
//     bc_Neum[ivar]=0;    bc_flag[ivar]=0;
//   }
  }

double MGSolDA::MGFunctional ( double, double & ) {
std:
  cout << "Not implemented in MsolverDA";
  return 0.;

  }

void MGSolDA::ic_read (
  int /*face_id_node*/,
  int /*mat_id_elem*/,     // If using GAMBIT_INTERFACE, mat_flag given in Gambit for this node
  double /*xp*/[],   // Coordinates of the point
  int /*iel*/,  // element
  double ic[]    // Initial value of the [ivar] variable of the system
) {
  for ( int ivar = 0; ivar < _n_vars; ivar++ ) {
      ic[ivar] = 0.;
      }

  //xp[0]*(1.-xp[0])*xp[1]*(1.-xp[1]);
  }


void MGSolDA::ActivateVectField ( int Order, int Field, std::string SystemFieldName, int & n_index, int coupled ) {

  std::string FieldX =  SystemFieldName + "X";
  std::string FieldY =  SystemFieldName + "Y";
  std::string FieldZ =  SystemFieldName + "Z";

  if ( coupled == 0 ) {
      ActivateEquation ( Order, Field,     FieldX,  n_index );
      ActivateEquation ( Order, Field + 1,   FieldY,  n_index );
#if DIMENSION==3
      ActivateEquation ( Order, Field + 2,   FieldZ,  n_index );
#endif
      }
  else {   // flag 1 in SimulationConfiguration -> coupled
      _data_eq[Order].tab_eqs[Field] = n_index;                                  // table
      _data_eq[Order].mg_eqs[n_index] = _mgeqnmap.get_eqs ( SystemFieldName );           // FSI equation pointer
      _data_eq[Order].indx_ub[n_index + 1] = _data_eq[Order].indx_ub[n_index] + DIMENSION; // _data_eq[2].ub index
      _data_eq[Order].n_eqs++;                                               // number of quadratic system
      n_index++;   // update counter
      }

  return;
  }


// ================================================================================================
void MGSolDA::ActivateControl (
  int Order,
  int Field,               // Field to activate
  std::string SystemFieldName, // system field name
  int & n_index,                // n_index (collecting index)
  int vector,
  int dimension
) { //===============================================================================================
// // flag 2 in SimulationConfiguration.in ->  split
// // flag 1 in SimulationConfiguration -> coupled
// -------------------------------------------------------------------------------------------------
  std::string FieldX = SystemFieldName + "X";
  std::string FieldY = SystemFieldName + "Y";
  std::string FieldZ = SystemFieldName + "Z";

  ActivateEquation ( Order, Field, FieldX, n_index );

  if ( vector == 1 ) {
      ActivateEquation ( Order, Field + 1, FieldY, n_index );

      if ( dimension == 3 )
        ActivateEquation ( Order, Field + 2, FieldZ, n_index );
      }

  return;
  }


void MGSolDA::ActivateScalar ( int Order, int Field, std::string SystemFieldName, int & n_index ) {
  ActivateEquation ( Order, Field, SystemFieldName, n_index );
  return;
  }

void MGSolDA::ActivateCoupled ( int Order, int Field, std::string SystemFieldName, std::string SystemFieldName2, int & n_index ) {
  ActivateEquation ( Order, Field,   SystemFieldName,  n_index );
  ActivateEquation ( Order, Field + 1, SystemFieldName2, n_index );
  return;
  }

void MGSolDA::ActivateEquation ( int Order, int Field, std::string SystemFieldName, int & n_index ) {
  _data_eq[Order].tab_eqs[Field] = n_index;                               // table
  _data_eq[Order].mg_eqs[n_index] = _mgeqnmap.get_eqs ( SystemFieldName ); // Navier-Stokes equation pointer
  _data_eq[Order].indx_ub[n_index + 1] = _data_eq[Order].indx_ub[n_index] + 1; // _data_eq[2].ub index
  _data_eq[Order].n_eqs++;                                                // number of quadratic system
  n_index++;  // update counter
  return;
  }

