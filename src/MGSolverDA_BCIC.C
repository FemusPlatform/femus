// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ----------
// #include "Printinfo_conf.h"

// class files ----------------
#include "MGSolverDA.h"
// #include "MGSclass_conf.h"
// #include "MGFE_conf.h"
// local includes --------------
// #include "EquationSystemsExtendedM.h"

#include "MGGeomEl.h"
// #include "MGSystem.h"
#include "MeshExtended.h"
// #include "dense_set.h"

#include "MGEquationsSystem.h"
#include "MGFE.h"
#include "MGFEMap.h"
#include "MGGraph.h"
#include "MGUtils.h"

// #ifdef HAVE_MED
// #include "InterfaceFunctionM.h"
// #include "MEDCouplingFieldDouble.hxx"
// #include "MEDCouplingUMesh.hxx"
// #include "MEDLoader.hxx"
// #endif

// algebric includes -----------
// #include "dense_matrixM.h"
// #include "dense_vectorM.h"
// #include "numeric_vectorM.h"
// #include "sparse_MmatrixM.h"
// #include "sparse_matrixM.h"

// ==========================================================
//               MGSolDA functions
// ==========================================================

// ****************************************************************************
// Boundary condition and INITIAL CONDITION
// ****************************************************************************

// ============================================================================
/// This function generates the initial conditions for a Base system:
void MGSolDA::GenIc() {
  std::string input_dir = _mgutils._inout_dir;
  std::string ibc = _mgutils.get_file("IBC");
  std::ostringstream ibc_file;
  ibc_file << input_dir << ibc << ".h5";
  std::ifstream in((ibc_file.str()).c_str());

  if (!in) {  // ------- function initial condition -> user fuinction -> ic_read

    // Get a constant reference to the mesh objects.
    const double* xyz_glob = _mgmesh._xyz;               // global coordinates
    const int* map_nodes = _mgmesh._el_map[0];           // node map (elem,loc_nodes)->gl_nodes
    const int offset = _mgmesh._NoNodes[_NoLevels - 1];  // offset variables (on level _NoLevels-1)
    const int* off_el = _mgmesh._off_el[0];              // offset element
    int face_id_node = 0;
    int mat_id_elem = 0;
    NumericVectorM& sol_top = *x[_NoLevels - 1];  // solution (top level)
    NumericVectorM& old_sol_top = *x_old[0][_NoLevels - 1];
    NumericVectorM& nl_sol_top = *x_nonl[_NoLevels - 1];
    NumericVectorM& oold_sol_top = *x_old[1][_NoLevels - 1];
    const int* node_dof_top = &_node_dof[_NoLevels - 1][0];
    int ntot_elements = 0;
    for (int ilev = 0; ilev < _NoLevels; ilev++) { ntot_elements += _mgmesh._NoElements[0][ilev]; }

    // temp vect
    double* u_value = new double[_n_vars];
    double xp[DIMENSION];

    // loop reading
    int off_proc = _iproc * _NoLevels;
    int ndof_lev = 0;

    // Reading  bc_id *********************************************************
    // Open an existing file -------------------------------------------------
    std::ostringstream file_bc;
    file_bc << _mgutils._inout_dir << _mgutils.get_file("INMESH");  //"/mesh.h5";
#ifdef PRINT_INFO
    std::cout << " Reading bc_id from= " << file_bc.str() << std::endl;
#endif
    hid_t file_id = H5Fopen(file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hsize_t dims[2];

    // face id vector ---------------------------------------------------------
    int* face_id_vect;
    face_id_vect = new int[offset];

    for (int i = 0; i < offset; i++) { face_id_vect[i] = 0; }

    // Getting dataset
    std::ostringstream Name;
    Name << "NODES/COORD/BC";
    hid_t dtset = H5Dopen(file_id, Name.str().c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dtset); /* Get filespace handle first. */
    hid_t status = H5Sget_simple_extent_dims(filespace, dims, NULL);

    if (status < 0)
      std::cerr << "GenIc::read dims not found";
    else {
      assert((int)dims[0] == offset);
      status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &face_id_vect[0]);
    }

    H5Dclose(dtset);
    H5Sclose(filespace);

    // Reading  mat_id ********************************************************
    // mat id vector
    int* mat_id_vect = new int[ntot_elements];
    for (int i = 0; i < ntot_elements; i++) { mat_id_vect[i] = 1; }

    // Getting dataset
    std::ostringstream Name1;
    Name1 << "ELEMS/SUB/MAT";
    dtset = H5Dopen(file_id, Name1.str().c_str(), H5P_DEFAULT);
    filespace = H5Dget_space(dtset); /* Get filespace handle first. */
    status = H5Sget_simple_extent_dims(filespace, dims, NULL);

    if (status < 0) {
      std::cerr << "GenIc::read mat dims not found";
    } else
      status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat_id_vect[0]);

    H5Sclose(filespace);
    H5Dclose(dtset);
    H5Fclose(file_id);
    //  end Reading  mat_id ******************************************************

    for (int pr = 0; pr < _mgmesh._iproc; pr++) {
      int delta = off_el[pr * _NoLevels + _NoLevels - 1 + 1] - off_el[pr * _NoLevels + _NoLevels - 1];
      ndof_lev += delta;
    }

    for (int iel = 0; iel < off_el[off_proc + _NoLevels] - off_el[off_proc + _NoLevels - 1]; iel++) {
      int elem_gidx = (iel + off_el[_iproc * _NoLevels + _NoLevels - 1]) * NDOF_FEM;
      int elem_indx = (iel + ndof_lev) * _el_dof[0];
      mat_id_elem = mat_id_vect[iel + ndof_lev];

      // the local nodes
      for (int i = 0; i < NDOF_FEM; i++) {
        int k = map_nodes[elem_gidx + i];
        for (int idim = 0; idim < DIMENSION; idim++) { xp[idim] = xyz_glob[k + idim * offset]; }
        face_id_node = face_id_vect[k];
        // ===================================================
        // user definition reading function ----------------==
        ic_read(face_id_node, mat_id_elem, xp, iel + off_el[off_proc + _NoLevels - 1], u_value);
        // -------------------------------------------------==
        // ===================================================

        // Set discontinuous fields
        if (i == NDOF_FEM - 1)
          for (int ivar = 0; ivar < _nvars[0]; ivar++) {
            sol_top.set(
                node_dof_top[elem_indx + (ivar + _nvars[2] + _nvars[1]) * offset],
                u_value[_nvars[2] + _nvars[1] + ivar]);

            for (int kdof0 = 1; kdof0 < _el_dof[0]; kdof0++) {
              sol_top.set(node_dof_top[kdof0 + elem_indx + (ivar + _nvars[2] + _nvars[1]) * offset], 0.);
            }
          }

        // Set the quadratic and linear fields
        if (i < _el_dof[1])
          for (int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++) {
            sol_top.set(node_dof_top[k + ivar * offset], u_value[ivar]);
            //             printf(" %d %d %d %20g  *************************** \n", _iproc,
            //             k,ivar,u_value[ivar]);
          }
        else
          for (int ivar = 0; ivar < _nvars[2]; ivar++) {
            int irrr = node_dof_top[k + ivar * offset];
            sol_top.set(irrr, u_value[ivar]);
            //               printf(" %d %d %d %20g  *************quad ************** \n", _iproc,
            //               k,ivar,u_value[ivar]);
          }
      }
    }  // end of element loop

    //     // delocalize
    sol_top.localize(old_sol_top);
    sol_top.localize(nl_sol_top);
    sol_top.localize(oold_sol_top);

    // clean
    delete[] u_value;
#ifdef PRINT_INFO
    std::cout << "\n GenIc(DA): Initial solution defined by ic_read(xp,u_value)"
              << "\n \n";
#endif
    delete[] face_id_vect;
    delete[] mat_id_vect;
  } else {  // -------------------- file reading --> data_in/case.h5
    const int restart_lev = stoi(_mgutils._sim_config["restart_lev"]);  // restart label
    read_u(ibc_file.str(), restart_lev);
#ifdef PRINT_INFO
    std::cout << "\n GenIc(DA): Initial solution defined by " << ibc_file.str() << "\n \n";
#endif
  }  // +++++++++++++++++++++++++++++++++++++++++++++

  in.close();
  return;
}

// ============================================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::GenBc() {  // ========================================================================
  /// A) Set up: mesh,dof, bc
  // mesh ----------------------------------------------------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int ntot_elements = 0;

  for (int ilev = 0; ilev < _NoLevels; ilev++) { ntot_elements += _mgmesh._NoElements[0][ilev]; }

  // Dof ----------------------------------------------------------------------
  const int n_kb_dofs = ((_nvars[0] > 0) ? DIMENSION + 1 : 0);  // surface dofs
  const int n_pb_dofs =
      ((_nvars[1] > 0) ? _fe[1]->_NoShape[DIMENSION - 2] : 0);  // get_n_shapes(DIMENSION-2);
  const int n_ub_dofs =
      ((_nvars[2] > 0) ? _fe[2]->_NoShape[DIMENSION - 2] : 0);  // get_n_shapes(DIMENSION-2);
  //   const int  n_dofs =  n_pb_dofs*_nvars[1] + n_ub_dofs*_nvars[2];
  const int n_k_dofs = ((_nvars[0] > 0) ? DIMENSION + 1 : 0);                    // volume dofs
  const int n_l_dofs = ((_nvars[1] > 0) ? _fe[1]->_NoShape[DIMENSION - 1] : 0);  // get_n_shapes(DIMENSION-1);
  const int n_u_dofs = ((_nvars[2] > 0) ? _fe[2]->_NoShape[DIMENSION - 1] : 0);  // get_n_shapes(DIMENSION-1);

  // set 1 all the points for  bc (boundary condition) ------------------------
  for (int i1 = 0; i1 < _Dim[_NoLevels - 1]; i1++) {
    _bc[0][i1] = 1;
    _bc[1][i1] = 1;
  }

  // **************************************************************************
  // B) Reading  face_id vector (boundary zones) if the dataset exists
  // Open an existing file ----------------------------------------------------
  std::ostringstream file_bc;
  file_bc << _mgutils._inout_dir << _mgutils.get_file("INMESH");  //"/mesh.h5";
#ifdef PRINT_INFO
  std::cout << " Reading bc_id from= " << file_bc.str() << std::endl;
#endif
  hid_t file_id = H5Fopen(file_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dims[2];
  // face id vector initial setting (set 0) -----------------------------------
  int* face_id_vect = new int[offset];

  for (int i = 0; i < offset; i++) { face_id_vect[i] = 0; }

  // Getting dataset ----------------------------------------------------------
  std::ostringstream Name;
  Name << "NODES/COORD/BC";
  hid_t dtset = H5Dopen(file_id, Name.str().c_str(), H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dtset); /* Get filespace handle first. */
  hid_t status = H5Sget_simple_extent_dims(filespace, dims, NULL);

  if (status < 0) {
    std::cerr << "GenIc::read dims not found";
  } else {  // reading (otherwise it stays 0)
    assert((int)dims[0] == offset);
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &face_id_vect[0]);
  }

  H5Dclose(dtset);
  H5Sclose(filespace);
  // **************************************************************************
  // C) Reading  mat_id (volume zones) if the dataset exists
  // mat id vector initialization ---------------------------------------------
  int* mat_id_vect = new int[ntot_elements];

  //   int *mat_id_vect_lev; int icount=0;
  for (int i = 0; i < ntot_elements; i++) { mat_id_vect[i] = 1; }

  // level loop (in the file are written for each level)
  // Getting dataset --------------------------------------------------------
  std::ostringstream Name1;
  Name1 << "ELEMS/SUB/MAT";
  dtset = H5Dopen(file_id, Name1.str().c_str(), H5P_DEFAULT);

  filespace = H5Dget_space(dtset); /* Get filespace handle first. */
  status = H5Sget_simple_extent_dims(filespace, dims, NULL);

  if (status < 0) {
    std::cerr << "GenIc::read mat dims not found";
  } else {  // reading if the dataset exists (otherwise it stays 1) ------------
    status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat_id_vect[0]);
  }  // end reading ---------------------------------------------------------

  // clean --------------------------------------------------------------------
  H5Sclose(filespace);
  H5Dclose(dtset);
  H5Fclose(file_id);
  // *************************************************************************
  /// C  reading bc from function --> bc_read

  GenBc_loop(0, NDOF_FEM, n_u_dofs, n_l_dofs, n_k_dofs, face_id_vect, mat_id_vect);

  //   GenBc_loop(1,NDOF_FEMB,n_ub_dofs,n_pb_dofs,n_kb_dofs,
  //              face_id_vect,mat_id_vect);
  // clean
  delete[] face_id_vect;
  delete[] mat_id_vect;
#ifdef PRINT_INFO
  std::cout << "\n GenBc(DA): boundary conditions defined by  bc_read(xp,normal,bc_value) \n";
#endif

  return;
}

// ========================================================
/// This function  defines the boundary conditions for  DA systems:
void MGSolDA::GenBc_loop(
    const int vvvb,          // 0-> volume   1-> boundary
    const int ndof_femv,     // number of elment nodes
    const int n_u_dofs,      // element quad dofs
    const int n_l_dofs,      // element linear dofs
    const int /*n_k_dofs*/,  // element konst dofs
    int face_id_vect[],      // group bc from gambit
    int mat_id_vect[]        // bc from gambit

) {  // ======================================================
  const int n_dofs = _nvars[0] + _nvars[1] + _nvars[2];
  /// A) Set up: mesh,dof, bc
  // mesh
  const double* xyzgl = _mgmesh._xyz;
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  //   double normal[DIMENSION];
  double xp[DIMENSION];
  double xxb_qnds[NDOF_FEM * DIMENSION];
  double xx_qnds[NDOF_FEM * DIMENSION];
  double normal[DIMENSION];
  int el_neigh[NDOF_FEM];
  int el_flag[NDOF_FEM];
  int el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];
  int sur_toply[NDOF_FEMB];  // boundary topology
  //   int  *bc_Neu  =new int[n_dofs]; // element bc
  //   int  *bc_value=new int[n_dofs]; // element bc
  int bc_Neu[DIMENSION];    // new int[n_dofs]; // element bc
  int bc_value[DIMENSION];  // new int[n_dofs]; // element bc
  int face_id_node = 0;
  int mat_id_elem = 0;
  int ndof_lev = 0;
  const int el_sides = _mgmesh._GeomEl._n_sides[0];

  // initial setting -----------------------------------------------------------------------------
  for (int isub = 0; isub < _mgmesh._n_subdom; ++isub) {
    int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
    int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
    int delta = ielf - iel0;
    for (int iel = 0; iel < delta; iel++) {  // element loop
      _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds, isub);
      _mgmesh.get_el_neighbor(el_sides, 0, _NoLevels - 1, iel, el_neigh, isub);
      for (int i = 0; i < ndof_femv; i++) {                              // node lement loop
        const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + i];  // global node
        // coordinates
        //         for(int idim = 0; idim < DIMENSION; idim++) { xp[idim] = xyzgl[k + idim * offset]; }
        if (_node_dof[_NoLevels - 1][k] > -1) { _bc[0][_node_dof[_NoLevels - 1][k]] = -1000; }
      }
    }
  }  // --------------------------------------------------------------------------------------------

  /// B) Element Loop to set  bc[] (which is a node vector)
  for (int isub = 0; isub < _mgmesh._n_subdom; ++isub) {
    int iel0 = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub];
    int ielf = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * isub + 1];
    int delta = ielf - iel0;

    for (int iel = 0; iel < delta; iel++) {  // element loop
      _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds, isub);
      _mgmesh.get_el_neighbor(el_sides, 0, _NoLevels - 1, iel, el_neigh, isub);
      mat_id_elem = 0;  // v==1 (boundari) no mat_id
      if (vvvb == 0) { mat_id_elem = mat_id_vect[iel + iel0]; }

      // volume -----------------------------------------------------------------------------------
      for (int i = 0; i < ndof_femv; i++) {                              // node lement loop
        const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + i];  // global node

        if (_node_dof[_NoLevels - 1][k] > -1) {
          for (int idim = 0; idim < DIMENSION; idim++) {
            xp[idim] = xyzgl[k + idim * offset];
          }  // coordinates
          bc_value[0] = 1;
          bc_Neu[0] = 11;  // variable loop
          face_id_node = face_id_vect[k];
          // boundary (face_id_vect) and volume zones (mat_id_elem)
          bc_intern_read(face_id_node, mat_id_elem, xp, bc_Neu, bc_value);  // (idgroup, mat,x,.,.)
          int old_val = _bc[0][_node_dof[_NoLevels - 1][k]];
          if (old_val == -1000) _bc[0][_node_dof[_NoLevels - 1][k]] = bc_Neu[0];
          //           _bc[0][_node_dof[_NoLevels - 1][k]] = (old_val == -1000) ? bc_Neu[0] : (old_val);

          // sharing boundary nodes on the same element
          for (int ivar = 0; ivar < n_dofs; ivar++) {
            int kdof = _node_dof[_NoLevels - 1][k + ivar * offset];  // kdof <-k
            if (kdof > -1) {
              int number = _bc[1][kdof] / 10000;  // number of old count for double pt in an element
                                                  // multiple corner
              if (abs(number) == 1) {
                _bc[1][kdof] = _bc[1][kdof] - _bc[1][kdof] * 10000 /
                                                  abs(_bc[1][kdof]);  // if the count is 1 set 0 if >1 leave
              }
            }
          }
        }
      }  // i-loop
         // --------------------------------------------------------------------------------------------

      // -----------------------------  -----------------------------------------------------------------
      // Boundary
      for (int iside = 0; iside < el_sides; iside++) {
        if (el_neigh[iside] == -1) {
          // setup boundary element -> connectivity+coordinates
          for (int lbnode = 0; lbnode < NDOF_FEMB; lbnode++) {
            int lnode = _mgmesh._GeomEl._surf_top[lbnode + NDOF_FEMB * iside];  // local nodes
            sur_toply[lbnode] = lnode;                                          // lbnode -> lnode
            const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[lbnode]];

            // read from user function -> bc_read(   ) ----------------------------------------------------
            if (_node_dof[_NoLevels - 1][k] > -1) {
              face_id_node = face_id_vect[k];  // boundary (face_id_vect)
              for (int idim = 0; idim < DIMENSION; idim++) {
                xp[idim] = xxb_qnds[idim * NDOF_FEMB + lbnode] = xx_qnds[idim * NDOF_FEM + lnode];
              }
              bc_value[0] = 1;
              bc_Neu[0] = 11;
              bc_read(face_id_node, mat_id_elem, xp, bc_Neu, bc_value);
              _bc[0][_node_dof[_NoLevels - 1][k]] = bc_Neu[0];
            }
          }  // ----------------------------------------------------------------------------------------------

          // normal
          _fe[2]->normal_g(xxb_qnds, normal);
          int dir_maxnormal = (fabs(normal[0]) > fabs(normal[1])) ? 0 : 1;
          dir_maxnormal =
              (fabs(normal[dir_maxnormal]) > fabs(normal[DIMENSION - 1])) ? dir_maxnormal : DIMENSION - 1;
          // computation of k_el and face_id_el
          int k_el = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[NDOF_FEMB - 1]];  // global node
          int k_el_dof = _node_dof[_NoLevels - 1][k_el];
          int bc_el = (k_el_dof > -1) ? (int)_bc[0][k_el_dof] % 100 : 0;

          // ++++++++++++++++++++++++++++++++++++++++++++++++++++
          const int k_face =
              _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[NDOF_FEMB - 1]];  // global node
          int k0dof = _node_dof[_NoLevels - 1][k_face];
          int face_id_mid_face = (k0dof > -1) ? face_id_vect[k0dof] : 0;
          int bc_face = (k0dof > -1) ? (int)_bc[0][k0dof] : 0;
          bc_face = bc_face % 100;

          for (int i = 0; i < NDOF_FEMB; i++) {                                         // node lement loop
            const int k = _mgmesh._el_map[0][(iel + iel0) * ndof_femv + sur_toply[i]];  // global node
            int k0dof = _node_dof[_NoLevels - 1][k];

            if (k0dof > -1) {
              face_id_node = face_id_vect[k0dof];
              int bc_id = (int)_bc[0][k0dof];  // label surface pt 00
              int mynormal = dir_maxnormal;    // dir_normal=geometrical normal
              int sign = (bc_id == 0) ? 1 : (bc_id / (abs(bc_id)));
              if ((bc_id % 1000) / 100 > 0) { mynormal = (bc_id % 1000) / 100 - 1; }
              if ((bc_id % 1000) / 100 >
                  3) {  //   pts label x00 ----------------------------------------------------------
                std::cout << "\n ... single pressure pts .. \n";
                if (i < NDOF_PB) { _bc[1][_node_dof[_NoLevels - 1][k + DIMENSION * offset]] = 4; }
                mynormal = abs((bc_id % 1000) / 100) % 4;  // force the normal
              }
              bc_id = bc_id % 100;  // pts label x00 forced normal

              //   pts label 00 ----------------------------------------------------------
              // set the local boundary conditions into global vector bc[]
              for (int ivar = 0; ivar < _nvars[2]; ivar++) {  // quad el --
                int kdof = _node_dof[_NoLevels - 1][k + (ivar)*offset];
                int number = abs(_bc[1][kdof] / 10000) + 1;
                _bc[1][kdof] += sign * 10000;  // updating number of common nodes
                if (abs(_bc[1][kdof]) < 10000 || bc_id == bc_face) {
                  _bc[1][kdof] = sign * (abs(bc_id) + (mynormal + 1) * 1000 + number * 10000);
                }
              }  // ----------------------------------------------------------------
            }
          }  // loop i +++++++++++++++++++++++++++++++++++++++++++++
        }    // iside -1
      }      // -----------------------------  End Boundary -------------------------------------
    }        // end of element loop

    ndof_lev += delta;
  }  // i-sub
  return;
}

// ========================================
/// This function  defines the boundary conditions for the system:
void MGSolDA::bc_intern_read(
    int /*face_id_node*/,  ///<  face identity           (in)
    int /*mat_flag*/,      ///<  volume identity         (in)
    double /*xp*/[],       ///< xp[] node coordinates    (in)
    int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
    int bc_flag[]          ///< boundary condition flag  (out)
) {                        // ===================================
  /// Default: all Neumann
  //   for(int ivar=0; ivar<_n_vars; ivar++) {
  bc_flag[0] = 1;
  bc_Neum[0] = 11;

  //   }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void MGSolDA::get_el(
    const int Level,                   // level
    const int ivar0,                   // initial variable  <-
    const int nvars,                   // # of variables to get  <-
    const int el_nds,                  // # of element nodes for this variable  <-
    const int el_conn[],               // connectivity <-
    const int offset,                  // offset for connectivity <-
    std::vector<int>& el_dof_indices,  // element DOFs->
    int bc_dofs[][NDOF_FEM],           // element boundary cond flags ->
    double uold[]                      // element node values ->
    ) const {                          // ==============================================================
  for (int id = 0; id < el_nds; id++) {
    // quadratic -------------------------------------------------
    for (int ivar = ivar0; ivar < ivar0 + nvars; ivar++) {       // ivarq is like idim
      const int indx_loc = id + ivar * NDOF_FEM;                 // local (element) index
      const int indx_glob = el_conn[id] + ivar * offset;         // global (mesh) index
      const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

      el_dof_indices[indx_loc] = _node_dof[Level][indx_glob];  // from mesh to dof
      bc_dofs[0][indx_loc] = _bc[0][kdof_top];                 // element bc
      bc_dofs[1][indx_loc] = _bc[1][kdof_top];                 // element bc
      uold[indx_loc] = (*x_old[0][_NoLevels - 1])(kdof_top);   // element sol
    }  // end quadratic ------------------------------------------------
  }

  return;
}

// ==============================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void MGSolDA::get_el_dof_bc(
    const int Level,  // level
    const int iel,    // element subdomain
    //   const int nvars[],       // # of variables to get  <-
    const int el_nds[],                // # of element nodes for this variable  <-
    const int el_conn[],               // connectivity <-
    const int offset,                  // offset for connectivity <-
    std::vector<int>& el_dof_indices,  // element connectivity ->
    int bc_vol[],                      // element boundary cond flags ->
    int bc_bd[]                        // element boundary cond flags ->
    ) const {                          // ==============================================================

  for (int id = 0; id < NDOF_FEM; id++) {
    // quadratic -------------------------------------------------
    if (id < el_nds[2])
      for (int ivar = 0; ivar < _nvars[2]; ivar++) {  // ivarq is like idim
        const int indx_loc = id + ivar * NDOF_FEM;
        const int indx_loc_ql = id + ivar * el_nds[2];
        const int indx_glob = el_conn[id] + ivar * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level
        //         std::cout << el_dof_indices[indx_loc_ql] << " " << _node_dof[Level][indx_glob] << "\n";
        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
        bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
        bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
      }  // end quadratic ------------------------------------------------

    //     // linear -----------------------------
    if (id < el_nds[1])
      for (int ivar = 0; ivar < _nvars[1]; ivar++) {  // ivarq is like idim
        //        const int  indx_loc_l = id +ivar*el_nds[1];
        const int indx_loc = id + (ivar + _nvars[2]) * NDOF_FEM;
        const int indx_loc_ql = id + ivar * el_nds[1] + _nvars[2] * el_nds[2];
        const int indx_glob = el_conn[id] + (ivar + _nvars[2]) * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
        bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
        bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
      }  // end quadratic ------------------------------------------------

    //     // piecewise -----------------------------
    if (id < el_nds[0])
      for (int ivar = 0; ivar < _nvars[0]; ivar++) {  // ivarq is like idim
        //        const int  indx_loc_l = id +ivar*el_nds[1];
        const int indx_loc = id + (ivar + _nvars[2] + _nvars[1]) * NDOF_FEM;
        const int indx_loc_ql = id + ivar * el_nds[0] + _nvars[2] * el_nds[2] + _nvars[1] * el_nds[1];
        const int indx_glob = id + iel * el_nds[0] + (ivar + _nvars[2] + _nvars[1]) * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
        bc_bd[indx_loc] = _bc[1][kdof_top];                         // element bc
        bc_vol[indx_loc] = _bc[0][kdof_top];                        // element bc
      }  // end piecewise ------------------------------------------------
  }

  return;
}

/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element
void MGSolDA::set_el_dof_bc(
    const int Level,  // level
    const int iel,    // element subdomain
    //   const int nvars[],       // # of variables to get  <-
    const int el_nds[],                // # of element nodes for this variable  <-
    const int el_conn[],               // connectivity <-
    const int offset,                  // offset for connectivity <-
    std::vector<int>& el_dof_indices,  // element connectivity ->
    int bc_vol[],                      // element boundary cond flags ->
    int bc_bd[]                        // element boundary cond flags ->
    ) const {                          // ==============================================================

  for (int id = 0; id < NDOF_FEM; id++) {
    // quadratic -------------------------------------------------
    if (id < el_nds[2])
      for (int ivar = 0; ivar < _nvars[2]; ivar++) {  // ivarq is like idim
        const int indx_loc = id + ivar * NDOF_FEM;
        const int indx_loc_ql = id + ivar * el_nds[2];
        const int indx_glob = el_conn[id] + ivar * offset;
        const int kdof_top = _node_dof[_NoLevels - 1][indx_glob];  // dof from top level

        el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];  // from mesh to dof
        _bc[1][kdof_top] = bc_bd[indx_loc];                         // element bc
        _bc[0][kdof_top] = bc_vol[indx_loc];                        // element bc
      }  // end quadratic ------------------------------------------------

    //     // linear -----------------------------
    //       if ( id < el_nds[1] )    for ( int ivar = 0; ivar < _nvars[1]; ivar++ ) { //ivarq is like idim
    // //        const int  indx_loc_l = id +ivar*el_nds[1];
    //             const int indx_loc = id + ( ivar + _nvars[2] ) * NDOF_FEM;
    //             const int indx_loc_ql = id + ivar * el_nds[1] + _nvars[2] * el_nds[2];
    //             const int indx_glob = el_conn[id] + ( ivar + _nvars[2] ) * offset;
    //             const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level
    //
    //             el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
    //             bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
    //             bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
    //             } // end quadratic ------------------------------------------------
    //
    //
    //       //     // piecewise -----------------------------
    //       if ( id < el_nds[0] )    for ( int ivar = 0; ivar < _nvars[0]; ivar++ ) { //ivarq is like idim
    // //        const int  indx_loc_l = id +ivar*el_nds[1];
    //             const int indx_loc = id + ( ivar + _nvars[2] + _nvars[1] ) * NDOF_FEM;
    //             const int indx_loc_ql = id + ivar * el_nds[0] + _nvars[2] * el_nds[2] + _nvars[1] *
    //             el_nds[1]; const int indx_glob = id + iel * el_nds[0] + ( ivar + _nvars[2] + _nvars[1] ) *
    //             offset; const int kdof_top = _node_dof[_NoLevels - 1][indx_glob]; // dof from top level
    //
    //             el_dof_indices[indx_loc_ql] = _node_dof[Level][indx_glob];    //from mesh to dof
    //             bc_bd[indx_loc]         = _bc[1][kdof_top];                    // element bc
    //             bc_vol[indx_loc]        = _bc[0][kdof_top];                    // element bc
    //             } // end piecewise ------------------------------------------------
  }

  return;
}

///---------USER SOURCE---------------
// ============================================================================
/// This function reads Boundary conditions  from function
void MGSolDA::bc_read(
    int /*face_id_node*/,  ///<  face identity          (in)
    int /*mat_flag*/,      ///<  volume identity         (in)
    double /*xp*/[],       ///< xp[] node coordinates    (in)
    int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
    int bc_flag[]          ///< boundary condition flag  (out)
) {                        // =========================================================================

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

// ===========================================================================================
void MGSolDA::ic_read(
    int /*face_id_node*/,
    int /*mat_id_elem*/,  // If using GAMBIT_INTERFACE, mat_flag given in Gambit for this node
    double /*xp*/[],      // Coordinates of the point
    int /*iel*/,          // element
    double ic[]           // Initial value of the [ivar] variable of the system
) {
  for (int ivar = 0; ivar < _n_vars; ivar++) { ic[ivar] = 0.; }
  // xp[0]*(1.-xp[0])*xp[1]*(1.-xp[1]);
}
