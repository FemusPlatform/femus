// #include "Equations_tab.h"
// #include "Equations_conf.h"
// ===============================
// ===============================
// std lib
// #include <cmath>
#include <iomanip>
// #include <limits>
// #include <sstream>

// configuration files ----------
#include "Printinfo_conf.h"

// class files ----------------
#include "MGSolverDA.h"
// #include "MGSclass_conf.h"

// local includes --------------
// #include "EquationSystemsExtendedM.h"
// #include "MGFE_conf.h"
#include "MGGeomEl.h"
#include "MGSystem.h"
#include "MeshExtended.h"
// #include "dense_set.h"

// #include "MGEquationsSystem.h"
#include "MGFE.h"
// #include "MGFEMap.h"
// #include "MGGraph.h"
// #include "MGUtils.h"

#ifdef HAVE_MED
// #include "InterfaceFunctionM.h"
// #include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"
#endif

// algebric includes -----------
// #include "dense_matrixM.h"
// #include "dense_vectorM.h"
#include "numeric_vectorM.h"
#include "sparse_MmatrixM.h"
#include "sparse_matrixM.h"

// ==========================================================
//               MGSolDA functions
// ==========================================================

//  ******************************************************************************************
//  PRINT AND READ
//  ******************************************************************************************

// ================================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_bc(std::string namefile, const int Level) {
  //  setup
  const int n_nodes = _mgmesh._NoNodes[Level];
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int* sol = new int[n_nodes];

  // file hdf5
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hsize_t dimsf[2];
  dimsf[0] = n_nodes;
  dimsf[1] = 1;

  // loop over all systems ---------------------------
  // print quad -------------------------------------
  for (int ivar = 0; ivar < _nvars[2]; ivar++) {
    std::string var_name = _var_names[ivar] + "bd";
    for (int i = 0; i < n_nodes; i++) { sol[i] = _bc[1][_node_dof[_NoLevels - 1][i + (ivar)*offset]]; }
    _mgutils.print_Ihdf5(file_id, var_name, dimsf, sol);
    std::string var_name2 = _var_names[ivar] + "vl";
    for (int i = 0; i < n_nodes; i++) { sol[i] = _bc[0][_node_dof[_NoLevels - 1][i + ivar * offset]]; }
    _mgutils.print_Ihdf5(file_id, var_name2, dimsf, sol);
  }

  // print linear -----------------------------------
  double* sol_c = new double[_mgmesh._GeomEl.n_l[0]];

  // put the attention to the end of the for cycle (_nvars[1]+_nvars[2])
  for (int ivar = _nvars[2]; ivar < _nvars[1] + _nvars[2]; ivar++) {
    std::string var_name = _var_names[ivar] + "bd";

    // class writing routine
    // 2bA projection of the fine nodes -------------------------
    for (int isubdom = 0; isubdom < _mgmesh._n_subdom; isubdom++) {
      for (int i = _mgmesh._off_nd[0][isubdom * _NoLevels];
           i < _mgmesh._off_nd[0][isubdom * _NoLevels] + _mgmesh._off_nd[1][Level + 1 + isubdom * _NoLevels] -
                   _mgmesh._off_nd[1][isubdom * _NoLevels];
           i++) {
        sol[i] = _bc[1][_node_dof[_NoLevels - 1][i + ivar * offset]];
      }
    }  // 2bA end proj fine grid ----------------------------------

    //  2bB element interpolation over the fine mesh -----------------------
    for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
      for (int iel = _mgmesh._off_el[0][_NoLevels * iproc + Level];
           iel < _mgmesh._off_el[0][Level + _NoLevels * iproc + 1]; iel++) {
        int indx = iel * NDOF_FEM;

        // element iel ----------------------------
        // vertices
        for (int in = 0; in < NDOF_P; in++) { sol_c[in] = sol[_mgmesh._el_map[0][indx + in]]; }

        for (int in = 0; in < NDOF_FEM; in++) {  // mid-points
          double sum = 0;

          for (int jn = 0; jn < NDOF_P; jn++) { sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn]; }

          sol[_mgmesh._el_map[0][indx + in]] = sum;
        }
      }  // ---- end iel -------------------------
    }    // 2bB end interpolation over the fine mesh ------------------------

    _mgutils.print_Ihdf5(file_id, var_name, dimsf, sol);

    std::string var_name2 = _var_names[ivar] + "vl";

    // class writing routine
    // 2bA projection of the fine nodes -------------------------
    for (int isubdom = 0; isubdom < _mgmesh._n_subdom; isubdom++) {
      for (int i = _mgmesh._off_nd[0][isubdom * _NoLevels];
           i < _mgmesh._off_nd[0][isubdom * _NoLevels] + _mgmesh._off_nd[1][Level + 1 + isubdom * _NoLevels] -
                   _mgmesh._off_nd[1][isubdom * _NoLevels];
           i++) {
        sol[i] = _bc[0][_node_dof[_NoLevels - 1][i + ivar * offset]];
      }
    }  // 2bA end proj fine grid ----------------------------------

    //  2bB element interpolation over the fine mesh -----------------------
    for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
      for (int iel = _mgmesh._off_el[0][_NoLevels * iproc + Level];
           iel < _mgmesh._off_el[0][Level + _NoLevels * iproc + 1]; iel++) {
        int indx = iel * NDOF_FEM;

        // element iel ----------------------------
        // vertices
        for (int in = 0; in < NDOF_P; in++) { sol_c[in] = sol[_mgmesh._el_map[0][indx + in]]; }

        for (int in = 0; in < NDOF_FEM; in++) {  // mid-points
          double sum = 0;

          for (int jn = 0; jn < NDOF_P; jn++) { sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn]; }

          sol[_mgmesh._el_map[0][indx + in]] = sum;
        }
      }  // ---- end iel -------------------------
    }    // 2bB end interpolation over the fine mesh ------------------------

    _mgutils.print_Ihdf5(file_id, var_name2, dimsf, sol);

  }  // ivar

  // clean -------------------------------
  delete[] sol_c;
  delete[] sol;
  H5Fclose(file_id);

  return;
}

// ===========================================================================================
/// This function reads the MSolDA system solution from namefile.h5
void MGSolDA::read_u(
    std::string namefile,  ///< filename (with path)
    int Level_restart      ///< restart Level-1 flag
) {  //========================================================================================

  /// A. Setup
  // mesh
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];                          // fine level # of nodes
  const int pie_offset = (NSUBDOM) * (_mgmesh._NoElements[0][_NoLevels - 1]);  // fine level # of nodes
  // field
  double* sol;
  double* sol_pie;
  hid_t file_sol;
  // ***********************************************************************************************
  // ***********************************************************************************************
  // reading from different level -> RESU_COARSE
  // ***********************************************************************************************
  // ***********************************************************************************************
  /// B. Read from  Level-1  (NoLevels-2)
  if (Level_restart != 0) {
    const int ndigits = stoi(_mgutils._sim_config["ndigits"]);  // digit for namefile
    const int restart = stoi(_mgutils._sim_config["restart"]);  // restart label

    std::string resu_dir = _mgutils._inout_dir;  // mesh aux dir
    std::ostringstream namefile_sol;             //  namefile_sol
    namefile_sol << resu_dir << "/sol.msh1." << setw(ndigits) << setfill('0') << restart << ".h5";
    file_sol = H5Fopen(namefile_sol.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);  // hdf5 file

    std::ostringstream namefile_mesh;  //   namefile_mesh
    namefile_mesh << "RESU_COARSE/mesh.msh1.h5";
    hid_t file_mesh = H5Fopen(namefile_mesh.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    int Level = _NoLevels - 2;                         // Reading from Level
    const int n_elem = _mgmesh._NoElements[0][Level];  // coarse level # of elements
    const int offset_lev = _mgmesh._NoNodes[Level];    // coarse level # of nodes
    int* map_f2c = new int[offset];                    // node map fine to coarse

    std::cout << " Reading Level-1 sol from file " << namefile_sol.str() << "\n";
    std::cout << " Reading Level-1 connectivity from file " << namefile_mesh.str() << "\n";

    // file to read
    double sol_c[NDOF_FEM];
    double sol_f[NDOF_FEM];        // element connectivity
    sol = new double[offset_lev];  // coarse solution     (sol.xxx.h5)
    sol_pie = new double[1];

    std::ostringstream namedir_conn;  // coarse connectivity (-> mesh.h5)
    namedir_conn << "/ELEMS/FEM0/MSH" << Level;
    int* map_coarse = new int[NDOF_FEM * n_elem];
    _mgutils.read_Ihdf5(file_mesh, namedir_conn.str().c_str(), map_coarse);

    int offel_coarse = 0;  // # of element at Level for proc < _iproc

    for (int jproc = 0; jproc < _iproc; jproc++) {
      offel_coarse +=
          _mgmesh._off_el[0][jproc * _NoLevels + Level + 1] - _mgmesh._off_el[0][jproc * _NoLevels + Level];
    }
    int end_elem = _mgmesh._off_el[0][_iproc * _NoLevels + Level + 1];
    int start_elem = _mgmesh._off_el[0][_iproc * _NoLevels + Level];
    for (int iel = 0; iel < end_elem - start_elem; iel++) {
      int elem_gidx = (iel + _mgmesh._off_el[0][_iproc * _NoLevels + Level]) * NDOF_FEM;
      int elem_coarse_gidx = (iel + offel_coarse) * NDOF_FEM;
      for (int in = 0; in < NDOF_FEM; in++) {
        map_f2c[_mgmesh._el_map[0][elem_gidx + in]] = map_coarse[elem_coarse_gidx + in];
      }
    }

    // reading loop over system varables
    for (int ivar = 0; ivar < _n_vars; ivar++) {
      int el_nds = NDOF_FEM;
      if (ivar >= _nvars[2]) el_nds = NDOF_P;  // quad and linear

      // reading ivar param
      _mgutils.read_Dhdf5(
          file_sol, "/" + _var_names[ivar], sol);  // reading coarse quad solution (sol.xxx.h5)
      double Irefval = 1. / _refvalue[ivar];       // units

      // storing  ivar variables (in parallell)
      int end_elem = _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels];
      int start_elem = _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1];
      for (int iel = 0; iel < end_elem - start_elem; iel++) {
        int elem_gidx = (iel + start_elem) * NDOF_FEM;
        // linear variables ----------------------------------------------------------------------
        // vertex storage
        for (int in = 0; in < NDOF_P; in++) {  //  only vertices
          int k_fine = _mgmesh._el_map[0][elem_gidx + in];
          sol_c[in] = sol[map_f2c[k_fine]];
        }
        // interpolation over quad field
        for (int in = 0; in < el_nds; in++) {  // all element points
          double sum = 0;
          for (int jn = 0; jn < NDOF_P; jn++) { sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn]; }
          sol_f[in] = sum;
        }  // -------------------------------------------------------------------------------------
        // quad fields ---------------------------------------------------------------------------
        for (int i = 0; i < el_nds; i++) {                               // linear and quad
          const int k = _mgmesh._el_map[0][elem_gidx + i];               // the global node
          const int kdof = _node_dof[_NoLevels - 1][k + ivar * offset];  // dof
          x[_NoLevels - 1]->set(kdof, sol_f[i] * Irefval);  // sol[k]*Irefval);    // set the field
        }  // --------------------------------------------------------------------------------------
      }    // int iel
    }      // int ivar

    x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);

    // clean
    H5Fclose(file_mesh);
    delete[] map_coarse;
  }
  // ***********************************************************************************************
  // ***********************************************************************************************
  /// C. Read from  the same Level (NoLevels-1) -> RESU dir
  // ***********************************************************************************************
  // ***********************************************************************************************
  else {
    int SolToRead = _NumRestartSol;
    if (_NumRestartSol > 3) {
      std::cout << " _NumRestartSol>3 \n";
      abort();
    }
    std::string SolSuffix[3];
    SolSuffix[0] = "";
    SolSuffix[1] = "_old";
    SolSuffix[2] = "_oold";

    // file to read
    sol = new double[offset];
    sol_pie = new double[pie_offset];  // temporary vector
    file_sol = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    const int end = _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels];
    const int start = _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1];

    for (int oldTstep = 0; oldTstep < SolToRead; oldTstep++) {
      int istep = SolToRead - oldTstep - 1;

      // linear and quad variable  over system varables --------------------------------------
      for (int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++) {
        int el_nds = NDOF_FEM;
        if (ivar >= _nvars[2]) el_nds = NDOF_P;  // quad and linear
        // reading
        _mgutils.read_Dhdf5(file_sol, "/" + _var_names[ivar] + SolSuffix[istep], sol);
        double Irefval = 1. / _refvalue[ivar];  // units
        // in parallel over elements
        for (int iel = 0; iel < end - start; iel++) {
          int elem_gidx = (iel + _mgmesh._off_el[0][_iproc * _NoLevels + _NoLevels - 1]) * NDOF_FEM;
          for (int i = 0; i < el_nds; i++) {            // linear and quad
            int k = _mgmesh._el_map[0][elem_gidx + i];  // the global node
            x[_NoLevels - 1]->set(
                _node_dof[_NoLevels - 1][k + ivar * offset], sol[k] * Irefval);  // set the field
          }
        }
      }  // ---------------------------------------------------------------------------------

      int ndof_lev = 0;  // first nodeo fthe level
      for (int pr = 0; pr < _mgmesh._iproc; pr++) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] -
                    _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
        ndof_lev += delta;
      }
      // Piecewise variables --------------------------------------------------------------
      for (int ivar = _nvars[2] + _nvars[1]; ivar < _n_vars; ivar++) {
        int el_nds = NDOF_FEM;
        el_nds = NSUBDOM;  // quad and linear
        // reading ivar param
        _mgutils.read_Dhdf5(file_sol, "/" + _var_names[ivar], sol_pie);
        double Irefval = 1. / _refvalue[ivar];  // units
        // storing  ivar variables (in parallell)
        for (int iel = 0; iel < end - start; iel++) {
          int k = (iel + ndof_lev);  // the global element
          x[_NoLevels - 1]->set(
              _node_dof[_NoLevels - 1][k + ivar * offset], sol_pie[k * el_nds] * Irefval);  // set the field
        }
      }  //-------------------------------------------------------------------------------

      x[_NoLevels - 1]->localize(*x_old[istep][_NoLevels - 1]);

    }  // end  else {
  }    // end if

  /// D. delocalization and clean

  H5Fclose(file_sol);
  delete[] sol;
  delete[] sol_pie;
  return;
}

/// ==================================================
///      ----------OPERATORS  READING-------------
///  -------------------------------------------------

/// This function read the Matrix linear-quad operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadMatrix(
    const int Level,              // Level
    const std::string& namefile,  // file name
    SparseMatrixM& Mat,           // Matrix to read
    const int nvars_in[]          // # quad variables
    //                            const int nvars_l           // # linear variables
) {  // -p----------------------------------------------------

#ifdef PRINT_INFO  // ------------------------------------------
  std::cout << " ReadMatrix(DA): start matrix reading    " << std::endl;
#endif
  int NoLevels = _mgmesh._NoLevels;  // max number of levels

  ///-------------------------------------------------------------------------------
  /// [1] Reading matrix from file hdf5---------------------------------------------
  ///-------------------------------------------------------------------------------

  // reading matrix dimensions
  // initialize dimension
  int dim_qlk[3][3], ldim[2];  // [*][*] row variable type column variable type

  for (int aa = 0; aa < 3; aa++)
    for (int bb = 0; bb < 3; bb++) { dim_qlk[aa][bb] = 0; }

  for (int iql = 0; iql < 3; iql++) {
    for (int jql = 0; jql < 3; jql++) {
      std::ostringstream mode;
      mode << (iql)*10 + (jql)*1;
      int status = Mat.read_dim_hdf5(namefile.c_str(), mode.str().c_str(), ldim);  // reading dimensions

      if (status == 0) {
        dim_qlk[iql][jql] = ldim[0];  // if "mode" found -> put the row number in dim_qlk[][]
      }
    }
  }

  // setup len and len_off vectors
  int** length_row = new int*[9];     // total length entry sparse struct
  int** length_offrow = new int*[9];  // offdiag entry sparse matrix struct

  // length and offlength structure
  //    [0]kk [1]kl [2]kq
  //    [3]lk [4]ll [5]lq
  //    [6]qk [7]ql [8]qq

#ifdef HAVE_LASPACKM
  int** pos_row = new int*[9];  // pos matrix  struct
#endif

  // reading diagonal row length and off row length (0=k,1=l,2=q)
  for (int iql = 0; iql < 3; iql++) {      // cycle on the row blocks
    if (nvars_in[iql] != 0) {              // check if the matrix is needed for this row block
      for (int jql = 0; jql < 3; jql++) {  // cycle on the column blocks
        if (nvars_in[jql] != 0) {          // check if the matrix is needed for this column block
          length_row[iql * 3 + jql] = new int[dim_qlk[iql][jql] + 1];  // prepare vectors of correct dimension
          length_offrow[iql * 3 + jql] = new int[dim_qlk[iql][jql] + 1];
          std::ostringstream mode;
          mode << (iql)*10 + (jql)*1;
          Mat.read_len_hdf5(
              namefile.c_str(), mode.str().c_str(), length_row[iql * 3 + jql],
              length_offrow[iql * 3 + jql]);  // reading sparse struct
#ifdef HAVE_LASPACKM
          // reading compressed column index
          pos_row[ql] = new int[length_row[ql][dim_ql[ql]]];     //  matrix pos
          Mat.read_pos_hdf5(namefile.c_str(), ql, pos_row[ql]);  // reading sparse pos
#endif
        } else {  // if there is no variable of jlq type assign -1 in the length vector
          length_row[iql * 3 + jql] = new int[1];
          length_row[iql * 3 + jql][0] = -1;
          length_offrow[iql * 3 + jql] = new int[1];
          length_offrow[iql * 3 + jql][0] = -1;
        }
      }  // end cycle on the column blocks
    }    // if there is no variable of ilq type assign -1 in the length vector for all the column blocks
    else {
      for (int kql = 0; kql < 3; kql++) {
        length_row[iql * 3 + kql] = new int[1];
        length_row[iql * 3 + kql][0] = -1;
        length_offrow[iql * 3 + kql] = new int[1];
        length_offrow[iql * 3 + kql][0] = -1;
      }
    }
  }  // end cycle on the row blocks

  ///-------------------------------------------------------------------------------------
  /// [2] Set up initial matrix positions for levels, processors and type of variables----
  ///-------------------------------------------------------------------------------------
  // mesh index vector
  int* off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0];  // quad offset index(subdomains and levels)
  off_nd[1] = _mgmesh._off_nd[1];  // linear offset index (subdomains and levels)

  // prepare row indeces
  int ml[3];
  int m = 0;

  for (int l = 0; l < 3; l++) {
    m += nvars_in[l] * dim_qlk[l][l];  // global last row
  }

  int n = m;  // global last column
  ml[0] =
      (_mgmesh._off_el[0][_iproc * NoLevels + Level + 1] - _mgmesh._off_el[0][_iproc * NoLevels + Level]) *
      _el_dof[0];                                                                   // const local row
  ml[1] = off_nd[1][_iproc * NoLevels + Level + 1] - off_nd[1][_iproc * NoLevels];  // linear local row
  ml[2] = off_nd[0][_iproc * NoLevels + Level + 1] - off_nd[0][_iproc * NoLevels];  // quad local row
  int m_l = ml[0] * nvars_in[0] + ml[1] * nvars_in[1] + ml[2] * nvars_in[2];  // quad+lin+const last local row
  int n_l = m_l;  // quad+lin+const last local column

  // setting levels
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev_c = (Level - 1 + _NoLevels + 1) % (_NoLevels + 1);
  int lev[3];
  lev[0] = Level;
  lev[1] = lev_c;
  lev[2] = Level;

  // initial local matrix positions
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0;  // initialize ml_init

  for (int isubd = 0; isubd < _iproc; isubd++) {
    ml_init[0] +=
        (_mgmesh._off_el[0][isubd * NoLevels + Level + 1] - _mgmesh._off_el[0][isubd * NoLevels + Level]) *
        _el_dof[0];                                                                         // const
    ml_init[1] += off_nd[1][Level + 1 + isubd * _NoLevels] - off_nd[1][isubd * _NoLevels];  // linear
    ml_init[2] += off_nd[0][Level + 1 + isubd * _NoLevels] - off_nd[0][isubd * _NoLevels];  // quad
  }

  int top_node_iproc = off_nd[0][_iproc * NoLevels];
  int ml_start = _node_dof[Level][top_node_iproc];

  ///-------------------------------------------------------------------------------------
  /// [3] Write matrix sparse structure---------------------------------------------------
  ///-------------------------------------------------------------------------------------

#ifdef HAVE_LASPACKM  // ------------------- only laspack ----------------

  // Graph dimension
  Graph graph;
  graph.resize(m);
  graph._m = m;
  graph._n = n;
  graph._ml = m_l;
  graph._nl = n_l;

  graph._ml_start = ml_start;

  // quadratic
  for (int ivar = 0; ivar < nvar[0]; ivar++) {
    for (int i = ml_init[0]; i < ml_init[0] + ml[2]; i++) {
      // row index
      len[0] = (length_row[0][i + 1] - length_row[0][i]);
      int irowl =
          _node_dof[Level][_mgmesh._node_map[lev[2]][i] + (ivar + 0 * nvars_in[2]) * offset] - ml_start;
      len[1] = (length_row[1][i + 1] - length_row[1][i]);

      graph[irowl].resize(nvars_in[2] * len[0] + nvars_in[1] * len[1] + 1);

      //       int  jlq=0;
      for (int jlq = 0; jlq < 2; jlq++) {  // ilq=0 quadratic-quadratic ilq=1 quadratic-linear
        for (int jvar = 0; jvar < nvars_in[jlq]; jvar++) {
          for (int j = 0; j < len[jlq]; j++) {
            graph[irowl][j + jvar * len[jlq] + jlq * nvar[0] * len[0]] =
                _node_dof[Level]
                         [_mgmesh._node_map[lev[jlq]][pos_row[jlq][j + length_row[jlq][i]]] +
                          (jvar + jlq * nvar[0]) * offset];
          }
        }

        // last stored value is the number of in-matrix nonzero values
        graph[irowl][nvar[0] * len[0] + nvar[1] * len[1]] =
            nvar[0] * (length_offrow[0][i + 1] - length_offrow[0][i]) +
            nvar[1] * (length_offrow[1][i + 1] - length_offrow[1][i]);
      }
    }  // end quadratic ++++++++++++++++++++++++++++++++++++
  }

  // linear part ++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int ivar = 0; ivar < nvar[1]; ivar++) {
    for (int i = ml_init[1]; i < ml_init[1] + ml[1]; i++) {
      // row index
      len[3] = (length_row[3][i + 1] - length_row[3][i]);
      int irowl = _node_dof[Level][_mgmesh._node_map[lev_c][i] + (ivar + nvar[0]) * offset];
      len[2] = (length_row[2][i + 1] - length_row[2][i]);

      graph[irowl].resize(nvar[1] * len[3] + nvar[0] * len[2] + 1);

      for (int ilq = 2; ilq < 4; ilq++) {  // ilq=3 linear-linear ilq=2 linear-quadratic
        for (int jvar = 0; jvar < nvar[ilq - 2]; jvar++) {
          for (int j = 0; j < len[ilq]; j++) {
            graph[irowl][j + jvar * len[ilq] + (ilq - 2) * nvar[0] * len[2]] =
                _node_dof[Level]
                         [_mgmesh._node_map[lev[ilq - 2]][pos_row[ilq][j + length_row[ilq][i]]] +
                          (jvar + (ilq - 2) * nvar[0]) * offset];
          }
        }
      }

      graph[irowl][nvar[1] * len[3] + nvar[0] * len[2]] =
          nvar[1] * (length_offrow[3][i + 1] - length_offrow[3][i]) +
          nvar[0] * (length_offrow[2][i + 1] - length_offrow[2][i]);
    }
  }  // end linear

  // Update sparsity pattern -----------------------------
  Mat.update_sparsity_pattern(graph);
  graph.clear();
#endif              // ------------------ ----------------------------------------
#ifdef HAVE_PETSCM  // -------------------  only petsc ----------------

  int cont_var = 0;            // shifter for variables of all type
  std::vector<int> n_dz(m_l);  // # diagonal entries
  std::vector<int> n_oz(m_l);  // # offset entries

  for (int i = 0; i < m_l; i++) {
    n_dz[i] = 0;  // initialize to zero
    n_oz[i] = 0;
  }

  for (int ilq = 2; ilq >= 0; ilq--) {                  // cycle on type of variable (2=q,1=l,0=k) - row block
    for (int ivar = 0; ivar < nvars_in[ilq]; ivar++) {  // cycle on variables of ilq type
      for (int i1 = 0; i1 < ml[ilq]; i1++) {            // cycle on the rows
        int i = ml_init[ilq] + i1;
        int i_top;

        if (ilq != 0) {
          i_top = _mgmesh._node_map[lev[ilq]][i];
        } else {
          i_top = i;  // for piecewise only
        }

        int ind = _node_dof[Level][i_top + (ivar + cont_var) * offset] -
                  ml_start;  // ind is local (per proc) row index

        for (int jlq = 2; jlq >= 0; jlq--) {  // cycle on type of variable (2=q,1=l,0=k) - column block
          int jlq1 = (ilq + jlq) % 3;

          for (int jvar = 0; jvar < nvars_in[jlq1]; jvar++) {  // cycle on variables of jlq1 type
            n_dz[ind] += (length_row[ilq * 3 + jlq1][i + 1] - length_row[ilq * 3 + jlq1][i]);
            n_oz[ind] += (length_offrow[ilq * 3 + jlq1][i + 1] - length_offrow[ilq * 3 + jlq1][i]);
          }  // jlq1 variables loop
        }    // column block loop
      }      // row loop
    }        // ilq variables loop

    cont_var +=
        nvars_in[ilq];  // shift in the _node_dof for all the variables written till now (q, then l, then k)
  }

  for (int i = 0; i < m_l; i++) {
    n_dz[i] -= n_oz[i];  // count correctly diagonal entries
  }

  // Update sparsity pattern -----------------------------
  Mat.update_sparsity_pattern(m, n, m_l, n_l, n_dz, n_oz);

#endif

  //  clean ----------------------------------------------------
  n_dz.clear();
  n_oz.clear();

  for (int ql = 0; ql < 9; ql++) {
#ifdef HAVE_LASPACKM
    delete[] pos_row[ql];  // sparse compressed postion
#endif
    delete[] length_row[ql];
    delete[] length_offrow[ql];
  }

#ifdef HAVE_LASPACKM
  delete[] pos_row;  // sparse compressed postion
#endif
  delete[] length_row;
  delete[] length_offrow;

#ifdef PRINT_INFO  // ------------------------------------------
  std::cout << " ReadMatrix(DA): end matrix reading    " << std::endl;
#endif
}

// =================================================================
/// This function read the Prolongation operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadProl(
    const int Level,          // Level
    const std::string& name,  // file name (reading from)
    SparseMMatrixM& Prol,     // Prolongation matrix
    const int nvars_in[],     // # quad, linear and const variables
    int node_dof_f[], int node_dof_c[]) {
#ifdef PRINT_INFO  // ------------------------------------------
  std::cout << " ReadProl(DA): start prolongator reading    " << std::endl;
#endif
  ///-------------------------------------------------------------------------------
  /// [1] Reading prolongator from file hdf5----------------------------------------
  ///-------------------------------------------------------------------------------

  // reading dimensions
  int ldim[2];
  int dim_qlk[3][2];  // [0][0] piecewise fine, [2][1] quad coarse

  for (int aa = 0; aa < 3; aa++)
    for (int bb = 0; bb < 2; bb++) { dim_qlk[aa][bb] = 0; }

  for (int jql = 0; jql < 3; jql++) {
    std::ostringstream mode;
    mode << jql;
    int status = Prol.read_dim_hdf5(name.c_str(), mode.str().c_str(), ldim);

    if (status == 0) {
      dim_qlk[jql][0] = ldim[0];
      dim_qlk[jql][1] = ldim[1];
    }
  }

  // setup length (0=k,1=l,2=q)
  int** length_row = new int*[3];      // total length entry sparse struct
  int** length_off_row = new int*[3];  // offdiag entry sparse matrix struct
  int** pos_row = new int*[3];         // compressed matrix  positions
  double** val_row = new double*[3];   // compressed matrix  values

  for (int ql1 = 0; ql1 < 3; ql1++)  // type of variable loop (0=k,1=l,2=q)
    if (nvars_in[ql1] != 0) {        // check if the prolongator is needed for this variable
      std::ostringstream mode;
      mode << ql1;
      // reading diagonal row length and off row length
      length_row[ql1] = new int[dim_qlk[ql1][0] + 1];  // prepare vectors of correct dimension
      length_off_row[ql1] = new int[dim_qlk[ql1][0] + 1];
      Prol.read_len_hdf5(name.c_str(), mode.str().c_str(), length_row[ql1], length_off_row[ql1]);
      // reading position and values of the prologator
      pos_row[ql1] = new int[length_row[ql1][dim_qlk[ql1][0]]];  // prepare vectors of correct dimension
      val_row[ql1] = new double[length_row[ql1][dim_qlk[ql1][0]]];
      Prol.read_pos_hdf5(name.c_str(), mode.str().c_str(), pos_row[ql1], val_row[ql1]);
    } else {  // if there is no variable of ql1 type assign -1 in the length vector
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
  int* off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0];  // quad offset index   (subdomains and levels)
  off_nd[1] = _mgmesh._off_nd[1];  // linear offset index (subdomains and levels)

  // Set up Projector pattern -------------------------------------------
  int off_proc = _NoLevels * _iproc;
  int ml[3];  // row processor dofs
  ml[0] = (_mgmesh._off_el[0][off_proc + Level + 1] - _mgmesh._off_el[0][off_proc + Level]) *
          _el_dof[0];                                             // row constant processor dofs
  ml[1] = off_nd[1][off_proc + Level + 1] - off_nd[1][off_proc];  // row linear processor dofs
  ml[2] = off_nd[0][off_proc + Level + 1] - off_nd[0][off_proc];  // row quadratic processor dofs

  int nl[3];  // column processor dofs
  nl[0] = (_mgmesh._off_el[0][off_proc + Level] - _mgmesh._off_el[0][off_proc + Level - 1]) *
          _el_dof[0];                                         // column constant processor dofs
  nl[1] = off_nd[1][off_proc + Level] - off_nd[1][off_proc];  // column linear processor dofs
  nl[2] = off_nd[0][off_proc + Level] - off_nd[0][off_proc];  // column quadratic processor dofs

  // initial prolongator entry (parallel)
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0;  // proc initial dofs [0]=const [1]=linear [2] quad

  for (int isubd = 0; isubd < _iproc; isubd++) {
    ml_init[0] +=
        (_mgmesh._off_el[0][Level + 1 + isubd * _NoLevels] - _mgmesh._off_el[0][Level + isubd * _NoLevels]) *
        _el_dof[0];
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
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev[3];
  lev[0] = Level - 1;
  lev[1] = _NoLevels;
  lev[2] = Level - 1;

  if (Level > 1) { lev[1] = Level - 2; }

    ///-------------------------------------------------------------------------------------
    /// [3] Write prolongator sparse structure----------------------------------------------
    ///-------------------------------------------------------------------------------------

#ifdef HAVE_LASPACKM  // ------------ only laspack -----------------
  // pattern dimension
  Graph pattern;
  pattern.resize(m);
  pattern._n = n;
  pattern._m = m;  // global
  pattern._ml = m_l;
  pattern._nl = n_l;  // local
  int ml_start = _node_dof[Level][off_nd[0][off_proc]];
  pattern._ml_start = ml_start;

  // pattern structure
  for (int iql = 0; iql < 2; iql++) {  // [0]=quad [1]=linear
    for (int ivar = 0; ivar < nvar[iql]; ivar++) {
      int off_varl = (ivar + iql * nvar[0]) * offset;

      for (int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++) {  // element loop
        int irow = _node_dof[Level][_mgmesh._node_map[Level - iql][i] + off_varl];
        int ncol = length_row[iql][i + 1] - length_row[iql][i];
        // pattern dimensions
        pattern[irow].resize(ncol + 1);
        pattern[irow][ncol] = length_off_row[iql][i + 1] - length_off_row[iql][i];

        // fill the pattern
        for (int j = 0; j < ncol; j++) {
          pattern[irow][j] =
              _node_dof[Level - 1]
                       [_mgmesh._node_map[lev[iql] - 1 + iql][pos_row[iql][j + length_row[iql][i]]] +
                        off_varl];
        }
      }
    }  // end var loop -------------------------------------
  }

  // Update sparsity pattern
  Prol.update_sparsity_pattern(pattern);
  // clean
  pattern.clear();
#endif  // ------------ end only laspack -----------------
#ifdef HAVE_PETSCM
  std::vector<int> n_nz(m_l);  // # diagonal entries
  std::vector<int> n_oz(m_l);  // # offset entries

  for (int ilq = 0; ilq < 3; ilq++) {  // type of variable loop (0=k,1=l,2=q)
    int mll_var = 0;                   // counter for preceding variables

    for (int count = ilq + 1; count < 3; count++) { mll_var += nvars_in[count] * ml[count]; }

    for (int ivar = 0; ivar < nvars_in[ilq]; ivar++) {                        // loop on variables of ilq type
      for (int i = ml_init[ilq]; i < ml_init[ilq] + ml[ilq]; i++) {           // row loop
        int len = (length_row[ilq][i + 1] - length_row[ilq][i]);              // diag row length
        int len_off = (length_off_row[ilq][i + 1] - length_off_row[ilq][i]);  // off row length
        n_oz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len_off;          // set petsc off-row length
        n_nz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len - len_off;    // set petsc diag-row length
      }                                                                       // row loop
    }                                                                         // variables loop of ilq type
  }                                                                           // type of variable loop

  Prol.update_sparsity_pattern(m, n, m_l, n_l, n_oz, n_nz);
#endif

  ///-------------------------------------------------------------------------------------
  /// [4] Assemblying the Prolongator-----------------------------------------------------
  ///-------------------------------------------------------------------------------------

  DenseMatrixM* valmat;
  std::vector<int> tmp(1);

  for (int iql = 0; iql < 3; iql++) {  // type of variable loop (0=k,1=l,2=q)
    int iql1 = 2 - iql;                // to set the correct level (0 for quad, 1 for linear)
    int mll_var = 0;                   // counter for preceding variables

    for (int count = iql + 1; count < 3; count++) { mll_var += nvars_in[count]; }

    for (int ivar = 0; ivar < nvars_in[iql]; ivar++) {  // loop on variables of iql type
      int off_varl = (ivar + mll_var) * offset;

      for (int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++) {  // row loop
        int ind_irow = i;

        if (iql > 0) { ind_irow = _mgmesh._node_map[Level - iql1][i]; }

        int irow = node_dof_f[ind_irow + off_varl];              // dof of the row
        int ncol = length_row[iql][i + 1] - length_row[iql][i];  // row length

        valmat = new DenseMatrixM(1, ncol);  // row storage
        tmp[0] = irow;
        std::vector<int> ind(ncol);  // row-col indeces

        for (int j = 0; j < ncol; j++) {  // fill the row with values and the indices - column loop
          int jcol = j + length_row[iql][i];
          int jpos = pos_row[iql][jcol];  // set the position
          int ind_jpos = jpos;

          if (iql > 0) { ind_jpos = _mgmesh._node_map[lev[iql]][jpos]; }

          ind[j] = node_dof_c[ind_jpos + off_varl];
          //      if (iql==1) std::cout<<"  bc[0]  "<< bc[0][ind[j]]<<"  pos  "<< ind[j]<<"\n";
        }  // column loop

        for (int j = 0; j < ncol; j++) {  // fill the row with values and the indices - column loop
          int jcol = j + length_row[iql][i];
          double tmp = 1.;

          if (iql == 1 && _bc[0][irow] > 1.5 && _bc[0][irow] < 2.5) { tmp = 1.; }

          (*valmat)(0, j) = tmp * val_row[iql][jcol];  // assign the value in the matrix
          int jpos = pos_row[iql][jcol];               // set the position
          int ind_jpos = jpos;

          if (iql > 0) { ind_jpos = _mgmesh._node_map[lev[iql]][jpos]; }

          ind[j] = node_dof_c[ind_jpos + off_varl];
        }  // column loop

        Prol.add_matrix(*valmat, tmp, ind);  // insert row in the Prolongator
        delete valmat;                       // clean the row
      }                                      // row loop
    }                                        // loop on variables of iql type
  }                                          // type of variable loop

  //  clean ----------------------------------------------
  for (int ql = 0; ql < 3; ql++) {
    delete[] pos_row[ql];
    delete[] val_row[ql];
    delete[] length_off_row[ql];
    delete[] length_row[ql];
  }

  delete[] pos_row;
  delete[] val_row;
  delete[] length_off_row;
  delete[] length_row;

#ifdef PRINT_INFO
  std::cout << " ReadProl(DA): end reading Prol " << name.c_str() << std::endl;
#endif
}

// =================================================================
/// This function read the Restriction linear-quad operators from file (with name)
/// and assemble the operator with (nvars_q) quad and (nvars_l) linear variables
/// for each level (Level)
void MGSolDA::ReadRest(
    const int Level,          // Level
    const std::string& name,  // file name (reading from)
    SparseMMatrixM& Rest,     // Restriction Matrix
    const int nvars_in[],     // # const,linear quad variables
    int node_dof_f[],         // dof map fine mesh
    int node_dof_c[],         // dof map coarse
    int /*_node_dof_top*/[]   // dof map top level
) {                           // -----------------------------------------------------------------

#ifdef PRINT_INFO  // ------------------------------------------
  std::cout << " ReadRest(DA): start restrictor reading    " << std::endl;
#endif
  ///-------------------------------------------------------------------------------
  /// [1] Reading restrictor from file hdf5-----------------------------------------
  ///-------------------------------------------------------------------------------

  // reading dimensions
  int ldim[2];
  int dim_qlk[3][2];  // [0][0] piecewise fine, [2][1] quad coarse

  for (int aa = 0; aa < 3; aa++)
    for (int bb = 0; bb < 2; bb++) { dim_qlk[aa][bb] = 0; }

  for (int jql = 0; jql < 3; jql++) {
    std::ostringstream mode;
    mode << jql;
    int status = Rest.read_dim_hdf5(name.c_str(), mode.str().c_str(), ldim);

    if (status == 0) {
      dim_qlk[jql][0] = ldim[0];
      dim_qlk[jql][1] = ldim[1];
    } else {
      std::cout << "  MGSolDA::ReadRest: error reading dimension ";
      abort();
    }
  }

  // setup length (0=k,1=l,2=q)
  int** length_row = new int*[3];      // total length entry sparse struct
  int** length_off_row = new int*[3];  // offdiag entry sparse matrix struct
  int** pos_row = new int*[3];         // compressed matrix  positions
  double** val_row = new double*[3];   // compressed matrix  values

  // reading diagonal row length and off row length (0=qq.1=ql,3=lq,3=ll)
  for (int ql1 = 0; ql1 < 3; ql1++) {  // type of variable loop (0=k,1=l,2=q)
    if (nvars_in[ql1] != 0) {          // check if the restrictor is needed for this variable
      std::ostringstream mode;
      mode << ql1;
      // reading diagonal row length and off row length
      length_row[ql1] = new int[dim_qlk[ql1][0] + 1];  // prepare vectors of correct dimension
      length_off_row[ql1] = new int[dim_qlk[ql1][0] + 1];
      Rest.read_len_hdf5(name.c_str(), mode.str().c_str(), length_row[ql1], length_off_row[ql1]);
      // reading position and values of the restrictor
      pos_row[ql1] = new int[length_row[ql1][dim_qlk[ql1][0]]];  // prepare vectors of correct dimension
      val_row[ql1] = new double[length_row[ql1][dim_qlk[ql1][0]]];
      Rest.read_pos_hdf5(name.c_str(), mode.str().c_str(), pos_row[ql1], val_row[ql1]);
    } else {  // if there is no variable of ql1 type assign -1 in the length vector
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
  int* off_nd[2];
  off_nd[0] = _mgmesh._off_nd[0];  // quad offset index(subdomainsand levels)
  off_nd[1] = _mgmesh._off_nd[1];  // linear offset index (subdomains and levels)

  // Set up Projector pattern -------------------------------------------
  // Matrix dimension
  // local indeces (parallel matrix)
  int off_proc = _NoLevels * _iproc;
  int ml[3];
  int nl[3];  //([2]=quad [1]=linear [0]=const)

  ml[0] = (_mgmesh._off_el[0][off_proc + Level + 1] - _mgmesh._off_el[0][off_proc + Level]) *
          _el_dof[0];                                             // local const m
  ml[1] = off_nd[1][off_proc + Level + 1] - off_nd[1][off_proc];  // local linear m
  ml[2] = off_nd[0][off_proc + Level + 1] - off_nd[0][off_proc];  // local quadratic m

  nl[0] = (_mgmesh._off_el[0][off_proc + Level + 2] - _mgmesh._off_el[0][off_proc + Level + 1]) *
          _el_dof[0];                                             // local const n
  nl[1] = off_nd[1][off_proc + Level + 2] - off_nd[1][off_proc];  // local linear n
  nl[2] = off_nd[0][off_proc + Level + 2] - off_nd[0][off_proc];  // local quadratic n

  // initial matrix entry (parallel)
  int ml_init[3];
  ml_init[0] = 0;
  ml_init[1] = 0;
  ml_init[2] = 0;  // proc initial dofs [0]=const [1]=linear [2] quad

  for (int isubd = 0; isubd < _iproc; isubd++) {
    ml_init[0] +=
        (_mgmesh._off_el[0][Level + 1 + isubd * _NoLevels] - _mgmesh._off_el[0][Level + isubd * _NoLevels]) *
        _el_dof[0];
    ml_init[1] += off_nd[1][Level + 1 + isubd * _NoLevels] - off_nd[1][isubd * _NoLevels];
    ml_init[2] += off_nd[0][Level + 1 + isubd * _NoLevels] - off_nd[0][isubd * _NoLevels];
  }

  // Local restrictor dimension (parallel)
  int m_l = ml[0] * nvars_in[0] + ml[1] * nvars_in[1] + ml[2] * nvars_in[2];  //  local m
  int n_l = nvars_in[0] * nl[0] + nvars_in[1] * nl[1] + nvars_in[2] * nl[2];  //  local n

  // global indices
  int m = nvars_in[0] * dim_qlk[0][0] + nvars_in[1] * dim_qlk[1][0] + nvars_in[2] * dim_qlk[2][0];
  int n = nvars_in[0] * dim_qlk[0][1] + nvars_in[1] * dim_qlk[1][1] + nvars_in[2] * dim_qlk[2][1];

  // pattern operator storage -----------------------------
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  int lev_c = (Level - 1 + _NoLevels + 1) % (_NoLevels + 1);  // coarse linear level
  int lev[3];
  lev[0] = Level;
  lev[1] = lev_c;
  lev[2] = Level;

  ///-------------------------------------------------------------------------------------
  /// [3] Write restrictor sparse structure-----------------------------------------------
  ///-------------------------------------------------------------------------------------

#ifdef HAVE_PETSCM             // ----------- only petsc -----------------------------
  std::vector<int> n_nz(m_l);  // # diagonal entries
  std::vector<int> n_oz(m_l);  // # offset entries

  //   if (nvars_in[0]==0) {
  // linear + quadratic  operator
  for (int ilq = 0; ilq < 3; ilq++) {  // type of variable loop (0=k,1=l,2=q)
    int mll_var = 0;                   // counter for preceding variables

    for (int count = ilq + 1; count < 3; count++) { mll_var += nvars_in[count] * ml[count]; }

    for (int ivar = 0; ivar < nvars_in[ilq]; ivar++) {                        // loop on variables of ilq type
      for (int i = ml_init[ilq]; i < ml_init[ilq] + ml[ilq]; i++) {           // row loop
        int len = (length_row[ilq][i + 1] - length_row[ilq][i]);              // diag row length
        int len_off = (length_off_row[ilq][i + 1] - length_off_row[ilq][i]);  // off row length
        n_oz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len_off;
        n_nz[i - ml_init[ilq] + mll_var + ivar * ml[ilq]] = len - len_off;
      }  // row loop
    }    // variables loop of ilq type
  }      // type of variable loop

  Rest.update_sparsity_pattern(m, n, m_l, n_l, n_oz, n_nz);
#endif  // ----------- end only petsc --------------------------------------

#ifdef HAVE_LASPACKM  // ----------- only laspack  -----------------------------
  // pattern structure dimension
  Graph pattern;
  pattern.resize(m);  // values.resize(nrowt);
  pattern._m = m;
  pattern._n = n;                                        // global dim _m x _n
  pattern._nl = n_l;                                     //  local _n
  pattern._ml = m_l;                                     //  local _m
  int ml_start = _node_dof[Level][off_nd[0][off_proc]];  //  offset proc nodes
  pattern._ml_start = ml_start;                          // starting indices for local matrix

  // quadratic matrix dof (mlq_init-mlq_init+ml_q)
  for (int iql = 0; iql < 2; iql++) {               // [0]=quad [1]=linear
    for (int ivar = 0; ivar < nvar[iql]; ivar++) {  // variable loop
      int off_var_q = (ivar + iql * nvar[0]) * offset;

      for (int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++) {               // element loop
        int ncol = length_row[iql][i + 1] - length_row[iql][i];                   // # columns
        int irow = _node_dof[Level][_mgmesh._node_map[lev[iql]][i] + off_var_q];  // row
        // pattern data
        pattern[irow].resize(ncol + 1);                                             // dimension
        pattern[irow][ncol] = length_off_row[iql][i + 1] - length_off_row[iql][i];  // # entry off diag

        for (int j = 0; j < ncol; j++) {  // inserting  row data
          pattern[irow][j] =
              _node_dof[Level + 1]
                       [_mgmesh._node_map[Level + 1][pos_row[iql][j + length_row[iql][i]]] + off_var_q];
        }  //  row data
      }    // element loop
    }      // variable loop
  }        // [0]=quad [1]=linear

  //  update sparsity pattern for Rst ------------------
  Rst[Level]->update_sparsity_pattern(pattern);
  pattern.clear();
#endif  // ----------- end only laspack  -----------------------------

  ///-------------------------------------------------------------------------------------
  /// [4] Assemblying the Restrictor------------------------------------------------------
  ///-------------------------------------------------------------------------------------

  DenseMatrixM* valmat;
  std::vector<int> tmp(1);

  for (int iql = 0; iql < 3; iql++) {  // type of variable loop (0=k,1=l,2=q)
    int iql1 = 2 - iql;                // to set the correct level (0 for quad, 1 for linear)
    int mll_var = 0;                   // counter for preceding variables

    for (int count = iql + 1; count < 3; count++) { mll_var += nvars_in[count]; }

    for (int ivar = 0; ivar < nvars_in[iql]; ivar++) {  // loop on variables of iql type
      int off_val_l = (ivar + mll_var) * offset;

      for (int i = ml_init[iql]; i < ml_init[iql] + ml[iql]; i++) {  // row loop
        int top_node = i;

        if (iql > 0) { top_node = _mgmesh._node_map[lev[iql]][i]; }

        int irow = node_dof_c[top_node + off_val_l];             // dof of the row
        int ncol = length_row[iql][i + 1] - length_row[iql][i];  // row length - # colunms
        int irow_top = _node_dof[_NoLevels - 1][top_node + off_val_l];
        valmat = new DenseMatrixM(1, ncol);  // row storage
        tmp[0] = irow;
        std::vector<int> ind(ncol);  // row-col indeces

        for (int j = 0; j < ncol; j++) {  // fill the row with values and the indices - column loop
          int jcol = j + length_row[iql][i];
          int jpos = pos_row[iql][jcol];  // set the position
          int ind_jpos = jpos;

          if (iql > 0) { ind_jpos = _mgmesh._node_map[Level + 1 - iql1][jpos]; }

          ind[j] = node_dof_f[ind_jpos + off_val_l];
          (*valmat)(0, j) = (_bc[0][irow_top] % 2) *
                            val_row[iql][j + length_row[iql][i]];  // assign the value in the matrix
        }                                                          // column loop

        Rest.add_matrix(*valmat, tmp, ind);  // insert row in the Restrictor
        delete valmat;                       // clean the row
      }                                      // row loop
    }                                        // loop on variables of iql type
  }                                          // type of variable loop

  //  clean   -----------------------------------
  for (int ql = 0; ql < 3; ql++) {
    delete[] pos_row[ql];
    delete[] val_row[ql];
    delete[] length_off_row[ql];
    delete[] length_row[ql];
  }

  delete[] pos_row;
  delete[] val_row;
  delete[] length_off_row;
  delete[] length_row;

#ifdef PRINT_INFO
  std::cout << " ReadRest(DA): end reading Rest " << name.c_str() << std::endl;
#endif
  return;
}

// ============================================================
void MGSolDA::print_u(
    std::string namefile,
    const int Level) {  // ============================================================
  int length = namefile.length();
  if (namefile.substr(length - 3, 3) == "med") { print_u_med(namefile, Level); }
  if (namefile.substr(length - 2, 2) == "h5") { print_u_hdf5(namefile, Level); }
  return;
}

// ============================================================
/// This function prints the solution for quad fem
void MGSolDA::print_u_hdf5_quad(
    const int Level, const int istep, const int offset, hid_t file_id,
    hsize_t dimsf[]) {  // ===============================================
  double* sol = new double[dimsf[0]];

  //     for ( int ivar = 0; ivar < dimsf[0]; ivar++ ) sol[ivar]  = 0.123;
  // print quad -------------------------------------
  for (int ivar = 0; ivar < _nvars[2]; ivar++) {
    std::string var_name = _var_names[ivar];
    for (int i = 0; i < dimsf[0]; i++) {
      sol[i] = (*x_old[0][Level])(_node_dof[Level][i + ivar * offset]) * _refvalue[ivar];
    }
    _mgutils.print_Dhdf5(file_id, var_name, dimsf, sol);
  }
  delete[] sol;
  //     printf("MGSolDA::print_u_hdf5_quad done %s \n", _var_names[0].c_str());
  return;
}

// ============================================================
/// This function prints the solution for  linear fem
void MGSolDA::print_u_hdf5_lin(
    const int Level, const int istep, const int offset, hid_t file_id,
    hsize_t dimsf[]) {  // ===============================================
                        // //       std::cout<<  " node " << in << " "; ;
                        //     for(int in = 0; in < 9; in++) {
                        //       int idofu = 25*2+in;
                        //  std::cout<<  " node " << in << " "; ;
  // //       for(int ivar = 0; ivar < _nNSdim; ivar++) {  // always oold for direction -> get_oold_sol
  //         double val=    _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F] ]->get_sol(0,idofu);
  //         std::cout<<     val<< " ";
  // //       }// int tn
  //   std::cout<<  " \n " ;
  //
  //     }
  double* sol_c = new double[_mgmesh._GeomEl.n_l[0]];
  double* sol = new double[dimsf[0]];
  for (int ivar = _nvars[2]; ivar < _nvars[2] + _nvars[1]; ivar++) {
    std::string var_name = _var_names[ivar];
    //  2bB element interpolation over the fine mesh ------------------------------------------
    for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
      int start = _mgmesh._off_el[0][_NoLevels * iproc + Level];
      int end = _mgmesh._off_el[0][Level + _NoLevels * iproc + 1];
      for (int iel = start; iel < end; iel++) {  // element iel ----------------------------
        // printf(" %d %d %d %d \n",iproc, iel, start, end);
        int indx = iel * NDOF_FEM;
        for (int in = 0; in < NDOF_P; in++) {  // vertices
          int gl_i = _mgmesh._el_map[0][indx + in];
          double val = (*x_old[istep][Level])(_node_dof[Level][gl_i + ivar * offset]);
          sol_c[in] = val * _refvalue[ivar];
          sol[gl_i] = sol_c[in];
        }
        for (int in = 0; in < NDOF_FEM; in++) {  // mid-points
          double sum = 0;
          for (int jn = 0; jn < NDOF_P; jn++) { sum += _mgmesh._GeomEl.Prol[in * NDOF_P + jn] * sol_c[jn]; }
          int k = _mgmesh._el_map[0][iel * NDOF_FEM + in];
          sol[k] = sum;
        }
      }  // ---- end iel -----------------------------------------------------------------
    }    // 2bB end interpolation over the fine mesh ------------------------

    _mgutils.print_Dhdf5(file_id, var_name, dimsf, sol);
  }  // ivar
  delete[] sol;
  delete[] sol_c;
  //     printf("MGSolDA::print_u_hdf5_lin done %s \n", _var_names[0].c_str());
  return;
}

// ======================================
// This function prints the solution for constant fem
void MGSolDA::print_u_hdf5_const(
    const int Level, const int istep, const int offset, hid_t file_id,
    hsize_t dimsf[]) {  // ===============================================
  // print constant function -----------------------------------
  double* sol_p = new double[dimsf[0]];
  for (int ivar = 0; ivar < _nvars[0]; ivar++) {
    std::string var_name = _var_names[ivar + _nvars[2] + _nvars[1]];
    //  2bB element interpolation over the fine mesh -----------------------
    int eldof_lev = 0;
    for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
      int delta =
          _mgmesh._off_el[0][Level + 1 + _NoLevels * iproc] - _mgmesh._off_el[0][_NoLevels * iproc + Level];
      for (int iel = 0; iel < delta; iel++) {
        int indx = iel + eldof_lev;  // print only the first dof
        double val = (*x_old[0][Level])(
                         _node_dof[Level][indx * _el_dof[0] + (ivar + _nvars[2] + _nvars[1]) * offset]) *
                     _refvalue[_nvars[0]];
        for (int isubcell = 0; isubcell < NSUBDOM; isubcell++) {
          sol_p[indx * NSUBDOM + isubcell] = val * _refvalue[DIMENSION];
        }
      }  // ---- end iel -------------------------
      eldof_lev += delta;
    }  // 2bB end interpolation over the fine mesh ------------------------
    _mgutils.print_Dhdf5(file_id, var_name, dimsf, sol_p);
  }  // ivar
  delete[] sol_p;
  return;
}

// =======================================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_u_hdf5(
    std::string namefile,
    const int Level) {  // ====================================================================
  //  set up
  int SolToPrint = 0;  // print of x_old[SolToPrint] (default=0)
  hid_t file_id = H5Fopen(namefile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  const int n_nodes = _mgmesh._NoNodes[Level];  // count of points
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  hsize_t dimsf[2];
  dimsf[0] = n_nodes;
  dimsf[1] = 1;
  // print quad -----------------------------------
  if (_nvars[2] > 0) print_u_hdf5_quad(Level, SolToPrint, offset, file_id, dimsf);
  // print linear -----------------------------------
  if (_nvars[1] > 0) print_u_hdf5_lin(Level, SolToPrint, offset, file_id, dimsf);

  // count of elements
  if (_nvars[1] > 0) {
    int size = 0;
    for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
      size +=
          _mgmesh._off_el[0][Level + 1 + _NoLevels * iproc] - _mgmesh._off_el[0][_NoLevels * iproc + Level];
    }
    dimsf[0] = size * NSUBDOM;
    dimsf[1] = 1;
    // print constant -----------------------------------
    print_u_hdf5_const(Level, SolToPrint, offset, file_id, dimsf);
  }
  // clean and close -----------------
  H5Fclose(file_id);

  return;
}

// =============================================
/// Print xdmf attrib
void MGSolDA::print_u_xdmf(
    std::ofstream& out,    ///<  xdmf file
    int nodes,             ///<  number of nodes
    int nelems,            ///<  number of elements
    std::string file_name  ///<  xdmf file name
    ) const {              // ================================
  // print xml quad and linear vars over pts
  for (int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++) {
    std::string var_name = _var_names[ivar];
    print_u_xdmf_pts(out, nodes, var_name, file_name);
  }
  // print xml constant vars over elements
  for (int ivar = 0; ivar < _nvars[0]; ivar++) {
    std::string var_name = _var_names[ivar + _nvars[2] + _nvars[1]];
    print_u_xdmf_elem(out, nelems, var_name, file_name);
  }
  return;
}

// =============================================
/// Print xml attrib
void MGSolDA::print_u_xdmf_pts(
    std::ofstream& out,    ///<  xdmf file//  file xdmf
    int nodes,             ///<  number of nodes
    std::string var_name,  ///<  var names
    std::string file_name  ///<  xdmf file name

    ) const {  // ================================

  out << "<Attribute Name=\"" << var_name << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << nodes << "  " << 1
      << "\" Format=\"HDF\">  \n";
  out << file_name << ":" << var_name << "\n";
  out << "</DataItem>\n"
      << "</Attribute>";

  return;
}
// =============================================
/// Print xml attrib
void MGSolDA::print_u_xdmf_elem(
    std::ofstream& out,    ///<  xdmf file//  file xdmf
    int nelems,            ///<  number of elements
    std::string var_name,  ///<  var names
    std::string file_name  ///<  xdmf file name

    ) const {  // ================================
  out << "<Attribute Name=\"" << var_name << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << nelems << "  " << 1
      << "\" Format=\"HDF\">  \n";
  out << file_name << ":" << var_name << "\n";
  out << "</DataItem>\n"
      << "</Attribute>";

  return;
}
  //**********************************************************************************
  //   MED LIBRARY
  //**********************************************************************************
  //**********************************************************************************
  //   MED LIBRARY
  //**********************************************************************************
  //**********************************************************************************
  //   MED LIBRARY
  //**********************************************************************************

#ifdef HAVE_MED
// ============================================================
/// This function prints the controlled domain in med format
void MGSolDA::print_weight_med(
    std::string namefile, const int /* Level1*/
) {                                 // ===============================================

  int n_nodes = _mgmesh._NoNodes[_NoLevels - 1];
  int n_elements = _mgmesh._NoElements[0][_NoLevels - 1];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * (_mgmesh._n_subdom - 1)];
  int nodes_el = _mgmesh._type_FEM[0];

#if ELTYPE == 27
#if DIMENSION == 3
  const unsigned int nodesinv[] = {7, 4,  5,  6,  3,  0,  1,  2,  19, 16, 17, 18, 11, 8,
                                   9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};
  const unsigned int nodesinvbd[] = {3, 0, 1, 2, 7, 4, 5, 6, 8};
#else
  const unsigned int nodesinv[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const unsigned int nodesinvbd[] = {0, 1, 2};
#endif
#endif
#if ELTYPE == 10
#if DIMENSION == 3
  const unsigned int nodesinv[] = {2, 3, 1, 0, 9, 8, 5, 6, 7, 4};
  const unsigned int nodesinvbd[] = {1, 2, 0, 4, 5, 3};
#else
  const unsigned int nodesinv[] = {1, 2, 0, 4, 5, 3};
  const unsigned int nodesinvbd[] = {0, 2, 1};
#endif
#endif

  double* coord;
  coord = new double[n_nodes * 3];

  for (int i = 0; i < n_nodes; i++) {
    coord[i * 3 + 1] = 0.;
    coord[i * 3 + 2] = 0.;

    for (int idim = 0; idim < _mgmesh._dim; idim++) {
      coord[i * 3 + idim] = _mgmesh._xyz[i + idim * n_nodes];
    }
  }

  int Level = _mgmesh._NoLevels - 1;
  int icount = 0;
  int* conn;
  conn = new int[n_elements * nodes_el];

  for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
    for (int el = _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level];
         el < _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level + 1]; el++) {
      for (int i = 0; i < nodes_el; i++) {
        conn[icount] = _mgmesh._el_map[0][el * nodes_el + nodesinv[i]];
        icount++;
      }
    }
  }

  MEDCoupling::MEDCouplingUMesh* mesh = MEDCoupling::MEDCouplingUMesh::New("Mesh_1", _mgmesh._dim);
  mesh->allocateCells(n_elements);

  for (int i = 0; i < n_elements; i++) { mesh->insertNextCell(MED_EL_TYPE, nodes_el, conn + i * nodes_el); }

  mesh->finishInsertingCells();

  MEDCoupling::DataArrayDouble* coordarr = MEDCoupling::DataArrayDouble::New();
  coordarr->alloc(n_nodes, 3);
  std::copy(coord, coord + n_nodes * 3, coordarr->getPointer());
  mesh->setCoords(coordarr);
  //     MEDLoader::WriteUMesh(namefile.c_str(), mesh, true);

  delete[] conn;
  delete[] coord;
  //  set up
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  double* sol = new double[n_nodes + 1];

  MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
  f->setMesh(mesh);
  f->setName("weight");
  MEDCoupling::DataArrayDouble* array = MEDCoupling::DataArrayDouble::New();
  array->alloc(n_elements, 1);
  const int nel_e = _mgmesh._off_el[0][_NoLevels];      // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1];  // stop element
  int i = 0;

  // print quad -------------------------------------
  for (int ivar = 0; ivar < 1; ivar++) {
    std::string var_name = "weight";

    for (int iel = 0; iel < (nel_e - nel_b); iel++) {
      sol[i] = _weight_ctrl[iel + nel_b];
      array->setIJ(i, 0, sol[i]);
      i++;
    }
  }

  f->setArray(array);
  std::string s = "weight_control";
  s += ".med";
  MEDCoupling::WriteField(s.c_str(), f, true);
  coordarr->decrRef();
  mesh->decrRef();
  f->decrRef();
  array->decrRef();
  delete[] sol;

  return;
}

// ============================================================
/// This function prints the solution: for quad and linear fem
void MGSolDA::print_u_med(
    std::string namefile, const int /* Level1*/
) {                                 // ===============================================

  int n_nodes = _mgmesh._NoNodes[_NoLevels - 1];
  int n_elements = _mgmesh._NoElements[0][_NoLevels - 1];
  int nodes_el = _mgmesh._type_FEM[0];

  double* coord;
  coord = new double[n_nodes * 3];

  for (int i = 0; i < n_nodes; i++) {
    coord[i * 3 + 1] = 0.;
    coord[i * 3 + 2] = 0.;

    for (int idim = 0; idim < _mgmesh._dim; idim++) {
      coord[i * 3 + idim] = _mgmesh._xyz[i + idim * n_nodes];
      //        std::cout << " "<< coord[i*3+idim];
    }

    //      std::cout<< " "<< coord[i*3+2] << std::endl;
  }

  std::cout << " " << n_nodes << " " << n_elements << " " << nodes_el << std::endl;
  int Level = _mgmesh._NoLevels - 1;
  int icount = 0;

  int* conn;
  conn = new int[n_elements * nodes_el];

  for (int iproc = 0; iproc < _mgmesh._n_subdom; iproc++) {
    for (int el = _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level];
         el < _mgmesh._off_el[0][iproc * _mgmesh._NoLevels + Level + 1]; el++) {
      for (int i = 0; i < nodes_el; i++) {
        conn[icount] = _mgmesh._el_map[0][el * nodes_el + i];
        icount++;
      }
    }
  }

  MEDCoupling::MEDCouplingUMesh* mesh = MEDCoupling::MEDCouplingUMesh::New("Mesh_1", _mgmesh._dim);
  mesh->allocateCells(n_elements);

  for (int i = 0; i < n_elements; i++) { mesh->insertNextCell(MED_EL_TYPE, nodes_el, conn + i * nodes_el); }

  mesh->finishInsertingCells();

  MEDCoupling::DataArrayDouble* coordarr = MEDCoupling::DataArrayDouble::New();
  coordarr->alloc(n_nodes, 3);
  std::copy(coord, coord + n_nodes * 3, coordarr->getPointer());
  mesh->setCoords(coordarr);
  MEDCoupling::WriteUMesh(namefile.c_str(), mesh, true);
  //   coordarr->decRef();
  //   mesh->decRef();

  delete[] coord;

  //  set up
  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  double* sol = new double[n_nodes + 1];

  MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
  f->setMesh(mesh);
  f->setName(_var_names[0].c_str());
  MEDCoupling::DataArrayDouble* array = MEDCoupling::DataArrayDouble::New();
  array->alloc(n_nodes, 1);

  // print quad -------------------------------------
  for (int ivar = 0; ivar < 1; ivar++) {
    std::string var_name = _var_names[ivar];

    for (int i = 0; i < n_nodes; i++) {
      sol[i] = (*x_old[0][Level])(_node_dof[Level][i + ivar * offset]) * _refvalue[ivar];
      array->setIJ(i, 0, sol[i]);
    }
  }

  f->setArray(array);
  std::string s = "my_test3_sol_1";
  s += ".med";
  MEDCoupling::WriteField(s.c_str(), f, true);

  delete[] conn;
  delete[] sol;
  return;
}
#endif
