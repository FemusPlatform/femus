// std libraries ------------------------------------------
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>

// LibMesh includes ---------------------------------------
#include "enum_elem_type.h"
#include "boundary_mesh.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "elem.h"
#include "mesh_generation.h"
#include "boundary_info.h"
#include "MGFE_conf.h"

// configure includes -------------------------------------
#include "MGGenCase.h"
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"
#include "Domain_conf.h"	//  domain dimensions
#include "MGFE_conf.h"		//


// MED includes -------------------------------------------
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCoupling.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#endif

// local includes -----------------------------------------
#include "MGGenCase_conf.h"
#include "MGUtils.h"
#include "MGGeomEl.h"


// =======================================
/// Constructor
MGGenCase::MGGenCase(const Parallel::Communicator & comm_in,
                     MGUtils & mgutils_in, MGGeomEl & geomel_in):
    _comm(comm_in),
    _geomel(geomel_in),
    _mgutils(mgutils_in) {
    // =================================
    const int NoLevels = (int) _mgutils._geometry["nolevels"];
    _bcmat = 0;
    _dim = DIMENSION;
    _n_levels = NoLevels;
    _N_CHILD = ((DIMENSION > 1) ? 4 * (DIMENSION - 1) : 2);
    _dcl_nel =
        ((DIMENSION >
          1) ? 1 + NDOF_FEM + 3 + 1 + _N_CHILD : 1 + NDOF_FEM + 3 + 1 + _N_CHILD);
    _N_NDV = ((DIMENSION > 1) ? 4 : 2);
    _dclb_nel = NDOF_FEMB + 5;
    _n_subdomains = libMesh::global_n_processors();

    return;
    }

// =======================================================
/// This function is the main mesh generation function
void
MGGenCase::GenCase() {
    // =======================================================

    int dim = DIMENSION;
    std::cout << " Dimension= " << dim << "\n";
#ifdef PRINT_TIME		// ------------------------------------------------
    std::cout << " \n ================================================ \n";
    std::clock_t start_timeA = std::clock();
#endif // -----------------------------------------------------------

    // *************************************************************
    //                   Coarse mesh generator
    // *************************************************************
    // coarse mesh (msh0) generation case 0) from file 1) Internal
    Mesh *msh0 = new Mesh(_comm, _dim);
    // fine mesh (msht) generation case   0) from file 1) Internal
    Mesh *msht = new Mesh(_comm, dim);
    const int libmesh_gen = (int) _mgutils._geometry["libmesh_gen"];	// gen param flag
    switch (libmesh_gen) {
    // mesh generator cases 0) File 1) Internal

    case 0: {
        // mesh from file at level 0 -----------------------
        std::cout << " Reading coarse Mesh File at level 0 \n";
        std::ostringstream mesh_infile;
        mesh_infile << _mgutils.
                    _mesh_dir << _mgutils.get_file("F_MESH_READ");
        std::cout << mesh_infile.str() << std::endl;
        msh0->allow_renumbering(false);
        msh0->read(mesh_infile.str().c_str());
        msht->allow_renumbering(false);
        msht->read(mesh_infile.str().c_str());
        //From HEX8 or HEX20 to HEX27 etc..
        if ((int) _mgutils._geometry["second_order"] == 1) {
            msh0->all_second_order(true);
            msht->all_second_order(true);
            }
        break;
        }				// ----------------------------------------------------------

    case 1: {
        // Internal mesh generator at level 0 ---------------
        std::cout << " Internal mesh generator at level 0 \n";
        libMesh::ElemType libmname;	// libmesh name elem type
        int nintervx = 0;
        int nintervy = 0;
        int nintervz = 0;	// rectangular dimension

        const double LXB = _mgutils._geometry["LXB"];
        const double LXE = _mgutils._geometry["LXE"];
        const double LYB = _mgutils._geometry["LYB"];
        const double LYE = _mgutils._geometry["LYE"];
        const double LZB = _mgutils._geometry["LZB"];
        const double LZE = _mgutils._geometry["LZE"];

        switch (dim) {
        // dimension cases 1,2,3 dim

        case 3: {
            // dimension 3 --------------------------------------
            nintervx = (int) _mgutils._geometry["nintervx"];
            nintervy = (int) _mgutils._geometry["nintervy"];
            nintervz = (int) _mgutils._geometry["nintervz"];
            if (_geomel.name[0] == "Hex_27") {
                libmname = HEX27;
                }
            else if (_geomel.name[0] == "Tet_10") {
                libmname = TET10;
                }
            // else if (_geomel.name[0] == "Wedge_18")  libmname = PRISM18;
            MeshTools::Generation::build_cube(*msh0, nintervx, nintervy, nintervz,	// number of intervals
                                              LXB, LXE, LYB, LYE, LZB, LZE,	// rectangular dimensions
                                              libmname);	// element type
            MeshTools::Generation::build_cube(*msht, nintervx, nintervy, nintervz,	// number of intervals
                                              LXB, LXE, LYB, LYE, LZB, LZE,	// rectangular dimensions
                                              libmname);	// element type
            break;		// end case 3 ----------------------------------------
            }
        case 2: {
            // dimension 2 --------------------------------------
            nintervx = (int) _mgutils._geometry["nintervx"];
            nintervy = (int) _mgutils._geometry["nintervy"];
            if (_geomel.name[0] == "Quad_9") {
                libmname = QUAD9;
                }
            else if (_geomel.name[0] == "Tri_6") {
                libmname = TRI6;
                }
            MeshTools::Generation::build_square(*msh0, nintervx, nintervy, LXB, LXE, LYB, LYE, libmname);	//AX0,AX1,AY0,AY1
            MeshTools::Generation::build_square(*msht, nintervx, nintervy, LXB, LXE, LYB, LYE, libmname);	//AX0,AX1,AY0,AY1
            break;		//  end case 2 ---------------------------------------
            }
        default: {
            // dimension 1 -------------------------------------
            nintervx = (int) _mgutils._geometry["nintervx"];
            libmname = EDGE3;
            MeshTools::Generation::build_square(*msh0, nintervx, 0, LXB,
                                                LXE, 0., 0., libmname);
            MeshTools::Generation::build_square(*msht, nintervx, 0, LXB,
                                                LXE, 0., 0., libmname);
            }			// end default ----------------------------------------------
            }			// switch dim
        break;			// ------------------------------------------------------
        }
    default:
        std::cout << "MGGenCase: Error mesh coarse generation" << std::endl;
        abort();

        }				//--------------------------------------------------------------------------
    std::cout <<
              "\n ========================================================== \n" <<
              "  ================= Coarse Mesh information ============= \n";
    msh0->print_info();		// print the mesh at coarse level

#ifdef PRINT_TIME		// ------------------------------------------------
    std::clock_t end_timeA = std::clock();
    std::cout << " *+* Generation/Reading coarse mesh time ="
              << double(end_timeA - start_timeA) / CLOCKS_PER_SEC << std::endl;
    std::cout << " \n ================================================= \n ";
    std::clock_t start_timeB = std::clock();
#endif // -------------------------------------------------------------

    // ********************************************************
    /// Refine the coarse mesh
    // ********************************************************
    const int NoLevels = (int) _mgutils._geometry["nolevels"];	// param
    // type mesh map -----------------------------------
    int *ttype_FEM;
    ttype_FEM = new int[2];
    ttype_FEM[0] = _geomel.n_q[0];
    ttype_FEM[1] = _geomel.n_q[1];

    // mesh_map_in ----------------------------------
    int *mesh_map_in;
    mesh_map_in = new int[2 * NoLevels];
    for (int itp = 0; itp < 2; itp++) {
        for (int ilev = 0; ilev < NoLevels; ilev++) {
            mesh_map_in[ilev + itp * NoLevels] = ttype_FEM[itp];
            }
        }

    // LibMesh refinement msht[0]-> msht[NoLevels-1]-----------
    if (NoLevels > 1) {
        std::cout << "\n LibMesh Mesh Refinement ---------  \n";
        MeshRefinement mesh_refinement(*msht);
        mesh_refinement.uniformly_refine(NoLevels - 1);
        }

    // *********************************************************
    //  Generating Boundary Mesh from top level
    // *********************************************************
    std::cout << " LibMesh BOUNDARY generation --------- \n";
    // boundary of the coarse  mesh
    BoundaryMesh *bd_msh0 =
        new BoundaryMesh(_comm, msh0->mesh_dimension() - 1);
    msh0->boundary_info->sync(*bd_msh0);
    BoundaryMesh *bd_msht =
        new BoundaryMesh(_comm, msht->mesh_dimension() - 1);
    msht->boundary_info->sync(*bd_msht);

#ifdef PRINT_TIME		// ------------------------------------------
    std::clock_t end_timeB = std::clock();
    std::cout << " *+* Generation refined mesh time ="
              << double(end_timeB - start_timeB) / CLOCKS_PER_SEC << std::endl;
    std::cout << " \n ==================================================  \n ";
    std::clock_t start_timeC = std::clock();
#endif // -------------------------------------------------------

    std::cout << "\n =============== Fine Mesh info ==================== \n";
    msht->print_info();

    // ********************************************************
    ///              PrintMesh
    // ********************************************************
    if (libMesh::global_processor_id() == 0) {
        printMesh(*bd_msht, *msht, *bd_msh0, *msh0, mesh_map_in);
        }

#ifdef PRINT_TIME		// ---------------------------------------------------------
    std::clock_t end_timeC = std::clock();
    std::cout << " Print and Operators time =" << double(end_timeC -
              start_timeC) /
              CLOCKS_PER_SEC << std::endl;
#endif // -----------------------------------------------------------------------

    // ********************************************************
    /// clean and stop
    // ********************************************************
    delete[]ttype_FEM;
    delete[]mesh_map_in;
    delete bd_msht;
    delete msht;
    delete bd_msh0;
    delete msh0;

#ifdef PRINT_TIME		// --------------------------------------------------------------------------------------------
    std::cout << " \n ====================================   \n";
    std::cout << " +*+*+* Total time =" << double(end_timeC -
              start_timeA) /
              CLOCKS_PER_SEC << std::endl;
#endif // --------------------------------------------------------------------------------------------

    std::cout << "\n +++++++++ End GenCase ++++++++++++++++ " << std::endl;
    return;
    }

// ==============================================
/// This function print the mesh
// ----------------------------------------------
void
MGGenCase::printMesh(BoundaryMesh & bd_msht,	// boundary mesh
                     Mesh & msht,	// Volume mesh -> refined one
                     BoundaryMesh & bd_msh0,	// boundary mesh
                     Mesh & msh0,	// coarse mesh
                     const int map_mesh_in[]	// mesh map
                    ) {
    // ============================================

    // set up
    int n_nodes = msht.n_nodes();	//from mesh
    int n_elements = msht.n_elem();	//from mesh
    int n_elements_b = bd_msht.n_elem() * _n_levels;	//from mesh

    Mesh::const_element_iterator it_t00 = msh0.elements_begin();
    const Mesh::const_element_iterator end_t00 = msh0.elements_end();
    Mesh::const_element_iterator it_tr = msht.elements_begin();
    const Mesh::const_element_iterator end_tr = msht.elements_end();
// // //

// // // // active elemnt iterator
//   Mesh::const_element_iterator        it_t00 = msh0.active_elements_begin();
//   const Mesh::const_element_iterator  end_t00 = msh0.active_elements_end();
//   Mesh::const_element_iterator           it_tr = msht.active_elements_begin();
//   const Mesh::const_element_iterator   end_tr = msht.active_elements_end();


    const int n_faces = (*it_t00)->n_sides();
//   const int  n_nodes_lev0=msh0.n_nodes();
    int n_variables = 1;
    // counting
    int *n_elements_lev = new int[2 * _n_levels];
    int *n_nodes_lev = new int[2 * _n_levels + 1];
    int *bd_n_elements_sl = new int[_n_subdomains * _n_levels + 1];
    int *n_nodes_slp = new int[_n_subdomains * _n_levels + 1];
    int *n_elements_sl = new int[_n_subdomains * _n_levels + 1];
    int *n_nodes_sl = new int[_n_subdomains * _n_levels + 1];
    int *elxnode = new int[n_nodes];

    // structures with original mesh -------
    double *nod_val;
    nod_val = new double[n_nodes * 3];
    int *nod_flag;
    nod_flag = new int[n_nodes];
    int *mat_flag;
    mat_flag = new int[n_elements];
    int *map;
    map = new int[n_elements * NDOF_FEM];
    int **elem_sto;
    elem_sto = new int *[n_elements];
    int **elem_conn;
    elem_conn = new int *[n_elements];
    int **bd_elem_sto;
    bd_elem_sto = new int *[n_elements_b];
    int **nod_sto;
    nod_sto = new int *[n_nodes];
    // order vectors -------------------
    std::vector < std::pair < int, int >>v_el(n_elements);
    int *v_inv_el = new int[n_elements];
    std::vector < std::pair < int, int >>v(n_nodes);
    int *v_inv_nd = new int[n_nodes];
    //  local bc and material for each element
    int bc_tmp[NDOF_FEM];
    int maptmp[NDOF_FEM];
    // zeroes ------------------
    for (int i = 0; i < n_nodes; i++) {
        nod_flag[i] = 0;		// bc condition
        }
    for (int i = 0; i < n_nodes * 3; i++) {
        nod_val[i] = 0.;		// nodes
        }
    for (int i = 0; i < n_elements * NDOF_FEM; i++) {
        map[i] = -1;		// nodes
        }
    for (int i = 0; i < _n_levels * _n_subdomains; i++) {
        n_elements_sl[i] = 0;
        }
    for (int i = 0; i < _n_levels * _n_subdomains; i++) {
        n_nodes_sl[i] = 0;
        }
    for (int i = 0; i < _n_levels * _n_subdomains; i++) {
        bd_n_elements_sl[i] = 0;
        }
    for (int ilev = 0; ilev < _n_levels * 2; ilev++) {
        n_elements_lev[ilev] = 0;	// nodes
        }
    for (int ilev = 0; ilev < _n_levels * 2 + 1; ilev++) {
        n_nodes_lev[ilev] = 0;	// nodes
        }
    for (int is = 0; is <= _n_subdomains * _n_levels; is++) {
        n_nodes_slp[is] = 0;
        }
    for (int i = 0; i < n_elements; i++) {
        // elements
        elem_sto[i] = new int[_dcl_nel];
        mat_flag[i] = 1;		// material flag (default->1)
        for (int k = 0; k < _dcl_nel; k++) {
            elem_sto[i][k] = -1;
            }
        elem_conn[i] = new int[n_faces + 1];
        elem_conn[i][0] = n_faces;
        for (int k = 1; k < n_faces + 1; k++) {
            elem_conn[i][k] = -1;
            }
        }
    for (int i = 0; i < n_elements_b; i++) {
        // bd elements
        bd_elem_sto[i] = new int[_dclb_nel];
        for (int j = 0; j < _dclb_nel; j++) {
            bd_elem_sto[i][j] = 0;
            }
        }
//   ================================================
//    ELEMENTS (in elem_sto and bd_elem_sto)
// ---------------------------------------------
//    elem_sto[Id(0), nodes (NDOF_FEM),
//           lev(NDOF_FEM+1),pr(NDOF_FEM+2),parent (NDOF_FEM+3),
//           Nchildren(NDOF_FEM+4), CHILDREN]
// --------------------------------------------------
//    bd_elem_sto[Id(0), elem_id (1),
//           lev(2),side(3), n_nodes (4),
//           side-nodes]
// ================================================
    // Storing the mesh information
    int count_e = 0;
    int count_eb = 0;
    int n_nodes_el;

//   Mesh::const_element_iterator       it_tr = msht->elements_begin();
//   const Mesh::const_element_iterator end_tr= msht->elements_end();
    for (; it_tr != end_tr; ++it_tr) {
        Elem *elem = *it_tr;
        // element id
        int id_el = elem->id();
        elem_sto[count_e][0] = id_el;

        // element nodes
        n_nodes_el = elem->n_nodes();
//     int a= NDOF_FEM;
        for (int inode = 0; inode < n_nodes_el; inode++) {
            int knode = elem->node(inode);
            elem_sto[count_e][1 + inode] = knode;
            // coordinates storage
            for (int idim = 0; idim < DIMENSION; idim++) {
                double xyz = msht.point(knode)(idim);
                nod_val[knode + idim * n_nodes] = xyz;
                }
            }
        // level and subdomain flag
        int lev = elem->level();
        elem_sto[count_e][NDOF_FEM + 1] = lev;
        int subdom_id = elem->processor_id();
        elem_sto[count_e][NDOF_FEM + 2] = subdom_id;
        // parent flag
        Elem *parent = elem->parent();
        if (parent != NULL) {
            elem_sto[count_e][NDOF_FEM + 3] = parent->id();
            }
        // children
        int n_childs = elem->n_children();
        if (elem->has_children()) {
            elem_sto[count_e][NDOF_FEM + 4] = n_childs;
            // children in the element
            for (int i_ch = 0; i_ch < n_childs; i_ch++) {
                Elem *child = (*it_tr)->child(i_ch);
                elem_sto[count_e][NDOF_FEM + 5 + i_ch] = child->id();
                }
            }

        //  boundary ------------------------------------------
        for (int s = 0; s < (int) elem->n_sides(); s++) {

            if (elem->neighbor(s) == NULL) {
                bd_elem_sto[count_eb][0] = count_eb;
                bd_elem_sto[count_eb][NDOF_FEMB + 1] = elem->id();
                bd_elem_sto[count_eb][NDOF_FEMB + 2] = lev;
                n_elements_lev[lev + _n_levels]++;
                bd_elem_sto[count_eb][NDOF_FEMB + 3] = s;
                UniquePtr < Elem > side(elem->build_side(s));
                bd_elem_sto[count_eb][NDOF_FEMB + 4] = (int) side->n_nodes();
                for (int ns = 0; ns < (int) side->n_nodes(); ns++) {
                    bd_elem_sto[count_eb][1 + ns] = side->node(ns);
                    }
                count_eb++;
                }
            else {
                elem_conn[count_e][s + 1] = (elem->neighbor(s))->id();
                }
            }
        count_e++;
        }

//    ================================================================

    int bd_n_elements = count_eb;
    assert(n_elements == count_e);
    int count_lev0 = 0;
    // ================================================
    // ELEMENT GROUPING  (lev(0)>lev(n))
    // ================================================
    // mesh0 sets subdomain (prodcessor id)
//   Mesh::const_element_iterator       it_t00 = msh0->elements_begin();
//   const Mesh::const_element_iterator end_t00= msh0->elements_end();

    for (; it_t00 != end_t00; ++it_t00) {
        // mesh00 processors
        Elem *elem = *it_t00;	// element pt
        int el_id = elem->id();	// el identity
        elem_sto[el_id][NDOF_FEM + 2] = elem->processor_id();	// processor
        elem_sto[el_id][NDOF_FEM + 1] = 0;	// father (coarse=0)
        count_lev0++;		// count
        }

    int n_groups_names = 1;
    int *group_id_names = new int[1];
#ifdef MATBC_INTERFACE

    const int libmesh_gen = (int) _mgutils._geometry["libmesh_gen"];	// gen param flag
    if (libmesh_gen == 0) {
        // ****************************************************************************
        // Boundary condition and material from file
        // file with mat bc info ------------------------------------------------------
        hid_t status = 0;		// read error flag
        hsize_t dims10[2];	// dimension vector for dset (two index)
        // file with bc and mat cond --------------------------------------------------
        std::ostringstream order_bc;	// file name to read
        order_bc << _mgutils._inout_dir << _mgutils.get_file("BASEBC");
        hid_t file_id =
            H5Fopen(order_bc.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

        // Reading bc_flag ----------------------------------------------------------
        hid_t dtset = H5Dopen(file_id, "/bc_flag"
#if HDF5_VERSIONM != 1808
                              , H5P_DEFAULT
#endif
                             );			// get dset
        hid_t filespace = H5Dget_space(dtset);	// get filespace
        status = H5Sget_simple_extent_dims(filespace, dims10, NULL);	// get dim
        int n_nodes_lev00 = dims10[0];	// level zero original element -> 00
        int *bc_condition_lev0 = new int[dims10[0]];	// coarse bc
        status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, bc_condition_lev0);	// read
        assert(status == 0);	// check

// Reading mat_condition ----------------------------------------------------
        dtset = H5Dopen(file_id, "/material"
#if HDF5_VERSIONM != 1808
                        , H5P_DEFAULT
#endif
                       );			// get dset
        filespace = H5Dget_space(dtset);	// get filespace
        status = H5Sget_simple_extent_dims(filespace, dims10, NULL);	// get dim
        int n_elements_lev0 = dims10[0];	// number of elem c level
        int *mat_condition_lev0 = new int[n_elements_lev0];	// coarse mat c
        status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat_condition_lev0);	// read
        assert(status == 0);	// check

//   int n_groups_names=1;
//   int * group_id_names =new int[1];

#ifdef HAVE_GROUP
        // Reading group name  ----------------------------------------------------------
        dtset = H5Dopen(file_id, "/group_names"
#if HDF5_VERSIONM != 1808
                        , H5P_DEFAULT
#endif
                       );			// get dset
        filespace = H5Dget_space(dtset);	// get filespace
        status = H5Sget_simple_extent_dims(filespace, dims10, NULL);	// get dim
        n_groups_names = dims10[0];	// level zero original element -> 00
        delete[]group_id_names;
        group_id_names = new int[dims10[0]];	// coarse bc
        status = H5Dread(dtset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_id_names);	// read
        assert(status == 0);
#endif // end HAVE_GROU  ----------------------------------------------------------
        // creating map file-mesh  to libmesh-mesh ----------------------------------
        for (int i = 0; i < n_nodes_lev00; i++) {
            nod_flag[i] = bc_condition_lev0[i];
            }
        const int order_quad = (int) _mgutils._geometry["second_order"];	// param
        for (int ielem = 0; ielem < n_elements_lev0; ielem++) {
            mat_flag[ielem] = mat_condition_lev0[ielem];
            if (order_quad != 0) {
                for (int inode = 0; inode < NDOF_FEM; inode++) {
                    maptmp[inode] = elem_sto[ielem][1 + inode];
                    }

                //writing bc conditions in the elements (min wins)
                for (int i = 0; i < NDOF_FEM; i++) {
                    // loop over all child nodes
                    int min = 100000000;
                    int count = 0;
                    for (int j = 0; j < NDOF_P; j++) {
                        // loop over father nodes
                        int bc_tmp = nod_flag[maptmp[j]];

                        if (MGGeomEl::Prol[j + i * NDOF_P] != 0) {
                            // nodes with proj
                            if (min > bc_tmp) {
                                min = bc_tmp;	// min wins
                                }
                            count++;
                            }
                        }
                    if (count < 5) {
                        nod_flag[maptmp[i]] = min;	// store the value
                        }
                    }
// +++++++++++++++++++++++++++++

                }			// ibc_flag
            }
//     // bc and material *****************************************


        //clean
        delete[]bc_condition_lev0;
        delete[]mat_condition_lev0;
        }
#endif

    // ****************************************************************************
    //  From coarse level to level n (ordering and filling)
    // subdomain ordering from level 0
    for (int ilev = 0; ilev < _n_levels - 1; ilev++) {
        // level loop
        for (int ielem = 0; ielem < n_elements; ielem++) {
            if (elem_sto[ielem][NDOF_FEM + 1] == ilev) {
                // only el with same level
                for (int ich = 0; ich < elem_sto[ielem][NDOF_FEM + 4]; ich++) {
                    // child loop

                    elem_sto[elem_sto[ielem][NDOF_FEM + 5 + ich]][NDOF_FEM + 2] = elem_sto[ielem][NDOF_FEM + 2];	// writing into my child the rank

#ifdef MATBC_INTERFACE
                    // MATBC_INTERFACE ++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // child material flag --------------------
                    mat_flag[elem_sto[ielem][NDOF_FEM + 5 + ich]] =
                        mat_flag[ielem];

                    // child node flag -----------------------------
                    // storing  the  NDOF_FEM  nod_flag nodes of the child (ich)
                    for (int inode = 0; inode < NDOF_FEM; inode++) {
                        int icnode =
                            elem_sto[elem_sto[ielem][NDOF_FEM + 5 + ich]][1 +
                                    inode];
                        bc_tmp[inode] = nod_flag[icnode];	// child node flag
                        maptmp[inode] = icnode;	// child nodes
                        }

                    //writing bc conditions in the elements (min wins)
                    for (int i = 0; i < NDOF_FEM; i++) {
                        // loop over all child nodes
                        int min = 100000000;
                        for (int j = 0; j < NDOF_P; j++) {
                            // loop over father nodes
                            if (MGGeomEl::Prol[j + i * NDOF_P] != 0) {
                                // nodes with proj
                                if (min > bc_tmp[j]) {
                                    min = bc_tmp[j];	// min wins
                                    }
                                }
                            }
                        nod_flag[maptmp[i]] = min;	// store the value in the child node
                        }
#endif
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    }		// end for ich
                }			//  ==ilev
            }			// end for ielem
        }				// end for (ilev)
    // *************************************************************

    // ================================================
    // REORDERING ELEMENTS (lev>subdomain)
    // ================================================

    // v_el[i].first  : level at which cell i belongs
    // v_el[i].second : cell id, with global numbering
    // level 0: nl0 cells, level 1: nl1 cells, level 2: nl2 cells
    // then global cells numbering is from 0 to (nl0 + nl1 + nl2)

    std::vector < int >ElementsPerLevel(_n_levels);
    for (int k = 0; k < _n_levels; k++) {
        ElementsPerLevel[k] = 0;
        }
    for (int kel = 0; kel < n_elements; kel++) {
        // counting (level and subdom)
        n_elements_lev[elem_sto[kel][NDOF_FEM + 1]]++;
        n_elements_sl[elem_sto[kel][NDOF_FEM + 1] +
                      elem_sto[kel][NDOF_FEM + 2] * _n_levels]++;
        //setup ordering vector
        v_el[kel].first =
            elem_sto[kel][NDOF_FEM + 1] + _n_levels * elem_sto[kel][NDOF_FEM + 2];
        v_el[kel].second = kel;
        ElementsPerLevel[elem_sto[kel][NDOF_FEM + 1]]++;
        }

    for (int k = 0; k < _n_levels; k++) {
        std::cout << "Elements in level " << k << ": " << ElementsPerLevel[k] <<
                  std::endl;
        }

    // std reordering on pairs NDOF_FEM
    std::sort(v_el.begin(), v_el.end());
    for (int i = 0; i < n_elements; i++) {
        int kel = v_el[i].second;
        v_inv_el[kel] = i;
        if (kel == 344 || i == 344)
            std::cout << "kel " << kel << " i " << i << "   " << v_el[i].first <<
                      std::endl;
        }
    // boundary ---------------------------------
    std::vector < std::pair < int, int >>v_elb(bd_n_elements);
    int *v_inv_elb = new int[bd_n_elements];
    for (int kel = 0; kel < bd_n_elements; kel++) {
        int isubdom = elem_sto[bd_elem_sto[kel][NDOF_FEMB + 1]][NDOF_FEM + 2];
        int ilev = bd_elem_sto[kel][NDOF_FEMB + 2];

        bd_n_elements_sl[ilev + isubdom * _n_levels]++;
        v_elb[kel].first = ilev + _n_levels * isubdom;
        v_elb[kel].second = kel;
        }
    // std reordering on pairs
    std::sort(v_elb.begin(), v_elb.end());
    //  print
    for (int i = 0; i < bd_n_elements; i++) {
        int kel = v_elb[i].second;
        v_inv_elb[i] = kel;
        }

// ================================================
//    NODES (in nd_sto and nod_val)
// ---------------------------------------------
//    nd_sto[Index(0), pr(1),lev (2),levp(3),var(4)]
// --------------------------------------------------
//    nod_val[x,y,z]
// ================================================

    // Set up -------------------------------------
    for (int i = 0; i < n_nodes; i++) {
        nod_sto[i] = new int[_N_NDV + 2];
        for (int k = 0; k < _N_NDV + 2; k++) {
            nod_sto[i][k] = -1;
            }
        nod_sto[i][2] = 100000;
        nod_sto[i][3] = _n_levels;
        }

    // Storing the node information--------------
    for (int i = 0; i < n_elements; i++) {
        for (int k = 0; k < n_nodes_el; k++) {

            int knode = elem_sto[i][k + 1];
            nod_sto[knode][0] = knode;
            if (elem_sto[i][NDOF_FEM + 2] > nod_sto[knode][1]) {
                nod_sto[knode][1] = elem_sto[i][NDOF_FEM + 2];	// subdomain
                }
            if (elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][2]) {
                nod_sto[knode][2] = elem_sto[i][NDOF_FEM + 1];	// level
                }
            if (k < NDOF_P && elem_sto[i][NDOF_FEM + 1] < nod_sto[knode][3]) {
                nod_sto[knode][3] = elem_sto[i][NDOF_FEM + 1];
                }
            nod_sto[knode][4] = 1;
            }
        }

    // -----------------------------------------------
    // REORDERING NODES (subdom>lev>levP>var)
    int n_nodes_0p = 0;
    for (int knode = 0; knode < n_nodes; knode++) {
        n_nodes_sl[_n_levels * nod_sto[knode][1] + nod_sto[knode][2]]++;
        if (nod_sto[knode][3] < _n_levels) {
            n_nodes_slp[_n_levels * nod_sto[knode][1] + nod_sto[knode][3]]++;
            }
        if (nod_sto[knode][3] == 0) {
            n_nodes_0p++;
            }
        for (int klev = nod_sto[knode][2]; klev < _n_levels; klev++) {
            n_nodes_lev[klev]++;
            }
        v[knode].first = nod_sto[knode][4] + n_variables *
                         (nod_sto[knode][3] + (_n_levels + 1) * (nod_sto[knode][2] +
                                 _n_levels *
                                 nod_sto[knode][1]));
        v[knode].second = knode;
        }
    n_nodes_lev[_n_levels * 2] = n_nodes_0p;
    sort(v.begin(), v.end());
    //  print
    for (int i = 0; i < (int) v.size(); i++) {
        v_inv_nd[v[i].second] = i;
        }

    for (int ind = 0; ind < n_nodes; ind++) {
        delete[]nod_sto[ind];
        }
    delete[]nod_sto;


// ============================================
// OFFSET
// ============================================
    int *off_nd[2];
    int *off_el = new int[_n_subdomains * _n_levels + 1];
    off_nd[0] = new int[_n_subdomains * _n_levels + 1];
    int *bd_off_el = new int[_n_subdomains * _n_levels + 1];
    off_nd[1] = new int[_n_subdomains * _n_levels + 1];

    off_el[0] = 0;
    int sum = 0;
    off_nd[0][0] = 0;
    int sum_nd = 0;
    bd_off_el[0] = 0;
    int sumb = 0;
    off_nd[1][0] = 0;
    int sum_nd1 = 0;

    for (int isub = 0; isub < _n_subdomains * _n_levels; isub++) {
        sum += n_elements_sl[isub];
        sum_nd += n_nodes_sl[isub];
        sumb += bd_n_elements_sl[isub];
        sum_nd1 += n_nodes_slp[isub];
        off_el[isub + 1] = sum;
        bd_off_el[isub + 1] = sumb;
        off_nd[0][isub + 1] = sum_nd;
        off_nd[1][isub + 1] = sum_nd1;
        }
// max number of element on each node
    for (int inode = 0; inode < n_nodes; inode++) {
        elxnode[inode] = 0;
        }
    for (int i = 0; i < n_elements; i++)	// ALERT
        for (int k = 0; k < n_nodes_el; k++) {
            if (elem_sto[i][1 + NDOF_FEM] == _n_levels - 1) {
                elxnode[elem_sto[i][1 + k]]++;
                }
            }


    int maxelxnode = 0;
    for (int inode = 0; inode < n_nodes; inode++)
        if (elxnode[inode] > maxelxnode) {
            maxelxnode = elxnode[inode];
            }
    std::cout << " Max number of element on each node = " << maxelxnode <<
              " \n";

// clean
    delete[]elxnode;
    delete[]bd_n_elements_sl;
    delete[]n_nodes_slp;
    delete[]n_elements_sl;
    delete[]n_nodes_sl;


    // =================
    // node map -------------------------------------------
    // packaging
    int **g_indexL = new int *[_n_levels + 1];
    for (int ilev = 0; ilev < _n_levels; ilev++) {
        g_indexL[ilev] = new int[n_nodes];
        int count = 0;
        for (int inode = 0; inode < n_nodes; inode++) {
            g_indexL[ilev][inode] = -1;
            }
        for (int iproc = 0; iproc < _n_subdomains; iproc++) {
            for (int inode = off_nd[0][iproc * _n_levels];
                    inode < off_nd[0][iproc * _n_levels + ilev + 1]; inode++) {
                g_indexL[ilev][inode] = count++;
                }
            }
        }
    g_indexL[_n_levels] = new int[n_nodes];
    int count = 0;
    for (int inode = 0; inode < n_nodes; inode++) {
        g_indexL[_n_levels][inode] = -1;
        }
    for (int iproc = 0; iproc < _n_subdomains; iproc++) {
        for (int inode = 0;
                inode <
                off_nd[1][iproc * _n_levels + 1] - off_nd[1][iproc * _n_levels];
                inode++) {
            g_indexL[_n_levels][off_nd[0][iproc * _n_levels] + inode] = count++;
            }
        }

//   print mesh
    print_mesh_h5(n_nodes_lev, map_mesh_in, n_nodes, nod_val,
                  off_nd, n_elements, n_elements_lev, bd_n_elements,
                  off_el, v_inv_nd, v_inv_el, elem_sto, elem_conn, v_el,
                  bd_elem_sto, v_elb, bd_off_el, g_indexL, v, nod_flag,
                  mat_flag);


    for (int iel = 0; iel < n_elements_b; iel++) {
        delete[]bd_elem_sto[iel];
        }
    delete[]bd_elem_sto;
    delete[]nod_val;
    delete[]bd_off_el;

    //data_in/mesh.xmf
    print_multimesh(n_elements_lev, n_nodes_lev);
    // print med

#ifdef HAVE_MED
    const int ibcflag = (int) _mgutils._geometry["ibc_gen"];	/// parameter defined in parameters.in


    if (ibcflag != 0) {
        // name file med-mesh
        std::string mesh_name = _mgutils.get_file("F_MESH_READ");
        unsigned pos = mesh_name.find(".");	// position of "live" in str
        std::string str3 = mesh_name.substr(0, pos);	// get from "live" to the end
        print_med(_n_levels - 1, bd_msht, msht, n_groups_names, group_id_names,
                  v, v_elb, nod_flag, mat_flag, str3.c_str());
        print_MedToMg(_n_levels - 1, bd_msht, msht, n_groups_names,
                      group_id_names, v, v_inv_nd, v_elb, nod_flag, mat_flag,
                      off_el, v_inv_el, v_el, ElementsPerLevel, elem_sto,
                      str3.c_str());
        }
#endif


    delete[]n_elements_lev;	//


// *********************************
//          MGOperators
// *********************************
// #if EL_TYPE==18
// return;
// #else
    const int mgops_gen = (int) _mgutils._geometry["mgops_gen"];

    if (mgops_gen) {
        compute_and_print_MGOps(n_nodes_lev, maxelxnode, off_el, v_el,
                                elem_sto, v_inv_nd, v_inv_el, g_indexL,
                                off_nd);
        }
    // clean
    delete[]v_inv_elb;
    delete[]v_inv_el;
    delete[]v_inv_nd;
    delete[]off_el;
    delete[]off_nd[0];
    for (int iel = 0; iel < n_elements; iel++) {
        delete[]elem_sto[iel];
        }
    delete[]elem_sto;
    for (int ind = 0; ind < _n_levels + 1; ind++) {
        delete[]g_indexL[ind];
        }
    delete[]g_indexL;
    delete[]off_nd[1];
    delete[]n_nodes_lev;
// #endif
#ifdef PRINT_INFO
    std::cout << " MGGenCase::printMesh: Operators  printed \n";
#endif
    return;
    }				//end print_mesh


// ====================================================
/// This function computes and prints MGOps
void
MGGenCase::compute_and_print_MGOps(int *n_nodes_lev,	// nodes for level
                                   int maxelxnode,	// max # of elements for node
                                   int *off_el,	// offset element
                                   std::vector < std::pair < int, int > >v_el,	// element ordering
                                   int **elem_sto,	// element storage
                                   int *v_inv_nd,	// node ordering
                                   int *v_inv_el, int **g_indexL,	// node maps
                                   int *off_nd[]	// node subdomain  offset
                                  ) {


    int max_elnd = maxelxnode;
#ifdef PRINT_INFO
    std::cout << " MGGenCase::compute_and_print_MGOps:  start \n";
#endif


    std::cout << " MGGenCase::compute_and_print_MGOps:  start \n";
    compute_matrix(n_nodes_lev, max_elnd, off_el, v_el, elem_sto, v_inv_nd,
                   g_indexL, off_nd);
#ifdef PRINT_INFO
    std::cout << " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
#endif
    compute_prol(n_nodes_lev, off_el, v_el, elem_sto, v_inv_nd, v_inv_el,
                 g_indexL, off_nd);
#ifdef PRINT_INFO
    std::cout << " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
#endif
    compute_rest(n_nodes_lev, max_elnd, off_el, v_el, elem_sto, v_inv_nd,
                 v_inv_el, g_indexL, off_nd);
#ifdef PRINT_INFO
    std::cout << " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
#endif



// //===========================================
// #ifdef Q2Q0                  //quadratic-piecewise
//     compute_matrix_Q2Q0(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
// #endif
//     compute_prol_Q2Q0(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd, v_inv_el,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
// #endif
//     compute_rest_Q2Q0(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,v_inv_el,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
// #endif
//
// #endif
//
// //==========================================
// #ifdef Q2P1                 //quadratic-3piecewise
//     compute_matrix_Q2P1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
// #endif
//     compute_prol_Q2P1(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd, v_inv_el,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
// #endif
//     compute_rest_Q2P1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,v_inv_el,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
// #endif
// #endif
// //=======================================
//
// #ifdef Q2Q1                 //quadratic-linear
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps:  start \n";
// #endif
//     compute_matrix_Q2Q1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,  g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_matrix end \n";
// #endif
//     compute_prol_Q2Q1(n_nodes_lev,off_el,v_el,elem_sto,v_inv_nd,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_prol end  \n";
// #endif
//     compute_rest_Q2Q1(n_nodes_lev, max_elnd,off_el,v_el,elem_sto,v_inv_nd,g_indexL,off_nd);
// #ifdef PRINT_INFO
//     std::cout<< " MGGenCase::compute_and_print_MGOps: compute_rest  \n";
// #endif
// #endif
// //==========================


//end switch


    return;
    }

// ==============================================================
/// This function prints the mesh structure elem_sto of the MGGenCase class
void
MGGenCase::print_lib_mesh(MGGeomEl & /*Eltype */ , int **elem_sto,
                          int n_elements, int dcl_nel) {

//  const int NDOF_FEM = Eltype.n_q[0]; const int max_elnd = Eltype.n_l[0];
    std::cerr << " N : Id ,  NDOF_P , NDOF_FEM; (lev,pr,parent) ; "
              << " Nchildren: CHILDREN \n ";
    for (int i = 0; i < n_elements; i++) {
        std::cerr << i << " : " << elem_sto[i][0] << ", ";
        for (int k = 1; k < 1 + NDOF_P; k++) {
            std::cerr << elem_sto[i][k] << " ";	// nodes
            }
        std::cerr << ", ";
        for (int k = NDOF_P + 1; k < 1 + NDOF_FEM; k++) {
            std::cerr << elem_sto[i][k] << " ";
            }
        std::cerr << "; ";
        std::cerr << "(" << elem_sto[i][1 + NDOF_FEM] << ",";	//  level
        std::cerr << elem_sto[i][NDOF_FEM + 2] << ",";	//  lproc
        std::cerr << elem_sto[i][NDOF_FEM + 3] << ");  ";	//  parent
        std::cerr << " " << elem_sto[i][NDOF_FEM + 4] << " :  ";	//  n children
        for (int k = 1 + NDOF_FEM + 4; k < dcl_nel; k++) {
            std::cerr << elem_sto[i][k] << " ";
            }
        std::cerr << " \n ";
        }
    return;
    }

// ================================================================
/// This function prints the mesh structure nod_sto of the MGGenCase class
void
MGGenCase::print_lib_node(int **nod_sto, int n_nodes, int /*n_ndv */) {

    std::cerr << " N : Id ,  pr,lev (P); var \n ";
    for (int i = 0; i < n_nodes; i++) {
        std::cerr << i << " : " << nod_sto[i][0] << "; ";	// id
        std::cerr << nod_sto[i][1] << ",";	//  lproc
        std::cerr << nod_sto[i][2] << " (";	//  level
        std::cerr << nod_sto[i][3] << ");  ";	//  level P
        std::cerr << nod_sto[i][4] << " \n ";	//  variable
        }
    return;
    }


void
MGGenCase::print_Mat(int *Mat, int nrow, int ncln) {
    std::cerr << " Matrix \n ";
    for (int i = 0; i < nrow; i++) {
        std::cerr << i << " - ";
        for (int j = 0; j < ncln; j++) {
            std::cerr << Mat[i * ncln + j] << " ";
            }
        std::cerr << "\n ";
        }
    return;
    }

// ======================================================
/// This function print the mesh xdmf format
void
MGGenCase::print_multimesh(const int *n_elements_lev,	// elements for level <-
                           const int *n_nodes_lev	// nodes for level <-
                          ) {
    // ==================================================

    std::string multimesh = _mgutils.get_file("MULTIMESH");
    std::string basemesh = _mgutils.get_file("BASEMESH");
    std::ostringstream inmesh_xmf;
    inmesh_xmf << _mgutils._inout_dir << multimesh << ".xmf";
    std::ofstream out(inmesh_xmf.str().c_str());

    out << "<?xml version=\"1.0\" ?> \n";
    out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
//   out << "[ <!ENTITY HeavyData \"mesh.h5\"> ] ";
    out << "> \n" << " \n";
    out << "<Xdmf> \n";
    out << "<Domain> \n";
    // Volume mesh ----------------------------------------------
    for (int ilev = 0; ilev < _n_levels; ilev++) {
        out << "<Grid Name=\"Mesh" << ilev << "\"> \n";
        out << "<Topology Type=\"" << _geomel.name[0] << "\" Dimensions=\"" <<
            n_elements_lev[ilev] << "\"> \n";
        out << "<DataStructure DataType=\"Int\" Dimensions=\"" <<
            n_elements_lev[ilev] << " " << NDOF_FEM << "\" Format=\"HDF\">  \n";
        out << basemesh << ".h5"	//"mesh.h5"
            << ":/ELEMS/FEM0/MSH" << ilev << "\n";
        out << "</DataStructure> \n";
        out << "</Topology> \n";
        out << "<Geometry Type=\"X_Y_Z\"> \n";
        for (int idim = 0; idim < 3; idim++) {
            out <<
                "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
                << n_nodes_lev[_n_levels - 1] << " 1 \" Format=\"HDF\">  \n";
            out << basemesh << ".h5"	//"mesh.h5"
                << ":/NODES/COORD/X" << idim + 1 << " \n";
            out << "</DataStructure> \n";
            }
        out << " </Geometry>\n";
        // #ifdef MATBC_INTERFACE
        out <<
            " <Attribute Name=\"Material\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        out << "<DataItem DataType=\"Int\" Precision=\"8\" Dimensions=\"" <<
            n_elements_lev[ilev] << " 1 \" Format=\"HDF\">  \n";
        out << basemesh << ".h5" << ":/ELEMS/SUB/MAT" << ilev << "\n";
        out << "</DataItem> \n";
        out << " </Attribute> \n";
        // end MATBC_INTERFACE
// #endif
        out << "</Grid> \n";
        }				// ---------- end volume mesh -----------------------------
    // Boundary mesh ----------------------------------------------
    for (int ilev = 0; ilev < _n_levels; ilev++) {
        out << "<Grid Name=\"Meshb" << ilev << "\"> \n";
        out << "<Topology Type=\"" << _geomel.name[1] << "\" Dimensions=\"" <<
            n_elements_lev[_n_levels + ilev] << "\"> \n";
        out << "<DataStructure DataType=\"Int\" Dimensions=\"" <<
            n_elements_lev[_n_levels +
                           ilev] << " " << NDOF_FEMB << "\" Format=\"HDF\">  \n";
        out << basemesh << ".h5"	//"mesh.h5"
            << ":/ELEMS/FEM1/MSH" << ilev << "\n";
        out << "</DataStructure> \n";
        out << "</Topology> \n";
        out << "<Geometry Type=\"X_Y_Z\"> \n";
        for (int idim = 0; idim < 3; idim++) {
            out <<
                "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
                << n_nodes_lev[_n_levels - 1] << " 1 \" Format=\"HDF\">  \n";
            out << basemesh << ".h5"	//"mesh.h5"
                << ":/NODES/COORD/X" << idim + 1 << " \n";
            out << "</DataStructure> \n";
            }
        out << " </Geometry>\n";
        // #ifdef MATBC_INTERFACE
        out << " <Attribute Name=\"BC\" > \n";
        out << "<DataStructure DataType=\"Int\" Precision=\"8\" Dimensions=\""
            << n_nodes_lev[_n_levels - 1] << " 1 \" Format=\"HDF\">  \n";
        out << basemesh << ".h5" << ":/NODES/COORD/BC \n";
        out << "</DataStructure> \n";
        out << " </Attribute> \n";
        // end MATBC_INTERFACE
        out << "</Grid> \n";
        }
    // ---------- end boundary mesh -----------------------------

    out << "</Domain> \n";
    out << "</Xdmf> \n";
    out.close();

    return;
    }


// =========================================================
/// This function print the mesh in hdf5 format (hdf5 version  1.8.8)
void
MGGenCase::print_mesh_h5(int *n_nodes_lev,
                         const int *map_mesh_in,
                         int n_nodes,
                         double *nod_val,
                         int *off_nd[],
                         int n_elements, int *n_elements_lev,
                         int bd_n_elements, int *off_el, int *v_inv_nd,
                         int *v_inv_el, int **elem_sto, int **elem_conn,
                         std::vector < std::pair < int, int > >v_el,
                         int **bd_elem_sto, std::vector < std::pair < int,
                         int > >v_elb, int *bd_off_el, int **g_indexL,
                         std::vector < std::pair < int, int > >v,
                         int *nod_flag, int *mat_flag) {
    // ===================================================

    hid_t status = 0;
    // file h5 to print ------------------------------------
    std::ostringstream name;
    std::string basemesh = _mgutils.get_file("BASEMESH");
    std::string appl_dir = _mgutils.get_file("APPL_DIR");

    std::ostringstream inmesh;

    inmesh << _mgutils._inout_dir << basemesh << ".h5";
    hid_t file = H5Fcreate(inmesh.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           H5P_DEFAULT);
    int n_meshes = 2 * _n_levels;

//==================================
// ROOT
// =====================
// _dim,_NoFamFEM,_NoLevels,Subdom ------------------
    // packaging data
    int *tdata;
    tdata = new int[4];
    tdata[0] = DIMENSION;
    tdata[1] = 2;
    tdata[2] = _n_levels;
    tdata[3] = _n_subdomains;
    // DFLS vector (hdf5 sorage)
    hsize_t dimsf[2];
    dimsf[0] = 4;
    dimsf[1] = 1;
    status = _mgutils.print_Ihdf5(file, "DFLS", dimsf, tdata);
    // clean
    assert(status == 0);
    delete[]tdata;
    //  _type_FEM -------------------------------
    //  data packaging
    tdata = new int[n_meshes];
    for (int itp = 0; itp < n_meshes; itp++) {
        tdata[itp] = map_mesh_in[itp];
        }
    // hdf5 sorage
    dimsf[0] = n_meshes;
    dimsf[1] = 1;
    status = _mgutils.print_Ihdf5(file, "FEM", dimsf, tdata);
    // clean
    assert(status == 0);
    delete[]tdata;
// // ===========================================
// //  NODE COORDINATES  (COORD)
// // ===========================================
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    hid_t group_id = H5Gcreate(file, "/NODES", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                               , H5P_DEFAULT, H5P_DEFAULT
#endif
                              );
// // ++++++++++++++++++++++++++++++++++++++++++++++++++
// //  COORDINATES
    hid_t subgroup_id = H5Gcreate(file, "/NODES/COORD", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                  , H5P_DEFAULT, H5P_DEFAULT
#endif
                                 );
    // _xyz ---------------------------------------------------------------------
    // packaging data
    double *xcoord;
    xcoord = new double[n_nodes];
    for (int kc = 0; kc < 3; kc++) {
        for (int inode = 0; inode < n_nodes; inode++) {
            xcoord[inode] = nod_val[v[inode].second + kc * n_nodes];
            }
        // hdf5 print
        dimsf[0] = n_nodes;
        dimsf[1] = 1;
        name.str("");
        name << "/NODES/COORD/X" << kc + 1;
        status = _mgutils.print_Dhdf5(file, name.str(), dimsf, xcoord);
        }
    assert(status == 0);

    //  boundary conditions -----------------------------------------------------
    int *bc_con2;
    bc_con2 = new int[n_nodes];
    for (int inode = 0; inode < n_nodes; inode++) {
        bc_con2[inode] = nod_flag[v[inode].second];
        }

    dimsf[0] = n_nodes;
    dimsf[1] = 1;
    name.str("");
    name << "/NODES/COORD/BC";
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, bc_con2);
    assert(status == 0);
    delete[]bc_con2;

    // clean --------------------------------------------------------------------
    delete[]xcoord;
    H5Gclose(subgroup_id);

    // ++++++++++++++++++++++++++++++++++++++++++++++++++
    //  /NODES/MAP
    subgroup_id = H5Gcreate(file, "/NODES/MAP", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                           );

    // +++++++++++++++++++++++++++++++++++++++++++++++++
    // node offset
    dimsf[0] = _n_subdomains * _n_levels + 1;
    dimsf[1] = 1;
    status = _mgutils.print_Ihdf5(file, "/NODES/MAP/OFF_ND", dimsf, off_nd[0]);
    std::cout << "\n ------ OFF_ND ----- \n ";
    for (int i = 0; i < _n_subdomains; i++) {
        std::cout << " Proc " << i << "\n ";
        for (int j = 0; j < _n_levels; j++) {
            std::cout << " Level Nodes " << j << " from " << off_nd[0][j +
                      i *
                      _n_levels]
                      << " to " << off_nd[0][j + 1 + i * _n_levels] << "\n ";
            std::cout << " Level Cells " << j << " from " << off_el[j +
                      i *
                      _n_levels]
                      << " to " << off_el[j + 1 + i * _n_levels] << "\n ";
            }
        }
    std::cout << " \n ";
    // node offset linear
    status =
        _mgutils.print_Ihdf5(file, "/NODES/MAP/OFF_ND1", dimsf, off_nd[1]);
    dimsf[0] = n_meshes + 1;
    dimsf[1] = 1;
    status =
        _mgutils.print_Ihdf5(file, "/NODES/MAP/NDxLEV", dimsf, n_nodes_lev);

    // node map -------------------------------------------
    // level map (hdf5 print)
    dimsf[0] = n_nodes;
    dimsf[1] = 1;
    for (int ilev = 0; ilev < _n_levels; ilev++) {
        name.str("");
        name << "/NODES/MAP/MAP" << ilev;
        status =
            _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[ilev]);
        }
    // linear coarse level map (hdf5 print)
    name.str("");
    name << "/NODES/MAP/MAP" << _n_levels;
    status =
        _mgutils.print_Ihdf5(file, name.str(), dimsf, g_indexL[_n_levels]);
    H5Gclose(subgroup_id);
    H5Gclose(group_id);


    // ===========================================
    //   /ELEMS (CONNECTIVITY)
    // ===========================================
    group_id = H5Gcreate(file, "/ELEMS", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                         , H5P_DEFAULT, H5P_DEFAULT
#endif
                        );
    // +++++++++++++++++++++++++++++++++++
    //  /ELEMS/FEM1  (volume mesh)
    // +++++++++++++++++++++++++++++++++++
    subgroup_id = H5Gcreate(file, "/ELEMS/FEM0", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                           );
    hid_t subgroup_id2 = H5Gcreate(file, "/ELEMS/SUB", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                                   , H5P_DEFAULT, H5P_DEFAULT
#endif
                                  );
    dimsf[0] = 3;
    dimsf[1] = 1;
    int ndofm[3];
    ndofm[0] = NDOF_FEM;
    ndofm[2] = NDOF_P;
    ndofm[1] = 0;
    status = _mgutils.print_Ihdf5(file, "/ELEMS/FEM0/NDOFEM", dimsf, ndofm);

    // NoElements ------------------------------------
    dimsf[0] = _n_levels;
    dimsf[1] = 1;
    status =
        _mgutils.print_Ihdf5(file, "/ELEMS/FEM0/NExLEV", dimsf,
                             &n_elements_lev[0]);

    // offset
    dimsf[0] = _n_subdomains * _n_levels + 1;
    dimsf[1] = 1;
    status = _mgutils.print_Ihdf5(file, "/ELEMS/FEM0/OFF_EL", dimsf, off_el);

    // Connectivity element-node all levels  packaging data  (volume) ------
    int *tempconn;
    tempconn = new int[n_elements * NDOF_FEM];	// coonectivity all mesh
    int *temp_material = new int[n_elements];
    for (int ielem = 0; ielem < n_elements; ielem++) {
        for (int inode = 0; inode < NDOF_FEM; inode++) {
            tempconn[inode + ielem * NDOF_FEM] = v_inv_nd[elem_sto[v_el[ielem].second][inode + 1]];	// swap for proc and lev
            }
        temp_material[ielem] = mat_flag[v_el[ielem].second];
        }
    dimsf[0] = n_elements * NDOF_FEM;
    dimsf[1] = 1;			// global mesh hdf5 storage
    status = _mgutils.print_Ihdf5(file, "/ELEMS/FEM0/MSH", dimsf, tempconn);

    // material condition on lev (ilev)    ------------------------------------
    dimsf[0] = n_elements;
    dimsf[1] = 1;
    name.str("");
    name << "/ELEMS/SUB/MAT";
    status = _mgutils.print_Ihdf5(file, name.str(), dimsf, temp_material);

    // -----------------------------------------------------------------------------

    // Connectivity element-element (neighbrs) at all levels (volume) -------
    const int n_faces = elem_conn[0][0];
    int *tempconn_el;
    tempconn_el = new int[n_elements * n_faces];
    for (int ielem = 0; ielem < n_elements; ielem++) {
        for (int iside = 0; iside < n_faces; iside++) {
            int el_side = elem_conn[v_el[ielem].second][iside + 1];
            // internal element ->  v_inv_el[el_side]; boundary element -1
            //   std::cout << "ielem new: " << ielem << "        ielem old: "
//       <<  v_el[ielem].second << "        el_side: " << el_side <<  "       v_inv_el:  "
//       <<   v_inv_el[el_side] << std::endl;
            if (el_side != -1) {
                tempconn_el[iside + ielem * n_faces] = v_inv_el[el_side];
                }
            else {
                tempconn_el[iside + ielem * n_faces] = -1;	// boundary element -1
                }
            }
        }
    dimsf[0] = n_elements * n_faces;
    dimsf[1] = 1;			// global mesh hdf5 storage
    status =
        _mgutils.print_Ihdf5(file, "/ELEMS/FEM0/EL_NEIG", dimsf, tempconn_el);
    // --------------------------------------------------------------------------
    // level connectivity (for each level) (start from zero)
#if ELTYPE == 18
    int wedge_18[18] =
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17 };
#endif

    for (int ilev = 0; ilev < _n_levels; ilev++) {
        // mesh level vector for connectivity and material
        int *tempconnf = new int[n_elements_lev[ilev] * NDOF_FEM];
        int *temp_material2 = new int[n_elements_lev[ilev]];
        int ltot = 0;
        for (int iproc = 0; iproc < _n_subdomains; iproc++) {
            for (int iel = off_el[iproc * _n_levels + ilev];
                    iel < off_el[iproc * _n_levels + ilev + 1]; iel++) {
                temp_material2[ltot] = temp_material[iel];
                for (int inode = 0; inode < NDOF_FEM; inode++) {
                    tempconnf[ltot * NDOF_FEM + inode] =
                        tempconn[iel * NDOF_FEM +
#if ELTYPE == 18
                                 wedge_18[inode]
#else
                                 inode
#endif
                                ];
                    }		// end for inode
                ltot++;
                }			// end for iel
            }			// end for iproc
        // mesh connectivity on lev (ilev)  ---------------------------------------
        dimsf[0] = n_elements_lev[ilev] * NDOF_FEM;
        dimsf[1] = 1;
        name.str("");
        name << "/ELEMS/FEM0/MSH" << ilev;
        status = _mgutils.print_Ihdf5(file, name.str(), dimsf, tempconnf);

        // material condition on lev (ilev)    ------------------------------------
        dimsf[0] = n_elements_lev[ilev];
        dimsf[1] = 1;
        name.str("");
        name << "/ELEMS/SUB/MAT" << ilev;
        status =
            _mgutils.print_Ihdf5(file, name.str(), dimsf, temp_material2);



        //clean lev structure
        delete[]temp_material2;
        delete[]tempconnf;

        }				// end for ilev




    // clean
    delete[]tempconn;
    delete[]tempconn_el;
    delete[]temp_material;
    H5Gclose(subgroup_id2);
    H5Gclose(subgroup_id);
    // =+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  CONN/FEM2 --> BOUNDARY MESH
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subgroup_id = H5Gcreate(file, "/ELEMS/FEM1", H5P_DEFAULT
#if HDF5_VERSIONM != 1808
                            , H5P_DEFAULT, H5P_DEFAULT
#endif
                           );
    dimsf[0] = 3;
    dimsf[1] = 1;
    ndofm[0] = NDOF_FEMB;
    ndofm[1] = 0;
    ndofm[2] = NDOF_PB;
    status = _mgutils.print_Ihdf5(file, "/ELEMS/FEM1/NDOFEM", dimsf, ndofm);
    // No elements
    dimsf[0] = _n_levels;
    dimsf[1] = 1;
    status =
        _mgutils.print_Ihdf5(file, "/ELEMS/FEM1/NExLEV", dimsf,
                             &n_elements_lev[_n_levels]);

    // Global Surface mesh at all levels  packaging data  (volume) ----------
    tempconn = new int[bd_n_elements * NDOF_FEMB];
    for (int ielem = 0; ielem < bd_n_elements; ielem++) {
        for (int inode = 0; inode < NDOF_FEMB; inode++) {
            tempconn[inode + ielem * NDOF_FEMB] =
                v_inv_nd[bd_elem_sto[v_elb[ielem].second][inode + 1]];
            }
        }
    // global boundary mesh hdf5 sorage ---------------------
    dimsf[0] = bd_n_elements * NDOF_FEMB;
    dimsf[1] = 1;
    status = _mgutils.print_Ihdf5(file, "/ELEMS/FEM1/MSH", dimsf, tempconn);

    dimsf[0] = _n_subdomains * _n_levels + 1;
    dimsf[1] = 1;
    status =
        _mgutils.print_Ihdf5(file, "/ELEMS/FEM1/OFF_EL", dimsf, bd_off_el);
    //  Level surface mesh connectivity ---------------------
    int tot_el = 0;
    for (int ilev = 0; ilev < _n_levels; ilev++) {
        tot_el += n_elements_lev[ilev + _n_levels];
        int *tempconnf = new int[tot_el * NDOF_FEMB];
        int ltot = 0;
        for (int isubdom = 0; isubdom < _n_subdomains; isubdom++) {
            for (int iel = bd_off_el[ilev + isubdom * _n_levels];
                    iel < bd_off_el[ilev + 1 + isubdom * _n_levels]; iel++) {
                for (int inode = 0; inode < NDOF_FEMB; inode++) {
                    tempconnf[ltot * NDOF_FEMB + inode] =
                        tempconn[iel * NDOF_FEMB + inode];
                    }
                ltot++;
                }
            }
        // partial print -------------------------
        dimsf[0] = ltot * NDOF_FEMB;
        dimsf[1] = 1;
        name.str("");
        name << "/ELEMS/FEM1/MSH" << ilev;
        status = _mgutils.print_Ihdf5(file, name.str(), dimsf, tempconnf);
        //clean
        delete[]tempconnf;
        }
    delete[]tempconn;
    H5Gclose(subgroup_id);
    H5Gclose(group_id);
// file closure
    H5Fclose(file);

    return;
    }

// =============================================================
/// This function prints (hdf5 format) the prolongation matrix structures
void
MGGenCase::print_op_h5(  //hdf5 version  1.8.8
    std::string name,	// file name
    int n_nodes_row,	// # of fine grid nodes
    int n_nodes_cln,	// # of coarse grid nodes
    int count_q,	// # of entries in Prol_q
    int *Res_q,	// compressed pos for prolongation Operator
    double *values_q,	// compressed prolongation Operator
    int *len_q,	// row length
    int *len_qoff,	// offset row len
    int qq	// mode linear-quadratic
) {
    // ====================================================

    std::ostringstream mode;
    if (qq == 0 || qq == 1 || qq == 2) {
        mode << "0" << qq;
        }
    else {
        mode << qq;
        }

    hid_t status = 0;
    // file
    hid_t fileR = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Prolongation dimensions
    std::ostringstream name0;
    name0 << "DIM" << mode.str().c_str();
    hsize_t dimsf[2];
    dimsf[0] = 2;
    dimsf[1] = 1;
    int rowcln[2];
    rowcln[0] = n_nodes_row;
    rowcln[1] = n_nodes_cln;
    status = _mgutils.print_Ihdf5(fileR, name0.str().c_str(), dimsf, rowcln);
    assert(status == 0);
    // node offset
    if (Res_q != NULL) {
        // print matrix position of P
        std::ostringstream name1;
        name1 << "POS" << mode.str().c_str();
        dimsf[0] = count_q;
        status =
            _mgutils.print_Ihdf5(fileR, name1.str().c_str(), dimsf, Res_q);
        }
    if (values_q != NULL) {
        // print values
        std::ostringstream name2;
        name2 << "VAL" << mode.str().c_str();
        status =
            _mgutils.print_Dhdf5(fileR, name2.str().c_str(), dimsf, values_q);
        }
    // print row length of R
    std::ostringstream name3;
    name3 << "LEN" << mode.str().c_str();
    dimsf[0] = n_nodes_row + 1;
    status = _mgutils.print_Ihdf5(fileR, name3.str().c_str(), dimsf, len_q);
    // print row length of R
    std::ostringstream name4;
    name4 << "OFFLEN" << mode.str().c_str();
    status =
        _mgutils.print_Ihdf5(fileR, name4.str().c_str(), dimsf, len_qoff);
    //   clean
    H5Fclose(fileR);

    return;
    }



//-----------------------------------------------------------------
/// This function computes the sparse matrix patern
void
MGGenCase::compute_matrix(int *n_nodes_lev,	// nodes for level
                          int max_elnd,	// max # of elements at each point
                          int *off_el,	// element offset
                          std::vector < std::pair < int, int > >v_el,	// element ordering
                          int **elem_sto,	// element structure
                          int *v_inv_nd,	// node inverse ordering
                          int **g_indexL,	// node map
                          int *off_nd[]	// node offset
                         ) {
    // ====================================================

    std::cout << " Computing matrix :\n";

    // structures
    int n_nodes = n_nodes_lev[_n_levels - 1];	// quadratic nodes
    int n_nodesl = n_nodes_lev[_n_levels - 1];	// linear nodes
    int *len_q = new int[n_nodes + 1];	// row length
    int *len_qoff = new int[n_nodes + 1];	// row offdiagonal length
    int *mem = new int[n_nodes];	// row memorization index
    int *Mat_q = new int[n_nodes * NDOF_FEM * max_elnd];
    int *aux = new int[NDOF_FEM * max_elnd * n_nodes];

    // directory
    std::string input_dir = _mgutils._inout_dir;
    std::string f_matrix = _mgutils.get_file("F_MATRIX");

    int level_row[2];
    int level_cln[2];		// level row x cln
    int nodes_top[3];		// nodes row x cln
    int el_nodes[3];		// nodes for elements
    int *offM[2];			// offnode limits

    for (int Level1 = 0; Level1 < _n_levels; Level1++) {

        std::ostringstream name;
        name << input_dir << f_matrix << Level1 << ".h5";
        hid_t fileM =
            H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT);

        // quad matrix dimensions
        // Set up --------------------------
#ifdef PRINT_TIME		//  TC +++++++++++++++
        std::clock_t start_time = std::clock();
#endif
        // nodes
        n_nodes = n_nodes_lev[Level1];	// quad
        int Level_l = (Level1 + 2 * _n_levels) % (2 * _n_levels + 1);
        n_nodesl = n_nodes_lev[Level_l];	// linear
        // levels
        std::cout << " Level =" << Level1 << " \n ";	// quad level print
        Level_l = (Level1 + _n_levels) % (_n_levels + 1);	// linear=Level-1

        // processor  limit for counting -----------------------------
        //  limits quad-linear
        for (int iql = 0; iql < 2; iql++) {
            offM[iql] = new int[_n_subdomains + 1];
            offM[iql][0] = 0;
            for (int idom = 0; idom < _n_subdomains; idom++) {
                offM[iql][idom + 1] =
                    offM[iql][idom] +
                    (off_nd[iql][idom * _n_levels + Level1 + 1] -
                     off_nd[iql][idom * _n_levels]);
                }
            }

#ifdef PRINT_INFO		// ---------------------------------------
        printf
        (" \n Off processor limits quad and linear for level = %d and proc %d \n ---------------------------  \n",
         Level1, _n_subdomains);
        for (int idom = 0; idom <= _n_subdomains; idom++) {
            printf("  %d %d %d \n", idom, offM[0][idom], offM[1][idom]);
            }
#endif
        offM[0][_n_subdomains] += 1;
        offM[1][_n_subdomains] += 1;	// for interval bound

        level_row[0] = Level1;
        level_cln[0] = Level1;
        level_row[1] = Level_l;
        level_cln[1] = Level_l;


        // element nodes  [0]=quad  [1]=linear [2]=piecewise
        el_nodes[0] = NDOF_FEM;
        el_nodes[1] = NDOF_P;
        el_nodes[2] = NDOF_K;

        // nodes [0]=quad  [1]=linear
        nodes_top[0] = n_nodes_lev[level_row[0]];	// linear
        nodes_top[1] = n_nodes_lev[2 * _n_levels];	// quadratic
        nodes_top[2] = 0;
        for (int kdom = 0; kdom < _n_subdomains; kdom++) {
            // piecewise
            nodes_top[2] +=
                (off_el[kdom * _n_levels + Level1 + 1] -
                 off_el[kdom * _n_levels + Level1]) * NDOF_K;
            }
        if (Level1 > 0) {
            nodes_top[1] = n_nodes_lev[Level1 - 1];
            }



        // Generating matrix operators ----------------------------------------

        // matrices qq-ql-lq-ll  ===============================================
        for (int iql = 0; iql < 2; iql++) {
            for (int jql = 0; jql < 2; jql++) {
                int qlt = iql * 2 + jql;
                int name_label = (2 - iql) * 10 + (2 - jql);
                int count_q =
                    compute_mat_qq(Level1, level_row[iql], level_cln[jql],
                                   el_nodes[iql], el_nodes[jql],
                                   n_nodes, max_elnd, off_el,
                                   v_el, elem_sto, v_inv_nd,
                                   g_indexL, off_nd[iql], offM[qlt % 2], aux,
                                   mem, Mat_q, len_q, len_qoff);
                print_op_h5(name.str(), nodes_top[iql], nodes_top[jql],
                            count_q, Mat_q, NULL, len_q, len_qoff, name_label);
                }
            }

        //  matrix   qk  ====================================================
        int iql = 0;
        int jql = 0;
        int qlt = 20;
        int count_q =
            compute_mat_qk(Level1, level_row[iql], level_cln[jql], el_nodes[iql],
                           el_nodes[2],
                           n_nodes, max_elnd, off_el,
                           v_el, elem_sto, v_inv_nd,
                           g_indexL, off_nd[iql], offM[iql], aux, mem, Mat_q,
                           len_q, len_qoff);
        print_op_h5(name.str(), nodes_top[iql], nodes_top[2] /* *NDOF_K */ ,
                    count_q, Mat_q, NULL, len_q, len_qoff, qlt);



        //  matrix   kq ====================================================
        iql = 0;
        jql = 0;
        qlt = 02;
        count_q =
            compute_mat_kq(Level1, level_row[iql], level_cln[jql], el_nodes[2],
                           el_nodes[jql], n_nodes, max_elnd, off_el, v_el,
                           elem_sto, v_inv_nd, g_indexL, off_nd[iql], offM[iql],
                           aux, mem, Mat_q, len_q, len_qoff);
        print_op_h5(name.str(), nodes_top[2] /* *NDOF_K */ , nodes_top[iql],
                    count_q, Mat_q, NULL, len_q, len_qoff, qlt);

        //  matrix   lk ====================================================
        iql = 1;
        jql = 1;
        qlt = 10;
        count_q =
            compute_mat_lk(Level1, level_row[iql], level_cln[jql], el_nodes[jql],
                           el_nodes[2], n_nodes, max_elnd, off_el, v_el,
                           elem_sto, v_inv_nd, g_indexL, off_nd[1], offM[iql],
                           aux, mem, Mat_q, len_q, len_qoff);
        print_op_h5(name.str(), nodes_top[iql], nodes_top[2] /* *NDOF_K */ ,
                    count_q, Mat_q, NULL, len_q, len_qoff, qlt);


        //  matrix   kl ====================================================
        iql = 1;
        jql = 1;
        qlt = 01;
        count_q =
            compute_mat_kl(Level1, level_row[iql], level_cln[jql], el_nodes[2],
                           el_nodes[jql], n_nodes, max_elnd, off_el, v_el,
                           elem_sto, v_inv_nd, g_indexL, off_nd[iql], offM[iql],
                           aux, mem, Mat_q, len_q, len_qoff);
        print_op_h5(name.str(), nodes_top[2] /* *NDOF_K */ , nodes_top[iql],
                    count_q, Mat_q, NULL, len_q, len_qoff, qlt);

        //  matrix   kk ====================================================
        iql = 0;
        jql = 0;
        qlt = 0;
        count_q =
            compute_mat_kk(Level1, level_row[iql], level_cln[jql], el_nodes[2],
                           el_nodes[2], n_nodes, max_elnd, off_el, v_el,
                           elem_sto, v_inv_nd, g_indexL, off_nd[iql], offM[iql],
                           aux, mem, Mat_q, len_q, len_qoff);
        print_op_h5(name.str(), nodes_top[2] /* *NDOF_K */ ,
                    nodes_top[2] /* *NDOF_K */ , count_q, Mat_q, NULL,
                    len_q, len_qoff, qlt);




#ifdef PRINT_TIME		//  TC +++++++++++++++
        std::clock_t end_time = std::clock();
        std::cout << " quad Matrix compute time =" << double(end_time -
                  start_time) /
                  CLOCKS_PER_SEC << std::endl;
#endif
        delete[]offM[0];
        delete[]offM[1];
        H5Fclose(fileM);
        }

    // cleaning -------------------------------------------------
    delete[]len_q;
    delete[]len_qoff;
    delete[]mem;
    delete[]Mat_q;
    delete[]aux;

    return;
    }


// =====================================================================
/// This function computes the prologator operator
void
MGGenCase::compute_prol(int *n_nodes_lev, int *off_el,
                        std::vector < std::pair < int, int > >v_el,
                        int **elem_sto, int *v_inv_nd, int *v_inv_el,
                        int **g_indexL, int *off_nd[]) {

    // name
    std::string input_dir = _mgutils._inout_dir;
    std::string f_prol = _mgutils.get_file("F_PROL");

    int level_row[3];
    int level_cln[3];		// level row x cln
    int nodes_row[3];
    int nodes_cln[3];		// nodes row x cln
    int el_nodes[3];		// nodes for elements

    for (int Level1 = 1; Level1 < _n_levels; Level1++) {

        // Set up --------------------------
        //  file  -------------------------------------------------------
        std::ostringstream name;
        name << input_dir << f_prol << Level1 - 1 << "_" << Level1 << ".h5";
        // Create file (data_in/Prol*.h5)
        hid_t fileP =
            H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT);

        // set up ----------------------------------------
        // levels [0]=quad  [1]=linear [2]=piecewise
        level_row[0] = Level1;
        level_cln[0] = Level1 - 1;	// quadratic
        level_row[1] = Level1 - 1;
        level_cln[1] = _n_levels;
        if (Level1 > 1) {
            level_cln[1] = Level1 - 2;	// linear
            }
        level_row[2] = Level1;
        level_cln[2] = Level1 - 1;	// piecewise

        // nodes [0]=quad  [1]=linear [2]=piecewise
        nodes_row[0] = n_nodes_lev[level_row[0]];
        nodes_cln[0] = n_nodes_lev[level_cln[0]];	// quadratic
        nodes_cln[1] = n_nodes_lev[2 * _n_levels];
        if (Level1 > 1) {
            nodes_cln[1] = n_nodes_lev[Level1 - 2];	// linear
            }
        nodes_row[1] = n_nodes_lev[level_row[1]];

        int sum_lev = 0;
        int sum_levm1 = 0;
        for (int jpr = 0; jpr < _n_subdomains; jpr++) {
            sum_levm1 +=
                off_el[jpr * _n_levels + Level1] - off_el[jpr * _n_levels +
                        Level1 - 1];
            sum_lev +=
                off_el[jpr * _n_levels + Level1 + 1] - off_el[jpr * _n_levels +
                        Level1 + 0];
            }
        nodes_row[2] = sum_lev;
        nodes_cln[2] = sum_levm1;	// piecewise

        // element nodes  [0]=quad  [1]=linear [2]=piecewise
        el_nodes[0] = NDOF_FEM;
        el_nodes[1] = NDOF_P;
        el_nodes[2] = NDOF_K;

        // matrix structures --------------------------
        int *Prol_q;		// compressed sparse prolongation
        double *values_q;		// compressed prolongation values
        int *len_q;		// # of diagonal row entries
        int *len_qoff;		// # of offdiagonal row entries

        // Prolongation Operator   (quad and linear)
        for (int iql = 0; iql < 3; iql++) {
            //[0]=quad  [1]=linear
            int count_q = 0;
            if (iql < 2) {
                len_q = new int[nodes_row[iql] + 1];
                len_qoff = new int[nodes_row[iql] + 1];
                Prol_q = new int[nodes_row[iql] * NDOF_FEM];
                values_q = new double[nodes_row[iql] * NDOF_FEM];

                // compute quadratic and linear prolongator
                count_q =
                    compute_prol_qq(Level1, level_row[iql], level_cln[iql],
                                    el_nodes[iql], nodes_row[iql], off_el, v_el,
                                    elem_sto, v_inv_nd, g_indexL, off_nd[iql],
                                    values_q, Prol_q, len_q, len_qoff);


                }
            else {
                len_q = new int[nodes_row[iql] + 1];
                len_qoff = new int[nodes_row[iql] + 1];
                Prol_q = new int[nodes_row[iql] /* *NDOF_K */ ];
                values_q = new double[nodes_row[iql] /* *NDOF_K */ ];

                // compute piecewise prolongator
                count_q =
                    compute_prol_kk(Level1, level_row[iql], level_cln[iql],
                                    el_nodes[iql], nodes_row[iql], off_el, v_el,
                                    elem_sto, v_inv_nd, v_inv_el, g_indexL,
                                    off_nd[iql], values_q, Prol_q, len_q,
                                    len_qoff);
                }

            print_op_h5(name.str(), nodes_row[iql], nodes_cln[iql], count_q,
                        Prol_q, values_q, len_q, len_qoff, 2 - iql);
            // clean
            delete[]Prol_q;
            delete[]values_q;
            delete[]len_q;
            delete[]len_qoff;
            }

        H5Fclose(fileP);

        }				//end Level1

    return;
    }


// =====================================================================
/// This function computes the restrictor operator
void
MGGenCase::compute_rest(int *n_nodes_lev,	// node level vector
                        int max_elnd,	// max # of elements at each point
                        int *off_el,	// offset element vectors
                        std::vector < std::pair < int, int > >v_el	// element  ordering
                        , int **elem_sto,	// storage element structure
                        int *v_inv_nd,	// node inverse ordering
                        int *v_inv_el, int **g_indexL,	//  map nodes (level)
                        int *off_nd[]	//  node offset vector
                       ) {
    // =====================================================
    // hdf5 print Restrictor    data_in/Rest*.h5
    std::string input_dir = _mgutils._inout_dir;
    std::string f_rest = _mgutils.get_file("F_REST");

    int level_row[3];
    int level_cln[3];		// level row x cln
    int nodes_row[3];
    int nodes_cln[3];		// nodes row x cln
    int el_nodes[3];		// nodes for elements

    // level loop
    for (int Level1 = 0; Level1 < _n_levels - 1; Level1++) {

        // file name and creation  -----------------
        std::ostringstream name;
        name << input_dir << f_rest << Level1 + 1 << "_" << Level1 << ".h5";
        hid_t fileR =
            H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT);

        // Set up --------------------------
        // levels
        level_row[0] = Level1;
        level_cln[0] = Level1 + 1;	// quadratic
        level_row[1] = _n_levels;
        if (Level1 > 0) {
            level_row[1] = Level1 - 1;
            }
        level_cln[1] = Level1;	// linear
        level_row[2] = Level1;
        level_cln[2] = Level1 + 1;	// piecewise
        // nodes
        nodes_row[0] = n_nodes_lev[Level1];
        nodes_cln[0] = n_nodes_lev[Level1 + 1];	// quadratic
        nodes_row[1] = n_nodes_lev[2 * _n_levels];
        if (Level1 > 0) {
            nodes_row[1] = n_nodes_lev[Level1 - 1];	// linear
            }
        nodes_cln[1] = n_nodes_lev[Level1];
        int sum_lev = 0;
        int sum_levp1 = 0;
        for (int jpr = 0; jpr < _n_subdomains; jpr++) {
            sum_levp1 +=
                off_el[jpr * _n_levels + Level1 + 2] - off_el[jpr * _n_levels +
                        Level1 + 1];
            sum_lev +=
                off_el[jpr * _n_levels + Level1 + 1] - off_el[jpr * _n_levels +
                        Level1];
            }
        nodes_row[2] = sum_lev;
        nodes_cln[2] = sum_levp1;	// piecewise
        // element nodes
        el_nodes[0] = NDOF_FEM;
        el_nodes[1] = NDOF_P;
        el_nodes[2] = NDOF_K;

        // matrix structures --------------------------
        int *Rest_q;		// pos
        double *values_q;		// values
        int *mem;			// mem row
        int *len_q;		// # of diagonal row entries
        int *len_qoff;		// # of offdiagonal row entries

        for (int iql = 0; iql < 3; iql++) {
            // [0]=quad [1]=linear [2]=piecewise
            int count_q = 0;
            if (iql < 2) {
                // set up
                Rest_q = new int[nodes_row[iql] * NDOF_FEM * max_elnd];
                values_q = new double[nodes_row[iql] * NDOF_FEM * max_elnd];
                mem = new int[nodes_row[iql] * NDOF_FEM * max_elnd];
                len_qoff = new int[nodes_row[iql] + 1];
                len_q = new int[nodes_row[iql] + 1];

                //  Restrictor quadratic and linear--------------------------------------
                count_q =
                    compute_res_qq(Level1, level_row[iql], level_cln[iql],
                                   el_nodes[iql], nodes_row[iql], nodes_cln[iql],
                                   max_elnd, off_el, v_el, elem_sto, v_inv_nd,
                                   g_indexL, off_nd[iql], mem, values_q, Rest_q,
                                   len_q, len_qoff);
                }
            else {
                // set up
                max_elnd = elem_sto[0][NDOF_FEM + 4];	//<-pay attention to this instruction, do first quad and linear!!
                Rest_q = new int[nodes_row[iql] * max_elnd];
                values_q = new double[nodes_row[iql] * max_elnd];
                mem = new int[nodes_row[iql] * max_elnd];
                len_qoff = new int[nodes_row[iql] + 1];
                len_q = new int[nodes_row[iql] + 1];

                //  Restrictor constant --------------------------------------
                count_q =
                    compute_res_kk(Level1, level_row[iql], level_cln[iql],
                                   el_nodes[iql], nodes_row[iql], nodes_cln[iql],
                                   max_elnd, off_el, v_el, elem_sto, v_inv_nd,
                                   v_inv_el, g_indexL, off_nd[iql], mem,
                                   values_q, Rest_q, len_q, len_qoff);
                }
            // print ----------------------------------------------------
            print_op_h5(name.str(), nodes_row[iql], nodes_cln[iql], count_q,
                        Rest_q, values_q, len_q, len_qoff, 2 - iql);
            // clean ----------------------------------------------------
            delete[]Rest_q;
            delete[]values_q;
            delete[]mem;		// Restrictor structures
            delete[]len_qoff;
            delete[]len_q;	// compressed row length
            }

        H5Fclose(fileR);

        }				// ----- end Level loop  ------

    return;
    }


// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_qq(int Level1,	// level
                          int Level_row,	// row level
                          int Level_clmn,	// column level
                          int ndof_row_ql,	// row nodes
                          int ndof_clmn_ql,	// column nodes
                          int n_nodes,	// top nodes
                          int max_elnd, int *off_el, std::vector < std::pair < int, int > >v_el, int **elem_sto, int *v_inv_nd, int **g_indexL, int *off_nd, int offM[], int aux[], int mem[],	//number elemets row
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    // zeroing matrix
    for (int im = 0; im < n_nodes; im++) {
        mem[im] = 0;
        }
    for (int im = 0; im < n_nodes * NDOF_FEM * max_elnd; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int op = v_inv_nd[elem_sto[el_lib][1 + kk]];
                int i_row = g_indexL[Level_row][op];
                for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                    int opj = v_inv_nd[elem_sto[el_lib][1 + kj]];
                    int ind_row = i_row * NDOF_FEM * max_elnd;
                    Mat_q[ind_row + mem[i_row] + kj] =
                        g_indexL[Level_clmn][opj];
                    }
                mem[i_row] += NDOF_FEM;
                }
            }
        }

    // compressing zeros and counting -------------------------------------
    // counters
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    int ki = 0;
    // zeroing aux
    for (int i = 0; i < NDOF_FEM * max_elnd * n_nodes; i++) {
        aux[i] = 0;
        }
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                (int)(off_nd[kdom * _n_levels + Level1 + 1] -
                      off_nd[kdom * _n_levels]); iki++) {
            // setting aux
            for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
                int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
                if (knod_q >= 0) {
                    aux[knod_q] = 1;
                    }
                }
            // computing various flags and lengths
            for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
                int knod_q = Mat_q[NDOF_FEM * max_elnd * ki + kj];
                if (knod_q >= 0 && aux[knod_q] == 1) {
                    aux[knod_q] = 0;
                    Mat_q[count_q] = knod_q;
                    count_q++;
                    if (knod_q < (int) offM[kdom]
                            || knod_q >= (int) offM[kdom + 1]) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }

    return count_q;
    }





// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_kq(int Level1,	// level
                          int /*Level_row */ ,	// row level
                          int Level_clmn,	// column level
                          int ndof_row_ql,	// row nodes
                          int /*ndof_clmn_ql */ ,	// column nodes
                          int /*n_nodes */ ,	// top nodes
                          int /*max_elnd */ ,
                          int *off_el,
                          std::vector < std::pair < int, int > >v_el,
                          int **elem_sto, int *v_inv_nd,
                          int **g_indexL, int * /*off_nd */ ,
                          int offM[],
                          int aux[],
                          int mem[],
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    int n_element_lev = 0;
    for (int pr = 0; pr < _n_subdomains; pr++) {
        n_element_lev +=
            off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
        }
    // zeroing matrix
    for (int im = 0; im < n_element_lev * ndof_row_ql; im++) {
        mem[im] = 0;		// row index mem position
        }
    for (int im = 0; im < n_element_lev * NDOF_FEM; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;

            int sum_lev = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int i_row =
                    (iel - off_el[pr * _n_levels + Level1] +
                     sum_lev) * ndof_row_ql + kk;
                int ind_row = i_row * NDOF_FEM;
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int opj = v_inv_nd[elem_sto[el_lib][1 + kj]];	// node top
                    Mat_q[ind_row + mem[i_row] + kj] = g_indexL[Level_clmn][opj];	// pos level Level_clmn
                    }
                }
            }
        }

    // compressing zeros and counting -------------------------------------
    // counters

    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    int ki = 0;
    // zeroing aux
    for (int i = 0; i < NDOF_FEM * n_element_lev; i++) {
        aux[i] = 0;
        }
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                off_el[kdom * _n_levels + Level1 + 1] - off_el[kdom * _n_levels +
                        Level1]; iki++) {
            for (int kii = 0; kii < ndof_row_ql; kii++) {
                // setting aux
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int knod_q = Mat_q[NDOF_FEM * ki + kj];
                    if (knod_q >= 0) {
                        aux[knod_q] = 1;
                        }
                    }
                // computing various flags and lengths
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int knod_q = Mat_q[NDOF_FEM * ki + kj];
                    if (knod_q >= 0 && aux[knod_q] == 1) {
                        aux[knod_q] = 0;
                        Mat_q[count_q] = knod_q;
                        count_q++;
                        if (knod_q < (int) offM[kdom]
                                || knod_q >= (int) offM[kdom + 1]) {
                            count_qoff++;
                            }
                        }
                    }
                len_q[ki + 1] = count_q;
                len_qoff[ki + 1] = count_qoff;
                ki++;
                }
            }
        }
    return count_q;
    }




// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_kk(int Level1,	// level
                          int /*Level_row */ ,	// row level
                          int /*Level_clmn */ ,	// column level
                          int ndof_row_ql,	// row nodes
                          int ndof_clmn_ql,	// column nodes
                          int /*n_nodes */ ,	// top nodes
                          int /*max_elnd */ ,
                          int *off_el,
                          std::vector < std::pair < int, int > > /*v_el */ ,
                          int ** /*elem_sto */ , int * /*v_inv_nd */ ,
                          int ** /*g_indexL */ , int * /*off_nd */ ,
                          int /*offM */ [],
                          int aux[],
                          int mem[],
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    int n_element_lev = 0;
    for (int pr = 0; pr < _n_subdomains; pr++) {
        n_element_lev +=
            off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
        }
    // zeroing matrix
    for (int im = 0; im < n_element_lev * ndof_row_ql; im++) {
        mem[im] = 0;		// row index mem position
        }
    for (int im = 0; im < n_element_lev * ndof_clmn_ql; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
//       int el_lib=v_el[iel].second; // libmesh element
            // sum_lev= number of elements (at Level1) with proc < actual proc
            int sum_lev = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }

            int iel_lev = iel - off_el[pr * _n_levels + Level1] + sum_lev;	// element at Level1
            int idof_lev = iel_lev * ndof_clmn_ql;	// element dofs at Level1

            for (int ki = 0; ki < ndof_row_ql; ki++) {
                int row_indx = (idof_lev + ki) * ndof_clmn_ql;	// matrix row
                for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                    Mat_q[row_indx + kj] = idof_lev + kj;
                    }
                }
            }
        }

    // compressing zeros and counting -------------------------------------
    // counters
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    int ki = 0;
    // zeroing aux
    for (int i = 0; i < ndof_clmn_ql * n_element_lev; i++) {
        aux[i] = 0;
        }
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        int delta_dof =
            (off_el[kdom * _n_levels + Level1 + 1] -
             off_el[kdom * _n_levels + Level1]) * ndof_row_ql;
        for (int iki = 0; iki < delta_dof; iki++) {
            // setting aux
            for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[ndof_clmn_ql * ki + kj];
                if (knod_q >= 0) {
                    aux[knod_q] = 1;
                    }
                }
            int dof_start_lev = 0;	// dof start for proc (kdom)
            for (int jpr = 0; jpr < kdom; jpr++) {
                dof_start_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            dof_start_lev *= ndof_row_ql;

            // computing various flags and lengths
            for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[ndof_clmn_ql * ki + kj];
                if (knod_q >= 0 && aux[knod_q] == 1) {
                    aux[knod_q] = 0;
                    Mat_q[count_q] = knod_q;
                    count_q++;
                    if (knod_q < dof_start_lev
                            || knod_q >= delta_dof + dof_start_lev) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }
    return count_q;
    }

// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_qk(int Level1,	// level
                          int Level_row,	// row level
                          int /*Level_clmn */ ,	// column level
                          int ndof_row_ql,	// row nodes
                          int ndof_clmn_ql,	// column nodes
                          int n_nodes,	// top nodes
                          int max_elnd,
                          int *off_el,
                          std::vector < std::pair < int, int > >v_el,
                          int **elem_sto, int *v_inv_nd,
                          int **g_indexL, int *off_nd, int /*offM */ [],
                          int aux[],
                          int mem[],
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    // zeroing matrix
    for (int im = 0; im < n_nodes; im++) {
        mem[im] = 0;
        }
    for (int im = 0; im < n_nodes * max_elnd * ndof_clmn_ql; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int sum_lev = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int op = v_inv_nd[elem_sto[el_lib][1 + kk]];	// libmesh number node
                int i_row = g_indexL[Level_row][op];	// my  number node
                for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                    int opj = (iel - off_el[pr * _n_levels + Level1] + sum_lev) * ndof_clmn_ql + kj;	// dof in each elem (pie const)
                    int ind_row = i_row * max_elnd * ndof_clmn_ql;
                    Mat_q[ind_row + mem[i_row] + kj] = opj;
                    }
                mem[i_row] += ndof_clmn_ql;
                }
            }
        }
    // compressing zeros and counting -------------------------------------
    // counters
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    // zeroing aux
    for (int i = 0; i < max_elnd * n_nodes * ndof_clmn_ql; i++) {
        aux[i] = 0;
        }
    int ki = 0;			//  ki=global matrix row (level)
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                (int)(off_nd[kdom * _n_levels + Level1 + 1] -
                      off_nd[kdom * _n_levels]); iki++) {
            // setting aux
            for (int kj = 0; kj < max_elnd * ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[max_elnd * ndof_clmn_ql * ki + kj];
                if (knod_q >= 0) {
                    aux[knod_q] = 1;
                    }
                }
            int sum_lev = 0;
            for (int jpr = 0; jpr < kdom; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            // computing various flags and lengths
            for (int kj = 0; kj < max_elnd * ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[max_elnd * ndof_clmn_ql * ki + kj];
                if (knod_q >= 0 && aux[knod_q] == 1) {
                    aux[knod_q] = 0;
                    Mat_q[count_q] = knod_q;
                    count_q++;
                    if (knod_q < sum_lev * ndof_clmn_ql
                            || knod_q >=
                            (off_el[kdom * _n_levels + Level1 + 1] -
                             off_el[kdom * _n_levels + Level1] +
                             sum_lev) * ndof_clmn_ql) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }

    return count_q;
    }

// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_lk(int Level1,	// level
                          int Level_row,	// row level
                          int /*Level_clmn */ ,	// column level
                          int ndof_row_ql,	// row nodes
                          int ndof_clmn_ql,	// column nodes
                          int n_nodes,	// top nodes
                          int max_elnd,
                          int *off_el,
                          std::vector < std::pair < int, int > >v_el,
                          int **elem_sto, int *v_inv_nd,
                          int **g_indexL, int *off_nd, int /*offM */ [],
                          int aux[],
                          int mem[],
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    // zeroing matrix
    for (int im = 0; im < n_nodes; im++) {
        mem[im] = 0;
        }
    for (int im = 0; im < n_nodes * max_elnd * ndof_clmn_ql; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int sum_lev = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int op = v_inv_nd[elem_sto[el_lib][1 + kk]];	// libmesh number node
                int i_row = g_indexL[Level_row][op];	// my  number node
                for (int kj = 0; kj < ndof_clmn_ql; kj++) {
                    int opj = (iel - off_el[pr * _n_levels + Level1] + sum_lev) * ndof_clmn_ql + kj;	// dof in each elem (pie const)
                    int ind_row = i_row * max_elnd * ndof_clmn_ql;
                    Mat_q[ind_row + mem[i_row] + kj] = opj;
                    }
                mem[i_row] += ndof_clmn_ql;
                }
            }
        }
// //   printf("*************************************************************************** %d  ", ndof_clmn_ql);
//   for (int i=0;i < ndof_row_ql; i++){
//     for (int j=0;j < mem[i]; j++) {
//      printf("%d ",Mat_q[i*max_elnd*ndof_clmn_ql+j]);
//     }
//     printf("\n");
// }
    // compressing zeros and counting -------------------------------------
    // counters
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    // zeroing aux
    for (int i = 0; i < max_elnd * n_nodes * ndof_clmn_ql; i++) {
        aux[i] = 0;
        }
    int ki = 0;			//  ki=global matrix row (level)
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                (int)(off_nd[kdom * _n_levels + Level1 + 1] -
                      off_nd[kdom * _n_levels]); iki++) {
            // setting aux
            for (int kj = 0; kj < max_elnd * ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[max_elnd * ndof_clmn_ql * ki + kj];
                if (knod_q >= 0) {
                    aux[knod_q] = 1;
                    }
                }
            int sum_lev = 0;
            for (int jpr = 0; jpr < kdom; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            // computing various flags and lengths
            for (int kj = 0; kj < max_elnd * ndof_clmn_ql; kj++) {
                int knod_q = Mat_q[max_elnd * ndof_clmn_ql * ki + kj];
                if (knod_q >= 0 && aux[knod_q] == 1) {
                    aux[knod_q] = 0;
                    Mat_q[count_q] = knod_q;
                    count_q++;
                    if (knod_q < sum_lev * ndof_clmn_ql
                            || knod_q >=
                            (off_el[kdom * _n_levels + Level1 + 1] -
                             off_el[kdom * _n_levels + Level1] +
                             sum_lev) * ndof_clmn_ql) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }

    return count_q;
    }

// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_mat_kl(int Level1,	// level
                          int Level_row,	// row level
                          int Level_clmn,	// column level
                          int ndof_row_ql,	// row nodes
                          int ndof_clmn_ql,	// column nodes
                          int n_nodes,	// top nodes
                          int max_elnd,
                          int *off_el,
                          std::vector < std::pair < int, int > >v_el,
                          int **elem_sto, int *v_inv_nd,
                          int **g_indexL, int *off_nd,
                          int offM[],
                          int aux[],
                          int mem[],
                          int Mat_q[], int len_q[], int len_qoff[]) {
    // ================================================

    int n_element_lev = 0;
    for (int pr = 0; pr < _n_subdomains; pr++) {
        n_element_lev +=
            off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
        }
    // zeroing matrix
    for (int im = 0; im < n_element_lev * ndof_row_ql; im++) {
        mem[im] = 0;		// row index mem position
        }
    for (int im = 0; im < n_element_lev * NDOF_FEM; im++) {
        Mat_q[im] = -1;
        }

    // generation matrix entries
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;

            int sum_lev = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_lev +=
                    off_el[jpr * _n_levels + Level1 + 1] -
                    off_el[jpr * _n_levels + Level1];
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int i_row =
                    (iel - off_el[pr * _n_levels + Level1] +
                     sum_lev) * ndof_row_ql + kk;
                int ind_row = i_row * NDOF_FEM;
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int opj = v_inv_nd[elem_sto[el_lib][1 + kj]];	// node top
                    Mat_q[ind_row + mem[i_row] + kj] = g_indexL[Level_clmn][opj];	// pos level Level_clmn
                    }
                }
            }
        }

    // compressing zeros and counting -------------------------------------
    // counters

    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    int ki = 0;
    // zeroing aux
    for (int i = 0; i < NDOF_FEM * n_element_lev; i++) {
        aux[i] = 0;
        }
    // compressing
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                off_el[kdom * _n_levels + Level1 + 1] - off_el[kdom * _n_levels +
                        Level1]; iki++) {
            for (int kii = 0; kii < ndof_row_ql; kii++) {
                // setting aux
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int knod_q = Mat_q[NDOF_FEM * ki + kj];
                    if (knod_q >= 0) {
                        aux[knod_q] = 1;
                        }
                    }
                // computing various flags and lengths
                for (int kj = 0; kj < NDOF_FEM; kj++) {
                    int knod_q = Mat_q[NDOF_FEM * ki + kj];
                    if (knod_q >= 0 && aux[knod_q] == 1) {
                        aux[knod_q] = 0;
                        Mat_q[count_q] = knod_q;
                        count_q++;
                        if (knod_q < (int) offM[kdom]
                                || knod_q >= (int) offM[kdom + 1]) {
                            count_qoff++;
                            }
                        }
                    }
                len_q[ki + 1] = count_q;
                len_qoff[ki + 1] = count_qoff;
                ki++;
                }
            }
        }
    return count_q;
    }

// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_prol_qq(int Level1,	// level
                           int Level_row,	// row level
                           int Level_clmn,	// column level
                           int ndof_row_ql,	// element dofs
                           int n_nodes_f, int *off_el,	// offset elements vector
                           std::vector < std::pair < int, int > >v_el,	// element ordering
                           int **elem_sto,	// element storage
                           int *v_inv_nd,	// node ordering
                           int **g_indexL,	// map nodes
                           int *off_nd,	// offset nodes vector
                           double values_q[],	// prolongation nonzero values
                           int Prol_q[],	// prolongation compressed pos
                           int len_q[],	// # of diagonal nonzero row entries
                           int len_qoff[]	// # of offdiagonal nonzero row entries
                          ) {
    // ======================================================

    // zeroing matrix
    for (int im = 0; im < n_nodes_f * NDOF_FEM; im++) {
        Prol_q[im] = -1;
        values_q[im] = 0.;
        }

    // Prolongation Operator linear ---------------------------------------------------
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int par = elem_sto[el_lib][NDOF_FEM + 3];
            int ch = 0;
            int ch_id = 0;
            while (el_lib != elem_sto[par][NDOF_FEM + 5 + ch]) {
                ch_id = ++ch;
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int el_n = v_inv_nd[elem_sto[el_lib][1 + kk]];
                for (int kj = 0; kj < ndof_row_ql; kj++) {
                    int par_n = v_inv_nd[elem_sto[par][1 + kj]];
                    Prol_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj] =
                        g_indexL[Level_clmn][par_n];
                    if (ndof_row_ql == NDOF_P)
                        values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj] =
                            _geomel._embedding_matrix_l[ch_id][kk][kj];
                    else
                        values_q[g_indexL[Level_row][el_n] * NDOF_FEM + kj] =
                            _geomel._embedding_matrix_q[ch_id][kk][kj];
                    }
                }
            }
        }

    //Setting boundary limits ---------------------------------
    int *offPc = new int[_n_subdomains + 1];
    offPc[0] = 0;
    for (int idom = 0; idom < _n_subdomains; idom++) {
        offPc[idom + 1] =
            offPc[idom] + (off_nd[idom * _n_levels + Level1] -
                           off_nd[idom * _n_levels]);
        }
    offPc[_n_subdomains] += 1;

    // Compressing prolongation ---------------------------
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;

    int ki = 0;
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0;
                iki <
                off_nd[kdom * _n_levels + Level1 + 1] - off_nd[kdom * _n_levels];
                iki++) {
            for (int kj = 0; kj < NDOF_FEM; kj++) {
                int kk = NDOF_FEM * ki + kj;
                if (fabs(values_q[kk]) > 1.e-8) {
                    Prol_q[count_q] = Prol_q[kk];
                    values_q[count_q] = values_q[kk];
                    count_q++;
                    if (Prol_q[kk] >= (int) offPc[kdom + 1]
                            || Prol_q[kk] < (int) offPc[kdom]) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }

    // clean -------------------------------------------------
    delete[]offPc;
    return count_q;
    }


// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_prol_kk(int Level1,	// level
                           int Level_row,	// row level
                           int Level_clmn,	// column level
                           int ndof_row_ql,	// element dofs
                           int n_nodes_f, int *off_el,	// offset elements vector
                           std::vector < std::pair < int, int > >v_el,	// element ordering
                           int **elem_sto,	// element storage
                           int *v_inv_nd,	// node ordering
                           int *v_inv_el, int **g_indexL,	// map nodes
                           int *off_nd,	// offset nodes vector
                           double values_q[],	// prolongation nonzero values
                           int Prol_q[],	// prolongation compressed pos
                           int len_q[],	// # of diagonal nonzero row entries
                           int len_qoff[]	// # of offdiagonal nonzero row entries
                          ) {
    // ======================================================

    // zeroing matrix
    for (int im = 0; im < Level_row * 1; im++) {
        Prol_q[im] = -1;
        values_q[im] = 0.;
        }

    // Simple !!! Prolongation Operator linear ---------------------------------------------------
    int sum_lev = 0;
    int count_q = 0;
    len_q[0] = 0;
    len_qoff[0] = 0;
    for (int pr = 0; pr < _n_subdomains; pr++) {
        int delta =
            off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int i_row = iel - off_el[pr * _n_levels + Level1] + sum_lev;
            int par = elem_sto[el_lib][NDOF_FEM + 3];
            int j_par = v_inv_el[par];
            int sum_levm1 = 0;
            for (int jpr = 0; jpr < pr; jpr++) {
                sum_levm1 +=
                    off_el[jpr * _n_levels + Level1] - off_el[jpr * _n_levels +
                            Level1 - 1];
                }
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int ind_row = i_row * ndof_row_ql + kk;
                int j_column =
                    (j_par - off_el[pr * _n_levels + Level1 - 1] +
                     sum_levm1) * ndof_row_ql + kk;
                Prol_q[ind_row] = j_column;
                values_q[ind_row] = 1.;
                len_q[ind_row + 1] = count_q + 1;
                len_qoff[ind_row + 1] = 0;
                count_q++;
                }
            }
        sum_lev += delta;
        }
    return count_q;
    }

// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_res_qq(int Level1,	// level
                          int Level_row,	// row level
                          int Level_cln,	// colum level
                          int ndof_row_ql,	// element dofs
                          int n_nodes_row,	// # of row nodes
                          int n_nodes_cln,	// # of column nodes
                          int max_elnd,	// max # of elements at each point
                          int *off_el,	// offset elements vector
                          std::vector < std::pair < int, int > >v_el,	// element ordering
                          int **elem_sto,	// element storage
                          int *v_inv_nd,	// node ordering
                          int **g_indexL,	// map nodes
                          int *off_nd,	// offset nodes vector
                          int mem[],	//  memory index
                          double values_q[],	// restriction nonzero values
                          int Rest_q[],	// restriction compressed pos
                          int len_q[],	// # of diagonal nonzero row entries
                          int len_qoff[]	// # of offdiagonal nonzero row entries
                         ) {
    // ====================================================

    // zeroing the structures
    for (int i = 0; i < n_nodes_row * NDOF_FEM * max_elnd; i++) {
        Rest_q[i] = -1;
        values_q[i] = 0.;
        }
    for (int im = 0; im < n_nodes_row; im++) {
        mem[im] = 0;
        }

    // Restriction operator
    for (int pr = 0; pr < _n_subdomains; pr++) {
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int n_childs = elem_sto[el_lib][NDOF_FEM + 4];	//(*it_n)->n_children();
            for (int i_ch = 0; i_ch < n_childs; i_ch++) {
                int ch_el = elem_sto[el_lib][NDOF_FEM + 5 + i_ch];	// Elem* child=(*it_n)->child(i_ch);
                for (int kk = 0; kk < ndof_row_ql; kk++) {
                    int op = v_inv_nd[elem_sto[el_lib][1 + kk]];
                    int irow = g_indexL[Level_row][op];
                    for (int kj = 0; kj < ndof_row_ql; kj++) {
                        int ch_op = v_inv_nd[elem_sto[ch_el][1 + kj]];

#ifdef REST_SIMPLE
                        if (op == ch_op) {
                            Rest_q[irow * NDOF_FEM * max_elnd] =
                                g_indexL[Level_cln][ch_op];
                            values_q[irow * NDOF_FEM * max_elnd] = NDOF_P;
                            if (ndof_row_ql == NDOF_P) {
                                values_q[irow * NDOF_FEM * max_elnd] = NDOF_P;
                                }
                            mem[irow] = 1;
                            }
#else
                        double val_qq = 0;
                        if (ndof_row_ql == NDOF_P) {
                            val_qq = _geomel._embedding_matrix_l[i_ch][kj][kk];
                            }
                        else {
                            val_qq = _geomel._embedding_matrix_q[i_ch][kj][kk];
                            }
                        if (fabs(val_qq) > 1.e-6) {

                            int jcln = g_indexL[Level_cln][ch_op];
                            int indx0 = irow * NDOF_FEM * max_elnd;
                            int indx_mem = mem[irow];
                            int flag = 0;
                            for (int ks = 0; ks < indx_mem; ks++)
                                if (Rest_q[indx0 + ks] == jcln) {
                                    flag = 1;
                                    }

                            if (flag == 0) {
                                Rest_q[indx0 + indx_mem] = jcln;
                                values_q[indx0 + indx_mem] = val_qq;
                                mem[irow]++;
                                }

                            }
#endif
                        }

                    }
                }
            }
        }

    // compressing zeros
    int *offPc = new int[_n_subdomains + 1];
    offPc[0] = 0;
    for (int idom = 0; idom < _n_subdomains; idom++) {
        offPc[idom + 1] =
            offPc[idom] + (off_nd[idom * _n_levels + Level1 + 2] -
                           off_nd[idom * _n_levels]);
        }
    offPc[_n_subdomains] += 1;
    //  Compressing quadratic restrictor
    int count_q = 0;
    len_q[0] = count_q;
    int count_qoff = 0;
    len_qoff[0] = count_qoff;
    int ki = 0;
    for (int kdom = 0; kdom < _n_subdomains; kdom++) {
        for (int iki = 0; iki < (int)(off_nd[kdom * _n_levels + Level1 + 1]
                                      - off_nd[kdom * _n_levels]); iki++) {
            // computing various flags and lengths
            for (int kj = 0; kj < NDOF_FEM * max_elnd; kj++) {
                int kk = NDOF_FEM * max_elnd * ki + kj;
                if (fabs(values_q[kk]) > 1.e-8) {
                    Rest_q[count_q] = Rest_q[kk];
                    values_q[count_q] = values_q[kk];
                    count_q++;
                    if (Rest_q[kk] >= (int) offPc[kdom + 1]
                            || Rest_q[kk] < (int) offPc[kdom]) {
                        count_qoff++;
                        }
                    }
                }
            len_q[ki + 1] = count_q;
            len_qoff[ki + 1] = count_qoff;
            ki++;
            }
        }

    delete[]offPc;
    return count_q;
    }


// =======================================================
/// This function computes the quadratic matrix
int
MGGenCase::compute_res_kk(int Level1,	// level
                          int Level_row,	// row level
                          int Level_cln,	// colum level
                          int ndof_row_ql,	// element dofs
                          int n_nodes_row,	// # of row nodes
                          int n_nodes_cln,	// # of column nodes
                          int max_elnd,	// max # of elements at each point
                          int *off_el,	// offset elements vector
                          std::vector < std::pair < int, int > >v_el,	// element ordering
                          int **elem_sto,	// element storage
                          int *v_inv_nd,	// node ordering
                          int *v_inv_el, int **g_indexL,	// map nodes
                          int *off_nd,	// offset nodes vector
                          int mem[],	//  memory index
                          double values_q[],	// restriction nonzero values
                          int Rest_q[],	// restriction compressed pos
                          int len_q[],	// # of diagonal nonzero row entries
                          int len_qoff[]	// # of offdiagonal nonzero row entries
                         ) {
    // ====================================================

    // zeroing the structures
    for (int i = 0; i < n_nodes_row * max_elnd; i++) {
        Rest_q[i] = -1;
        values_q[i] = 0.;
        }
    for (int im = 0; im < n_nodes_row; im++) {
        mem[im] = 0;
        }

    // Restriction operator
    int sum_lev = 0;
    int count_q = 0;
    len_q[0] = 0;
    len_qoff[0] = 0;
    for (int pr = 0; pr < _n_subdomains; pr++) {
        int delta =
            off_el[pr * _n_levels + Level1 + 1] - off_el[pr * _n_levels + Level1];
        for (int iel = off_el[pr * _n_levels + Level1];
                iel < off_el[pr * _n_levels + Level1 + 1]; iel++) {
            int el_lib = v_el[iel].second;
            int n_childs = elem_sto[el_lib][NDOF_FEM + 4];	//(*it_n)->n_children();
            for (int kk = 0; kk < ndof_row_ql; kk++) {
                int i_row =
                    kk + (iel - off_el[pr * _n_levels + Level1] +
                          sum_lev) * ndof_row_ql;
                int ind_row = i_row * n_childs;
                for (int i_ch = 0; i_ch < n_childs; i_ch++) {
                    int ch_el = elem_sto[el_lib][NDOF_FEM + 5 + i_ch];	// Elem* child=(*it_n)->child(i_ch);
                    int j_par = v_inv_el[ch_el];
                    int sum_levp1 = 0;
                    for (int jpr = 0; jpr < pr; jpr++) {
                        sum_levp1 +=
                            off_el[jpr * _n_levels + Level1 + 2] -
                            off_el[jpr * _n_levels + Level1 + 1];
                        }
                    int j_column =
                        kk + (j_par - off_el[pr * _n_levels + Level1 + 1] +
                              sum_levp1) * ndof_row_ql;
                    double val_qq = 1.;
                    Rest_q[ind_row + i_ch] = j_column;
                    values_q[ind_row + i_ch] = val_qq;
                    count_q++;
                    }
                len_q[i_row + 1] = count_q;
                len_qoff[i_row + 1] = 0;
                }
            }
        sum_lev += delta;
        }
    return count_q;
    }




// ====================================================
// Print mesh in med format:
// conversion form libmesh to med format
void
MGGenCase::print_med(int Level,	// Level
                     BoundaryMesh & bd_msh0,	// boundary coarse mesh
                     Mesh & msh0,	// coarse mesh
                     int n_groups, int *group_id_names, std::vector < std::pair < int, int > >v,	// node order
                     std::vector < std::pair < int, int > >v_elb, int *nod_flag,	// bc node
                     int *mat_flag,	// mat element
                     std::string filename	// filename
                    ) {
    // ==================================================
#ifdef HAVE_MED
    // name file directory
    std::string input_dir = _mgutils._mesh_dir;
    std::ostringstream name;
    name << input_dir << filename << "_gen.med";

    // coarse  mesh
    int n_nodes = msh0.n_nodes();	// from mesh
    int n_elements = msh0.n_elem();	// from mesh
    //  boundary coarse mesh
    int bd_n_nodes = bd_msh0.n_nodes();	// from boundary mesh
    int bd_n_elements = bd_msh0.n_elem();	// from boundary mesh

    // cordinates
    double *coord;
    coord = new double[n_nodes * _dim];
    // connectivity (vol+boundary)
    int *conn;
    conn = new int[n_elements * NDOF_FEM];
    int *conn_bd;
    conn_bd = new int[bd_n_elements * NDOF_FEMB];
    int *elem_bd_id2;
    elem_bd_id2 = new int[bd_n_elements];
    // element (nodes for element)
    int n_nodes_el;
    int bd_n_nodes_el;
    int count_eb = 0;

// med->libmesh map  second order only
#if ELTYPE==27
#if DIMENSION==3
    const unsigned int nodesinv[] = {
        4, 7, 3, 0, 5, 6, 2, 1, 19, 15, 11, 12, 17, 14, 9, 13,
        16, 18, 10, 8, 24, 25, 23, 20, 21, 22, 26
        };
    const unsigned int nodesinvbd[] = { 3, 0, 1, 2, 7, 4, 5, 6, 8 };
#else
    const unsigned int nodesinv[] = { 3, 0, 1, 2, 7, 4, 5, 6, 8 };

    const unsigned int nodesinvbd[] = { 0, 2, 1 };
#endif
#endif
#if ELTYPE==10
#if DIMENSION==3
    const unsigned int nodesinv[] = { 2, 3, 1, 0, 9, 8, 5, 6, 7, 4 };
    const unsigned int nodesinvbd[] = { 1, 2, 0, 4, 5, 3 };
#else
    const unsigned int nodesinv[] = { 1, 2, 0, 4, 5, 3 };
    const unsigned int nodesinvbd[] = { 0, 2, 1 };
#endif
#endif
    //  mesh (volume) -----------------------------------------------------
    int *map_elem = new int[n_elements];
    int n_element_top = 0;
    Mesh::const_element_iterator it_tr = msh0.elements_begin();
    const Mesh::const_element_iterator end_tr = msh0.elements_end();
    for (; it_tr != end_tr; ++it_tr) {
        Elem *elem = *it_tr;	// element
        int lev = elem->level();
        if (lev == Level) {
            int id_el = n_element_top;
            map_elem[id_el] = elem->id();	// element id
            n_nodes_el = elem->n_nodes();	// number of element nodes
            for (int inode = 0; inode < n_nodes_el; inode++) {
                int knode = elem->node(nodesinv[inode]);	// global node through map
                conn[id_el * n_nodes_el + inode] = knode;	// connectivity

                // coordinates storage
                for (int idim = 0; idim < _dim; idim++) {
                    double xyz = msh0.point(knode)(idim);
                    coord[knode * _dim + idim] = xyz;	// med cordinate tuple
                    }
                }
            n_element_top++;
            }
        //  mesh (boundary) ----------------------------------------
        if (lev == Level) {
            for (int s = 0; s < (int) elem->n_sides(); s++) {
                if (elem->neighbor(s) == NULL) {
                    UniquePtr < Elem > side(elem->build_side(s));	// face element
                    bd_n_nodes_el = (int) side->n_nodes();	// face nodes
                    int min = 100000;
                    for (int ns = 0; ns < bd_n_nodes_el; ns++) {

                        conn_bd[count_eb * bd_n_nodes_el + ns] =
                            side->node(nodesinvbd[ns]);
                        if (min >
                                nod_flag[conn_bd[count_eb * bd_n_nodes_el + ns]]) {
                            min =
                                nod_flag[conn_bd[count_eb * bd_n_nodes_el + ns]];
                            }

                        }
                    elem_bd_id2[count_eb] = min;
                    count_eb++;	// counter
                    }
                }
            }
        }
//   int *conn_med= new int[n_element_top];
//   for(int ie=0; ie<n_element_top; ie++) {
//         for(int inode=0; inode<n_nodes_el; inode++) {
//            conn_med[ie+inode*n_element_top]= conn[ie*n_nodes_el+inode];
//   }
    bd_n_elements = count_eb;
    assert(bd_n_elements == count_eb);	// check
    std::cout << " Printing med file " << name.
              str().c_str() << "\n for bc and mat: nodes =" << n_nodes <<
              "; elements =" << n_elements << "; n_nodes_el= " << n_nodes_el << std::
              endl;

    // MED mesh *************************************************

    // MEDCouplingUMesh mesh connectivity (volume)
    MEDCoupling::MEDCouplingUMesh * mesh1 =
        MEDCoupling::MEDCouplingUMesh::New("Mesh_1", _dim);
    mesh1->allocateCells(n_element_top);
    for (int i = 0; i < n_element_top; i++) {
        if (MED_EL_TYPE!=n_nodes_el) {
            std::cout<<"Attention!!! MED_EL_TYPE!=n_nodes_el in MGGenCase::print_med"<<std::endl;
            std::cout<<"Usually this means that your mesh is not linear or bi-quadratic"<<std::endl;
            abort();
        }
        mesh1->insertNextCell(MED_EL_TYPE, n_nodes_el, conn + i * n_nodes_el);
        }
    mesh1->finishInsertingCells();

    // MEDCouplingUMesh Mesh connectivity (boundary)
    MEDCoupling::MEDCouplingUMesh * mesh2 =
        MEDCoupling::MEDCouplingUMesh::New("Mesh_1", DIMENSION - 1);
    mesh2->allocateCells(bd_n_elements);
    for (int i = 0; i < bd_n_elements; i++) {
        mesh2->insertNextCell(MED_EL_BDTYPE, bd_n_nodes_el,
                              conn_bd + i * bd_n_nodes_el);
        }
    mesh2->finishInsertingCells();

    // coord (same node set for both meshes)
    MEDCoupling::DataArrayDouble * coordarr =
        MEDCoupling::DataArrayDouble::New();
    coordarr->alloc(n_nodes, _dim);
    std::copy(coord, coord + n_nodes * _dim, coordarr->getPointer());
    mesh1->setCoords(coordarr);
    mesh2->setCoords(coordarr);

    // Setting MEDCouplingUMesh into MEDFileUMesh
    MEDCoupling::MEDFileUMesh * mm = MEDCoupling::MEDFileUMesh::New();
    mm->setName("Mesh_1");	//name needed to be non empty
    mm->setDescription("Description Mesh_1");
    mm->setCoords(mesh1->getCoords());
    mm->setMeshAtLevel(0, mesh1, false);
    mm->setMeshAtLevel(-1, mesh2, false);

    // Volume Groups
//   int n_vol_group_nodes=0;
//   for(int i=0;i<n_groups_names;i++) if(group_id_names[i]<10) n_vol_group_nodes;

    std::map < int, std::vector < int >>vol_group_nodes;
    std::map < int, int >vol_group;
    for (int i = 0; i < n_element_top; i++) {
        vol_group[mat_flag[map_elem[i]]]++;
        vol_group_nodes[mat_flag[map_elem[i]]].push_back(i);
        }

    int n_vol_group = vol_group.size();
    std::vector < const MEDCoupling::DataArrayInt * >gr_vol(n_vol_group);
    MEDCoupling::DataArrayInt ** g_vol =
        new MEDCoupling::DataArrayInt *[n_vol_group];

    int js = 0;
    // defining the vol group data to store
    std::map < int, std::vector < int >>::iterator it_vol;
    for (it_vol = vol_group_nodes.begin(); it_vol != vol_group_nodes.end();
            ++it_vol) {
        int igroup = it_vol->first;
        int is = it_vol->second.size();
        g_vol[js] = MEDCoupling::DataArrayInt::New();
        g_vol[js]->alloc(is, 1);
        std::ostringstream name_p;
        name_p << igroup;
        g_vol[js]->setName(name_p.str().c_str());
        int *val1 = new int[is];
        for (int iv = 0; iv < is; iv++) {
            val1[iv] = it_vol->second[iv];
            }
        std::copy(val1, val1 + is, g_vol[js]->getPointer());
        delete[]val1;
        gr_vol[js] = g_vol[js];
        js++;
        }
    // inserting  the volume groups into the med-file
    mm->setGroupsAtLevel(0, gr_vol, false);


    // Boundary group ******************************************
    // Finding bd_group and  bd_group_nodes map
    std::map < int, std::vector < int >>bd_group_nodes;
    std::map < int, int >bd_group;
    for (int i = 0; i < bd_n_elements; i++) {
        bd_group[elem_bd_id2[i]]++;
        bd_group_nodes[elem_bd_id2[i]].push_back(i);
        }
    int n_bd_group = bd_group.size();
    // group vector
    std::vector < const MEDCoupling::DataArrayInt * >gr_bd(n_bd_group);
    MEDCoupling::DataArrayInt ** g_bd =
        new MEDCoupling::DataArrayInt *[n_bd_group];

    js = 0;
    // defining the  group data to store
    std::map < int, std::vector < int >>::iterator it;
    for (it = bd_group_nodes.begin(); it != bd_group_nodes.end(); ++it) {
        int igroup = it->first;
        int is = it->second.size();
//     std::cout << igroup << "\n";
        g_bd[js] = MEDCoupling::DataArrayInt::New();
        g_bd[js]->alloc(is, 1);
        std::ostringstream name_p;
        name_p << igroup;
        g_bd[js]->setName(name_p.str().c_str());
        int *valb1 = new int[is];	/*int icount=0; */
        for (int iv = 0; iv < is; iv++) {
            valb1[iv] = it->second[iv];
            }
        std::copy(valb1, valb1 + is, g_bd[js]->getPointer());
        delete[]valb1;
        gr_bd[js] = g_bd[js];
        js++;
        }
    // insert the boundary groups into the med-file
    mm->setGroupsAtLevel(-1, gr_bd, false);

    // Printing Group and Family
    std::cout << "\n \n =================================== ";
    std::cout << "\n Group Names -> Families : \n";
    std::map < std::string, std::vector < std::string > >a =
        (mm->getGroupInfo());
    std::map < std::string, std::vector < std::string > >::iterator ita;
    for (ita = a.begin(); ita != a.end(); ++ita) {
        std::string igroup = ita->first;
        int is = ita->second.size();
        std::cout << "\n " << igroup << "-> ";
        for (int i = 0; i < is; i++) {
            std::string a2 = ita->second[i];
            std::cout << a2 << "  ";
            }
        }
    std::cout << "\n \n =================================== ";
    std::cout << " \n Families -> Groups id: \n";
    std::map < std::string, int >fama = mm->getFamilyInfo();
    std::map < std::string, int >::iterator itfa;
    for (itfa = fama.begin(); itfa != fama.end(); ++itfa) {
        std::string igroup = itfa->first;
        std::cout << "\n " << igroup << "-> ";
        int a2 = itfa->second;
        std::cout << a2 << "  ";
        }
    std::cout << "\n\n ";
    mm->write(name.str().c_str(), 2);

    // clean
    delete[]coord;
    fama.clear();
    a.clear();
    delete[]conn;
    coordarr->decrRef();
    mesh1->decrRef();
    mesh2->decrRef();
    mm->decrRef();

#else
    std::cout <<
              "\n \n MGGenCase::Print_med you don't have MED library included in your project\n\n";
#endif




    return;
    }



void
MGGenCase::print_MedToMg(int Level,	// Level
                         BoundaryMesh & bd_msh0,	// boundary coarse mesh
                         Mesh & msh0,	// coarse mesh
                         int n_groups, int *group_id_names, std::vector < std::pair < int, int > >v,	// node order
                         int *v_inv_nd, std::vector < std::pair < int, int > >v_elb, int *nod_flag,	// bc node
                         int *mat_flag,	// mat element
                         int *off_el, int *v_inv_el, std::vector < std::pair < int, int >>v_el, std::vector < int >ElementsPerLevel, int **elem_sto, std::string filename	// filename
                        ) {
    // ==================================================
    // THIS ROUTINE IS USED TO PRINT MED MESH FILES WITH FIELDS THAT HELP CREATING INTERFACES
    // IN PARTICULAR WE SAVE THE PROC BY PROC MESH PARTITIONING AND THE NODE RENUMBERING FOR PARALLEL DISTRIBUTION

    std::cout << "\033[1;31m\n............................................... \
                        \n PRINTING THE MED FILE FOR COUPLING INTERFACES \
                        \n...............................................\n\033[0m\n";

#ifdef HAVE_MED
    // name file directory
    std::string input_dir = _mgutils._mesh_dir;
    std::string info_name = input_dir + filename + "_MedToMg" + ".med";

    //  boundary coarse mesh
    int bd_n_nodes = bd_msh0.n_nodes();	// from boundary mesh
    int bd_n_elements = bd_msh0.n_elem();	// from boundary mesh
    int *conn_bd;
    conn_bd = new int[bd_n_elements * NDOF_FEMB];
    int *elem_bd_id2;
    elem_bd_id2 = new int[bd_n_elements];
    // element (nodes for element)
    int n_nodes_el;
    int NumBdNodes;
    int count_eb = 0;

// MAP BETWEEN MED AND LIBMESH ELEMENT NODE NUMBERING
#if ELTYPE==27
#if DIMENSION==3
    const unsigned int LibToMed[] = {
        7, 4, 5, 6, 3, 0, 1, 2, 19, 16, 17, 18, 11, 8, 9, 10, 15, 12, 13, 14,
        25, 24, 21, 22, 23, 20, 26
        };
    const unsigned int LibToMed_bd[] = { 3, 0, 1, 2, 7, 4, 5, 6, 8 };
#else
    const unsigned int LibToMed[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    const unsigned int LibToMed_bd[] = { 0, 1, 2 };
#endif
#endif
#if ELTYPE==10
#if DIMENSION==3
    const unsigned int LibToMed[] = { 2, 3, 1, 0, 9, 8, 5, 6, 7, 4 };
    const unsigned int LibToMed_bd[] = { 1, 2, 0, 4, 5, 3 };
#else
    const unsigned int LibToMed[] = { 1, 2, 0, 4, 5, 3 };
    const unsigned int LibToMed_bd[] = { 0, 2, 1 };
#endif
#endif


    //  mesh (volume) -----------------------------------------------------
    std::vector < int >ElemPerLevel(Level + 1);
    ElemPerLevel[0] = 0;
    if (Level > 0)
        for (int LEVEL = 1; LEVEL <= Level; LEVEL++) {
            ElemPerLevel[LEVEL] =
                ElemPerLevel[LEVEL - 1] + ElementsPerLevel[LEVEL - 1];
            }


    for (int LEVEL = 0; LEVEL <= Level; LEVEL++) {
        Mesh::const_element_iterator it_tr = msh0.level_elements_begin(LEVEL);
        const Mesh::const_element_iterator end_tr =
            msh0.level_elements_end(LEVEL);
        int n_element_top = 0;

        MEDCoupling::DataArrayDouble * Volcells =
            MEDCoupling::DataArrayDouble::New();
        Volcells->alloc(ElementsPerLevel[LEVEL], 1);
        double *CellMap = const_cast < double *>(Volcells->getPointer());
        MEDCoupling::DataArrayDouble * ProcArray =
            MEDCoupling::DataArrayDouble::New();
        ProcArray->alloc(ElementsPerLevel[LEVEL], 1);
        double *ProcMap = const_cast < double *>(ProcArray->getPointer());

        int *conn;
        conn = new int[ElementsPerLevel[LEVEL] * NDOF_FEM];
        std::vector < int >Nodes;

        for (; it_tr != end_tr; ++it_tr) {
            // LOOP OVER THE MESH ELEMENTS ------------------------
            Elem *elem = *it_tr;	// element
            int lev = elem->level();	// actual level of the mesh element
            int id_el = n_element_top;

            n_nodes_el = elem->n_nodes();	// number of element nodes
//             CellMap[id_el] = v_inv_el[elem->id()];
            CellMap[id_el] = v_el[id_el].second;
            ProcMap[id_el] =
                elem_sto[id_el + ElemPerLevel[LEVEL]][NDOF_FEM + 2];

            // VOLUME
            for (int inode = 0; inode < n_nodes_el; inode++) {
                int knode = elem->node(LibToMed[inode]);	// global node through map
                conn[id_el * n_nodes_el + inode] = knode;	// element connectivity
                Nodes.push_back(knode);
                }
            n_element_top++;

            // BOUNDARY
            for (int s = 0; s < (int) elem->n_sides(); s++) {
                // loop over element sides
                if (elem->neighbor(s) == NULL) {
                    // if neighbor == null then the side is on boundary
                    UniquePtr < Elem > side(elem->build_side(s));	// face element
                    NumBdNodes = (int) side->n_nodes();	// face nodes
                    int min = 100000;
                    for (int ns = 0; ns < NumBdNodes; ns++) {
                        conn_bd[count_eb * NumBdNodes + ns] = side->node(LibToMed_bd[ns]);	// connectivity of boundary elements
                        if (min > nod_flag[conn_bd[count_eb * NumBdNodes + ns]]) {
                            min = nod_flag[conn_bd[count_eb * NumBdNodes + ns]];	// we set the boundary flag of the node -> group id
                            }
                        }
                    elem_bd_id2[count_eb] = min;	// boundary flag of the element -> it is equal to the lowest bd node flag of the side
                    count_eb++;	// counter
                    }
                }
            }			// END LOOP OVER MESH ELEMENTS -------------------------------------------------------
        bd_n_elements = count_eb;	// total number of boundary elements (sides)
        assert(bd_n_elements == count_eb);	// check

        // MED mesh *************************************************
        MEDCoupling::MEDCouplingUMesh * VolumeMesh =
            MEDCoupling::MEDCouplingUMesh::New("Mesh_Lev_" +
                                               std::to_string(LEVEL), _dim);
        VolumeMesh->allocateCells(n_element_top);
        for (int i = 0; i < n_element_top; i++) {
            VolumeMesh->insertNextCell(MED_EL_TYPE, n_nodes_el,
                                       conn + i * n_nodes_el);
            }
        VolumeMesh->finishInsertingCells();
        // MEDCouplingUMesh Mesh connectivity (boundary)
        MEDCoupling::MEDCouplingUMesh * BoundaryMesh =
            MEDCoupling::MEDCouplingUMesh::New("Mesh_Lev_" +
                                               std::to_string(LEVEL), _dim - 1);
        BoundaryMesh->allocateCells(bd_n_elements);
        for (int i = 0; i < bd_n_elements; i++) {
            BoundaryMesh->insertNextCell(MED_EL_BDTYPE, NumBdNodes,
                                         conn_bd + i * NumBdNodes);
            }
        BoundaryMesh->finishInsertingCells();

        std::set < int >s(Nodes.begin(), Nodes.end());
        Nodes.clear();
        int LevelMeshNodes = s.size();
        double *Coords = new double[_dim * LevelMeshNodes];

        MEDCoupling::DataArrayDouble * FinerConn =
            MEDCoupling::DataArrayDouble::New();
        FinerConn->alloc(LevelMeshNodes, 1);
        double *MapArray = const_cast < double *>(FinerConn->getPointer());

        for (int node = 0; node < LevelMeshNodes; node++) {
            MapArray[node] = v_inv_nd[node];
            for (int dim = 0; dim < _dim; dim++) {
                Coords[node * _dim + dim] = msh0.point(node)(dim);
                }
            }

        MEDCoupling::DataArrayDouble * coordarr =
            MEDCoupling::DataArrayDouble::New();
        coordarr->alloc(LevelMeshNodes, _dim);
        std::copy(Coords, Coords + LevelMeshNodes * _dim,
                  coordarr->getPointer());
        VolumeMesh->setCoords(coordarr);
        BoundaryMesh->setCoords(coordarr);

        // Setting MEDCouplingUMesh into MEDFileUMesh
        MEDCoupling::MEDFileUMesh * mm = MEDCoupling::MEDFileUMesh::New();
        mm->setName("Mesh_Lev_" + std::to_string(LEVEL));	//name needed to be non empty
        mm->setDescription("Description Mesh_1");
        mm->setCoords(VolumeMesh->getCoords());
        mm->setMeshAtLevel(0, VolumeMesh, false);
        mm->setMeshAtLevel(-1, BoundaryMesh, false);

        if (LEVEL == Level) {
            // boundary and volume groups available only at finer level
            std::map < int, std::vector < int >>vol_group_nodes;
            std::map < int, int >vol_group;
            for (int i = 0; i < n_element_top; i++) {
                vol_group[mat_flag[(int) CellMap[i]+ElemPerLevel[LEVEL]]]++;    //ElemPerLevel[LEVEL] is the cumulative number of elems of the coarser Levels
                vol_group_nodes[mat_flag[(int) CellMap[i]+ElemPerLevel[LEVEL]]].push_back(i);
                }

            int n_vol_group = vol_group.size();
            std::vector <
            const MEDCoupling::DataArrayInt * >gr_vol(n_vol_group);
            MEDCoupling::DataArrayInt ** g_vol =
                new MEDCoupling::DataArrayInt *[n_vol_group];

            int js = 0;
            // defining the vol group data to store
            std::map < int, std::vector < int >>::iterator it_vol;
            for (it_vol = vol_group_nodes.begin();
                    it_vol != vol_group_nodes.end(); ++it_vol) {
                int igroup = it_vol->first;
                int is = it_vol->second.size();
                g_vol[js] = MEDCoupling::DataArrayInt::New();
                g_vol[js]->alloc(is, 1);
                std::ostringstream name_p;
                name_p << igroup;
                g_vol[js]->setName(name_p.str().c_str());
                int *val1 = new int[is];
                for (int iv = 0; iv < is; iv++) {
                    val1[iv] = it_vol->second[iv];
                    }
                std::copy(val1, val1 + is, g_vol[js]->getPointer());
                delete[]val1;
                gr_vol[js] = g_vol[js];
                js++;
                }
            // inserting  the volume groups into the med-file
            mm->setGroupsAtLevel(0, gr_vol, false);

            // Boundary group ******************************************
            // Finding bd_group and  bd_group_nodes map
            std::map < int, std::vector < int >>bd_group_nodes;
            std::map < int, int >bd_group;
            for (int i = 0; i < bd_n_elements; i++) {
                bd_group[elem_bd_id2[i]]++;
                bd_group_nodes[elem_bd_id2[i]].push_back(i);
                }
            int n_bd_group = bd_group.size();
            // group vector
            std::vector < const MEDCoupling::DataArrayInt * >gr_bd(n_bd_group);
            MEDCoupling::DataArrayInt ** g_bd =
                new MEDCoupling::DataArrayInt *[n_bd_group];

            js = 0;
            // defining the  group data to store
            std::map < int, std::vector < int >>::iterator it;
            for (it = bd_group_nodes.begin(); it != bd_group_nodes.end();
                    ++it) {
                int igroup = it->first;
                int is = it->second.size();
                g_bd[js] = MEDCoupling::DataArrayInt::New();
                g_bd[js]->alloc(is, 1);
                std::ostringstream name_p;
                name_p << igroup;
                g_bd[js]->setName(name_p.str().c_str());
                int *valb1 = new int[is];	/*int icount=0; */
                for (int iv = 0; iv < is; iv++) {
                    valb1[iv] = it->second[iv];
                    }
                std::copy(valb1, valb1 + is, g_bd[js]->getPointer());
                delete[]valb1;
                gr_bd[js] = g_bd[js];
                js++;
                }
            // insert the boundary groups into the med-file
            mm->setGroupsAtLevel(-1, gr_bd, false);
            }

        if (LEVEL == 0) {
            mm->write(info_name, 2);
            }
        else {
            mm->write(info_name, 1);
            }

        // CELL FIELD STORING THE MED-TO-MG MESH CELL NUMBERING
        MEDCoupling::MEDCouplingFieldDouble * CellField =
            MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
        CellField->setMesh(VolumeMesh);
        CellField->setArray(Volcells);
        CellField->setName("MG_cell_id_Lev_" + std::to_string(LEVEL));
        MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, CellField);
//         MEDCoupling::WriteField ( "MESH/ciao.med", CellField, true );

        // PROC FIELD STORING PROC ID
        MEDCoupling::MEDCouplingFieldDouble * ProcField =
            MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
        ProcField->setMesh(VolumeMesh);
        ProcField->setArray(ProcArray);
        ProcField->setName("Proc_Lev_" + std::to_string(LEVEL));
        MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, ProcField);

        // MAP TO MG NODE NUMBERING
        MEDCoupling::MEDCouplingFieldDouble * FinerLevelIDS =
            MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
        FinerLevelIDS->setMesh(VolumeMesh);
        FinerLevelIDS->setArray(FinerConn);
        FinerLevelIDS->setName("FinerLevelNodeIDS_Lev_" +
                               std::to_string(LEVEL));
        MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name,
                FinerLevelIDS);

        CellField->decrRef();
        ProcField->decrRef();
        FinerLevelIDS->decrRef();
        ProcArray->decrRef();
        FinerConn->decrRef();
        Volcells->decrRef();

        delete[]Coords;
        delete[]conn;
        coordarr->decrRef();
        VolumeMesh->decrRef();
        BoundaryMesh->decrRef();
        mm->decrRef();
        }
#else
    std::cout <<
              "\n \n MGGenCase::Print_med you don't have MED library included in your project\n\n";
#endif
    std::cout <<
              "\033[1;32m ------------------------------------------------------------------------------------------------------- \n";
    std::cout << " Mesh for interface creation printed in file " << info_name <<
              std::endl;
    std::cout << " The domain is defined on " << _n_levels << " levels \n";
    std::cout << " For each level the following maps are printed: \n\
    Proc_Lev_<level>               ->  cell-wise field showing the processor id for each cell \n\
    MG_cell_id_Lev_<level>         ->  cell-wise field the cell-id within the global cell numbering (from level 0 to highest level)\n\
    FinerLevelNodeIDS_Lev_<level>  ->  point-wise field: it's a map from local (level) node numbering to highest level node numbering\n";
    std::cout <<
              " ------------------------------------------------------------------------------------------------------- \033[0m\n";

    return;
    }


void
MGGenCase::GenLevelMed(int Level,
                       int celle,
                       int CumulativeElementsOtherLevels,
                       int *v_inv_el,
                       int *v_inv_nd,
                       int **elem_sto,
                       double *coord,
                       std::string info_name, const unsigned int LibToMed[]) {

    MEDCoupling::DataArrayDouble * Volcells2 =
        MEDCoupling::DataArrayDouble::New();
    Volcells2->alloc(celle, 1);

    MEDCoupling::DataArrayDouble * proc_num2 =
        MEDCoupling::DataArrayDouble::New();
    proc_num2->alloc(celle, 1);

    int *NodeConnLev1 = new int[celle * NDOF_FEM];
    int *NodeConnFinerLev = new int[celle * NDOF_FEM];
    double *coords2 = new double[celle * NDOF_FEM * _dim];

    for (int l = 0; l < celle; l++) {
        int reverse_element = l + CumulativeElementsOtherLevels;
        int level = elem_sto[reverse_element][NDOF_FEM + 1];
        int proc = elem_sto[reverse_element][NDOF_FEM + 2];

        std::cout << elem_sto[reverse_element][0] << std::endl;

        for (int k = 0; k < NDOF_FEM; k++) {
            NodeConnLev1[l * NDOF_FEM + k] = l * NDOF_FEM + k;
            NodeConnFinerLev[l * NDOF_FEM + k] =
                elem_sto[reverse_element][1 + LibToMed[k]];
            for (int c = 0; c < _dim; c++) {
                coords2[(l * NDOF_FEM + k) * _dim + c] =
                    coord[elem_sto[reverse_element][1 + LibToMed[k]] * _dim + c];
                }
            }
        Volcells2->setIJ(l, 0, v_inv_el[reverse_element]);
        proc_num2->setIJ(l, 0, (double) proc);
        }

    MEDCoupling::MEDCouplingUMesh * VolumeMesh2 =
        MEDCoupling::MEDCouplingUMesh::New("Mesh_Lev_" + std::to_string(Level),
                                           _dim);
    VolumeMesh2->allocateCells(celle);
    for (int i = 0; i < celle; i++) {
        VolumeMesh2->insertNextCell(MED_EL_TYPE, NDOF_FEM,
                                    NodeConnLev1 + i * NDOF_FEM);
        }
    VolumeMesh2->finishInsertingCells();
    MEDCoupling::DataArrayDouble * coordarr2 =
        MEDCoupling::DataArrayDouble::New();
    coordarr2->alloc(celle * NDOF_FEM, _dim);
    std::copy(coords2, coords2 + celle * NDOF_FEM * _dim,
              coordarr2->getPointer());
    VolumeMesh2->setCoords(coordarr2);
    std::cout << "               " << VolumeMesh2->getNumberOfNodes() << std::
              endl;

    int node_number = VolumeMesh2->getCoords()->getNumberOfTuples();
    bool areNodesMerged;
    int newNbOfNodes;
    VolumeMesh2->mergeNodes(1.e-13, areNodesMerged, newNbOfNodes);	// removing double nodes
    if (areNodesMerged) {
        std::cout << "Level " << Level <<
                  ": mesh nodes merged from initial value " << node_number << " to " <<
                  newNbOfNodes << std::endl;
        }
    VolumeMesh2->zipCoords();

    MEDCoupling::DataArrayDouble * FinerConn =
        MEDCoupling::DataArrayDouble::New();
    FinerConn->alloc(VolumeMesh2->getNumberOfNodes(), 1);
    double *MapArray = const_cast < double *>(FinerConn->getPointer());

    for (int cel = 0; cel < celle; cel++) {
        std::vector < int >CellConn;
        VolumeMesh2->getNodeIdsOfCell(cel, CellConn);
        for (int nod = 0; nod < NDOF_FEM; nod++) {
            MapArray[CellConn[nod]] =
                v_inv_nd[NodeConnFinerLev[cel * NDOF_FEM + nod]];
            }
        }

    MEDCoupling::WriteUMesh(info_name, VolumeMesh2, false);

    // CELL FIELD STORING THE PROC ID FOR EVERY MESH CELL
    MEDCoupling::MEDCouplingFieldDouble * proc_field2 =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
    proc_field2->setMesh(VolumeMesh2);
    proc_field2->setArray(proc_num2);
    proc_field2->setName("Proc_Lev_" + std::to_string(Level));

    MEDCoupling::MEDCouplingFieldDouble * FinerLevelIDS =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
    FinerLevelIDS->setMesh(VolumeMesh2);
    FinerLevelIDS->setArray(FinerConn);
    FinerLevelIDS->setName("FinerLevelNodeIDS_Lev_" + std::to_string(Level));

    // CELL FIELD STORING THE MED-TO-MG MESH CELL NUMBERING
    MEDCoupling::MEDCouplingFieldDouble * cell_field2 =
        MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
    cell_field2->setMesh(VolumeMesh2);
    cell_field2->setArray(Volcells2);
    cell_field2->setName("MG_cell_id_Lev_" + std::to_string(Level));
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, proc_field2);
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, cell_field2);
    MEDCoupling::WriteFieldUsingAlreadyWrittenMesh(info_name, FinerLevelIDS);
//     CumulativeElementsOtherLevels += celle;

    VolumeMesh2->decrRef();
    cell_field2->decrRef();
    proc_field2->decrRef();
    FinerConn->decrRef();
    delete[]coords2;
    delete[]NodeConnLev1;
    delete[]NodeConnFinerLev;
    return;
    }









// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
