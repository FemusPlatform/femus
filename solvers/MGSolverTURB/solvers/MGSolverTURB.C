#include "Equations_conf.h"


// ============================================
#ifdef _TURBULENCE_ // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverTURB.h"
#include "Printinfo_conf.h"
#include <sstream>
#include "MGGeomEl.h"
#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"

#include "MGFE.h"
#include "MGUtils.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"
#include "parallelM.h"



// ======================================================
/// This function constructs the 3d-2D MGSolTURB class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
* This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
* equation name "T", basic variable name "T"
*/
MGSolTURB::MGSolTURB (
    MGEquationsSystem & mg_equations_map_in, ///<  mg_equations_map_in pointer
    const int nvars_in[],                   ///< KLQ number of variables
    std::string eqname_in,                  ///< equation name
    std::string varname_in                  ///< basic variable name
) :
    MGSolDA ( mg_equations_map_in, nvars_in, eqname_in, varname_in ),
    _offset ( _mgmesh._NoNodes[_NoLevels - 1] ), // mesh nodes
    _uref ( _mgutils._mat_prop["Uref"] ),      // parameter  u reference
    _lref ( _mgutils._mat_prop["Lref"] ),      // parameter  l reference
    _Tref ( _mgutils._mat_prop["Tref"] ),      // parameter  l reference
    _rhof ( _mgutils._mat_prop["rho0"] ),      // parameter density
    _muf ( _mgutils._mat_prop["mu0"] )   // parameter  viscosity reference
{
    //  =========================================================================

    _var_names[0] = varname_in;

    if ( !varname_in.compare ( "dist" ) ) {
        _dir = 0;  // kappa
        _refvalue[0] = _uref * _uref;
    }

    if ( !varname_in.compare ( "muT" ) ) {
        _dir = 1;  // omega
        _refvalue[0] = _uref * _uref * _uref / _lref;
    }

    if ( !varname_in.compare ( "alphaT" ) ) {
        _dir = 2;  // omega
        _refvalue[0] = _uref * _uref * _uref / _lref;
    }

    _FirstAssembly = 1;
    _AssembleOnce  = 0;

//      if(_dir==1) _AssembleOnce = 1;

    std::cout << _dir << "\n";
    /// A) reading parameters  for field coupling (in _FF_idx[])
    _nTdim = DIMENSION;

    for ( int k_index = 0; k_index < 30; k_index++ ) {
        _FF_idx[k_index] = -1;
    }

    /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
    _var_names[0] = varname_in;
    _refvalue[0] = _Tref;

    /// C ) Setting the  solver type (with _solver[l]->set_solver_type(SOLVERT))
    for ( int l = 0; l < _NoLevels; l++ ) {
        _solver[l]->set_solver_type ( GMRESM );
    }

    _IRe = _muf / ( _rhof * _uref * _lref );

    std::string GroupString = _mgutils._sim_config["WallGroups"];
    std::string temps;
    int Length = GroupString.length();
    int count, pos1;
    count = pos1 = 0;

    while ( count < Length ) {
        if ( GroupString.at ( count ) == ',' ) {
            temps = GroupString.substr ( pos1, count - pos1 );
            _WallGroupsIDs.push_back ( stoi ( temps ) );
            pos1 = count + 1;
        }

        count++;
    }

    _WallGroupsIDs.push_back ( stoi ( GroupString.substr ( pos1, Length - pos1 ) ) );
    _FirstCall = 1;

    _DiffusionCoefficient = 1.e-6;
    if ( _dir == 0 ) {
        _DiffusionCoefficient = 1.e-8;
    }

    for ( int i = 0; i < 2 * NDOF_FEM; i++ ) {
        _KW[i] = _KHWH[i] = 1.;
    }

    if ( _dir == 0 ) {
        _Solve = 1;
    } else if ( _dir == 1 ) {
        _Solve = _mgutils._TurbParameters->_SolveMuT;
    } else if ( _dir == 2 ) {
        _Solve = _mgutils._TurbParameters->_SolveAlphaT;
    }

    return;
}



//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void  MGSolTURB::GenMatRhs (
    const double    /**< time (in) */,
    const int Level /**< discretization Level (in) */,
    const int mode  /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
)    // ===============================================
{
//   double Crank_Nicolson =1.;
    /// a) Set up
    const int unsteady_flag = _C_parameter.UNSTEADY;
    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels - 1]; // mesh nodes
    const int  el_sides = _mgmesh._GeomEl._n_sides[0];   // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];                       // bd element connectivity
    int        sur_toply[NDOF_FEMB];                     // boundary topology

    // gauss integration  -----------------------------------------------------------------------------


    const int el_ngauss = _fe[2]->_NoGauss1[ _nTdim - 1];                 // elem gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[ _nTdim - 2];           // bd elem gauss points

    // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
    int el_ndof[3];
    el_ndof[0] = 1;
    int elb_ndof[3];
    elb_ndof[0] = 1; // number of el dofs
    int el_mat_nrows = 0;                                              // number of mat rows (dofs)

    for ( int ideg = 1; ideg < 3; ideg++ ) {
        el_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 1];
        elb_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 2];
        el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
    };

    const int el_ndof2 = _fe[2]->_NoShape[_nTdim - 1];

    int el_mat_ncols = el_mat_nrows;   // square matrix

    std::vector<int> el_dof_indices ( el_mat_ncols );   // element dof vector

    // coupling  fields -------------------------------------------------------------------------------
    for ( int k = 0; k < 30; k++ ) { // coupling  basic system fields
        const int idx = _data_eq[2].tab_eqs[k];
        _FF_idx[k] = ( idx >= 0 ) ? _data_eq[2].indx_ub[idx] : -1;
    }

    double vel_g[DIMENSION];

    for ( int idim = 0; idim < _nTdim; idim++ ) {
        vel_g[idim] = 0.;    // velocity not coupled
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
    A[Level]->zero();

    if ( mode == 1 ) {
        b[Level]->zero();     // global matrix+rhs
    }

    DenseMatrixM KeM;
    DenseVectorM FeM;                // local  matrix+rhs
    KeM.resize ( el_mat_nrows, el_mat_ncols );
    FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

    int ndof_lev = 0;

    for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
        ndof_lev += delta;
    }


    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

    _res = 0.;
    
    for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {
        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        /// 1. geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
        _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, _xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );
        // fill the node data vector

        if ( _FF_idx[K_F] >= 0 ) {
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _KW ); // pressure
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F] + 1]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 1, _KW ); // pressure
        }

        if ( _FF_idx[KTT_F] >= 0 ) {
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[KTT_F]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _KHWH ); // pressure
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[KTT_F] + 1]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 1, _KHWH ); // pressure
        }

        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DIST]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _DIST ); // pressure

        for ( int deg = 0; deg < 3; deg++ ) { // OLD SOLUTION
            for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                                       el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
            }
        }

        if ( _Solve == 0 ) {

            double kw[2], khwh[2];

            for ( int i = 0; i < el_ndof2; i++ ) {

                const int  kdof_top = _node_dof[_NoLevels - 1][ el_conn[i]]; // dof from top level
                kw[0] = _KW[i];
                kw[1] = _KW[i + el_ndof2];

                if ( _dir == 1 )  {
                    double mu_turb[1];
                    mu_turb[0] = _mgutils._TurbParameters->CalcMuTurb ( kw, _DIST[i] );
                    if ( _bc_vol[i] == 0 ) {
                        mu_turb[0] = 0.;
                    }

                    x[_NoLevels - 1]->set ( kdof_top, mu_turb[0] );
                }

                if ( _dir == 2 )  {
                    khwh[0] = _KHWH[i];
                    khwh[1] = _KHWH[i + el_ndof2];
                    double alpha_turb[1];
                    alpha_turb[0] = _mgutils._TurbParameters->CalcAlphaTurb ( kw, khwh, _DIST[i] );
                    if ( _bc_vol[i] == 0 ) {
                        alpha_turb[0] = 0.;
                    }
                    x[_NoLevels - 1]->set ( kdof_top, alpha_turb[0] );
                }
            }
        }

        if ( _Solve == 1 ) {
            // ----------------------------------------------------------------------------------
            /// 2. Boundary integration  (bc)
            // ----------------------------------------------------------------------------------

            for ( int i = 0; i < NDOF_FEM; i++ ) {
                _bc_el[i] = 1;
            }

            for ( int iside=0; iside< el_sides; iside++ )  {
                if ( el_neigh[iside] == -1 ) {
                    for ( int idof=0; idof<elb_ndof[2]; idof++ ) {
                        sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                        int idofb=sur_toply[idof];
                        if ( _bc_bd[idofb]%2 == 0 ) {
                            _bc_el[idofb] = 0;
                        }
//                          for ( int idim=0; idim< _nTdim; idim++ ) _xxb_qnds[idim*NDOF_FEMB+idof]=_xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                    }
//                     bc_set(KeM, FeM, sur_toply, el_ndof2, elb_ndof[2], elb_ngauss,0);
                }
            }
            // ----------------------------------------------------------------------------------
            //   3. Volume integration
            // ----------------------------------------------------------------------------------
            //  external cell properties -------------------------------------

            // volume integral
            _y_dist = _mgmesh._dist[ iel + nel_b];
            vol_integral ( KeM, FeM, el_ndof2, el_ngauss, _xx_qnds, unsteady_flag, mode, el_conn );

            // ----------------------------------------------------------------------------------
            //  4. add local to global
            // ----------------------------------------------------------------------------------
            A[Level]->add_matrix ( KeM, el_dof_indices );              // global matrix

            if ( mode == 1 ) {
                b[Level]->add_vector ( FeM, el_dof_indices );    // global rhs
            }

        }
        
        if(_dir==1){
            int WallEl = 0;
            for ( int iside = 0; iside < el_sides; iside++ )  {
                if ( el_neigh[iside] == -1 ) {
                        int nodeid = _mgmesh._GeomEl._surf_top[elb_ndof[2] - 1 + NDOF_FEMB * iside]; // local nodes
                        if ( _bc_bd[nodeid] % 2 == 0 ) {
                            _bc_el[nodeid] = 0;
                            WallEl = 1;
                        }
                }
            }
            if(WallEl == 1){
                compute_y_plus( el_ndof2, el_ngauss, _xx_qnds, unsteady_flag, mode, el_conn, iel+nel_b );               
            }                       
        }
        

    } // end of element loop
    
    /// 5. clean
    el_dof_indices.clear();
    A[Level]->close();

    if ( mode == 1 ) {
        b[Level]->close();
    }

    //   A[Level]->print(); b[Level]->print();
#ifdef PRINT_INFO
    std::cout << " Matrix Assembled(T)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

    if(_dir==1) {
        int nProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
        int *shifts = new int[nProcs];
        int *counts = new int[nProcs];
        
        for(int np=0; np<nProcs; np++){
            const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * np + 1]; // start element
            const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * np]; // stop element
            shifts[np] = nel_b;
            counts[np] = nel_e - nel_b;
        }
//         MPI_Gather(_mgmesh._yplus + nel_b, nel_e-nel_b, MPI_DOUBLE, _mgmesh._yplus, nel_e-nel_b, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Gathering yplus values on proc 0 for printing 
        MPI_Gatherv(_mgmesh._yplus + nel_b, nel_e-nel_b, MPI_DOUBLE, _mgmesh._yplus, counts, shifts, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        delete [] shifts;
        delete [] counts;
    }
    
    return;
}

//  ===============================================================================================
/// This function assembles the rhs:
//  ===============================================================================================
void  MGSolTURB::GenRhs (
    const double    /**< time (in) */,
    const int Level /**< discretization Level (in) */,
    const int mode  /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
)    // ===============================================
{
//   double Crank_Nicolson =1.;
    /// a) Set up
    const int unsteady_flag = _C_parameter.UNSTEADY;
    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels - 1]; // mesh nodes
    const int  el_sides = _mgmesh._GeomEl._n_sides[0];   // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];                       // bd element connectivity
    int        sur_toply[NDOF_FEMB];                     // boundary topology

    // gauss integration  -----------------------------------------------------------------------------


    const int el_ngauss = _fe[2]->_NoGauss1[ _nTdim - 1];                 // elem gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[ _nTdim - 2];           // bd elem gauss points

    // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
    int el_ndof[3];
    el_ndof[0] = 1;
    int elb_ndof[3];
    elb_ndof[0] = 1; // number of el dofs
    int el_mat_nrows = 0;                                              // number of mat rows (dofs)

    for ( int ideg = 1; ideg < 3; ideg++ ) {
        el_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 1];
        elb_ndof[ideg] = _fe[ideg]->_NoShape[_nTdim - 2];
        el_mat_nrows += _nvars[ideg] * el_ndof[ideg];
    };

    const int el_ndof2 = _fe[2]->_NoShape[_nTdim - 1];

    int el_mat_ncols = el_mat_nrows;   // square matrix

    std::vector<int> el_dof_indices ( el_mat_ncols );   // element dof vector

    // coupling  fields -------------------------------------------------------------------------------
    for ( int k = 0; k < 30; k++ ) { // coupling  basic system fields
        const int idx = _data_eq[2].tab_eqs[k];
        _FF_idx[k] = ( idx >= 0 ) ? _data_eq[2].indx_ub[idx] : -1;
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
    if ( mode == 1 ) {
        b[Level]->zero();     // global matrix+rhs
    }

    DenseVectorM FeM;                // local  matrix+rhs
    FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

    int ndof_lev = 0;

    for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
        ndof_lev += delta;
    }

    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

    for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) {

        // set to zero matrix and rhs and center
        FeM.zero();

        /// 1. geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
        _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, _xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

        // fill the node data vectors
        if ( _FF_idx[K_F] >= 0 ) {
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _KW ); // pressure
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F] + 1]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 1, _KW ); // pressure
        }

        if ( _FF_idx[KTT_F] >= 0 ) {
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[KTT_F]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _KHWH ); // pressure
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[KTT_F] + 1]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 1, _KHWH ); // pressure
        }

        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[DIST]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _DIST ); // pressure

        // ----------------------------------------------------------------------------------
        /// 2. Boundary integration  (bc)
        // ----------------------------------------------------------------------------------

        for ( int i = 0; i < NDOF_FEM; i++ ) {
            _bc_el[i] = 1;
        }

        for ( int iside = 0; iside < el_sides; iside++ )  {
            if ( el_neigh[iside] == -1 ) {
                for ( int idof = 0; idof < elb_ndof[2]; idof++ ) {
                    sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside]; // local nodes
                    int idofb = sur_toply[idof];

                    if ( _bc_bd[idofb] % 2 == 0 ) {
                        _bc_el[idofb] = 0;
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------------
        //   3. Volume integration
        // ----------------------------------------------------------------------------------
        _y_dist = _mgmesh._dist[ iel + nel_b];
        rhs_integral ( FeM, el_ndof2, el_ngauss, _xx_qnds, unsteady_flag, mode, el_conn );

        // ----------------------------------------------------------------------------------
        //  4. add local to global
        // ----------------------------------------------------------------------------------
        if ( mode == 1 ) {
            b[Level]->add_vector ( FeM, el_dof_indices );    // global rhs
        }

    } // end of element loop

    /// 5. clean
    el_dof_indices.clear();

    if ( mode == 1 ) {
        b[Level]->close();
    }

#ifdef PRINT_INFO
    std::cout << " Matrix Assembled(T)  for  Level " << Level << " dofs " << A[Level]->n() << "\n";
#endif

    return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolTURB::MGTimeStep (
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
)
{
// =========================================================================================

/// A) Set up the time step
    std::cout  << std::endl << "\033[038;5;" << DIST + 50 << ";1m "
               << "--------------------------------------------------- \n\t"
               <<  _eqname.c_str()
               << " solution of problem " << _mgutils.get_name()
               << "\n ---------------------------------------------------\n\033[0m";



    /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
    std::clock_t start_time = std::clock();
#endif

    _res = 0.;

    if ( _FirstAssembly == 1 || _AssembleOnce == 0 ) {
        GenMatRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs

        for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
            GenMatRhs ( time, Level, 0 );   // matrix
        }
    } else {
        GenRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs

        for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
            GenRhs ( time, Level, 0 );   // matrix
        }
    }

#if PRINT_TIME==1
    std::clock_t end_time = std::clock();
    std::cout << "  Assembly time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
    /// C) Solution of the linear MGsystem (MGSolTURB::MGSolve).

    if ( _Solve == 1 ) {
        MGSolve ( 1.e-6, 40 );
    }

#if PRINT_TIME==1
    end_time = std::clock();
    std::cout << " Assembly+solution time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC
              << "s " << std::endl;
#endif
    /// D) Update of the old solution at the top Level  (MGSolTURB::OldSol_update),
    x[_NoLevels - 1]->localize ( *_x_olds[_NoLevels - 1][0] );


    _FirstAssembly = 0;


    return;
}// =======================================================================================


#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
