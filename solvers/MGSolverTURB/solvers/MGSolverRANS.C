#include "Equations_conf.h"

// ============================================
#ifdef RANS_EQUATIONS // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverRANS.h"
#include "Printinfo_conf.h"

// local include -------------------
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
/// This function constructs the 3d-2D MGSolRANS class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
* This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
* equation name "T", basic variable name "T"
*/

MGSolRANS::MGSolRANS (
    MGEquationsSystem & mg_equations_map_in, ///<  mg_equations_map_in pointer
    const int nvars_in[],                   ///< KLQ number of variables
    std::string eqname_in,                  ///< equation name
    std::string varname_in                  ///< basic variable name
) :
    MGSolDA ( mg_equations_map_in, nvars_in, eqname_in, varname_in ),
    _offset ( _mgmesh._NoNodes[_NoLevels - 1] ),        // mesh nodes (top level)
    _dt ( stod ( _mgutils._sim_config["dt"] ) ),                  // parameter  dt
    _uref ( _mgutils._mat_prop["Uref"] ),      // parameter  u reference
    _lref ( _mgutils._mat_prop["Lref"] ),      // parameter  l reference
    _rhof ( _mgutils._mat_prop["rho0"] ),      // parameter density
    _muf ( _mgutils._mat_prop["mu0"] ),
    _Tref ( _mgutils._mat_prop["Tref"] ), // parameter  temperature reference
    _cp0 ( _mgutils._mat_prop["cp0"] ),   // parameter  Cp reference
    _kappa0 ( _mgutils._mat_prop["kappa0"] )   // parameter  conductivity reference
{
    //  =========================================================================

    /// A) reading parameters  for field coupling (in _FF_idx[])
    _RANS_parameter.read_param ( _mgutils );
    _TimeDer = _RANS_parameter._TimeDer;

    
    FillBoundaryMap();
    _nTdim = DIMENSION;
    _euler_impl = 1.;
    _IRe = _muf / ( _rhof * _uref * _lref );

    for ( int k_index = 0; k_index < 30; k_index++ ) {
        _FF_idx[k_index] = -1;
    }

    for ( int l = 0; l < _NoLevels; l++ ) {
        _solver[l]->set_solver_type ( _RANS_parameter._SolverType );
    }

    // SETTING TURBULENCE MODEL PARAMETERS -> READING FROM TurbUtils CLASS
    _WallDist = _mgutils._geometry["Wall_dist"];
    _InvSigma = 1. / 1.4;
    _SolveRANS = true;

//   _SUPG, _SCPG, _UPWIND, _RNSTAB, _RESSTAB
    std::map<std::string, bool> YesNo;
    YesNo["yes"] = true;
    YesNo["no"]  = false;


    _SolveRANS = YesNo[_mgutils._sim_config["SolveDynamicalTurbulence"]];

    _SUPG   = _RANS_parameter._Supg;
    _UPWIND = _RANS_parameter._Upwind;
    _SCPG   = _RANS_parameter._Scpg;
    _NLITER = _RANS_parameter._MaxNonLinearIterations;
    _FlatProfile = _RANS_parameter._FlatProfile;
    _SolveSteady = _RANS_parameter._SolveSteady;
    _AxiSym      = ( int ) ( _mgutils._geometry["Axysim"] );

    _InterpolatedMuTurb = _RANS_parameter._InterpolatedMuTurb;

    std::map<std::string, int> StabMap;
    StabMap["gls"]   = 1;
    StabMap["supgc"] = 0;
    StabMap["sgs"]   = -1;

    if ( _RANS_parameter._ModifiedSupg == "no" ) {
        _ModifiedSupg = false;
        _theta = 0;
    } else {
        _ModifiedSupg = true;
        _theta = StabMap[_RANS_parameter._ModifiedSupg];
    }

    _Restart = 0;


    return;
}



//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void  MGSolRANS::GenMatRhs (
    const double    /**< time (in) */,
    const int Level /**< discretization Level (in) */,
    const int mode  /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
)    // ===============================================
{
    SetTurbPar();
    /// a) Set up
    const int unsteady_flag = _RANS_parameter._SolveSteady;
    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels - 1]; // mesh nodes
    const int  el_sides = _mgmesh._GeomEl._n_sides[0];   // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];                       // bd element connectivity
    int        sur_toply[NDOF_FEMB];                     // boundary topology

    // gauss integration  -----------------------------------------------------------------------------
    double     u_old[NDOF_FEM];
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
    A[Level]->zero();

    if ( mode == 1 ) {
        b[Level]->zero();  // global matrix+rhs
    }

    _KeM.resize ( el_mat_nrows, el_mat_ncols );
    _FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

    int ndof_lev = 0;

    for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + Level + 1] - _mgmesh._off_el[0][pr * _NoLevels + Level];
        ndof_lev += delta;
    }

    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level + _NoLevels * _iproc + 1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level + _NoLevels * _iproc]; // stop element

    const int BoundLoop = ( ( _RANS_parameter._FractionalStep == 1 && _NonLinearIt == 0 ) ||  _RANS_parameter._FractionalStep == 0 ) ? 1 : 0;

    double GlobResidual = 0.;

    for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) { // LOOP OVER MESH ELEMENTS
        // set to zero matrix and rhs and center
        _KeM.zero();
        _FeM.zero();

        // ----------------------------------------------------------------------------------
        /// 1. Geometry and element  fields
        // ----------------------------------------------------------------------------------
        _mgmesh.get_el_nod_conn ( 0, Level, iel, el_conn, _xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides, 0, Level, iel, el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

        for ( int deg = 0; deg < 3; deg++ ) { // OLD SOLUTION
            for ( int eq = 0; eq < _data_eq[deg].n_eqs; eq++ ) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0, _data_eq[deg].indx_ub[eq + 1] - _data_eq[deg].indx_ub[eq],
                                                       el_ndof[deg], el_conn, offset, _data_eq[deg].indx_ub[eq], _data_eq[deg].ub );
            }
        }

        for ( int dim = 0; dim < 2; dim++ ) { // SOLUTION OF NON LINEAR ITERATIONS
            _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F + dim]]->get_el_nonl_sol ( 0, 1, el_ndof[2], el_conn, offset, dim, _tur_nl );
        }

        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F + _dir]]->get_el_sol ( 0, 1, el_ndof[2], el_conn, offset, 0, _x_1ts );
        _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[K_F + _dir]]->get_el_oldsol ( 0, 1, el_ndof[2], el_conn, offset, 0, _x_2ts );

        // ----------------------------------------------------------------------------------
        /// 2. Boundary integration  (bc)
        // ----------------------------------------------------------------------------------

        int ElemOnBoundary_flags = 0;
        for ( int k = 0; k < el_ndof[2]; k++ ) {
            int bc_node = _bc_vol[k] % 100;
            _bc_el[k] = ( bc_node / 10 == 0 ) ? 0 : 1;
            _bc_bound[k] = 1;
            if ( bc_node != 11 ) {
                ElemOnBoundary_flags=1;
            }
        }

        _BoundElem = false;
        _WallElement = 0;
        _y_dist = _mgmesh._dist[ iel + nel_b];

        for ( int dir = 0; dir < _nTdim; dir++ ) {
            _NormMidCell[dir] = 0.;
        }

        for ( int iside = 0; iside < el_sides; iside++ ) // CHECK IF THE ELEMENT HAS BOUNDARY SIDES =======
            if ( el_neigh[iside] == -1 ) {
                _BoundElem = true;
            }

        // ===========================================================================================

        int NumOfWallSides = 0;
        if ( _BoundElem ) { // INTEGRATION OVER BOUNDARY ==================================================

            for ( int j=0; j<el_ndof[2]; j++ ) {
                _NodeOnNwallSides[j]=0;
            }

            for ( int iside = 0; iside < el_sides; iside++ )  {
                if ( el_neigh[iside] == -1 ) {

                    int mid_node = _mgmesh._GeomEl._surf_top[ ( elb_ndof[2]-1 ) + NDOF_FEMB * iside];
                    const int  kdof_top = _node_dof[_NoLevels - 1][ el_conn[mid_node]]; // dof from top level
                    if ( _BoundaryMap[_mgmesh._NodeBDgroup[kdof_top]] == WALL ) {
                        NumOfWallSides ++;
                        for ( int idof = 0; idof < elb_ndof[2]; idof++ ) {
                            int pos = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside]; // local nodes
                            _NodeOnNwallSides[pos] ++;
                        }
                    }
                }
            }

            for ( int iside = 0; iside < el_sides; iside++ )  {
                if ( el_neigh[iside] == -1 ) {
                    for ( int idof = 0; idof < elb_ndof[2]; idof++ ) {
                        sur_toply[idof] = _mgmesh._GeomEl._surf_top[idof + NDOF_FEMB * iside]; // local nodes
                        int idofb = sur_toply[idof];                                    // connectivity vector
                        _bc_bound[idofb] = 0;

                        for ( int idim = 0; idim < _nTdim; idim++ ) {
                            _xxb_qnds[idim * NDOF_FEMB + idof] = _xx_qnds[idim * NDOF_FEM + idofb];  // coordinates
                        }
                    }
//                     int before = _WallElement;
                    if ( BoundLoop == 1 ) {
                        bc_set ( sur_toply, el_ndof[2], elb_ndof[2], elb_ngauss, el_conn );
                    }
//                     int after  = _WallElement;
//                     if ( after * before == 0 && after + before > 0 ) {
//                         _WallElement = ( after > before ) ? after : before;
//                     }
//                     if ( _RANS_parameter._WallFunctionApproach == 0 ) {
//                         _WallElement = 0;
//                     }
                }
            }// =========================================================================================
        }

        // Nodes on boundary but no element side on it
//       if (ElemOnBoundary_flags == 1 && !_BoundElem ){
//         for ( int k = 0; k < el_ndof[2]; k++ ) {
//           int bc_node = _bc_vol[k] % 100;
//           if(bc_node != 11) {_bc_el[k]=0;
//           std::cout<<"iel "<<iel<<" node "<<el_conn[k]<<" bc_node "<<bc_node<<std::endl;}
//           }
//       }

        // ----------------------------------------------------------------------------------
        //   3. Volume integration
        // ----------------------------------------------------------------------------------


//         _bc_el[i] != 0
        _EquationNodes.clear();
        for ( int i=0; i<el_ndof2; i++ )
            if ( _bc_el[i] != 0 ) {
                _EquationNodes.push_back ( i );
            }

        // VOLUME INTEGRATION
        if ( _WallElement == 1 && _RANS_parameter._WallFunctionApproach == 1 ) {
            wall_element ( el_ndof2, el_ngauss, _xx_qnds, iel, el_conn );
        } else {
            vol_integral ( el_ndof2, el_ngauss, _xx_qnds, iel, el_conn );
        }

        // COMPUTE RESIDUAL
        double ElemResidual = 0.;
        double RowResidual  = 0.;
        for ( int i=0; i<el_mat_nrows; i++ ) {
            RowResidual = 0.;
            if ( _bc_el[i]!=0 ) {
                for ( int j=0; j<el_mat_ncols; j++ ) {
                    RowResidual += _KeM ( i,j ) * _x_1ts[j];
                }
                RowResidual -= _FeM ( i );
            }
            ElemResidual += RowResidual*RowResidual;
        }

        GlobResidual += ElemResidual;
//       if (fabs(ElemResidual) > 1.e-3)
//         vol_stab ( el_ndof2, el_ngauss, mode, el_conn );
//
        // ----------------------------------------------------------------------------------
        //  4. add local to global
        // ----------------------------------------------------------------------------------
        A[Level]->add_matrix ( _KeM, el_dof_indices );               // global matrix

        if ( mode == 1 ) {
            b[Level]->add_vector ( _FeM, el_dof_indices );  // global rhs
        }

    }// END LOOP OVER MESH ELEMENTS

    std::cout<<"Cumulative Residual "<<sqrt ( GlobResidual ) <<std::endl;

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

    return;
}


//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void  MGSolRANS::SetUpFlags()
{
    /// a) Set up
    const int unsteady_flag = _RANS_parameter._SolveSteady;
    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels - 1]; // mesh nodes
    const int  el_sides = _mgmesh._GeomEl._n_sides[0];   // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];                       // bd element connectivity
    int        sur_toply[NDOF_FEMB];                     // boundary topology

    // gauss integration  -----------------------------------------------------------------------------
    double     u_old[NDOF_FEM];
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

    _KeM.resize ( el_mat_nrows, el_mat_ncols );
    _FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

    int ndof_lev = 0;

    for ( int pr = 0; pr < _mgmesh._iproc; pr++ ) {
        int delta = _mgmesh._off_el[0][pr * _NoLevels + _NoLevels + 1] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels];
        ndof_lev += delta;
    }

    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][_NoLevels + _NoLevels * _iproc + 1]; // start element
    const int nel_b = _mgmesh._off_el[0][_NoLevels + _NoLevels * _iproc]; // stop element

    for ( int iel = 0; iel < ( nel_e - nel_b ); iel++ ) { // LOOP OVER MESH ELEMENTS
        // set to zero matrix and rhs and center
        _KeM.zero();
        _FeM.zero();

        // ----------------------------------------------------------------------------------
        /// 1. Geometry and element  fields
        // ----------------------------------------------------------------------------------
        _mgmesh.get_el_nod_conn ( 0, _NoLevels, iel, el_conn, _xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides, 0, _NoLevels, iel, el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( _NoLevels, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

        _WallElement = 0;

        for ( int iside = 0; iside < el_sides; iside++ )  {
            if ( el_neigh[iside] == -1 ) {
                const int node_id = _mgmesh._GeomEl._surf_top[NDOF_FEMB + NDOF_FEMB * iside]; // local nodes
                const int FaceBD = _bc_vol[node_id] % 100;
                switch ( FaceBD ) {
                case 1:
                case 12:
                    _WallElement = 1;
                    break;
                default:
                    break;
                }
            }
        }

        // VOLUME INTEGRATION
        if ( _WallElement == 1 ) {
            for ( int i=0; i<el_ndof2; i++ ) {
                _bc_vol[i] = 8;
            }

            set_el_dof_bc ( _NoLevels, iel + ndof_lev, el_ndof, el_conn, offset, el_dof_indices, _bc_vol, _bc_bd );

        }

    }// END LOOP OVER MESH ELEMENTS

    return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolRANS::MGTimeStep (
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
)
{
// =========================================================================================

/// A) Set up the time step
    if ( _SolveRANS ) {
        std::cout  << std::endl << "\033[038;5;" << 217 << ";1m "
                   << "--------------------------------------------------- \n\t"
                   <<  _eqname.c_str() << "   " << _var_names[0].c_str()
                   << " solution of problem " << _mgutils.get_name() << " with dir " << _dir
                   << "\n --------------------------------------------------- \n \033[0m";

        // SET UP XOOLD AND XNONL VECTORS AFTER RESTART
        if ( _Restart == 0 ) {
            _TimeDer = 1;
            x_old[_NoLevels - 1]->localize ( *x_oold[_NoLevels - 1] );
            x_old[_NoLevels - 1]->localize ( *x_nonl[_NoLevels - 1] );
        }

        _Restart = 1;

        // EQUATION ASSEMBLY + SOLUTION
        if ( _RANS_parameter._FractionalStep == 1 ) {
            FractionalTimeStep ( time );
        }

        if ( _RANS_parameter._FractionalStep == 0 ) {
            StandardTimeStep ( time );
        }

        // UPDATE XOLD AND XOOLD VECTORS
        MGUpdateStep();

        _TimeDer = _RANS_parameter._TimeDer;
    }

    return;
}// =======================================================================================

void MGSolRANS::FractionalTimeStep (
    const double time
)
{
    /// B) Assemblying of the Matrix-Rhs

    _NonLinearIt = 1;

    std::clock_t start_time = std::clock();
    GenMatRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs
    for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
        GenMatRhs ( time, Level, 0 );  // matrix
    }
    std::clock_t end_time = std::clock();
    MGSolve ( 1.e-6, 40 );
    std::clock_t end_time2 = std::clock();

#if PRINT_TIME==1
    std::cout <<" Ass. time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC << "s "
              <<" Ass. and sol. time: =" << double ( end_time2 - start_time ) / CLOCKS_PER_SEC << "s " << std::endl;
#endif

    _NonLinearIt = 0;
    x[_NoLevels - 1]->localize ( *x_nonl[_NoLevels - 1] );
    GenMatRhs ( time, _NoLevels - 1, 1 );                                                // matrix and rhs
    for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
        GenMatRhs ( time, Level, 0 );  // matri
    }
    MGSolve ( 1.e-5, 40 );                                                               // solve

    return;
}

void MGSolRANS::StandardTimeStep (
    const double time
)
{
    /// B) Assemblying of the Matrix-Rhs
    _NonLinearIt = 0;

    std::clock_t start_time = std::clock();
    GenMatRhs ( time, _NoLevels - 1, 1 );                                         // matrix and rhs
    for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
        GenMatRhs ( time, Level, 0 );  // matrix
    }
    std::clock_t end_time = std::clock();
    MGSolve ( 1.e-6, 40 );
    std::clock_t end_time2 = std::clock();

#if PRINT_TIME==1
    std::cout <<" Ass. time -----> =" << double ( end_time - start_time ) / CLOCKS_PER_SEC << "s "
              <<" Ass. and sol. time: =" << double ( end_time2 - start_time ) / CLOCKS_PER_SEC << "s " << std::endl;
#endif

    for ( int iter = 0; iter < _NLITER; iter++ ) {
        _NonLinearIt = 1;
        x[_NoLevels - 1]->localize ( *x_nonl[_NoLevels - 1] );

        GenMatRhs ( time, _NoLevels - 1, 1 );                                                // matrix and rhs

        for ( int Level = 0 ; Level < _NoLevels - 1; Level++ ) {
            GenMatRhs ( time, Level, 0 );  // matri
        }

        MGSolve ( 1.e-5, 40 );                                                               // solve
    }

    return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolRANS::MGTimeStep_no_up (
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
)
{
// =========================================================================================

    if ( _SolveRANS ) {
        std::cout  << std::endl << "\033[038;5;" << 217 << ";1m "
                   << "--------------------------------------------------- \n\t"
                   <<  _eqname.c_str() << "   " << _var_names[0].c_str()
                   << " solution of problem " << _mgutils.get_name()
                   << "\n --------------------------------------------------- \n \033[0m";

        // SET UP XOOLD AND XNONL VECTORS AFTER RESTART
        if ( _Restart == 0 ) {
            x_old[_NoLevels - 1]->localize ( *x_oold[_NoLevels - 1] );
            x_old[_NoLevels - 1]->localize ( *x_nonl[_NoLevels - 1] );
        }

        _Restart = 1;

        // EQUATION ASSEMBLY + SOLUTION
        if ( _RANS_parameter._FractionalStep == 1 ) {
            FractionalTimeStep ( time );
        }

        if ( _RANS_parameter._FractionalStep == 0 ) {
            StandardTimeStep ( time );
        }
    }

    return;
}// =======================================================================================

void MGSolRANS::MGUpdateStep()
{
    x_old[_NoLevels - 1]->localize ( *x_oold[_NoLevels - 1] );

    int size = x[_NoLevels-1]->size();
    for ( int i=0; i<size; i++ ) {
       double k[1] = { ( * ( x[_NoLevels -1] ) ) ( i ) };
       if ( k[0] < _TurLowerLim[_dir] ) {
          k[0] = _TurLowerLim[_dir];
          x[_NoLevels-1]->set ( i, k[0] );
       }
    }

    x[_NoLevels - 1]->localize ( *x_old[_NoLevels - 1] );

    return;
}

void MGSolRANS::SetTurbPar()
{
    if ( _mgutils._TurbParameters->_IsFilled ) {
        _InvSigma = _mgutils._TurbParameters->_Nagano / 1.4 + _mgutils._TurbParameters->_Wilcox / 2.;
        _WallDist = _mgutils._TurbParameters->_BoundWallDist;
    }
}

void MGSolRANS::FillBoundaryMap()
{

    int nBound = _RANS_parameter._BoundaryGroupsIDs.size();

    for ( int i=0; i<nBound; i++ ) {
//     std::pair<int, BoundaryType> coup = std::make_pair( _RANS_parameter._BoundaryGroupsIDs[i], OTHER );
//     _BoundaryMap.insert(coup) ;
        _BoundaryMap[_RANS_parameter._BoundaryGroupsIDs[i]] = OTHER;
    }

    std::vector<int> WallGroups;

    std::string GroupString = _mgutils._sim_config["WallGroups"];
    std::string temps;
    int Length = GroupString.length();
    int count, pos1;
    count = pos1 = 0;

    while ( count < Length ) {
        if ( GroupString.at ( count ) == ',' ) {
            temps = GroupString.substr ( pos1, count - pos1 );
            WallGroups.push_back ( stoi ( temps ) );
            pos1 = count + 1;
        }
        count++;
    }
    WallGroups.push_back ( stod ( GroupString.substr ( pos1, Length - pos1 ) ) );

    int nWalls = WallGroups.size();
    std::cout<<"Num Wall Sides "<<nWalls<<" string sides "<<GroupString<<std::endl;
    for ( int i=0; i<nWalls; i++ ) {
        _BoundaryMap[WallGroups[i]] = WALL;
    }

    for ( int i=0; i<nBound; i++ ) {
        std::cout<<"Group Num: "<<_RANS_parameter._BoundaryGroupsIDs[i]<<" Bound Type "<<_BoundaryMap[_RANS_parameter._BoundaryGroupsIDs[i]]<<std::endl;
    }

    return;
}

#endif
// #endif // personal application

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
