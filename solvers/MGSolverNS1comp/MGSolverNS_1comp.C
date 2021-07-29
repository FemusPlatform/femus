// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#include "MGSolverNS_1comp.h"       // Navier-Stokes class header file

#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "EquationSystemsExtendedM.h"  // Equation map class

// local alg lib -----------------------------------------------
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers



// ==================================================================
/// This routine constructs the FSI class:
MGSolNS_1comp::MGSolNS_1comp (
    MGEquationsSystem &mg_equations_map_in,
    int             nvars_in[],
    std::string     eqname_in,
    std::string     varname_in
) :  MGSolDA ( mg_equations_map_in,nvars_in,eqname_in,varname_in ),
//
//===================================================================================================//
//                                    A) Reading parameters                                          //
//===================================================================================================//
//
_offset ( _mgmesh._NoNodes[-1] ),          // mesh nodes (top level)
_dt ( stod ( _mgutils._sim_config["dt"] ) ),                  // parameter  dt
_uref ( _mgutils._mat_prop["Uref"] ),      // parameter  u reference
_lref ( _mgutils._mat_prop["Lref"] ),      // parameter  l reference
_rhof ( _mgutils._mat_prop["rho0"] ),      // parameter density
_muf ( _mgutils._mat_prop["mu0"] ) {      // parameter viscosity
    //===================================================================================================//
    //                                    B) Setting class variables                                     //
    //===================================================================================================//

    // READ PARAMETERS FROM CLASS NS_parameter
    _NS_parameter.read_param ( _mgutils );    
//     _NumRestartSol = _NS_parameter._NumRestartSol;
    _nNSdim=DIMENSION;
    // class equation ---------------------------------------------------------------------
    for ( int k_index=0; k_index<40; k_index++ ) {
        _FF_idx[k_index]=-1;
    }
    _pres_order= ( _nvars[0]>0 ) ? 0:1;
    // class variable names ------------------------------------------------------------
    _dir=0;
    _var_names[0]="u";
    _refvalue[0]=_uref; // velocity 3D

    //===================================================================================================//
    //                                    C) Setting solver type                                         //
    //===================================================================================================//
    for ( int l=0; l<_NoLevels; l++ ) {
        _solver[l]->set_solver_type ( _NS_parameter._SolverType );
    }
    //===================================================================================================//
    //                                    D) Setting non_dimensional parameters  01                        //
    //===================================================================================================//

    _IRe=_muf/ ( _rhof*_lref*_uref );       // Reynolds number
    _IFr=9.81*_lref/ ( _uref*_uref );       // Froud number
    _dirg[0] = _mgutils._mat_prop["dirgx"];      // x-gravity
    _dirg[1] = _mgutils._mat_prop["dirgy"];      // y-gravity
    _dirg[2] = _mgutils._mat_prop["dirgz"];      // z-gravity
    for ( int i=0; i<_nNSdim; i++ ) {
        _xyzg[i]=0.;
    }
    _y_bcout=sqrt ( 2*0.09*_IRe/80. );

    std::map<std::string, bool> YesNo;
    YesNo["yes"] = true;
    YesNo["no"]  = false;

    _SolveNS  = YesNo[_mgutils._sim_config["SolveNavierStokes"]];
    _Wall_dist =_mgutils._geometry["Wall_dist"];
    _AxiSym  = ( int ) ( _mgutils._geometry["Axisym"] );

    if ( _AxiSym == 1 && _nNSdim == 3 ) {
        std::cout<<"\033[1;31m MGSolNS: CANNOT USE AXISYM OPTION WITH 3D GEOMETRY \n \033[0m";
        std::abort();
    }

    return;
} //****************************************************************************************************//
//
//

// #include "MGsolverNS_vol.h"01

// ====================================================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================================================
void  MGSolNS_1comp::GenMatRhs ( 
const double/* time*/, 
const int Level,
const  int mode )
{
    // ===================================================================================================
    //                                    A) Set up
    // ===================================================================================================

    // NS parameters

    const double les=_NS_parameter._Les;
    // -------------------------- Geometry -------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
    int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
    int        el_neigh[NDOF_FEM];                                         // element connectivity
    int        sur_toply[NDOF_FEMB];                                       // boundary topology
    double     normal[DIMENSION];                                          // normal to the boundary
    // -------------------------- Gauss integration -----------------------------------------------------

    const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points

    // Number of  element dof: constant[0]-linear[1]-quadratic[2] ---------------------------------------
    int el_ndof[3];
    el_ndof[0]=NDOF_K;
    int elb_ndof[3];
    elb_ndof[0]=0; // number of el dofs
    int el_mat_nrows =0;                                            // number of mat rows (dofs)
    for ( int ideg=1; ideg<3; ideg++ ) {                            //     ...
        el_ndof[ideg]= ( ( _nvars[ideg]>0 ) ?    _fe[ideg]->_NoShape[ _nNSdim-1]:0 );              //   computing
        elb_ndof[ideg]= ( ( _nvars[ideg]>0 ) ?_fe[ideg]->_NoShape[ _nNSdim-2]:0 );            //     ...
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    }

    el_mat_nrows +=  el_ndof[0]*_nvars[0];
    int el_mat_ncols = el_mat_nrows;                                //square matrix
    std::vector<int> el_dof_indices ( el_mat_ncols );               // element dof vector

    // fields -> Navier-Stokes ----------------------------------------------------------------------
//   double u_nlg[DIMENSION];
    double u_old[DIMENSION*NDOF_FEM];
    double u_oold[DIMENSION*NDOF_FEM];
    double u_nl[DIMENSION*NDOF_FEM];       // velocity vector for non linear terms -> it contains all the velocity components //
    double x_m[DIMENSION];
    int WallNodes[NDOF_FEMB];


    for ( int k=0; k<40; k++ ) { // coupling  basic system fields
        const int idx= _data_eq[2].tab_eqs[k];
        _FF_idx[k]= ( idx>=0 ) ?_data_eq[2].indx_ub[idx]:-1;
    }
    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
    A[Level]->zero();
    if ( mode ==1 ) {
        b[Level]->zero();     // global matrix+rhs
    }
    DenseMatrixM KeM;
    DenseVectorM FeM;                  // local  matrix+rhs
    KeM.resize ( el_mat_nrows,el_mat_ncols );
    FeM.resize ( el_mat_nrows );       // resize  local  matrix+rhs

    int ndof_lev=0;   // for multilevel multi proc point ordering ------------------------------------------
    for ( int pr=0; pr <_mgmesh._iproc; pr++ ) {
        ndof_lev +=_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    }

    // ===================================================================================================
    //                                    B) Element  Loop over the volume (n_elem)
    // ===================================================================================================

    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    
    for ( int iel=0; iel < ( nel_e - nel_b ); iel++ ) {
        // set to zero matrix and rhs and center01
        KeM.zero();
        FeM.zero();

        // ---------------------------------------------------------------------------------
        // Volume  geometry and element  fields
        // ---------------------------------------------------------------------------------
        // Element Connectivity (el_conn)  and coordinates (xx_qnds)
        _mgmesh.get_el_nod_conn ( 0,Level,iel,el_conn,_xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides,0,Level,iel,el_neigh );
        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd );
        // field data  ------------------------------------------------------
        get_el_field_data ( iel, Level,el_conn, offset,el_ndof,ndof_lev, u_old,u_oold, u_nl );

        // initializing  volume quantities
        for ( int idim=0; idim< 1; idim++ ) { // quad loop entities (vector)
            x_m[idim]=0.;
            for ( int d=0; d<el_ndof[2]; d++ ) {
                const int  dnode=idim*el_ndof[2]+d;    // index local points
                x_m[idim] +=_xx_qnds[dnode]/el_ndof[2];
                _normal_pt[dnode]=0.;// normal boundary
                _bc_el[dnode]=0;
                _AlreadyWrittenDirBC[dnode] = false;
            }
        }

        // SETTING _bc_el = -8 FOR EVERY ELEMENT BOUNDARY NODE
        for ( int  iside=0; iside< el_sides; iside++ ) {
            if ( el_neigh[iside] == -1 ) {
                for ( int  lbnode=0; lbnode<elb_ndof[2]; lbnode++ ) {
                    int lnode=_mgmesh._GeomEl._surf_top[lbnode+elb_ndof[2]*iside];// local nodes
                    for ( int idim=0; idim< 1; idim++ ) {
                        _bc_el[lnode + idim*el_ndof[2]] = _BdFlagId;
                    }
                }
            }
        }
        // element fields ----------------------------------
        if ( _FF_idx[K_F]>=0 ) 
            _y_bcout= _mgmesh._dist[ iel+nel_b];     //  _y_bcout=_mgmesh._dist[ iel+nel_b];
        
        _WallElement=0;
        // --------------------------------------------------------------------------
        //  Boundary and boundary bc
        // --------------------------------------------------------------------------
        //  New bc flags for Navier Stokes equation:
        //  _bc_el = -1 -2 -3 -> Neumann bc along tangential 1, tangential 2 (only 3d), normal directions
        //  _bc_el =  1  2  3 -> Dirichlet bc along tangential 1, tangential 2 (only 3d), normal directions
        //  _bc_el =  0 -> normal interior point
        //  Integration is performed for _bc_el <=0
        //  Test functions are set to zero when _bc_el >0
        for ( int  iside=0; iside< el_sides; iside++ ) {
            if ( el_neigh[iside] == -1 ) {
                // setup boundary element  ----------------------------------------------------------------
                for ( int  lbnode=0; lbnode<elb_ndof[2]; lbnode++ ) { // quad quantities
                    int lnode=_mgmesh._GeomEl._surf_top[lbnode+elb_ndof[2]*iside];// local nodes
                    sur_toply[lbnode]=lnode;          // lbnode -> lnode
                    elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
                    for ( int idim=0; idim< _nNSdim; idim++ ) { // coordinates
                        _xxb_qnds[idim*elb_ndof[2]+lbnode]=_xx_qnds[idim*el_ndof[2]+lnode];
                    }
                }
                _fe[2]->normal_g ( _xxb_qnds,x_m,normal );
                int dir_maxnormal = ( fabs ( normal[0] ) >fabs ( normal[1] ) ) ?0:1 ;
                dir_maxnormal= ( fabs ( normal[dir_maxnormal] ) >fabs ( normal[DIMENSION-1] ) ) ? dir_maxnormal:DIMENSION-1;
                //         std::cout<<" .........................................................\n Element "<<iel<<" side "<< iside <<std::endl;
                int before  = _WallElement;
                set_bc_matrix ( KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof,elb_ngauss,normal,u_old, el_conn );
                int after  = _WallElement;
                if ( after * before == 0 && after + before == 1 )   for ( int idof=0; idof<elb_ndof[2]; idof++ ) {
                        WallNodes[idof] = sur_toply[idof];
                    }
//                 std::cout << " bc after : " << iel <<"\t"<<iside<< " x "<< x_m[0] << " y" << x_m[1]<<std::endl;
//             std::cout <<  "  KeM \n" << std::setprecision(10)<<KeM << " rhs \n" <<  FeM <<"\n\n" <<std::endl;
            } // iside -1
        }  // -----------------------------  End Boundary -------------------------------------
        
        matrixrhsvol ( KeM, FeM, el_ndof, u_old, u_nl, 0, mode, el_conn );

        // ----------------------------------------------------------------------------------
        //   E) Add them to the global matrix
        A[Level]->add_matrix ( KeM,el_dof_indices );
        if ( mode == 1 ) 
            b[Level]->add_vector ( FeM,el_dof_indices );
        
    } //  =============== End of element loop =============================================

// ===============================================================================
    el_dof_indices.clear();
    A[Level]->close();
    if ( mode == 1 ) 
        b[Level]->close();
    
// ----------------------------------------------------------
#ifdef PRINT_INFO
    //   A[Level]->print();  b[Level]->print();
    std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
} //****************************************************************************************************//
// 
// 
// // =========================================================================================
// /// This function controls the assembly and the solution of the NS_equation system:
void MGSolNS_1comp::MGTimeStep (
    const double time,  // time
    const int max_iter  // Number of max inter
) {
// =========================================================================================

    // ========================================================================================= //
    //              A) Set up the time step                                                      //
    // ========================================================================================= //
    if ( _SolveNS ) {
        std::cout  << std::endl << "\033[038;5;"<<NS_F + 50<<";1m "
        << "--------------------------------------------------- \n\t"
        <<  _eqname.c_str()
        << " solution of problem " << _mgutils.get_name()
        << "\n ---------------------------------------------------\n\033[0m";

        x[_NoLevels-1]->localize ( *x_nonl[_NoLevels-1] );
        // ========================================================================================= //
        //              B) Assemblying of the Matrix-Rhs                                             //
        // ========================================================================================= //

#if PRINT_TIME==1
        std::clock_t start_time=std::clock();
#endif
        GenMatRhs ( time,_NoLevels-1,1 );                                           // matrix and rhs
        for ( int Level = 0 ; Level < _NoLevels-1; Level++ ) {
            GenMatRhs ( time,Level,0 );     // matrix
        }
#if PRINT_TIME==1
        std::clock_t end_time=std::clock();
        std::cout << "  Assembly time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
        // ========================================================================================= //
        //              C) Solution of the linear MGsystem (MGSolFSI::MGSolve)                       //
        // ========================================================================================= //
        MGSolve ( 1.e-6,15 );
#if PRINT_TIME==1
        end_time=std::clock();
        std::cout << " Assembly+solution time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC
        << "s "<< std::endl;
#endif

        // ========================================================================================= //
        //              D) Non Linear Iterations                                                     //
        // ========================================================================================= //

        int iter_nl=_NS_parameter._MaxNonLinearIterations;
        for ( int iter=0; iter<iter_nl; iter++ ) {
            x[_NoLevels-1]->localize ( *x_nonl[_NoLevels-1] );

            std::cout  << std::endl << " === NON LINEAR ITERATOR: "<< iter+1 << "- " << _eqname.c_str() << " solution "  << std::endl;

            GenMatRhs ( time,_NoLevels-1,1 );                                             // matrix and rhs
            for ( int Level = 0 ; Level < _NoLevels-1; Level++ ) {
                GenMatRhs ( time,Level,0 );     // matrix
            }
            MGSolve ( 1.e-6,15 );                                                         // solve
//
            // ----------- Check error --------------------------------------------------------------- //

        }
        x_old[0][_NoLevels-1]->localize ( *x_old[1][_NoLevels-1] );
        x[_NoLevels-1]->localize ( *x_old[0][_NoLevels-1] );

    }
    return;
} //****************************************************************************************************//

void MGSolNS_1comp::SetBCFlags (
    double normal[], // normal unit
    int sur_toply[], // surface connectivity (respect to volume)
    int el_ndof,     // number of dofs
    int el_conn[]    // volume connectivity
) {// -----------------------------------------------------------------------------------------------------
    // bc conditions
    // Dirichlet >5  (surface point setting ( no volume integration))
    // 6-> set zero component 7-> set component       ------------------------> _bc_el[] = 3,0,-3; _bc_el[] =1,2,-1,-2;
    // 8-> zero along dir (normal or tg) 9-> value along dir  ----------------> _bc_el[] = 3,0,-3; _bc_el[] =1,2,-1,-2;
    //  Neumann (volume integration (+surface integration))   ---------------->
    // 1-> no surface integration                             ---------------->
    // 2-> pressure integration (normal only)  3-> pressure setting    ------->
    // 4-> tau=a.velocity                      4-> stau=a.(velocity-velocity_0)------>
    const int bd_face = abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 100; // element normal-tg bc flag
    const int norm_face = ( abs ( _bc_bd[sur_toply[NDOF_FEMB - 1]] ) % 10000 ) / 1000 - 1 ; // element max-normal-dir bc flag

    for ( int  lbnode = 0; lbnode < NDOF_FEMB; lbnode++ ) { // loop on surface nodes
        const int bd_node = abs ( _bc_bd[sur_toply[lbnode]] ) % 100; // point normal-tg  bc flag
        const int bc_var_check = abs ( _bc_bd[sur_toply[lbnode]] ) % 10000; // bc_var
        int bc_var_normal = bc_var_check / 1000 - 1; // point max-normal-dir  bc flag
        const int bc_n_flag = bd_node / 10;         // normal bc flag
        const int bc_tg_flag = bd_node % 10;         // tg bc flag

        bool Norm_Dirichlet_Bound_Cond = ( bc_n_flag > 5 ) ? true : false;
        bool Tang_Dirichlet_Bound_Cond = ( bc_tg_flag > 5 ) ? true : false;
        int NormDir = bc_var_normal;

        // BC are imposed only if node bc equal to face bc
        if ( bd_face == bd_node || ( bd_node == 88 || bd_node == 66 ) ) {
            int NormRow = sur_toply[lbnode] + NormDir * el_ndof;

            // Setting bc flags for normal direction --------------------------
            if ( _bc_el[NormRow] == _BdFlagId ) { // only _BdFlagId=-8 (first time)
                if ( Norm_Dirichlet_Bound_Cond ) {
                    _bc_el[NormRow] = 3;
                    }
                else {// Neumann bc (surface integration)
                    if ( bc_n_flag == 1 || bd_node == 31 ) {
                        _bc_el[NormRow] = 0;
                        }
                    else {
                        _bc_el[NormRow] = -3;
                        }
                    }
                }//  norm  ----------------------------------------------------------

            // Setting bc flags for tangential directions -------------------
            int dir_tg0 = 0, NeuTg = 0;

            for ( int kdir = NormDir + 1; kdir < NormDir + _nNSdim; kdir ++ ) {
                int idir = kdir % _nNSdim;
                int tg_dirRow = sur_toply[lbnode] + idir * el_ndof;

                if ( _bc_el[tg_dirRow] == _BdFlagId ) { // only _BdFlagId=-8 (first time)
                    if ( Tang_Dirichlet_Bound_Cond ) {
                        dir_tg0 ++;  // Dirichlet bc
                        _bc_el[tg_dirRow] = dir_tg0;
                        }
                    else { // Neumann bc (surface integration)
                        if ( bd_node == 11 || bd_node == 31 ) {
                            _bc_el[tg_dirRow] = 0; // for outflow bc we don't project the equation along normal-tangential directions
                            }
                        else {
                            NeuTg++;
                            _bc_el[tg_dirRow] = -NeuTg;
                            }
                        }
                    }
                }//  Tan ----------------------------------------------------------
            }// if (bd_face==bd_node)
        }

    return;
    }
    


/******************************************************************************************************/
#endif
// #endif  // NS_equation is personal

