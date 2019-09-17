// class header
#include "FEMUS.h"

// forwarded classes
#include "MeshExtended.h"
#include "MGFemusInit.h"
#include "EquationSystemsExtendedM.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "MGTimeLoop.h"
#include "EquationsMap.h"
// additional local includes
#include "Printinfo_conf.h"
#include "MGFE.h"

// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// c++ libraries
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <assert.h>

#include "EquationsMap.h"

// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
//  This function is the Basic constructor
FEMUS::FEMUS() :
    _comm ( MPI_COMM_WORLD )  // MPI_COMM_WORLD communicator
{
    // ==========================================================================
    init_femus();
    return;
}

// ============================================================================
// This function is a constructor with  communicator
FEMUS::FEMUS ( MPI_Comm comm ) :
    _comm ( comm ) // use communicator
{
    // ==========================================================================
    init_femus();// Call init class variables
    return;
}
// ============================================================================
// This function is a constructor with  fem and mesh classes
FEMUS::FEMUS ( MGUtils & mgutils )  :
    _comm ( MPI_COMM_WORLD )  // MPI_COMM_WORLD communicator
{
    // ==========================================================================
    EquationsMap FieldClass;
    mgutils.AddFieldClass ( &FieldClass );

    init_femus();           // Init class variables
    init_param ( mgutils ); // Init parameters
    init_fem();             // Init finite element
    set_mesh();              // Set mesh
    init_equation_system(); // Init equation system
    FieldClass.FillEquationMap ( get_MGExtSystem() ); // Init the equation class objects and set systems for each class
    init_systems();         // Init system data

    return;
}
// ============================================================================
// ============================================================================
// This function initializes the class
void FEMUS::init_femus (
)  // ====================================================================
{
    int flag=0;
    MPI_Initialized ( &flag );
    if ( flag ) {
        _local_MPI_Init = false;   //set _local_MPI_Init
    } else {
        _local_MPI_Init = true;    //_local_MPI_Init
    }
    // femus init
    int argc = 0;
    char ** argv = NULL;
    _start= new MGFemusInit ( argc,argv );

    _MgEquationMapInitialized = false;
    _MgMeshInitialized = false;
    return;
}

// =======================================================================
// This function
void FEMUS::init_param (
    MGUtils  &  mgutils,  ///< utils class
    int name              ///< utils class name
)   // ====================================================================
{
    _mg_utils=&mgutils;
    _mg_utils->set_name ( name );  // set ultils class
    return;
}

// ============================================================================
// This class set _mg_geomel and _mg_femap private variable
void FEMUS::init_fem (
    MGGeomEl & mggeomel, ///< geom class
    MGFEMap & mgfemap   ///< fem class
)   // ========================================================================
{
// A) setting MGGeomEl
    _mg_geomel=&mggeomel;  // ***************************************************
    if ( _mg_geomel == NULL ) {
        std::cout<< "FEMUS::init_fem: no _mg_geomel";
        abort();
    }
    /// B) setting MGFEMap (fem)
    _mg_femap=&mgfemap;  // *****************************************************
    if ( _mg_femap == NULL ) {
        std::cout<< "FEMUS::init_fem: no _mg_femap";
        abort();
    }

    return;
}

// ============================================================================
void FEMUS::init_fem ()
{
    // ========================================================================
    /// A) setting MGGeomEl
    _mg_geomel=new  MGGeomEl();
    if ( _mg_geomel == NULL ) {
        std::cout<< "FEMUS::init_fem: no _mg_geomel";
        abort();
    }

    /// B) setting MGFEMap (fem)
    _mg_femap=new MGFEMap();
    MGFE * dfe_q;
    dfe_q=new MGFE ( 2,ELTYPE );
    dfe_q->init_qua();
    MGFE * dfe_l;
    dfe_l=new MGFE ( 1,ELTYPE );
    dfe_l->init_lin();
    MGFE * dfe_k;
    dfe_k=new MGFE ( 0,ELTYPE );
    dfe_k->init_pie();

    _mg_femap->set_FE ( dfe_q ); // quadratic fem
    _mg_femap->set_FE ( dfe_l ); // linear fem
    _mg_femap->set_FE ( dfe_k ); // piecewise fem

    if ( _mg_femap == NULL ) {
        std::cout<< "FEMUS::init_fem: no _mg_femap";
        abort();
    }
    return;
}
// ============================================================================
// This function is the destructor
FEMUS::~FEMUS (
)  // ==========================================================================
{
// DO NOT TOUCH ================
    delete _mg_time_loop;
    delete _start;
//==============================

//   delete _mg_equations_map;
//   delete _mg_utils;
//   delete _mg_mesh;
    delete _mg_geomel;
    delete _mg_femap;
#ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif

}

// ============================================================================
// This function is the problem destructor
void FEMUS::terminate (
)   // =========================================================================
{
}

// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************
#ifdef   TWO_PHASE
void FEMUS::set_mgcc (
    MGSolCC & cc
)
{
    _mg_equations_map->set_mgcc ( cc );
    return;
}
#endif
//=============================================================================
// This function sets the controlled domain
void FEMUS::setCtrlDomain (
    const double xMin,
    const double xMax,
    const double yMin,
    const double yMax,
    const double zMin,
    const double zMax
)   // ========================================================================
{
    _mg_equations_map->eqnmap_ctrl_domain ( xMin,xMax,yMin,yMax,zMin,zMax );
    return;
}
//=============================================================================
/// This function returns a value from the systems
double FEMUS::GetValue ( const int  & ff,int flag )
{
    _mg_equations_map->GetValue ( ff,flag );
}
//=============================================================================
/// This function sets a value in the systems
void FEMUS::SetValue ( const int  & ff,double value )
{
    _mg_equations_map->SetValue ( ff,value );
}
//=============================================================================
/// This function sets a set of values in the systems
void FEMUS::SetValueVector ( const int  & ff,std::vector<double> value )
{
    _mg_equations_map->SetValueVector ( ff,value );
}
// =============================================================================

/// This function sets the type of problem
/*MGEquationsSystem &*/ void FEMUS::init_equation_system (
    int n_data_points,
    int n_data_cell
)   // ==========================================================================
{
    _mg_equations_map=new EquationSystemsExtendedM ( *_mg_utils,*_mg_mesh,*_mg_femap,n_data_points,n_data_cell ); // MGEquationsMap class
//     return *_mg_equations_map;
    return;
}

void FEMUS::init_systems()
{

    _mg_equations_map->init_data ( 0 );
    _mg_equations_map->setDofBcOpIc();                            // set operators
    _mg_equations_map->set_mesh_mg ( *_mg_mesh );
#ifdef HAVE_MED
    _mg_equations_map->set_mesh_med ( *_med_mesh );
#endif
    if ( _mg_geomel == NULL ) {
        std::cout<< "FEMUS::setSystem: no _mg_equations_map";
        abort();
    }
    //time loop
    _mg_time_loop=new  MGTimeLoop ( *_mg_utils,*_mg_equations_map );
    if ( _mg_time_loop == NULL ) {
        std::cout<< "FEMUS::setSystem: no _mg_time_loop";
        abort();
    }
    _MgEquationMapInitialized = true;
    //     std::cout<<"Creating interface for mesh "<<_mg_utils->_interface_mesh.c_str()<<std::endl;
//     init_interface ( _GlobInterfaceId, 2, _mg_utils->_interface_mesh.c_str() );
//     init_par_interface ( 2,true );
}




// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::set_mesh (
)   // ==========================================================================
{

    const int NoLevels= ( int ) _mg_utils->_geometry["nolevels"];
    _mg_mesh=new MeshExtended ( _start->comm(), *_mg_utils,*_mg_geomel );
    // check insanity
    if ( _mg_mesh == NULL ) {
        std::cout<< "FEMUS::setMesh: no _mg_mesh";
        abort();
    }
    if ( NoLevels != _mg_mesh->_NoLevels ) {
        std::cout << "Inconsistent Number of Levels between Mesh and SolBase"
                  << std::endl;
        abort();
    }
    // print mesh at level NoLevels-1 (linear connectivity)
    _mg_mesh->print ( NoLevels-1,0 );

#ifdef HAVE_MED
    // prind mesh at level NoLevels-1 (med format)
    std::string mesh_name= _mg_utils->get_file ( "F_MESH_READ" ); //://= _mg_utils.get_file("F_MESH_READ");
    unsigned pos = mesh_name.find ( "." );      // position of "live" in str
    std::ostringstream name;
    if ( _mg_mesh->_iproc==0 ) {
        name << _mg_utils->_mesh_dir <<  mesh_name.substr ( 0,pos ) << "_fine.med" ;
        _mg_mesh->print_med ( NoLevels-1,name.str().c_str() );
    }
#endif
    _MgMeshInitialized = true;
    return;
}


void FEMUS::setMeshTurbCase ( )
{
    set_mesh();
    setMedMesh ();
    return;
}
// *******************************************************************
// *******************************************************************
int  FEMUS::get_proc() const
{
    return _mg_mesh->_iproc;
}
// *******************************************************************



// *******************************************************************
// **************** Solve  *******************************************
// *******************************************************************
/// This function sets up the intial set
void FEMUS::solve_setup (
    int     &    t_in,                 ///< initial time iteration
    double    &   time                 ///< actual time
)
{
    const int restart      = stoi ( _mg_utils->_sim_config["restart"] ); // restart or not
    _mg_time_loop->transient_setup ( restart,t_in,time );  //  MGTimeLoop: setup

    return;
}

//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_onestep (
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double    &   time,                ///< actual time
    double    &   dt,                   ///< step time
    const int  & eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int  & eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)   // ========================================================================
{
    _mg_time_loop->transient_onestep ( t_in,t_step,print_step,time,dt,eq_min,eq_max ); ///< step time
    return;
}
// ========================================================================
void FEMUS::solve_and_update ( const int & t_in,        ///< initial time iteration
                               const int & t_step,       ///< actual time iteration
                               const int & print_step,   ///< print every
                               double & time,            ///< actual time
                               double & dt,              ///< step time
                               const int & eq_min,   ///< eq min to solve -> enum  FIELDS (equations_conf.h)
                               const int & eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
                             )
{
    _mg_time_loop->transient_solve_and_update ( t_in,t_step,print_step,time,dt,eq_min,eq_max ); ///< step time
    return;
}

// ========================================================================
// This function solves one step  for transient problems
void  FEMUS::solve_steady (
    const int & nmax_step,  ///< number max of steps
    const double & toll,  ///< tolerance
    const int  & it_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double    &   dt,         ///< inial time step
    const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int    &   eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
)   // ========================================================================
{
    _mg_time_loop->steady ( nmax_step,toll,it_step,print_step,dt,eq_min,eq_max ); ///< step time
    return;
}

// This function solves one step  for transient problems
void  FEMUS::set_uooold (
    const int & vec_from,  ///< source vector to be copied
    const int & vec_to,    ///< target vector
    const double & toll,   ///< tolerance
    const double delta_t_step_in,  //   (in)
    const int  & eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int    &   eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
)   // ========================================================================
{
    _mg_time_loop->set_uooold ( vec_from, vec_to, toll,delta_t_step_in,eq_min,eq_max ); ///< step time
    return;
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep (
    const int  & nmax_step,       ///< number max of steps         (in)
    const int  & it,              ///< initial time iteration
    const int  & t_step,          ///< actual time iteration
    const int  & print_step,      ///< print every
    double   &   time,            ///< actual time
    double   &   dt,              ///< step time
    const int  & eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int  & eq_max,          ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    bool    &    converged        ///< check if the solution converged (1->converged)     (out)
)   // ========================================================================
{
    _mg_time_loop->transient_control_onestep ( nmax_step, it,t_step,print_step,time,dt,eq_min,eq_max,converged ); ///< step time
    return;
}
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep (
    const int  & nmax_step,       ///< number max of steps         (in)
    const int  & it,              ///< initial time iteration
    const int  & t_step,          ///< actual time iteration
    const int  & print_step,      ///< print every
    double   &   time,            ///< actual time
    double   &   dt,              ///< step time
    const int  & eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int  & eq_max,          ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    std::vector<double>    controlled_eq,  ///< vector containing numbers of controlled equations
    bool    &    converged,       ///< check if the solution converged (1->converged)     (out)
    const double   &   toll             ///< tolerance
)   // ========================================================================
{
    _mg_time_loop->transient_control_onestep ( nmax_step, it,t_step,print_step,time,dt,eq_min,eq_max,controlled_eq,converged,toll ); ///< step time
    return;
}
// ========================================================================
double  FEMUS::System_functional (
    const int  & ff,                 ///< initial time iteration
    double      parameter,             ///< functional parameter
    double   &   control                   ///< step control
)
{
    return _mg_equations_map->System_functional ( ff,parameter,control );
}

//=============================================================================
// This function solves one step  for transient problems
void FEMUS::dummy_step (
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double    &   time,                ///< actual time
    double    &   dt                   ///< step time
)   // ========================================================================
{
    _mg_time_loop->dummy_step ( t_in,t_step,print_step,time,dt ); ///< step time
    return;
}
//=============================================================================
// This function solves one step for under relaxations problems (no loop inside)
void FEMUS::solve_underrelaxed_onestep  (
    const int  & it,              ///< initial time iteration
    const int  & t_step,          ///< actual time iteration
    const int  & print_step,      ///< print every
    double   &   time,            ///< actual time
    double   &   dt,              ///< step time
    const int  & eq_min,          ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int  & eq_max,          ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    std::vector<double>    controlled_eq,  ///< vector containing numbers of controlled equations
    bool    &    converged,       ///< check if the solution converged (1->converged)     (out)
    const double   &   toll             ///< tolerance
)   // ========================================================================
{
    _mg_time_loop->transient_underrelaxed_onestep ( it,t_step,print_step,time,dt,eq_min,eq_max,controlled_eq,converged,toll ); ///< step time
    return;
}
