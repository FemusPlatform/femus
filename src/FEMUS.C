#include <iostream>
#include <cstdlib>
#include <sstream>
#include <assert.h>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
// #include "Equations_conf.h"
#include "Equations_tab.h"
#include "MGTimeLoop.h"


// class include
#include "FEMUS.h"

#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionM.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"
#endif

#include "MeshExtended.h"
#include "EquationSystemsExtendedM.h"

// #include "TurbUtils.h"
// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
FEMUS::FEMUS()  :
_comm ( MPI_COMM_WORLD ) { // communicator
    // Init MPI flag --------------------------------
    int flag=0;
    MPI_Initialized ( &flag );
    if ( flag ) {
        _local_MPI_Init = false;
    } else {
        _local_MPI_Init = true;
    }
    // femus init -----------------------------------
    int argc = 0;
    char **argv = NULL;
    _start=new  MGFemusInit ( argc,argv,_comm );

    return;
}

// ============================================================================
// This function is a constructor with  communicator
FEMUS::FEMUS (
    MPI_Comm comm
) :
_comm ( comm ) // communicator
//   _num_mgmesh(0)
{
    // n of femus_meshes
    // transient system
    // Init MPI flag
    int flag=0;
    MPI_Initialized ( &flag );
    if ( flag ) {
        _local_MPI_Init = false;
    } else {
        _local_MPI_Init = true;
    }

    // femus init
    int argc = 0;
    char **argv = NULL;
    _start= new MGFemusInit ( argc,argv );

    return;
}


// void FEMUS::init_param(MGUtils &mgutils, TurbUtils &Parameters,int name){
//   _mg_utils=&mgutils;
//   _mg_utils->set_name(name);
//   _mg_utils->set_turbulence_info(Parameters);
//   return;
// }


// =======================================================================
void FEMUS::init_param (
    MGUtils    &mgutils,
    int name
) { // ====================================================================
    _mg_utils=&mgutils;
    _mg_utils->set_name ( name );
    return;
}
// ============================================================================
void FEMUS::init_fem (
    MGGeomEl &mggeomel,
    MGFEMap &mgfemap
) { // ========================================================================
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
// This function is the destructor
FEMUS::~FEMUS() {
    // ==========================================================================
//   delete _start;
    delete _mg_time_loop;delete _start;
//   delete _mg_equations_map;
//   delete _mg_utils;
//   delete _mg_mesh;
//   delete _mg_femap;
#ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif

}

// ============================================================================
// This function is the problem destructor
void FEMUS::terminate (
) { // =========================================================================

}

// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************
#ifdef   TWO_PHASE
void FEMUS::set_mgcc ( MGSolCC &cc

                     ) {
    _mg_equations_map->set_mgcc ( cc );
    return;
}
#endif
//=============================================================================
// This function sets the controlled domain 
void FEMUS::setCtrlDomain(
     const double xMin,
     const double xMax,
     const double yMin,
     const double yMax,
     const double zMin,
     const double zMax
) { // ========================================================================
  _mg_equations_map->eqnmap_ctrl_domain(xMin,xMax,yMin,yMax,zMin,zMax);
  return;
}
// // // =============================================================================
// // // This function sets the type of problem
void FEMUS::setSystem (
    const std::vector<FIELDS> &pbName,
    int n_data_points,
    int n_data_cell
) { // ==========================================================================
    _mg_equations_map=new EquationSystemsExtendedM ( *_mg_utils,*_mg_mesh,*_mg_femap,n_data_points,n_data_cell ); // MGEquationsMap class
//     _mg_equations_map->read_par();
#ifdef PRINT_INFO  // ---- info ---------------
//     _mg_equations_map->print_par();       // print parameters
#endif
    _mg_equations_map->init_data ( 0 );
    _mg_equations_map->init ( pbName );                           // adds the equations to the map
    _mg_equations_map->setDofBcOpIc();                            // set operators
    _mg_equations_map->set_mesh_mg ( *_mg_mesh );
#ifdef HAVE_MED
    _mg_equations_map->set_mesh_med ( *_med_mesh );
#endif
//   }
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

    init_interface ( _GlobInterfaceId, 2, _mg_utils->_interface_mesh.c_str() );
    init_par_interface ( 2,true );
    return;
}


// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::setMesh (
) { // ==========================================================================

//   const int NoLevels= _mg_utils->get_par("nolevels");  // numb of Level
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

void FEMUS::setMeshTurbCase ( ) { // ==========================================================================
    setMesh();
    setMedMesh ();
    return;
}
// // // *******************************************************************
// // // *******************************************************************
int  FEMUS::get_proc() const {
    return _mg_mesh->_iproc;
}
// // // *******************************************************************

void FEMUS::InitTurbulence() {
    bool DynTurb, TherTurb;

    ( stoi ( _mg_utils->_sim_config["DynamicalTurbulence"] ) > 0 ) ? DynTurb = true:false;
    ( stoi ( _mg_utils->_sim_config["ThermalTurbulence"] )   > 0 ) ? TherTurb = true:false;

//      MyAssert ( _NodeWallDist.size() != 0, "FEMUS::InitTurbulence() _NodeWallDist not computed! \n" );

    TurbUtils *Parameters = new TurbUtils ( get_proc(), _mg_mesh->_NoLevels, _NodeMap, DynTurb, TherTurb );
    _mg_utils->set_turbulence_info ( Parameters );
    return;
}

void FEMUS::CalcTurbulence() {
    int FinerLevel = _mg_mesh->_NoLevels;
    if ( _mg_utils->_TurbParameters->_IsWallDistSet == false ) {
        ParaMEDMEM::MEDCouplingFieldDouble *Dist = getProcSolution ( "C",1,0,FinerLevel-1 );
        MEDLoader::WriteField ( "RESU_MED/wd"+to_string ( get_proc() ) +".med",Dist,1 );
        _mg_utils->_TurbParameters->SetWallDistAtLevel ( FinerLevel-1, Dist );
        _mg_utils->_TurbParameters->_IsWallDistSet = true;
        Dist->decrRef();
    }
    ParaMEDMEM::MEDCouplingFieldDouble *K1W = getProcSolution ( "K1W",1,0,FinerLevel-1 ); // CONTROLLARE CHE LA SOLUZIONE VENGA PRESA DA OGNI SINGOLO LIVELLO
    ParaMEDMEM::MEDCouplingFieldDouble *K2K = getProcSolution ( "K2K",1,0,FinerLevel-1 ); //
    _mg_utils->_TurbParameters->CalcMuTurb ( K2K,K1W,FinerLevel -1 );
    if ( stoi ( _mg_utils->_sim_config["ThermalTurbulence"] )   > 0 ) {
        ParaMEDMEM::MEDCouplingFieldDouble *TK = getProcSolution ( "TK",1,0,FinerLevel-1 ); // CONTROLLARE CHE LA SOLUZIONE VENGA PRESA DA OGNI SINGOLO LIVELLO
        ParaMEDMEM::MEDCouplingFieldDouble *TW = getProcSolution ( "TK2",1,0,FinerLevel-1 ); //
        _mg_utils->_TurbParameters->CalcAlphaTurb ( K2K,K1W,TK,TW,FinerLevel -1 );
        TW->decrRef();
        TK->decrRef();
    }
    K1W->decrRef();
    K2K->decrRef();

    return;
}

// // // *******************************************************************
// // // **************** Solve  *******************************************
// // // *******************************************************************
/// This function sets up the intial set
void FEMUS::solve_setup (
    int         &t_in,                 ///< initial time iteration
    double       &time                 ///< actual time
) {
    const int restart      = stoi ( _mg_utils->_sim_config["restart"] ); // restart or not
    _mg_time_loop->transient_setup ( restart,t_in,time );  //  MGTimeLoop: setup

//   if(_mg_utils->_TurbParameters != NULL){
//     CalcTurbulence();
//   }

    return;
}

//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_onestep (
    const int   &t_in,                 ///< initial time iteration
    const int   &t_step,               ///< actual time iteration
    const int   &print_step,            ///< print every
    double       &time,                ///< actual time
    double       &dt,                   ///< step time
    const int   &eq_min, ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int   &eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) { // ========================================================================
    _mg_time_loop->transient_onestep ( t_in,t_step,print_step,time,dt,eq_min,eq_max ); ///< step time
    return;
}

void FEMUS::solve_and_update ( const int &t_in,         ///< initial time iteration
                               const int &t_step,        ///< actual time iteration
                               const int &print_step,    ///< print every
                               double &time,             ///< actual time
                               double &dt,               ///< step time
                               const int &eq_min,    ///< eq min to solve -> enum  FIELDS (equations_conf.h)
                               const int &eq_max    ///< eq max to solve -> enum  FIELDS (equations_conf.h)
                             ) {
    _mg_time_loop->transient_solve_and_update ( t_in,t_step,print_step,time,dt,eq_min,eq_max ); ///< step time
    return;
}


// This function solves one step  for transient problems
void  FEMUS::solve_steady (
    const int &nmax_step,   ///< number max of steps
    const double &toll,   ///< tolerance
    const int   &it_step,               ///< actual time iteration
    const int   &print_step,            ///< print every
    double       &dt,         ///< inial time step
    const int   &eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h)
    const int       &eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h)
) { // ========================================================================
    _mg_time_loop->steady ( nmax_step,toll,it_step,print_step,dt,eq_min,eq_max ); ///< step time
    return;
}

// This function solves one step  for transient problems
void  FEMUS::set_uooold (
    const int &flag,   ///<  0 xold-> x_ooold   1 x_ooold-> xold
    const double &toll,   ///< tolerance
    const double delta_t_step_in,  //   (in)
    const int   &eq_min,     ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
    const int       &eq_max ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
) { // ========================================================================
    _mg_time_loop->set_uooold ( flag, toll,delta_t_step_in,eq_min,eq_max ); ///< step time
    return;
}
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_control_onestep (
    const int   &t_in,                 ///< initial time iteration
    const int   &t_step,               ///< actual time iteration
    const int   &print_step,            ///< print every
    double       &time,                ///< actual time
    double       &dt                   ///< step time
) { // ========================================================================
    _mg_time_loop->transient_control_onestep ( t_in,t_step,print_step,time,dt ); ///< step time
    return;
}

double  FEMUS::System_functional (
    const int   &ff,                 ///< initial time iteration
    double      parameter,             ///< functional parameter
    double      &control                   ///< step control
) {
    return _mg_equations_map->System_functional ( ff,parameter,control );
}

#ifdef HAVE_MED

// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
void FEMUS::init_interface (
    const int interface_name,
    std::vector<int> IDSvec,
    int order_cmp,
    const std::string &medfile_name,  // medfile name    (in)
    bool on_nodes,
    const int index_medmesh                       // med-mesh index  (in)
) { // =========================================================================

    std::vector<std::string> vG ( IDSvec.size() );
    std::cout<<" Creating an interface with id "<<interface_name<< " containing groups ";
    for ( int j=0; j<IDSvec.size(); j++ ) {
        vG[j] = to_string ( IDSvec[j] );
        std::cout<<IDSvec[j]<<" ";
    }
    std::cout<<std::endl;

    std::string mesh_dir, localFile, filename;
    std::vector<std::string> meshNames, MeshNames, FieldNames;

    GetInfo ( medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );

    ParaMEDMEM::MEDCouplingFieldDouble *FieldContainingMap;
    ParaMEDMEM::DataArrayDouble *MedToMgMapArray;

    int FinerLevel = _mg_utils->_geometry["nolevels"] -1;

    if ( on_nodes ) {
        FieldContainingMap = MEDLoader::ReadField ( ParaMEDMEM::ON_NODES, filename.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), 0, "FinerLevelNodeIDS_Lev_"+to_string ( FinerLevel ), -1,-1 );
    } else {
        FieldContainingMap = MEDLoader::ReadField ( ParaMEDMEM::ON_CELLS, filename.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), 0, "MG_cell_id_Lev_"+to_string ( FinerLevel ), -1,-1 );
    }

    MedToMgMapArray = FieldContainingMap->getArray();

    ParaMEDMEM::MEDCouplingUMesh *support;
    int id_level=-1;
    if ( IDSvec[0]<10 ) {
        id_level=0;
    }

    support = MEDLoader::ReadUMeshFromGroups ( localFile.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), id_level,vG );
    std::cout<<" NUMBER OF NODES "<<support->getNumberOfNodes() <<std::endl;

    MEDLoader::WriteUMesh ( "RESU_MED/Interface_mesh.med",support,true );

    init_interface ( interface_name, order_cmp, support, MedToMgMapArray );
    FieldContainingMap->decrRef();

    return;
}


void FEMUS::init_interface (
    const int interface_name,
    int order_cmp,
    const std::string &medfile_name,  // medfile name    (in)
    bool on_nodes
) { // =========================================================================

    // Reading mesh Names from med file  ------------------------------

    std::string mesh_dir, localFile, filename;
    std::vector<std::string> meshNames, MeshNames, FieldNames;

    GetInfo ( medfile_name, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );
    ParaMEDMEM::MEDCouplingFieldDouble *FieldContainingMap;
    ParaMEDMEM::DataArrayDouble *MedToMgMapArray;

    int FinerLevel = _mg_utils->_geometry["nolevels"] -1;
    if ( on_nodes ) {
        FieldContainingMap = MEDLoader::ReadField ( ParaMEDMEM::ON_NODES, filename.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), 0, "FinerLevelNodeIDS_Lev_"+to_string ( FinerLevel ), -1,-1 );
    } else {
        FieldContainingMap = MEDLoader::ReadField ( ParaMEDMEM::ON_CELLS, filename.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), 0, "MG_cell_id_Lev_"+to_string ( FinerLevel ), -1,-1 );
    }

    MedToMgMapArray = FieldContainingMap->getArray();
    int id_level=0;
    ParaMEDMEM::MEDCouplingUMesh *support = MEDLoader::ReadUMeshFromFile ( localFile.c_str(), "Mesh_Lev_"+to_string ( FinerLevel ), id_level );

    init_interface ( interface_name, order_cmp, support, MedToMgMapArray );

    FieldContainingMap->decrRef();

    return;
}

void FEMUS::init_par_interface (
    int order_cmp,                    /**< Order (piecewise, linear, quadratic) */
    bool on_nodes
) {
    // Reading mesh Names from med file  ------------------------------
    const int levels = _mg_mesh->_NoLevels;
    ParaMEDMEM::MEDCouplingFieldDouble *NodeMap = NULL;
    ParaMEDMEM::MEDCouplingUMesh *support = NULL;
    getNodeMapAndProcMeshAtLevel ( levels-1, support, NodeMap );
    ParaMEDMEM::DataArrayDouble *MedToMgMapArray = NodeMap->getArray();

    _ParallelInterfaceId = 1235 + get_proc();
    init_interface ( _ParallelInterfaceId, order_cmp, support, MedToMgMapArray );

    return;
}


void FEMUS::getNodeMapAndProcMeshAtLevel ( int level,  ParaMEDMEM::MEDCouplingUMesh *&FemusPar, ParaMEDMEM::MEDCouplingFieldDouble *&NodeMap ) {
    const int proc = get_proc();
    std::string mesh_dir, localFile, filename;
    std::vector<std::string> meshNames, MeshNames, FieldNames;
    GetInfo ( _mg_utils->_interface_mesh, mesh_dir, localFile, filename, meshNames, MeshNames, FieldNames );

    ParaMEDMEM::MEDCouplingUMesh *FemusPart       = MEDLoader::ReadUMeshFromFile ( filename.c_str(),"Mesh_Lev_"+to_string ( level ), 0 ); // ORIGINAL FEMUS FROM MED FILE

    ParaMEDMEM::MEDCouplingFieldDouble *ProcField = MEDLoader::ReadField (
                ParaMEDMEM::ON_CELLS,
                filename.c_str(),
                "Mesh_Lev_"+to_string ( level ),
                0,
                "Proc_Lev_"+to_string ( level ),
                -1,-1
            );             // PIECEWISE FIELD CONTAINING PROC IDS
    ParaMEDMEM::MEDCouplingFieldDouble *FinerLevelIDS = MEDLoader::ReadField ( ParaMEDMEM::ON_NODES,
            filename.c_str(),
            "Mesh_Lev_"+to_string ( level ), 0, "FinerLevelNodeIDS_Lev_"+to_string ( level ), -1,-1 );
    double *FinLevId = const_cast<double *> ( FinerLevelIDS->getArray()->getPointer() );
    std::vector<int> CellIdProc;
    for ( int i=0; i<ProcField->getNumberOfTuples(); i++ ) {
        if ( ( int ) ProcField->getIJ ( i,0 ) == proc ) {
            CellIdProc.push_back ( i );
        }
    }

    int *CellId = new int[CellIdProc.size()];                                                                      // ARRAY CONTAINING CELL IDS OF SINGLE PROC MESH
    for ( int i=0; i<CellIdProc.size(); i++ ) {
        CellId[i] = CellIdProc[i];
    }

    FemusPar = FemusPart->buildPartOfMySelf ( CellId, CellId + CellIdProc.size() );  // EXTRACTION OF PROC MESH FROM GLOBAL MESH
    FemusPar->zipCoords();                                                                                         // It changes nodes numbering -> new numbering is from 0 to FemusPar->getNumberOfNodes
    FemusPar->setName ( "FemusProcMesh_Lev_"+to_string ( level ) );

    // CREATING A MAP FROM LOCAL NODE NUMBERING TO GLOBAL NODE NUMBERING - RELATIVE TO MED FILE
    double *NodesId = NULL;
    ParaMEDMEM::DataArrayDouble *NodeArray = ParaMEDMEM::DataArrayDouble::New();
    NodeArray->alloc ( FemusPar->getNumberOfNodes(),1 );
    NodesId = const_cast<double *> ( NodeArray->getPointer() );

    std::vector<int> FeParCellConn, FeTotCellConn;       // VECTOR FOR PARALLEL AND GLOBAL MESH CONNECTIVITIES

    int FemusCells        = FemusPar->getNumberOfCells();
    int FemusNodesPerCell = FemusPar->getNumberOfNodesInCell ( 0 );

    for ( int i_cell=0; i_cell<FemusCells; i_cell++ ) {
        FemusPar->getNodeIdsOfCell ( i_cell,FeParCellConn );
        FemusPart->getNodeIdsOfCell ( CellId[i_cell],FeTotCellConn );
        for ( int i_cnode=0; i_cnode<FemusNodesPerCell; i_cnode++ ) {
            NodesId[FeParCellConn[i_cnode]] = FinLevId[FeTotCellConn[i_cnode]];
        }
        FeParCellConn.clear();
        FeTotCellConn.clear();
    }// END MAP ---------------------------------------------------------------------------------

    NodeMap = ParaMEDMEM::MEDCouplingFieldDouble::New ( ParaMEDMEM::ON_NODES ) ;
    NodeMap->setMesh ( FemusPar );
    NodeMap->setArray ( NodeArray );
    NodeMap->setName ( "NodesID_Lev_"+to_string ( level ) );

    return ;
}



void FEMUS::init_interface (
    const int interface_name,
    int order_cmp,
    ParaMEDMEM::MEDCouplingUMesh *InterfaceMesh,
    ParaMEDMEM::DataArrayDouble *MEDToMgMapArray,
    int interface_id
) {

    double *MapArray = const_cast<double *> ( MEDToMgMapArray->getPointer() );
    std::map<int, int> MedToMg;

    int celle = InterfaceMesh->getNumberOfCells();
    int NodesPerCell = InterfaceMesh->getNumberOfNodesInCell ( 0 );

    std::vector<int> MgNodesIds, MedNodesIds;

    // Storing cell connectivities inside vector MgNodesIds -> numeration of mesh mg
    for ( int i =0; i<celle; i++ ) {
        InterfaceMesh->getNodeIdsOfCell ( i,MgNodesIds );
    }

    // Renumbering nodes and cells -> new numeration in accordance with mesh total number of nodes and cells
    InterfaceMesh->zipCoords();

    // Storing cell connectivities inside vector MedNodesIds -> numeration of mesh med
    for ( int i =0; i<celle; i++ ) {
        InterfaceMesh->getNodeIdsOfCell ( i,MedNodesIds );
    }

    // Creating the map
    for ( int i=0; i<celle*NodesPerCell; i++ ) {
        MedToMg[MedNodesIds[i]] = ( int ) MapArray[MgNodesIds[i]];
    }

    const int NODI = InterfaceMesh->getNumberOfNodes();
    int *map_med = new int[NODI];
    int *map_mg  = new int[NODI];

    for ( int i = 0; i<NODI; i++ ) {
        map_med[i] = i;
        map_mg[i]  = MedToMg[i];
    }

    MgNodesIds.clear();
    MedNodesIds.clear();

    if ( interface_id==0 ) printf ( "\n \033[1;31m Interface %d \
                                          support set to mesh without volume group and order %d \
                                          \033[0m\n",interface_name,order_cmp );
    else std::cout << "\n \033[1;31m Interface "<<interface_name
        <<" support set to mesh with interface " <<interface_id
        <<" and order "<<order_cmp
        << "\033[0m\n";

    InterfaceFunctionM *fun = new InterfaceFunctionM;
    fun->set_maps ( map_med,map_mg, NODI );
    fun->set_order ( order_cmp );
    fun->set_mg_mesh ( _mg_mesh );
    fun->set_support_med ( InterfaceMesh );
    fun->set_NumberOfNodes ( NODI );

    _mg_equations_map->add_interface_fun ( interface_name, fun );
    MapArray = NULL;
    MedToMg.clear();
    return;
}


void FEMUS::write_Boundary_value (
    int id_boundary_name ,      ///< identity interface name (in)
    std::string mgsystem_name, ///< system name          (in)
    int n_cmp,             ///< from variable system (in)
    int first_cmp               ///< to variable system   (in)

) {
    _mg_equations_map->write_Boundary_value ( id_boundary_name, mgsystem_name,n_cmp,first_cmp );
    return;
}



// =============================================================================
//This function reads and sets the med-mesh
// (also the libmesh calling the other setMesh function)
void FEMUS::setMedMesh ( ) {
    std::string dataFile = _mg_utils->_mesh_dir + _mg_utils->_interface_mesh;
    int l = dataFile.size();
    int proc = get_proc();

    // CHECK THAT MGMESH AND EQUATION MAP ARE INITIALIZED
    MyAssert ( _MgMeshInitialized, "FEMUS::setMedMesh MGMESH not initialized" );
    MyAssert ( dataFile.substr ( l-4 ) == ".med", "FEMUS::setMesh: " + dataFile +" does not exist!" );

    std::string localFile = dataFile;
    std::vector<std::string> meshNames = MEDLoader::GetMeshNames ( localFile.c_str() );

    MyAssert ( meshNames.size() > 0.5, " FEMUS::setMesh : no meshes in the file'" );

    _med_mesh = MEDLoader::ReadUMeshFromFile ( localFile.c_str(), meshNames[0].c_str(), 0 );
    if ( _med_mesh == NULL ) {
        std::cout<<  " FEMUS::setMesh : unable to read the med-mesh'";
        abort();
    }

    // BUILDING LOCAL PROC MESH FROM GLOBAL MED MESH ================================================================
    int FinerLevel = _mg_utils->_geometry["nolevels"]-1;

    std::cout<<" Creating proc mesh for level "<<FinerLevel<<std::endl;

    for ( int lev=0; lev<FinerLevel+1; lev++ ) {
        ParaMEDMEM::MEDCouplingUMesh *FemusPart       = MEDLoader::ReadUMeshFromFile ( dataFile,"Mesh_Lev_"+to_string ( lev ), 0 ); // ORIGINAL FEMUS FROM MED FILE
        ParaMEDMEM::MEDCouplingFieldDouble *ProcField = MEDLoader::ReadField ( ParaMEDMEM::ON_CELLS,
                dataFile,
                "Mesh_Lev_"+to_string ( lev ), 0, "Proc_Lev_"+to_string ( lev ), -1,-1 );             // PIECEWISE FIELD CONTAINING PROC IDS
        ParaMEDMEM::MEDCouplingFieldDouble *CellField = MEDLoader::ReadField ( ParaMEDMEM::ON_CELLS,
                dataFile,
                "Mesh_Lev_"+to_string ( lev ), 0, "MG_cell_id_Lev_"+to_string ( lev ), -1,-1 );
        double *MgCellId = const_cast<double *> ( CellField->getArray()->getPointer() );

        ParaMEDMEM::MEDCouplingFieldDouble *FinerLevelIDS = MEDLoader::ReadField ( ParaMEDMEM::ON_NODES,
                dataFile,
                "Mesh_Lev_"+to_string ( lev ), 0, "FinerLevelNodeIDS_Lev_"+to_string ( lev ), -1,-1 );
        double *FinLevId = const_cast<double *> ( FinerLevelIDS->getArray()->getPointer() );

        std::vector<int> CellIdProc;
        for ( int i=0; i<ProcField->getNumberOfTuples(); i++ ) {
            if ( ( int ) ProcField->getIJ ( i,0 ) == proc ) {
                CellIdProc.push_back ( i );
            }
        }

        int *CellId = new int[CellIdProc.size()];                                                                      // ARRAY CONTAINING CELL IDS OF SINGLE PROC MESH
        for ( int i=0; i<CellIdProc.size(); i++ ) {
            CellId[i] = CellIdProc[i];
        }

        ParaMEDMEM::MEDCouplingUMesh *FemusPar = FemusPart->buildPartOfMySelf ( CellId, CellId + CellIdProc.size() );  // EXTRACTION OF PROC MESH FROM GLOBAL MESH
        FemusPar->zipCoords();                                                                                         // It changes nodes numbering -> new numbering is from 0 to FemusPar->getNumberOfNodes
        FemusPar->setName ( "FemusProcMesh_Lev_"+to_string ( lev ) );

        // CREATING A MAP FROM LOCAL NODE NUMBERING TO GLOBAL NODE NUMBERING - RELATIVE TO MED FILE
        double *NodesId = NULL;
        ParaMEDMEM::DataArrayDouble *NodeArray = ParaMEDMEM::DataArrayDouble::New();
        NodeArray->alloc ( FemusPar->getNumberOfNodes(),1 );
        NodesId = const_cast<double *> ( NodeArray->getPointer() );

        std::vector<int> FeParCellConn, FeTotCellConn;       // VECTOR FOR PARALLEL AND GLOBAL MESH CONNECTIVITIES

        int FemusCells        = FemusPar->getNumberOfCells();
        int FemusNodesPerCell = FemusPar->getNumberOfNodesInCell ( 0 );

        for ( int i_cell=0; i_cell<FemusCells; i_cell++ ) {
            FemusPar->getNodeIdsOfCell ( i_cell,FeParCellConn );
            FemusPart->getNodeIdsOfCell ( CellId[i_cell],FeTotCellConn );
            for ( int i_cnode=0; i_cnode<FemusNodesPerCell; i_cnode++ ) {
                NodesId[FeParCellConn[i_cnode]] = FinLevId[FeTotCellConn[i_cnode]];
            }
            FeParCellConn.clear();
            FeTotCellConn.clear();
        }// END MAP ---------------------------------------------------------------------------------

        _NodeMap.push_back ( ParaMEDMEM::MEDCouplingFieldDouble::New ( ParaMEDMEM::ON_NODES ) );
        _NodeMap[lev]->setMesh ( FemusPar );
        _NodeMap[lev]->setArray ( NodeArray );
        _NodeMap[lev]->setName ( "NodesID_Lev_"+to_string ( lev ) );

        MEDLoader::WriteUMesh ( "RESU_MED/FEMUS_PAR_"+to_string ( proc ) +".med", FemusPar, false );
        MEDLoader::WriteFieldUsingAlreadyWrittenMesh ( "RESU_MED/FEMUS_PAR_"+to_string ( proc ) +".med",_NodeMap[lev] );

        FemusPar->decrRef();
        FemusPart->decrRef();
        NodeArray->decrRef();
        CellField->decrRef();
        FinerLevelIDS->decrRef();
        ProcField->decrRef();
        CellId  = NULL;
        MgCellId = NULL;
        CellIdProc.clear();
    }

    return;
}

void FEMUS::setExtField ( const std::string &systemName,
                          ParaMEDMEM::MEDCouplingFieldDouble *bcField ) {
    _mg_equations_map->write_Boundary_value ( systemName,bcField );
    return;
}

 ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::GetExtField ( const std::string &systemName) {
     ParaMEDMEM::MEDCouplingFieldDouble * ExtField = _mg_equations_map->GetField ( systemName);
    return ExtField;
}

// ===================================================================
void FEMUS::setAnalyticSource (
//   const std::string & bcName,           // boundary name
    int interface_name,
    int n_cmp,
    const std::string &bcExpression       // boundary symbolic expr
) {
    _mg_equations_map->setBC ( interface_name,n_cmp , bcExpression.c_str() );
    return;
}



void FEMUS::setFieldSource (
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble *srcField ) {

    if ( srcField ==NULL ) {
        return;
    }
    _mg_equations_map->setBC ( interface_name,n_cmp, srcField );

    return;
}



// ============================================================================
/// This function gets all the values on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble *FEMUS::getValuesOnBoundary (  // field (out)
    int  interface_name,                     // boundary name (char*) (in)
    const std::string &systemName,                   // system name           (in)
    int n_cmp,                                        // component             (in)
    int first_cmp                                        // component             (in)
) {
    return _mg_equations_map->getValuesOnBoundary_nodes ( interface_name,systemName.c_str(),n_cmp,first_cmp );
}


ParaMEDMEM::MEDCouplingFieldDouble *FEMUS::getProcSolution (  // field (out)
    const std::string &systemName,                   // system name           (in)
    int n_cmp,                                        // component             (in)
    int first_cmp,                                        // component             (in)
    int Level
) {

    MyAssert ( _NodeMap.size() >Level, "FEMUS::getProcSolution _NodeMap not initialized for level " +to_string ( Level ) + " maximum levels "+to_string ( _NodeMap.size() ) );
    return _mg_equations_map->getProcValues ( _NodeMap[Level],_GlobInterfaceId,systemName.c_str(),n_cmp,first_cmp, Level );
}

// ============================================================================
/// This function gets the actual mesh of FEMUS problem
const ParaMEDMEM::MEDCouplingUMesh *FEMUS::getUMesh (
    int name
) {
    return _mg_equations_map->getUMeshCoupling ( name );
}
// // // ===========================================================================
/// This function gets the original mesh of FEMUS problem
const ParaMEDMEM::MEDCouplingUMesh *FEMUS::getUMesh_orig (
    int name
) {
    return _mg_equations_map->getUMeshCoupling_orig ( name );
}
// ------------------------------------------------
/// This function moves the FEMUS interface according to a given displacement field
void FEMUS::update_interface (
    const int interface_name,
    int n_cmp,
    const   std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> &srcField
) { // =========================================================================
    ParaMEDMEM::DataArrayDouble   *mg_disp= ParaMEDMEM::DataArrayDouble::Meld ( srcField[0]->getArray(),
                                            srcField[1]->getArray() );
    if ( n_cmp==3 ) {
        mg_disp= ParaMEDMEM::DataArrayDouble::Meld ( mg_disp,srcField[2]->getArray() );
    }
    std::cout<< "FEMUS::UPDATE: Tuples  src_field  "<< mg_disp ->getNumberOfTuples() <<
    "  Comp    " << mg_disp ->getNumberOfComponents()
    << " and x value is " <<  mg_disp->getIJ ( 0,0 ) << " and y value is " <<  mg_disp->getIJ ( 0,1 )
    <<std::endl;
    std::cout << "FEMUS::UPDATE: support  set to boundary with name "<<
    interface_name  << "\n";


    ParaMEDMEM::MEDCouplingUMesh *support;
    ParaMEDMEM::MEDCouplingUMesh *support_up;
//   int id_level=-1;  if(interface_id_1<10) {id_level=0;}
//   support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
    support= ( getUMesh_orig ( interface_name ) )->clone ( 1 );
    support_up= ( getUMesh ( interface_name ) )->clone ( 1 );
// support->zipCoords();

    ParaMEDMEM::DataArrayDouble *coord;
    ParaMEDMEM::DataArrayDouble *coord_up;

    coord= support->getCoords();
    coord_up= support_up->getCoords();
    int npt=coord->getNumberOfTuples();
    int ncomp= coord->getNumberOfComponents();

    ParaMEDMEM::DataArrayDouble *new_cord=ParaMEDMEM::DataArrayDouble::Add ( coord,mg_disp );

    std::cout<< "Tuples    "<< npt << "Comp    " << ncomp<<std::endl;
    support->setCoords ( new_cord );

    InterfaceFunctionM *fct = _mg_equations_map->get_interface_fun ( interface_name );
    fct->update_support ( support );


    return;
}

void FEMUS::GetInfo (
    string medfile_name,
    std::string &mesh_dir,
    std::string &localFile,
    std::string &filename,
    std::vector<std::string> &meshNames,
    std::vector<std::string> &MeshNames,
    std::vector<std::string> &FieldNames,
    const int index_medmesh
) {
    mesh_dir=_mg_utils->_mesh_dir;
    localFile=mesh_dir+medfile_name;

    meshNames = MEDLoader::GetMeshNames ( localFile.c_str() );
    // Reading mesh from the index_medmesh med-mesh name
    if ( meshNames.size() < 1 ) {
        std::cout<<  " FEMUS::setMesh : no meshes in the file'";
    }
//   std::string MMMM = medfile_name;
//   int posP = MMMM.find("_");  // position of "live" in str
//
//   filename =   mesh_dir+MMMM.substr(0,posP)  + "_MedToMg.med" ;
//   std::cout<<"\n 1:"<<localFile<<"\n 2:"<<filename<<std::endl;
    filename =localFile;

    FieldNames = MEDLoader::GetAllFieldNames ( filename.c_str() );
    MeshNames = MEDLoader::GetMeshNames ( filename.c_str() );

    return;
};

#endif
//=============================================================================
// This function solves one step  for transient problems
void FEMUS::dummy_step (
    const int   &t_in,                 ///< initial time iteration
    const int   &t_step,               ///< actual time iteration
    const int   &print_step,            ///< print every
    double       &time,                ///< actual time
    double       &dt                   ///< step time
) { // ========================================================================
    _mg_time_loop->dummy_step ( t_in,t_step,print_step,time,dt ); ///< step time
    return;
}


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 

