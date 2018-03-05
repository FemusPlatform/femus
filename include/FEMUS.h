#ifndef __FEMUS__
#define __FEMUS__

#include <vector>
#include "Solverlib_conf.h"

#ifdef HAVE_MED
namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}
#endif

class MeshExtended;
class MGFemusInit;
class EquationSystemsExtendedM;
class MGSystem;
class MGUtils;
class MGGeomEl;
class MGFEMap;
class MGTimeLoop;
class BoundInterp;

#ifdef   TWO_PHASE_LIB 
class MGSolCC;
#endif

class TurbUtils;
class IbUtils;

class FEMUS
{

  // PROTECTED VARIABLES ================================================================
protected:

  //data communication (defined in Constructor) 
  MGFemusInit * _start;		// start function
  MPI_Comm _comm;		// communicator
  bool _local_MPI_Init;		// initial mpi flag

  // param and file data (defined in init(.,.)) 
  MGUtils *_mg_utils;		// param and files

  // fem ------- (defined in init(.,.))  
  MGGeomEl *_mg_geomel;		// set element type
  MGFEMap *_mg_femap;		// fem map

  MGTimeLoop *_mg_time_loop;	// transient
  // data meshes
  MeshExtended *_mg_mesh;	// FEMus-mesh
  
#ifdef HAVE_MED
  ParaMEDMEM::MEDCouplingUMesh *         _med_mesh;     // Med-mesh
#endif  
  
  bool _MgMeshInitialized = false;

  EquationSystemsExtendedM *_mg_equations_map;	// system
  bool _MgEquationMapInitialized = false;

// #ifdef HAVE_MED
//     ParaMEDMEM::MEDCouplingUMesh * _med_mesh;	// Med-mesh
// #endif

  // PUBLIC VARIABLES ===================================================================
public:

  // Constructor-Destructor
    FEMUS ();
    FEMUS (MPI_Comm comm);

  const int _GlobInterfaceId = 1234; // default interface id for interface covering entire domain
  int _ParallelInterfaceId; 
    
  void init_param (MGUtils & mgutils, TurbUtils & Parameters, int name = 0);

  void init_param (MGUtils & mgutils, int name = 0);
  void init_fem (MGGeomEl & mggeomel, MGFEMap & mgfemap);
   ~FEMUS ();
  void terminate ();

#ifdef HAVE_MED
  //! MED field containing node interpolated wall distance
    std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> _NodeWallDist;
  //! MED field containing node numbering map between med local parallel mesh and med global mesh
    std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> _NodeMap;
#endif
 #ifdef  TWO_PHASE_LIB 
void set_mgcc(MGSolCC & cc);
#endif
  /// MESH RELATED FUNCTIONS -------------------------------------------------
  //! Default setMesh function
  void setMesh ();
  //! setMesh function and computation of parallel mesh
  void setMeshTurbCase ();
  //! String for assigning a mesh name label
    std::string MeshName;
  //! Function assigning mesh name label
  inline void setMeshName (std::string meshname)
  {
    MeshName = meshname;
  };
  

  inline int IsVelFieldActive(){
    int a=0;
    if(stoi(_mg_utils->_sim_config["NavierStokes"])!=0 ||  stoi(_mg_utils->_sim_config["FluidStructure"])!=0) a = 1;
    return a;
  }
  
  
  //! Function getting mesh name label
  inline std::string getMeshName ()
  {
    return MeshName;
  };
  //! Function getting mesh dir 
  inline string GetMeshDir ()
  {
    return _mg_utils->_mesh_dir;
  };
  //! Function returning _mg_mesh pointer
  const MeshExtended & get_MGMesh ()
  {
    return *_mg_mesh;
  };
  
  
  //! Initialization of equations to solve
  void setSystem (const std::vector < FIELDS > &pbName,
		  int n_data_points = 0, int n_data_cell = 0);

  void MyAssert (bool cond, std::string message)
  {
    if (cond == false)
      {
	std::cout << "\033[1;31m " << message << " \033[0m" << std::endl;
	abort ();
      }
  }

  //! Function returning the proc id
  int get_proc () const;



//=============================================================================
//   FUNCTIONS RELATED TO PROBLEM SOLUTION
//=============================================================================

  //! This function sets up the intial set
  void solve_setup (int &t_in,	///< initial time iteration
		    double &time	///< actual time
                   );
  
// This function solves one step  for transient problems
  void solve_steady (const int &nmax_step,	///< number max of steps
		     const double &toll,	///< tolerance
		     const int &it_step,	///< actual num iteration
		     const int &print_step,	///< print every
		     double &dt,	        ///< inial time step 
		     const int &eq_min = 0,	///< eq min to solve -> enum  FIELDS (equations_conf.h) 
		     const int &eq_max = 30	///< eq max to solve -> enum  FIELDS (equations_conf.h)
    );
  
  //! This function solves one step  for transient problems
  void solve_and_update (const int &t_in,          ///< initial time iteration
                      const int &t_step,        ///< actual time iteration
                      const int &print_step,    ///< print every
                      double &time,             ///< actual time
                      double &dt,               ///< step time
                      const int &eq_min = 0,    ///< eq min to solve -> enum  FIELDS (equations_conf.h) 
                      const int &eq_max = 30    ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );  
  
  //! This function solves one step  for transient problems
  void solve_onestep (const int &t_in,	        ///< initial time iteration
		      const int &t_step,	///< actual time iteration
		      const int &print_step,	///< print every
		      double &time,	        ///< actual time
		      double &dt,	        ///< step time
		      const int &eq_min = 0,	///< eq min to solve -> enum  FIELDS (equations_conf.h) 
		      const int &eq_max = 30	///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );

  //! This function solves the problem
  void dummy_step (const int &t_in,	        ///< initial time iteration
		   const int &t_step,	        ///< actual time iteration
		   const int &print_step,	///< print every
		   double &time,	        ///< actual time
		   double &dt	                ///< step time
    );

  //! This function solves one step  for transient problems
  void solve_non_linear_onestep (const int &t_in,	///< initial time iteration
				 const int &t_step,	///< actual time iteration
				 const int &print_step,	///< print every
				 double &time,	        ///< actual time
				 double &dt	        ///< step time
    );

  //! This function solves one step  for transient problems
  void solve_control_onestep (const int &t_in,	        ///< initial time iteration
			      const int &t_step,	///< actual time iteration
			      const int &print_step,	///< print every
			      double &time,	        ///< actual time
			      double &dt	        ///< step time
    );

  //! This function write solution to/from x_ooold vector
  void set_uooold (const int &flag,	                ///<  0 xold-> x_ooold   1 x_ooold-> xold
		   const double &toll,	                ///< tolerance
		   const double delta_t_step_in,	//   (in)
		   const int &eq_min,	                ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
		   const int &eq_max	                ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
    );

  EquationSystemsExtendedM & get_MGExtSystem ()
  {
    return *_mg_equations_map;
  };

  double System_functional (const int &ff,	///< initial time iteration
			    double parameter,	///< functional parameter
			    double &control	///< step control
    );
  //! This function sets the controlled domain
   void setCtrlDomain(
     const double xMin,
     const double xMax,
     const double yMin,
     const double yMax,
     const double zMin,
     const double zMax
   );





  
  
//=============================================================================
//                              FUNCTIONS USING
//                                    MED
//                                 LIBRARIES
//=============================================================================  
#ifdef HAVE_MED
  //! This function sets the Med-Mesh
  void setMedMesh ( );
  
  //! This function returns the Med-Mesh
  const ParaMEDMEM::MEDCouplingUMesh & getMedMesh () {return *_med_mesh;};
  
  //! Function getting node interpolated wall distance
  inline ParaMEDMEM::MEDCouplingFieldDouble * GetWallDistField (int Level)
  {
    MyAssert (_NodeWallDist.size() != 0,
	      "FEMUS::GetWallDistField NodeWallDist not calculated yet! ");
    return _NodeWallDist[Level];
  };

  //! Function for setting a MED field as external source field in <systemName> equation solver
  void setExtField (const std::string & systemName,
		    ParaMEDMEM::MEDCouplingFieldDouble * bcField);
  ParaMEDMEM::MEDCouplingFieldDouble * GetExtField (const std::string & systemName);
//=============================================================================
//   FUNCTIONS RELATED TO INTERFACE CREATION AND UPDATE
//=============================================================================
  
  //! This function sets up an interface with all group ids contained in IDSvec
  void init_interface (const int interface_name,              /**< Interface name */
                       std::vector<int> IDSvec,               /**< Vector containing group IDS */
		       int order_cmp,                         /**< Order of the numerical field */
		       const std::string & medfile_name,      /**< Name of the MED file containing the mesh */
		       bool on_nodes = true,                  /**< Flag to indicate if the field is ON_NODES or ON_CELLS */
		       const int index_medmesh = 0	      /**< med-mesh index  */
                       );

  //! This function sets up an interface for the whole mesh -> interface_id not required
  void init_interface (const int interface_name,             /**< Interface name */
		       int order_cmp,                        /**< Order (piecewise, linear, quadratic) */
		       const std::string & medfile_name,     /**< Name of the MED file containing the mesh */
		       bool on_nodes = true                  /**< Flag to indicate if the field is ON_NODES or ON_CELLS */
		      );
  
  //! This functions actually create the interface 
  /*!
   * The interface is created and added to #EquationSystemsExtendedM::_interfaceFunMap
   */
  void init_interface (const int interface_name,                        /**< Interface name */
		       int order_cmp,                                   /**< Order (piecewise, linear, quadratic) */
		       ParaMEDMEM::MEDCouplingUMesh * InterfaceMesh,    /**< MED interface mesh */
		       ParaMEDMEM::DataArrayDouble * MEDToMgMapArray,   /**< Array containing the map between MED and MG numeration */
		       int interface_id = 0                             /**< Default interface_id */
		      );
  
  void init_par_interface (int order_cmp,                    /**< Order (piecewise, linear, quadratic) */
                       bool on_nodes = true                  /**< Flag to indicate if the field is ON_NODES or ON_CELLS */
                      );  
    
  void update_interface (const int interface_name,	                                        /**< Interface name */
			 int n_cmp,	                                                        /**< Number of components of displacement field */
			 const std::vector < ParaMEDMEM::MEDCouplingFieldDouble * >&srcField	/**< Vector containing displacement MED field */
                        );

  //! This function sets the value from first_cmp to end_cmp of the field on the old solution x_old of the interface id. From Interface to System solution
  void write_Boundary_value (int id_boundary_name,	  ///< identity interface name (in)
			     std::string mgsystem_name,	  ///< system name          (in)
			     int n_cmp,	                  ///< from variable system (in)
			     int first_cmp = 0	          ///< to variable system   (in)
                            );

  //! This function sets the field of interface function (interface_id) with an analitic expression
  void setAnalyticSource (int interface_name,                   /**< Interface name */
			  int n_cmp,                            /**< Number of components of displacement field */
			  const std::string & bcExpression	/**< String containing analytical expression */
                         );

  //! This function sets the field of interface function (interface_id)
  void setFieldSource (int interface_name,                                    /**< Interface name */
		       int n_cmp,                                             /**< Number of components of displacement field */
		       const ParaMEDMEM::MEDCouplingFieldDouble * srcField    /**< MED field to set */
		      );  

  //! This function gets all the values on boundary with identity id
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary (int interface_name,	                ///< boundary name (char*) (in)
							    const std::string & systemName,	///< system name           (in)
							    int n_cmp,	                        ///< component             (in)
							    int first_cmp = 0	                ///< component             (in)
                                                           );

  ParaMEDMEM::MEDCouplingFieldDouble * getProcSolution (const std::string & system_name,	///< system name             (in)
							int n_cmp,	                        ///<  first variable         (in)
							int first_cmp = 0,	                ///< n variables             (in)
                            int Level = 0);

  void getNodeMapAndProcMeshAtLevel ( int level, ParaMEDMEM::MEDCouplingUMesh *&FemusPar, ParaMEDMEM::MEDCouplingFieldDouble *&NodeMap );
  
  //! This function gets actual support of the interface 
  const ParaMEDMEM::MEDCouplingUMesh * getUMesh (int name);
  
  //! This routin gets original support of the interface 
  const ParaMEDMEM::MEDCouplingUMesh * getUMesh_orig (int name);

  void GetInfo (string medfile_name,
		std::string & mesh_dir,
		std::string & localFile,
		std::string & filename,
		std::vector < std::string > &meshNames,
		std::vector < std::string > &MeshNames,
		std::vector < std::string > &FieldNames,
		const int index_medmesh = 0);

  void InitTurbulence ();
  void CalcTurbulence ();

#endif

};

#endif
