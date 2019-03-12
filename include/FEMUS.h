#ifndef __FEMUS__
#define __FEMUS__

#include <vector>
#include "Solverlib_conf.h"
#include "mpi.h"
#include "MGUtils.h"


#ifdef HAVE_MED
namespace MEDCoupling
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
    MGFemusInit                *_start;	              // start function
    MGUtils                    *_mg_utils;	      // Parameters and files
    MGGeomEl                   *_mg_geomel;	      // Element type
    MGFEMap                    *_mg_femap;	      // Finite Element Map
    MGTimeLoop                 *_mg_time_loop;	      // Time control
    MeshExtended               *_mg_mesh;	      // FEMus-mesh
    EquationSystemsExtendedM   *_mg_equations_map;    // Equation Map
    
    MPI_Comm _comm;		// Communicator
    
    bool _local_MPI_Init;		// initial mpi flag
    bool _MgEquationMapInitialized;
    bool _MgMeshInitialized;

#ifdef HAVE_MED
    MEDCoupling::MEDCouplingUMesh *         _med_mesh;     // Med-mesh
#endif
    // PUBLIC VARIABLES ===================================================================
public:

#ifdef HAVE_MED
    const int _GlobInterfaceId = 1234; // default interface id for interface covering entire domain
    int _ParallelInterfaceId;
#endif    
    
    std::vector<FIELDS> _myproblemP;
    
    // Constructor-Destructor
    FEMUS ();
    FEMUS (MPI_Comm comm);
    FEMUS (MGUtils & mgutils);
    ~FEMUS ();

    // Init structures
    //! This function initialize FEMUS class object
    void init_femus();
    //! This function initialize parameters with turbulence
    void init_param (MGUtils & mgutils, TurbUtils & Parameters, int name = 0);
    //! This function initialize parameters 
    void init_param (MGUtils & mgutils, int name = 0);
    //! This function initialize the finite element and functions
    void init_fem (MGGeomEl & mggeomel, MGFEMap & mgfemap);
    //! This function initialize the finite element and functions
    void init_fem ();
    //! This function initialize the map of equations to solve
    /*MGEquationsSystem &*/ void init_equation_system(int n_data_points = 0, int n_data_cell = 0);

    void terminate ();

#ifdef HAVE_MED
    //! MED field containing node interpolated wall distance
    std::vector<MEDCoupling::MEDCouplingFieldDouble *> _NodeWallDist;
    //! MED field containing node numbering map between med local parallel mesh and med global mesh
    std::vector<MEDCoupling::MEDCouplingFieldDouble *> _NodeMap;
#endif
#ifdef  TWO_PHASE_LIB
    void set_mgcc(MGSolCC & cc);
#endif
    /// MESH RELATED FUNCTIONS -------------------------------------------------
    //! Default setMesh function
    void set_mesh();
    //! setMesh function and computation of parallel mesh
    void setMeshTurbCase ();
    //! String for assigning a mesh name label
    std::string MeshName;
    
    
    //! Initialization of equations to solve
    void init_systems ();
    

    inline EquationSystemsExtendedM & get_MGExtSystem ()
    {
        return *_mg_equations_map;
    };
    
    //! Function assigning mesh name label
    inline void setMeshName (std::string meshname)
    {
        MeshName = meshname;
    };


    inline int IsVelFieldActive() {
        int a=0;
        if(stoi(_mg_utils->_sim_config["MG_NavierStokes"])!=0 ||  stoi(_mg_utils->_sim_config["MG_FluidStructure"])!=0) a = 1;
        return a;
    }


    //! Function getting mesh name label
    inline std::string getMeshName ()
    {
        return MeshName;
    };
    //! Function getting mesh dir
    inline std::string GetMeshDir ()
    {
        return _mg_utils->_mesh_dir;
    };
    //! Function returning _mg_mesh pointer
    const MeshExtended & get_MGMesh ()
    {
        return *_mg_mesh;
    };

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
    void solve_control_onestep (const int &nmax_step,     ///< number max of steps         (in)
                                const int &it,	        ///< iteration
                                const int &t_step,	///< actual time iteration
                                const int &print_step,	///< print every
                                double    &time,	        ///< actual time
                                double    &dt,	        ///< step time
                                const int &eq_min,        ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
                                const int &eq_max,        ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
                                bool      &converged      ///< check if the solution converged (1->converged) (out)
                               );
        //! This function solves one step  for transient problems
    void solve_control_onestep (const int               &nmax_step,         ///< number max of steps         (in)
                                const int               &it,	            ///< iteration
                                const int               &t_step,	        ///< actual time iteration
                                const int               &print_step,	    ///< print every
                                double                  &time,	            ///< actual time
                                double                  &dt,	            ///< step time
                                const int               &eq_min,            ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
                                const int               &eq_max,            ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
                                std::vector<double>     controlled_eq,
                                bool                    &converged,         ///< check if the solution converged (1->converged) (out)
                                const double            &toll = 1e-5        /// tolerance set by default 1e-5
                               );

    //! This function write solution to/from x_ooold vector
    void set_uooold (const int &vec_from,	        ///< source vector to be copied
                     const int &vec_to,	                ///< target vector 
                     const double &toll,	        ///< tolerance
                     const double delta_t_step_in,	//   (in)
                     const int &eq_min,	                ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
                     const int &eq_max	                ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
                    );

    

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
//=====================================================================================
#ifdef HAVE_MED
    //  MED-mesh function =============================================================
    //! This function sets the Med-Mesh
    void setMedMesh ( );
    //! This function returns the Med-Mesh
    const MEDCoupling::MEDCouplingUMesh & getMedMesh () {return *_med_mesh;};
    //! Function getting node interpolated wall distance
    inline MEDCoupling::MEDCouplingFieldDouble * GetWallDistField (int Level)
    { std::cout<< "FEMUS::GetWallDistField NodeWallDist not calculated yet! ";abort();
//         MyAssert (_NodeWallDist.size() != 0,
//                   "FEMUS::GetWallDistField NodeWallDist not calculated yet! ");
        return _NodeWallDist[Level];
    };

 
//=============================================================================
//   FUNCTIONS RELATED TO INTERFACE CREATION AND UPDATE
//=============================================================================
   // ========================================================================= 
   //  external functions 
     // ================================================================================================
    //! GetInfo med from file    medfile_name   
    void GetInfo (                                    
      std::string medfile_name,  ///<  medfile_name (in)
      std::string & mesh_dir,    ///< mesh_dir (out std::string )
      std::string & localFile,   ///< localFile  (out std::string )
      std::string & filename,    ///< filename   (out std::string  )
      std::vector < std::string > &meshNames,  ///< meshNames  (out std::vector<std::string>)
      std::vector < std::string > &MeshNames,  ///< MeshNames  (out std::vector<std::string>)
      std::vector < std::string > &FieldNames, ///< FieldNames  (out std::vector<std::string>)
      const int index_medmesh = 0
    );
    
// ================================================================================================
   void getNodeMapAndProcMeshAtLevel ( 
   int level, 
   MEDCoupling::MEDCouplingUMesh *&FemusPar, 
   MEDCoupling::MEDCouplingFieldDouble *&NodeMap 
   );   
    
   // ========================================================================= 
   //! Function for setting a MED field as external source field in <systemName> equation solver
    void setExtField (const std::string & systemName,
                      MEDCoupling::MEDCouplingFieldDouble * bcField);
    //! Function for getting a MED field as external source field in <systemName> equation solver
    MEDCoupling::MEDCouplingFieldDouble * GetExtField (const std::string & systemName);
    // ========================================================================= 
    //  interfaces functions 
    // ========================================================================= 
    // ================================================================================================
    //! This function sets up an interface with all group ids contained in IDSvec
    void init_interface (
      const int interface_name,              ///< Interface name 
      std::vector<int> IDSvec,               ///< Vector containing group IDS 
      int order_cmp,                         ///< Order of the numerical field 
      const std::string & medfile_name,      ///< Name of the MED file containing the mesh 
      bool on_nodes = true,                  ///< Flag to indicate if the field is ON_NODES or ON_CELLS 
      const int index_medmesh = 0	     ///< med-mesh index  
      );
    // ------------------------------------------------------------------------------------
    //! This function sets up an interface for the whole mesh -> interface_id not required
    void init_interface (
      const int interface_name,             ///< Interface name 
      int order_cmp,                        ///< Order (piecewise, linear, quadratic) 
      const std::string & medfile_name,     ///< Name of the MED file containing the mesh 
      bool on_nodes = true                  ///< Flag to indicate if the field is ON_NODES or ON_CELLS 
      );
    // ------------------------------------------------------------------------------------
    /// Basic interface construction with support (InterfaceMesh) and  MEDToMg map (MEDToMgMapArray) 
    /// The interface is created and added to #EquationSystemsExtendedM::_interfaceFunMap
    // ------------------------------------------------------------------------------------
    void init_interface (
      const int interface_name,                        ///< Interface name 
      int order_cmp,                                   ///< Order (piecewise, linear, quadratic) 
      MEDCoupling::MEDCouplingUMesh * InterfaceMesh,   ///< MED interface mesh 
      MEDCoupling::DataArrayDouble * MEDToMgMapArray,  ///< Array containing the map with MED/MG numeration 
      int interface_id = 0                             ///< Default interface_id 
     );
    // ------------------------------------------------------------------------------------
   /// This function gets the volume mesh from the med file as default
     // ------------------------------------------------------------------------------------
    void init_par_interface (
      int order_cmp,        ///< Order (piecewise, linear, quadratic) 
      bool on_nodes = true  ///< Flag to indicate if the field is ON_NODES or ON_CELLS 
     );
    
    
 
    // ================================================================================================
    //! This function sets the value from first_cmp to end_cmp of the field on the old solution x_old of the interface id. From Interface to System solution
    void write_Boundary_value (int id_boundary_name,	  ///< identity interface name (in)
                               std::string mgsystem_name,	  ///< system name          (in)
                               int n_cmp,	                  ///< from variable system (in)
                               int first_cmp = 0	          ///< to variable system   (in)
                              );
    // ================================================================================================
    //! This function sets the field of interface function (interface_id) with an analitic expression
    void setAnalyticSource (int interface_name,                   /**< Interface name */
                            int n_cmp,                            /**< Number of components of displacement field */
                            const std::string & bcExpression	/**< String containing analytical expression */
                           );
    // ================================================================================================
    //! This function sets the field of interface function (interface_id)
    void setFieldSource (
      int interface_name,                                  ///< Interface name 
      int n_cmp,                                           ///< Number of components of displacement field 
      const MEDCoupling::MEDCouplingFieldDouble * srcField ///< MED field to set 
                        );
    // ================================================================================================
    // ============================================================================
/// This function gets all the values on boundary with identity id
MEDCoupling::MEDCouplingFieldDouble *getValuesOnInterface_from_node(   // field (out)
  int  interface_name,           // boundary name (char*) (in)
  const std::string &systemName, // system name           (in)
  int n_cmp,                     // component             (in)
  int order=2,                   // return order quad=2 lin=1 (in)
  int first_cmp =0               // component             (in)
) ;
// ============================================================================
/// This function gets all the values on boundary with identity id
MEDCoupling::MEDCouplingFieldDouble *getValuesOnInterface_from_cell(   // field (out)
  int  interface_name,           // boundary name (char*) (in)
  const std::string &systemName, // system name           (in)
  int n_cmp,                     // component             (in) 
  int order =2,                  // return order quad=2 lin=1 (in)
  int first_cmp=0                // component             (in)
) ;

    
    
    
    //! This function gets all the values on boundary with identity id old !!!!!!!!!!!!!!!!!!!!
    MEDCoupling::MEDCouplingFieldDouble * getValuesOnInterface (
      int interface_name,	        ///< boundary name (char*) (in)
      const std::string & systemName,	///< system name           (in)
      int n_cmp,	                ///< component             (in)
      int first_cmp = 0	                ///< component             (in)
     );
    // ================================================================================================
    MEDCoupling::MEDCouplingFieldDouble * getProcSolution (
      const std::string & system_name,	///< system name             (in)
      int n_cmp,	                ///<  first variable         (in)
      int first_cmp = 0,	        ///< n variables             (in)
      int Level = 0
    );   
    
   // ==================================================================================
   // MESH movement
   // ==================================================================================
    // ================================================================================================                                                          );
    //! This function gets all the values on boundary with identity id
    MEDCoupling::MEDCouplingFieldDouble * getDisplacement (
      int interface_name,	        ///< boundary name (char*) (in)
      const std::string & systemName,	///< system name           (in)
      int n_cmp,	                ///< component             (in)
      int first_cmp = 0	                ///< component             (in)
       );
    // ================================================================================================
    /// This function moves the FEMUS interface according to a given displacement field
    void update_interface (
      const int interface_name,	                   ///< Interface name 
      int n_cmp,	                           ///< Number of components of displacement field
      const std::vector< MEDCoupling::MEDCouplingFieldDouble * >&srcField ///< Vector containing displacement MED field */
     );
 
 
 
 
// ================================================================================================
    //! This function gets actual support of the interface
    const MEDCoupling::MEDCouplingUMesh * getUMesh (
      int name
    );
    // ================================================================================================
    //! This routine gets original support of the interface
    const MEDCoupling::MEDCouplingUMesh * getUMesh_orig (
      int name
    );
   
    // ================================================================================================
    void InitTurbulence ();
    void InitTurbulence ( int MeshID );
    void CalcTurbulence ();

    // ================================================================================================
    double GetValue(
      const int & ff,
      int flag
    );
    void   SetValue(
      const int & ff,
      double value
    );
    void   SetValueVector(
      const int & ff,
      std::vector<double> value
    );
    // ================================================================================================
    void SetPieceFieldOnYdist(
      MEDCoupling::MEDCouplingFieldDouble * Field
    );
    
//     MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ReadField(
//        MEDCoupling::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order
//                                                   );
    

#endif

};

#endif
