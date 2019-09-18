#ifndef __mgsolverda_h__
#define __mgsolverda_h__

// conf includes
#include "MGFE_conf.h"
// algebra class
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "numeric_vectorM.h"

// local class
#include "MGFEMap.h"
#include "MGSolverBase.h"  //for the inherited class
// #include "MGFEMap.h"

// Forward declarations --------------------------
class MGUtils;
class MGSystem;
class MGMesh;
class MGEquationsMap;

// ================================================================================================
// ================================================================================================
/// This Class contains all the pointer table of the possible external fields for the system
class external_field
/// This class is based on enum FIELDS ( to define in Equation_conf.h)
/// FIELDS: equation symbol (NS_F) -> preassigned order (0)
/// enum  FIELDS{
///    NS_F  =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSX_F =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSY_F =1,     // [1] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSZ_F =2,     // [2] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
///  };
///
/// This class defines:
/// 1)  tab_eqs : FIELDS -> Eqs order (system appearance order)
///               (NS_F=0) ->  0
///               ( P_F=3) ->  1
/// 2)  indx_ub : Eqs order -> data pos in  ub
/// 3)  ub      : I -> data
/// The  tab_eqs + indx_ub are defined in the
///                    MGSolverDA function   set_ext_fields(const std::vector<FIELDS> &pbName)
/// The vector ub is filled during assemblying by
/// MGSolverDA function get_el_sol(#block data,KLQ,el_ndof,el_conn,offset,starting,_data_eq[1].ub)
/// Each  MGSolverXX class has a _data  external_field in its data
// ================================================================================================
{  // ===============================================================================================
 public:
  ///@{ \name EXTERNAL FIELDS INFORMATION
  int n_eqs;                ///< number of equations in the system
  const int max_neqs = 40;  ///< max number of equations (can be changed)
  MGSolBase* mg_eqs[40];    ///< equation system pointer (mg_eqs)
  int tab_eqs[40];          ///< map external system and index

  int indx_ub[40];           ///< index of external equations(const/linear/quad):
  double ub[40 * NDOF_FEM];  ///< element external field old solution= max_neqs*NDOF_FEM
  ///@}
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  external_field() {}   //< Empty Constructor
  ~external_field() {}  //< Empty destructor
                        ///@}
};
// =====================================================

// ====================================================
/// Class for MG Diffusion-Advection equation solvers
// ====================================================

class MGSolDA : public MGSolBase {
 public:
  // -------------------------------------------------------------
  // -------------------- class field ----------------------------

  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
  double _phi_g[3][NDOF_FEM];                            ///< shape field (0-1-2 degree)
  double _dphi_g[3][NDOF_FEM * DIMENSION];               ///< shape derivative  (0-1-2 degree)
  double _ddphi_g[3][NDOF_FEM * DIMENSION * DIMENSION];  ///< shape second derivative  (0-1-2 degree)

  //  external field
  external_field _data_eq[3];  ///< external data structure

  // data -----------------------------------------
  int _el_dof[3];  ///< number of dof in each element  (piecewise, linear, quadratic)
  MGFE* _fe[3];    ///< fem  (piecewise linear,piecewise quadratic)
  double _dt = 0.1;
  int _dir = 0;

  //! Number of solutions needed for restart
  /* This number is used for printing and reading solutions for system restart.
   * Default value is equal to 1. With value _NumRestartSol=2 the solution
   * stored in _x_oold vector is printed with _old suffix, while with _NumRestartSol=3
   * also the _x_ooold stored solution is printed, with _oold suffix (visible inside .h5 file)
   */
  int _NumRestartSol;

  // ========================================================================
  ///@{ \name  Constructor destructor memory allocation
  // ========================================================================
  /// Level constructor  I
  MGSolDA(
      MGEquationsSystem& mg_equations_map,
      const int nvars_in[],  ///< number of piecewise[0], linear[1], quadratic[2] variables
      std::string eq_name_in = "DA", std::string varname_in = "u");

  /// Destructor (level structure)
  ~MGSolDA();
  /// Clean all substructures
  void clean();
  // Setting --------------------------------------------------
  virtual void init_dof(  ///< Setting dof
      const int Level,    // MG Level
      const int vb_0 = 0, const int n_vb = 1);
  virtual void init(   ///< Setting Memmory alloc
      const int Level  // MG Level
  );
  ///@}
  // ========================================================================
  ///@{ \name BOUNDARY/INITIAL CONDITIONS
  // ========================================================================
  void GenBc();  ///< Setting boundary conditions
  void GenBc_loop(
      const int vb, const int ndof_femv, const int n_ub_dofs, const int n_pb_dofs, const int n_kb_dofs,
      int bc_id[], int mat_id[]);
  virtual void GenIc();  ///< Element Initial conditions
                         //   virtual void GenIc0(  ///< Setting Initial conditions
                         //     const std::string ib_name,
                         //     const int node_dof_top[],
                         //     NumericVectorM &sol_top,
                         //     NumericVectorM &old_sol_top
                         //   );
                         ///@}

  // ========================================================================
  ///@{ \name MULTILEVEL OPERATORS
  // ========================================================================
  // Reading/Writing MG operators ------------------------------------------
  void ReadProl(                ///< Read Prolongation Op
      const int Level,          // MG Level <-
      const std::string& name,  // file name for P <-
      SparseMMatrixM& Mat,      //  Matrix for P <-
      const int nvars_in[], int node_dof_c[], int node_dof_f[]);

  void ReadRest(                ///< Restriction Op.
      const int Level,          // MG Level
      const std::string& name,  // file name (reading from)
      SparseMMatrixM& Rest,     // Restriction Matrix
      const int nvars_in[],     // # cost,linear quad variables
      int node_dof_f[],         // dof map fine mesh
      int node_dof_c[],         // dof map coarse
      int _node_dof_top[]       // dof map top level
  );                            // ----

  void ReadMatrix(              ///< Reading Matrix
      const int Level,          // MG Level <-
      const std::string& name,  // file name for M <-
      SparseMatrixM& Mat, const int* nvars_in);
  ///@}

  // ========================================================================
  ///@{ \name    MultiGrid Solver
  // ========================================================================
  virtual void GenMatRhs(  ///< Volume Assemblying  matrix-rhs
      const double time,   // time
      const int Lev,       // Level
      const int m          // rhs assembly control
  );
  /// MG time step solver (backward Euler).
  virtual void MGTimeStep(
      const double time,  // time               <-
      const int mode      // rhs assembler flag <-
  );

  ///@}

  // ========================================================================
  ///@{ \name SOLUTION WRITE AND READ
  // ========================================================================
  // Print
  virtual void print_bc(     ///< Print boundary conditions to a xdmf file.
      std::string namefile,  // filename <-
      const int Level        // MGLevel  <-
  );
  virtual void print_u(      /// Print solution to a xdmf file.
      std::string namefile,  // filename <-
      const int Level);      // MGLevel  <-

#ifdef HAVE_MED
  virtual void print_u_med(  /// Print solution to a med file.
      std::string namefile,  // filename <-
      const int Level);      // MGLevel  <-
#endif
#ifdef HAVE_MED
  virtual void print_weight_med(  /// Print weight for control problems to a med file.
      std::string namefile,       // filename <-
      const int Level);           // MGLevel  <-
#endif
  virtual void print_xml_attrib(
      std::ofstream& out,  // file stream to print ->
      int nodes, int nelems, std::string file_name) const;
  // Read
  virtual void read_u(   /// Read solution from a xdmf file
      std::string name,  // filename <-
      int Level          // restart level <-
  );

  // ============================================================================
  /// This function  defines the boundary conditions for the system:
  virtual void bc_intern_read(
      int /*face_id_node*/,  ///<  face identity           (in)
      int /*mat_flag*/,      ///<  volume identity         (in)
      double /*xp*/[],       ///< xp[] node coordinates    (in)
      int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
      int bc_flag[]          ///< boundary condition flag  (out)
  );
  // ============================================================================
  /// This function reads Boundary conditions  from function
  virtual void bc_read(
      int /*face_id_node*/,  ///<  face identity          (in)
      int /*mat_flag*/,      ///<  volume identity         (in)
      double /*xp*/[],       ///< xp[] node coordinates    (in)
      int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
      int bc_flag[]          ///< boundary condition flag  (out)
  );

  virtual void ic_read(  /// initial conditions  from function
      int k,             // bc from gambit    <-
      int m,             // material from gambit    <-
      double xp[],       // point coordinates <-
      int iel,           // element           <-
      double value[]     //  point values     ->
  );
  virtual void set_vector(const int&, const int&);
  virtual void set_xooold2x();
  virtual double MGFunctional(double, double&);
  ///@}

  // Return ----------------------------------------------------
  // ========================================================================
  ///@{ \name INTERPOLATIONS FUNCTIONS
  // ========================================================================
  // void get_el_lq(const uint Level,const int nvars[],const int el_nds[],
  // 			 const uint el_conn[],const uint offset,
  // 			 std::vector<uint>   & el_dof_indices,
  //                           std::vector<uint>   bc_bd[],  std::vector<uint>   bc_vol[],
  //                          std::vector<double> uold[])  const;
  // void  get_el_lq(
  //   const uint Level,       // level
  //   const int nvars[],       // # of variables to get  <-
  //   const int el_nds[],      // # of element nodes for this variable  <-
  //   const uint el_conn[],   // connectivity <-
  //   const uint offset,      // offset for connectivity <-
  //   std::vector<uint>   & el_dof_indices, // element connectivity ->
  //   uint  bc_vol[],        // element boundary cond flags ->
  //   uint bc_bd[],        // element boundary cond flags ->
  //   double uold[]           // element node values ->
  // )  const;

  void interp_el_sol(
      const double uold_b[],  // node values <-
      const int ivar0,        // init variable  <-
      const int nvars,        // # of variables  <-
      const double phi[],     // shape functions  <-
      const int n_shape,      // # of shape functions  <-
      double uold[]           // interpolated function ->
      ) const;

  void interp_el_bd_sol(
      const double uold_b[],  // node values <-
      const int sur_tpgly[],  // surface nodes topology <-
      const int el_ndof,      // surface nodes topology <-
      const int ivar0,        // init variable  <-
      const int nvars,        // # of variables  <-
      const double phi[],     // shape functions  <-
      const int n_shape,      // # of shape functions  <-
      double uold[]           // interpolated function ->
      ) const;                // =======================================

  void interp_el_gdx(
      double uold_b[],      // node values <-
      const int ivar0,      // init variable  <-
      const int nvars,      // # of variables  <-
      const double dphi[],  // derivatives of the shape functions  <-
      const int n_shape,    // # of shape functions  <-
      double uold_dx[]      // interpolated derivatives ->
      ) const;

  void interp_el_bd_gdx(
      const double uold_b[],  // node values <-
      const int sur_tpgly[],  // surface nodes topology <-
      const int el_ndof,      // surface nodes topology <-
      const int ivar0,        // init variable  <-
      const int nvars,        // # of variables  <-
      const double dphi[],    // derivatives of the shape functions  <-
      const int n_shape,      // # of shape functions  <-
      double uold_dx[]        // interpolated derivatives ->
      ) const;                // =======================================

  void interp_el_gddx(
      double uold_b[],      // node values <-
      const int ivar0,      // init variable  <-
      const int nvars,      // # of variables  <-
      const double dphi[],  // derivatives of the shape functions  <-
      const int n_shape,    // # of shape functions  <-
      double uold_dx[]      // interpolated derivatives ->
      ) const;

  // void  interp_el_dx (
  //   double uold_b[], // node values <-
  //   const uint ivar0,      // init variable  <-
  //   const uint nvars,      // # of variables  <-
  //   const double phi[],    // shape functions  <-
  //   const double dphi[],   // derivatives of the shape functions  <-
  //   const uint n_shape,    // # of shape functions  <-
  //   double uold[],         // interpolated function ->
  //   double uold_dx[]       // interpolated derivatives ->
  // )  const;

  void compute_jac(
      const int j, const int idim,
      double uold_b[],      // node values <-
      const int nvars,      // # of variables  <-
      const double phi[],   // shape functions  <-
      const double dphi[],  // derivatives of the shape functions  <-
      const int n_shape,    // # of shape functions  <-
      double u_forw[],      // interpolated function ->
      double u_back[],      // interpolated function ->
      double u_forw_dx[],   // interpolated derivatives ->
      double u_back_dx[]    // interpolated derivatives ->

      ) const;
  ///@}

  // ========================================================================
  ///@{ \name RETURN FUNCTIONS
  // ========================================================================
  void get_el_sol(
      const int ivar0,      // initial variable  <-
      const int nvars,      // # of variables to get  <-
      const int el_nds,     // # of element nodes for this variable  <-
      const int el_conn[],  // connectivity <-
      const int offset,     // offset for connectivity <-
      const int kvar0,      // offset  variable for  uold <-
      double uold[]         // element node values ->
      ) const;
  // virtual void  f_aux (double a1[]);

  //  void  get_el_q(const uint Level,const uint nvars, const uint el_nds,
  // 		const uint el_conn[], const uint offset,
  // 		std::vector<uint>& el_dof_indices,std::vector<uint> & bc_vol,
  // 		std::vector<uint>   & bc_bd,std::vector<double> & uold) const;

  /// Dof , the bc and the solution  vector at the nodes of  an element
  void get_el(
      const int Level,                   // level <-
      const int nvar0,                   // intial varables number based on NDOF_FEM <-
      const int nvar,                    // final varables number based on NDOF_FEM <-
      const int el_nds,                  // number of nodes for element
      const int el_conn[],               // connectivity <-
      const int offset,                  // offset <-
      std::vector<int>& el_dof_indices,  // DOF indices ->
      int bc_dofs[][NDOF_FEM],           // boudary conditions ->
      double uold[]                      // solution ->
      ) const;

  /// Dof , the bc and the solution  vector at the nodes of  an element
  /*  void  get_el(const uint Level,const uint el_off0,const uint nvar0,
             const uint nvar,const uint el_nds,
                 const uint el_conn[],const uint offset,
                 std::vector<uint>  &el_dof_indices,
             std::vector<uint>   bc_dofs[],
             std::vector<double>  & uold)  const;*/
  void get_el_dof_bc(
      const int Level,  // level
      const int iel,
      //   const int nvars[],       // # of variables to get  <-
      const int el_nds[],                // # of element nodes for this variable  <-
      const int el_conn[],               // connectivity <-
      const int offset,                  // offset for connectivity <-
      std::vector<int>& el_dof_indices,  // element connectivity ->
      int bc_vol[],                      // element boundary cond flags ->
      int bc_bd[]                        // element boundary cond flags ->
      ) const;                           // ==============================================================

  void set_el_dof_bc(
      const int Level,  // level
      const int iel,
      //   const int nvars[],       // # of variables to get  <-
      const int el_nds[],                // # of element nodes for this variable  <-
      const int el_conn[],               // connectivity <-
      const int offset,                  // offset for connectivity <-
      std::vector<int>& el_dof_indices,  // element connectivity ->
      int bc_vol[],                      // element boundary cond flags ->
      int bc_bd[]                        // element boundary cond flags ->
      ) const;                           // ==============================================================

  ///@}

  void set_dt(double dt) { _dt = dt; };  // computation val3=musker
  // ========================================================================
  ///@{ \name EXTERNAL FIELDS
  // ========================================================================
  //     void set_ext_fields(const std::vector<FIELDS> &pbName);
  virtual double eval_var1(double[]) { return 0.; };
  virtual double eval_var2(double[]) { return 0.; };
  virtual double eval_var3(double[]) { return 0.; };

  void setUpExtFieldData();
  void ActivateVectField(int Order, int Field, std::string SystemFieldName, int& n_index, int coupled);
  void ActivateScalar(int Order, int Field, std::string SystemFieldName, int& n_index);
  void ActivateControl(
      int Order, int Field, std::string SystemFieldName, int& n_index, int vector, int dimension);
  void ActivateCoupled(
      int Order, int Field, std::string SystemFieldName, std::string SystemFieldName2, int& n_index);
  void ActivateEquation(int Order, int Field, std::string SystemFieldName, int& n_index);
  ///@}  double _xxg[DIMENSION];
};

#endif

// #endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
