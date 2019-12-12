#ifndef __mgsolbase__
#define __mgsolbase__

// c++ class include files ----------------------------------------------------

// algebra  matrix-vector class files (contrib) -------------------------------
class SparseMatrixM;   // algebra sparse matrices
class SparseMMatrixM;  // algebra sparse quadrangolar matrices
class NumericVectorM;  // algebra numerical vector
class LinearSolverM;   // algebra linear solver

// FEMus configure include files ----------------------------------------------
#include "Printinfo_conf.h"  //
#include "Solverlib_conf.h"
#include "MGFE_conf.h"

// FEMus  include files -------------------------------------------------------
#include "EquationsMap.h"
class MGUtils;
class MGSystem;
class MGEquationsSystem;
class MGFEMap;
class MGMesh;

#ifdef HAVE_MED  // ---------------------------------------------------------
#include "MEDCouplingFieldDouble.hxx"
#endif
#ifdef TWO_PHASE_LIB
class MGSolCC;
#endif

// ================================================================================================
// ================================================================================================
/// This Class contains all the pointer table of the possible external fields for the system
class external_field
/// This class is based on enum FIELDS ( to define in Equation_conf.h)
/// FIELDS: equation symbol (*_F) -> preassigned order (0)
/// enum  FIELDS{
///    NS_F  =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSX_F =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSY_F =1,     // [1] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    NSZ_F =2,     // [2] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
///    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
///    .......................
///    ........................ see  Equation_conf.h
///  };
/// This class defines:
/// 1)  tab_eqs : FIELDS -> Eqs order (system appearance order)  (.)_F=0) ->  0
/// 2)  indx_ub : Eqs order -> data pos in  ub
/// 3)  ub      : I -> data
/// The  tab_eqs + indx_ub are defined in the
///                    MGSolverDA function   set_ext_fields(const std::vector<FIELDS> &pbName)
/// The vector ub is filled during assemblying by
/// MGSolverDA function get_el_sol(#block data,KLQ,el_ndof,el_conn,offset,starting,_data_eq[1].ub)
/// Each  MGSolverXX class has a _data  external_field in its data
// ================================================================================================
{
  // ===============================================================================================
 public:
  ///@{ \name EXTERNAL FIELDS INFORMATION
  int n_eqs;                 ///< number of equations in the system
  const int max_neqs = 40;   ///< max number of equations (can be changed)
  MGSolBase* mg_eqs[40];     ///< equation system pointer (mg_eqs)
  int tab_eqs[40];           ///< map external system and index
  int indx_ub[40];           ///< index of external equations(const/linear/quad):
  double ub[40 * NDOF_FEM];  ///< element external field old solution= max_neqs*NDOF_FEM

  /// CONSTRUCTOR-DESTRUCTOR
  external_field() {}   //< Empty Constructor
  ~external_field() {}  //< Empty destructor
  ///@}
};
// =====================================================

//   #define VANKA (0)

// ****************************************************************************
/// This class is the basic equation class
class MGSolBase
#ifdef LM_REFCOUNT
    : public ReferenceCountedObject<MGSolBase>
#endif
{
  // **************************************************************************
 protected:
  // --------------------------------------------------------------------------
  //   DOF MAP
  int* _Dim;  ///< dimension number of dofs per level
  //  -------------------------------------------------------------------------
  // DATA POINTER
  MGEquationsSystem& _mgeqnmap;  ///<  equation map  pointerz
  MGUtils& _mgutils;             ///<  utility class pointer
  MGFEMap& _mgfemap;             ///<  FEM class pointer
  MGMesh& _mgmesh;               ///<  mesh pointer

#ifdef TWO_PHASE_LIB
  MGSolCC* _msolcc;
#endif
  //  -------------------------------------------------------------------------
  // DATA PARALLEL
  int _iproc;               ///< processor
  LinearSolverM** _solver;  ///< linear system solver type (each level)

  //  -------------------------------------------------------------------------
  //  MULTIGRID MATRICES AND RHS'S
  std::vector<SparseMatrixM*> A;     ///< Matrix A
  std::vector<NumericVectorM*> b;    ///< rhs b
  std::vector<NumericVectorM*> x;    ///< solution x
  std::vector<NumericVectorM*> res;  ///< residual

  // --------------------------------------------------------------------------
  //   MULTIGRID OPERATORS
  std::vector<SparseMMatrixM*> Rst;  ///< Restrictor
  std::vector<SparseMMatrixM*> Prl;  ///< Prolongation

  // --------------------------------------------------------------------------
 public:
  // --------------------------------------------------------------------------
  const int _NoLevels;  ///< level number
  int** _node_dof;      ///< dof map

  double _dt;  ///< time step
  int _dir;    ///< dir= eq index for vector system of equations

  double _control;                   ///< control flag
  std::vector<double> _weight_ctrl;  ///< controlled region for optimal control
  // --------------------------------------------------------------------------
  // LABELING
  const std::string _eqname;  ///< equation name
  // --------------------------------------------------------------------------
  //    VARIABLES
  const int _n_vars;        ///< number of variables
  int _nvars[3];            ///< number of variables quadratic[2], linear[1] and piecewise [0]
  std::string* _var_names;  ///< variable names
  double* _refvalue;        /// reference values

  //  ----------------------------------------------------------------------------------------
  // OLD VECTOR see SOLUTION SECTION  MGSolBase_SOL.C  (see section  below)
  //   std::vector<NumericVectorM*> x_old[3];  ///< Multi level [l=0..NoL] old solution x old x[0] oold x[1]
  //   ooold x[2] std::vector<NumericVectorM*> x_nonl;    ///<  Multi level [l=0..NoL] non linear solution x
  //   std::vector<NumericVectorM*> x_aux;     ///< vector for multiple uses on top
  //   std::vector<NumericVectorM*> d_aux;     ///< vector for multiple uses
  //   std::vector<double> _weight_ctrl;      ///< controlled region for optimal control problems
  // -----------------------------------------------------------------------------------------
  // EXTERNAL FIELD to rhe MGSolverBase  (see section  below)
  //   external_field _data_eq[3];  ///< external data structure
  // -----------------------------------------------------------------------------------------
  // BOUNDARY CONDITIONS see MGSolDA_BCIC.C  (see section  below)
  // int  *_bc[2];  ///< boundary conditions map (top level) see section below

  // ========================================================================
  // CONSTRUCTOR-DESTRUCTOR
  // ========================================================================
  MGSolBase(
      MGEquationsSystem& mg_equations_map,  ///< \param[in] <> MG equation map
      const int nvars_in[],                 ///< \param[in] <> number of variables
      std::string eq_name_in = "Base"       ///< \param[in] <> equation name
  );                                        ///< Level constructor
  MGSolBase(
      MGUtils& mgutils_in, MGFEMap& mgfemap_in, MGMesh& mgmesh_in,
      const int nvars_in[],  // # variables
      std::string eqname_in  // equation name
  );
  //-------------------------------------------------------------------------
  ~MGSolBase();  ///< Destructor (level structure)
  void clear();  ///< Substructure destructor
  //-------------------------------------------------------------------------
  /// Distributing dof function (virtual MGSolDA)
  virtual void init_dof(
      const int Level,     ///< \param[in]
      const int vb_0 = 0,  ///< type of mesh
      const int n_vb = 1   ///< number type of meshes
      ) = 0;
  //-------------------------------------------------------------------------
  ///  Initializing memory allocation (virtual MGSolDA)
  virtual void init(const int Level  ///<  Level
                    ) = 0;
  // =======================================================================
// set functions
// =======================================================================
#ifdef TWO_PHASE_LIB
  void set_mgcc(MGSolCC& cc);
#endif
  //-------------------------------------------------------------------------
  virtual void set_dt(double dt);  ///< MG time step solver (backward Euler)
                                   //-------------------------------------------------------------------------
  void set_ctrl_dom(
      const double x_min, const double x_max, const double y_min, const double y_max, const double z_min,
      const double z_max);

  // *******************************************************************************
  //  SOLUTION OPERATIONS ON SOLUTIONS
  //  defined in MGSolBase_SOL.C
  // *******************************************************************************
  // =======================================================================
  std::vector<NumericVectorM*>
      x_old[3];  ///< Multi level [l=0..NoL] old solution x old x[0] oold x[1]  ooold x[2]
  std::vector<NumericVectorM*> x_nonl;  ///<  Multi level [l=0..NoL] non linear solution x
  std::vector<NumericVectorM*> x_aux;   ///< vector for multiple uses on top
  std::vector<NumericVectorM*> d_aux;   ///< vector for multiple uses

  // ========================================================================
  //  SET FUNCTIONS
  // ========================================================================
  ///@{ \name RETURN FUNCTIONS
  ///< \param[in]   <ivar0>   initial variable
  ///< \param[in]   <nvars>   number of variables to get
  ///< \param[in]   <el_nds>  number  of element nodes for this variable
  ///< \param[in]  <el_conn> connectivity
  ///< \param[in]  <offset>  offset for connectivity
  ///< \param[in]  <kvar0>   offset  variable for  uold
  ///< \param[out]  <uold>   solution
  // ----------------------------------------------------------------

  // defined in MGSolBase
  void set_sol(int i, int kdofs, double value);
  void set_x_aux(int i, int kdofs, double value);
  void set_d_aux(int i, int kdofs, double value);

  //  -------------------------------------------------------------
  virtual void set_xooold2x() = 0;
  //  -------------------------------------------------------------
  virtual void SetValue(double value);  ///< This function sets a value in the solvers
  //  -------------------------------------------------------------
  virtual void SetValueVector(std::vector<double> value);  ///< This function sets a value in the solvers
  // ========================================================================
  //  GET FUNCTIONS
  // ========================================================================
  ///@{ \name RETURN FUNCTIONS
  ///< \param[in]   <ivar0>   initial variable
  ///< \param[in]   <nvars>   number of variables to get
  ///< \param[in]   <el_nds>  number  of element nodes for this variable
  ///< \param[in]  <el_conn> connectivity
  ///< \param[in]  <offset>  offset for connectivity
  ///< \param[in]  <kvar0>   offset  variable for  uold
  ///< \param[out]  <uold>   solution
  // ----------------------------------------------------------------
  double get_sol(int i, int kdofs);
  double get_x_aux(int i, int kdofs);
  double get_d_aux(int i, int kdofs);

  // ----------------------------------------------------------------
  void get_el_sol(
      const int i_step,     ///<  i-step
      const int ivar0,      ///< \param[in]   <ivar0>   initial variable
      const int nvars,      ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,     ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],  ///< \param[in]  <el_conn> connectivity
      const int offset,     ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,      ///< \param[in]  <kvar0>   offset  variable for  uold
      double uold[]         ///< \param[out]  <uold>   solution
      ///<
      ) const;  ///< Return element solution
  // ----------------------------------------------------------------
  void get_el_sol_piece(
      const int ivar0,   ///< \param[in]   <ivar0>   initial variable
      const int nvars,   ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,  ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int iel,     ///< \param[in]  <iel> element id (ofset)
      const int offset,  ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,   ///< \param[in]  <kvar0>   offset  variable for  uold
      double uold[]      ///< \param[out]  <uold>   solution
      ///<
      ) const;  ///< Return element solution
  // --------------------------------------
  /// Return no linear  solution on element
  void get_el_nonl_sol(
      const int ivar0,      ///< \param[in]   <ivar0>   initial variable
      const int nvars,      ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,     ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],  ///< \param[in]  <el_conn> connectivity
      const int offset,     ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,      ///< \param[in]  <kvar0>   offset  variable for  uold
      double uold[]         ///< \param[out]  <uold>   solution
      ) const;
  // --------------------------------------
  /// Return d_aux  solution on element
  void get_el_d_aux(
      const int istep,
      const int ivar0,      // initial variable  <-
      const int nvars,      // # of variables to get  <-
      const int el_nds,     // # of element nodes for this variable  <-
      const int el_conn[],  // connectivity <-
      const int offset,     // offset for connectivity <-
      const int kvar0,      // offset  variable for  uold <-
      double uold[]         // element node values ->
      ) const;
  // -------------------------------------------------------------------------------
  /// Return element dof indices using local-to-global map
  void get_el_dof_indices(
      const int Level,      ///< \param[in] <Level>   level
      const int iel,        ///< \param[in] <iel>     eLement number
      const int el_conn[],  ///< \param[in] <el_conn> connectivity
      const int el_dof[],   ///< \param[in] <el_dof>  quadratic[2] linear[1] const[0] dofs
      const int offset,     ///< \param[in] <offset>  offset for connectivity
      std::map<int, std::vector<int>>& el_dof_indices  ///< \param[out]<el_dof-indices>  dof indices
      ) const;
  // -----------------------------------------------------------------------------------
  virtual double GetValue(int flag);  ///< This function returns a value from the solvers

  //  ------------------------------------------------------
  /// This function interpolates a vector field over the fem element
  void interp_el(
      const double uold[],    // node values <-
      const int ivar0,        // init variable  <-
      const int nvars,        // # of variables  <-
      const double phi[],     // shape functions  <-
      const int n_shape,      // # of shape functions  <-
      double u_int[],         // interpolated function ->
      int dim,                // deriv or
      const int sur_tpgly[],  // surface nodes topology <-
      const int el_ndof       // surface nodes topology <-
      ) const;

  // -------------------------------------------------------------------------
  // This function interpolates the solution at gaussian point (phi[])
  void interp_el_sol(
      const double uold_b[],  // node values <-
      const int ivar0,        // init variable  <-
      const int nvars,        // # of variables  <-
      const double phi[],     // shape functions  <-
      const int n_shape,      // # of shape functions  <-
      double uold[]           // interpolated function ->
      ) const;
  // -------------------------------------------------------------------------
  // This function interpolates the 1derivative  at gaussian point (dphi[])
  void interp_el_gdx(
      double uold_b[],      // node values <-
      const int ivar0,      // init variable  <-
      const int nvars,      // # of variables  <-
      const double dphi[],  // derivatives of the shape functions  <-
      const int n_shape,    // # of shape functions  <-
      double uold_dx[]      // interpolated derivatives ->
      ) const;
  // -------------------------------------------------------------------------
  // This function interpolates the  2derivative at gaussian point (dphi[])
  void interp_el_gddx(
      double uold_b[],      // node values <-
      const int ivar0,      // init variable  <-
      const int nvars,      // # of variables  <-
      const double dphi[],  // derivatives of the shape functions  <-
      const int n_shape,    // # of shape functions  <-
      double uold_dx[]      // interpolated derivatives ->
      ) const;
  // -------------------------------------------------------------------------
  // This function interpolates the  1derivative
  // on the boundary (-> sur_tpgly[]) at gaussian point (dphi[])
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
  // -------------------------------------------------------------------------
  // This function interpolates the  solution
  // on the boundary (-> sur_tpgly[]) at gaussian point (dphi[])
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
                              // -------------------------------------------------------------------------
  // This function copy  the sol x_old[i1_old_sol] into   x_old[i2_old_sol]
  // on the boundary (-> sur_tpgly[]) at gaussian point (dphi[])
  virtual void set_cp_vector(const int& i1_old_sol, const int& i2_old_sol);
  // -----------------------------------------------------------

  void localize_xooold();

  // ************************************************************************************
  // BOUNDARY CONDITION and INITIAL CONDITIONS
  // all virtual defined in MGSOLDA_BCIC.C
  // ************************************************************************************

  // Data ===============================================================================
  // Boundary Condition Data in the vector *_bc[2]
  // _bc[0]=volume
  //
  int* _bc[2];  ///< boundary conditions map (top level)
                // set get ============================================================================

  // Functions =========================================================================
  // reading print ic
  virtual void GenIc() = 0;  ///< Reading Initial conditions (IC)
  // ------------------------------------------------------------------------------------
  virtual void ic_read(  /// initial conditions  from function
      int k,             // bc from gambit    <-
      int m,             // material from gambit    <-
      double xp[],       // point coordinates <-
      int iel,           // element  <-
      double valueic[]   //  point values     ->
      ) = 0;             ///< Reading IC  function (element)
              // reading print bc ==================================================================
  virtual void
  GenBc() = 0;  ///< Reading boundary conditions (BC)
                // ----------------------------------------------------------------------------------
  virtual void bc_read(
      int k,
      int m,         ///< global node id
      double x[],    ///< point vector
      int bc_vol[],  ///< Volume flag
      int bc_sur[]   ///< Surface flag
      ) = 0;         ///< Reading sur BC function (element)
  // ------------------------------------------------------------------------------------
  /// Reading vol BC function (element)
  virtual void bc_intern_read(
      int k, int m,
      double x[],    // point vector
      int bc_vol[],  // Volume flag
      int u[]        // value vector
      ) = 0;
  // -------------------------------------------------------------------------------------
  virtual void print_bc(
      std::string fname,  ///< \param[in] <>  file name
      const int Level     ///< \param[in] <> MG level
      ) = 0;              ///< print boundary conditions

  // ====================================================================================
  // Print/Read function (All virtual in MGSOLDA_PR.C
  // ====================================================================================
  // ------------------------------------------------------------------------------------
  //       Print/Read soution
  // ------------------------------------------------------------------------------------
  // void print_ext_data() is not implemented (to do)
  // void print_xml_attrib() is defined
  // the other are defined in MGSolverDA
  // Print
  // fname=  file name
  // Level=  MG level
  // of_out=  out stream file
  // n_nodes=  number of nodes
  // n_elems=  number of elements
  // Data ==============================================================================

  // Functions =========================================================================
  // --------------------------------------------------------------------
  virtual void print_ext_data(double /*vect_data*/
                                  []){}
  /*printf("\n \n Wrong use of function print_ext_data in SolverBase for coupled mesh \n \n") =0*/
  ;
  // --------------------------------------------------------------------
  virtual void print_u_xdmf(
      std::ofstream& of_out,  ///< \param[in] <>  out stream file
      int n_nodes,            ///< \param[in] <>  number of nodes
      int n_elems,            ///< \param[in] <>  number of elements
      std::string fname       ///< \param[in] <>  file name
      ) const;                ///< print xml file
  // --------------------------------------------------------------------
  virtual void print_u(
      std::string fname,  ///< \param[in] <>  file name
      const int Level     ///< \param[in] <> MG level
      ) = 0;              ///< print solution
  // --------------------------------------------------------------------
  virtual void read_u(
      std::string fname,  ///< \param[in] <>  file name
      int Level           ///< \param[in] <>  restart MG level
      ) = 0;              ///< Read solution

  // ---------------------------------------------------------------------
  //       Print/Read    multilevel operator
  //------------------------------------------------------------------------
  //  MGDofBcOp() is defined also in  MGSolverBase
  //  all the other are defined in MGSolverDA_PR.C
  // --------------------------------------------------------------------
  // This function read the Operators (defined in the MGSolverDA)
  virtual void MGDofBcOp();
  // This function read the Matrix structure  (defined in the MGSolverDA)
  virtual void ReadMatrix(  ///< Reading matrix A
      const int Level,      //  MG level <-
      const std::string& name, SparseMatrixM& Mat, const int* nvars_all) = 0;
  // This function read the Prolongation Operator   (defined in the MGSolverDA)
  virtual void ReadProl(  ///< Reading Prolongation
      const int Level,    //  MG level <-
      const std::string& name, SparseMMatrixM& Mat, const int nvars_in[], int node_dof_c[],
      int node_dof_f[]) = 0;
  // This function read the Restriction Operator   (defined in the MGSolverDA)
  virtual void ReadRest(
      const int Level,  //  MG level <-
      const std::string& name, SparseMMatrixM& Mat, const int nvars[], int node_dof_c[], int node_dof_f[],
      int _node_dof_top[]) = 0;

  // ==================================================
  ///  MULTILEVEL SOLUTION/ASSEMBLYING  (MGSolverBase.C)
  // ==================================================
  // Data =========================================================================================
  // Functions ===================================================================================
  // --------------------------------------------------------------------
  /// This function performes the V,W and Fast MultiGrid scheme
  virtual void MGSolve(
      double Eps,   ///< tolerance
      int MaxIter,  ///< n iterations
      // -------------------------------------------------------------
      const int clearing = 0,    ///< 1= clean the init matrix and precond
      const int Gamma = 1,       ///< Control V W cycle
      const int Nc_pre = 8,      ///< n pre-smoothing cycles
      const int Nc_coarse = 40,  ///< n coarse cycles
      const int Nc_post = 8      ///< n post-smoothing cycles
  );
  // ------------------------------------------------------------------
  /// This function solves a MultiGrid step
  virtual double MGStep(
      int Level,            ///<  MG level <-
      double Eps1,          ///< tolerance
      int MaxIter,          ///< n iteratio
      const int Gamma,      ///< Control V W cycle
      const int Nc_pre,     ///< n pre-smoothing cycles
      const int Nc_coarse,  ///< n coarse cycles
      const int Nc_post,    ///< n post-smoothing cycles
      const int clearing    ///< 1= clean the init matrix and precond
  );
  // --------------------------------------------------------------------
  virtual void MGCheck(int Level  /// \param[in] <>   MG level
                       ) const;   ///< Check Operators
  // -----------------------------------------------------------
  virtual void MGTimeStep(
      const double time,  /// \param[in] <>  time
      const int mode      /// \param[in] <>   rhs assembler flag
  );                      ///< MG time step solver (backward Euler)
  virtual void MGTimeStep_no_up(
      const double time,  /// \param[in] <>  time
      const int mode      /// \param[in] <>   rhs assembler flag
  );
  virtual void MGUpdateStep();
  virtual void MGUndo_disp();  ///< MG time step solver (backward Euler)

  // -----------------------------------------------------------
  // This function assembles the matrix   (defined in the MGSolver__)
  virtual void GenMatRhs(
      const double time,  /// \param[in] <>  time
      const int Level,    /// \param[in] <>   MG level
      const int mode      /// \param[in] <>   rhs assembler flag
      ) = 0;              /// Assemblying A matrix function

  // -----------------------------------------------------------
  virtual double MGFunctional(
      double parameter,  /// Use of the function: (0) compute functional OR (1) set _eta
      double& control    /// \param[in] <>  eta multiplier for optimal method
  );

  // ************************************************************************
  ///    EXTERNAL FIELDS (MGSolverDA.C)
  // ************************************************************************
  // virtual to be defined in MGSolverDA
  //    int Order                       = pol interpolation order
  //    int Field,                      = Field
  //    std::string SystemFieldName,    = Field name
  //    std::string SystemFieldName2,   = Coupled Field name
  //    int& n_index,                   = index
  //    int coupled                     = coupled or segregated
  //    int vector
  //   int dimension
  // ------------------------------------------------------------------------
  // Data =========================================================================================
  external_field _data_eq[3];  ///< external data structure
  // Functions ===================================================================================
  //   virtual void set_ext_fields(const std::vector<FIELDS> & pbName)=0;   ///< set external fields
  virtual void setUpExtFieldData();
  virtual void ActivateControl(
      int Order, int Field, std::string SystemFieldName, int& n_index, int vector, int dimension) = 0;
  virtual void ActivateEquation(int Order, int Field, std::string SystemFieldName, int& n_index) = 0;
  virtual void ActivateVectField(
      int Order, int Field, std::string SystemFieldName, int& n_index, int coupled) = 0;  // quad vector
  virtual void ActivateScalar(
      int Order, int Field, std::string SystemFieldName, int& n_index) = 0;  // quad scalar
  virtual void ActivateCoupled(
      int Order, int Field, std::string SystemFieldName, int& n_index, std::string SystemFieldName2) = 0;
  // ************************************************************************
  ///    MED
  // ****************** ******************************************************
  // Data =========================================================================================
  // Functions ===================================================================================
#ifdef HAVE_MED
  MEDCoupling::MEDCouplingFieldDouble* _ExtField;
  // --------------------------------------------------------------------
  virtual void print_u_med(  /// Print solution to a med file.
      std::string namefile,  // filename <-
      const int Level) = 0;  // MGLevel  <-
  // --------------------------------------------------------------------
  virtual void print_weight_med(  /// Print weight for control problems to a med file.
      std::string namefile,       // filename <-
      const int Level) = 0;       // MGLevel  <-
#endif
  // ****************** ******************************************************

#ifdef VANKA

  // ****************** ******************************************************
  // * (MGSolverBase_VANKA.C)
  // ******************  *****************************************************
  // Data =========================================================================================
  // Functions ===================================================================================
  void Vanka_solve(
      int Level, SparseMatrixM& matrix_in, NumericVectorM& solution_in, NumericVectorM& rhs_in,
      const double tol, const int m_its);
  // --------------------------------------------------------------------
  double Vanka_test(int Level  // Level
  );
  // --------------------------------------------------------------------
  virtual double MGStep_Vanka(  ///< MultiGrid Step
      int Level, double Eps1, int MaxIter, const int Gamma, const int Nc_pre, const int Nc_coarse,
      const int Nc_post);

#endif  // ****************************************************************
};
#endif
