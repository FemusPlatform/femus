#ifndef __mgsolverda_h__
#define __mgsolverda_h__

// conf includes --------------------------------------------------------------
#include "MGFE_conf.h"

// algebra class --------------------------------------------------------------
#include "dense_matrixM.h"
#include "dense_vectorM.h"
class NumericVectorM;

// local class ----------------------------------------------------------------
#include "MGFEMap.h"
#include "MGSolverBase.h"  //for the inherited class
class MGUtils;
class MGSystem;
class MGMesh;
class MGEquationsMap;

#include "hdf5.h"

// ================================================================================================
/// Class for MG Diffusion-Advection equation generic solvers
// ================================================================================================
class MGSolDA : public MGSolBase {
 public:
  // ************************************************************************************************
  //                   DATA
  // ************************************************************************************************
  // 1. Common data MGSolDA -------------------------------------------------------------
  // _dt = time step
  // _dir= eq index for vector equation, neq=3 _dir=0 (x) _dir=1 (y)  _dir=2 (z)
  // _NumRestartSol= Number of solutions needed for restart
  //                 This number is used for printing and reading solutions for system restart.
  //                 Default value is equal to 1. With value _NumRestartSol=2 the solution
  //                 stored in _x_oold vector is printed with _old suffix, while with _NumRestartSol=3
  //             also the _x_ooold stored solution is printed, with _oold suffix (visible inside .h5 file)

  int _NumRestartSol;  ///< Number of solutions needed for restart

  // 2. FEM  ----------------------------------------------------------------------------------------
  // int _el_dof[3]= number of dof in each element  (piecewise, linear, quadratic)
  // MGFE* _fe[3]  = fem  (piecewise linear,piecewise quadratic)
  int _el_dof[3];  ///< number of dof in each element  (piecewise, linear, quadratic)
  MGFE* _fe[3];    ///< fem  (piecewise linear,piecewise quadratic)
  // shape and derivative function at a  gaussian point
  double _phi_g[3][NDOF_FEM];                            ///< shape field (0-1-2 degree)
  double _dphi_g[3][NDOF_FEM * DIMENSION];               ///< shape derivative  (0-1-2 degree)
  double _ddphi_g[3][NDOF_FEM * DIMENSION * DIMENSION];  ///< shape second derivative  (0-1-2 degree)

  // 3. Local MATRIX and VECTOR  ---------------------------------------------------------------
  // DenseMatrixM _KeM  = local  matrix
  //  DenseVectorM _FeM = local  rhs
  DenseMatrixM _KeM;  ///< local  matrix
  DenseVectorM _FeM;  ///< local  rhs

  // ========================================================================
  // 1. Constructor destructor memory allocation
  // ========================================================================
  // const int nvars_in[],    = number of piecewise[0], linear[1], quadratic[2] variables
  //  std::string eq_name_in  = Equation name  ("DA")
  //  std::string varname_in  = Variable suffix name (u)
  //  const int Level,        = MG Level
  //  const int vb_0 = 0,     = mesh element fem type (vb_0=0 volume  (vb_0=1 surface)
  //  const int n_vb = 1      = mesh element groups (n_vb=1 only one fem type)
  // ========================================================================
  /// Level constructor  I
  MGSolDA(
      MGEquationsSystem& mg_equations_map,
      const int nvars_in[],           ///< number of piecewise[0], linear[1], quadratic[2] variables
      std::string eq_name_in = "DA",  ///< Equation name  ("DA")
      std::string varname_in = "u"    ///< Variable suffix name (u)
  );
  // Destructor------------------------------------------------------------------
  ~MGSolDA();    ///< Destructor (level structure)
  void clean();  ///< Clean all substructures
  // Initial Setting -----------------------------------------------------------
  void init_dof(        ///< Setting dof
      const int Level,  // MG Level
      const int vb_0 = 0, const int n_vb = 1);
  void init(           ///< Setting Memmory alloc
      const int Level  // MG Level
  );
  // ==========================================================================
  // BOUNDARY/INITIAL CONDITIONS
  // ==========================================================================
  //   _bc[2];  ///< boundary conditions map (top level
  // -------------------------------------------------------------------------
  virtual void GenBc();  ///< Setting boundary conditions
  // -------------------------------------------------------------------------
  virtual void GenBc_loop(
      const int vb, const int ndof_femv, const int n_ub_dofs, const int n_pb_dofs, const int n_kb_dofs,
      int bc_id[], int mat_id[]);
  // -------------------------------------------------------------------------
  virtual void GenIc();  ///< Element Initial conditions

  // ========================================================================
  /// MULTILEVEL OPERATORS
  // ========================================================================
  // Reading/Writing MG operators ------------------------------------------
  virtual void ReadProl(        ///< Read Prolongation Op
      const int Level,          // MG Level <-
      const std::string& name,  // file name for P <-
      SparseMMatrixM& Mat,      //  Matrix for P <-
      const int nvars_in[], int node_dof_c[], int node_dof_f[]);
  // -------------------------------------------------------------------------
  virtual void ReadRest(        ///< Restriction Op.
      const int Level,          // MG Level
      const std::string& name,  // file name (reading from)
      SparseMMatrixM& Rest,     // Restriction Matrix
      const int nvars_in[],     // # cost,linear quad variables
      int node_dof_f[],         // dof map fine mesh
      int node_dof_c[],         // dof map coarse
      int _node_dof_top[]       // dof map top level
  );                            // ----
                                // -------------------------------------------------------------------------
  virtual void ReadMatrix(      ///< Reading Matrix
      const int Level,          // MG Level <-
      const std::string& name,  // file name for M <-
      SparseMatrixM& Mat, const int* nvars_in);

  // ========================================================================
  //    MultiGrid Solver (defined in the MGSolver__)
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

  ///
  virtual double CalcFUpwind(
      double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes);

  // ========================================================================
  /// RETURN FUNCTIONS
  // ========================================================================
  // get_sol  retrun sol and old sol on top level --------------------------
  //   const int ivar0,     = initial variable  <-
  //    const int nvars,    = number of variables to get  <-
  //    const int el_nds,   =   # of element nodes for this variable  <-
  //    const int el_conn[] =  element  connectivity <-
  //    const int offset,   =   offset for connectivity <-
  //    const int kvar0,    =   offset  variable for  uold <-
  //    double uold[]       =   element node values ->
  // -----------------------------------------------------------------------
  void get_el_sol(
      const int ivar0,      // initial variable  <-
      const int nvars,      // # of variables to get  <-
      const int el_nds,     // # of element nodes for this variable  <-
      const int el_conn[],  // connectivity <-
      const int offset,     // offset for connectivity <-
      const int kvar0,      // offset  variable for  uold <-
      double uold[]         // element node values ->
      ) const;
  // -----------------------------------------------------------------------
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
  // -----------------------------------------------------------------------
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
      ) const;

  // ========================================================================
  /// SET FUNCTIONS
  // ========================================================================
  //    virtual void set_cp_vector(const int&, const int&);
  // -----------------------------------------------------------------------
  virtual void set_xooold2x();
  // -----------------------------------------------------------------------
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
      ) const;
  // -----------------------------------------------------------------------

  void set_dt(double dt) { _dt = dt; };  // computation val3=musker

  // ========================================================================
  /// COMPUTE FUNCTIONS
  // ========================================================================
  // -------------------------------------------------------------------------
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

  // // -------------------------------------------------------------------------
  virtual double MGFunctional(double, double&);

  // ========================================================================
  ///  EXTERNAL FIELDS
  // ========================================================================
  //    int Order                       = pol interpolation order
  //    int Field,                      = Field
  //    std::string SystemFieldName,    = Field name
  //    std::string SystemFieldName2,   = Coupled Field name
  //    int& n_index,                   = equation index
  //    int coupled                     = coupled (1) or segregated (0)  solver
  //    int vector                      = vector
  //    int neq                         = number of equations
  //  =========================================================================================
  // ------------------------------------------------------------------------
  //  virtual
  // auxiliary function fro user purpose
  virtual double eval_var1(double[]) { return 0.; };
  virtual double eval_var2(double[]) { return 0.; };
  virtual double eval_var3(double[]) { return 0.; };

  // ------------------------------------------------------------------------
  void ActivateVectField(int Order, int Field, std::string SystemFieldName, int& n_index, int coupled);
  void ActivateScalar(int Order, int Field, std::string SystemFieldName, int& n_index);
  void ActivateControl(int Order, int Field, std::string SystemFieldName, int& n_index, int vector, int neq);
  void ActivateCoupled(
      int Order, int Field, std::string SystemFieldName, int& n_index, std::string SystemFieldName2);
  void ActivateEquation(int Order, int Field, std::string SystemFieldName, int& n_index);

  // ========================================================================
  /// SOLUTION WRITE AND READ
  // ========================================================================
  // Print
  virtual void print_bc(     ///< Print boundary conditions to a xdmf file.
      std::string namefile,  // filename <-
      const int Level        // MGLevel  <-
  );

  // -------------------------------------------------------------------------
  // Read
  virtual void read_u(   /// Read solution from a xdmf file
      std::string name,  // filename <-
      int Level          // restart level <-
  );
  // -------------------------------------------------------------------------
  /// This function  defines the boundary conditions for the system:
  virtual void bc_intern_read(
      int /*face_id_node*/,  ///<  face identity           (in)
      int /*mat_flag*/,      ///<  volume identity         (in)
      double /*xp*/[],       ///< xp[] node coordinates    (in)
      int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
      int bc_flag[]          ///< boundary condition flag  (out)
  );
  // -------------------------------------------------------------------------
  /// This function reads Boundary conditions  from function
  virtual void bc_read(
      int /*face_id_node*/,  ///<  face identity          (in)
      int /*mat_flag*/,      ///<  volume identity         (in)
      double /*xp*/[],       ///< xp[] node coordinates    (in)
      int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
      int bc_flag[]          ///< boundary condition flag  (out)
  );
  // -------------------------------------------------------------------------
  virtual void ic_read(  /// initial conditions  from function
      int k,             // bc from gambit    <-
      int m,             // material from gambit    <-
      double xp[],       // point coordinates <-
      int iel,           // element           <-
      double value[]     //  point values     ->
  );

  // -------------------------------------------------------------------------
  void print_u(              /// Print solution to a xdmf file.
      std::string namefile,  // filename <-
      const int Level        // MGLevel  <-
  );
  // -------------------------------------------------------------------------
  virtual void print_u_xdmf(
      std::ofstream& out,  // file stream to print ->
      int nodes, int nelems, std::string file_name) const;

  // =============================================
  /// Print xml attrib
  void print_u_xdmf_pts(
      std::ofstream& out,    ///< xdmf file
      int nodes,             ///< number of nodes
      std::string var_name,  ///< quad linear var name
      std::string file_name  ///<  xdmf file name

      ) const;
  // =============================================
  /// Print xml attrib
  void print_u_xdmf_elem(
      std::ofstream& out,  //  file xdmf
      int nelems, std::string var_name, std::string file_name

      ) const;

  /// Print solution in hdf5
  virtual void print_u_hdf5(std::string namefile, const int Level);
  /// This function prints the solution for quad in hdf5
  virtual void print_u_hdf5_quad(
      const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[]);
  /// This function prints the solution for  linear fem in  hdf5
  virtual void print_u_hdf5_lin(
      const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[]);
  /// This function prints the solution: for quad and linear fem
  virtual void print_u_hdf5_const(
      const int Level, const int istep, const int offset, hid_t file_id, hsize_t dimsf[]);

  // ************************************************************************
  //                                  MED
  // ************************************************************************
#ifdef HAVE_MED
  // ------------------------------------------------------------------------
  virtual void print_u_med(       /// Print solution to a med file.
      std::string namefile,       // filename <-
      const int Level);           // MGLevel  <-
                                  // ------------------------------------------------------------------------
  virtual void print_weight_med(  /// Print weight for control problems to a med file.
      std::string namefile,       // filename <-
      const int Level);           // MGLevel  <-
#endif
};

#endif
