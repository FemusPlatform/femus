#ifndef __mgsolbase__
#define __mgsolbase__

#include "Printinfo_conf.h"   //
#include "EquationsMap.h"
#include "Solverlib_conf.h"
#include <vector>
#include <string>
#include <stdio.h>
#include "MEDCouplingFieldDouble.hxx"
class SparseMatrixM;
class SparseMMatrixM;
class NumericVectorM;
class LinearSolverM;

class MGUtils        ;
class MGSystem       ;
class MGEquationsSystem ;
class MGFEMap        ;
class MGMesh         ;
#ifdef    TWO_PHASE_LIB
class MGSolCC;
#endif

///This class is the basic equation class
class MGSolBase
#ifdef LM_REFCOUNT
  : public ReferenceCountedObject<MGSolBase>
#endif
  {


  protected:
    //  ----------------------------------
    ///@{ \name DATA POINTER
//   MGUtils          & _mgutils;  ///<  utility class pointer
// //   MGSystem         & _mgphys;   ///<  parameter class pointer
//   MGFEMap          & _mgfemap;  ///<  FEM class pointer
    MGMesh      &      _mgmesh;   ///<  mesh pointer
    MGEquationsSystem  &  _mgeqnmap; ///<  equation map  pointerz
    ///<
    ///@}
    //  ----------------------------------
    ///@{ \name DATA PARALLEL
    int _iproc;                  ///< processor
    LinearSolverM ** _solver;    ///< linear system solver type (each level)
    ///<
    ///@}
    //  ---------------------
    ///@{ \name  MULTIGRID MATRICES AND RHS'S
    std::vector<SparseMatrixM *> A;     ///< Matrix A
    std::vector<NumericVectorM *> b;    ///< rhs b
    std::vector<NumericVectorM *> x;    ///< solution x
    std::vector<NumericVectorM *> res;  ///< residual
    ///@}
    // Multigrid operators ------------------------------
    ///@{  \name MULTIGRID OPERATORS
    std::vector<SparseMMatrixM *> Rst;  ///< Restrictor
    std::vector<SparseMMatrixM *> Prl;  ///< Prolongation
    ///<
    ///@}
  public:
#ifdef    TWO_PHASE_LIB
    MGSolCC * _msolcc;
#endif
    double  _control;
    MGUtils      &     _mgutils;  ///<  utility class pointer
//   MGSystem         & _mgphys;   ///<  parameter class pointer
    MGFEMap       &      _mgfemap;  ///<  FEM class pointer
    //  --------------------------
    ///@{  \name OLD VECTOR

    std::vector<NumericVectorM *> x_old;    ///< old solution x
    std::vector<NumericVectorM *> x_oold;   ///< oold solution x
    std::vector<NumericVectorM *> x_nonl;   ///< non linear solution x
    std::vector<NumericVectorM *> disp;     ///< displacement solution x
    std::vector<NumericVectorM *> disp_old; ///< displacement old solution x
    std::vector<NumericVectorM *> disp_oold; ///< displacement old solution x
    std::vector<NumericVectorM *> x_ooold;  ///< vector for multiple uses
    std::vector<double> _weight_ctrl;       ///<controlled region for optimal control problems
#ifdef HAVE_MED
    MEDCoupling::MEDCouplingFieldDouble * _ExtField;
#endif
    ///<
    ///@}
    // ---------------------------------------
    ///@{  \name  ATTRIBUTES
    const int _NoLevels;        ///< level number
    int * _Dim;                 ///< dimension number of dofs per level
    ///<
    ///@}
    // ---------------------------------------
    ///@{  \name LABELING
    const std::string _eqname;  ///< equation name
    ///<
    ///@}
    // ---------------------------------------
    ///@{  \name  VARIABLES
    const int _n_vars;         ///< number of variables
    int _nvars[3];        ///< number of variables quadratic[2], linear[1] and piecewise [0]
    std::string * _var_names;  ///< variable names
    double * _refvalue;        ///reference values
    ///<
    ///@}
    // ---------------------------------------
///@{  \name BC AND DOF MAP
    int * bc[2];              ///< boundary conditions map (top level)
    int ** _node_dof;         ///< dof map
///<
    ///@}

    // ======================================
    ///@{ \name CONSTRUCTOR-DESTRUCTOR

    MGSolBase (
      MGEquationsSystem & mg_equations_map, ///< \param[in] <> MG equation map
      const int nvars_in[],              ///< \param[in] <> number of variables
      std::string eq_name_in = "Base"    ///< \param[in] <> equation name
    );                          ///< Level constructor
    MGSolBase (
      MGUtils & mgutils_in,
      MGFEMap & mgfemap_in,
      MGMesh & mgmesh_in,
      const int nvars_in[],            // # variables
      std::string eqname_in     // equation name
    );

//-------------------------------------------------------------------------
    /// Destructor (level structure)
    ~MGSolBase();
//-------------------------------------------------------------------------
    /// Substructure destructor
    void clear();
//-------------------------------------------------------------------------
    /// Distributing dof function
    virtual void init_dof (
      const int Level                ///< \param[in]
    ) = 0;
//-------------------------------------------------------------------------
    ///  Initializing memory allocation
    virtual void     init (
      const int Level               ///< \param[in] <> Level
    ) = 0;
    ///
    ///@}

    //----------------------------------------------------------------------
    ///@{ \name RETURN FUNCTIONS
    ///< \param[in]   <ivar0>   initial variable
    ///< \param[in]   <nvars>   number of variables to get
    ///< \param[in]   <el_nds>  number  of element nodes for this variable
    ///< \param[in]  <el_conn> connectivity
    ///< \param[in]  <offset>  offset for connectivity
    ///< \param[in]  <kvar0>   offset  variable for  uold
    ///< \param[out]  <uold>   solution

#ifdef    TWO_PHASE_LIB
    void set_mgcc ( MGSolCC & cc );
#endif
    // ----------------------------------------------------------------
    void set_ctrl_dom (
      const double x_min,
      const double x_max,
      const double y_min,
      const double y_max,
      const double z_min,
      const double z_max
    );
    // ----------------------------------------------------------------
    void  get_el_sol (
      const int ivar0,              ///< \param[in]   <ivar0>   initial variable
      const int nvars,              ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,             ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],          ///< \param[in]  <el_conn> connectivity
      const int offset,             ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,              ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                ///< \param[out]  <uold>   solution
      ///<
    )  const  ;               ///< Return element solution
// ----------------------------------------------------------------
    void  get_el_sol_piece (
      const int ivar0,              ///< \param[in]   <ivar0>   initial variable
      const int nvars,              ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,             ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int iel,                ///< \param[in]  <iel> element id (ofset)
      const int offset,             ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,              ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                ///< \param[out]  <uold>   solution
      ///<
    )  const  ;               ///< Return element solution

// --------------------------------------
/// Return old element solution
    void  get_el_oldsol (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
    ///@}
// --------------------------------------
/// Return no linear  solution on element
    void  get_el_nonl_sol (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
    ///@}

/// Return ooold element solution
    void  get_el_oooldsol (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
///@}
/// Return disp element solution
    void  get_el_disp (
      const int Level,                ///< \param[in]   <Level>   Level
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
/// Return disp element solution
    void  get_el_disp (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
///@}
/// Return disp element solution
    void  get_el_new_disp (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
///@}
/// Return disp element solution
    void  get_el_oold_disp (
      const int ivar0,                ///< \param[in]   <ivar0>   initial variable
      const int nvars,                ///< \param[in]   <nvars>   number of variables to get  <-
      const int el_nds,               ///< \param[in]   <el_nds>  number  of element nodes for this variable
      const int el_conn[],            ///< \param[in]  <el_conn> connectivity
      const int offset,               ///< \param[in]  <offset>  offset for connectivity
      const int kvar0,                ///< \param[in]  <kvar0>   offset  variable for  uold
      double  uold[]                  ///< \param[out]  <uold>   solution
    )  const  ;
///@}
    // =============================
    virtual void set_vector ( const int &, const int & ) = 0;
    virtual void set_xooold2x() = 0;
// virtual void  f_aux (double a1[]) const;
    // ==================================================================
    // ===== INITIAL AND BOUNDARY CONDITIONS (all virtual)  =============
    // ==================================================================
///@{ \name INITIAL AND BOUNDARY CONDITIONS
    virtual void    GenIc() = 0;    ///< Reading Initial conditions (IC)
//   virtual void    GenIc0(
//     const std::string ibc_file,        // initial condition file
//     const int node_dof_top[],          // node map (top mesh)
//     NumericVectorM &sol_top,           // solution (top mesh)
//     NumericVectorM &old_sol_top        // old solution (top mesh)
//   )=0;                            ///<  Reading IC  function


    virtual void  ic_read ( /// initial conditions  from function
      int    k,           // bc from gambit    <-
      int    m,           // material from gambit    <-
      double xp[],        // point coordinates <-
      int iel,           // element  <-
      double valueic[]        //  point values     ->
    ) = 0;                          ///< Reading IC  function (element)
    virtual void  GenBc() = 0;      ///< Reading boundary conditions (BC)
    virtual void  bc_read (
      int k,
      int m, ///< global node id
      double x[],                       ///< point vector
      int bc_vol[],                     ///< Volume flag
      int bc_sur[]                      ///< Surface flag
    ) = 0;                          ///< Reading sur BC function (element)
    /// Reading vol BC function (element)
    virtual void  bc_intern_read (
      int k,
      int m,
      double x[],                       // point vector
      int bc_vol[],                     // Volume flag
      int u[]                           // value vector
    ) = 0; // --------------------------------------------------------------
    ///@}

    double CalcFUpwind ( double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes );

//   double CalcPhiSupgContribution(int i,
//         double vel_g[],
//         double ParVel[],
//         double tauc,
//                     double f_upwind,
//         double implicit_diss,
//         int el_ndof2);

    // ===================================================================
    // ================== Print/Read function (All virtual)====================
    // ===================================================================
    // Print
    ///@{ \name PRINT/READ
    ///< \param[in] <fname>  file name
    ///< \param[in] <Level>  MG level
    ///< \param[in] <of_out>  out stream file
    ///< \param[in] <n_nodes>  number of nodes
    ///< \param[in] <n_elems>  number of elements


    // --------------------------------------------------------------------
    virtual void  print_u (
      std::string fname,          ///< \param[in] <>  file name
      const int Level             ///< \param[in] <> MG level
    ) = 0;                    ///< print solution
#ifdef HAVE_MED
    // --------------------------------------------------------------------
    virtual void print_u_med ( /// Print solution to a med file.
      std::string namefile,    // filename <-
      const int Level ) = 0;     // MGLevel  <-
#endif
#ifdef HAVE_MED
    virtual void print_weight_med ( /// Print weight for control problems to a med file.
      std::string namefile,    // filename <-
      const int Level ) = 0;     // MGLevel  <-
#endif
    // --------------------------------------------------------------------
    virtual void  print_bc (
      std::string fname,          ///< \param[in] <>  file name
      const int Level             ///< \param[in] <> MG level
    ) = 0;                    ///< print boundary conditions
    // --------------------------------------------------------------------
    void print_ext_data (
      double /*vect_data*/[]
    ) {}/*printf("\n \n Wrong use of function print_ext_data in SolverBase for coupled mesh \n \n") =0*/;
    // --------------------------------------------------------------------
    virtual void print_xml_attrib (
      std::ofstream & of_out,     ///< \param[in] <>  out stream file
      int n_nodes,                ///< \param[in] <>  number of nodes
      int n_elems,                ///< \param[in] <>  number of elements
      std::string fname           ///< \param[in] <>  file name
    ) const;                 ///< print xml file
    // --------------------------------------------------------------------
    virtual void read_u (
      std::string fname,         ///< \param[in] <>  file name
      int Level                  ///< \param[in] <>  restart MG level
    ) = 0;                   ///< Read solution
///@}

    // ==================================================
///@{ \name  MULTILEVEL OPERATOR
    //======================================================
    virtual void MGDofBcOp();    ///< Reading Operators
    virtual void ReadMatrix (    ///< Reading matrix A
      const int Level,           //  MG level <-
      const std::string & name,
      SparseMatrixM   &  Mat,
      const int     *    nvars_all
    ) = 0;
    virtual void ReadProl (       ///< Reading Prolongation
      const int Level,            //  MG level <-
      const std::string & name,
      SparseMMatrixM  &  Mat,
      const int          nvars_in [],
      int                node_dof_c[],
      int                node_dof_f[]
    ) = 0;
    virtual void ReadRest (      ///< Reading Restriction
      const int Level,           //  MG level <-
      const std::string & name,
      SparseMMatrixM   &  Mat,
      const int           nvars[],
      int                 node_dof_c[],
      int                 node_dof_f[],
      int                _node_dof_top[]
    ) = 0;
///@}
// ==================================================
    ///@{ \name MULTILEVEL SOLUTION/ASSEMBLYING
// ==================================================
    virtual void GenMatRhs (
      const double time,        /// \param[in] <>  time
      const int    Level,       /// \param[in] <>   MG level
      const int    mode         /// \param[in] <>   rhs assembler flag
    ) = 0;              /// Assemblying A matrix function
    // -----------------------------------------------------------
    virtual  void  MGTimeStep (
      const double time,        /// \param[in] <>  time
      const int    mode         /// \param[in] <>   rhs assembler flag
    );                ///< MG time step solver (backward Euler)

    virtual  void  MGTimeStep_no_up (
      const double time,        /// \param[in] <>  time
      const int    mode         /// \param[in] <>   rhs assembler flag
    );
    virtual  void MGUpdateStep ();

    // -----------------------------------------------------------
    virtual  void  MGUndo_disp();                ///< MG time step solver (backward Euler)
    // -----------------------------------------------------------
    virtual  void  set_dt ( double dt );             ///< MG time step solver (backward Euler)
    // -----------------------------------------------------------
    virtual  double  GetValue ( int flag );             ///< This function returns a value from the solvers
    // -----------------------------------------------------------
    virtual  void  SetValue ( double value );             ///< This function sets a value in the solvers
    // -----------------------------------------------------------
    virtual  void  SetValueVector ( std::vector<double> value );             ///< This function sets a value in the solvers
    // -----------------------------------------------------------
    virtual  double  MGFunctional (
      double parameter,     /// Use of the function: (0) compute functional OR (1) set _eta
      double & control /// \param[in] <>  eta multiplier for optimal method
    );
// --------------------------------------------------------------------
    /// This function performes the V,W and Fast MultiGrid scheme
    virtual  void  MGSolve (
      double    Eps,          ///< tolerance
      int       MaxIter,      ///< n iterations
      // -------------------------------------------------------------
      const int clearing = 0, ///< 1= clean the init matrix and precond
      const int Gamma = 1,    ///< Control V W cycle
      const int Nc_pre = 8,   ///< n pre-smoothing cycles
      const int Nc_coarse = 40, ///< n coarse cycles
      const int Nc_post = 8   ///< n post-smoothing cycles
    );
    // ------------------------------------------------------------------
    /// This function solves a MultiGrid step
    virtual double MGStep (
      int       Level,         ///<  MG level <-
      double    Eps1,          ///< tolerance
      int       MaxIter,       ///< n iteratio
      const int Gamma,         ///< Control V W cycle
      const int Nc_pre,        ///< n pre-smoothing cycles
      const int Nc_coarse,     ///< n coarse cycles
      const int Nc_post,       ///< n post-smoothing cycles
      const int clearing       ///< 1= clean the init matrix and precond
    );
    // --------------------------------------------------------------------
    virtual  void  MGCheck (
      int Level               /// \param[in] <>   MG level
    ) const;             ///< Check Operators
///@}

    void   localize_xooold();
// ==================================================
    ///@{ \name external fields
    // ==================================================
//   virtual void set_ext_fields(const std::vector<FIELDS> & pbName)=0;   ///< set external fields
    virtual void setUpExtFieldData() = 0;
    virtual void ActivateControl   ( int Order, int Field, std::string SystemFieldName, int & n_index, int vector, int dimension ) = 0;
    virtual void ActivateEquation  ( int Order, int Field, std::string SystemFieldName, int & n_index ) = 0;
    virtual void ActivateVectField ( int Order, int Field, std::string SystemFieldName, int & n_index, int coupled ) = 0; // quad vector
    virtual void ActivateScalar    ( int Order, int Field, std::string SystemFieldName, int & n_index ) = 0; // quad scalar
    virtual void ActivateCoupled   ( int Order, int Field, std::string SystemFieldName, std::string SystemFieldName2, int & n_index ) = 0;



///@}



// ========================================================================
// ****************** HERE STARTS VANKA SECTION ***************************
// ========================================================================
void Vanka_solve(
    int Level,
  SparseMatrixM&  matrix_in,
  NumericVectorM& solution_in,
  NumericVectorM& rhs_in,
  const double tol,
  const  int m_its
);
//
  double Vanka_test(
    int Level           // Level
  );
//
  virtual double MGStep_Vanka(    ///< MultiGrid Step
    int       Level,       //  MG level <-
    double    Eps1,
    int       MaxIter,
    const int Gamma,
    const int Nc_pre,
    const int Nc_coarse,
    const int Nc_post
  );

// #endif // ****************************************************************
  };
#endif
