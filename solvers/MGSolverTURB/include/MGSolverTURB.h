#ifndef __mgsolverTURB_h__
#define __mgsolverTURB_h__

#include "Equations_conf.h"

// ===================================
#ifdef _TURBULENCE_
// ==================================

// config files ------------------
#include "MGFE_conf.h"
#include "MGSolverDA.h"
#include "UserTURB.h"

// Forward declarations ----------
class MGEquationsSystem;

// =================================================
/// Class for mg energy equation solvers with name T_EQUATIONS. Multilevel and mulitporcessor class (see <a
/// href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolTURB : public MGSolDA {
  // ======================================================================================================
  // ============= MGSolTURB class data ======================================================================
  // ======================================================================================================
 public:
  TURB_param _C_parameter;

 private:
  int _nTdim;       //<dimension
  int _FF_idx[40];  //< field equation flag
  // parameters--------------------------------------------------------------------------------------------
  // constant reference parameters (From paramater.in file _mgphys.get_par("Uref")) -----------------------
  const double _uref; /**< "Uref" */
  const double _lref; /**< "lref" */
  const double _Tref; /**< "Tref" */

  // constant fluid properties (From paramater.in file _mgphys.get_par("Uref")) ---------------------------
  const double _rhof; /**< "rho0" */
  const double _muf;
  // nondimensional numbers -------------------------------------------------------------------------------
  double _IRe; /**< Reynolds*/

  // turbulence
  double _y_dist;  ///< distance from the wall

  // mesh -------------------------------------------------------------------------------------------------
  const int _offset;                     ///< = _mgmesh._NoNodes[_NoLevels-1]= mesh nodes
  double _xx_qnds[NDOF_FEM * DIMENSION]; /**< elem coords */
  double _xxb_qnds[NDOF_FEMB * DIMENSION];

  // -------------------- class field ---------------------------------------------------------------------
  // element boundary conditions
  int _bc_vol[NDOF_FEM]; /**<  b.cond from function */
  int _bc_bd[NDOF_FEM];  /**< b.cond flags */
  int _bc_el[NDOF_FEM];  /**<  b.cond in matrix assemblying */
  // ------------------ integration -----------------------
  //  fields at gaussian points
  double _ub_g[3][40];    /**< external field  (0-1-2 degree)*/
  double _xxg[DIMENSION]; /**< gauss pts*/
  double _InvJac2[DIMENSION * DIMENSION];
  double _ub_dxg[2 * DIMENSION];  ///< external field derivative  (0-1-2 degree)

  double _KW[2 * NDOF_FEM], _KHWH[2 * NDOF_FEM], _DIST[NDOF_FEM];

  std::vector<int> _WallGroupsIDs;
  int _FirstCall;
  int _FirstAssembly;
  int _AssembleOnce;
  int _Solve;
  double _DiffusionCoefficient;
  double _res;

 public:
  // ==========================================================================
  // =========         Constructor - Destructor  ==============================
  // ==========================================================================
  /// This function constructs the 3d-2D MGSolTURB class
  MGSolTURB(                                   ///< Constructor
      MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
      const int nvars_in[],                    ///< KLQ number of variables
      std::string eqname_in = "C",             ///< equation name
      std::string varname_in = "c"             ///< basic variable name
  );
  // ==========================================================================
  /// This function destructs the 3d-2D MGSolTURB class
  ~MGSolTURB() {}  ///< Destructor

  //===========================================================================
  // =========== Read MGSolTURB functions  =======================================
  // ==========================================================================

  // This function reads the boundary conditions
  void bc_read(  ///< Read bc
      int bc_gam /**< group idx  */, int bc_mat /**< mat idx      */, double xp[] /**< pts coords */,
      int bc_Neu[] /**< bc flags   */, int bc_value[]
      /**< aux bc flags*/
  );
  // ==========================================================================
  // This function reads the initial solution
  void ic_read(  ///< Read bc
      int bc_gam /**< group idx    */, int bc_mat /**< mat id     */, double xp[] /**< pts coords   */,
      int iel /**< element id */, double u_value[]
      /**< value vector */
  );  ///< Read ic
  /// d)  Assemblying MGSolNS Operators

  // ==============================================================================
  // This function assembles the Volume integral
  void GenMatRhs(
      const double time, /**< time*/
      const int Level /**< discrtization Level*/, const int mode
      /**< y/n assembnle rhs */
  );

  void GenRhs(
      const double time, /**< time*/
      const int Level /**< discrtization Level*/, const int mode
      /**< y/n assembnle rhs */
  );

  // ==========================================================================
  /// This function  computes a time step
  void MGTimeStep(
      const double time,  ///< Time-step manager function
      const int           /*iter */
  );
  // ==========================================================================
  /// This function  computes the  functional defined (user write}
  void MGFunctional(
      const double /*time */, double /*starting_distance */
  ) {}
  // ==========================================================================
  /// This function sets the _bc_el with boundary condition flags.
  /// This function also assembles the surface integral to obtain the algebraic sytem.
  /// It is called by MGSolTURB::GenMatRhs in MGSolverT.C
  /// This function sets  the  functional defined (user write}
  void bc_set(
      DenseMatrixM& KeM /**< Local matrix  */, DenseVectorM& FeM /**< Local rhs */,
      int sur_toply[] /**< Local matrix  */, int el_ndof2 /**< el dofs*/, int elb_ndof2 /**< el bd dofs*/,
      int elb_ngauss /**< #pt gauss  */, int sign_normal
      /**< old velocity field  */
  );
  // ==========================================================================
  /// This function assembles the volume integral to obtain the algebraic sytem.
  /// It is called by MGSolTURB::GenMatRhs in MGSolverT.C
  void vol_integral(
      DenseMatrixM& KeM /**< Local matrix  */, DenseVectorM& FeM /**< Local rhs */,
      const int el_ndof2 /**< el dofs*/, const int el_ngauss /**< #pt gauss  */,
      double xx_qnds[] /**< el coords  */, const int unsteady /**< un/steady flag  */,
      const int mode /**<  rhs int flag */, int el_conn[]);

  void compute_y_plus(
      const int el_ndof2 /**< el dofs*/, const int el_ngauss /**< #pt gauss  */,
      double xx_qnds[] /**< el coords  */, const int unsteady /**< un/steady flag  */,
      const int mode /**<  rhs int flag */, int el_conn[], int iel);

  void rhs_integral(
      DenseVectorM& FeM /**< Local rhs */, const int el_ndof2 /**< el dofs*/,
      const int el_ngauss /**< #pt gauss  */, double xx_qnds[] /**< el coords  */,
      const int unsteady /**< un/steady flag  */, const int mode /**<  rhs int flag */, int el_conn[]);
};

#endif  // define T_EQUATIONS
#endif  //__mgsolverT_h__

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
