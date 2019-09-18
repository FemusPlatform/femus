#ifndef __mgsolver_RANS_thermal_h__
#define __mgsolver_RANS_thermal_h__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_THERMAL_EQUATIONS
// ==================================

#include "User_RANS_thermal.h"
#include "RANS_thermal_parameters.h"
#include "MGFE_conf.h"
#include "MGSolverDA.h"

class MGEquationsSystem;

// =================================================
/// Class for mg energy equation solvers with name TBK_EQUATIONS. Multilevel and mulitporcessor class (see <a
/// href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolRANS_thermal : public MGSolDA {
  // ======================================================================================================
  // ============= MGSolTBK class data ======================================================================
  // ======================================================================================================
 public:
  RANS_thermal_param _RANS_t_parameter;

 protected:
  int _nTKdim;      //<dimension
  double _dt;       ///< =_mgutils.get_par("dt");
  int _FF_idx[30];  //< field equation flag

  // parameters--------------------------------------------------------------------------------------------
  // constant reference parameters (From paramater.in file _mgphys.get_par("Uref")) -----------------------
  const double _uref; /**< "Uref" */
  const double _lref; /**< "lref" */
  const double _Tref; /**< "Tref" */

  // constant fluid properties (From paramater.in file _mgphys.get_par("Uref")) ---------------------------
  const double _rhof;   /**< "rho0" */
  const double _muf;    /**< "mu0" */
  const double _cp0;    /**< "cp0" */
  const double _kappa0; /**< "kappa0" see paramater.in file */

  // nondimensional numbers -------------------------------------------------------------------------------
  double _alpha; /**< conducibility*/
  double _IPrdl; /**< Prandl*/
  double _IRe;   /**< Reynolds*/
  double _qheat; /** vol source*/
  double _qs;    /**< surface flux*/

  // turbulence
  double _alpha_eff;
  const double _InvSigma = 0.714285714286;
  double _alpha_turb;               /**< turb conducibility*/
  double _IPrdl_turb;               /**< turb Prandl number*/
  double _Tkappa_g[2], _kappa_g[2]; /**< reference kappa*/
  double _y_dist;                   ///< distance from the wall - interior points
  double _WallDist;                 ///< first mesh point
  double _sP;                       ///< turbulent tensor modulus
  double _sT;                       ///< temperature gradient modulus
  int _WallElement;
  double _NormMidCell[DIMENSION];
  // mesh -------------------------------------------------------------------------------------------------
  const int _offset;                     ///< = _mgmesh._NoNodes[_NoLevels-1]= mesh nodes
  double _xx_qnds[NDOF_FEM * DIMENSION]; /**< elem coords */
  double _xxb_qnds[NDOF_FEMB * DIMENSION];
  // -------------------- class field ---------------------------------------------------------------------

  // ------------------ integration -----------------------
  //  fields at gaussian points
  double _xxg[DIMENSION];         /**< gauss pts*/
  double _ub_g[3][30];            /**< external field  (0-1-2 degree)*/
  double _ub_dxg[2 * DIMENSION];  ///< external field derivative  (0-1-2 degree)
  double _InvJac2[DIMENSION * DIMENSION];

  // EXPLICIT AND IMPLICIT SOURCE AND DISSIPATION
  double _FemDiss[2], _KemDiss[2], _FemSource[2], _KemSource[2], _diss[2], _source[2];
  double _LowerLimit;

  // TURBULENCE MODEL
  double _KH_der[DIMENSION], _WH_der[DIMENSION];
  double _KH_2der[DIMENSION * DIMENSION], _WH_2der[DIMENSION * DIMENSION];

  // STABILIZATION PARAMETERS
  int _theta;
  int _NLITER, _RNSTAB, _Restart;
  bool _SUPG, _UPWIND, _SCPG, _ModifiedSupg, _UnderRelCorr;
  double _UnderRel;

  // BOUNDARY AND INITIAL CONDITION
  bool _BoundElem;
  bool _FlatProfile;

  std::vector<int> _BoundaryGroupsIDs;
  std::map<int, bound_condRANS_T> _BoundaryConditionMap;
  int _bc_vol[NDOF_FEM]; /**<  b.cond from function */
  int _bc_bd[NDOF_FEM];  /**< b.cond flags */
  int _bc_el[NDOF_FEM];  /**<  b.cond in matrix assemblying */

  int _SolveSteady;
  int _SolveTTBK;
  int _AxiSym;

  std::vector<int> _EquationNodes;

  int _ExplicitNearWallDer[2];  // 1: yes, 0:no
  double _explicit_diss[2], _explicit_source[2], _implicit_diss[2], _implicit_source[2], _mecc_term[2];
  double _Cross[2], _Log_Cross[2];
  double _Adv, _Lap, _Vel_g[DIMENSION], _AlphaTurbDxg[DIMENSION], _LapSupg, _LapMuTurb;

 public:
  // ==========================================================================
  // =========         Constructor - Destructor  ==============================
  // ==========================================================================
  /// This function constructs the 3d-2D MGSolTBK class
  MGSolRANS_thermal(                           ///< Constructor
      MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
      const int nvars_in[],                    ///< KLQ number of variables
      std::string eqname_in = "TK",            ///< equation name
      std::string varname_in = "kh"            ///< basic variable name
  );
  // ==========================================================================
  /// This function destructs the 3d-2D MGSolTBK class
  //     ~MGSolTBK() {}                              ///< Destructor
  ~MGSolRANS_thermal() {
    _BoundaryGroupsIDs.clear();
    _BoundaryConditionMap.clear();
  }
  //===========================================================================
  // =========== Read MGSolTBK functions  =======================================
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
  // ==========================================================================
  /// This function  computes a time step
  void MGTimeStep(
      const double time,  ///< Time-step manager function
      const int           /*iter */
  );
  void MGTimeStep_no_up(
      const double time,  ///< Time-step manager function
      const int           /*iter */
  );
  void MGUpdateStep();
  // ==========================================================================
  /// This function  computes the  functional defined (user write}
  void MGFunctional(
      const double /*time */, double /*starting_distance */
  ) {}
  // ==========================================================================
  /// This function sets the _bc_el with boundary condition flags.
  /// This function also assembles the surface integral to obtain the algebraic sytem.
  /// It is called by MGSolTBK::GenMatRhs in MGSolverT.C
  /// This function sets  the  functional defined (user write}
  void bc_set(
      int sur_toply[] /**< Local matrix  */, int el_ndof2 /**< el dofs*/, int elb_ndof2 /**< el bd dofs*/,
      int elb_ngauss
      /**< #pt gauss  */
  );
  // ==========================================================================
  /// This function assembles the volume integral to obtain the algebraic sytem.
  /// It is called by MGSolTBK::GenMatRhs in MGSolverT.C
  void vol_integral(
      const int el_ndof2,
      /**< el dofs*/
      const int el_ngauss,
      /**< #pt gauss  */
      const int mode, int el_conn[]);

  void turb_utils();
  void GetSolution(int iel, int Level, int el_conn[], int offset, int el_ndof[], int ndof_lev);
  void GaussInterpolationOoldSolution();

  void CalcSUPG(double& h_eff, double& f_upwind, double mod2_vel, double vel_g[]);

  double CalcPhiSupg(
      int i, double vel_g[], double ParVel[], double tauc, double f_upwind, double implicit_diss,
      int el_ndof2);

  void StandardTimeStep(const double time);

  virtual void inline CalcSourceAndDiss() {
    _explicit_source[0] = _source[0];
    _explicit_diss[0] = _diss[0];
    _implicit_diss[0] = 0.;
    // second equation source and diss -> eh or wh or log_wh
    _explicit_source[1] = _source[1] + _mecc_term[1];
    _explicit_diss[1] = _diss[1] + _mecc_term[0];
    _implicit_diss[1] = 0.;
  };

  virtual void CalcAdvectiveAndDiffusiveTerms(
      int TestFunctionID, int InterpolationNode, int NumOfNodes, double f_upwind);
  virtual void VelocityForSUPG(double& mod2_vel, double vel_g[], double VEL[]);
  void VelOnGaussAndTensorModulus(double& mod2_vel, int NumOfNodes);
  void TempGradientModulus(double& mod2_vel, int NumOfNodes);
};

#endif  // define TTBK_EQUATIONS
#endif  //__mgsolverTK_h__
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
