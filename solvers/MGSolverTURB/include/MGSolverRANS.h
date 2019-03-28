#ifndef __mgsolverRANS_h__
#define __mgsolverRANS_h__


#include "Equations_conf.h"
// ===================================
#ifdef RANS_EQUATIONS
// ==================================

// config files ------------------
#include "User_RANS.h"
#include "RANS_parameters.h"
#include "MGFE_conf.h"
// class files ---------------------
#include "MGSolverDA.h"
#include "TurbUtils.h"


class MGEquationsSystem;



enum BoundaryType
{
  WALL = 0,
  OTHER = 1
};



// =================================================
/// Class for mg energy equation solvers with name TBK_EQUATIONS. Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolRANS:public MGSolDA
{

  // ======================================================================================================
  // ============= MGSolTBK class data ======================================================================
  // ======================================================================================================
public:
  RANS_param _RANS_parameter;
  DYNturModels *_MUTmodel;
protected:
  int _nTdim;			//<dimension
  double _dt;			///< =_mgutils.get_par("dt");
  int _FF_idx[30];		//< field equation flag
  double _euler_impl;
  // parameters--------------------------------------------------------------------------------------------
  // constant reference parameters (From paramater.in file _mgphys.get_par("Uref")) -----------------------
  const double _uref;	  /**< "Uref" */
  const double _lref;			    /**< "lref" */
  const double _Tref;					      /**< "Tref" */
  // constant fluid properties (From paramater.in file _mgphys.get_par("Uref")) ---------------------------
  const double _rhof;	   /**< "rho0" */
  const double _muf;				/**< "mu0" */
  const double _cp0;						 /**< "cp0" */
  const double _kappa0;	       /**< "kappa0" see paramater.in file */
  const int _offset;		///< = _mgmesh._NoNodes[_NoLevels-1]= mesh nodes

  // nondimensional numbers -------------------------------------------------------------------------------
  double _IRe;						  /**< Reynolds*/
  double _explicit_diss[2], _explicit_source[2], _implicit_diss[2],
    _implicit_source[2];
  double _Cross[2], _Log_Cross[2];
  double _Adv, _Lap, _Vel_g[DIMENSION], _MuTurbDxg[DIMENSION], _LapSupg,
    _LapMuTurb;
  double _NodeMuTurb[NDOF_FEM], _NodeWallDist[NDOF_FEM];

  int _TimeDer;
  double _x_1ts[NDOF_FEM], _x_2ts[NDOF_FEM];

  // turbulence
  double _nueff;
  double _kappa_g[2];	   /**< reference kappa*/
  double _y_dist;		///< distance from the wall
  double _vol_frac;
  double _sP;			///< turbulent tensor modulus
  double _mu_turb;

  // mesh -------------------------------------------------------------------------------------------------
  double _xx_qnds[NDOF_FEM * DIMENSION];    /**< elem coords */
  double _xxb_qnds[NDOF_FEMB * DIMENSION];
  // -------------------- class field ---------------------------------------------------------------------
  // element boundary conditions

  // ------------------ integration -----------------------
  //  fields at gaussian points
  double _ub_g[3][30];	   /**< external field  (0-1-2 degree)*/
  double _xxg[DIMENSION];				    /**< gauss pts*/
  double _InvJac2[DIMENSION * DIMENSION];
  double _ub_dxg[2 * DIMENSION];	///< external field derivative  (0-1-2 degree)

  double _WallDist;

  double _InvSigma;
  bool _SolveRANS;
  double _tur_nl[2 * NDOF_FEM];
    std::vector < int >_EquationNodes;

  // STABILIZATION
  bool _SUPG, _ModifiedSupg, _SCPG;
  int _theta;
  double _UPWIND;
  int _NonLinearIt;
  int _NLITER;

  int _InterpolatedMuTurb;
  double _kappa_nl_g[2];
  int _SolveSteady;

  int _WallElement = 0;
  double _NormMidCell[DIMENSION];

  // TURBULENCE MODEL
  double _source[2], _diss[2];
  double _T_nl_dxg[2][DIMENSION], _T_dxg[2][DIMENSION],
    _T_2dxg[2][DIMENSION * DIMENSION];


  int _NodeOnNwallSides[NDOF_FEM];
    std::map < int, BoundaryType > _BoundaryMap;
  int _Restart;

  // BOUNDARY AND INITIAL CONDITION
    std::vector < int >_BoundaryGroupsIDs;	///< vector containing boundary group ids
    std::map < int, bound_condK > _BoundaryConditionMap;	///< map containing the b.c. for every group
  int _bc_vol[NDOF_FEM];	 /**<  b.cond from function */
  int _bc_bd[NDOF_FEM];		 /**< b.cond flags */
  int _bc_el[NDOF_FEM];		 /**<  b.cond in matrix assemblying */
  int _bc_bound[NDOF_FEM];
  bool _BoundElem;		///< flag indicating if an element has at least one boundary side
  bool _FlatProfile;
  int _AxiSym;

  int _ExplicitNearWallDer[2];	// 1: yes, 0:no
public:
  // ==========================================================================
  // =========         Constructor - Destructor  ==============================
  // ==========================================================================
  /// This function constructs the 3d-2D MGSolTBK class
    MGSolRANS (			///< Constructor
		MGEquationsSystem & mg_equations_map_in,	///<  mg_equations_map_in pointer
		const int nvars_in[],	///< KLQ number of variables
		std::string eqname_in = "K",	///< equation name
		std::string varname_in = "k"	///< basic variable name
    );


  // ==========================================================================
  /// This function destructs the 3d-2D MGSolTBK class
//     ~MGSolTBK() {}                              ///< Destructor
   ~MGSolRANS ()
  {
    _BoundaryGroupsIDs.clear ();
    _BoundaryConditionMap.clear ();
  }
  //===========================================================================
  // =========== Read MGSolTBK functions  =======================================
  // ==========================================================================

  // This function reads the boundary conditions
  void bc_read (		///< Read bc
		 int bc_gam /**< group idx  */ ,
		 int bc_mat /**< mat idx      */ ,
		 double xp[] /**< pts coords */ ,
		 int bc_Neu[] /**< bc flags   */ , int bc_value[]
	/**< aux bc flags*/
    );
  // ==========================================================================
  // This function reads the initial solution
  void ic_read (		///< Read bc
		 int bc_gam /**< group idx    */ ,
		 int bc_mat /**< mat id     */ ,
		 double xp[] /**< pts coords   */ ,
		 int iel /**< element id */ ,
		 double u_value[]
	/**< value vector */
    );				///< Read ic


  /// d)  Assemblying MGSolNS Operators

  // ==============================================================================
  // This function assembles the Volume integral
  void GenMatRhs (const double time,	/**< time*/
		  const int Level /**< discrtization Level*/ ,
		  const int mode
		     /**< y/n assembnle rhs */
    );
  // ==========================================================================
  /// This function  computes a time step
  void MGTimeStep (const double time,	///< Time-step manager function
		   const int	/*iter */
    );
  void MGTimeStep_no_up (const double time,	///< Time-step manager function
			 const int	/*iter */
    );
  void MGUpdateStep ();

  void vol_stab (const int el_ndof2,
		 const int el_ngauss, const int mode, int el_conn[]);

  // ==========================================================================
  /// This function  computes the  functional defined (user write}
  void MGFunctional (const double /*time */ ,
		     double	/*starting_distance */
    )
  {
  }

  void SetTurbPar ();
  // ==========================================================================
  /// This function sets the _bc_el with boundary condition flags.
  /// This function also assembles the surface integral to obtain the algebraic sytem.
  /// It is called by MGSolTBK::GenMatRhs in MGSolverT.C
  /// This function sets  the  functional defined (user write}
  void bc_set (int sur_toply[] /**< Local matrix  */ ,
	       int el_ndof2 /**< el dofs*/ , int elb_ndof2 /**< el bd dofs*/ ,
	       int elb_ngauss, int el_conn[]
		  /**< #pt gauss  */
    );
  // ==========================================================================
  /// This function assembles the volume integral to obtain the algebraic sytem.
  /// It is called by MGSolTBK::GenMatRhs in MGSolverT.C
  void vol_integral (const int el_ndof2 /**< el dofs*/ ,
		     const int el_ngauss /**< #pt gauss  */ ,
		     double xx_qnds[] /**< el coords  */ ,
		     int iel, int el_conn[]);

  void wall_element (const int el_ndof2 /**< el dofs*/ ,
		     const int el_ngauss /**< #pt gauss  */ ,
		     double xx_qnds[] /**< el coords  */ ,
		     int iel, int el_conn[]);

  void CalcSCPG (double ParVel[], double &tauc, double DefSource,
		 double vel_g[], double mod2_vel, double h_eff);


  double CalcPhiSupg (int i, double vel_g[], double ParVel[], double tauc,
		      double f_upwind, double implicit_diss, int el_ndof2);


  void FractionalTimeStep (const double time);
  void StandardTimeStep (const double time);

  virtual void CalcSourceAndDiss (int el_ndof2);
  virtual void CalcAdvectiveAndDiffusiveTerms (int TestFunctionID,
					       int InterpolationNode,
					       int NumOfNodes,
					       double f_upwind);
  virtual void VelocityForSUPG (double &mod2_vel, double vel_g[],
				double VEL[]);

  void CalcSUPG (double &h_eff, double &f_upwind, double mod2_vel,
		 double vel_g[]);

  void VelOnGaussAndTensorModulus (double &mod2_vel, int NumOfNodes);
  void FillBoundaryMap ();
  void SetUpFlags ();
};




#endif // define TBK_EQUATIONS
#endif //__mgsolverT_h__
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
