#ifndef __mgsolverCTRL_h__
#define __mgsolverCTRL_h__

#include "Equations_conf.h"
// ===================================
#ifdef CTRL_EQUATIONS
// ==================================

// config files ------------------
#include "UserCTRL.h"
#include "MGFE_conf.h"
// class files ---------------------
#include "MGSclass_conf.h"
// Local Includes -----------------
#include "MGSolverDA.h"


// Forward declarations ----------
class MGEquationsSystem;




// =================================================
/// Class for mg energy equation solvers with name T_EQUATIONS. Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolCTRL: public MGSolDA
{

    // ======================================================================================================
    // ============= MGSolCTRL class data ======================================================================
    // ======================================================================================================
public:
    CTRL_param   _CTRL_parameter;

private:
    int   _nDSdim;  //<dimension
    double       _dt;            ///< =_mgutils.get_par("dt");
    int _FF_idx[30];    //< field equation flag
    // parameters--------------------------------------------------------------------------------------------
    // constant reference parameters (From paramater.in file _mgphys.get_par("Uref")) -----------------------
    const double _uref; /**< "Uref" */ const double _lref;/**< "lref" */ const double _Tref;/**< "Tref" */
    // constant fluid properties (From paramater.in file _mgphys.get_par("Uref")) ---------------------------
    const double _rhof;  /**< "rho0" */   const double _muf;  /**< "mu0" */ 
    const double _rhos;   ///< solid density 
    // nondimensional numbers -------------------------------------------------------------------------------
    
    int    _dir;          ///< direction
    double _disp_d;
    double _alpha; /**< conducibility*/   
    double _IPrdl;/**< Prandl*/             double _IRe; /**< Reynolds*/
    
    double _qheat;  /** vol source*/        double    _qs;/**< surface flux*/

    // turbulence
    double _alpha_turb;  /**< turb conducibility*/ double _IPrdl_turb; /**< turb Prandl number*/
    double _kappa_g[2];  /**< reference kappa*/    double _kappaT_g[2];    /**< reference omega */
    double _y_dist;     ///< distance from the wall
    double _sP;         ///< turbulent tensor modulus
    double _nut_ratio;  ///< effective turbulent viscosity
    double _ctrl_alpha;
    double _ctrl_beta;
      
    // mesh -------------------------------------------------------------------------------------------------
    const int    _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1]= mesh nodes
    double _xx_qnds[NDOF_FEM *DIMENSION]; /**< elem coords */
    double _xxb_qnds[NDOF_FEMB *DIMENSION];
    double _xyz_g[DIMENSION]; 
    double _dphijdx_g2[DIMENSION];
    double _dphijdx_g1[DIMENSION];
    double _dphiidx_g2[DIMENSION];
    double _dphiidx_g1[DIMENSION];
    // -------------------- class field ---------------------------------------------------------------------
    // element boundary conditions
    int   _bc_vol[NDOF_FEM]; /**<  b.cond from function */ 
    int   _bc_bd[NDOF_FEM]; /**< b.cond flags */
    int   _bc_el[NDOF_FEM];  /**<  b.cond in matrix assemblying */
    // ------------------ integration -----------------------
    //  fields at gaussian points
    double  _ub_g[3][30];/**< external field  (0-1-2 degree)*/ 
    double _xxg[DIMENSION];/**< gauss pts*/
 
    double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree)
    double _InvJac2[DIMENSION *DIMENSION];
    bool _SolveCTRL = true;

public:
    // ==========================================================================
    // =========         Constructor - Destructor  ==============================
    // ==========================================================================
    /// This function constructs the 3d-2D MGSolCTRL class
    MGSolCTRL (                                  ///< Constructor
        MGEquationsSystem &mg_equations_map_in, ///<  mg_equations_map_in pointer
        const int nvars_in[],                   ///< KLQ number of variables
        std::string eqname_in="CTRL",              ///< equation name
        std::string varname_in="c"              ///< basic variable name
    );
    // ==========================================================================
    /// This function destructs the 3d-2D MGSolCTRL class
    ~MGSolCTRL() {}                              ///< Destructor

    //===========================================================================
    // =========== Read MGSolCTRL functions  =======================================
    // ==========================================================================
    
    // This function reads the boundary conditions 
    void bc_read ( ///< Read bc
        int    bc_gam/**< group idx  */, int bc_mat    /**< mat idx      */,
        double xp[]  /**< pts coords */,
        int bc_Neu[] /**< bc flags   */, int bc_value[]/**< aux bc flags*/
    );
    // ==========================================================================
    // This function reads the initial solution
    void ic_read ( ///< Read bc
        int    bc_gam    /**< group idx    */, int bc_mat  /**< mat id     */,
        double xp[]      /**< pts coords   */, int iel     /**< element id */,
        double u_value[] /**< value vector */
    );    ///< Read ic
    /// d)  Assemblying MGSolNS Operators
    
    // ==============================================================================
    // This function assembles the Volume integral
    void GenMatRhs (
        const double time, /**< time*/ const int Level /**< discrtization Level*/,
        const int mode     /**< y/n assembnle rhs */
    );
    // ==========================================================================
     /// This function  computes a time step 
    void MGTimeStep (
        const double time,           ///< Time-step manager function
        const int /*iter*/
    );
    // ==========================================================================
    /// This function  computes the  functional defined (user write} 
    double  MGFunctional(double p1,  double & p2);
    // ==========================================================================
    /// This function sets the _bc_el with boundary condition flags.
    /// This function also assembles the surface integral to obtain the algebraic sytem.
    /// It is called by MGSolCTRL::GenMatRhs in MGSolverDS.C
    /// This function sets  the  functional defined (user write}
    void set_bc_matrix(
        DenseMatrixM &KeM, DenseVectorM &FeM,       ///< local Matrix and rhs
        int dir_maxnormal,                          ///<  normal dir
        int sur_toply[],                            ///< boundary topology map
        int el_ndof[],                              ///< number of volume dofs
        int elb_ndof[],                             ///< number of boundary dofs
        int elb_ngauss,                             ///<  number of surface gaussian points
        double normal[],                            ///< normal
        const int iaxis                              ///< axisymmetric
    );
    // ==============================================================================================
void matrixrhsvol_liq_ds(
  DenseMatrixM &KeM,
  DenseVectorM &FeM,
  int el_ndof[],
  int flag_group[]
);// ======

#ifdef TBK_EQUATIONS
    void  f_mu ( double val[] );
#endif

};

#endif  // define CTRL_EQUATIONS
#endif //__mgsolverCTRL_h__

