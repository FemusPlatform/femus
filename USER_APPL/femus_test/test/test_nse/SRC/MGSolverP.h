#ifndef __mgsolverPb_h__
#define __mgsolverPb_h__

#include "Equations_conf.h"
// ===================================
#ifdef NS_EQUATIONS
// ==================================
#if (NS_EQUATIONS%2==0)
// config files ------------------
#include "UserP.h"
#include "MGFE_conf.h"
// class files ---------------------
#include "MGSclass_conf.h"
// Local Includes -----------------
#include "MGSolverDA.h"


// Forward declarations ----------
class MGEquationsSystem;


// ======================================================================================
// ====================== MGSolFSIP  =================================== =================
// ======================================================================================


// =================================================
/// Class for mg energy equation solvers with name T_EQUATIONS. 
// Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolP: public MGSolDA
{

    // ======================================================================================================
    // ============= MGSolP class data ======================================================================
    // ======================================================================================================
public:
    P_param   _P_parameter;

private:
    int   _nPdim;  //<dimension
    double       _dt;            ///< =_mgutils.get_par("dt");
    int _FF_idx[30];    //< field equation flag
    double _euler_impl;
    // parameters--------------------------------------------------------------------------------------------
    // constant reference parameters (From paramater.in file _mgphys.get_par("Uref")) -----------------------
    const double _uref; /**< "Uref" */ const double _lref;/**< "lref" */ // const double _Tref;/**< "Tref" */

    // constant fluid properties (From paramater.in file _mgphys.get_par("Uref")) ---------------------------
    const double _rhof;  /**< "rho0" */   const double _muf;  /**< "mu0" */ // const double _cp0;/**< "cp0" */
//     const double _kappa0;    /**< "kappa0" see paramater.in file */

    // nondimensional numbers -------------------------------------------------------------------------------
    double _alpha; /**< conducibility*/     double _IPrdl;/**< Prandl*/        double _IRe; /**< Reynolds*/
    double _qheat;  /** vol source*/        double    _qs;/**< surface flux*/

    // turbulence
    double _alpha_turb;  /**< turb conducibility*/ double _IPrdl_turb; /**< turb Prandl number*/
    double _kappa_g[2];  /**< reference kappa*/    double _kappaT_g[2];    /**< reference omega */
    double _y_dist;     ///< distance from the wall
    double _sP;         ///< turbulent tensor modulus
    double _nut_ratio;  ///< effective turbulent viscosity

    // mesh -------------------------------------------------------------------------------------------------
    const int    _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1]= mesh nodes
    double _xx_qnds[NDOF_FEM *DIMENSION]; /**< elem coords */ double _xxb_qnds[NDOF_FEMB *DIMENSION];
    double _xyzg[DIMENSION];
    // -------------------- class field ---------------------------------------------------------------------
    // element boundary conditions
    int   _bc_vol[NDOF_FEM]; /**<  b.cond from function */  int   _bc_bd[NDOF_FEM]; /**< b.cond flags */
    int   _bc_el[NDOF_FEM];  /**<  b.cond in matrix assemblying */
    // ------------------ integration -----------------------
    //  fields at gaussian points
    double  _ub_g[3][30];/**< external field  (0-1-2 degree)*/ double _xxg[DIMENSION];/**< gauss pts*/
    double _InvJac2[DIMENSION *DIMENSION];
    double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree)

    bool _AxiSym;
    bool _SolveP;
public:
    // ==========================================================================
    // =========         Constructor - Destructor  ==============================
    // ==========================================================================
    /// This function constructs the 3d-2D MGSolP class
    MGSolP (                                  ///< Constructor
        MGEquationsSystem &mg_equations_map_in, ///<  mg_equations_map_in pointer
        const int nvars_in[],                   ///< KLQ number of variables
        std::string eqname_in="NS2P",              ///< equation name
        std::string varname_in="p"              ///< basic variable name
    );
    // ==========================================================================
    /// This function destructs the 3d-2D MGSolP class
    ~MGSolP() {}                              ///< Destructor

    //===========================================================================
    // =========== Read MGSolP functions  =======================================
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
    void MGFunctional (
        const double /*time*/,           
        double /*starting_distance*/
    ) {}
    // ==========================================================================
    /// This function sets the _bc_el with boundary condition flags.
    /// This function also assembles the surface integral to obtain the algebraic sytem.
    /// It is called by MGSolP::GenMatRhs in MGSolverT.C
    /// This function sets  the  functional defined (user write} 
//     void bc_set (
//         DenseMatrixM &KeM /**< Local matrix  */,DenseVectorM &FeM/**< Local rhs */,
//         int sur_toply[] /**< Local matrix  */,
//         int el_ndof2 /**< el dofs*/,int elb_ndof2 /**< el bd dofs*/,int elb_ngauss /**< #pt gauss  */,
//         double u_old[] /**< old velocity field  */
//     );
    // ==========================================================================
    /// This function assembles the volume integral to obtain the algebraic sytem.
    /// It is called by MGSolP::GenMatRhs in MGSolverT.C
//     void vol_integral (
//         DenseMatrixM &KeM/**< Local matrix  */, DenseVectorM &FeM/**< Local rhs */,
//         const int el_ndof2/**< el dofs*/,const int el_ngauss/**< #pt gauss  */,
//         double xx_qnds[]/**< el coords  */,
//         const int unsteady/**< un/steady flag  */,const int mode/**<  rhs int flag */
//     );
    // ==========================================================================
    /// This function computes the field values of all fields in the data structure.
    /// The number _FF_index[]*NDOF_FEM sets the location of the fields in the data
    /// vector _data.ub[]
//     void el_vol_data (
//         double xx_qnds[]/**< el coords  */,int el_conn[]/**< connectivity  */,
//         int el_ndof[]/**< el dofs*/,int offset /**< number of nodes on top level */
//     );
    // ==========================================================================
    /// This function computes the effective element measure heff if velocity field vel[] is
    /// given (average on the element)
//     double   heff ( double vel[] );
    
#ifdef TBK_EQUATIONS
    void  f_mu ( double val[] );
#endif


};


#endif //  NS_EQUATIONS%2=0

#endif  // define NS_EQUATIONS
#endif //__mgsolverP_h__

