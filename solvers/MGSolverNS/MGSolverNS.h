#ifndef __mgsolverns_h__
#define __mgsolverns_h__

#include "Equations_conf.h"
// =================================
#ifdef NS_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"
// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;

/// @image html swirl2.png

#include "UserNS.h"

enum TANG_CALC_TYPE{
  GEOM_BASED = 0,  
  VEL_BASED  = 1  
};

enum TIME_DER_TYPE{
  OLDER_STEPS = 0,
  ACTUAL_STEP = 1  
};

// ===========================================================================
//                                Navier-Stokes equation class
//=============================================================================
/// Navier-Stokes equation  (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolNS:public MGSolDA
{
/// Class for mg Navier-Stokes equation  with name NS_EQUATIONS.
/// Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)



// ===========================================================================
protected:
    // =====================  start data =======================================

    int _nNSdim;          //<dimension
    int _FF_idx[30];      //< field equation flag
    const double _dt;     ///< =_mgutils.get_par("dt");
    //! Flag to highlight boundary nodes
    /*!
     * For each element, before performing the loop over boundary sides with #MGSolNS::set_bc_matrix,
     * we set the values of #MGSolNS::_bc_el equal to _BdFlagId in order to highlight boundary nodes
     * and to see if an equation has already been assigned to a matrix row. After looping over boundary sides
     * no row can have a value of _bc_el equal to _BdFlagId. See #MGSolNS::_bc_el for admissible values
     */
    const int _BdFlagId = -8;
    /// a) MGSolNS element data
    // mesh
    const int _offset;        ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
    // constant reference parameters and  fluid properties (=_mgutils.get_par(" "))
    const double _uref /**< "Uref" */ ;
    const double _lref /**< "Lref" */ ;
    const double _rhof /**< "rhof" */ ;
    const double _muf /**< "muf"  */ ;
    // nondimensional numbers
    double _IRe, _IReEff /**< Reynolds number */ ;
    double _IFr /**<  =Froud number */ ;
    double _dirg[3] /**< =gravity */ ;
    // turbulence
    double _kappa_g[2] /**< reference kappa */ ;
    double _Wall_dist;    /**< fixed wall distance of first mesh points */
    double _mu_turb /**< turb eff */ ;
    double _sP /**< tensor modulus*/ ;
//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
    int _bc_vol[NDOF_FEM * ( DIMENSION + 1 )]; ///<  element  b.cond flags (vol int)
    int _bc_bd[NDOF_FEM * ( DIMENSION + 1 )]; ///<  element  b.cond flags  (bd int)

    double _u_1ts[NDOF_FEM*DIMENSION], _u_2ts[NDOF_FEM*DIMENSION], _u_nl[NDOF_FEM*DIMENSION];
    double _p_1ts[NDOF_P], _p_2ts[NDOF_P];
    
    double _Pressure[NDOF_P],   _P_gdx[DIMENSION];
    double _u_1ts_g[DIMENSION], _u_2ts_g[DIMENSION], _u_nl_g[DIMENSION], _vel_gdx[DIMENSION*DIMENSION], _vel_gddx[DIMENSION * DIMENSION * DIMENSION];
    
    double _Phi_supg, _f_upwind;
    double _ElemVolume;
    int _Coupled;
    
    //! Flag for equation assembly
    /*! This flag sets the equation that needs to be assembled on each element node.
     * The possible values are:<br>
     * -1, -2 -> projection along tangent directions t1 and t2 with stress calculation <br>
     * -3     -> projection along normal direction with stress calculation             <br>
     *  0     -> usual equation discretization, interior nodes                         <br>
     *  1,  2 -> projection along tangent directions t1 and t2 for Dirichlet bc        <br>
     *  3     -> projection along normal direction for Dirichlet bc
     */
    int _bc_el[NDOF_FEM * ( DIMENSION + 1 )];

    int _bc_bound[NDOF_FEM];
    bool _SolveNS;
    double _InvJac2[DIMENSION * DIMENSION];
    double _InvJac1[DIMENSION * DIMENSION];
    double _xx_qnds[NDOF_FEM * DIMENSION];
    double _xxb_qnds[NDOF_FEMB * DIMENSION];
    double _xyzg[DIMENSION];

    int _Restart = 0;
    int _TimeStep=0;
    //! Node normal vectors
    /*!
     * Array containing normal vectors for each element node. It is used to calculate
     * tangential vectors and to project the equation along normal or tangential directions
     * on nodes laying on boundary sides
     */
    double _normal_pt[NDOF_FEM * DIMENSION];
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
    double _ub_g[3][12];      ///< external field  (0-1-2 degree)
    double _ub_dxg[DIMENSION * DIMENSION];    ///< external field derivative  (0-1-2 degree)
    bool _BoundElem;
    int _AxiSym;
    int _pres_order;
    int _NComp;
    int _WallElement;
    
    int _BoundEquation, _NormalEquation;
    double _ProjDir[3][DIMENSION];

    bool _NormalFlag[NDOF_FEM * DIMENSION];
    //! Flag for Dirichlet boundary conditions
    /*!
     * This flag is used in order to see if a Dirichlet boundary condition has already been set
     * on a particular row.
     */
    bool _AlreadyWrittenDirBC[NDOF_FEM * DIMENSION];

    //! Option for tangential direction calculation
    /*!
     * Tangent vector can be calculated by geometric construction, namely, given normal vector \f$ \bf{n} \f$,<br>
     * - Dimension = 3, calculation of tangential vectors \f$ \bf{t_1}, \bf{t_2}\f$: <br>
     *   -# \f$ t_1[0] = -(n[1] + n[2]), t_1[1] = n[0], t_1[2]=n[0] \f$  <br>
     *   -# \f$ \bf{t_2} = \bf{n} \times \bf{t_1} \f$ <br>
     * - Dimension = 2, calculation of tangential vector\f$ \bf{t}\f$: <br>
     *   -# \f$ t[0] = -(n[1]), t[1] = n[0] \f$ <br>
     *
     * or it can be calculated following flow main direction: <br>
     * - Dimension = 3, calculation of tangential vectors \f$ \bf{t_1}, \bf{t_2}\f$: <br>
     *   -# \f$  \bf{t_1} = \frac{[\bf{u} -(\bf{u} \cdot \bf{n})\bf{n}]}{|\bf{u}|}  \f$<br>
     *   -# \f$  \bf{t_2} = \bf{n} \times \bf{t_1}\f$
     */

// ParaMEDMEM::MEDCouplingUMesh *_PalaMesh = NULL;
// ParaMEDMEM::DataArrayDouble  *_CoordArray = NULL;
// ParaMEDMEM::DataArrayDouble  *_RefinedConn = NULL;
    std::vector<int> _AdditionalCellsConn;
    std::vector<int> _CellsToRemove;
    std::vector<double> _AdditionalCellsCoords;
    std::vector<double> _AdditionalCellsVelField;
    int _InitNumberCells, _InitNumberNodes;


    // ======================  end data ==========================================

    // ======================  start functions ===================================
public:
    NS_param _NS_parameter;
//  DenseMatrixM _KeM;
//      DenseVectorM _FeM;
    /// b)   Init MGSolNS functions (constructor,destructor,external fields):
    MGSolNS (   ///< Constructor
        MGEquationsSystem & mg_equations_map, ///< equation map class (Mesh and parameters)
        int nvars_in[],   ///< KLQ number of variables
        std::string eqname_in = "NS0",    ///< base name system
        std::string varname_in = "u"  ///< base name variable
    );
    ~MGSolNS ()
    {
    };                ///< Destructor

    /// c)          Read MGSolNS functions:
    /// This function  read bc
    void ic_read ( int bc_gam, int bc_mat, double xp[],   // point coordinates
                   int iel, double u_value[]   // point field values
                 );

    /// This function reads ic
    void bc_read ( int bc_gam, int bc_mat, double x[], // point coordinates
                   int bc_Neu[],   //  bc volume integral flag
                   int bc_bd[] //  bc surface integral flag
                 );

    /// d)  Assemblying MGSolNS Operators
    /// This function computes volume and surface assemblying
    void GenMatRhs ( const double time,   // time
                     const int Lev,    // Level
                     const int m   // rhs assembly control
                   );
// Multigrid function ------------------------------------------------------------------
    /// This function is the  time-step manager function
    void MGTimeStep ( const double time,  ///< time
                      const int    /*iter *////< number max of iterations
                    );
    
    void MGTimeStep_no_up ( const double time,  ///< time
                      const int    /*iter *////< number max of iterations
                    );
    
    void MGUpdateStep();

    // ==============================================================================================
    virtual void set_bc_matrix ( int dir_maxnormal,    ///<  normal dir
                         int sur_toply[],  ///< boundary topology map
                         int el_ndof[],    ///< number of volume dofs
                         int elb_ndof[],   ///< number of boundary dofs
                         int elb_ngauss,   ///<  number of surface gaussian points
                         double normal[],  ///< normal
                         int el_conn[] ){};

    // ==============================================================================================

    virtual void matrixrhsvol (
        int el_ndof[],
        const int mode,
        int el_conn[]
    );

    // ==============================================================================================

    //! Function for the calculation of tangential vectors
    /*! The function sets the tangential vectors for boundary elements.
     * The normal vector is taken from #MGSolNS::_normal_pt .
     * This function is used for the projection of Navier Stokes equation along
     * tangential directions, in particular for assemblying the volume part
     */
    void CalcTangDir ( int PhiNumber, ///< [in] ID of node for which we are writing the equation - local element numbering
                       int el_dof, ///< [in] Number of element dof (global element, not boundary)
                       int bc_el,  ///< [in] Boundary condition flag
                       int qp, ///< [in] Number of quadrature point
                       double Tang1[], ///< [in,out] Array for storing tangential vector along principal flow direction
                       double Tang2[], ///< [in,out] Array for storing tangential vector along secondary flow direction
                       TANG_CALC_TYPE TangType   ///< [in] Method for calculation of tangent vectors, see #MGSolNS::TangType
                     );

    //! Function for the calculation of tangential vectors
    /*! The function sets the tangential vectors for boundary elements.
     * This function is used for the projection of Navier Stokes equation along
     * tangential directions, in particular for assemblying the boundary part
     */
    void CalcTangDir ( int PhiNumber, ///< [in] ID of node for which we are writing the equation - local element numbering
                       int el_dof, ///< [in] Number of element dof (global element, not boundary)
                       int bc_el,  ///< [in] Boundary condition flag
                       int qp, ///< [in] Number of quadrature point
                       double Tang1[], ///< [in,out] Array for storing tangential vector along principal flow direction
                       double Tang2[], ///< [in,out] Array for storing tangential vector along secondary flow direction
                       double Normal[],    ///< [in] Array containing boundary element normal vector
                       TANG_CALC_TYPE TangType    ///< [in] Method for calculation of tangent vectors, see #MGSolNS::TangType
                     );

    //! Function for the calculation of tangential vectors
    /*! The function sets the tangential vectors for boundary elements.
     * This function is used for the projection of Navier Stokes equation along
     * tangential directions, in particular for assemblying the boundary part
     */
    void CalcTangDir ( double Tang1[], ///< [in,out] Array for storing tangential vector along principal flow direction
                       double Tang2[], ///< [in,out] Array for storing tangential vector along secondary flow direction
                       double Normal[],    ///< [in] Array containing boundary element normal vector
                       TANG_CALC_TYPE TangType    ///< [in] Method for calculation of tangent vectors, see #MGSolNS::TangType
                     );

    //! Function that sets boundary flags
    /*!
     * This function reads the boundary conditions for each element boundary node and sets the condition
     * for normal and tangential directions. The flags are written into #MGSolNS::_bc_el
     */
    void SetBCFlags ( double normal[], int sur_toply[],
                      int el_ndof, int el_conn[] );

    //! Correction of BC flags
    /*!
     * This function corrects the values of #MGSolNS::_bc_el <br>
     * It is used when #MGSolNS::_nNSdim == 3 in order to switch the rows taking the equations along
     * directions \f$ \bf{t_1} \f$ and \f$ \bf{t_2} \f$ for corvengence issues.     <br>
     * \f{bmatrix}{  \dots row = i  \dots \\ \\ \dots row = i + el\_ dof \dots \\ \\ \dots row = i + 2 el\_ dof \dots  \f}
     * Equation for normal direction is written in \f$ row = i + j*el\_ dof \f$ where j is the maximum normal vector component.
     * In a similar way equations for directions \f$ \bf{t_1} \f$ and \f$ \bf{t_2} \f$ are written in \f$ row = i + k*el\_ dof \f$ and
     * \f$ row = i + l*el\_ dof \f$ where k and l are the maximum vector components of \f$ \bf{t_1} \f$ and \f$ \bf{t_2} \f$, respectively.
     */
    void CorrectBCFlags ( double Tang1[], ///< [in, out] Array containing components of vector \f$ \bf{t_1} \f$
                          double Tang2[],  ///< [in, out] Array containing components of vector \f$ \bf{t_2} \f$
                          int MaxNormal,   ///< [in] Maximum normal vector component
                          int ElementNode, ///< [in] Node id, local numbering
                          int ElDof,   ///< [in] Number of element nodes
                          int Dirichlet    ///< [in] Flag for Dirichlet bc (0) or Neumann bc (1)
                        );
    void CalcMuTurbAndWallDist ( double MuTurb[], double WallDist[], int WallNodes[],int ElId, int NumOfNodes, int Level );
    
    void CalcVolume ();
    
    void GetVolumeTestFunctions( int qp );
    virtual void SetVariableNames(std::string varname_in){};
    virtual void get_el_field_data ( int iel, int Level, int el_conn[], int offset, int el_dof[], int ndof_lev){};
    virtual void InterpolateSolutions( int el_dof[] );
    virtual void RowSetUp( int nPhi, int indx_eq, int qp, int el_ndof[]);
    virtual inline int GetIndxEquation( int nPhi, int rowShift, int el_ndof[] ){};
    
    // MOMENTUM EQUATION DISCRETIZATION
    virtual void CalcAdv_Lap_LapSupg( double & Adv, double & Lap, double & LapSupg, int nPhi, int j, int el_ndof[], int dir );
    virtual void CalcPresGrad( double &PresGrad, int nPhi, int interpNodeID, int spaceDir, int el_ndof[] );
    virtual void CalcTimeDer ( double &TimeDerivative, TIME_DER_TYPE Type, int SpaceDir, int nPhi, int interpNodeID=0 );
    
    // CONTINUITY EQUATION DISCRETIZATION
    virtual void ContinuityEquation( double JxW_g2, int el_ndof[] ){};
    virtual void CalcVelDivergence( double &Div, int PresOrder, int nPhi, int InterpNodeID, int SpaceDir, int el_ndof[] );
    
    // ADD TO KEM OR FEM
    virtual void MomentumEquation ( double JxW_g2, int el_ndof[], int qp ){};
    
    
    // BOUNDARY FUNCTIONS 
    void NodeBCFlags (int NodeElemNumber, int &bc_node, int &bc_tang, int& bc_norm, int &bc_var_normal );

    
};


#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;
