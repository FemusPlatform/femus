
#ifndef __userCTRL_h__
#define __userCTRL_h__

#include "Equations_conf.h"
// ===================================
#ifdef CTRL_EQUATIONS
// ==================================
#include <string>
#include <map>
#include <vector>
#include "Solvertype_enum.h"
class MGUtils;

  enum bound_condCTRL {
    fix_in0_ctrl=      1,     ///<  1  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$
    fix_tg0_ctrl=      2,     ///<  2  normal velocity inlet  \f$ {\bf v}.{\bf t}={ 0}{\bf t}  \f$
    fix_disp0_ctrl=              3,     ///<  3  normal velocity inlet  \f$ {\bf v}={\bf 0}  \f$
    fix_in_ctrl=       5,     ///<  5  normal velocity inlet  \f$ {\bf v}.{\bf n}={\bf v}_0.{\bf n}  \f$
    fix_tg_ctrl=       6,     ///<  6  normal velocity inlet  \f$ {\bf v}.{\bf t}={\bf v}_0.{\bf n}  \f$
    fix_dis_ctrl=          7,     ///<  7  normal velocity inlet  \f$ {\bf v}={\bf v}_0  \f$
    free_disp_ctrl=          10,     ///< 10  Neuman homogeneus         (dT.n=0)
    inter_ctrl=         11,     ///<       
    free_disp_outlet_ctrl=  12,     ///< 12  Neuman homogeneus         (dT.n=0)
     free_wall_turb_ctrl=        13,     ///< 13   Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
     free_disp_p_ctrl=        14,     ///< 14  Neuman nonhomogeneus      (dT.n=q_0)
     free_disp_inlet_ctrl=   16,     ///< 16  Neuman nonhomogeneus      (dT.n=q_0)
    
    simm_dispx_ctrl=            21,     ///< 21  simmetry      \f$ (p{\bf n}+{\bf tau} \cdot \widehat{n})\cdot \widehat{i}_x=0 \f$   
    simm_dispy_ctrl=            22,     ///< 22  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simm_dispz_ctrl=            23,     ///< 23  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simm_dispxy_ctrl=           24,     ///< 24  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simm_dispxz_ctrl=           25,     ///< 25  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
    simm_dispyz_ctrl=           26      ///< 26  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)     
};


// ============================================================================
/*!  \defgroup CTRL_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup CTRL_param
// ============================================================================
class CTRL_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    int _SolveSteady;                    /// Flag for steady state solution - value: 0 or 1
    int _Supg;                           /// Flag for SUPG stabilization of advection term - value: 0 or 1
    double   _Upwind;                         /// Flag for normal upwind stabilization - value: 0 <Up <1
    double _UnderRelaxation;             /// Flag for stabilization based on UnderRelaxation - value between 0 and 1
    int _ReactionNumberBased;            /// Flag for stabilization based on reaction number - value: 0 or 1
    int _FlatProfile;                    /// Flag for initial solution: flat profile or modulated as function of wall distance - value: 0 or 1
    bool _Solve;                         /// Flag for solution of temperature equation
    double _Prt;
    double _ctrl_alpha;
    double _ctrl_beta;
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_condCTRL>  _map_CTRLgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_condCTRL> _BoundMap;  /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    
   
public:
    // constructor --------------------------------------------------------------
    CTRL_param() {// SETTING DEFAULT VALUES
        _SolverType                  =GMRESM;
        _SolveSteady                 =0;
        _Supg                        =1;
        _Upwind                      =0;
        _ReactionNumberBased         =0;
        _FlatProfile                 =1;
        _UnderRelaxation             =0.;
        _Prt                         =0.85;
//         _BoundMap["fix_in0"]          =   Twall0    ;
//         _BoundMap["Twall"]           =   Twall     ;
//         _BoundMap["insulation"]      =   insulation;
//         _BoundMap["heat_flux"]       =   heat_flux ;
//         _BoundMap["robinT"]          =   robinT    ;
//         _BoundMap["ext_conv"]        =   ext_conv  ;
//         _BoundMap["simmetry"]        =   simmetry  ;
        
        
        
                                         
        _BoundMap["fix_in0_ctrl"]              =    fix_in0_ctrl        ; 
        _BoundMap["fix_tg0_ctrl"]              =    fix_tg0_ctrl        ;
        _BoundMap["fix_disp0_ctrl"]            =    fix_disp0_ctrl      ;
        _BoundMap["fix_in_ctrl"]               =    fix_in_ctrl         ;
        _BoundMap["fix_tg_ctrl"]               =    fix_tg_ctrl         ;
        _BoundMap["fix_dis_ctrl"]              =    fix_dis_ctrl        ;
        _BoundMap["free_disp_ctrl"]             =   free_disp_ctrl        ;
        _BoundMap["inter_ctrl"]                =   inter_ctrl            ;
        _BoundMap["free_disp_outlet_ctrl"]     =   free_disp_outlet_ctrl; 
        _BoundMap["free_wall_turb_ctrl"]       =   free_wall_turb_ctrl  ;
        _BoundMap["free_disp_p_ctrl"]         =   free_disp_p_ctrl      ;
        _BoundMap["free_disp_inlet_ctrl"]      =   free_disp_inlet_ctrl ; 
        _BoundMap["simm_dispx_ctrl"]           =   simm_dispx_ctrl      ;
        _BoundMap["simm_dispy_ctrl"]           =   simm_dispy_ctrl      ;
        _BoundMap["simm_dispz_ctrl"]           =   simm_dispz_ctrl     ; 
        _BoundMap["simm_dispxy_ctrl"]          =   simm_dispxy_ctrl    ;
        _BoundMap["simm_dispxz_ctrl"]          =   simm_dispxz_ctrl    ;
        _BoundMap["simm_dispyz_ctrl"]          =   simm_dispyz_ctrl    ;       
    }// END T_param CONSTRUCTOR

    ~CTRL_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_CTRLgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

 


// class DS_param
// {//< This class defines the physical and numerical  energy equation parameters
// public:
//     
//     
//     int AXISYM  /**< axisymmetry  */;
//     // stabilization NS --------------------------------------------------------
//      double SUPG/**< supg */; double UPWIND;/**< Normal Upwind */ double UPWIND2;/**< Transv Upwind */
//      
//     // turbulence
//     double LES/**< les */; double DIST_FIX/**< distance from wall */;
//     
//    // non linear -------------------------------------------------------------
//     int    NL_ITER;/**< NON LIN IT */    int    NL_ITER0;/**< INIT NON LIN IT */
//     double NL_TIME0; // initial time for  non linear regime
//     double H_EXT/**< h convective in bc  */ ;double V_EXT/**< v convective in bc  */;
//     
//     // time discretization
//     double CRANK_NICK/**< implicit (1) explicit (0) Crank-Nicolson (0.5)*/; 
//     int UNSTEADY;    /**< un(1)/steady(0) flag  */
//   int COMPRESSIBLE;
//      
//     DS_param(){
//         AXISYM=0;
//          // stabilization NS -----------------------------------------------------
//         SUPG=1.0/**< SUPG */; UPWIND=.0/**< Normal Upwind */; UPWIND2=.0/**< Transv Upwind */;
//         
//         // turbulence
//         LES=.0/**< LES */;DIST_FIX=1.e-3/**< distance from wall */;
//         // non linear ------------------------------------------------------------
//         
//         NL_ITER=0;/**< NON LIN IT */ NL_ITER0=0;/**< INIT NON LIN IT */
//         NL_TIME0=-.0001;
//         
//         // bc  -------------------------------------------------------------
//         H_EXT=1.;/**< h convective in bc  */ V_EXT=1.;/**< v convective in bc  */
//         
//           // time discretization
//         CRANK_NICK=1.;  /**< implicit (1) explicit (0) Crank-Nicolson (0.5) */
//         UNSTEADY=1;     /**< un(1)/steady(0) flag  */
//         COMPRESSIBLE=0;
//         
//     }
//     inline void set_AXISYM(int val){AXISYM=val;}
//     inline void set_SUPG(double    val) {SUPG=val;}
//     inline void set_UPWIND(double  val) {UPWIND=val;}
//     inline void set_UPWIND2(double val) {UPWIND2=val;}
//     inline void set_UNSTEADY(int val){UNSTEADY=val;}
//     
// };



// // NSboundary conditions===============================================================================================
// enum bound_condDS {
//     simm_disp=0, disp_in0=1,disp_tg0=2,wall_fix=3,
//     disp_in=5,disp_tg=6,wall_disp=7,
//     free_disp=10, free_dispn=12,int_disp=11,slip=13,
// };
// NS boundary conditions=======================================================
/*!     \defgroup CTRL_Boundary_conditions     Enum Table: NS Boundary conditions  */ 
/// \ingroup CTRL_Boundary_conditions
    // ========================================================================
    /// Navier-Stokes boundary conditions
//     enum bound_condDS {
//     fix_in0=      1,     ///<  1  normal velocity inlet  \f$ {\bf v}.{\bf n}= 0  \f$
//     fix_tg0=      2,     ///<  2  normal velocity inlet  \f$ {\bf v}.{\bf t}={ 0}{\bf t}  \f$
//     fix_disp0=              3,     ///<  3  normal velocity inlet  \f$ {\bf v}={\bf 0}  \f$
//     fix_in=       5,     ///<  5  normal velocity inlet  \f$ {\bf v}.{\bf n}={\bf v}_0.{\bf n}  \f$
//     fix_tg=       6,     ///<  6  normal velocity inlet  \f$ {\bf v}.{\bf t}={\bf v}_0.{\bf n}  \f$
//     fix_dis=          7,     ///<  7  normal velocity inlet  \f$ {\bf v}={\bf v}_0  \f$
//     
//     free_disp=          10,     ///< 10  Neuman homogeneus         (dT.n=0)
//     inter=         11,     ///<       
//     free_disp_outlet=  12,     ///< 12  Neuman homogeneus         (dT.n=0)
//  
//      free_wall_turb=        13,     ///< 13   Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
//      free_disp_p=        14,     ///< 14  Neuman nonhomogeneus      (dT.n=q_0)
//      free_disp_inlet=   16,     ///< 16  Neuman nonhomogeneus      (dT.n=q_0)
//     
//     simm_dispx=            21,     ///< 21  simmetry      \f$ (p{\bf n}+{\bf tau} \cdot \widehat{n})\cdot \widehat{i}_x=0 \f$   
//     simm_dispy=            22,     ///< 22  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispz=            23,     ///< 23  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispxy=           24,     ///< 24  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispxz=           25,     ///< 25  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
//     simm_dispyz=           26      ///< 26  simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)     
// };
    
  
  
#endif
#endif
