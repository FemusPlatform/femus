
#ifndef __userT_h__
#define __userT_h__

#include "Equations_conf.h"
// ===================================
#ifdef T_EQUATIONS
// ==================================
#include <string>
#include <map>
#include <vector>
#include "Solvertype_enum.h"

class MGUtils;

// T boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
enum bound_condT {
    // ========================================================================
    Twall0=     0, ///< 0= Dirichlet homogeneus     (T=0)
    Twall=      1, ///< 1= Dirichlet nonhomogeneus  (T=T_0)
    insulation=10, ///< 10= Neuman homogeneus         (dT.n=0)
    heat_flux= 11, ///< 11= Neuman nonhomogeneus      (dT.n=q_0)
    robinT=    12, ///< 12= Robin    \f$ \nabla T \cdot \widehat{n}= \beta T \f$  (dT.n=beta*T)
    ext_conv=  13, ///< 13= Turbulence      (dT.n=tau_0^2=beta*ug(u))
    simmetry=  20  ///< 20= Simmetry      \f$ \nabla T \cdot \widehat{n}=0 \f$    (dT.n=0)
               // ========================================================================

};

/*!  \defgroup T_param   Class Table:  energy equation parameters (T_param) */
/// \ingroup T_param
// ============================================================================
class T_param
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
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_condT>  _map_Tgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_condT> _BoundMap;  /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    
   
public:
    // constructor --------------------------------------------------------------
    T_param() {// SETTING DEFAULT VALUES
        _SolverType                  =GMRESM;
        _SolveSteady                 =0;
        _Supg                        =1;
        _Upwind                      =0;
        _ReactionNumberBased         =0;
        _FlatProfile                 =1;
        _UnderRelaxation             =0.;
        _Prt                         =0.85;
        _BoundMap["Twall0"]          =   Twall0    ;
        _BoundMap["Twall"]           =   Twall     ;
        _BoundMap["insulation"]      =   insulation;
        _BoundMap["heat_flux"]       =   heat_flux ;
        _BoundMap["robinT"]          =   robinT    ;
        _BoundMap["ext_conv"]        =   ext_conv  ;
        _BoundMap["simmetry"]        =   simmetry  ;
    }// END T_param CONSTRUCTOR

    ~T_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_Tgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};

#endif
#endif
