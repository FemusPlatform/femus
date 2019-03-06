#ifndef __userP_h__
#define __userP_h__
// ============================================================================
#include "Equations_conf.h"
#include "Solvertype_enum.h"
// ============================================================================
// P boundary conditions=======================================================
/*!     \defgroup Boundary_conditions     Enum Table: Boundary conditions  */
/// \ingroup Boundary_conditions
// P boundary conditions===========================================================================
enum bound_cond_p {

// ================================================================================================
// outflowp0 = 0= Dirichlet homogeneus     \f$ p= 0   \f$    (p=0)
// outflowp  = 4= Dirichlet nonhomogeneus  \f$ p= p_0 \f$   (p=p_0)
//
// vel_fix   =10= Neuman homogeneus or simmetry \f$ \nabla p \cdot \widehat{n}= 0 \f$  (dp.n=0)
// interiorp =11= Neuman nonhomogeneus \f$ \nabla p \cdot \widehat{n}= \Delta p_0 \f$   (dp.n=dp0)
// ================================================================================================
// Dirichlet
    outflowp0 =0,///< 0=Dirichlet homogeneus     \f$ p= 0   \f$    (p=0)
    outflowp  =4,///< 1=Dirichlet nonhomogeneus  \f$ p= p_0 \f$    (p=p_0)
// Neuman
    vel_fix  =10,///< 10= Neuman homogeneus or simmetry \f$ \nabla p \cdot \widehat{n}= 0 \f$ (dp.n=0)
    interiorp=11 ///< 11= Neuman nonhomogeneus \f$ \nabla p \cdot \widehat{n}=\Delta p_0 \f$  (dp.n=dp0)
};


class P_param
{
    //< This class defines the physical and numerical  energy equation parameters
public:
    std::vector<int>  _BoundaryGroupsIDs ;         /// Vector containing boundary group ids
    std::map<int,  bound_cond_p>  _map_NSgroup;      /// Map containing boundary group ids and their relative boundary condition
    std::map<std::string, bound_cond_p> _BoundMap;   /// Map that associates a bound_condT condition to the relative string
    std::map<std::string,std::string> _FileMap;    /// String map containing Tproperties.in parameters
    SolverTypeM _SolverType;             // for other solver types see Solvertype_enum.h
    int _TimeDisc;
    int _AssembleOnce;
    int _NodeIDrefPressure;
    P_param() {
        _SolverType              =GMRESM;
        _BoundMap["interior"]        = interiorp;
        _BoundMap["nostress"]        = interiorp;
        _BoundMap["outflow"]         = interiorp;
        _BoundMap["pressure_outlet"] = outflowp;
        _BoundMap["outflow_p"]       = outflowp;
        _BoundMap["pressure_inlet"]  = outflowp0;
        _BoundMap["slip"]            = vel_fix;
        _BoundMap["wall"]            = vel_fix;
        _BoundMap["penalty_turb"]    = vel_fix;
        _BoundMap["velocity"]        = vel_fix;
        _BoundMap["velocity_norm"]   = vel_fix;
        _BoundMap["velocity_tang"]   = vel_fix; 
        _BoundMap["accelerating_swirl"]     = vel_fix;
        _BoundMap["decelerating_swirl"]     = vel_fix;
        _BoundMap["accelerating_stress"]    = vel_fix;
        _BoundMap["decelerating_stress"]    = vel_fix;
        _BoundMap["stress"]                 = vel_fix;
        _BoundMap["swirl"]                  = vel_fix;
    }
    
   ~P_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_NSgroup.clear();
        _BoundaryGroupsIDs.clear();
    };
    
    void  read_param ( MGUtils &mgutils, int proc=0 ); /// This function sets all the T_param parameters
    void  read_file();                     /// This function reads the parameters from Tproperties.in file
    void  print_par();                     /// This function prints the parameters contained in Tproperties.in
};



#endif   //    __userP_h__  <-----------------------------------------------------------------------------------------
