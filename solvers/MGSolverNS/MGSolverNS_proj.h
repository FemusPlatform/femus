#ifndef __mgsolverns_proj_h__
#define __mgsolverns_proj_h__

#include "Equations_conf.h"
// =================================
#ifdef NS_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"
#include "MGSolverNS.h"
#include "UserNS.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;

// ===========================================================================
//                                Navier-Stokes equation class
//=============================================================================
/// Navier-Stokes equation  (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolNS_proj:public MGSolNS
{
/// Class for mg Navier-Stokes equation  with name NS_EQUATIONS.
/// Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)

// ===========================================================================
private:


    // ======================  start functions ===================================
public:

    /// b)   Init MGSolNS functions (constructor,destructor,external fields):
    MGSolNS_proj (   ///< Constructor
        MGEquationsSystem & mg_equations_map, ///< equation map class (Mesh and parameters)
        int nvars_in[],   ///< KLQ number of variables
        std::string eqname_in = "NS0",    ///< base name system
        std::string varname_in = "u"  ///< base name variable
    );
    ~MGSolNS_proj () {};               ///< Destructor

    void SetVariableNames ( std::string varname_in );
    void get_el_field_data ( int iel, int Level, int el_conn[], int offset, int el_dof[], int ndof_lev );

    inline int GetIndxEquation ( int nPhi, int rowShift, int el_ndof[] )
    {
        return nPhi + _dir * el_ndof[2];
    };

    // ADD TO KEM OR FEM
    void MomentumEquation ( double JxW_g2, int el_ndof[], int qp );

    // BOUNDARY
    void set_bc_matrix ( int dir_maxnormal,    ///<  normal dir
                         int sur_toply[],  ///< boundary topology map
                         int el_ndof[],    ///< number of volume dofs
                         int elb_ndof[],   ///< number of boundary dofs
                         int elb_ngauss,   ///<  number of surface gaussian points
                         double normal[],  ///< normal
                         int el_conn[] );
};


#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; ;
