#ifndef __mgsolver_wilcox__
#define __mgsolver_wilcox__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_EQUATIONS
// ==================================

// config files ------------------
#include "MGSolverRANS.h"
#include "MGSolverDA.h"
class MGEquationsSystem;




// =================================================
/// Class for mg energy equation solvers with name TBK_EQUATIONS. Multilevel and mulitporcessor class (see <a href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolWilcox:public MGSolRANS
{
public:

  MGSolWilcox (			///< Constructor
		MGEquationsSystem & mg_equations_map_in,	///<  mg_equations_map_in pointer
		const int nvars_in[],	///< KLQ number of variables
		std::string eqname_in = "K",	///< equation name
		std::string varname_in = "k"	///< basic variable name
    );

   ~MGSolWilcox ()
  {
  }
  void CalcSourceAndDiss (int el_ndof2);
  void CalcAdvectiveAndDiffusiveTerms (int TestFunctionID,
				       int InterpolationNode, int NumOfNodes,
				       double f_upwind);
  void VelocityForSUPG (double &mod2_vel, double vel_g[], double VEL[]);
};

#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
