#ifndef __mgsolver_nagano_keT__
#define __mgsolver_nagano_keT__

#include "Equations_conf.h"
// ===================================
#ifdef RANS_THERMAL_EQUATIONS
// ==================================

// config files ------------------
#include "MGSolverDA.h"
#include "MGSolverRANS_thermal.h"
class MGEquationsSystem;

// =================================================
/// Class for mg energy equation solvers with name TBK_EQUATIONS. Multilevel and mulitporcessor class (see <a
/// href="ns_discretization.pdf"  target="_blank"><b>Overview</b></a>)
class MGSolNaganoKET : public MGSolRANS_thermal {
 public:
  MGSolNaganoKET(                              ///< Constructor
      MGEquationsSystem& mg_equations_map_in,  ///<  mg_equations_map_in pointer
      const int nvars_in[],                    ///< KLQ number of variables
      std::string eqname_in = "K",             ///< equation name
      std::string varname_in = "k"             ///< basic variable name
  );

  ~MGSolNaganoKET() {}

  void CalcAdvectiveAndDiffusiveTerms(
      int TestFunctionID, int InterpolationNode, int NumOfNodes, double f_upwind);
  void VelocityForSUPG(double& mod2_vel, double vel_g[], double VEL[]);

  void inline CalcSourceAndDiss() {
    // first equation source and diss -> kh or log_kh
    _explicit_source[0] = _source[0];
    _explicit_diss[0] = 0.;
    _implicit_diss[0] = _diss[0];
    // second equation source and diss -> eh or wh or log_wh
    _explicit_source[1] = _source[1] + _mecc_term[1];
    _explicit_diss[1] = 0.;
    _implicit_diss[1] = (_diss[1] + _mecc_term[0]);
  };
};

#endif
#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
