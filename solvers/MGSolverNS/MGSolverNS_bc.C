// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS

#include <iomanip>       // std::setprecision
#include "MGFE.h"        // Mesh class
#include "MGSolverNS.h"  // Navier-Stokes class header file

void MGSolNS::NodeBCFlags(int NodeElemNumber, int& bc_node, int& bc_tang, int& bc_norm, int& bc_var_normal) {
  bc_node = abs(_bc_bd[NodeElemNumber]) % 100;
  const int bc_var_check2 = abs(_bc_bd[NodeElemNumber]);  // total  bc_var
  const int bc_var_npt = bc_var_check2 / 10000;           // recurrent  points in one element
  const int bc_var_check = bc_var_check2 % 10000;         // bc_var
  bc_var_normal = bc_var_check / 1000 - 1;                //
  const int bc_var = bc_var_check % 1000;
  bc_norm = bc_var / 10;
  bc_tang = bc_var % 10;

  return;
}

#endif
