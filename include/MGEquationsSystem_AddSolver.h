#ifndef __mgequationsmap_add_solver_h__
#define __mgequationsmap_add_solver_h__

template <class SolverClass>
void MGEquationsSystem::AddSolver(
    std::string SystemName, int nSys, int nPieceWise, int nLinear, int nQuadratic, std::string VarName) {
  int nvars_in[3];
  nvars_in[0] = nPieceWise;
  nvars_in[1] = nLinear;
  nvars_in[2] = nQuadratic;
  SolverClass* mgs = new SolverClass(*this, nvars_in, SystemName, VarName);
  set_eqs(mgs);
  set_num_eqs(mgs->_eqname, nSys);
  return;
}

template <class SolverClass>
void MGEquationsSystem::AddSolver(
    std::string SystemName, int nSys, int nPieceWise, int nLinear, int nQuadratic, std::string VarName,
    std::vector<FIELDS> PBname) {
  int nvars_in[3];
  nvars_in[0] = nPieceWise;
  nvars_in[1] = nLinear;
  nvars_in[2] = nQuadratic;
  SolverClass* mgs = new SolverClass(*this, nvars_in, SystemName, VarName);
  set_eqs(mgs);
  set_num_eqs(mgs->_eqname, nSys);
  mgs->set_ext_fields(PBname);
  return;
}

#endif
