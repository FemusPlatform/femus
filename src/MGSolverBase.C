// lib include ------------------
// #include <limits>
// #include <sstream>

// configure file -------------
// #include "MGFE_conf.h"
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"
// // class include
#include "MGSolverBase.h"
#ifdef TWO_PHASE
#include "MGSolverCC.h"
#endif
// local inlude -----------------
#include "MGMesh.h"
#include "MGSystem.h"
#include "MGUtils.h"
#include "MGEquationsSystem.h"
// #include "MGFEMap.h"

// alg include-------------------
#include "linear_solverM.h"
#include "numeric_vectorM.h"
#include "sparse_MmatrixM.h"
#include "sparse_matrixM.h"

// #include <mpi.h>  //for MPI_COMM_WORLD
// #include "petsc_macroM.h"

//  #define VANKA (0)
//
// ============================================================================
/// This function  is the MGSolBase constructor :
MGSolBase::MGSolBase(
    MGEquationsSystem& e_map_in,  // equation map
    const int nvars_in[],         // # variables
    std::string eqname_in         // equation name
    )
    : _mgutils(e_map_in._mgutils),  // mgutils pointer from equation map pointer
      _mgfemap(e_map_in._mgfemap),  // mgfemap pointer from equation map pointer
      _mgmesh(e_map_in._mgmesh),    // mgmesh pointer from equation map pointer
      _mgeqnmap(e_map_in),          // equation map pointer
      _NoLevels((int)(e_map_in._mgutils._geometry["nolevels"])),
      _eqname(eqname_in),
      _n_vars(nvars_in[0] + nvars_in[1] + nvars_in[2]),
      _var_names(NULL),
      _refvalue(NULL)
#ifdef TWO_PHASE_LIB
      ,
      _msolcc(NULL)
#endif
{
  // processor number ---------------------------------------------------------
  _iproc = _mgmesh._iproc;  // > 1CPU
  //  Lagrange piecewise  variables (const,linear,quad) -----------------------
  _nvars[0] = nvars_in[0];  // Lagrange piecewise constant variables
  _nvars[1] = nvars_in[1];  // Lagrange piecewise linear variables
  _nvars[2] = nvars_in[2];  // Lagrange piecewise linear variables

  // allocation of dynamic system ---------------------------------------------
  _Dim = new int[_NoLevels];   // matrix and vect  dim
  A.resize(_NoLevels);         // matrix
  x.resize(_NoLevels);         // matrix vect sol
  b.resize(_NoLevels);         // rhs
  res.resize(_NoLevels);       // residual
  x_old[0].resize(_NoLevels);  // old solution 1 step (old)
  x_old[1].resize(_NoLevels);  // old solution 2 step (oold)
  x_old[2].resize(_NoLevels);  // old solution 3 step (ooold)
  x_nonl.resize(_NoLevels);    // non linear solution
  x_aux.resize(0);
  d_aux.resize(0);
  // dof info -----------------------------------------------------------------
  _node_dof = new int*[_NoLevels + 1];  // dof (+1)
  // restr and prol operators -------------------------------------------------
  Rst.resize(_NoLevels);
  Prl.resize(_NoLevels);  // Projector and restrictor
  // solver -------------------------------------------------------------------
  _solver = new LinearSolverM*[_NoLevels];  // one solver for each level
  for (int l = 0; l < _NoLevels; l++)
    _solver[l] = LinearSolverM::build(_mgmesh._comm.comm(), LSOLVER).release();
  // class parameters ----------------------------------------------------------
  _control = 0.;
  _dt = 0.1;  // time step
  _dir = 0;   // number equation for segregated system

  return;
}

// ===================================================
/// This function  is the MGSolBase destructor.
MGSolBase::~MGSolBase() {  // ===================================================
  // clear substructrures
  clear();
  A.clear();
  x.clear();
  b.clear();    //  A and x and b
  res.clear();  //  rhs and residual vector

  x_nonl.clear();  // nonlinear solution tmp
  x_old[0].clear();
  x_old[1].clear();
  x_old[2].clear();  // old solutions
  x_aux.clear();
  d_aux.clear();
  Rst.clear();
  Prl.clear();  // Restrictor and projector
  //   _attrib.clear();                // Cell properties
  delete[] _Dim;  // dimension system Ax=b
  delete[] _bc[0];
  delete[] _bc[1];      // boundary condition flag
  delete[] _node_dof;   // dof distribution
  delete[] _var_names;  // variables names
  delete[] _refvalue;   // destroy variable unit vec
  delete[] _solver;     // delete solver;
}

// =================================================================
/// This function  is the substructure MGSolBase destructor.
void MGSolBase::clear() {
  // ==============================================================

  for (int Level = 0; Level < _NoLevels; Level++) {
    delete A[Level];
    delete x[Level];  //  A and x  at Level
    delete b[Level];
    delete res[Level];  //  old solutions  at Level
    delete x_old[0][Level];
    delete x_old[1][Level];
    delete x_old[2][Level];
    delete x_nonl[Level];                          //  rhs and residual vector
    delete[] _node_dof[Level];                     // dof distribution at Level
    delete _solver[Level];                         // delete solver  at Level
    if (Level < _NoLevels - 1) delete Rst[Level];  // Restrictor
    if (Level > 0) delete Prl[Level];              // projector
  }
}

// ==================================================
//               MULTILEVEL OPERATOR
//======================================================
//  MGDofBcOp() is defined also in  MGSolverBase
//  all the other are defined in MGSolverDA

// ===============================================================
/// This  function reads all the Operators from files
void MGSolBase::MGDofBcOp() {
// =============================================================
//  initialization of all levels: dofs and matrices;
#ifdef PRINT_INFO  // -------  info --------------------
  std::cout << "\n Reading " << _eqname << " Dof, Bc, Op \n";
#endif  // -------  end      info --------------------
  for (int Level = 0; Level < _NoLevels; Level++) init_dof(Level);  // init the dofmap
  GenBc();                                                          // init the bc-map
  for (int Level = 0; Level < _NoLevels; Level++) init(Level);      // allocation struct
#ifdef PRINT_INFO
  std::cout << " MGDofBcOp(MGSolverBase): Dof, Bc, Op settled \n";
#endif
  return;
}

// ====================================================================
/// This function solves the discrete problem with multigrid solver
void MGSolBase::MGSolve(
    double Eps1,  // tolerance
    int MaxIter,  // n iterations
    // -----------------------------------------------------------
    const int clearing,   ///< 1= clean the init matrix and precond
    const int Gamma,      ///< Control V W cycle
    const int Nc_pre,     ///< n pre-smoothing cycles
    const int Nc_coarse,  ///< n coarse cycles
    const int Nc_post     ///< n post-smoothing cycles
) {
  // ===================================================================
  double rest = 0.;
  // rhs norm ----------------------------------------
  b[_NoLevels - 1]->close();
  double bNorm = b[_NoLevels - 1]->l2_norm();

#ifdef PRINT_CONV
  std::cout << " bNorm " << bNorm << std::endl;
  if (bNorm != bNorm) {
    std::cout << "\033[38;5;196m  --- bNorm is Nan ABORT!! ---- \033[0m \n";
    return;
  }
  if (bNorm > 1.e15) {
    std::cout << "\033[38;5;196m  --- bNorm is too high!!---  \033[0m \n";
    return;
  }
  if (bNorm < 1.e-12) {
    x[_NoLevels - 1]->close();
    x[_NoLevels - 1]->zero();
    std::cout << "\033[38;5;196m  ----------bNorm <1.e-8 !!----------  \033[0m \n";
    return;
  }
#endif
  // FAS Multigrid (Nested) ------------------------------------------------------
  int NestedMG = 1;
  if (NestedMG == 0) {
    x[0]->close();
    x[0]->zero();
    rest = MGStep(0, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post, clearing);
    for (int Level = 1; Level < _NoLevels; Level++) {
      x[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);                                 // projection
      rest = MGStep(Level, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post, clearing);  // MGStepsolution
    }
  }  // -------------------------------------------------------------------------

  // V or W cycle
  int cycle = 0;
  int exit_mg = 0;
  //   MaxIter=10;//ATTENTION
  while (exit_mg == 0 && cycle < MaxIter) {  // multigrid step start -------
#ifndef VANKA
    if (rest > 1.e+20) {
      std::cout << "\033[38;5;196m  ----------residual is too high----------  \033[0m" << endl;
      return;
    }
    rest = MGStep(_NoLevels - 1, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post, clearing);  // MGStep
//       rest=0.; MGCheck(_NoLevels-1); // check projection-restriction
#else
    rest = MGStep_Vanka(
        _NoLevels - 1, 1.e-20, MaxIter, Gamma, Nc_pre, Nc_coarse,
        Nc_post);  // MGStep_Vanka
                   // Vanka_test(0); // check projection-restriction
#endif
    if (rest < Eps1 * (1. + bNorm)) exit_mg = 1;  // exit test
    cycle++;
    // #ifdef PRINT_CONV
    std::cout << " cycle= " << cycle << " residual= " << rest << " \n";
    // #endif
  }  // --------------------- end multigrid step  ------------------
  return;
}

// ====================================================================
/// This function does one multigrid step
double MGSolBase::MGStep(
    int Level,            // Level
    double Eps1,          // Tolerance
    int MaxIter,          // n iterations
    const int Gamma,      // Control V W cycle
    const int Nc_pre,     // n pre-smoothing cycles
    const int Nc_coarse,  // n coarse cycles
    const int Nc_post,    // n post-smoothing cycles
    const int clearing    // clean the init matrix and precond
) {
  // ====================================================================
  std::pair<int, double> rest(0, 0.);
  if (Level == 0) {  // coarse level ----------------------------------
    rest = _solver[Level]->solve(*A[Level], *x[Level], *b[Level], Eps1, 200, clearing);

#ifdef PRINT_CONV
    std::cout << " Coarse res " << rest.second << " " << rest.first << std::endl;
#endif
    // coarse residual
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
  }       // --------------------------------------------------------------
  else {  // fine levels

    // presmoothing (Nu1) ---------------------------------
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t start_time = std::clock();
#endif  //  TC +++++++++++++++
    int Nc_pre1 = Nc_pre;
    if (Level < _NoLevels - 1) Nc_pre1 *= 2;
    rest = _solver[Level]->solve(*A[Level], *x[Level], *b[Level], Eps1, Nc_pre1, clearing);

#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Pre Lev " << Level << " res " << rest.second << " " << rest.first;
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    std::clock_t end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
    // presmoothing residual
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
    // --------- end presmoothing (Nc_pre) ------------------------
    // restriction
    b[Level - 1]->matrix_mult(*res[Level], *Rst[Level - 1]);
    //  solving of system of equations for the residual on the coarser grid
    x[Level - 1]->close();
    x[Level - 1]->zero();
    double coarser_rest;
    for (int g = 1; g <= Gamma; g++)
      coarser_rest = MGStep(Level - 1, Eps1, MaxIter, Gamma, Nc_pre, Nc_coarse, Nc_post, clearing);
    // interpolation of the solution from the coarser grid (projection)
    res[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);
    // adding the coarse solution
    if (coarser_rest < rest.second)
      x[Level]->add(0.95, *res[Level]);  // add the coarser solution only if it helps
      // postsmoothing (Nc_post) --------------------------------------------
#ifdef PRINT_TIME               //  TC +++++++++++++++
    start_time = std::clock();  //   initial set
#endif                          //  TC +++++++++++++++
    //  postsmooting
    int Nc_post1 = Nc_post;
    if (Level < _NoLevels - 1) Nc_post1 *= 2;
    rest = _solver[Level]->solve(*A[Level], *x[Level], *b[Level], Eps1, Nc_post1, clearing);

#ifdef PRINT_CONV  //  CC +++++++++++++++
    std::cout << " Post Lev " << Level << " res " << rest.second << " " << rest.first;
#endif             //  CC +++++++++++++++
#ifdef PRINT_TIME  //  TC +++++++++++++++
    end_time = std::clock();
    std::cout << " time =" << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
#endif  //  TC +++++++++++++++
    //  postsmoothing residual
    res[Level]->resid(*b[Level], *x[Level], *A[Level]);
    // ----------------  end postsmoothing ---------------------------
  }
  // end cycle -------------------------------------
  res[Level]->close();
  return rest.second;
}

// ============================================================================================
/// Check for Prolong and Restr Operators
void MGSolBase::MGCheck(int Level) const {
  std::cout << "\nxlevel-1 before rest\n";
  x[Level - 1]->print();
  std::cout << "\n x level before rest\n";
  x[Level]->print();
  x[Level - 1]->matrix_mult(*x[Level], *Rst[Level - 1]);
  std::cout << "\nxlevel-1 after rest\n";
  x[Level - 1]->print();
  x[Level]->matrix_mult(*x[Level - 1], *Prl[Level]);
  std::cout << "\n x level after prol\n";
  x[Level]->print();
  //   x[Level-1]->matrix_mult(*x[Level],*Rst[Level-1]);
  //   x[Level]  ->matrix_mult(*x[Level-1],*Prl[Level]);
  return;
}

// *************************************************************************
//  SET FUNCTIONS
// ***********************************************************************=

#ifdef TWO_PHASE_LIB
void MGSolBase::set_mgcc(MGSolCC& cc) { _msolcc = &cc; }
#endif

// =========================================
/// This function sets up data structures for each problem class
// =========================================
void MGSolBase::setUpExtFieldData() {
  /// A) set up _mg_eqs
  // external system and index vectors
  for (int deg = 0; deg < 3; deg++) {
    for (int kl = 0; kl < _data_eq[deg].max_neqs; kl++) {
      _data_eq[deg].n_eqs = 0;
      _data_eq[deg].mg_eqs[kl] = NULL;
      _data_eq[deg].indx_ub[kl] = -1;
      _data_eq[deg].tab_eqs[kl] = -1;

      for (int kk = 0; kk < NDOF_FEM; ++kk) {
        _data_eq[deg].ub[kk + kl * NDOF_FEM] = 0.;  // data
      }
    }
  }
  // start index K from 0, L from 0, Q from DIMENSION (coordinates+ q variable)
  _data_eq[0].indx_ub[0] = 0;  //_data_eq[0].n_eqs=0; // piecewice constant  (0)
  _data_eq[1].indx_ub[0] = 0;  //_data_eq[1].n_eqs=0; // piecewice linear    (1)
  _data_eq[2].indx_ub[0] = 0;  //_data_eq[2].n_eqs=0; // piecewice quadratic (2)
  return;
}
// *****************************************************************************
//  SOLVER
// *****************************************************************************

/// ======================================================
/// This function controls the time step operations:
/// ======================================================
void MGSolBase::MGTimeStep(const double time, const int) {
  std::cout << std::endl << "  " << _eqname.c_str() << " solution " << std::endl;

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if PRINT_TIME == 1
  std::clock_t start_time = std::clock();
#endif
  GenMatRhs(time, _NoLevels - 1, 1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for (int Level = 0; Level < _NoLevels - 1; Level++) GenMatRhs(time, Level, 0);

#if PRINT_TIME == 1
  std::clock_t end_time = std::clock();
  std::cout << " Assembly time =" << double(end_time - start_time) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve(1.e-6, 40);
#if PRINT_TIME == 1
  std::clock_t end_timef = std::clock();
  std::cout << " Assembly+solution time =" << double(end_timef - start_time) / CLOCKS_PER_SEC << "s "
            << std::endl;
#endif
  /// [d] Update of the old solution at the top Level
  x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);
  return;
}

/// ======================================================
/// This function controls the time step operations with no update for old solution:
/// ======================================================
void MGSolBase::MGTimeStep_no_up(const double time, const int) {
  std::cout << std::endl << "  " << _eqname.c_str() << " solution " << std::endl;

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if PRINT_TIME == 1
  std::clock_t start_time = std::clock();
#endif
  GenMatRhs(time, _NoLevels - 1, 1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for (int Level = 0; Level < _NoLevels - 1; Level++) GenMatRhs(time, Level, 0);

#if PRINT_TIME == 1
  std::clock_t end_time = std::clock();
  std::cout << " Assembly time =" << double(end_time - start_time) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve(1.e-6, 40);
#if PRINT_TIME == 1
  std::clock_t end_timef = std::clock();
  std::cout << " Assembly+solution time =" << double(end_timef - start_time) / CLOCKS_PER_SEC << "s "
            << std::endl;
#endif

  return;
}
/// ======================================================
/// This function  updates the solution x_old from x
/// ======================================================
void MGSolBase::MGUpdateStep() {
  std::cout << std::endl << "  " << _eqname.c_str() << " update solution " << std::endl;
  x[_NoLevels - 1]->localize(*x_old[0][_NoLevels - 1]);
  return;
}
/// ======================================================
/// This function displace back the mesh
/// ======================================================
void MGSolBase::MGUndo_disp() {
  std::cout << "Wrong use of MGUndo_disp from MGSolBase.C, aborting";
  abort();
  return;
}

// *****************************************************************************
//  COMPUTE fUNCTIONS
// *****************************************************************************

/// ======================================================
/// This function  computes the  functional
/// ======================================================
double MGSolBase::MGFunctional(
    double /*parameter*/,  /// Use of the function: (0) compute functional OR (1) set _eta
    double& /*control*/    /// \param[in] <>  eta multiplier for optimal method
) {
  std::cout << "Wrong use of MGFunctional from MGSolBase.C, aborting";
  abort();
  return 1;
}
/// =======================================================================
/// This function sets the controlled domain for optimal control problems
/// =======================================================================
void MGSolBase::set_ctrl_dom(
    const double xMin, const double xMax, const double yMin, const double yMax, const double zMin,
    const double zMax) {
  int el_conn[NDOF_FEM];
  double x_m[DIMENSION];
  double xx_qnds[NDOF_FEM * DIMENSION];
  const int n_elem = _mgmesh._off_el[0][_NoLevels + _NoLevels * (_mgmesh._n_subdom - 1)];
  double eps = 1.e-6;                                                            // tolerance for coordinates
  const int nel_e = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc + 1];  // start element
  const int nel_b = _mgmesh._off_el[0][_NoLevels - 1 + _NoLevels * _iproc];      // stop element

  _weight_ctrl.resize(n_elem);
  for (int iel = 0; iel < (nel_e - nel_b); iel++) {
    _mgmesh.get_el_nod_conn(0, _NoLevels - 1, iel, el_conn, xx_qnds);  // gets element coordinates
    for (int idim = 0; idim < DIMENSION; idim++) {
      x_m[idim] = 0;
      for (int d = 0; d < NDOF_FEM; d++) {
        x_m[idim] += xx_qnds[idim * NDOF_FEM + d] / NDOF_FEM;
      }                             // end d loop
    }                               // end idim loop
    _weight_ctrl[iel + nel_b] = 0;  // default is 0
    if (x_m[0] > xMin - eps && x_m[0] < xMax + eps && x_m[1] > yMin - eps && x_m[1] < yMax + eps
#if DIMENSION == 3
        && x_m[2] > zMin - eps && x_m[2] < zMax + eps
#endif
    ) {
      _weight_ctrl[iel + nel_b] = 1;  // 1 inside control region
      //  std::printf("iel %4d wieght %4f xm0 %4.5f xm1 %4.5f xm2 %4.5f\n",iel+nel_b,
      // _weight_ctrl[iel+nel_b],x_m[0],x_m[1],x_m[2]);
    }
  }  // end iel loop
  return;
}

/// ======================================================
/// This function change current dt
/// ======================================================

void MGSolBase::set_dt(double /*dt*/) {}

// ===================================================================
// ================== Print/Read function (All virtual)===============
// ===================================================================

// =============================================
/// Print xml attrib
void MGSolBase::print_u_xdmf(
    std::ofstream& out,  //  file xdmf
    int nodes, int /*nelems*/,
    std::string file_name) const {  // ================================

  for (int ivar = 0; ivar < _n_vars; ivar++) {
    std::string var_name = _var_names[ivar];
    out << "<Attribute Name=\"" << var_name << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << nodes << "  " << 1
        << "\" Format=\"HDF\">  \n";
    out << file_name
        //       femus_dir << "/" << output_dir << basesol << "."
        //       << setw(ndigits) << setfill('0') << t_step << ".h5"
        << ":" << var_name << "\n";
    out << "</DataItem>\n"
        << "</Attribute>";
  }
  return;
}
