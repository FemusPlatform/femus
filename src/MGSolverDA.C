// #include "Equations_tab.h"
// #include "Equations_conf.h"
// ===============================
// ===============================
// std lib
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

// configuration files ----------
// #include "Printinfo_conf.h"

// class files ----------------
#include "MGSolverDA.h"
// #include "MGSclass_conf.h"

// local includes --------------
#include "EquationSystemsExtendedM.h"
#include "MGFE_conf.h"
#include "MGGeomEl.h"
#include "MGSystem.h"
#include "MeshExtended.h"
// #include "dense_set.h"

#include "MGEquationsSystem.h"
#include "MGFE.h"
#include "MGFEMap.h"
#include "MGGraph.h"
#include "MGUtils.h"

// #ifdef HAVE_MED
// #include "InterfaceFunctionM.h"
// #include "MEDCouplingFieldDouble.hxx"
// #include "MEDCouplingUMesh.hxx"
// #include "MEDLoader.hxx"
// #endif

// algebric includes -----------
#include "dense_matrixM.h"
#include "dense_vectorM.h"
#include "numeric_vectorM.h"
#include "sparse_MmatrixM.h"
#include "sparse_matrixM.h"

// ==========================================================
//               MGSolDA functions
// ==========================================================

/// This function is the MGSolDA constructor I
MGSolDA::MGSolDA(
    MGEquationsSystem& mg_equations_map_in,  //
    const int nvars_in[],                    // # of quad variables
    //   const int nvarsl_in,                // # of linear variables
    std::string eqname_in,       // equation name
    std::string /*varname_in*/)  // basic variable name
    : MGSolBase(
          mg_equations_map_in, nvars_in, eqname_in) {  // ===============================================
  // piecewise constant linear, quadratic and (...cubic)
  _fe[0] = _mgfemap.get_FE(0);  ///> Lagrange piecewise constant
  _fe[1] = _mgfemap.get_FE(1);  ///> Lagrange piecewise linear
  _fe[2] = _mgfemap.get_FE(2);  ///> Lagrange piecewise quadratic

  _el_dof[0] = ((nvars_in[0] > 0) ? NDOF_K : 0);    // Lagrange piecewise constant variables
  _el_dof[1] = ((nvars_in[1] > 0) ? NDOF_P : 0);    // Lagrange piecewise linear variables
  _el_dof[2] = ((nvars_in[2] > 0) ? NDOF_FEM : 0);  // Lagrange piecewise linear variables

  // System names and units
  _var_names = new std::string[_n_vars];  // names
  _refvalue = new double[_n_vars];

  // output string stream
  int iname = 0;
  for (; iname < nvars_in[2]; iname++) {  // 2 quadratic
    std::ostringstream ostr;
    ostr << "qua" << iname + 1;  // use the string stream just like cout,
    _var_names[iname] = ostr.str();
    _refvalue[iname] = 1;
  }

  for (; iname < (nvars_in[1] + nvars_in[2]); iname++) {  // 1 linear
    std::ostringstream ostr;
    ostr << "lin" << iname - nvars_in[2] + 1;  // use the string stream just like cout,
    _var_names[iname] = ostr.str();
    _refvalue[iname] = 1;
  }

  for (; iname < (nvars_in[0] + nvars_in[1] + nvars_in[2]); iname++) {  // 0 piecewise
    std::ostringstream ostr;
    ostr << "pie" << iname - nvars_in[2] - nvars_in[1] + 1;  // use the string stream just like cout,
    _var_names[iname] = ostr.str();
    _refvalue[iname] = 1;
  }

  _NumRestartSol = 1;
  return;
}

// =========================================
/// This function is the destructor
// =========================================
MGSolDA::~MGSolDA() {
  //     delete []_data_eq;
}

// =========================================
/// This function is the main initialization function:
// =========================================
void MGSolDA::init(const int Level) {
  // name directories to read
  std::string f_matrix = _mgutils.get_file("F_MATRIX");
  std::string f_rest = _mgutils.get_file("F_REST");
  std::string f_prol = _mgutils.get_file("F_PROL");
  const ParallelM::Communicator& comm1 = _mgmesh._comm.comm();
  // number of partial dofs ---------------------------------------------------------
  int off_proc = _NoLevels * _iproc;
  int n_local_q = (_el_dof[2] > 0) ? _mgmesh._off_nd[0][Level + 1 + off_proc] - _mgmesh._off_nd[0][off_proc]
                                   : 0;  // 0=QUADRATIC
  int n_local_l = (_el_dof[1] > 0) ? _mgmesh._off_nd[1][Level + 1 + off_proc] - _mgmesh._off_nd[1][off_proc]
                                   : 0;  // 1=LINEAR
  int n_local_k =
      (_el_dof[0] > 0)
          ? (_mgmesh._off_el[0][Level + 1 + off_proc] - _mgmesh._off_el[0][Level + off_proc]) * _el_dof[0]
          : 0;  // 0=Volume
  //   int nel_local=_mgmesh._off_el[0][Level+off_proc+1]-_mgmesh._off_el[0][Level+off_proc];
  // global dofs  ---------------------------------------------------------------------
  int n_glob = _Dim[Level];
  int n_local = _nvars[1] * n_local_l + _nvars[2] * n_local_q + _nvars[0] * n_local_k;
  //  int  nel_glob=_mgmesh._NoElements[0][Level];

  //   matrix ------------------------------------------------------------------------------
  std::ostringstream filename("");
  filename << _mgutils._inout_dir << f_matrix << Level << ".h5";
  A[Level] = SparseMatrixM::build(_mgmesh._comm.comm()).release();
  A[Level]->init(n_glob, n_glob, n_local, n_local);
  ReadMatrix(Level, filename.str(), *A[Level], _nvars);

  // vectors ------------------------------------------------------------------------------
  b[Level] = NumericVectorM::build(comm1).release();
  b[Level]->init(n_glob, n_local, false, AUTOMATICM);
  res[Level] = NumericVectorM::build(comm1).release();
  res[Level]->init(n_glob, n_local, false, AUTOMATICM);
  x[Level] = NumericVectorM::build(comm1).release();
  x[Level]->init(n_glob, n_local, false, AUTOMATICM);
  for (int ieq = 0; ieq < 3; ieq++) {
    x_old[ieq][Level] = NumericVectorM::build(comm1).release();
    x_old[ieq][Level]->init(n_glob, false, SERIALM);
  }
  x_nonl[Level] = NumericVectorM::build(comm1).release();
  x_nonl[Level]->init(n_glob, false, SERIALM);

  // Restrictor -------------------------------------------------------------------------
  if (Level < _NoLevels - 1) {
    Rst[Level] = SparseMMatrixM::build(comm1).release();
    Rst[Level]->init(n_glob, _Dim[Level + 1], n_glob, _Dim[Level + 1]);
    //     int n_nodes_qp1= _mgmesh._NoNodes[Level+1];
    filename.str("");
    filename << _mgutils._inout_dir << f_rest << Level + 1 << "_" << Level << ".h5";
    ReadRest(
        Level, filename.str(), *Rst[Level], _nvars, _node_dof[Level + 1], _node_dof[Level],
        _node_dof[_NoLevels - 1]);
    Rst[Level]->close();
  }
  // Prolongation ----------------------------------------------------------------------
  if (Level > 0) {
    Prl[Level] = SparseMMatrixM::build(comm1).release();
    Prl[Level]->init(n_glob, _Dim[Level - 1], _Dim[Level], _Dim[Level - 1]);
    filename.str("");
    filename << _mgutils._inout_dir << f_prol << Level - 1 << "_" << Level << ".h5";
    ReadProl(Level, filename.str(), *Prl[Level], _nvars, _node_dof[Level], _node_dof[Level - 1]);
    Prl[Level]->close();
    //  Prl[Level]->print_personal(filename);//print on screen Prol
  }

  return;
}

// ===========================================================
/// This function initializes the system degrees of freedom (dof)
void MGSolDA::init_dof(
    const int Level,    //< Level
    const int vb_0,     //< type of element
    const int /*n_vb*/  //< number of type element
) {                     // ========================================================

  // Set up from mesh -----------------------------
  const int n_subdom = _mgmesh._n_subdom;
  const int n_nodes = _mgmesh._NoNodes[Level];
  const int n_elem = _mgmesh._NoElements[vb_0][Level];  // 0=volume

  const int offset = _mgmesh._NoNodes[_NoLevels - 1];
  const int* off_nd_q = _mgmesh._off_nd[0];
  const int* off_nd_l = _mgmesh._off_nd[1];
  const int* off_el = _mgmesh._off_el[0];  // volume elements

  // number of total dofs
  int n_nodes_l = _mgmesh._NoNodes[_mgmesh._NoFamFEM * _NoLevels];

  if (Level > 0) { n_nodes_l = _mgmesh._NoNodes[Level - 1]; }

  _Dim[Level] = _nvars[0] * n_elem * _el_dof[0] + _nvars[1] * n_nodes_l + _nvars[2] * n_nodes;

#ifdef PRINT_INFO
  std::cout << Level << " node " << n_nodes << " dof " << _Dim[Level] << " pres " << n_nodes_l << std::endl;
#endif

  // Set up boundary conditions ++++++++++++++++++++++++++++
  if (Level == _NoLevels - 1) {
    _bc[0] = new int[_Dim[Level]];
    _bc[1] = new int[_Dim[Level]];

    for (int k1 = 0; k1 < _Dim[Level]; k1++) {
      _bc[0][k1] = 1;
      _bc[1][k1] = 1;
    }
  }

  // construction dof node vector(_node_dof) +++++++++++++++++++++
  _node_dof[Level] = new int[_n_vars * offset];

  for (int k1 = 0; k1 < offset * _n_vars; k1++) { _node_dof[Level][k1] = -1; }

  int count = 0;
  int ndof_lev = 0;

  for (int isubdom = 0; isubdom < n_subdom; isubdom++) {
    // quadratic -----------------------------------
    int off_proc = isubdom * _NoLevels;
    for (int ivar = 0; ivar < _nvars[2]; ivar++) {
      for (int k1 = off_nd_q[off_proc]; k1 < off_nd_q[off_proc + Level + 1]; k1++) {
        _node_dof[Level][k1 + ivar * offset] = count;
        count++;
      }
    }

    // linear -----------------------------------
    for (int ivar = 0; ivar < _nvars[1]; ivar++) {
      for (int k1 = off_nd_q[off_proc];
           k1 < off_nd_q[off_proc] + off_nd_l[Level + 1 + off_proc] - off_nd_l[off_proc]; k1++) {
        _node_dof[Level][k1 + (_nvars[2] + ivar) * offset] = count;
        count++;
      }
    }

    // constant polynomial of order _el_dof[0] -----------------------------------
    int delta_el = off_el[off_proc + Level + 1] - off_el[off_proc + Level];
    for (int ivar = 0; ivar < _nvars[0]; ivar++) {
      for (int iel = 0; iel < delta_el; iel++) {
        for (int idof = 0; idof < _el_dof[0]; idof++) {
          _node_dof[Level][idof + (iel + ndof_lev) * _el_dof[0] + (_nvars[2] + _nvars[1] + ivar) * offset] =
              count;
          count++;
        }
      }
    }
    ndof_lev += delta_el;
  }

#ifdef PRINT_INFO
  std::cout << "MGSol::init_dof(D)   Level= " << Level << std::endl;
#endif
  return;
}

//  ******************************************************************************************
//  COMPUTE FUNCTION
//  ******************************************************************************************
// =====================================================================
//-----------------------------------------------------------------------------------
double MGSolDA::CalcFUpwind(
    double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes) {
  double vel_modulus = 1.e-10, h_eff = 1.e-20, f_upwind = 0.;

  for (int i = 0; i < Dim; i++) vel_modulus += VelOnGauss[i] * VelOnGauss[i];
  vel_modulus = sqrt(vel_modulus);

  for (int i = 0; i < NbOfNodes; i++) {
    double hh = 1.e-20;
    for (int idim = 0; idim < Dim; idim++)
      hh += VelOnGauss[idim] * PhiDer[i + idim * NbOfNodes] / vel_modulus;
    h_eff += fabs(hh);
  }
  h_eff = 2. / h_eff;
  if (h_eff < 1.e-10) {
    h_eff = 1.;
    std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";
  }
  // STANDARD SUPG
  const double Pe_h = 0.5 * vel_modulus * h_eff / (Diffusivity);
  const double a_opt = (1. / tanh(Pe_h) - 1. / Pe_h);
  if (a_opt > 1.) { std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n"; }
  //   f_upwind = 0.5*a_opt*h_eff/ (vel_modulus);
  f_upwind = 0.5 * a_opt * h_eff / (vel_modulus);
  return f_upwind;
}
// =====================================================================
double MGSolDA::MGFunctional(double, double&) {
  std::cout << "Not implemented in MsolverDA";
  return 0.;
}

// =====================================================================
void MGSolDA::compute_jac(
    const int j, const int idim,
    double uold_b[],      // node values <-
    const int nvars,      // # of variables  <-
    const double phi[],   // shape functions  <-
    const double dphi[],  // derivatives of the shape functions  <-
    const int n_shape,    // # of shape functions  <-
    double u_forw[],      // interpolated function ->
    double u_back[],      // interpolated function ->
    double u_forw_dx[],   // interpolated derivatives ->
    double u_back_dx[]    // interpolated derivatives ->
    ) const {             // =================================================================
  // All data vectors are NDOF_FEM long
  const double alfa = 1.e-08;
  // variable loop
  for (int ivar = 0; ivar < nvars; ivar++) {
    u_forw[ivar] = 0.;
    u_back[ivar] = 0.;  // set zero
    for (int jdim = 0; jdim < DIMENSION; jdim++) {
      u_forw_dx[ivar * DIMENSION + jdim] = 0.;  // set zero
      u_back_dx[ivar * DIMENSION + jdim] = 0.;  // set zero
    }
    // interpolation with shape functions
    for (int eln = 0; eln < n_shape; eln++) {
      const int indx = eln + ivar * NDOF_FEM;
      if (indx == j + idim * NDOF_FEM) {
        u_forw[ivar] += phi[eln] * (uold_b[indx] + alfa);
        u_back[ivar] += phi[eln] * (uold_b[indx] - alfa);
      } else {
        u_forw[ivar] += phi[eln] * uold_b[indx];
        u_back[ivar] += phi[eln] * uold_b[indx];
      }

      for (int jdim = 0; jdim < DIMENSION; jdim++) {
        if (indx == j + idim * NDOF_FEM) {
          u_forw_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * (uold_b[indx] + alfa);
          u_back_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * (uold_b[indx] - alfa);
        } else {
          u_forw_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
          u_back_dx[ivar * DIMENSION + jdim] += dphi[eln + jdim * n_shape] * uold_b[indx];
        }
      }
    }
  }

  return;
}

// =========================================
void MGSolDA::set_xooold2x() {
  for (int Level = 1; Level <= _NoLevels; Level++) {
    /// A. Setup
    const int offset = _mgmesh._NoNodes[Level - 1];  // fine level # of nodes
    //     const int pie_offset = (NSUBDOM) * (_mgmesh._NoElements[0][Level - 1]);  // fine level # of nodes

    // field
    // file to read
    // reading loop over system varables
    for (int ivar = 0; ivar < _nvars[2] + _nvars[1]; ivar++) {
      int el_nds = NDOF_FEM;

      if (ivar >= _nvars[2]) {
        el_nds = NDOF_P;  // quad and linear
      }

      // reading ivar param
      //       double Irefval = 1. / _refvalue[ivar];  // units

      // storing  ivar variables (in parallell)
      for (int iel = 0; iel < _mgmesh._off_el[0][_iproc * _NoLevels + Level] -
                                  _mgmesh._off_el[0][_iproc * _NoLevels + Level - 1];
           iel++) {
        int elem_gidx = (iel + _mgmesh._off_el[0][_iproc * _NoLevels + Level - 1]) * NDOF_FEM;

        for (int i = 0; i < el_nds; i++) {            // linear and quad
          int k = _mgmesh._el_map[0][elem_gidx + i];  // the global node
          const double value = (*x_old[2][Level - 1])(_node_dof[_NoLevels - 1][k + ivar * offset]);
          x[Level - 1]->set(_node_dof[_NoLevels - 1][k + ivar * offset], value);  // set the field
        }
      }
    }

    int ndof_lev = 0;

    for (int pr = 0; pr < _mgmesh._iproc; pr++) {
      int delta =
          _mgmesh._off_el[0][pr * _NoLevels + _NoLevels] - _mgmesh._off_el[0][pr * _NoLevels + _NoLevels - 1];
      ndof_lev += delta;
    }

    /// D. delocalization and clean
    x[Level - 1]->localize(*x_old[0][Level - 1]);
    x[Level - 1]->localize(*x_old[1][Level - 1]);
  }

  return;
}

/// ==================================================
///   ==============SOLVERS  SOURCE================
/// ==================================================

/// ======================================================
/// This function controls the time step operations:
/// ======================================================
void MGSolDA::MGTimeStep(const double time, const int) {
  std::cout << std::endl << "  " << _eqname.c_str() << " solution " << std::endl;
  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if PRINT_TIME == 1
  std::clock_t start_time = std::clock();
#endif
  GenMatRhs(time, _NoLevels - 1, 1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for (int Level = 0; Level < _NoLevels - 1; Level++) { GenMatRhs(time, Level, 0); }

#if PRINT_TIME == 1
  std::clock_t end_time = std::clock();
  std::cout << " Assembly time =" << double(end_time - start_time) / CLOCKS_PER_SEC << " s " << std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve(1.e-8, 40);
#if PRINT_TIME == 1
  std::clock_t end_timef = std::clock();
  std::cout << " Assembly+solution time =" << double(end_timef - start_time) / CLOCKS_PER_SEC << "s "
            << std::endl;
#endif
  /// [d] Update of the old solution at the top Level
  //   x[_NoLevels-1]->localize(*x_old[0][_NoLevels-1]);
  return;
}

//  =========================================================================================
// EQUATION ACTIVATION (START THE CLASS POINTER DURING EXECUTION
//  =========================================================================================
// int Order,      = order: const,linear,quad
// int Field,      = Field to activate
// std::string sfn,= system field name
// int& n_index,  = n_index (collecting index)
// int vector,    = vector
// int neqs       = dimension (number of eqs) neqs
// int coupled    = coupled (1) or segregated (0)  solver
// ==========================================================================================

void MGSolDA::ActivateVectField(
    int Order,        ///< order: const,linear,quad
    int Field,        ///< Field to activate
    std::string sfn,  ///< SystemFieldName
    int& n_index,     ///< n_index (collecting index)
    int coupled       ///< coupled (1) or segregated (0)  solver
) {
  std::string FieldX = sfn + "X";
  std::string FieldY = sfn + "Y";
  std::string FieldZ = sfn + "Z";

  if (coupled == 0) {  // flag 0 in SimulationConfiguration -> UNcoupled -------------------
    if (sfn.compare("TAU") == 0) {
      std::string FieldX = sfn + "XX";
      std::string FieldY = sfn + "XY";
      std::string FieldZ = sfn + "YY";
      ActivateEquation(Order, Field, FieldX, n_index);
      ActivateEquation(Order, Field + 1, FieldY, n_index);
      ActivateEquation(Order, Field + 2, FieldZ, n_index);
    } else {  //  if (sfn.compare("TAU") != 0)
      ActivateEquation(Order, Field, FieldX, n_index);
      ActivateEquation(Order, Field + 1, FieldY, n_index);
#if DIMENSION == 3
      ActivateEquation(Order, Field + 2, FieldZ, n_index);
#endif
    }
  }  // -------------------------------------------------------------------------------------
  else {
    // flag 1 in SimulationConfiguration -> coupled ----------------------------------------
    _data_eq[Order].tab_eqs[Field] = n_index;                  // table
    _data_eq[Order].mg_eqs[n_index] = _mgeqnmap.get_eqs(sfn);  // FSI equation pointer
    _data_eq[Order].indx_ub[n_index + 1] =
        _data_eq[Order].indx_ub[n_index] + DIMENSION;  // _data_eq[2].ub index
    _data_eq[Order].n_eqs++;                           // number of quadratic system
    n_index++;                                         // update counter
  }  // -------------------------------------------------------------------------------------

  return;
}

// ================================================================================================
void MGSolDA::ActivateControl(
    int Order,        ///< order: const,linear,quad
    int Field,        ///< Field to activate             // Field to activate
    std::string sfn,  ///< SystemFieldName // system field name
    int& n_index,     ///< n_index (collecting index)             // n_index (collecting index)
    int vector,       ///< coupled (1) or segregated (0)  solver
    int neqs          ///< dimension (number of eqs)
) {  //===============================================================================================
  // // flag 2 in SimulationConfiguration.in ->  split
  // // flag 1 in SimulationConfiguration -> coupled
  // -------------------------------------------------------------------------------------------------
  std::string FieldX = sfn + "X";
  std::string FieldY = sfn + "Y";
  std::string FieldZ = sfn + "Z";
  ActivateEquation(Order, Field, FieldX, n_index);
  if (vector == 1) {
    ActivateEquation(Order, Field + 1, FieldY, n_index);
    if (neqs == 3) ActivateEquation(Order, Field + 2, FieldZ, n_index);
  }

  return;
}
//--------------------------------------------------------------------------------------------------
// This function activates a scalar equation
void MGSolDA::ActivateScalar(
    int Order,        ///< order: const,linear,quad
    int Field,        ///< Field to activate
    std::string sfn,  ///< System Field Name
    int& n_index      ///< n_index (collecting index)
) {
  ActivateEquation(Order, Field, sfn, n_index);
  return;
}
//--------------------------------------------------------------------------------------------------
// This function activates coupled equations
void MGSolDA::ActivateCoupled(
    int Order,        ///< order: const,linear,quad
    int Field,        ///< Field 1 to activate
    std::string sfn,  ///< System Field Name
    int& n_index,     ///< n_index (collecting index)
    std::string sfn2  ///< Field 2 to activate
) {
  ActivateEquation(Order, Field, sfn, n_index);
  ActivateEquation(Order, Field + 1, sfn2, n_index);
  return;
}
//--------------------------------------------------------------------------------------------------
// This function activates a single  equation
void MGSolDA::ActivateEquation(
    int Order,        ///< order: const,linear,quad
    int Field,        ///< Field to activate
    std::string sfn,  ///< System Field Name
    int& n_index      ///< n_index (collecting index)
) {
  _data_eq[Order].tab_eqs[Field] = n_index;                  // table
  _data_eq[Order].mg_eqs[n_index] = _mgeqnmap.get_eqs(sfn);  // Navier-Stokes equation pointer
  _data_eq[Order].indx_ub[n_index + 1] = _data_eq[Order].indx_ub[n_index] + 1;  // _data_eq[2].ub index
  _data_eq[Order].n_eqs++;                                                      // number of quadratic system
  n_index++;                                                                    // update counter
  return;
}
