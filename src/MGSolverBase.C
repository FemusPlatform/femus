// lib include ------------------
#include <sstream>
#include <limits>

// configure file -------------
#include "Solverlib_conf.h"
#include "Printinfo_conf.h"
#include "MGFE_conf.h"
// class include
#include "MGSolverBase.h"
#ifdef   TWO_PHASE
#include "MGSolverCC.h"
#endif
// local inlude -----------------
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGMesh.h"

#include "MGEquationsSystem.h"
#include "MGFEMap.h"

// alg include-------------------
#include "sparse_matrixM.h"
#include "sparse_MmatrixM.h"
#include "numeric_vectorM.h"
#include "linear_solverM.h"


#include "petsc_macroM.h"
#include <mpi.h>  //for MPI_COMM_WORLD

//  #define VANKA (0)
//
// ========================================================
/// This function  is the MGSolBase constructor :
MGSolBase::MGSolBase (MGEquationsSystem &e_map_in, // equation map
                      const int nvars_in[],            // # variables
                      std::string eqname_in     // equation name
                     ) :
  _mgutils (e_map_in._mgutils),
//     _mgphys(e_map_in._mgphys),
  _mgfemap (e_map_in._mgfemap),
  _mgmesh (e_map_in._mgmesh),
  _mgeqnmap (e_map_in),
  _NoLevels ( (int) (e_map_in._mgutils._geometry["nolevels"])), //you can do that
  _eqname (eqname_in),
  _n_vars (nvars_in[0]+nvars_in[1]+nvars_in[2]), _var_names (NULL), _refvalue (NULL)
#ifdef  TWO_PHASE_LIB
  , _msolcc (NULL)
#endif
{

  _nvars[0]= nvars_in[0]; //Lagrange piecewise constant variables
  _nvars[1]= nvars_in[1]; //Lagrange piecewise linear variables
  _nvars[2]= nvars_in[2]; //Lagrange piecewise linear variables
  // processor number ---------------------------
#ifndef HAVE_LASPACKM
  _iproc=_mgmesh._iproc;         // > 1CPU
#else
  _iproc=0;                      // 1 CPU
#endif

  //allocation of dynamic system ---------------
  _Dim =new int[_NoLevels];                          // matrix and vect  dim
  A.resize (_NoLevels);
  x.resize (_NoLevels);    // matrix vect sol
  x_old.resize (_NoLevels);
  x_oold.resize (_NoLevels); // old solution
  x_nonl.resize (_NoLevels); // non linear solution
  disp.resize (_NoLevels); // displacement for mesh
  disp_old.resize (_NoLevels); // displacement for mesh
  disp_oold.resize (_NoLevels); // displacement for mesh
  x_ooold.resize (_NoLevels); // vector for multiple uses
  b.resize (_NoLevels);
  res.resize (_NoLevels);  // rhs

  // restr and prol operators -----------------
  Rst.resize (_NoLevels);
  Prl.resize (_NoLevels);   // Projector and restrictor

  // dof info ---------------
  _node_dof= new int*[_NoLevels+1];  // dof (+1)
//   bc= new int*[_NoLevels];           // bc
//   _attrib.resize(_NoLevels);         // Attribute vector (cell)

  //solver ---------------------------------------
  _solver=new LinearSolverM*[_NoLevels];     // one solver for each level
  // get the communicator from the PETSc object
//   ParallelM::communicator comm;
//   PetscObjectGetComm((PetscObject)pc, &comm);
//   const ParallelM::Communicator communicator(comm);

  for (int l=0; l<_NoLevels; l++) _solver[l]=
      LinearSolverM::build (_mgmesh._comm.comm(),LSOLVER).release();
  _control=0.;
#ifdef   TWO_PHASE

#endif


  return;
}



// ===================================================
/// This function  is the MGSolBase destructor.
MGSolBase::~MGSolBase (
) {
// ===================================================
  // clear substructrures
   clear(); 
  A.clear();
  x.clear();        //  A and x
  x_old.clear();
  x_oold.clear();  //  old solutions
  x_nonl.clear(); // nonlinear solution tmp
  disp.clear(); // displacement for mesh
  disp_old.clear(); // displacement for mesh
  disp_oold.clear(); // displacement for mesh
  b.clear();
  res.clear();      //  rhs and residual vector
  Rst.clear();
  Prl.clear();    // Restrictor and projector
//   _attrib.clear();                // Cell properties
  delete [] _Dim;                 // dimension system Ax=b
  delete[] bc[0];
  delete[] bc[1];                    // boundary condition flag
  delete [] _node_dof;            // dof distribution
  delete [] _var_names;           // variables names
  delete []_refvalue;             // destroy variable unit vec
  delete []_solver;               //delete solver;
}

// =================================================================
/// This function  is the substructure MGSolBase destructor.
void MGSolBase::clear (
) {
// ==============================================================

  for (int Level =0; Level<_NoLevels; Level++) {
    delete A[Level];
    delete x[Level];             //  A and x  at Level
    delete b[Level];
    delete res[Level];             //  old solutions  at Level
    delete x_old[Level];
    delete x_oold[Level];
    delete x_ooold[Level];
    delete x_nonl[Level];    //  rhs and residual vector
    delete disp[Level];
    delete disp_old[Level];
    delete disp_oold[Level];
//     delete _attrib[Level];                         // Cell properties  at Level
    delete [] _node_dof[ Level];                   // dof distribution at Level
    delete _solver[Level];                         //delete solver  at Level
    if (Level < _NoLevels - 1) delete Rst[Level];  // Restrictor
    if (Level > 0) delete Prl[Level];              // projector
  }
}

#ifdef    TWO_PHASE_LIB
void MGSolBase::set_mgcc (MGSolCC &cc) {
  _msolcc=&cc;
}
#endif
// ===============================================================
/// This  function reads all the Operators from files
void MGSolBase::MGDofBcOp (
) {
// =============================================================
  //  initialization of all levels: dofs and matrices;
#ifdef PRINT_INFO  // -------  info --------------------
  std::cout << "\n Reading "  <<  _eqname << " Dof, Bc, Op \n";
#endif   // -------  end      info --------------------
  for (int Level = 0; Level< _NoLevels; Level++)  init_dof (Level); //init the dofmap

  GenBc();   //init the bc-map

  for (int Level = 0; Level< _NoLevels; Level++)  {
    init (Level); // allocation struct

  }

#ifdef PRINT_INFO
  std::cout << " MGDofBcOp(B): Dof, Bc, Op settled \n";
#endif
  return;
}

// ====================================================================
/// This function solves the discrete problem with multigrid solver
void MGSolBase::MGSolve (
  double Eps1,         // tolerance
  int MaxIter,         // n iterations
  // -----------------------------------------------------------
  const int clearing,  // 1= clean the init matrix and precond 
  const int Gamma,     // Control V W cycle
  const int Nc_pre,    // n pre-smoothing cycles
  const int Nc_coarse, // n coarse cycles
  const int Nc_post    // n post-smoothing cycles                      
                        ) {
// ===================================================================
  double rest=0.;
  // rhs norm ----------------------------------------
  b[_NoLevels-1]->close();
  double bNorm=b[_NoLevels-1]->l2_norm();

#ifdef PRINT_CONV
  std::cout << " bNorm " << bNorm  << std::endl;
    if (bNorm!=bNorm) {
     std::cout << "\033[38;5;196m  ----------bNorm is Nan ABORT!!----------  \033[0m" <<  endl;
    return;}
    if (bNorm>1.e15) {
     std::cout << "\033[38;5;196m  ----------bNorm is too high!!----------  \033[0m" <<  endl;
    return;} 
#endif
  // FAS Multigrid (Nested) ---------------------------
  int NestedMG=1;
  if (NestedMG==0) {
    x[0]->zero();
    MGStep (0,1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);
    for (int Level = 1; Level < _NoLevels; Level++) {
      x[Level]->matrix_mult (*x[Level-1],*Prl[Level]); // projection
      rest= MGStep (Level,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);             // MGStepsolution
    }
  } // ------------------------------------------------

  // V or W cycle
  int cycle=0;
  int exit_mg=0;
//   MaxIter=10;//ATTENTION
  while (exit_mg==0 && cycle<MaxIter) { // multigrid step start -------
#ifndef VANKA
    if (rest>1.e+20) {std::cout << "\033[38;5;196m  ----------residual is too high----------  \033[0m" <<  endl; return;}
    rest= MGStep (_NoLevels-1,1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post, clearing);  // MGStep
//       rest=0.; MGCheck(_NoLevels-1); // check projection-restriction
#else
    rest= MGStep_Vanka (_NoLevels-1,1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);  // MGStep_Vanka
    //Vanka_test(0); // check projection-restriction
#endif
    if (rest <       Eps1* (1.+bNorm))   exit_mg=1;     // exit test
    cycle++;
#ifdef PRINT_CONV
    std::cout<<" cycle= "<< cycle  <<" residual= "<<rest<<" \n";
#endif
  }  // --------------------- end multigrid step  ------------------
  return;
}

// ====================================================================
/// This function does one multigrid step
double MGSolBase::MGStep (int Level,           // Level
                          double Eps1,          // Tolerance
                          int MaxIter,          // n iterations
                          const int Gamma,     // Control V W cycle
                          const int Nc_pre,    // n pre-smoothing cycles
                          const int Nc_coarse, // n coarse cycles
                          const int Nc_post,    // n post-smoothing cycles
                          const int clearing   // clean the init matrix and precond
                         ) {
// ====================================================================
  std::pair<int,double> rest (0,0.);
  if (Level == 0) {   // coarse level ----------------------------------
    rest=_solver[Level]->solve(*A[Level],*x[Level],*b[Level],Eps1,200,clearing);

#ifdef PRINT_CONV
    std::cout<<" Coarse res " << rest.second << " " << rest.first << std::endl;
#endif
    // coarse residual
    res[Level]->resid (*b[Level],*x[Level],*A[Level]);
  } // --------------------------------------------------------------
  else {// fine levels

    // presmoothing (Nu1) ---------------------------------
#ifdef PRINT_TIME  //  TC +++++++++++++++ 
    std::clock_t start_time=std::clock();
#endif             //  TC +++++++++++++++ 
    int  Nc_pre1=Nc_pre;
    if (Level<_NoLevels-1) Nc_pre1 *=2;
    rest=_solver[Level]->solve(*A[Level],*x[Level],*b[Level],Eps1, Nc_pre1,clearing);

#ifdef PRINT_CONV      //  CC +++++++++++++++ 
    std::cout<<" Pre Lev " << Level << " res " << rest.second << " " << rest.first;
#endif                 //  CC +++++++++++++++ 
#ifdef PRINT_TIME      //  TC +++++++++++++++      
    std::clock_t end_time=std::clock();
    std::cout<< " time ="<< double (end_time- start_time) / CLOCKS_PER_SEC << std::endl;
#endif                //  TC +++++++++++++++      
    // presmoothing residual
    res[Level]->resid (*b[Level],*x[Level],*A[Level]);
    // --------- end presmoothing (Nc_pre) ------------------------
    // restriction
    b[Level-1]->matrix_mult (*res[Level],*Rst[Level-1]);
    //  solving of system of equations for the residual on the coarser grid
    x[Level-1]->zero();
    double coarser_rest;
    for (int g=1; g <= Gamma; g++)
      coarser_rest=MGStep (Level-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post,clearing);
    // interpolation of the solution from the coarser grid (projection)
    res[Level]->matrix_mult (*x[Level-1],*Prl[Level]);
    // adding the coarse solution
    if (coarser_rest < rest.second) x[Level]->add (0.95,*res[Level]); // add the coarser solution only if it helps
    // postsmoothing (Nc_post) --------------------------------------------
#ifdef PRINT_TIME              //  TC +++++++++++++++ 
    start_time=std::clock();   //   initial set
#endif                         //  TC +++++++++++++++
    //  postsmooting
    int Nc_post1= Nc_post ;
    if (Level<_NoLevels-1) Nc_post1 *=2;
    rest=_solver[Level]->solve (*A[Level],*x[Level],*b[Level],Eps1,Nc_post1,clearing);

#ifdef PRINT_CONV              //  CC +++++++++++++++ 
    std::cout<<" Post Lev " << Level << " res " << rest.second << " " << rest.first;
#endif                         //  CC +++++++++++++++ 
#ifdef PRINT_TIME              //  TC +++++++++++++++
    end_time=std::clock();
    std::cout<< " time ="<< double (end_time- start_time) / CLOCKS_PER_SEC << std::endl;
#endif                         //  TC +++++++++++++++
    //  postsmoothing residual
    res[Level]->resid (*b[Level],*x[Level],*A[Level]);
    // ----------------  end postsmoothing ---------------------------
  }
  // end cycle -------------------------------------
  res[Level]->close();
  return  rest.second;

}

// =============================================
/// Check for Prolong and Restr Operators
void MGSolBase::MGCheck (int Level) const {
  std::cout<<"\nxlevel-1 before rest\n";
  x[Level-1]->print();
  std::cout<<"\n x level before rest\n";
  x[Level]->print();
  x[Level-1]->matrix_mult (*x[Level],*Rst[Level-1]);
  std::cout<<"\nxlevel-1 after rest\n";
  x[Level-1]->print();
  x[Level]  ->matrix_mult (*x[Level-1],*Prl[Level]);
  std::cout<<"\n x level after prol\n";
  x[Level]->print();
  //   x[Level-1]->matrix_mult(*x[Level],*Rst[Level-1]);
//   x[Level]  ->matrix_mult(*x[Level-1],*Prl[Level]);
  return;
}
// =============================================
/// Print xml attrib
void MGSolBase::print_xml_attrib (
  std::ofstream &out,  //  file xdmf
  int nodes,
  int /*nelems*/,
  std::string file_name
) const { // ================================

  for (int ivar=0; ivar<_n_vars; ivar++)   {
    std::string var_name = _var_names[ivar];
    out << "<Attribute Name=\""<< var_name <<"\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
        << nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
    out << file_name
//       femus_dir << "/" << output_dir << basesol << "."
//       << setw(ndigits) << setfill('0') << t_step << ".h5"
        << ":" << var_name << "\n";
    out << "</DataItem>\n" << "</Attribute>";
  }
  return;
}


// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_sol (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int ivar=0; ivar<nvars; ivar++) {  //ivar is like idim
        
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+ (ivar+ivar0) *offset]; // dof from top level
      uold[id + (ivar+kvar0) *NDOF_FEM]= ( (*x_old[_NoLevels-1]) (kdof_top)); // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_sol_piece (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int iel,  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int
       id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int ivar=0; ivar<nvars; ivar++) {  //ivar is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][iel*el_nds+id+ (ivar+ivar0) *offset]; // dof from top level
      uold[ id + (ivar+kvar0) *NDOF_FEM]= ( (*x_old[_NoLevels-1]) (kdof_top)); // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_oldsol (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+ (ivar+ivar0) *offset]; // dof from top level
      uold[ id + (kvar0+ivar) *NDOF_FEM]= ( (*x_oold[_NoLevels-1]) (kdof_top)); // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}
// ==========================================================================================
/// This function gets  the solution  vector at the nodes of  an element.
void  MGSolBase::get_el_nonl_sol (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int
       id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+ (ivar+ivar0) *offset]; // dof from top level
      double val= ( (*x_nonl[_NoLevels-1]) (kdof_top));
      uold[ id + (kvar0+ivar) *NDOF_FEM]=val;   // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_disp(
  const int Level,      // Level  <-
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(ivar+ivar0)*offset]; // dof from top level
      uold[ id +(kvar0+ivar)*NDOF_FEM]= ((*disp_old[Level])(kdof_top));     // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_disp(
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(ivar+ivar0)*offset]; // dof from top level
      uold[ id +(kvar0+ivar)*NDOF_FEM]= ((*disp_old[_NoLevels-1])(kdof_top));     // element sol
    } // end quadratic ------------------------------------------------
//     disp_old[_NoLevels-1]->print();
  }
  return;
}

// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_oooldsol (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+ (ivar+ivar0) *offset]; // dof from top level
      uold[ id + (kvar0+ivar) *NDOF_FEM]= ( (*x_ooold[_NoLevels-1]) (kdof_top)); // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_new_disp(
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(ivar+ivar0)*offset]; // dof from top level
      uold[ id +(kvar0+ivar)*NDOF_FEM]= ((*disp[_NoLevels-1])(kdof_top));     // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}
// ==========================================================================================
/// This function gets  the dof , the bc and the solution  vector at the nodes of  an element.
/// Note that indx_loc = id +ivar*NDOF_FEM with NDOF_FEM max dof (quad)
void  MGSolBase::get_el_oold_disp(
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const { // ==============================================================
  for (int id=0; id<el_nds; id++)    {
    // quadratic -------------------------------------------------
    for (int
         ivar=0; ivar<nvars; ivar++) {  //ivarq is like idim
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(ivar+ivar0)*offset]; // dof from top level
      uold[ id +(kvar0+ivar)*NDOF_FEM]= ((*disp_oold[_NoLevels-1])(kdof_top));     // element sol
    } // end quadratic ------------------------------------------------
  }
  return;
}

/// ======================================================
/// This function controls the time step operations:
/// ======================================================
void MGSolBase::MGTimeStep (const double time, const int) {

  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if  PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs (time,_NoLevels-1,1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for (int Level = 0 ; Level < _NoLevels-1; Level++) GenMatRhs (time,Level,0);

#if    PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << " Assembly time ="<< double (end_time- start_time) / CLOCKS_PER_SEC
            << " s "<< std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve (1.e-6,40);
#if    PRINT_TIME==1
  std::clock_t end_timef=std::clock();
  std::cout << " Assembly+solution time ="<< double (end_timef- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
  /// [d] Update of the old solution at the top Level
  x[_NoLevels-1]->localize (*x_old[_NoLevels-1]);
  return;
}


void MGSolBase::MGTimeStep_no_up (const double time, const int) {

  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if  PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs (time,_NoLevels-1,1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for (int Level = 0 ; Level < _NoLevels-1; Level++) GenMatRhs (time,Level,0);

#if    PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << " Assembly time ="<< double (end_time- start_time) / CLOCKS_PER_SEC
            << " s "<< std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve (1.e-6,40);
#if    PRINT_TIME==1
  std::clock_t end_timef=std::clock();
  std::cout << " Assembly+solution time ="<< double (end_timef- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif

  return;
}

void MGSolBase::MGUpdateStep () {

  std::cout  << std::endl << "  " << _eqname.c_str() << " update solution "  << std::endl;
  x[_NoLevels-1]->localize (*x_old[_NoLevels-1]);
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

double  MGSolBase::MGFunctional (
  double /*parameter*/,			/// Use of the function: (0) compute functional OR (1) set _eta
  double &/*control*/  /// \param[in] <>  eta multiplier for optimal method
) {
  std::cout << "Wrong use of MGFunctional from MGSolBase.C, aborting";
  abort();
  return 1;
}
  /// =======================================================================
  /// This function sets the controlled domain for optimal control problems
  /// =======================================================================
  void MGSolBase::set_ctrl_dom(
    const double xMin,const double xMax,
    const double yMin,const double yMax,
    const double zMin,const double zMax
  ){
    int    el_conn[NDOF_FEM];
    double x_m[DIMENSION];
    double xx_qnds[NDOF_FEM*DIMENSION];
    const int n_elem=_mgmesh._off_el[0][_NoLevels+_NoLevels*(_mgmesh._n_subdom-1)];
    double eps=1.e-6;  //tolerance for coordinates
    const int nel_e = _mgmesh._off_el[0][_NoLevels-1+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][_NoLevels-1+_NoLevels*_iproc];   // stop element

    _weight_ctrl.resize(n_elem);
    for(int iel=0; iel < (nel_e - nel_b); iel++) {
      _mgmesh.get_el_nod_conn(0,_NoLevels-1,iel,el_conn,xx_qnds);         //gets element coordinates
      for(int idim=0; idim< DIMENSION; idim++) {
	x_m[idim]=0;
	for(int d=0; d< NDOF_FEM; d++) {x_m[idim] +=xx_qnds[idim*NDOF_FEM+d]/NDOF_FEM;}//end d loop
      }//end idim loop
      _weight_ctrl[iel+nel_b]=0;   //default is 0
      if (x_m[0]>xMin-eps && x_m[0]<xMax+eps &&x_m[1]>yMin-eps && x_m[1]<yMax+eps
#if DIMENSION==3
	 &&x_m[2]>zMin-eps && x_m[2]<zMax+eps
#endif
        ){
        _weight_ctrl[iel+nel_b]=1;  //1 inside control region
// 	std::printf("iel %4d wieght %4f xm0 %4.5f xm1 %4.5f xm2 %4.5f\n",iel+nel_b, _weight_ctrl[iel+nel_b],x_m[0],x_m[1],x_m[2]);
	}
    } //end iel loop
    return;
  }

/// ======================================================
/// This function change current dt
/// ======================================================

double MGSolBase::GetValue (int flag) {}

void   MGSolBase::SetValue (double value) {}

void   MGSolBase::SetValueVector (std::vector<double> value) {}

void   MGSolBase::set_dt (double /*dt*/) {}

double MGSolBase::CalcFUpwind (double VelOnGauss[], double PhiDer[], double Diffusivity, int Dim, int NbOfNodes) {

  double vel_modulus =1.e-10, h_eff = 1.e-20, f_upwind = 0.;

  for (int i = 0; i<Dim; i++)  vel_modulus += VelOnGauss[i]*VelOnGauss[i];
  vel_modulus = sqrt (vel_modulus);

  for (int i=0; i<NbOfNodes; i++) {
    double hh=1.e-20;
    for (int idim=0; idim< Dim; idim++) hh += VelOnGauss[idim]*PhiDer[i+idim*NbOfNodes]/vel_modulus;
    h_eff += fabs (hh);
  }
  h_eff=2./h_eff;
  if (h_eff<1.e-10) {
    h_eff=1. ;
    std::cout << h_eff << " <1.e-10 in SUPG !!!!!!!!!\n";
  }
  // STANDARD SUPG
  const double Pe_h     = 0.5*vel_modulus*h_eff/ (Diffusivity);
  const double a_opt    = (1./tanh (Pe_h)-1./Pe_h);
  if (a_opt >1.) {
    std::cout << a_opt << " a_opt >1 in SUPG !!!!!!!!!\n";
  }
//   f_upwind = 0.5*a_opt*h_eff/ (vel_modulus);
  f_upwind = 0.5*a_opt*h_eff/ (vel_modulus);
  return f_upwind;
}


// ========================================================================

// ========================================================================
// ****************** HERE STARTS VANKA SECTION ***************************
// ========================================================================
#include "numeric_vectorM.h"
#include "sparse_matrixM.h"
#include "sparse_MmatrixM.h"
#include "dense_vectorM.h"
#include "dense_matrixM.h"
#include "petsc_matrixM.h"
#include "petsc_macroM.h"
#include "petsc_linear_solverM.h"
#include "petsc_preconditionerM.h"
#include "petsc_vectorM.h"
#include "petsc_matrixM.h"
//
// ====================================================================
/// This function does one multigrid step
// double MGSolBase::Vanka_test(int Level
//                             ) {
// // ====================================================================
//   std::pair<int,double> rest(0,0.);
//   if (Level == 0) {   // coarse level ----------------------------------
//     Vanka_solve(Level,*A[Level],*x[Level],*b[Level],1.e-6,40);
// // #ifdef PRINT_CONV
// //     std::cout<<" Coarse res " << rest.second << " " << rest.first << std::endl;
// // #endif
//   } // --------------------------------------------------------------
//   ;
//   return  0;
// }

// ====================================================================
/// This function does one multigrid step with Vanka
// double MGSolBase::MGStep_Vanka(
//   int Level,            // Level
//   double Eps1,          // Tolerance
//   int MaxIter,          // n iterations
//   const int Gamma,     // Control V W cycle
//   const int Nc_pre,    // n pre-smoothing cycles
//   const int Nc_coarse, // n coarse cycles
//   const int Nc_post    // n post-smoothing cycles
//   ) {
// // ====================================================================
//   std::pair<int,double> rest(0,0.);
//   if (Level == 0) {   // coarse level ----------------------------------
// //       b[Level]->close(); double bNorm=b[Level]->l2_norm();
// //       x[Level]->close(); double xNorm=x[Level]->l2_norm();
// //       A[Level]->close(); double aNorm=A[Level]->linfty_norm();
//     // coarse solution
// //    std::pair<int,double> rest1(0.,0);
//     Vanka_solve(Level,*A[Level],*x[Level],*b[Level],Eps1,40);
//     // coarse residual
//     res[Level]->resid(*b[Level],*x[Level],*A[Level]);
// #ifdef PRINT_CONV
//     std::cout<<" Coarse res " << 200 << " " <<  res[Level]->l2_norm() << std::endl;
// #endif
//
//   } // --------------------------------------------------------------
//   else {// fine levels
//
//     // presmoothing (Nu1) ---------------------------------
// #ifdef PRINT_TIME  //  TC +++++++++++++++
//     std::clock_t start_time=std::clock();
// #endif             //  TC +++++++++++++++
//     int  Nc_pre1=Nc_pre;
//     if (Level<_NoLevels-1) Nc_pre1 *=2;
//     Vanka_solve(Level,*A[Level],*x[Level],*b[Level],Eps1, Nc_pre1);
// #ifdef PRINT_CONV      //  CC +++++++++++++++
//     std::cout<<" Pre Lev " << Level << " res " << rest.second << " " << rest.first;
// #endif                 //  CC +++++++++++++++
// #ifdef PRINT_TIME      //  TC +++++++++++++++
//     std::clock_t end_time=std::clock();
//     std::cout<< " time ="<< double(end_time- start_time) / CLOCKS_PER_SEC << std::endl;
// #endif                //  TC +++++++++++++++
//     // presmoothing residual
//     res[Level]->resid(*b[Level],*x[Level],*A[Level]);
//     // --------- end presmoothing (Nc_pre) ------------------------
//
//     // restriction
//     b[Level-1]->matrix_mult(*res[Level],*Rst[Level-1]);
//
//     //  solving of system of equations for the residual on the coarser grid
//     x[Level-1]->zero();
//     for (int g=1; g <= Gamma; g++)
//       MGStep_Vanka(Level-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);
//
//     // interpolation of the solution from the coarser grid (projection)
//     res[Level]->matrix_mult(*x[Level-1],*Prl[Level]);
//     // adding the coarse solution
//     x[Level]->add(*res[Level]);//  *x[Level] +=*res[Level];
//
//     // postsmooting (Nc_post) --------------------------------------------
// #ifdef PRINT_TIME              //  TC +++++++++++++++
//     start_time=std::clock();   //   initial set
// #endif                         //  TC +++++++++++++++
//     //  postsmooting
//     int Nc_post1= Nc_post ;
//     if (Level<_NoLevels-1) Nc_post1 *=2;
//     Vanka_solve(Level,*A[Level],*x[Level],*b[Level],Eps1,Nc_post1);
// #ifdef PRINT_CONV              //  CC +++++++++++++++
//     std::cout<<" Post Lev " << Level << " res " << rest.second
//              << " " << rest.first;
// #endif                         //  CC +++++++++++++++
// #ifdef PRINT_TIME              //  TC +++++++++++++++
//     end_time=std::clock();
//     std::cout<< " time ="<< double(end_time- start_time) / CLOCKS_PER_SEC << std::endl;
// #endif                         //  TC +++++++++++++++
//     //  postsmooting residual
//     res[Level]->resid(*b[Level],*x[Level],*A[Level]);
//     // ----------------  end postsmoothing ---------------------------
//   }
//   // end cycle -------------------------------------
//   res[Level]->close();
//    double norm2= res[Level]->l2_norm();
// //    std::cout<< " True res l2norm (not prec) " << norm2<< std::endl;
//   return   norm2;
//
// }



// ========================================================
void MGSolBase::Vanka_solve(
  int Level,
  SparseMatrixM&  matrix_in,
  NumericVectorM& solution_in,
  NumericVectorM& rhs_in,
  const double tol,
  const  int m_its
) {

//   START_LOG("solve()", "PetscLinearSolverM");
  // Make sure the data passed in are really of Petsc types
  PetscMatrixM* matrix   = libmeshM_cast_ptr<PetscMatrixM*>(&matrix_in);
  PetscVectorM* solution = libmeshM_cast_ptr<PetscVectorM*>(&solution_in);
  PetscVectorM* rhs      = libmeshM_cast_ptr<PetscVectorM*>(&rhs_in);
//   this->init(matrix);
//   int Level=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;
  // Close the matrices and vectors in case this wasn't already done.
  matrix->close();    solution->close(); rhs->close();

  PetscErrorCode ierr;
  PetscInt m=1;
  PetscInt n=1;
//   const PetscInt idxm[1]={9};
//   const PetscInt idxn[1]={9};
  PetscScalar v[1];


  /// a) Set up
  // geometry -----------------------------------------------------------------------------------
  const int  ndim = DIMENSION;                                           //dimension
  const int  offset = _mgmesh._NoNodes[0];                     // mesh nodes
//   const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  double     normal[DIMENSION]; double    mu_m;                                          // normal to the boundary

  // Gauss integration ---------------------------------------------------------------------------
//   const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points
//   const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           // Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];  double dphiidx_g[3][DIMENSION]; // global derivatives at g point

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] -----------------------------------

  PetscInt ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt idxm[1];

//   DenseMatrixM KeM;    DenseVectorM FeM;                              // local  matrix+rhs
//   KeM.resize(el_mat_nrows,el_mat_ncols);    FeM.resize(el_mat_nrows); // resize  local  matrix+rhs

  vector<double> x_loc[offset];
  solution->localize(*x_loc);

  vector<double> b_loc[offset];
  rhs->localize(*b_loc);


  int nsm=2;
  for (int ismooth=0; ismooth<= nsm; ismooth++) {
    //element type loop
    const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1];
    const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];
    for (int iel=0; iel < (nel_e - nel_b); iel++) {

      /// b) Element  Loop over the volume (n_elem)
      // set to zero matrix and rhs
//     KeM.zero();         FeM.zero();

      // geometry element quantities ---------------------------
      // Element connectivity  and coordinates (xx_qnds)
      _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);  // get connectivity and coord
//     _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);// get neighboring element
//       std::cout <<  " element nodes \n" ;
      for (int inode=0; inode <NDOF_FEM; inode++) {
//         std::cout << el_conn[inode] << " " ;


        idxm[0]=el_conn[inode];
        ierr=MatGetValues(matrix->mat(),1, idxm,1, idxm,v);

        double sum= (*b_loc)[el_conn[inode]]+v[0]*(*x_loc)[el_conn[inode]];
        ierr=MatGetRow(matrix->mat(),el_conn[inode],&ncols,&cols,&vals);
        for (int i=0; i <ncols; i++)  sum -=vals[i]*(*x_loc)[cols[i]];
        sum /=v[0];

        ierr=MatRestoreRow(matrix->mat(),el_conn[inode],&ncols,&cols,&vals);
        (*x_loc)[el_conn[inode]]=sum;
	if(ismooth == nsm) {
// 	  std::cout << " sol " << sum << " node  " <<  el_conn[inode] << "  \n";
	  (*x[Level]).set(el_conn[inode],sum);
	}
      }  // inode

//      std::cout <<  " matrix  val \n" ;
    }
  }




  return;
}
