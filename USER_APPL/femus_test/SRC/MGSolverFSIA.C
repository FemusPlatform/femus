// ===============================================================
// --------------   NAVIER-STOKES system [FSA_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef FSIA_EQUATIONS
// ==============================================================
// FSIA_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSIA_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSIA_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Fluid Structure class conf file
#include "MGSolverFSIA.h"       // Fluid Structure class header file
#ifdef TBK_EQUATIONS
#include "MGSolverTBK.h"
#endif
#include <iomanip>      // std::setprecision
// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif



// ==================================================================
/// This routine constructs the FSI class:
MGSolFSIA::MGSolFSIA (
  MGEquationsSystem &mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
) :  MGSolDA (mg_equations_map_in,nvars_in,eqname_in,varname_in),
//
//===================================================================================================//
//                                    A) Reading parameters                                          //
//===================================================================================================//
//
  _offset (_mgmesh._NoNodes[-1]),            // mesh nodes (top level)
  _dt (stod (_mgutils._sim_config["dt"])),                      // parameter  dt
  _uref (_mgutils._mat_prop["Uref"]),        // parameter  u reference
  _lref (_mgutils._mat_prop["Lref"]),        // parameter  l reference
  _rhof (_mgutils._mat_prop["rho0"]),        // parameter density
  _muf (_mgutils._mat_prop["mu0"]) ,          // parameter viscosity
  _rhos(_mgutils._mat_prop["rhos"]),    // parameter density
  _ni(_mgutils._mat_prop["nis"]),
  _Emod(_mgutils._mat_prop["Es"]),
  _hs(_mgutils._mat_prop["hs"]) 
  {      
  //===================================================================================================//
  //                                    B) Setting class variables                                     //
  //===================================================================================================//

  // READ PARAMETERS FROM CLASS FSIA_parameter
  _FSIA_parameter.read_param (_mgutils);
  _nNSdim=DIMENSION;
  // class equation ---------------------------------------------------------------------
  for (int k_index=0; k_index<30; k_index++)     _FF_idx[k_index]=-1;
  _pres_order= (_nvars[0]>0) ? 0:1;
  // class variable names ------------------------------------------------------------
  _dir=0; 
#if FSIA_EQUATIONS==2       //   segregated ( P in NSP_EQUATIONS)
  if (!varname_in.compare ("ua"))     _dir=0;  // u-equation 
  if (!varname_in.compare ("va"))     _dir=1;  // v-equation  
  if (!varname_in.compare ("wa"))     _dir=2;  // w-equation  
  _var_names[0]=varname_in;  _refvalue[0]=_uref;
#else
  _var_names[_nNSdim-1]="wa";  _refvalue[_nNSdim-1]=_uref; // velocity 3D
  _var_names[0]="ua";          _refvalue[0]=_uref; // velocity 2D
  _var_names[1]="va";          _refvalue[1]=_uref; // velocity 2D
#if FSIA_EQUATIONS==1       //  coupled  (+P)
  _var_names[_nNSdim]="pa";   _refvalue[_nNSdim]=_rhof*_uref*_uref;  // pressure
#endif
#endif
  //===================================================================================================//
  //                                    C) Setting solver type                                         //
  //===================================================================================================//
  for (int l=0; l<_NoLevels; l++) {
    _solver[l]->set_solver_type (_FSIA_parameter._SolverType);
  }
  //===================================================================================================//
  //                                    D) Setting non_dimensional parameters  01                        //
  //===================================================================================================//
  double nu=_muf/_rhof;  if (nu<1.e-7) nu=1.e-7;
  _IRe=_muf/ (_rhof*_lref*_uref);         // Reynolds number
  _IFr=9.81*_lref/ (_uref*_uref);         // Froud number
  _dirg[0] = _mgutils._mat_prop["dirgx"];      // x-gravity
  _dirg[1] = _mgutils._mat_prop["dirgy"];      // y-gravity
  _dirg[2] = _mgutils._mat_prop["dirgz"];      // z-gravity
  for (int i=0; i<_nNSdim; i++) _xyzg[i]=0.;
  _y_bcout=sqrt (2*0.09*_IRe/80.);
  _lambda=(_ni*_Emod)/(_rhos*(1+_ni)*(1-2*_ni));//.001;
  _mus=_Emod/(_rhos*2*(1+_ni)*_lref*_uref);
  
  _SolveFSIA=(_mgutils._sim_config["SolveAdjointFluidStructure"].compare("yes") == 0) ? true: false;
  _Wall_dist =_mgutils._geometry["Wall_dist"];
  _AxiSym  = (int) (_mgutils._geometry["Axisym"]);

  if (_AxiSym == 1 && _nNSdim == 3) {
    std::cout<<"\033[1;31m MGSolFSIA: CANNOT USE AXISYM OPTION WITH 3D GEOMETRY \n \033[0m";
    std::abort();
  }
  return;
} //****************************************************************************************************//

// ====================================================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================================================
void  MGSolFSIA::GenMatRhs (const double/* time*/, const int
                          Level,const  int mode) {
// FSIA_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSIA_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSIA_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// =====================================================================================================

  // ===================================================================================================
  //                                    A) Set up
  // ===================================================================================================


  const double les=_FSIA_parameter._Les;
  // -------------------------- Geometry -------------------------------------------------------------
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     normal[DIMENSION];                                          // normal to the boundary
  double NodeMuTurb[NDOF_FEM], WallDist[NDOF_FEM];
  // -------------------------- Gauss integration -----------------------------------------------------

  const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] ---------------------------------------
  int el_ndof[3];  el_ndof[0]=NDOF_K;
  int elb_ndof[3];  elb_ndof[0]=0; // number of el dofs
  int el_mat_nrows =0;                                            // number of mat rows (dofs)
  for (int ideg=1; ideg<3; ideg++) {                              //     ...
    el_ndof[ideg]= ( (_nvars[ideg]>0) ?    _fe[ideg]->_NoShape[ _nNSdim-1]:0);                 //   computing
    elb_ndof[ideg]= ( (_nvars[ideg]>0) ?_fe[ideg]->_NoShape[ _nNSdim-2]:0);               //     ...
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  }
#if FSIA_EQUATIONS%2==0
  el_ndof[1]= _fe[1]->_NoShape[ _nNSdim-1];
#endif
  el_mat_nrows +=  el_ndof[0]*_nvars[0];
  int el_mat_ncols = el_mat_nrows;                                //square matrix
  std::vector<int> el_dof_indices (el_mat_ncols);                 // element dof vector

  // fields -> Fluid Structure ----------------------------------------------------------------------
//   double u_nlg[DIMENSION];
  double u_old[DIMENSION*NDOF_FEM];
  double u_oold[DIMENSION*NDOF_FEM];
  double u_nl[DIMENSION*NDOF_FEM];       // velocity vector for non linear terms -> it contains all the velocity components //
  double p_proj[DIMENSION*NDOF_FEM];
  double dp_proj[DIMENSION*NDOF_FEM];
  double x_m[DIMENSION];    int WallNodes[NDOF_FEMB];
  int flag_group[NDOF_FEM];
  MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
  for (int i0=0; i0<DIMENSION*NDOF_FEM; i0++) {
    p_proj[i0]=0.;    dp_proj[i0]=0.;
  }
// coupling  basic system fields -------------------------------------------------------------------------------
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for (int k=0; k<30; k++) { // coupling  basic system fields
    const int idx= _data_eq[2].tab_eqs[k];
    _FF_idx[k]= (idx>=0) ?_data_eq[2].indx_ub[idx]:-1;
  }
  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
  A[Level]->zero();  if (mode ==1)  b[Level]->zero();  // global matrix+rhs
  DenseMatrixM KeM; DenseVectorM FeM;                  // local  matrix+rhs
  KeM.resize (el_mat_nrows,el_mat_ncols);  FeM.resize (el_mat_nrows);         // resize  local  matrix+rhs

   int ndof_lev=0;   // for multilevel multi proc point ordering ------------------------------------------
   for (int pr=0; pr <_mgmesh._iproc; pr++)   
    ndof_lev +=_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];

  // ===================================================================================================
  //                                    B) Element  Loop over the volume (n_elem)
  // ===================================================================================================

  const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for (int iel=0; iel < (nel_e - nel_b); iel++) {
    // set to zero matrix and rhs and center01
    KeM.zero();    FeM.zero();
    
    // ---------------------------------------------------------------------------------
    // Volume  geometry and element  fields 
    // ---------------------------------------------------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn (0,Level,iel,el_conn,_xx_qnds);
    _mgmesh.get_el_neighbor (el_sides,0,Level,iel,el_neigh);
    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc (Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
    // field data  ------------------------------------------------------
    get_el_field_data (iel, Level,el_conn, offset,el_ndof,ndof_lev, u_old,u_oold, u_nl,p_proj,dp_proj);

   // initializing  volume quantities
    for (int idim=0; idim< _nNSdim; idim++) { // quad loop entities (vector)
       x_m[idim]=0.;
      for (int d=0; d<el_ndof[2]; d++) {
        const int  dnode=idim*el_ndof[2]+d;    // index local points
        x_m[idim] +=_xx_qnds[dnode]/el_ndof[2];
	flag_group[d]=fabs(ext_mesh->_bc_id[el_conn[d]]);
        _normal_pt[dnode]=0.;// normal boundary
        _bc_el[dnode]=0;    _AlreadyWrittenDirBC[dnode] = false;
      }
    }
#if FSIA_EQUATIONS==1
    for (int d=0; d< el_ndof[1]; d++) { // pressure  loop entities (scalar)
      if (_bc_bd[_nNSdim*el_ndof[2]+d]%10 >3)    _bc_el[_nNSdim*el_ndof[2]+d]=0;  
      else _bc_el[_nNSdim*el_ndof[2]+d]=1;
    }
#endif

   
      // SETTING _bc_el = -8 FOR EVERY ELEMENT BOUNDARY NODE
    for (int  iside=0; iside< el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        for (int  lbnode=0; lbnode<elb_ndof[2]; lbnode++) {
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+elb_ndof[2]*iside];// local nodes
          for (int idim=0; idim< _nNSdim; idim++) _bc_el[lnode + idim*el_ndof[2]] = _BdFlagId;
        }
      }
    }
    // element fields ----------------------------------
    if (_FF_idx[K_F]>=0) _y_bcout= _Wall_dist + _mgmesh._dist[ iel+nel_b]; //  _y_bcout=_mgmesh._dist[ iel+nel_b];
     _WallElement=0;
    
    // --------------------------------------------------------------------------
    //  Boundary and boundary bc 
    // --------------------------------------------------------------------------
    //  New bc flags for Navier Stokes equation:
    //  _bc_el = -1 -2 -3 -> Neumann bc along tangential 1, tangential 2 (only 3d), normal directions
    //  _bc_el =  1  2  3 -> Dirichlet bc along tangential 1, tangential 2 (only 3d), normal directions
    //  _bc_el =  0 -> normal interior point
    //  Integration is performed for _bc_el <=0
    //  Test functions are set to zero when _bc_el >0
    for (int  iside=0; iside< el_sides; iside++) {
      if (el_neigh[iside] == -1) {
        // setup boundary element  ----------------------------------------------------------------
        for (int  lbnode=0; lbnode<elb_ndof[2]; lbnode++) { // quad quantities
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+elb_ndof[2]*iside];// local nodes
          sur_toply[lbnode]=lnode;          // lbnode -> lnode
          elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
          for (int idim=0; idim< _nNSdim; idim++) { // coordinates
            _xxb_qnds[idim*elb_ndof[2]+lbnode]=_xx_qnds[idim*el_ndof[2]+lnode];
          }
        }
        // normal --------------------------------------------------------------------------------------
        _fe[2]->normal_g (_xxb_qnds,x_m,normal);
        int dir_maxnormal = (fabs (normal[0]) >fabs (normal[1])) ?0:1 ;
        dir_maxnormal= (fabs (normal[dir_maxnormal]) >fabs (normal[DIMENSION-1])) ? dir_maxnormal:DIMENSION-1;
       //         std::cout<<" .........................................................\n Element "<<iel<<" side "<< iside <<std::endl;
	int before  = _WallElement;
        set_bc_matrix (KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof,elb_ngauss,normal,u_old, el_conn);	
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        int after  = _WallElement;
        if (after * before == 0 && after + before == 1)   for (int idof=0; idof<elb_ndof[2]; idof++)      WallNodes[idof] = sur_toply[idof];
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                 std::cout << " bc after : " << iel <<"\t"<<iside<< " x "<< x_m[0] << " y" << x_m[1]<<std::endl;
      } // iside -1
    }  // -----------------------------  End Boundary -------------------------------------
    
//   std::cout << "-------------------------------------" << std::endl;
//                    std::cout << " element: " << iel << std::endl;
//             std::cout << std::setprecision(10) << "  KeM \n" << KeM << " rhs \n" <<  FeM <<"\n\n" <<std::endl;
// std::cout<<"Element "<<iel<<std::endl;
    // ------------------------------ Volume --------------------------------------------
   if(ext_mesh->_mat_id[ iel+nel_b]  !=4){  // liquid region
      matrixrhs_liq_vol (KeM, FeM, el_ndof, u_old, u_nl, 0, mode, NodeMuTurb, el_conn);
   }
   else{
      matrixrhs_sol_vol(KeM, FeM, el_ndof, u_old, u_nl, 0, mode, NodeMuTurb, el_conn,flag_group,iel+nel_b);
   }
    
    // ---------------------- end volume (element) --------------------------------------

// NORMALIZATION BOUNDARY BC EQUATIONS
// if boundary stresses are imposed, then with matrix normalization symmetry is lost
//    double max_diag=0.; int count=0;
//      for(int id=0; id< el_ndof[2]*_nvars[2]; id++) {
//        if(_bc_el[id] <= 0) if(max_diag<fabs(KeM(id,id))) { max_diag=fabs(KeM(id,id)); }
//      }
//      for(int id=0; id< el_ndof[2]*_nvars[2]+el_ndof[1]; id++)  if(_bc_el[id] >0) {
//        const double pivot=fabs(KeM(id,id));
//        FeM(id) *=  max_diag/pivot;
//        for(int jd=0; jd<el_ndof[2]*_nvars[2]+el_ndof[1]; jd++) { KeM(id,jd) *= max_diag/pivot; }
//      }
//      for(int id=el_ndof[2]*_nvars[2]; id< el_ndof[2]*_nvars[2]+el_ndof[1]; id++) if(_bc_el[id]==0) {
//        const double pivot=fabs(KeM(id,id));
//        FeM(id) *=  max_diag/pivot;
//        for(int jd=0; jd<el_ndof[2]*_nvars[2]; jd++) { KeM(id,jd) *= max_diag/pivot; }
//      }

//     std::cout << " element: " << iel << std::endl;
//     std::cout << std::setprecision(10) << "  KeM \n" << KeM << " rhs \n" <<  FeM <<"\n\n" <<std::endl;

    // ----------------------------------------------------------------------------------
    //   E) Add them to the global matrix
    A[Level]->add_matrix (KeM,el_dof_indices);   if (mode == 1)   b[Level]->add_vector (FeM,el_dof_indices);
    
  } //  =============== End of element loop =============================================

// ===============================================================================
  el_dof_indices.clear();
  A[Level]->close();  if (mode == 1)     b[Level]->close();
// ----------------------------------------------------------
#ifdef PRINT_INFO
//     A[Level]->print();  b[Level]->print();
  std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

  return;
} //****************************************************************************************************//


// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolFSIA::MGTimeStep (
  const double time,  // time
  const int max_iter  // Number of max inter
) {
// =========================================================================================

  // ========================================================================================= //
  //              A) Set up the time step                                                      //
  // ========================================================================================= //
  if (_SolveFSIA) {
    std::cout  << std::endl << "\033[038;5;"<<FSA_F + 50<<";1m "
               << "--------------------------------------------------- \n\t"
               <<  _eqname.c_str()
               << " solution of problem " << _mgutils.get_name()
               << "\n ---------------------------------------------------\n\033[0m";
    for(int k=0; k<30; k++) {// coupling  basic system fields
    const int idx= _data_eq[2].tab_eqs[k];
    _FF_idx[k]=(idx>=0)?_data_eq[2].indx_ub[idx]:-1;
  }
  for(int kdim=0;kdim<_nNSdim;kdim++) {
     const int num=_data_eq[2].tab_eqs[SDSX_F+kdim];
   _mgmesh.Translate(kdim,(*_data_eq[2].mg_eqs[num]->disp[_NoLevels-1]));   //Moves the mesh (interpolation of mid-points)
  }
    x[_NoLevels-1]->localize (*x_nonl[_NoLevels-1]);
    // ========================================================================================= //
    //              B) Assemblying of the Matrix-Rhs                                             //
    // ========================================================================================= //

#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif
    GenMatRhs (time,_NoLevels-1,1);                                             // matrix and rhs
    for (int Level = 0 ; Level < _NoLevels-1; Level++) GenMatRhs (time,Level,0);  // matrix
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double (end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif
    // ========================================================================================= //
    //              C) Solution of the linear MGsystem (MGSolFSIA::MGSolve)                       //
    // ========================================================================================= //
    MGSolve (1.e-6,40);
#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double (end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
//   x_oold[_NoLevels-1]->close();
//   x_oold[_NoLevels-1]->zero();
//   double norm_check=  x_oold[_NoLevels-1]->l2_norm();
//   if(norm_check>2.e+5) { x_oold[_NoLevels-1]->scale(1.e-3); }
//   double norm_diff=1.;
//
    // ========================================================================================= //
    //              D) Non Linear Iterations                                                     //
    // ========================================================================================= //

    int iter_nl=_FSIA_parameter._MaxNonLinearIterations;
//   if(time < _FSIA_parameter.NL_TIME0) { iter_nl=_FSIA_parameter.NL_ITER0; }
//   double penalty =1.e+1; double factor =1.e+0;
    for (int iter=0; iter<iter_nl; iter++) {
      x[_NoLevels-1]->localize (*x_nonl[_NoLevels-1]);
//     double norm_in= x_nonl[_NoLevels-1]->l2_norm();
//     std::cout <<" norm xnl*********** " << norm_in  << std::endl;
// #if FSIA_EQUATIONS%2==0
//     x_oold[_NoLevels-1]->add(penalty,*x_nonl[_NoLevels-1]); // penalty
// //     penalty*=factor;
// #endif
//     std::cout <<" norm x_oold--------- " <<  x_oold[_NoLevels-1]->l2_norm()  << std::endl;
//
      // ----------- Non linear solution ------------------------------------------------------- //
//
      std::cout  << std::endl << " === NON LINEAR ITERATOR: "<< iter+1 << "- " << _eqname.c_str() << " solution "  << std::endl;
//     A[_NoLevels-1]->close();
      GenMatRhs (time,_NoLevels-1,1);                                               // matrix and rhs
      for (int Level = 0 ; Level < _NoLevels-1; Level++)      GenMatRhs (time,Level,0);  // matrix
      MGSolve (1.e-6,15);                                                           // solve
//
      // ----------- Check error --------------------------------------------------------------- //
//
//     x[_NoLevels-1]->localize(*disp[_NoLevels-1]);
//     disp[_NoLevels-1]->add(-1.e+0,*x_nonl[_NoLevels-1]);
//     norm_diff=  disp[_NoLevels-1]->l2_norm();
//     std::cout << " Check:  "<< " norm_in "<< norm_in <<
//               "  norm_diff/-norm_in " <<norm_diff/norm_in <<std::endl
//               <<" === END NON LINEAR ITERATOR"<<std::endl;
//     if(norm_diff/norm_in<1.e-6) { break; }

    }
    x_old[_NoLevels-1]->localize (*x_oold[_NoLevels-1]);
    x[_NoLevels-1]->localize (*x_old[_NoLevels-1]);
#if (FSIA_EQUATIONS%2==0)
    x[_NoLevels-1]->localize (*x_nonl[_NoLevels-1]);
#endif
  }
  return;
} //****************************************************************************************************//
//
//

void MGSolFSIA::CalcMuTurbAndWallDist (double MuTurb[], double WallDist[], int WallNodes[], int ElId, int NumOfNodes, int Level) {
#ifdef HAVE_MED
  std::vector<int> NodeConn;
  _mgutils._TurbParameters->_MuTurbField[Level]->getMesh()->getNodeIdsOfCell (ElId,NodeConn);

 // CALCULATION OF NODE WALL DIST
  for (int f=0; f<NDOF_FEM; f++)  WallDist[f] = _mgutils._TurbParameters->_NodeWallDist[Level]->getIJ (NodeConn[f],0);
 

  if (_WallElement==0) {// CALCULATION OF MU TURB NODE VALUES -> READING FROM TURBPARAMETERS MUTURB FIELD
    for (int f=0; f<NDOF_FEM; f++) MuTurb[f] = _mgutils._TurbParameters->_MuTurbField[Level]->getIJ (NodeConn[f],0);
  }
  else {
    for (int f=0; f<NDOF_FEM; f++) {
      double TmpDist = 100.;   double TmpDist2 = 0.;
      for (int idof=0; idof<NumOfNodes; idof++) {
        TmpDist2 = 0.;  for (int dir=0; dir<DIMENSION; dir++)
          TmpDist2 += (_xx_qnds[dir*NDOF_FEM+f] - _xx_qnds[dir*NDOF_FEM+WallNodes[idof]]) * (_xx_qnds[dir*NDOF_FEM+f] - _xx_qnds[dir*NDOF_FEM+WallNodes[idof]]);

        TmpDist2 = sqrt (TmpDist2);
        if (TmpDist > TmpDist2) TmpDist = TmpDist2;
      }
      WallDist[f] = TmpDist;
//           printf (" el %d node %d walldist %12.8f \n",iel,f,WallDist[f]);
    }

    double Vel[DIMENSION], VelNorm[DIMENSION], UMag = 1.e-20;
    double norm_mag = 1.e-20, UNormMag = 1.e-20;
    for (int dir=0; dir<DIMENSION; dir++) {
      Vel[dir] = _data_eq[2].ub[ (_FF_idx[FSA_F]+dir) *NDOF_FEM + (NDOF_FEM -1)];
      norm_mag +=  _normal_pt[ (NDOF_FEM -1) *_nNSdim+dir]* _normal_pt[ (NDOF_FEM -1) *_nNSdim+dir];
    }
    norm_mag = sqrt (norm_mag);

    for (int dir=0; dir<DIMENSION; dir++) {
      _normal_pt[ (NDOF_FEM -1) *_nNSdim+dir] /= norm_mag;
      VelNorm[dir] = Vel[dir]* _normal_pt[ (NDOF_FEM -1) *_nNSdim+dir];
      UNormMag    += VelNorm[dir]*VelNorm[dir];
      UMag        += Vel[dir]*Vel[dir];
    }

    double UTangMag = sqrt (UMag - UNormMag);
    double utau     = _mgutils._TurbParameters->CalcUtau (UTangMag,WallDist[NDOF_FEM -1]);
    double yplus    = WallDist[NDOF_FEM -1]*utau/_IRe;
    _mgutils._TurbParameters->_utau = utau;
    for (int f=0; f<NDOF_FEM; f++) { // CALCULATION OF MU TURB NODE VALUES -> READING FROM TURBPARAMETERS MUTURB FIELD
      yplus = WallDist[f]*utau/_IRe;
      double k = utau*utau* 1./(1./(0.041*yplus*yplus) + 1./(1./0.3));
      const double wlin  = 2.*_IRe/ (0.09*WallDist[f]*WallDist[f]);
      const double wlog  = utau/ (0.41 * WallDist[f] * sqrt (0.09));
      double w= wlog;
    
      double KandW[2];
      KandW[0] = (1.-_mgutils._TurbParameters->_klog) *k + _mgutils._TurbParameters->_klog*log (k);
      KandW[1] = (1.-_mgutils._TurbParameters->_wlog) *w + _mgutils._TurbParameters->_wlog*log (w);
      
      MuTurb[f] = 1./ (1./ (0.001093*yplus*yplus*yplus) + 1./ (0.41*yplus));
//       MuTurb[f] = _mgutils._TurbParameters->CalcMuTurb(KandW, WallDist[f]);
//       printf("yplus %12.8f     utau %12.8f      mu turb %12.8f\n",yplus,utau,MuTurb[f]);
    }
  }
  #endif
  return;
}

// ------------------------------------------------------------------------------
void FSIA_param::read_param(
  MGUtils &mgutils
) {
  //  Reading parameter  -> FSIAproperties.in ------------------------------------------------
  read_file();

  // Boundary condition block ------------------------------------------------------------
  std::cout<<"Adjoint Fluid Structure boundary condition block \n";
  // boundary group names  from
  std::string GroupString = mgutils._sim_config["BoundaryGroups"];
  int count=0; int pos1=0; int Length = GroupString.length();
  while(count<Length) {
    if(GroupString.at(count)==',') {
      std::string  temps = GroupString.substr(pos1,count-pos1);
      _BoundaryGroupsIDs.push_back(stoi(temps));
      pos1=count+1;
    }
    count++;
  }
  // boundary group values from
  _BoundaryGroupsIDs.push_back(stod(GroupString.substr(pos1,Length-pos1)));
  for(int i=0; i<_BoundaryGroupsIDs.size(); i++) {
    std::string BDcond = "FSIAgroup"+to_string(_BoundaryGroupsIDs[i]);
    _map_FSIAgroup[_BoundaryGroupsIDs[i]]=_BoundMap[_FileMap[BDcond]];
  }
  // ---------------------  end bc ------------------------------------
  // Fluid Structure parameters block ------------------------------------------------------------
  std::cout<<" Adjoint FLUID STRUCTURE PARAMETER map (_FSIA_parameter) in  FSIAproperties.in + UserFSIA.h :\n";
  if(_FileMap ["SolveSteady"]!="") { _SolveSteady  = stoi(_FileMap["SolveSteady"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._SolveSteady (in UserFSIA.h) \n"; }
  if(_FileMap ["Compressible"]!="") { _Compressible  = stoi(_FileMap["Compressible"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._Compressible (in UserFSIA.h) \n"; }
  if(_FileMap ["MaxNonLinearIterations"]!="") { _MaxNonLinearIterations   =stoi(_FileMap["MaxNonLinearIterations"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._MaxNonLinearIterations (in UserFSIA.h)\n"; }
  if(_FileMap ["DynamicUnderRelaxation"]!="") { _UnderRelaxation =stof(_FileMap["DynamicUnderRelaxation"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._DynamicUnderRelaxation (in UserFSIA.h)\n"; }
  std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._SolverType  (in UserFSIA.h) \n"; 
  
  if(_FileMap ["Supg"]!="") { _Supg= stoi(_FileMap["Supg"]); }
  else { std::cout <<" FSIAproperties.in:  default value for _FSIA_parameter.Supg (in UserFSIA.h)\n"; }
  if(_FileMap ["Upwind"]!="")  _Upwind=stof(_FileMap["Upwind"]);
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._Upwind (in UserFSIA.h)\n"; }
  if(_FileMap ["Les"]!="")  {_Les=stof(_FileMap["Les"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._Les (in UserFSIA.h)\n"; }
  if(_FileMap ["ReactionNumberBased"]!="") { _ReactionNumberBased     = stoi(_FileMap["ReactionNumberBased"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._ReactionNumberBased  (in UserFSIA.h)\n"; }

  if(_FileMap ["FlatProfile"]!="")  _FlatProfile = stoi(_FileMap["FlatProfile"]);
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._FlatProfile (in UserFSIA.h)\n"; }
  if(_FileMap ["InterpolatedMuTurb"]!="")  _InterpolatedMuTurb= stoi(_FileMap["InterpolatedMuTurb"]);
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._InterpolatedMuTurb (in UserFSIA.h)\n"; }

  if(_FileMap ["WallFunctionApproach"]!="") _WallFunctionApproach= stoi(_FileMap["WallFunctionApproach"]);
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter._WallFunctionApproach (in UserFSIA.h)\n"; }
  if(_FileMap ["Penalty_n"]!="") { _Penalty_n = stof(_FileMap["Penalty_n"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter.Penalty_n (in UserFSIA.h)\n"; }
  if(_FileMap ["Penalty_tg"]!="") { _Penalty_tg=stof(_FileMap["Penalty_tg"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter.Penalty_tg (in UserFSIA.h)\n"; }
  
    if(_FileMap ["Tg1_stress"]!="") { _Tg1_stress = stof(_FileMap["Tg1_stress"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter.Tg1_stress (in UserFSIA.h)\n"; }
    if(_FileMap ["Tg2_stress"]!="") { _Tg2_stress = stof(_FileMap["Tg2_stress"]); }
  else { std::cout <<" FSIAproperties.in: default value for _FSIA_parameter.Tg2_stress (in UserFSIA.h)\n"; }
// ---------------------  end parameter NS ------------------------------------
  _FileMap.clear();
  return;
}

// ============================================================
/// This function reads NS parameter file (NSparameter.in)
void FSIA_param::read_file(
) {//  ===========================================================

  //  getting file name --------------------------------------------------------------------
  std::ostringstream file;   file << getenv("FEMUS_DIR")<<"/USER_APPL/"
  <<getenv("FM_MYAPP")<<"/DATA/FSIAproperties.in";
  std::ifstream fin;  fin.open(file.str().c_str()); // stream file
#ifdef PRINT_INFO
  if(fin.is_open()) {  std::cout << "\nInit Reading = " << file.str() <<  std::endl; }
#endif
  //  reading param file -----------------------------------------------------------------
  if(fin.is_open()) {
    std::string string_value; std::string buf="";  // read double, string, dummy
    while(buf != "/") { fin >> buf; }  // find "/" file start
    fin >> buf;    while(buf != "/") {
      if(buf == "#") { getline(fin, buf); }  // comment line
      else { fin >> string_value; _FileMap[buf] = string_value; }
      fin >> buf;
    }
  } else { std::cerr << "FSIA_param.read_file(): no parameter file found" << std::endl;   abort(); }

// printing after reading ----------------------------------------------------------------------------------------
#ifdef PRINT_INFO
  std::cout << "\033[038;5;"<<"\n ADJOINT NAVIER-STOKES PARAMETER MAP \n \033[0m"  << std::endl;
  for(std::map<std::string, std::string >::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }
  std::cout << "\033[038;5;"<<FSA_F + 50<<";1m \
                \n----------------------------------------------\n\033[0m"  << std::endl;
#endif

  fin.close();
  return;
}

//
#endif  //ENDIF FSIA_EQUATIONS
// #endif  // NS_equation is personal

