#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverT.h"


// configuration files -----------
#include "Printinfo_conf.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif


// standard  library
#include <sstream>

// local include -------------------
#include "MGGeomEl.h"
// #include "MGMesh.h"
#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"
#include "MGSystem.h"
#include "MGFE.h"
#include "MGUtils.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"
#include "parallelM.h"



// ======================================================
/// This function constructs the 3d-2D MGSolT class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
* This equation has 1 quadratic variable (T) defined in nvars_in[]=(0,0,1),
* equation name "T", basic variable name "T"
*/
MGSolT::MGSolT (
    MGEquationsSystem &mg_equations_map_in, ///<  mg_equations_map_in pointer
    const int nvars_in[],                   ///< KLQ number of variables
    std::string eqname_in,                  ///< equation name
    std::string varname_in                  ///< basic variable name
) :
MGSolDA ( mg_equations_map_in,nvars_in,eqname_in,varname_in ),
        _offset ( _mgmesh._NoNodes[_NoLevels-1] ), // mesh nodes
        _dt ( stod ( _mgutils._sim_config["dt"] ) ),                  // parameter  dt
        _uref ( _mgutils._mat_prop["Uref"] ),      // parameter  u reference
        _lref ( _mgutils._mat_prop["Lref"] ),      // parameter  l reference
        _rhof ( _mgutils._mat_prop["rho0"] ),      // parameter density
        _muf ( _mgutils._mat_prop["mu0"] ) ,       // parameter viscosity
        _Tref ( _mgutils._mat_prop["Tref"] ), // parameter  temperature reference
        _cp0 ( _mgutils._mat_prop["cp0"] ),   // parameter  Cp reference
_kappa0 ( _mgutils._mat_prop["kappa0"] ) { // parameter  conductivity reference
    //  =========================================================================

    // READ PARAMETERS FROM CLASS T_parameter
    _T_parameter.read_param ( _mgutils );


    /// A) reading parameters  for field coupling (in _FF_idx[])
    _nTdim=DIMENSION;
    _euler_impl=1;
    for ( int k_index=0; k_index<30; k_index++ ) {
        _FF_idx[k_index]=-1;
    }
    /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
    _var_names[0]=varname_in;
    _refvalue[0]=_Tref;

    /// C ) Setting the  solver type (with _solver[l]->set_solver_type(SOLVERT))
    for ( int l=0; l<_NoLevels; l++ ) {
        _solver[l]->set_solver_type ( _T_parameter._SolverType );
    }

    /// D) Setting nondimensional parameters _alpha _IPrdl _IRe .....
    _alpha=_kappa0/ ( _rhof*_cp0 );
    _IPrdl=_rhof*_alpha/_muf;
    _IRe=_muf/ ( _rhof*_uref*_lref );

    //todo prt da definire in file di turbutils -> deve essere lo stesso anche per turbolenza termica
    _IPrdl_turb=1./_T_parameter._Prt;//1./PRT;
    _alpha_turb=0.;
    //   _Nusselt=(_h_conv*_lref)/_kappa0;

// // //   AGGIUNGERE A FILE DI INPUT!!!!!
    _qheat=_mgutils._mat_prop["qheat"]*_lref/ ( _rhof*_cp0*_Tref*_uref );
    _qs=_mgutils._mat_prop["qs"]/ ( _rhof*_cp0*_Tref*_uref );
    _Wall_dist =_mgutils._geometry["Wall_dist"];


    std::map<std::string, bool> YesNo;
    YesNo["yes"] = true;
    YesNo["no"]  = false;
    _SolveT = YesNo[_mgutils._sim_config["SolveTemperature"]];
    return;
}



//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void  MGSolT::GenMatRhs (
    const double    /**< time (in) */,
    const int Level /**< discretization Level (in) */,
    const int mode  /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
) {  // ===============================================
//   double Crank_Nicolson =1.;
    /// a) Set up
    const int unsteady_flag=_T_parameter._SolveSteady;
    const int iaxisim = ( int ) ( _mgutils._geometry["Axysim"] );

    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];   // mesh nodes
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];    // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];                       // bd element connectivity
    int        sur_toply[NDOF_FEMB];                     // boundary topology

    // gauss integration  -----------------------------------------------------------------------------
    double x_m[DIMENSION];
    double normal[DIMENSION];
    double     u_old[NDOF_FEM];
    const int el_ngauss = _fe[2]->_NoGauss1[ _nTdim-1];                   // elem gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[ _nTdim-2];             // bd elem gauss points

    double WallDist[NDOF_FEM], AlphaTurb[NDOF_FEM];
double src_value[4];
#ifdef  HAVE_MED
//  int imesh=atoi(_mgutils.get_file("MESHNUMBER").c_str());
 
  int interface_id=9;
  EquationSystemsExtendedM *ext_es=dynamic_cast<EquationSystemsExtendedM *>(&_mgeqnmap);
  MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
  std::vector<int> mat=ext_mesh->_mat_id;
  std::vector<int> bc_id=ext_mesh->_bc_id;
  double bc_value[DIMENSION];
  InterfaceFunctionM * source= ext_es->get_interface_fun(interface_id);
   int * map_mg ;
  int * map_med;
  int n_map;
//   if(imesh==1){
  map_mg  = source->get_map_mg();
  map_med = source->get_map_med();
  n_map = source->get_n();
//  }
  
#endif    
    // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
    int el_ndof[3];
    el_ndof[0]=1;
    int elb_ndof[3];
    elb_ndof[0]=1;   // number of el dofs
    int el_mat_nrows =0;                                               // number of mat rows (dofs)
    for ( int ideg=1; ideg<3; ideg++ ) {
        el_ndof[ideg]=_fe[ideg]->_NoShape[_nTdim-1];
        elb_ndof[ideg]=_fe[ideg]->_NoShape[_nTdim-2];
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    const int el_ndof2=_fe[2]->_NoShape[_nTdim-1];

    int el_mat_ncols = el_mat_nrows;   // square matrix
    std::vector<int> el_dof_indices ( el_mat_ncols );   // element dof vector

    // coupling  fields -------------------------------------------------------------------------------
    for ( int k=0; k<30; k++ ) { // coupling  basic system fields
        const int idx= _data_eq[2].tab_eqs[k];
        _FF_idx[k]= ( idx>=0 ) ?_data_eq[2].indx_ub[idx]:-1;
    }
    double vel_g[DIMENSION];
    for ( int idim=0; idim< _nTdim; idim++ ) {
        vel_g[idim] =0.;    // velocity not coupled
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
    A[Level]->zero();
    if ( mode ==1 ) {
        b[Level]->zero();    // global matrix+rhs
    }
    DenseMatrixM KeM;
    DenseVectorM FeM;                // local  matrix+rhs
    KeM.resize ( el_mat_nrows,el_mat_ncols );
    FeM.resize ( el_mat_nrows );     // resize  local  matrix+rhs

    int ndof_lev=0;
    for ( int pr=0; pr <_mgmesh._iproc; pr++ ) {
        int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
        ndof_lev +=delta;
    }

//     printf("MGSolverT: V mean = %f \n", _mgutils._TurbParameters->_vmid);
    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for ( int iel=0; iel < ( nel_e - nel_b ); iel++ ) {

        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        /// 1. geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
        _mgmesh.get_el_nod_conn ( 0,Level,iel,el_conn,_xx_qnds );
        _mgmesh.get_el_neighbor ( el_sides,0,Level,iel,el_neigh );

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc ( Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd );
        // fill the node data vectors
        for ( int deg=0; deg<3; deg++ ) {
            for ( int eq=0; eq<_data_eq[deg].n_eqs; eq++ ) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol ( 0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                                       el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub );
            }
        }

        // ----------------------------------------------------------------------------------
        /// 2. Boundary integration  (bc)
        // ----------------------------------------------------------------------------------

        for ( int k=0; k< el_ndof[2]; k++ ) {
            _bc_el[k]= ( _bc_vol[k]/10==0 ) ?0:1;   // boundary condition
            _bc_bd[k]=1;
        }
        for ( int idim=0; idim< _nTdim; idim++ ) {
            x_m[idim]=0.;
            for ( int idofk=0; idofk< el_ndof[2]; idofk++ )   x_m[idim] +=_xx_qnds[idim*NDOF_FEM+idofk]/NDOF_FEM;
        }

        for ( int iside=0; iside< el_sides; iside++ )  {
            if ( el_neigh[iside] == -1 ) {
                for ( int idof=0; idof<elb_ndof[2]; idof++ ) {                    
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    int idofb=sur_toply[idof];             // connectivity vector          
                    _bc_bd[idofb]=0;
                    for ( int idim=0; idim< _nTdim; idim++ ) {
                        _xxb_qnds[idim*NDOF_FEMB+idof]=_xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                    }
                }
                int sign_normal=1;
                _fe[2]->normal_g ( _xxb_qnds,x_m,normal,sign_normal );
                bc_set ( KeM,FeM,sur_toply,el_ndof[2],elb_ndof[2],elb_ngauss,sign_normal, iaxisim );
            }
        }
        // ----------------------------------------------------------------------------------
        //   3. Volume integration
        // ----------------------------------------------------------------------------------
            src_value[0]=0.;_qheat=0.;
#ifdef  HAVE_MED
//     if(imesh==1){
     if (mat[iel+nel_b]==2) {
      int i_mg =0; int i_med=-1; 
      for(int i_mg =0; i_mg<n_map; i_mg++){
        if (map_mg[i_mg] ==  iel+nel_b) 
        {i_med=i_mg; break;}
      }
      if(i_med >-1) {int elem_med= map_med[i_med];
      source->eval(elem_med, 1, src_value);
      }
    }
    _qheat=src_value[0];
//     }
#endif
        
        
        //  external cell properties -------------------------------------
        if ( _FF_idx[K_F]>=0 ) _y_dist=_mgmesh._dist[ iel+nel_b];
        // volume integral
        
//         _mgutils._TurbParameters->GetLevelElemNodeWallDist(iel,Level,WallDist);
//         _mgutils._TurbParameters->GetLevelElemAlphaTurb(iel,Level,AlphaTurb);
        vol_integral ( KeM,FeM,el_ndof2,el_ngauss,mode, WallDist, AlphaTurb );

        A[Level]->add_matrix ( KeM,el_dof_indices );               // global matrix
        if ( mode == 1 ) {
            b[Level]->add_vector ( FeM,el_dof_indices );    // global rhs
        }
    } // end of element loop

    /// 5. clean
    el_dof_indices.clear();
    A[Level]->close();
    if ( mode == 1 ) {
        b[Level]->close();
    }
    //   A[Level]->print(); b[Level]->print();
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif
    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT::MGTimeStep (
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
// =========================================================================================
    if ( _SolveT ) {
/// A) Set up the time step
        std::cout  << std::endl << "\033[038;5;"<<196<<";1m "
        << "--------------------------------------------------- \n\t"
        <<  _eqname.c_str()
        << " solution of problem " << _mgutils.get_name()
        << "\n ---------------------------------------------------\n  \033[0m";

        /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
        std::clock_t start_time=std::clock();
#endif
        GenMatRhs ( time,_NoLevels-1,1 );                                             // matrix and rhs
        for ( int Level = 0 ; Level < _NoLevels-1; Level++ ) {
            GenMatRhs ( time,Level,0 );    // matrix
        }
#if PRINT_TIME==1
        std::clock_t end_time=std::clock();
        std::cout << "  Assembly time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

        /// C) Solution of the linear MGsystem (MGSolT::MGSolve).
        if ( _mgutils.get_name() != 1 ) {
            MGSolve ( 1.e-6,40 );
        }
#if PRINT_TIME==1
        end_time=std::clock();
        std::cout << " Assembly+solution time -----> ="<< double ( end_time- start_time ) / CLOCKS_PER_SEC
        << "s "<< std::endl;
#endif

        /// D) Update of the old solution at the top Level  (MGSolT::OldSol_update),
        x[_NoLevels-1]->localize ( *x_old[_NoLevels-1] );
    }
    return;
}// =======================================================================================



// ===================================================
void T_param::read_param (
MGUtils &mgutils ) {
    // ==================================================

    read_file();  // Reading parameters from file Tproperties.in

    // Boundary condition block ------------------------------------------------------------
    std::cout<<"Temprature boundary condition block \n";
    std::string GroupString = mgutils._sim_config["BoundaryGroups"]; // from Simulation configuration (es "20,21,13")
    int Length = GroupString.length(); // string length from Simulation configuration (es "20,21,13")
    int count=0;
    int  pos1=0;
    std::string temps;
    while ( count<Length ) {
        if ( GroupString.at ( count ) ==',' ) {
            temps = GroupString.substr ( pos1,count-pos1 );
            _BoundaryGroupsIDs.push_back ( stoi ( temps ) );
            pos1=count+1;
        }
        count++;
    }
    _BoundaryGroupsIDs.push_back ( stod ( GroupString.substr ( pos1,Length-pos1 ) ) );
    const int numero = _BoundaryGroupsIDs.size();

    std::string BDcond;
    for ( int i=0; i<_BoundaryGroupsIDs.size(); i++ ) {
        BDcond = "Tgroup"+to_string ( _BoundaryGroupsIDs[i] );
        _map_Tgroup[_BoundaryGroupsIDs[i]]=_BoundMap[_FileMap[BDcond]];
    }
    // end boundary condition block --------------------------------------------------------

    // SETTING OTHER PARAMETERS ------------------------------------------------------------
    std::cout<<" TEMPERATURE PARAMETER map (_T_parameter) in  Tproperties.in + UserT.h :\n";
    std:cout << " solve"<< _FileMap["SolveSteady"];
    if ( _FileMap ["SolveSteady"]!="" ) {
        _SolveSteady  = stoi ( _FileMap["SolveSteady"] );
    } else {
        std::cout <<" Tproperties.in: default value for _T_parameter._SolveSteady (in UserT.h) \n";
    }
    std::cout <<" Tproperties.in: default value for _T_parameter._SolverType  (in UserT.h) \n";
    if ( _FileMap ["DynamicUnderRelaxation"]!="" ) {
        _UnderRelaxation =stof ( _FileMap["DynamicUnderRelaxation"] );
    } else {
        std::cout <<" Tproperties.in: default value for _T_parameter._DynamicUnderRelaxation (in UserT.h)\n";
    }
    if ( _FileMap ["Supg"]!="" )   _Supg= stoi ( _FileMap["Supg"] );
    else {
        std::cout <<" Tproperties.in:  default value for _T_parameter._Supg (in UserT.h)\n";
    }
    if ( _FileMap ["Upwind"]!="" )  _Upwind=stof ( _FileMap["Upwind"] );
    else {
        std::cout <<" Tproperties.in: default value for _T_parameter._Upwind (in UserT.h)\n";
    }
    if ( _FileMap ["ReactionNumberBased"]!="" ) {
        _ReactionNumberBased     = stoi ( _FileMap["ReactionNumberBased"] );
    } else {
        std::cout <<" Tproperties.in: default value for _T_parameter._ReactionNumberBased  (in UserT.h)\n";
    }

    if ( _FileMap ["FlatProfile"]!="" )  _FlatProfile = stoi ( _FileMap["FlatProfile"] );
    else {
        std::cout <<" Tproperties.in: default value for _T_parameter._FlatProfile (in UserT.h)\n";
    }
    if ( _FileMap ["Prt"]!="" ) {
        _Prt=stof ( _FileMap["Prt"] );
    } else {
        std::cout <<" Tproperties.in: default value for _T_parameter.Prt (in UserT.h)\n";
    }
    
    if ( _FileMap ["InterpolatedAlphaTurb"]!="" ) {
        _InterpolatedAlphaTurb=stof ( _FileMap["InterpolatedAlphaTurb"] );
    } else {
        std::cout <<" Tproperties.in: default value for _T_parameter._InterpolatedAlphaTUrb (in UserT.h)\n";
    }
    _FileMap.clear();
    return;
}// END read_param FUNCTION ==============================================================


// ============================================================
/// This function reads Tparameter file (Tparameter.in)
void T_param::read_file (
) {//  ===========================================================

    //  getting file name --------------------------------------------------------------------
    std::ostringstream file;
    file << getenv ( "APP_PATH" ) <<"/DATA/Tproperties.in";
    std::ifstream fin;
    fin.open ( file.str().c_str() ); // stream file
#ifdef PRINT_INFO
    if ( fin.is_open() ) {
        std::cout << "\nInit Reading = " << file.str() <<  std::endl;
    }
#endif
    //  reading param file -----------------------------------------------------------------
    if ( fin.is_open() ) {
        std::string string_value;
        std::string buf="";  // read double, string, dummy
        while ( buf != "/" ) {
            fin >> buf;    // find "/" file start
        }
        fin >> buf;
        while ( buf != "/" ) {
            if ( buf == "#" ) {
                getline ( fin, buf );    // comment line
            } else {
                fin >> string_value;
                _FileMap[buf] = string_value;
            }
            fin >> buf;
        }
    } else {
        std::cerr << "T_param.read_file(): no parameter file found" << std::endl;
        abort();
    }

// printing after reading ----------------------------------------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << "\033[038;5;"<<"\n  TEMPERATURE PARAMETER MAP \n \033[0m"  << std::endl;
    for ( std::map<std::string, std::string >::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it ) {
        std::cout << it->first << " " << it->second << "\n";
    }
    std::cout << "\033[038;5;"<<T_F + 50<<";1m \
                \n----------------------------------------------\n\033[0m"  << std::endl;
#endif

    fin.close();
    return;
}



// =====================================================
#ifdef TBK_EQUATIONS
// =====================================================
void  MGSolT::f_mu ( double val[] ) {

//   if(_kappa_g[0]< 1.e-20) { _kappa_g[0]= 1.e-20; }  // kappa
//   if(_kappa_g[1]< 1.e-20) { _kappa_g[1]= 1.e-20; }  // kappa
//   double tau_k=1.;  // turbulent time
//
// #if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
//   tau_k=_y_dist/sqrt(_kappa_g[0]);
// #endif  // end kappa -------------------------------------------------------
//
// #if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ----------
//   tau_k= CMU*_kappa_g[0]/_kappa_g[1];// tau=1/omega=CMU*kappa/epsilon
// #endif   // kappa-epsilon --------------------------------------------------
//
// #if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
//   tau_k=1./_kappa_g[1]; // tau=1/omega
//
// #endif  // -----------------------------------------------------------------
//   if(tau_k> MAX_TAU) { tau_k= MAX_TAU; }
// // turbulent viscosity
//   double Re_t= _kappa_g[0]*tau_k/_IRe;
//   if(Re_t > MU_TOP) { Re_t =MU_TOP; }
//   if(Re_t < MU_LOW) { Re_t =MU_LOW; }
//   _nut_ratio=Re_t;
//
// //     // Boundary corrections
//   double R_t=Re_t/CMU;                                             // turbulent Reynolds number
//   double R_eps = _y_dist*sqrt(_kappa_g[0]/sqrt(R_t))/_IRe;
// //    double R_eps =_y_dist/sqrt((_muf/_rhof)*sqrt((_muf/_rhof)/_kappa_g[1]));    // *sqrt(_kappa_g[0]/sqrt(R_t))/_IRe;
// //    double y_plus = _y_dist*0.547722557505*sqrt(_kappa_g[0])/_IRe;
//
//
//
// #ifdef LOWRE     /// B) Low-Reynolds  correction model
//   double f_mu =(1.-exp(-1.*R_eps/(14.)))*(1.-exp(-1.*R_eps/14.));
//   _nut_ratio *= f_mu * (1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));    // correction
// #endif
//
// #ifdef SST       ///C) SST k-w model
//
//   // F1 and F2 coeff calculation
//   double F1,F2;
//   double F1_first  = 500.*_IRe*tau_k/(_y_dist*_y_dist);
//   double F2_first= 2.*sqrt(_kappa_g[0])*tau_k/(BETASTAR*_y_dist);
//   if(F1_first > F2_first) { F2_first=F1_first; }
//   F2=tanh(F2_first*F2_first);
//
//   // nu_t calculation
// //       double alpha_star        = 1.;//(0.024+ Re_t/6.)/(1.+ Re_t/6.);
//   double alpha_starb = sqrt(_sP)*F2*tau_k/0.31;
//   if(alpha_starb < 1.) {alpha_starb=1.;}  /*printf("alpha star is %f and and %f and %f\n",alpha_starb, 1./tau_k, F2);*/
//   _nut_ratio /= alpha_starb;
//
//   if(_nut_ratio > MU_TOP) { _nut_ratio =MU_TOP; }
//   if(_nut_ratio < MU_LOW) { _nut_ratio =MU_LOW; }
// //       printf("mu_turb is %f \n",_mu_turb);
//
// #endif
//
//
//
// #ifdef TTBK_EQUATIONS
//   /// A) Energy  Turbulent viscosity
//   // Energy  Turbulent viscosity
// //     if(_kappaT_g[1]< 1.e-20) _kappaT_g[1]= 1.e-20;            // kappa
// //     if(_kappaT_g[0]< 1.e-20) _kappaT_g[0]= 1.e-20;           // kappa
//
// #if (TTBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ----------
//   double tauT_k= CMU*_kappaT_g[0]/_kappaT_g[1]; // tau=1/omega=CMU*kappa/epsilon
// #endif   // kappa-epsilon --------------------------------------------------
//
// #if (TTBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
//   double tauT_k=1./_kappaT_g[1]; // tau=1/omega
// #endif
//   if(tauT_k> MAX_TAU) { tauT_k= MAX_TAU; }  // tau=1/omega=CMU*kappa/epsilon
//   double rT=tauT_k/tau_k;
//
// #if (TTBK_EQUATIONS/2==1)     // kappa_theta-epsilon_theta
//   double f_alpha =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl))))*(1.-exp(-1.*R_eps/14.));
//   double f_d=exp(-1.*R_t*R_t/(200.*200.));
//   double f_asym =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl))))*(1.-exp(-1.*R_eps/14.));
//
//   double a_wall = sqrt(2.*rT)*1.3*_IPrdl/pow(R_t,0.75)*f_d*f_alpha;
//   double a_inter =2.*rT/(.3 +rT)*f_alpha*exp(-1.*R_t*R_t/(500.*500.));
//   double asymp = 0.9*f_asym;
//
//   double b_nu = f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));
//
//   _IPrdl_turb=1.11111*(a_inter+a_wall+asymp)/b_nu;
// //        _IPrdl_turb=1./(0.9+0.7*_IPrdl/(_nut_ratio)); // kays model
// //       _IPrdl_turb=1./2.3;                     // SED model
// //      std::cout <<  _y_dist<<  " " << 1./_IPrdl_turb << "\n ";
// #endif
//
// #if (TTBK_EQUATIONS/2==2)       // kappa_theta-omega_theta ------------------------------
//   double f_alpha =(1.-exp(-1.*R_eps/(15.*sqrt(_IPrdl))))*(1.-exp(-1.*R_eps/16.));
//   double f_d=exp(-1.*R_t*R_t/(200.*200.));
//
//   double a_wall = sqrt(2.*rT)*5.*_IPrdl/pow(R_t,0.75)*f_d;
//   double a_inter =2.*rT/(.3+rT)*exp(-1.*R_t*R_t/(500.*500.));
//   double asymp = 0.8;
//
// //       double b_nu = f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));
//
//   _IPrdl_turb=1.11111111111*(f_alpha*(a_inter+a_wall+asymp));
//
// //       _IPrdl_turb=1./2.4;
// //       _IPrdl_turb=1./(0.9+0.7*_IPrdl/(_nut_ratio)); //Kays correlation
// //     std::cout <<  _y_dist<<  " " << 1./_IPrdl_turb << "\n ";
// #endif
//
// #endif
    return;
}
// end ---   f_mu -------------------------
#endif // ----------  end TBK_EQUATIONS  

// =============================================================================================================
// double   MGSolT::heff(double u_nlg_av[]) {
// //  h_eff computation for SUPG, --------------------------------------------------------------------
//   double  h_eff=1.e-20; double vdothmax_old=1.e-20;
//   for(int lpt=0; lpt< NDOF_P; lpt++) {// -------->
//     double hcurr=0, vdothmax=0.;
//     for(int idim=0; idim<_nTdim; idim++) {// -------->
//       const double dist_idim=(_xx_qnds[idim*NDOF_FEM+(lpt+1)%NDOF_P]-_xx_qnds[idim*NDOF_FEM+lpt]);
//       hcurr += dist_idim*dist_idim;  vdothmax += abs(u_nlg_av[idim]*dist_idim);
//     }// <-------- idim for
//     if(vdothmax > vdothmax_old) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;} // max scalar product
//     else if(vdothmax > vdothmax_old-1.e-10 && h_eff>sqrt(hcurr)) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;}
//   } // <------------ lpt for
//   // end h_eff computation----------------------------------------------------------------------------
//   return h_eff;
// }

// =============================================================================================================
// double   MGSolT::heff(double u_nlg_av[]) {
// //  h_eff computation for SUPG, --------------------------------------------------------------------
//   double  h_eff=1.e-20; double vdothmax_old=1.e-20;
//   for(int lpt=0; lpt< NDOF_P; lpt++) {// -------->
//     double hcurr=0, vdothmax=0.;
//     for(int idim=0; idim<_nTdim; idim++) {// -------->
//       const double dist_idim=(_xx_qnds[idim*NDOF_FEM+(lpt+1)%NDOF_P]-_xx_qnds[idim*NDOF_FEM+lpt]);
//       hcurr += dist_idim*dist_idim;  vdothmax += abs(u_nlg_av[idim]*dist_idim);
//     }// <-------- idim for
//     if(vdothmax > vdothmax_old) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;} // max scalar product
//     else if(vdothmax > vdothmax_old-1.e-10 && h_eff>sqrt(hcurr)) {h_eff=sqrt(hcurr); vdothmax_old=vdothmax;}
//   } // <------------ lpt for
//   // end h_eff computation----------------------------------------------------------------------------
//   return h_eff;
// }
// ============================================================================================================


#endif
// #endif // personal application

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
