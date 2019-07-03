#include "Equations_conf.h"

// ============================================
#ifdef DS_EQUATIONS // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverDS.h"
#include "MGSclass_conf.h"

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
/// This function constructs the 3d-2D MGSolDS class
// ==========================================================================
/*! This constructor needs    MGEquationsSystem &mg_equations_map_in object to be constructed.
* This equation has 1 quadratic variable (dx) defined in nvars_in[]=(0,0,1),
* equation name "dx", basic variable name "dx"
*/
MGSolDS::MGSolDS(
    MGEquationsSystem &mg_equations_map_in, ///<  mg_equations_map_in pointer
    const int nvars_in[],                   ///< KLQ number of variables
    std::string eqname_in,                  ///< equation name
    std::string varname_in                  ///< basic variable name
):
    MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
    _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
    _dt(stod(_mgutils._sim_config["dt"])),                        // parameter  dt
    _uref(_mgutils._mat_prop["Uref"]),         // parameter  u reference
    _lref(_mgutils._mat_prop["Lref"]),         // parameter  l reference
    _rhof(_mgutils._mat_prop["rho0"]),         // parameter density
    _muf(_mgutils._mat_prop["mu0"]) ,          // parameter viscosity
    _Tref(_mgutils._mat_prop["Tref"]),         // parameter  temperature reference
    _rhos(mg_equations_map_in.get_par("rhos")),  // parameter solid density
    _disp_d(mg_equations_map_in.get_par("disp_d"))
{   //  =========================================================================

    // READ PARAMETERS FROM CLASS DS_parameter
    _DS_parameter.read_param(_mgutils);
    _NumRestartSol = 2;

    /// A) reading parameters  for field coupling (in _FF_idx[])
    _nDSdim=DIMENSION;
    _dir=0;
    if(!varname_in.compare("dx")) {
        _dir=0;
    }
    if(!varname_in.compare("dy")) {
        _dir=1;
    }
    if(!varname_in.compare("dz")) {
        _dir=2;
    }
#if DIMENSION ==2
    if(_dir==2) {
        std::cout<<"Too many Dimension!!\n";
    }
#endif
    _var_names[0]=varname_in;
    _refvalue[0]=_lref;


    for(int k_index=0; k_index<30; k_index++) {
        _FF_idx[k_index]=-1;
    }
    /// B) setting class variable name T (in _var_names[0]) and ref value T_ref (in _refvalue[0])
    _var_names[0]=varname_in;
    _refvalue[0]=_Tref;

    /// C ) Setting the  solver type (with _solver[l]->set_solver_type(SOLVERT))
    for(int l=0; l<_NoLevels; l++) {
        _solver[l]->set_solver_type(_DS_parameter._SolverType);
    }

    /// D) Setting nondimensional parameters _alpha _IPrdl _IRe .....

    _SolveDS=(_mgutils._sim_config["SolveFluidStructure"].compare("yes") == 0) ? true: false;
    return;
}



//  ===============================================================================================
/// This function assembles the matrix and the rhs:
//  ===============================================================================================
void  MGSolDS::GenMatRhs(
    const double    /**< time (in) */,
    const int Level /**< discretization Level (in) */,
    const int mode  /**< y/n assemble rhs  (1=rhs) (0=only matrix) (in)*/
) {  // ===============================================

    /// a) Set up
    const int unsteady_flag=stoi(_mgutils._sim_config["SolveSteady"]);
    const int iaxisim=(int) (_mgutils._geometry["Axysim"]);
    int fl_int;
    int phase;

    // geometry ---------------------------------------------------------------------------------------
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];   // mesh nodes
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];    // element sides
    int        el_conn[NDOF_FEM];   // element connectivity
    int        el_neigh[NDOF_FEM];   // bd element connectivity
    int        sur_toply[NDOF_FEMB];   // boundary topology
    int        flag_group[NDOF_FEM];
    int        flag_sur_group[NDOF_FEMB];
    double x_m[DIMENSION];
    double normal[DIMENSION];
    // gauss integration  -----------------------------------------------------------------------------

    double u_old[NDOF_FEM];
    double  p_old;
    double l_old[DIMENSION*NDOF_FEM];
    const int el_ngauss = _fe[2]->_NoGauss1[ _nDSdim-1];                   // elem gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[ _nDSdim-2];             // bd elem gauss points

    // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
    int el_ndof[3];
    el_ndof[0]=1;
    int elb_ndof[3];
    elb_ndof[0]=1;   // number of el dofs
    int el_mat_nrows =0;                                               // number of mat rows (dofs)
    for(int ideg=1; ideg<3; ideg++) {
        el_ndof[ideg]=_fe[ideg]->_NoShape[_nDSdim-1];
        elb_ndof[ideg]=_fe[ideg]->_NoShape[_nDSdim-2];
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    const int el_ndof2=_fe[2]->_NoShape[_nDSdim-1];

    int el_mat_ncols = el_mat_nrows;   // square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector
    MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
    // coupling  fields -------------------------------------------------------------------------------
    for(int k=0; k<30; k++) {// coupling  basic system fields
        const int idx= _data_eq[2].tab_eqs[k];
        _FF_idx[k]=(idx>=0)?_data_eq[2].indx_ub[idx]:-1;
    }
    double vel_g[DIMENSION];
    for(int idim=0; idim< _nDSdim; idim++) {
        vel_g[idim] =0.;    // velocity not coupled
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
    A[Level]->zero();
    if(mode ==1) {
        b[Level]->zero();    // global matrix+rhs
    }
    DenseMatrixM KeM;
    DenseVectorM FeM;                // local  matrix+rhs
    KeM.resize(el_mat_nrows,el_mat_ncols);
    FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs

    int ndof_lev=0;
    for(int pr=0; pr <_mgmesh._iproc; pr++) {
        int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
        ndof_lev +=delta;
    }

    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for(int iel=0; iel < (nel_e - nel_b); iel++) {

        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        /// 1. geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (_xx_qnds)
        _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,_xx_qnds);
        _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
        // fill the node data vectors
        for(int deg=0; deg<3; deg++) {
            for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                                     el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
            }
        }

        // ----------------------------------------------------------------------------------
        /// 2. Boundary integration  (bc)
        // ----------------------------------------------------------------------------------

        for(int idim=0; idim< _nDSdim; idim++) {
            x_m[idim]=0.;
            for(int idofk=0; idofk< el_ndof[2]; idofk++) {
                int d=idim*NDOF_FEM+idofk;
                x_m[idim] +=_xx_qnds[d]/NDOF_FEM;
                flag_group[idofk]=ext_mesh->_bc_id[el_conn[idofk]];
            }
        }
        //  external cell properties -------------------------------------
        for(int d=0; d< el_mat_nrows; d++) {
            _bc_el[d]=1;
            l_old[d]= _data_eq[2].ub[(_FF_idx[SDSX_F]+_dir)*NDOF_FEM+d];
        }

        phase=(ext_mesh->_mat_id[iel+nel_b]==2)?0:1;

        for(int iside=0; iside< el_sides; iside++)  {
            if(el_neigh[iside] == -1) {
                for(int idof=0; idof<elb_ndof[2]; idof++) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    int idofb=sur_toply[idof];                                      // connectivity vector
                    for(int idim=0; idim< _nDSdim; idim++) {
                        _xxb_qnds[idim*NDOF_FEMB+idof]=_xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                    }
                }

                int sign_normal=1;
                _fe[2]->normal_g(_xxb_qnds,x_m,normal,sign_normal);
                int dir_maxnormal = (fabs(normal[0])>fabs(normal[1]))?0:1 ;
                dir_maxnormal= (fabs(normal[dir_maxnormal])>fabs(normal[DIMENSION-1]))? dir_maxnormal:DIMENSION-1;
                set_bc_matrix(KeM,FeM,dir_maxnormal,sur_toply,el_ndof,elb_ndof, elb_ngauss,normal,0); //  boundary conditions
            }
        }
        // ----------------------------------------------------------------------------------
        //   3. Volume integration
        // ----------------------------------------------------------------------------------
        if(_FF_idx[NS_F]>=0) {
            for(int idim=0; idim< _nDSdim; idim++) {
                for(int d=0; d< NDOF_FEM; d++)  vel_g[idim] += _data_eq[2].ub[NDOF_FEM*(_FF_idx[NS_F]+idim)+d]/NDOF_FEM;
            }
        } else {
            vel_g[ _nDSdim-1] = 0.;
            vel_g[0] = 1.;
            vel_g[1] = 1.;
        }

        // volume integral
        if(phase==0) {
            matrixrhsvol_liq_ds(KeM,FeM, el_ndof,flag_group);
        }
        else { // ==========================  solid ===================================================
            matrixrhsvol_sol_ds(KeM,FeM, el_ndof);
        }

        // ----------------------------------------------------------------------------------
        //  4. add local to global
        // ----------------------------------------------------------------------------------

        A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
        if(mode == 1) {
            b[Level]->add_vector(FeM,el_dof_indices);    // global rhs
        }

    } // end of element loop

    /// 5. clean
    el_dof_indices.clear();
    A[Level]->close();
    if(mode == 1) {
        b[Level]->close();
    }
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(DS)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolDS::MGTimeStep(
    const double time,  ///< time
    const int /*iter*/  ///< Number of max inter
) {
// =========================================================================================
    if(_SolveDS) {
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
        GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
        for(int Level = 0 ; Level < _NoLevels-1; Level++) {
            GenMatRhs(time,Level,0);    // matrix
        }
#if PRINT_TIME==1
        std::clock_t end_time=std::clock();
        std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

        /// C) Solution of the linear MGsystem (MGSolDS::MGSolve).
        if(_mgutils.get_name() != 1) {
            MGSolve(1.e-6,40);
        }
#if PRINT_TIME==1
        end_time=std::clock();
        std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
                  << "s "<< std::endl;
#endif

        /// D) Update of the old solution at the top Level  (MGSolDS::OldSol_update),
        _x_olds[_NoLevels-1][0]->localize(*_x_olds[_NoLevels-1][1]);
        x[_NoLevels-1]->localize(*_x_olds[_NoLevels-1][0]);
        const int flag_moving_mesh = _mgutils._geometry["moving_mesh"];
        if (flag_moving_mesh==1)  MoveMesh(_NoLevels-1);
    }
    return;
}// =======================================================================================

//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolDS::MoveMesh(
    const int Level  // Level <-
) {  // ===============================================
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
    const ParallelM::Communicator &comm1=_mgmesh._comm.comm();
    MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
    // print linear -----------------------------------
    double *sol_c=new double[_mgmesh._GeomEl.n_l[0]];
    int ndof_lev=0;
    disp[Level]->close();
    disp[Level]->zero();
    for(int iproc=0; iproc<_mgmesh._n_subdom; iproc++) {

        int n_iel0=_mgmesh._off_el[0][_NoLevels*iproc+0];
        int n_ielf=_mgmesh._off_el[0][0+_NoLevels*iproc+1];
        for(int iel=0; iel <n_ielf-n_iel0; iel++) {
            int e_indx=(iel+n_iel0)*NDOF_FEM;
            if(ext_mesh->_mat_id[iel+n_iel0]==4) {   //solid disp
                const int ivar =0;
                for(int in=0; in<NDOF_FEM; in++) {
                    int gl_i=_mgmesh._el_map[0][ e_indx+in];
                    double val= (*_x_olds[Level][0])(_node_dof[Level][gl_i+ivar*offset]);
                    disp[Level]->set(_node_dof[Level][ _mgmesh._el_map[0][ e_indx+in]],val);
                }
            } //end solid if
            else { //liquid at coarse level
                const int ivar =0;
                for(int in=0; in<NDOF_P; in++) {
                    int gl_i=_mgmesh._el_map[0][ e_indx+in];
                    double val= (*_x_olds[Level][0])(_node_dof[Level][gl_i+ivar*offset]);
                    sol_c[in]= val;
                }
                for(int in=0; in<NDOF_FEM; in++) {    // mid-points
                    double sum=0.;
                    for(int jn=0; jn<NDOF_P; jn++) {
                        sum += _mgmesh._GeomEl.Prol[in*NDOF_P+jn]*sol_c[jn];
                    }
                    disp[Level]->set(_node_dof[Level][ _mgmesh._el_map[0][ e_indx+in]],sum);
                }
            }
        }
        for (int ilev=1; ilev<=Level; ilev++) { //finer levels
            int n_iel0=_mgmesh._off_el[0][_NoLevels*iproc+ilev];
            int n_ielf=_mgmesh._off_el[0][ilev+_NoLevels*iproc+1];
            for(int iel=0; iel <n_ielf-n_iel0; iel++) {
                int e_indx=(iel+n_iel0)*NDOF_FEM;
                if(ext_mesh->_mat_id[iel+n_iel0]==4) {   //solid disp
                    const int ivar =0;
                    for(int in=0; in<NDOF_FEM; in++) {
                        int gl_i=_mgmesh._el_map[0][ e_indx+in];
                        double val= (*_x_olds[Level][0])(_node_dof[Level][gl_i+ivar*offset]);
                        disp[Level]->set(_node_dof[Level][ _mgmesh._el_map[0][ e_indx+in]],val);
                    }
                } //end solid if
                else {  //liquid
                    const int ivar =0;
                    for(int in=0; in<NDOF_P; in++) {
                        int gl_i=_mgmesh._el_map[0][ e_indx+in];
                        double val= (*disp[Level])(_node_dof[Level][gl_i+ivar*offset]);
                        sol_c[in]= val;
                    }
                    for(int in=0; in<NDOF_FEM; in++) {    // mid-points
                        double sum=0.;
                        for(int jn=0; jn<NDOF_P; jn++) {
                            sum += _mgmesh._GeomEl.Prol[in*NDOF_P+jn]*sol_c[jn];
                        }
                        disp[Level]->set(_node_dof[Level][ _mgmesh._el_map[0][ e_indx+in]],sum);
                    }
                } //end liquid
            } // ---- end iel -------------------------
        } // 2bB end interpolation over the fine mesh ------------------------
    }
    disp[Level]->close();
    delete sol_c;
    return;
}

// ===================================================
void DS_param::read_param(
    MGUtils &mgutils)
{   // ==================================================

    read_file();  // Reading parameters from file DSproperties.in

    // Boundary condition block ------------------------------------------------------------
    std::cout<<"Displacement boundary condition block \n";
    std::string GroupString = mgutils._sim_config["BoundaryGroups"]; // from Simulation configuration (es "20,21,13")
    int Length = GroupString.length(); // string length from Simulation configuration (es "20,21,13")
    int count=0;
    int  pos1=0;
    std::string temps;
    while(count<Length) {
        if(GroupString.at(count)==',') {
            temps = GroupString.substr(pos1,count-pos1);
            _BoundaryGroupsIDs.push_back(stoi(temps));
            pos1=count+1;
        }
        count++;
    }
    _BoundaryGroupsIDs.push_back(stod(GroupString.substr(pos1,Length-pos1)));
    const int numero = _BoundaryGroupsIDs.size();

    std::string BDcond;
    for(int i=0; i<_BoundaryGroupsIDs.size(); i++) {
        BDcond = "DSgroup"+to_string(_BoundaryGroupsIDs[i]);
        _map_DSgroup[_BoundaryGroupsIDs[i]]=_BoundMap[_FileMap[BDcond]];
    }
    // end boundary condition block --------------------------------------------------------

    // SETTING OTHER PARAMETERS ------------------------------------------------------------
    std::cout<<" DISPLACEMENT PARAMETER map (_DS_parameter) in  DSproperties.in + UserDS.h :\n";
//     std:cout << " solve"<< _FileMap["SolveSteady"];
    if(_FileMap ["SolveSteady"]!="") {
        _SolveSteady  = stoi(_FileMap["SolveSteady"]);
    }
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter._SolveSteady (in UserDS.h) \n";
    }
    std::cout <<" DSproperties.in: default value for _DS_parameter._SolverType  (in UserDS.h) \n";
    if(_FileMap ["DynamicUnderRelaxation"]!="") {
        _UnderRelaxation =stof(_FileMap["DynamicUnderRelaxation"]);
    }
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter._DynamicUnderRelaxation (in UserDS.h)\n";
    }
    if(_FileMap ["Supg"]!="")   _Supg= stoi(_FileMap["Supg"]);
    else {
        std::cout <<" DSproperties.in:  default value for _DS_parameter._Supg (in UserDS.h)\n";
    }
    if(_FileMap ["Upwind"]!="")  _Upwind=stof(_FileMap["Upwind"]);
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter._Upwind (in UserDS.h)\n";
    }
    if(_FileMap ["ReactionNumberBased"]!="") {
        _ReactionNumberBased     = stoi(_FileMap["ReactionNumberBased"]);
    }
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter._ReactionNumberBased  (in UserDS.h)\n";
    }

    if(_FileMap ["FlatProfile"]!="")  _FlatProfile = stoi(_FileMap["FlatProfile"]);
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter._FlatProfile (in UserDS.h)\n";
    }
    if(_FileMap ["Prt"]!="") {
        _Prt=stof(_FileMap["Prt"]);
    }
    else {
        std::cout <<" DSproperties.in: default value for _DS_parameter.Prt (in UserDS.h)\n";
    }

    _FileMap.clear();
    return;
}// END read_param FUNCTION ==============================================================


// ============================================================
/// This function reads Tparameter file (Tparameter.in)
void DS_param::read_file(
) {//  ===========================================================

    //  getting file name --------------------------------------------------------------------
    std::ostringstream file;
    file << getenv("APP_PATH")<<"/DATA/DSproperties.in";
    std::ifstream fin;
    fin.open(file.str().c_str()); // stream file
#ifdef PRINT_INFO
    if(fin.is_open()) {
        std::cout << "\nInit Reading = " << file.str() <<  std::endl;
    }
#endif
    //  reading param file -----------------------------------------------------------------
    if(fin.is_open()) {
        std::string string_value;
        std::string buf="";  // read double, string, dummy
        while(buf != "/") {
            fin >> buf;    // find "/" file start
        }
        fin >> buf;
        while(buf != "/") {
            if(buf == "#") {
                getline(fin, buf);    // comment line
            }
            else {
                fin >> string_value;
                _FileMap[buf] = string_value;
            }
            fin >> buf;
        }
    } else {
        std::cerr << "DS_param.read_file(): no parameter file found" << std::endl;
        abort();
    }

// printing after reading ----------------------------------------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << "\033[038;5;"<<"\n  DISPLACEMENT PARAMETER MAP \n \033[0m"  << std::endl;
    for(std::map<std::string, std::string >::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it) {
        std::cout << it->first << " " << it->second << "\n";
    }
    std::cout << "\033[038;5;"<<SDS_F + 50<<";1m \
                \n----------------------------------------------\n\033[0m"  << std::endl;
#endif

    fin.close();
    return;
}

#endif
// #endif // personal application
