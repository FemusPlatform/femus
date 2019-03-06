// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"
#if (NS_EQUATIONS%2==0)
// ======================================================================================
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSolverP.h"
#include <fstream>
#include "MGUtils.h"


// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for 
/// the pressure equation in the pressure projection method
// =================================================
void MGSolP::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  int iel,
  double u_value[]) {
double mod = 0.;
 u_value[0] = mod;
  
// ============================================================================
// ============================================================================
}

// ============================================================================
/// This function  defines the boundary conditions  for 
/// the pressure equation in the pressure projection method
void MGSolP::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],int bc_Neum[],int bc_flag[]) {
  
  bc_Neum[0] = vel_fix;
  bc_Neum[0] = _P_parameter._map_NSgroup[bc_gam];
  
  return;
} // end boundary conditions ==========================


// ============================================================
/// This function reads NS parameter file (NSparameter.in)
void P_param::read_file (
) { //  ===========================================================

    //  getting file name --------------------------------------------------------------------
    std::ostringstream file;
    file << getenv ( "APP_PATH" ) <<"/DATA/NSproperties.in";
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
            fin >> buf;     // find "/" file start
        }
        fin >> buf;
        while ( buf != "/" ) {
            if ( buf == "#" ) {
                getline ( fin, buf );     // comment line
            } else {
                fin >> string_value;
                _FileMap[buf] = string_value;
            }
            fin >> buf;
        }
    } else {
        std::cerr << "NS_param.read_file(): no parameter file found" << std::endl;
        abort();
    }

// printing after reading ----------------------------------------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << "\033[038;5;"<<"\n  NAVIER-STOKES PARAMETER MAP \n \033[0m"  << std::endl;
    for ( std::map<std::string, std::string >::const_iterator it = _FileMap.begin(); it != _FileMap.end(); ++it ) {
        std::cout << it->first << " " << it->second << "\n";
    }
    std::cout << "\033[038;5;"<<NS_F + 50<<";1m \
                \n----------------------------------------------\n\033[0m"  << std::endl;
#endif

    fin.close();
    return;
}

void P_param::read_param (
    MGUtils &mgutils,
    int Proc
) {
    //  Reading parameter  -> NSproperties.in ------------------------------------------------
    read_file();

    // Boundary condition block ------------------------------------------------------------
    std::cout<<"Navier-Stokes - Pressure Equation boundary condition block \n";
    // boundary group names  from
    std::string GroupString = mgutils._sim_config["BoundaryGroups"];
    int count=0;
    int pos1=0;
    int Length = GroupString.length();
    while ( count<Length ) {
        if ( GroupString.at ( count ) ==',' ) {
            std::string  temps = GroupString.substr ( pos1,count-pos1 );
            _BoundaryGroupsIDs.push_back ( stoi ( temps ) );
            pos1=count+1;
        }
        count++;
    }
    // boundary group values from
    _BoundaryGroupsIDs.push_back ( stod ( GroupString.substr ( pos1,Length-pos1 ) ) );
    for ( int i=0; i<_BoundaryGroupsIDs.size(); i++ ) {
        std::string BDcond = "NSgroup"+to_string ( _BoundaryGroupsIDs[i] );
        _map_NSgroup[_BoundaryGroupsIDs[i]]=_BoundMap[_FileMap[BDcond]];
    }

    if ( Proc==0 ) {// PRINTING ON MESSAGE.LOG
        std::cerr<<"=============================================================================\n";
        std::cerr<<" NAVIER-STOKES - PRESSURE EQUATION - BOUNDARY CONDITIONS \n";
        for ( auto& x: _map_NSgroup ) {
	  std::cerr<<"Group n. "<<x.first<<"  Pressure condition "<<x.second<<std::endl;
        }
    }
    
    if ( _FileMap ["TimeDisc"] != "" ) {
      _TimeDisc = stoi ( _FileMap["TimeDisc"] );
      }
  else {
      _TimeDisc = 1;
      }
    if ( _FileMap ["AssembleOnce"] != "" ) {
      _AssembleOnce = stoi ( _FileMap["AssembleOnce"] );
      }
  else {
      _AssembleOnce = 0;
      }
      
      
    _FileMap.clear();
    return;
}

#endif // ENDIF NS_EQUATIONS==0

