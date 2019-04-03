#include "Tparameters.h"

#include "Equations_conf.h"

#ifdef T_EQUATIONS

#include "User_T.h"
#include "MGUtils.h"
#include <fstream>

T_param::T_param()   {// SETTING DEFAULT VALUES
        _SolverType                  =GMRESM;
        _SolveSteady                 =0;
        _Supg                        =1;
        _Prt                         =0.85;
        _InterpolatedAlphaTurb       =0;
        _NumRestartSol               =1;
        _TimeDer                     =1;
        _BoundMap["Twall0"]          =   Twall0    ;
        _BoundMap["Twall"]           =   Twall     ;
        _BoundMap["insulation"]      =   insulation;
        _BoundMap["heat_flux"]       =   heat_flux ;
        _BoundMap["robinT"]          =   robinT    ;
        _BoundMap["ext_conv"]        =   ext_conv  ;
        _BoundMap["simmetry"]        =   simmetry  ;
    }// END T_param CONSTRUCTOR

T_param::~T_param(){
        _BoundMap.clear();
        _FileMap.clear();
        _map_Tgroup.clear();
        _BoundaryGroupsIDs.clear();
    }
    
    
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
    
    std:cout << " solve "<< _FileMap["SolveSteady"];
    
    if ( _FileMap ["SolveSteady"]!="" )   _SolveSteady  = stoi ( _FileMap["SolveSteady"] );
    if ( _FileMap ["Supg"]!="" )          _Supg= stoi ( _FileMap["Supg"] );
    if ( _FileMap ["Upwind"]!="" )        _Upwind=stof ( _FileMap["Upwind"] );
    if ( _FileMap ["Prt"]!="" )           _Prt=stof ( _FileMap["Prt"] );
    if ( _FileMap ["TimeDer"]!="" )       _TimeDer = stoi(_FileMap["TimeDer"]);
    
    if ( _FileMap ["NumRestartSol"]!="" )         _NumRestartSol=stoi ( _FileMap["NumRestartSol"] );
    if ( _FileMap ["InterpolatedAlphaTurb"]!="" ) _InterpolatedAlphaTurb=stof ( _FileMap["InterpolatedAlphaTurb"] );
    if ( mgutils._SolverTypeMap[_FileMap["SolverType"]] != 0 )  _SolverType = mgutils._SolverTypeMap[_FileMap["SolverType"]];
    
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

#endif
