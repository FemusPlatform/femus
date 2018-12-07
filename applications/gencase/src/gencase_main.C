// std lib ------------------------------------------------
// #include <iostream>     // for input output

// Libmesh include ----------------------------------------
#include "libmesh.h"     //for libmeshinit

// local include ------------------------------------------
#include "MGUtils.h"
#include "MGGeomEl.h" 
#include "MGGenCase.h"
#include "MGGenCase_conf.h"
// ========================================================
// =================  GENCASE   ===========================
// ========================================================

int main(int argc, char** argv) {

  LibMeshInit init(argc, argv);   // Initialize libMesh

  // ======= MGUtils =====
  std::vector<MGUtils*> mg_utils;
  for(int mesh=0; mesh<NUM_MESH; mesh++) mg_utils.push_back(new MGUtils(mesh+1)); // MGUtils class: read and print parameters 

  // =======MGGeomEl =====
  MGGeomEl mg_geomel; // set the fem element

  // ========= MGGenCase =====
  for(int mesh=0; mesh<NUM_MESH; mesh++){
  // parameters
   const int restart     = stoi(mg_utils[mesh]->_sim_config["restart"]);
   const int libmesh_gen = (int) (mg_utils[mesh]->_geometry["libmesh_gen"]);
   const int mgops_gen   = (int) (mg_utils[mesh]->_geometry["mgops_gen"]);
  
  if (!restart) {
    if (mgops_gen)  {
      MGGenCase mg_gencase(init.comm(), *mg_utils[mesh],mg_geomel);  // case Constructor
      mg_gencase.GenCase();                             // main class function
    } 
    else {std::cout << "No GenCase because it was already run once" << std::endl;}
  } 
  else {std::cout << "No GenCase because of restart flag !=0 " << std::endl;return 0;}
  std::cout<<'\n' << "======= Mesh No. " << mesh+1 << " printed ========" << std::endl;
  }

  return 0;  
}
