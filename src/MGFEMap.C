#include "MGFE.h"
#include "MGFEMap.h"

void MGFEMap::set_FE(MGFE* mgfe) {
  FE_map.insert(make_pair(mgfe->_order, mgfe));
  return;
}  ///< Insert in the map a FEM

// ==============================
/// This function destroys the MGFE class
MGFEMap::~MGFEMap() {  // ==============================
  for (std::map<int, MGFE*>::iterator itr = FE_map.begin(); itr != FE_map.end(); itr++) {
    (itr->second)->~MGFE();
  }
  FE_map.clear();
#if PRINT_INFO == 2
  std::cout << "~MGFEMap() called" << std::endl;
#endif
}
