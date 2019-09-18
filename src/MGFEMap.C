#include "MGFE.h"
#include "MGFEMap.h"

void MGFEMap::set_FE(MGFE* mgfe) {
  FE_map.insert(make_pair(mgfe->_order, mgfe));
  return;
}  ///< Insert in the map a FEM
