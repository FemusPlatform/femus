#ifndef __mgfemap_h__
#define __mgfemap_h__

#include <map>
using namespace std;



/// Class containing the map of the FEMs
class MGFE;

// ============================================================================
// ============================ MGFEMap =======================================
// ============================================================================
class MGFEMap {

  // ===============  data ====================================================
  private:
      map<int,MGFE*>         FE_map; ///< Map of the FEMs
       
  public:
    // =======================  Constructor ===================================
    ///@{ \name CONSTRUCTOR-DESTRUCTOR
    /// Constructor
    MGFEMap(){}; 
    /// Destructor
    ~MGFEMap(){FE_map.clear();}
    ///@}
    
  // ====================  functions ==========================================  
  ///@{ \name FE GET/SET
  void  set_FE(MGFE *mgfe); //{FE_map.insert(make_pair(mgfe->_order,mgfe));} ///< Insert in the map a FEM
  
  /// Get a FEM from the map
  inline MGFE* get_FE(const int  label) const {return FE_map.find(label)->second;}
  ///@}
  
  ///@{ \name ITERATORS FOR FE MAP
  typedef std::map<int, MGFE*>::const_iterator const_giterator;     ///< Iterator for the map
  const_giterator gc_begin() const { return FE_map.begin();} ///< Return begin of the map
  
  ///< Return end of the map
  const_giterator gc_end() const { return FE_map.end();}
  ///@}
};

#endif