#ifndef __dense_set_h__
#define __dense_set_h__

// #include "compare_types.h"
#include "dense_set.h"

// C++ includes
#include <vector>
#include <cassert>
#include <cmath>


// ------------------------------------------------------------
// DenseVector class definition
class DenseSet {
 
public:
  /// The actual data values, stored as a 1D array.
  int *_set;
  int _n;

  /// Constructor.  Creates a dense vector of dimension \p n.
  DenseSet(const unsigned int n=0){ _n=n;   _set = new int[_n]; }
  
  /// Destructor.  Does nothing.     
  ~DenseSet() {delete []_set;}
  
  
   ///  Resize the vector
  void resize (const unsigned int n){delete [] _set; _set = new int[n];_n=n; };
  /// Set every element in the vector to 0.
  void zero(){for(int i=0;i<_n;i++) _set[i]=0; }
  
};
#endif // #ifndef __dense_vector_h__

