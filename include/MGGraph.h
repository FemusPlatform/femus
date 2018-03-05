#ifndef _mggraph_
#define _mggraph_

#include <vector>

typedef std::vector<int> Row;

/// Contains information for the matrix sparse pattern

class Graph : public std::vector<Row> {

  public:

  int _m;
  int _n; 
  int _ml;
  int _nl; 
  int _ndz;
  int _noz; 
  int _ml_start;
  
  Graph():_m(0),_n(0),_ml(0),_nl(0),_ndz(0),_noz(0),_ml_start(0){}

};

#endif