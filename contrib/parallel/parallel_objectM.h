// The libMesh Finite Element Library.
// Copyright (C) 2002-2013 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef PARALLELM_OBJECT_H
#define PARALLELM_OBJECT_H

#include "Solverlib_conf.h"


// Local includes
#include "parallelM.h"
// mydefine --------------------------------------------------------
 #ifndef processor_id_typeM
#define processor_id_typeM int
 #endif
// Macro to identify and debug functions which should only be called in
// parallel on every processor at once
 #undef parallel_object_onlyM
 #ifndef NDEBUG
 #define parallel_object_onlyM() parallel_onlyM(this->comm())
 #else
 #define parallel_object_onlyM()  ((void) 0)
 #endif

// {
// ======================================================================================
// This class forms the base class for all other classes
// that are expected to be implemented in paralel. Each
// \p ParalelObject *requires* a \p ParallelM::Communicator object
// for construction.
// \author Benjamin S. Kirk \date 2013 modified \date 2016 S.Manservisi 

class ParallelObjectM{
    
protected:
  const ParallelM::Communicator &_communicator;
  
public:
  // ==================================================================================== 
  // Constructor. Requires a reference to the communicator
  // that defines the object's parallel decomposition.
  ParallelObjectM (
      const ParallelM::Communicator &comm_in
  ) :
    _communicator(comm_in)
  {}
  
  // ==================================================================================== 
  // Copy Constructor.
  ParallelObjectM (
      const ParallelObjectM &other
  ) :
    _communicator(other._communicator)
  {}
  
  // ==================================================================================== 
  // "Assignment" operator.  Simply asserts our references
  // are identical because this is the only thing that makes sense
  ParallelObjectM & operator= (
      const ParallelObjectM & other
  ){
    assert(&_communicator== &other._communicator);
    return *this;
  }

  // ==================================================================================== 
  // Destructor.  Virtual because we are a base class.
  virtual ~ParallelObjectM () {}
  // ==================================================================================== 
  // @returns a reference to the \p ParallelM::Communicator object used by this mesh.
  const ParallelM::Communicator & comm () const{ return _communicator; }
  // ==================================================================================== 
  // @returns the number of processors in the group.
  processor_id_typeM n_processors()const{return static_cast<processor_id_typeM>(_communicator.size());}
  // ==================================================================================== 
  // @returns the rank of this processor in the group.
  processor_id_typeM processor_id () const {return static_cast<processor_id_typeM>(_communicator.rank()); }


};
// } // namespace libMesh

#endif // LIBMESH_PARALLEL_OBJECT_H
