// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs
#include <cassert>
#include <limits> // std::numeric_limits<T>::min()

// Local Includes
#include "distributed_vectorM.h"
#include "dense_vectorM.h"
#include "dense_subvectorM.h"
#include "parallelM.h"
// #include "tensor_toolsM.h"


//--------------------------------------------------------------------------
// DistributedVectorM methods

// * =====================================
double DistributedVectorM::sum () const
{  // This function must be run on all processors at once
  parallel_object_onlyM();
//   parallel_onlyM();

  assert (this->initialized());assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_sum = 0.;
  for ( int i=0; i<local_size(); i++)  local_sum += _values[i];
  this->comm().sum(local_sum);
//   ParallelM::sum(local_sum);
  return local_sum;
}

// * ========================================
double DistributedVectorM::l1_norm () const
{  // This function must be run on all processors at once
//   parallel_onlyM();
  parallel_object_onlyM();
  assert (this->initialized());  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_l1 = 0.;
  for ( int i=0; i<local_size(); i++)  local_l1 += std::abs(_values[i]);
  this->comm().sum(local_l1);
//   ParallelM::sum(local_l1);
  return local_l1;
}

// * ==========================================
double DistributedVectorM::l2_norm () const
{  // This function must be run on all processors at once
//   parallel_onlyM();
  parallel_object_onlyM();
  assert (this->initialized());  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_l2 = 0.;
  for (numeric_index_type i=0; i<local_size(); i++)  local_l2 += (_values[i]*_values[i]);
  //   ParallelM::sum(local_l2);
  this->comm().sum(local_l2);
  return std::sqrt(local_l2);
}

// * =============================================
double DistributedVectorM::linfty_norm () const
{  // This function must be run on all processors at once
//   parallel_onlyM();
  parallel_object_onlyM();
  assert (this->initialized());  assert ((int)(int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  double local_linfty = 0.;
  for (numeric_index_type i=0; i<local_size(); i++)
    local_linfty  = std::max(local_linfty,
			     static_cast<double>(std::abs(_values[i]))
			     ); // Note we static_cast so that both
                                // types are the same, as required
                                // by std::max
  this->comm().max(local_linfty);
//   ParallelM::max(local_linfty);
  return local_linfty;
}

// ===============================================
NumericVectorM& DistributedVectorM::operator += (
  const NumericVectorM& v
){
  assert (this->closed());  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add(1., v);
  return *this;
}

// ============================================
NumericVectorM& DistributedVectorM::operator -= (
  const NumericVectorM& v
){
  assert (this->closed());  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add(-1., v);
  return *this;
}

// * ============================================
void DistributedVectorM::add_vector (
  const std::vector<double>& v,
  const std::vector<numeric_index_type>& dof_indices
){
  assert (!v.empty());  assert (v.size() == dof_indices.size());
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  for (std::size_t i=0; i<v.size(); i++)    add (dof_indices[i], v[i]);
}

// * ============================================
void DistributedVectorM::add_vector (
  const NumericVectorM& V,				      
  const std::vector<numeric_index_type>& dof_indices
){
  assert ((int)V.size() == (int)dof_indices.size());  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  for (int i=0; i<(int)V.size(); i++)  add(dof_indices[i], V(i));
}

// * ============================================
void DistributedVectorM::add_vector (
  const DenseVectorM& V,				    
  const std::vector<numeric_index_type>& dof_indices)
{
  assert ((int)V.size() == (int)dof_indices.size());     assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (std::size_t i=0; i<V.size(); i++)    add (dof_indices[i], V(i));
}
// * ============================================
void DistributedVectorM::reciprocal()
{
  for (numeric_index_type i=0; i<local_size(); i++){
      // Don't divide by zero
      assert(_values[i]== double(0));
      _values[i] = 1. / _values[i];
    }
}
// ============================================
void DistributedVectorM::add (const double v)
{
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( numeric_index_type i=0; i<local_size(); i++)    _values[i] += v;
}

// ============================================
void DistributedVectorM::add (const NumericVectorM& v)
{
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  add (1., v);
}

// ============================================
void DistributedVectorM::add (
  const double a, 
  const NumericVectorM& v
){
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  add(a, v);
}

// ============================================
void DistributedVectorM::insert (
  const std::vector<double>& v,				   
  const std::vector<numeric_index_type>& dof_indices
){
  assert (!v.empty());  assert (v.size() == dof_indices.size());
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);
  
  for (std::size_t  i=0; i<v.size(); i++)    this->set (dof_indices[i], v[i]);
}

// ============================================
void DistributedVectorM::insert (
  const NumericVectorM& V,				   
  const std::vector<numeric_index_type>& dof_indices
){
  assert ((int)V.size() == (int)dof_indices.size());
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (int  i=0; i<(int)V.size(); i++) this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVectorM::insert (
  const DenseVectorM& V,				   
  const std::vector< int>& dof_indices
){
  assert ((int)V.size() == (int)dof_indices.size());  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (int i=0; i<(int)V.size(); i++) this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVectorM::insert (
  const DenseSubVectorM& V,				   
  const std::vector< numeric_index_type>& dof_indices)
{
  assert ((int)V.size() == (int)dof_indices.size());  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (int i=0; i<(int)V.size(); i++)  this->set (dof_indices[i], V(i));
}

// ============================================
void DistributedVectorM::scale (
  const double factor
){
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (int i=0; i<(int)local_size(); i++)   _values[i] *= factor;  
}

// ============================================
void DistributedVectorM::abs()
{
  assert (this->initialized());
  assert ((_last_local_index - _first_local_index) == _local_size);

  for (int i=0; i<local_size(); i++) this->set(i,std::abs(_values[i]));
}

// * ============================================
double DistributedVectorM::dot (const NumericVectorM& V) const
{
  // This function must be run on all processors at once
//   parallel_onlyM();
  parallel_object_onlyM();

  // Make sure the NumericVector passed in is really a DistributedVectorM
  const DistributedVectorM* v = dynamic_cast<const DistributedVectorM*>(&V);

  // Make sure that the two vectors are distributed in the same way.
  assert ( this->first_local_index() == v->first_local_index() );
  assert ( this->last_local_index()  == v->last_local_index()  );
  
  // The result of dotting together the local parts of the vector.
  double local_dot = 0;

  for (int i=0; i<this->local_size(); i++)
    local_dot += this->_values[i] * v->_values[i];

  // The local dot products are now summed via MPI
   this->comm().sum(local_dot);
  //   ParallelM::sum(local_dot);
  return local_dot;
} 

// ============================================
NumericVectorM& 
DistributedVectorM::operator = (const double s)
{
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  for ( int i=0; i<local_size(); i++) _values[i] = s;
  return *this;
}

// ============================================
NumericVectorM&
DistributedVectorM::operator = (const NumericVectorM& v_in)
{
  // Make sure the NumericVector passed in is really a DistributedVectorM
  const DistributedVectorM* v = static_cast<const DistributedVectorM*>(&v_in);
  *this = *v;
  return *this;
}

// ============================================
DistributedVectorM&
DistributedVectorM::operator = (const DistributedVectorM& v)
{
  this->_is_initialized    = v._is_initialized;
  this->_is_closed         = v._is_closed;

  _global_size       = v._global_size;
  _local_size        = v._local_size;
  _first_local_index = v._first_local_index;
  _last_local_index  = v._last_local_index;
  
  if (v.local_size() == this->local_size()) { _values = v._values;    }
  else  { abort();  }
  
  return *this;
}

// ============================================
NumericVectorM&
DistributedVectorM::operator = (
  const std::vector<double>& v
){
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  if ((int)v.size() == local_size())    _values = v;

  else if ((int)v.size() == size())
    for ( int i=first_local_index(); i<last_local_index(); i++)
      _values[i-first_local_index()] = v[i];
  else    {      abort();    }
 
  return *this;
}

// ============================================
void DistributedVectorM::localize (
  NumericVectorM& v_local_in
) const{
  assert (this->initialized());  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  DistributedVectorM* v_local = dynamic_cast<DistributedVectorM*>(&v_local_in);
  v_local->_first_local_index = 0;
  v_local->_global_size = v_local->_local_size = v_local->_last_local_index = size();
  v_local->_is_initialized =    v_local->_is_closed = true;
  
  // Call localize on the vector's values.  This will help
  // prevent code duplication
  localize (v_local->_values);    
  
#ifndef LIBMESH_HAVE_MPI
  assert (local_size() == size());
#endif
}

// ============================================
void DistributedVectorM::localize (
  NumericVectorM& v_local_in,				     
  const std::vector<numeric_index_type>&) const
{
  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  // TODO: We don't yet support the send list; this is inefficient:
  localize (v_local_in);
}

// ============================================
void DistributedVectorM::localize (
  const  numeric_index_type first_local_idx,
  const numeric_index_type last_local_idx,				     
  const std::vector<numeric_index_type>& send_list
){
  // Only good for serial vectors
  assert (this->size() == this->local_size());
  assert (last_local_idx > first_local_idx);
  assert ((int)send_list.size() <= this->size());
  assert (last_local_idx < this->size());
  
  const numeric_index_type size       = this->size();
  const numeric_index_type local_size = (last_local_idx - first_local_idx + 1);

    // Don't bother for serial cases
  if ((first_local_idx == 0) &&      (local_size == size))    return;
  
  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  
  DistributedVectorM parallel_vec(this->comm());

  parallel_vec.init (size, local_size, true, PARALLELM);

  // Copy part of *this into the parallel_vec
  for (numeric_index_type i=first_local_idx; i<=last_local_idx; i++)
    parallel_vec._values[i-first_local_idx] = _values[i];

  // localize like normal
  parallel_vec.localize (*this, send_list);    
}

// ============================================
void DistributedVectorM::localize (std::vector<double>& v_local) const
{
  // This function must be run on all processors at once
//   parallel_onlyM();
  parallel_object_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local = this->_values;
  this->comm().allgather (v_local);
//   ParallelM::allgather (v_local);

#ifndef HAVE_MPI
  assert (local_size() == size());
#endif  
}

// ============================================
void DistributedVectorM::localize_to_one (
  std::vector<double>& v_local,					    
  const  processor_id_type pid) const
{
  // This function must be run on all processors at once
   parallel_object_onlyM();
//    parallel_onlyM();

  assert (this->initialized());
  assert ((int)_values.size() == _local_size);
  assert ((_last_local_index - _first_local_index) == _local_size);

  v_local = this->_values;
  this->comm().gather (pid, v_local);
  //   ParallelM::gather (pid, v_local);

#ifndef LIBMESH_HAVE_MPI
  assert (local_size() == size());
#endif  
}

// ============================================
void DistributedVectorM::pointwise_mult (
  const NumericVectorM&,					   
  const NumericVectorM&)
//void DistributedVectorM::pointwise_mult (const NumericVectorM& vec1,
//					   const NumericVectorM& vec2)
{
  std::cout << " Not implemented"; abort();
}

