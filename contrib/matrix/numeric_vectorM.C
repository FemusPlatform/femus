#include "numeric_vectorM.h"

#include "Solverlib_conf.h"

// Local Includes

#include "distributed_vectorM.h"
#include "laspack_vectorM.h"
#include "petsc_vectorM.h"
// #include "trilinos_epetra_vector.h"
// #include "shell_matrix.h"

// C++ includes
#include <cmath> // for std::abs
#include <memory>

//------------------------------------------------------------------
// NumericVectorM methods

// Full specialization for double datatypes

std::unique_ptr<NumericVectorM >
NumericVectorM::build(
  const ParallelM::Communicator &comm,
  const SolverPackageM solver_package){
  // Build the appropriate vector
  switch (solver_package){
#ifdef HAVE_LASPACKM
    case LASPACK_SOLVERSM:{
	std::unique_ptr<NumericVectorM > ap(new LaspackVectorM(comm, AUTOMATICM));
	return ap;
      }
#endif
#ifdef HAVE_PETSCM
    case PETSC_SOLVERSM:{
	std::unique_ptr<NumericVectorM > ap(new PetscVectorM(comm, AUTOMATICM));
	return ap;
      }
#endif
// #ifdef LIBMESH_HAVE_TRILINOS
//     case TRILINOS_SOLVERSM:{
// 	std::unique_ptr<NumericVectorM > ap(new EpetraVector<double>(comm, AUTOMATICM));
// 	return ap;
//       }
// #endif
    default:
//       std::unique_ptr<NumericVectorM > ap(NULL);
      std::cout << "NumericVectorM::build: Library not defined"; abort();
//       AutoPtr<NumericVectorM<T> > ap(new DistributedVector<T>);
//       return ap;
    } 
//   std::unique_ptr<NumericVectorM > ap(NULL);
  std::cout << "NumericVectorM::build: Library not defined"; abort();
//   return ap;    
}



// Full specialization for float datatypes (DistributedVector wants this)

// template <>
// int NumericVectorM<float>::compare (const NumericVectorM<float> &other_vector,
// 				   const double threshold) const
// {
//   assert (this->initialized());
//   assert (other_vector.initialized());
//   assert (this->first_local_index() == other_vector.first_local_index());
//   assert (this->last_local_index()  == other_vector.last_local_index());
// 
//   int rvalue     = -1;
//    int i = first_local_index();
// 
//   do
//     {
//       if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
// 	rvalue = i;
//       else
// 	i++;
//     }
//   while (rvalue==-1 && i<last_local_index());
// 
//   return rvalue;
// }

// Full specialization for double datatypes

int NumericVectorM::compare (const NumericVectorM &other_vector,
				    const double threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )first_different_i = i;
      else i++;
    }
  while (first_different_i==std::numeric_limits<int>::max() && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())   return -1;

  return first_different_i;
  
//   int rvalue     = -1;
//    int i = first_local_index();
// 
//   do
//     {
//       if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
// 	rvalue = i;
//       else
// 	i++;
//     }
//   while (rvalue==-1 && i<last_local_index());
// 
//   return rvalue;
}


// =================================
int NumericVectorM::local_relative_compare (
  const NumericVectorM &other_vector,			                     
  const double threshold) const
{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());


  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > threshold *
           std::max(std::abs((*this)(i)), std::abs(other_vector(i))))
	first_different_i = i;
      else	i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()   && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())    return -1;

  return first_different_i;
}


// ===================================================
int NumericVectorM::global_relative_compare (
  const NumericVectorM &other_vector,			             
  const double threshold
) const{
  assert (this->initialized());
  assert (other_vector.initialized());
  assert (this->first_local_index() == other_vector.first_local_index());
  assert (this->last_local_index()  == other_vector.last_local_index());
  

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  const double my_norm = this->linfty_norm();
  const double other_norm = other_vector.linfty_norm();
  const double abs_threshold = std::max(my_norm, other_norm) * threshold;

  do    {
      if ( std::abs( (*this)(i) - other_vector(i) ) > abs_threshold )
	first_different_i = i;
      else	i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()         && i<last_local_index());
  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);
  if (first_different_i == std::numeric_limits<int>::max())    return -1;

  return first_different_i;
}
// #ifdef TRIPLE_PRECISION
// // Full specialization for long double datatypes
// template <>
// int NumericVectorM<long double>::compare (const NumericVectorM<long double> &other_vector,
// 				         const double threshold) const
// {
//   assert (this->initialized());
//   assert (other_vector.initialized());
//   assert (this->first_local_index() == other_vector.first_local_index());
//   assert (this->last_local_index()  == other_vector.last_local_index());
// 
//   int rvalue     = -1;
//    int i = first_local_index();
// 
//   do
//     {
//       if ( std::abs( (*this)(i) - other_vector(i) ) > threshold )
// 	rvalue = i;
//       else
// 	i++;
//     }
//   while (rvalue==-1 && i<last_local_index());
// 
//   return rvalue;
// }
// #endif



// ---------------------------------------------------------
double NumericVectorM::subset_l1_norm (const std::set< int> & indices){
  NumericVectorM & v = *this;
  
  std::set< int>::iterator it = indices.begin();
  const std::set< int>::iterator it_end = indices.end();
  double norm = 0;
  for(; it!=it_end; ++it)    norm += std::abs(v(*it));
  this->comm().sum(norm);
//   ParallelM::sum(norm);
  return norm;
}

// --------------------------------------------------------
double NumericVectorM::subset_l2_norm (const std::set< int> & indices){
  NumericVectorM & v = *this;
  std::set< int>::iterator it = indices.begin();
  const std::set< int>::iterator it_end = indices.end();
  double norm = 0;
  for(; it!=it_end; ++it)    norm += (v(*it)*v(*it));
  this->comm().sum(norm);
//   ParallelM::sum(norm);
  return std::sqrt(norm);
}


double NumericVectorM::subset_linfty_norm (const std::set< int> & indices){
  NumericVectorM & v = *this;
  std::set< int>::iterator it = indices.begin();
  const std::set< int>::iterator it_end = indices.end();
  double norm = 0;
  for(; it!=it_end; ++it)    {
      double value = std::abs(v(*it));
      if(value > norm)     norm = value;
    }
//   ParallelM::max(norm); 
  this->comm().max(norm);

  return norm;
}



// template <typename T>
// void NumericVectorM<T>::add_vector (const NumericVectorM<T>& v,
// 				   const ShellMatrix<T>& a)
// {
//   a.vector_mult_add(*this,v);
// }



//------------------------------------------------------------------
// Explicit instantiations
// template class NumericVectorM<Number>;
