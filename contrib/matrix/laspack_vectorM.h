#ifndef __laspack_vectorM_h__
#define __laspack_vectorM_h__

#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM

#include "Paralleltype_enum.h"
#include "Typedefs_conf.h"
// -----------------------------------
// C++ includes 
#include <cstdio> // for std::sprintf
#include <cassert> // for assert
// -----------------------------------
// Local includes 
#include "numeric_vectorM.h" // numerics
#include <qvector.h>  // laspack
#include <operats.h> // laspack
// -----------------------------------
// Forward declarations 
class LaspackLinearSolver;
class SparseMatrixM;
// -----------------------------------

// ==================================================
// Laspack vector.  Provides a nice interface to the
// Laspack C-based data structures for serial vectors.
// ===================================================

class LaspackVectorM : public NumericVectorM  
{
   // =====================================
  // DATA
  // =====================================
public:
  /// Actual Laspack vector datatype
  QVector _vec;  
  
   // =====================================
  // Constructor /Destructor
  // =====================================
 public:
  /// Dummy-Constructor. Dimension=0
  explicit
  LaspackVectorM (const ParallelTypeM = AUTOMATICM);
  
  /// Constructor. Set dimension to \p n and initialize all elements with zero.
  explicit
  LaspackVectorM (const unsigned int n,
                 const ParallelTypeM = AUTOMATICM);
    
  /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
  LaspackVectorM (const unsigned int n,
		 const unsigned int n_local,
                 const ParallelTypeM = AUTOMATICM);
  
  /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
  LaspackVectorM (const unsigned int N,
		 const unsigned int n_local,
		 const std::vector<unsigned int>& ghost,
                 const ParallelTypeM = AUTOMATICM);
    
  /// Destructor, deallocates memory. 
  ~LaspackVectorM ();
  
  /// Call the assemble functions
  void close (); 
  /// @returns the \p LaspackVectorM to a pristine state.
  void clear ();
  
  
  /// Set all entries to zero. Equivalent to \p v = 0, 
  void zero ();    
  /// Creates a copy of this vector and returns it in an \p AutoPtr.
  std::unique_ptr<NumericVectorM > clone () const;
  
  /// Change the dimension of the vector to \p N. 
  //  The reserved memory for
  // this vector remains unchanged if possible, to make things faster, but
  // this may waste some memory, so take this in the back of your head.
  // However, if \p N==0 all memory is freed, i.e. if you want to resize
  // the vector and release the memory not needed, you have to first call
  // \p init(0) and then \p init(N). 
  // On \p fast==false, the vector is filled by zeros.
  void init (const unsigned int N,
	     const unsigned int n_local,
	     const bool         fast=false,
             const ParallelTypeM type=AUTOMATICM);
    
  /// call init with n_local = N
  void init (const unsigned int N,
	     const bool         fast=false,
             const ParallelTypeM type=AUTOMATICM);
    
  /// Create a vector that holds tha local indices plus 
  /// those specified in the \p ghost argument.
  void init (const unsigned int /*N*/,
	     const unsigned int /*n_local*/,
	     const std::vector<unsigned int>& /*ghost*/,
	     const bool /*fast*/ = false,
             const ParallelTypeM = AUTOMATICM);

  /// Creates a vector that has the same dimension and storage type as
  /// p other, including ghost dofs.
  virtual void init (const NumericVectorM& other,
                     const bool fast = false){
  this->init(other.size(),other.local_size(),fast,other.type());
};
    
 

  // ===========================
  // SETTING FUNCTIONS
  // ===========================
    /// v(i) = value
  void set (const unsigned int i, const Real value);
  /// v(i) += value
  void add (const unsigned int i, const Real value);
  
   /// \f$U(0-N) = s\f$: fill all components.
  NumericVectorM & operator= (const Real s);
  ///  \f$U = V\f$: copy all components.
  NumericVectorM & operator= (const NumericVectorM &V);
  ///  \f$U = V\f$: copy all components.
  LaspackVectorM & operator= (const LaspackVectorM &V);
  ///  \f$U = V\f$: copy all components.
  NumericVectorM & operator= (const std::vector<Real>  &v);
  
  ///  insert the std::vect  V into the vector
  virtual void insert (const std::vector<Real> & v,
		       const std::vector<unsigned int>& dof_indices);
   ///  insert the numvect  V into the vector
  virtual void insert (const NumericVectorM& V,
		       const std::vector<unsigned int>& dof_indices); 
  ///  insert the vect  V into the vector
  virtual void insert (const DenseVectorM & V,
		       const std::vector<unsigned int>& dof_indices);
   ///  insert the subvect V into the vector
  virtual void insert (const DenseSubVectorM & V,
		       const std::vector<unsigned int>& dof_indices);
   
 // ===========================
  // return FUNCTIONS
  // ===========================  
  /// @returns the minimum element in the vector.
  Real min () const;
  /// @returns the maximum element in the vector.
  Real max () const;
  /// @returns the sum of values in a vector
  Real sum () const; 
  /// @returns the \f$l_1\f$-norm of the vector
  Real l1_norm () const;
  /// @returns the \f$l_2\f$-norm of the vector
  Real l2_norm () const;
  /// @returns the maximum absolute value of the elements of this vector
  Real linfty_norm () const;

  // ===========================  
  /// @returns dimension of the vector. 
  unsigned int size () const;
  /// @returns the local size of the vector
  unsigned int local_size() const;
  /// @returns the index of the first vector element
  unsigned int first_local_index() const;
  /// @returns the index of the last vector element
  unsigned int last_local_index() const;
     
  // ===========================
  // algebra FUNCTIONS
  // ===========================  
  /// Access components, returns \p U(i).
  Real operator() (const unsigned int i) const;  
  /// Addition operator. \p U.add(1, V).
  NumericVectorM & operator += (const NumericVectorM &V);
  /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
  NumericVectorM & operator -= (const NumericVectorM &V);
    
  // ===========================  
  ///  add the scalar s to the vector
  void add (const Real s); 
  ///  add the V to the vector
  void add (const NumericVectorM& V);
  ///  add the a*V to the vector
  void add (const Real a, const NumericVectorM& v);
  ///  add the v +=A*V to the vector
  void add_vector (const NumericVectorM &,const SparseMatrixM  &A);	
  void add_vector (const NumericVectorM& /*A*/,const SparseMMatrixM& /*V*/);
   ///  the residua  v =b-A*V to the vector
  void resid (const NumericVectorM &/*rhs_in*/,const NumericVectorM& /*x_in*/,const SparseMatrixM& /*A_in*/);
   ///  matrix multiplication  v =A*V 
  void matrix_mult (const NumericVectorM &vec_in,const SparseMMatrixM &mat_in);
  
  ///  add the std::vect  V to the vector
  void add_vector (const std::vector<Real> & v,
		   const std::vector<unsigned int>& dof_indices);
  ///  add the num vect  V to the vector
  void add_vector (const NumericVectorM& V,
		   const std::vector<unsigned int>& dof_indices);
  ///  add the vect  V to the vector
  void add_vector (const DenseVectorM & V,
		   const std::vector<unsigned int>& dof_indices);
  

  ///Scale each element
  void scale (const Real factor);

  /// v = abs(v)
  virtual void abs();

  /// Computes the dot product, p = U.V
  virtual Real dot(const NumericVectorM& V) const;

  // ===========================  
  // PARALLEL OPERATIONS
  // ===========================
  /// Creates a copy of the global vector in the local vector \p v_local.
  void localize (std::vector<Real> & v_local) const;
  /// Same, but fills a \p NumericVectorM instead of a \p std::vector.
  void localize (NumericVectorM& v_local) const;
  /// Creates a local vector \p v_local containing defined by the \p send_list.
  void localize (NumericVectorM& v_local,
		 const std::vector<unsigned int>& send_list) const;
  /// Updates a local vector with selected values from neighboring
  /// processors, as defined by \p send_list.
  void localize (const unsigned int first_local_idx,
		 const unsigned int last_local_idx,
		 const std::vector<unsigned int>& send_list);
  /// Creates a local copy of the global vector in
  /// \p v_local only on processor \p proc_id.  By
  /// default the data is sent to processor 0.  
  void localize_to_one (std::vector<Real> & v_local,
			const unsigned int proc_id=0) const;
    
  /// Computes the pointwise (i.e. component-wise) product of \p vec1
  /// and \p vec2 and stores the result in \p *this.
  virtual void pointwise_mult (const NumericVectorM& vec1,
			       const NumericVectorM& vec2);

  /// Swaps the raw QVector contents.
  virtual void swap (NumericVectorM &v);
  
  // PRINT
   /// Print the contents of the matrix in hdf5 sparse matrix format. 
 void print_hdf5(const std::string /*name="NULL"*/) const;
void print_personal(std::ostream& os) const;
 private:
  /// Make other Laspack datatypes friends
  friend class LaspackLinearSolverM;
};



// ======================================
// LaspackVectorM inline methods
// ====================================== 
inline LaspackVectorM::LaspackVectorM (const ParallelTypeM type){  this->_type = type;}
// ======================================
inline LaspackVectorM::LaspackVectorM (const unsigned int n,const ParallelTypeM type){  
  this->init(n, n, false, type);
}
// ======================================
inline LaspackVectorM::LaspackVectorM (const unsigned int n,
				 const unsigned int n_local,
                                 const ParallelTypeM type){
  this->init(n, n_local, false, type);
}
// ======================================  
inline LaspackVectorM::LaspackVectorM (const unsigned int N,
	                         const unsigned int n_local,
	                         const std::vector<unsigned int>& ghost,
                                 const ParallelTypeM type){
  this->init(N, n_local, ghost, false, type);
}
// ====================================== 
// ====================================== 
inline LaspackVectorM::~LaspackVectorM (){this->clear ();}
// ==================================   
inline void LaspackVectorM::clear (){
  if (this->initialized())  V_Destr (&_vec);
  this->_is_closed = this->_is_initialized = false;
}
// ======================================  
// ======================================  
inline void LaspackVectorM::init (const unsigned int n,
			     const unsigned int n_local,
			     const bool fast,
                             const ParallelTypeM){
  // Laspack vectors only for serial cases,
  assert (n == n_local);  this->_type = SERIALM;
  // Clear initialized vectors
  if (this->initialized()) {std::cout << "LaspackVectorM::init: already initialized";exit(0);}  // this->clear();
  // create a sequential vector
  static int cnt = 0;  char foo[80];  std::sprintf(foo,  "Vec-%d", cnt++); 
  V_Constr(&_vec, const_cast<char*>(foo), n, Normal, _LPTrue);
  this->_is_initialized = true;
  // Optionally zero out all components
  if (fast == false)    this->zero ();
  return;
}
// ==================================   
inline void LaspackVectorM::init (const unsigned int n,
			     const bool fast,
                             const ParallelTypeM type){
  this->init(n,n,fast,type);
}
// ==================================   
inline void LaspackVectorM::init (const unsigned int n,
			     const unsigned int n_local,
	                     const std::vector<unsigned int>& ghost,
			     const bool fast,
                             const ParallelTypeM type){
  assert(ghost.empty());  this->init(n,n_local,fast,type);
}
// ====================================== 


/* Default implementation for solver packages for which ghosted
   vectors are not yet implemented.  */
// void LaspackVectorM::init (const NumericVectorM& other,
//                              const bool fast)
// {
//   this->init(other.size(),other.local_size(),fast,other.type());
// }



// ==================================   
inline void LaspackVectorM::close (){assert (this->initialized()); this->_is_closed = true;}

// ================================== 
inline void LaspackVectorM::zero (){ // V_SetAllCmp (&_vec, 0.)
  assert (this->initialized());   for(uint i=1;i<= this->size();i++) V_SetCmp(&_vec,i,0.);
}

// ==================================  
inline std::unique_ptr<NumericVectorM > LaspackVectorM::clone () const{
  std::unique_ptr<NumericVectorM > cloned_vector (new LaspackVectorM);
  cloned_vector->init(*this, true);
  *cloned_vector = *this;
  return cloned_vector;
}

// ==================================   
inline unsigned int LaspackVectorM::size () const{
  assert (this->initialized());
  return static_cast<unsigned int>(V_GetDim(const_cast<QVector*>(&_vec)));
}
// ==================================   
inline unsigned int LaspackVectorM::local_size () const{
  assert (this->initialized());  return this->size();
}

// ================================== 
inline unsigned int LaspackVectorM::first_local_index () const{
  assert (this->initialized());
  return 0;
}

// ==================================    
inline unsigned int LaspackVectorM::last_local_index () const{
  assert (this->initialized());  return this->size();
}

// ==================================    
inline void LaspackVectorM::set (const unsigned int i, const Real value){
  assert (this->initialized());  assert (i < this->size());
  V__SetCmp(&_vec, i+1, value);
}

// ==================================  
inline void LaspackVectorM::add (const unsigned int i, const Real value){
  assert (this->initialized());  assert (i < this->size());
  V__AddCmp (&_vec, i+1, value);
}

// ==================================
inline Real LaspackVectorM::operator() (const unsigned int i) const{
  assert (this->initialized()); assert ( ((i >= this->first_local_index()) &&  (i <  this->last_local_index())) );
  return V__GetCmp(&_vec, i+1);
}

// ==================================
inline void LaspackVectorM::swap (NumericVectorM &other){
  LaspackVectorM& v = static_cast<LaspackVectorM&>(other);
  // This is all grossly dependent on Laspack version...
  std::swap(_vec.Name, v._vec.Name);
  std::swap(_vec.Dim, v._vec.Dim);
  std::swap(_vec.Instance, v._vec.Instance);
  std::swap(_vec.LockLevel, v._vec.LockLevel);
  std::swap(_vec.Multipl, v._vec.Multipl);
  std::swap(_vec.OwnData, v._vec.OwnData);
  // This should still be O(1), since _vec.Cmp is just a pointer to data on the heap
  std::swap(_vec.Cmp, v._vec.Cmp);
}



#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifdef __laspack_vector_h__
