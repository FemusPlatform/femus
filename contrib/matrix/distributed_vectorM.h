#ifndef __distributed_vectorMM_h__
#define __distributed_vectorMM_h__

#include "Solverlib_conf.h"
#include "Paralleltype_enum.h"
// ++++++++++++++++++++++++++++++++
#define numeric_index_type int
#define processor_id_type int
// ++++++++++++++++++++++++++++++++

// C++ includes
#include <vector>
#include <algorithm>
// #include <complex>
#include <limits>

// Local includes
#include "numeric_vectorM.h"
#include "parallelM.h"

// ======================================
// Distributed vector. Provides an interface for simple
// parallel, distributed vectors. Offers some collective
// communication capabilities.  Note that the class will
// sill function without MPI, but only on one processor.
// This lets us keep the parallel details behind the scenes.
// =====================================================
class DistributedVectorM : public NumericVectorM
{ // ==============================================
private:

    /// Actual vector datatype to hold vector entries
    std::vector<double> _values;
    /// The global vector size
    numeric_index_type _global_size;
    /// The local vector size
    numeric_index_type _local_size;

    /// The first component stored locally
    numeric_index_type _first_local_index;
    /// The last component (+1) stored locally
    numeric_index_type _last_local_index;


public:

    /// ** Dummy-Constructor. Dimension=0 
    explicit
    DistributedVectorM (
      const ParallelM::Communicator &comm,
      const ParallelTypeM = AUTOMATICM);

    /// ** Constructor. Set dimension to \p n and initialize all elements with zero.
    explicit
    DistributedVectorM (
      const ParallelM::Communicator &comm,
      const numeric_index_type n,
      const ParallelTypeM type = AUTOMATICM);

    /// Constructor. Set local dimension to \p n_local,
    DistributedVectorM (
      const ParallelM::Communicator &comm,
      const numeric_index_type n,
      const numeric_index_type n_local,
      const ParallelTypeM type = AUTOMATICM);

    /// Constructor. Set local dimension to \p n_local
    DistributedVectorM (
                     const ParallelM::Communicator &comm,
		     const numeric_index_type N,
		     const numeric_index_type n_local,
		     const std::vector<numeric_index_type>& ghost,
                     const ParallelTypeM ptype = AUTOMATICM); 
 
    /// Destructor, deallocates memory. Made virtual to allow
    ~DistributedVectorM ();

    /// Call the assemble functions
    void close ();

    /// @returns the \p DistributedVectorM to a pristine state.
    void clear ();

    /// Set all entries to zero. Equivalent to \p v = 0, but more obvious and faster.
    void zero ();

    /// Creates a copy of this vector and returns it in an \p AutoPtr.
    std::unique_ptr<NumericVectorM > clone () const;
    std::unique_ptr<NumericVectorM > zero_clone () const;
    
    /// init 1. Change the dimension of the vector to \p N. The reserved memory for
    // this vector remains unchanged if possible, to make things faster, but
    // this may waste some memory, so take this in the back of your head.
    // However, if \p N==0 all memory is freed, i.e. if you want to resize
    // the vector and release the memory not needed, you have to first call
    // \p init(0) and then \p init(N).  On \p fast==false, the vector is filled by zeros.
    void init (
      const numeric_index_type N,
      const numeric_index_type n_local,
      const bool         fast=false,
      const ParallelTypeM ptype=AUTOMATICM
    );


    /// init 2. call init with n_local = N,
    void init (
             const numeric_index_type N,
	     const bool         fast=false,
	     const ParallelTypeM ptype=AUTOMATICM
	      );      
 

    /// init 3. Create a vector that holds tha local indices plus those specified in the \p ghost argument.
    void init (
                     const numeric_index_type /*N*/,
		     const numeric_index_type /*n_local*/,
		     const std::vector<numeric_index_type>& /*ghost*/,
		     const bool /*fast*/ = false,
		     const ParallelTypeM = AUTOMATICM);     

    /// init4. Creates a vector that has the same dimension and storage type as \p other, including ghost dofs.
    void init (const NumericVectorM& other,
                       const bool fast = false);
//     this->init(other.size(),other.local_size(),fast,other.type());
//     }
    // ===================================
    //  SETTINGS
    // ===================================

    /// v(i) = value
    void set (
      const numeric_index_type i, 
      const double value
    );
    /// v(i) += value
    void add (
      const  numeric_index_type i, 
      const double value
    );

    /// \f$U(0-N) = s\f$: fill all components.
    NumericVectorM & operator= (const double s);
    ///  \f$U = V\f$: copy all components.
    NumericVectorM & operator= (const NumericVectorM &V);

    ///  \f$U = V\f$: copy all components.
    DistributedVectorM & operator= (const DistributedVectorM &V);
    ///  \f$U = V\f$: copy all components.
    NumericVectorM & operator= (const std::vector<double> &v);

    /// \f$ U=v \f$ where v is a std::vect
    virtual void insert (
      const std::vector<double>& v,
      const std::vector< numeric_index_type>& dof_indices);
    /// \f$ U=v \f$ where v is a NumericVect
    virtual void insert (
      const NumericVectorM& V,
      const std::vector<numeric_index_type>& dof_indices);
    /// \f$ U=V \f$ where V is type
    // DenseVectorM and you want to specify WHERE to insert it
    virtual void insert (
      const DenseVectorM& V,
      const std::vector< numeric_index_type>& dof_indices);
    /// \f$ U=V \f$ where V is type  DenseSubVectorM
    virtual void insert (
      const DenseSubVectorM& V,
      const std::vector<numeric_index_type>& dof_indices);

// ===================================
//  RETURN FUNCTIONS
// ===================================

    /// @returns the minimum element in the vector.
    double min () const;
    /// @returns the maximum element in the vector.
    double max () const;
    /// @returns the sum of all values in the vector
    double sum() const;
    /// @returns the \f$l_1\f$-norm of the vector
    double l1_norm () const;
    /// @returns the \f$l_2\f$-norm of the vector
    double l2_norm () const;
    /// @returns the maximum absolute value of the elements
    double linfty_norm () const;

    /// @returns dimension of the vector.
    numeric_index_type size () const;
    /// @returns the local size of the vector (index_stop-index_start)
    numeric_index_type local_size() const;
    /// @returns the index of the first vector element actually stored on this processor
    numeric_index_type first_local_index() const;
    /// @returns the index of the last vector element
    numeric_index_type last_local_index() const;

    /// Access components, returns \p U(i).
    double operator() (const  numeric_index_type i) const;

    // ===================================
    //  ALGEBRA
    // ===================================

    /// Addition operator. Fast equivalent to \p U.add(1, V).
    NumericVectorM & operator += (const NumericVectorM &V);
    /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
    NumericVectorM & operator -= (const NumericVectorM &V);
    
    /// ** Replace each entry v_i of this vector by its reciprocal, 1/v_i.
    virtual void reciprocal();

    /// Addition of \p s (scalar) to all components. Note
    void add (const double s);
    /// \f$U+=V\f$ Simple vector addition operator +=.
    void add (const NumericVectorM& V);
    /// \f$U+=a*V\f$. Simple vector addition
    void add (const double a, const NumericVectorM& v);

    /// \f$U+=v\f$ where v is a \p std::vector<double>
    void add_vector (
      const std::vector<double>& v,
      const std::vector<numeric_index_type>& dof_indices
    );
    /// \f$U+=V\f$ where U and V are type  NumericVectorM and you
    // want to specify WHERE to add the \p NumericVectorM V
    void add_vector (
      const NumericVectorM& V,
      const std::vector< numeric_index_type>& dof_indices
    );
    /// \f$U+=A*V\f$. Add the product A*V
    void add_vector (
      const NumericVectorM&/* V*/,
      const SparseMatrixM&/* A*/
    ) {std::cout <<"error";  exit(0);}
    /// \f$U+=A*V\f$. Add the product A*V
    void add_vector (
      const NumericVectorM& /*A*/,
      const SparseMMatrixM& /*V*/
    ) {std::cout <<"error";  exit(0);}
    /// \f$U+=V\f$ with DenseVectorM V
    void add_vector (
      const DenseVectorM& V,
      const std::vector< numeric_index_type>& dof_indices
    );   
   /// *** \f$U+=A^T*V\f$. Add the product of the transpose 
   /// of a Sparse matrix \p A_trasp  and a Numeric vector
   ///  \p V to this Numeric vector.
   void add_vector_transpose (
     const NumericVectorM&,
     const SparseMatrixM&
     )  { std::cout <<"error";  exit(0);}
    
    
    /// *** residual computation new function   
     void resid (
       const NumericVectorM &/*rhs_in*/,
       const NumericVectorM& /*x_in*/,
       const SparseMatrixM& /*A_in*/
     ) {        std::cout <<"error";  exit(0);    }
    /// *** matrix mult computation new function   
    void matrix_mult (
      const NumericVectorM &/*vec_in*/,
      const SparseMMatrixM &/*mat_in*/
    ){        std::cout <<"error";  exit(0);    }	
		  
 



    /// Scale each element of the vector by the given factor.
    void scale (const double factor);
    /// v = abs(v)... that is, each entry in v is replaced by its absolute value.
    virtual void abs();

    /// Computes the dot product, p = U.V
    virtual double dot(const NumericVectorM& V) const;

    /// Computes the pointwise (i.e. component-wise) product of \p vec1
    /// and \p vec2 and stores the result in \p *this.
    virtual void pointwise_mult (const NumericVectorM& vec1,
                                 const NumericVectorM& vec2);
    /// Swaps the vector data and metadata
    virtual void swap (NumericVectorM &v);

    // ===================================
    //  DISTRIBUTED FUNCTIONS
    // ===================================
    /// Creates a copy of the global vector in the local vector \p v_local.
    void localize (std::vector<double>& v_local) const;
    /// Same, but fills a \p NumericVectorM instead of a \p std::vector.
    void localize (NumericVectorM& v_local) const;
    /// Creates a local vector \p v_local containing only
    /// information relevant to this processor, as defined by the \p send_list.
    void localize (NumericVectorM& v_local,
                   const std::vector< numeric_index_type>& send_list) const;
    ///  Updates a local vector with selected values from neighboring
    /// processors, as defined by \p send_list.
    void localize (const  int first_local_idx,
                   const  int last_local_idx,
                   const std::vector<numeric_index_type>& send_list);

    /// Creates a local copy of the global vector in \p v_local only
    /// on processor \p proc_id.  By default the data is sent to processor 0.
    /// This method is useful for outputting data from one processor.
    void localize_to_one (std::vector<double>& v_local,
                          const  processor_id_type proc_id=0) const;

  // PRINT
 /// Print the contents of the matrix in hdf5 sparse matrix format. 
 void print_hdf5(const std::string /*name="NULL"*/) const{std::cout<<"Not implemented \n";abort();}
 void print_personal(std::ostream& os=std::cout) const{os<<"Not implemented \n";abort();}
};


//--------------------------------------------------------------------------
// DistributedVectorM inline methods
// ============================================
inline DistributedVectorM::DistributedVectorM (
  const ParallelM::Communicator &comm,
  const ParallelTypeM ptype
) :     NumericVectorM(comm, ptype),
        _global_size      (0),  _local_size       (0),
        _first_local_index(0),  _last_local_index (0) {
    this->_type = ptype;
}
// ============================================
inline DistributedVectorM::DistributedVectorM (
  const ParallelM::Communicator &comm,
  const  numeric_index_type n,
  const ParallelTypeM ptype)
  : NumericVectorM(comm, ptype)
{
    this->init(n, n, false, ptype);
}
// ============================================
inline DistributedVectorM::DistributedVectorM (
  const ParallelM::Communicator &comm,
  const numeric_index_type n,
  const numeric_index_type n_local,
  const ParallelTypeM ptype) 
 : NumericVectorM(comm, ptype){
  this->init(n, n_local, false, ptype);
}
// ============================================
inline DistributedVectorM::DistributedVectorM (
  const ParallelM::Communicator &comm,
  const   numeric_index_type n,
  const  numeric_index_type n_local,
  const std::vector< numeric_index_type>& ghost,
  const ParallelTypeM ptype)
 : NumericVectorM(comm, ptype){
    this->init(n, n_local, ghost, false, ptype);
}
// ============================================
inline DistributedVectorM::~DistributedVectorM () {  this->clear (); }

// *** ============================================
inline void DistributedVectorM::init (
  const  numeric_index_type n,
  const  numeric_index_type n_local,
  const bool fast,
  const ParallelTypeM type) {

    // This function must be run on all processors at once
    parallel_onlyM();
    assert (n_local <= n);

    if (type == AUTOMATICM)    {
        if (n == n_local)    this->_type = SERIALM;
        else  this->_type = PARALLELM;
    }
    else    this->_type = type;
    assert ((this->_type==SERIALM && n==n_local) || this->_type==PARALLELM);

    // Clear the data structures if already initialized
    if (this->initialized())    this->clear();
    // Initialize data structures
    _values.resize(n_local);    _local_size  = n_local;
    _global_size = n;    _first_local_index = 0;

#ifdef HAVE_MPI
    
    
   std::vector<numeric_index_type> local_sizes (this->n_processors(), 0);
   local_sizes[this->processor_id()] = n_local;
   this->comm().sum(local_sizes);

   // _first_local_index is the sum of _local_size
   // for all processor ids less than ours
   for (processor_id_type p=0; p!=this->processor_id(); p++)
                  _first_local_index += local_sizes[p];

// old_code --------------------
//     int n_proc=0, proc_id=0;
//     MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
//     MPI_Comm_size (MPI_COMM_WORLD, &n_proc);
// 
//     std::vector<int> local_sizes     (n_proc, 0);
//     local_sizes[proc_id] = n_local;
//     ParallelM::sum(local_sizes);
// 
//     // _first_local_index is the sum of _local_size
//     // for all processor ids less than ours
//     for (int p=0; p<proc_id; p++)    _first_local_index += local_sizes[p];

#  ifdef DEBUG
    // Make sure all the local sizes sum up to the global
    // size, otherwise there is big trouble!
    numeric_index_type dbg_sum=0;
    for (processor_id_type p=0;  p!=this->n_processors(); p++)   dbg_ sum += local_sizes[p];
    assert (dbg_sum == n);

#  endif

#else
    // No other options without MPI!
    if (n != n_local)    {
        std::cerr << "ERROR:  MPI is required for n != n_local!"
                  << std::endl;
        abort();
    }

#endif

    _last_local_index = _first_local_index + n_local;
    // Set the initialized flag
    this->_is_initialized = true;
    // Zero the components unless directed otherwise
    if (!fast)    this->zero();
}
// ============================================
inline
void DistributedVectorM::init (
  const  numeric_index_type n,
  const numeric_index_type n_local,
  const std::vector< int>& /*ghost*/,
  const bool fast,
  const ParallelTypeM type)
{
    // TODO: we shouldn't ignore the ghost sparsity pattern
    this->init(n, n_local, fast, type);
}
// ============================================
/// *** Default implementation for solver packages for which ghosted
///   vectors are not yet implemented. 
inline
void DistributedVectorM::init (
  const NumericVectorM& other,
  const bool fast)
{
    this->init(other.size(),other.local_size(),fast,other.type());
}
// ============================================
inline void DistributedVectorM::init (
  const numeric_index_type n,
  const bool fast,
  const ParallelTypeM type)
{
    this->init(n,n,fast,type);
}
// ============================================
inline void DistributedVectorM::close () {
    assert (this->initialized());
    this->_is_closed = true;
}
// ============================================
inline void DistributedVectorM::clear () {
    _values.clear();
    _global_size =_local_size =_first_local_index=_last_local_index = 0;
    this->_is_closed = this->_is_initialized = false;
}
// ============================================
inline void DistributedVectorM::zero () {
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    std::fill (_values.begin(),_values.end(),0.);
}
// *** ============================================
inline std::unique_ptr<NumericVectorM > DistributedVectorM::zero_clone () const {
    std::unique_ptr<NumericVectorM > cloned_vector (new DistributedVectorM(this->comm()));
    cloned_vector->init(*this);
    return cloned_vector;
}
// ============================================
inline std::unique_ptr<NumericVectorM > DistributedVectorM::clone () const {
    std::unique_ptr<NumericVectorM > cloned_vector (new DistributedVectorM(this->comm()));
    cloned_vector->init(*this, true);
    *cloned_vector = *this;
    return cloned_vector;
}
// ============================================
inline numeric_index_type DistributedVectorM::size () const {
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    return _global_size;
}
// ============================================
inline  int DistributedVectorM::local_size () const {
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    return _local_size;
}
// ============================================
inline  int DistributedVectorM::first_local_index () const
{
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);

    return _first_local_index;
}

// ============================================
inline  int DistributedVectorM::last_local_index () const
{
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);

    return _last_local_index;
}

// ============================================
inline double DistributedVectorM::operator() (
  const  numeric_index_type i
) const{
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    assert ( ((i >= first_local_index()) &&
                      (i <  last_local_index())) );

    return _values[i - _first_local_index];
}

// ============================================
inline void DistributedVectorM::set (
  const  numeric_index_type i, 
  const double value
){
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    assert (i<size());
    assert (i-first_local_index() < local_size());

    _values[i - _first_local_index] = value;
}

// ============================================
inline void DistributedVectorM::add (
  const numeric_index_type i,
  const double value
){
    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);
    assert (i<size());
    assert (i-first_local_index() < local_size());

    _values[i - _first_local_index] += value;
}

// ** ============================================
inline double DistributedVectorM::min () const{
    // This function must be run on all processors at once
    parallel_onlyM();

    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);

    double local_min = _values.size() ?
                     (double)(_values[0]) : std::numeric_limits<double>::max();
    for (numeric_index_type i = 1; i < (int)_values.size(); ++i)
        local_min = std::min((double)(_values[i]), local_min);

    this->comm().min(local_min);
//     ParallelM::min(local_min);

    return local_min;
}

// * ============================================
inline double DistributedVectorM::max() const
{
    // This function must be run on all processors at once
    parallel_onlyM();

    assert (this->initialized());
    assert ((int)_values.size() == _local_size);
    assert ((_last_local_index - _first_local_index) == _local_size);

    double local_max = _values.size() ?
                     (double)(_values[0]) : -std::numeric_limits<double>::max();
    for ( int i = 1; i < (int)_values.size(); ++i)
        local_max = std::max((double)(_values[i]), local_max);
    this->comm().max(local_max);
//     ParallelM::max(local_max);

    return local_max;
}

// ============================================
inline void DistributedVectorM::swap (NumericVectorM &other)
{
    DistributedVectorM& v = dynamic_cast<DistributedVectorM&>(other);
    std::swap(_global_size, v._global_size);
    std::swap(_local_size, v._local_size);
    std::swap(_first_local_index, v._first_local_index);
    std::swap(_last_local_index, v._last_local_index);
    // This should be O(1) with any reasonable STL implementation
    std::swap(_values, v._values);
}


#endif  // #ifdef __distributed_vector_h__
