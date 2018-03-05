#ifndef __numeric_vectorM_h__
#define __numeric_vectorM_h__

#include "Printinfo_conf.h"
#include "Solverlib_conf.h"
#include "SolverPackage_enum.h"
#include "Paralleltype_enum.h"
#include "parallel_objectM.h"


// C++ includes
#include <vector>
#include <set>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <iostream>
// Local includes
// #include "libmesh_common.h"
// #include "enum_parallel_type.h"
// #include "enum_solver_package.h"
#define TOLERANCEM (1.e-20)

// forward declarations
//  class NumericVectorM;
class DenseVectorM;
class DenseSubVectorM;
class SparseMatrixM;
class SparseMMatrixM;
// template <typename T> class ShellMatrix;


// ==========================================
// Numeric vector. Provides a uniform interface
// to vector storage schemes for different linear
// algebra libraries.
// ===============================================

class NumericVectorM:
#ifdef LM_REFCOUNT
    public ReferenceCountedObject<NumericVectorM>,
#endif
    public ParallelObjectM
{
    // =====================================
    // DATA
    // =====================================
protected:

    /// Flag to see if the Numeric assemble routines have been called yet
    bool _is_closed;
    /// Flag to tell if init  has been called yet
    bool _is_initialized;
    /// Type of vector
    ParallelTypeM _type;

    // =====================================
    // Constructor /Destructor
    // =====================================
public:
    /// Dummy-Constructor. Dimension=0
    explicit
    NumericVectorM (
        const ParallelM::Communicator &comm_in,     ///< parallel communicator <-  
        const ParallelTypeM = AUTOMATICM
		   );

    /// Constructor. Set dimension to \p n and initialize all elements with zero.
    explicit
    NumericVectorM (
        const ParallelM::Communicator &comm_in,         ///< parallel communicator <-
        const  int n,                                   ///< global index
        const ParallelTypeM = AUTOMATICM                ///< type
		   );

    /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
    NumericVectorM (
        const ParallelM::Communicator &comm_in,         ///< parallel communicator <-
        const  int n,                                   ///< global index
        const  int n_local,                             ///< lobal index
        const ParallelTypeM = AUTOMATICM                ///< type
		   );             

    /// Constructor. Set local dimension to \p n_local, the global dimension to \p n
    NumericVectorM (
        const ParallelM::Communicator &comm_in,            ///< parallel communicator <-
        const  int N,                                      ///< global index
        const  int n_local,                                ///< lobal index
        const std::vector< int>& ghost,                    ///< ghost vector
        const ParallelTypeM = AUTOMATICM                   ///< type
		   );

/// Builds a \p NumericVectorM using the linear solver package
/// specified by \p solver_package
    static std::unique_ptr<NumericVectorM>
    build(
        const ParallelM::Communicator &comm_in,              ///< parallel communicator <-
        const SolverPackageM solver_package = LSOLVER        ///< solver package <-
    );
    /// Creates a copy of this vector and returns it in an \p AutoPtr.
    virtual std::unique_ptr<NumericVectorM > clone () const = 0;
    virtual std::unique_ptr<NumericVectorM > zero_clone () const = 0;

    /// Destructor, deallocates memory.
    virtual ~NumericVectorM () {
        clear ();
    }
    /// @returns the \p NumericalVectorM to a pristine state.
    virtual void clear () {
        _is_closed= false;
        _is_initialized = false;
    }

/// Call the assemble functions
    virtual void close () = 0;
    /// Change the dimension of the vector to \p N. The reserved memory for
    /// this vector remains unchanged if possible, to make things faster, but
    /// this may waste some memory, so take this in the back of your head.
    /// However, if \p N==0 all memory is freed, i.e. if you want to resize
    /// the vector and release the memory not needed, you have to first call
    /// \p init(0) and then \p init(N).
    virtual void init (const  int,
                       const  int,
                       const bool = false,
                       const ParallelTypeM = AUTOMATICM) = 0;

    /// call init with n_local = N,
    virtual void init (const  int,
                       const bool = false,
                       const ParallelTypeM = AUTOMATICM) = 0;

    /// Create a vector that holds tha local indices plus those specified
    /// in the \p ghost argument.
    virtual void init (const  int /*N*/,
                       const  int /*n_local*/,
                       const std::vector< int>& /*ghost*/,
                       const bool /*fast*/ = false,
                       const ParallelTypeM = AUTOMATICM) = 0;

    /// Creates a vector that has the same dimension and storage type as
    /// \p other, including ghost dofs.
    virtual void init (const NumericVectorM& other,
                       const bool fast = false) = 0;


///  Creates the subvector "subvector" from the indices in the
/// "rows" array.  Similar to the create_submatrix routine for
/// the SparseMatrix class, it is currently only implemented for
/// PetscVectors.
    virtual void create_subvector(NumericVectorM&,const std::vector< int>&) const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        exit(0);
    }

    // =====================================
    // SETTING FUNCTIONS
    // =====================================
    /// v(i) = value
    virtual void set (const  int i, const double value) = 0;
    /// v(i) += value
    virtual void add (const  int i, const double value) = 0;

/// Set all entries to zero. Equivalent to \p v = 0
    virtual void zero () = 0;
    /// \f$U(0-N) = s\f$: fill all components.
    virtual NumericVectorM & operator= (const double s) = 0;
    ///  \f$U = V\f$: copy all components.
    virtual NumericVectorM & operator= (const NumericVectorM &V) = 0;
    ///  \f$U = V\f$: copy all components.
    virtual NumericVectorM & operator= (const std::vector<double> &v) = 0;

    /// \f$ U=v \f$ where v is a std::vector<double> andspecify WHERE to insert it
    virtual void insert (const std::vector<double>& v,
                         const std::vector< int>& dof_indices) = 0;
    /// \f$U=V\f$, and specify WHERE to insert
    virtual void insert (const NumericVectorM& V,
                         const std::vector< int>& dof_indices) = 0;
    /// \f$ U=V \f$ insert
    virtual void insert (const DenseVectorM& V,
                         const std::vector< int>& dof_indices) = 0;
    /// \f$ U=V \f$ and specify WHERE to insert it
    virtual void insert (const DenseSubVectorM& V,
                         const std::vector< int>& dof_indices) = 0;
    // =====================================
    // RETURN FUNCTIONS
    // =====================================
    /// @returns true if the vector has been initialized,
    virtual bool initialized() const {
        return _is_initialized;
    }
    /// @returns true if the vector is closed and ready
    virtual bool closed() const {
        return _is_closed;
    }

    /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
    ParallelTypeM type() const {
        return _type;
    }
    /// @returns the type (SERIAL, PARALLEL, GHOSTED) of the vector.
    ParallelTypeM & type() {
        return _type;
    }


    /// @returns the minimum element in the vector.
    virtual double min () const = 0;
    /// @returns the maximum element in the vector.
    virtual double max () const = 0;
    /// returns the sum of the elements in a vector
    virtual double sum() const = 0;

    /// @returns the \f$l_1\f$-norm of the vector, i.e.
    virtual double l1_norm () const = 0;
    /// @returns the \f$l_2\f$-norm of the vector, i.e.
    virtual double l2_norm () const = 0;
    /// @returns the maximum absolute value of the
    virtual double linfty_norm () const = 0;

    /// @returns the \f$l_1\f$-norm of the vector, i.e.
    virtual double subset_l1_norm (const std::set< int> & indices);
    /// @returns the \f$l_2\f$-norm of the vector, i.e.
    virtual double subset_l2_norm (const std::set< int> & indices);
    /// @returns the maximum absolute value of the
    virtual double subset_linfty_norm (const std::set< int> & indices);


    /// @returns dimension of the vector.
    virtual  int size () const = 0;
    /// @returns the local size of the vector (index_stop-index_start).
    virtual  int local_size() const = 0;
    /// @returns the index of the first vector element
    virtual  int first_local_index() const = 0;
    /// @returns the index+1 of the last vector element
    virtual  int last_local_index() const = 0;

    ///Access components, returns \p U(i).
    virtual double operator() (const  int i) const = 0;
    /// @returns the element \p U(i)
    virtual double el(const  int i) const {
        return (*this)(i);
    }

    /**
     * Access multiple components at once.  \p values will be resized,
     * if necessary, and filled.  The default implementation calls \p
     * operator() for each index, but some implementations may supply
     * faster methods here.
     */
    virtual void get(const std::vector< int>& index, std::vector<double>& values) const;

    // =====================================
    // algebra FUNCTIONS
    // =====================================


    /// Addition operator. Fast equivalent to \p U.add(1, V).
    virtual NumericVectorM & operator += (const NumericVectorM &V) = 0;
    /// Subtraction operator. Fast equivalent to \p U.add(-1, V).
    virtual NumericVectorM & operator -= (const NumericVectorM &V) = 0;
    /// Multiplication operator. Equivalent to \p U.scale(a)
    NumericVectorM & operator *= (const double a) {
        this->scale(a);
        return *this;
    }
    /// Division operator. Equivalent to \p U.scale(1./a)
    NumericVectorM & operator /= (const double a) {
        this->scale(1./a);
        return *this;
    }
    /// Replace each entry v_i of this vector by its reciprocal, 1/v_i.
    virtual void reciprocal() = 0;


    /// \f$U(0-LIBMESH_DIM)+=s\f$. Addition of \p s to all components. Note
    virtual void add (const double s) = 0;
    /// \f$U+=V\f$: Simple vector addition, equal to the
    virtual void add (const NumericVectorM& V) = 0;
    /// \f$U+=a*V\f$. Simple vector addition, equal to the
    virtual void add (const double a, const NumericVectorM& v) = 0;

    /// \f$ U+=v \f$ where v is a DenseVectorM
    virtual void add_vector (const std::vector<double>& v,
                             const std::vector< int>& dof_indices) = 0;
    /// \f$U+=V\f$, where U and V are type
    virtual void add_vector (const NumericVectorM& V,
                             const std::vector< int>& dof_indices) = 0;
    /// \f$U+=A*V\f$, add the product A*v
    virtual void add_vector (const NumericVectorM& /*A*/,const SparseMatrixM& /*V*/) = 0;
    virtual void add_vector (const NumericVectorM& /*A*/,const SparseMMatrixM& /*V*/) = 0;
    virtual void resid (const NumericVectorM &/*rhs_in*/,const NumericVectorM& /*A*/,const SparseMatrixM& /*V*/) = 0;
    virtual void matrix_mult (const NumericVectorM &vec_in,const SparseMMatrixM &mat_in) = 0;
    /// \f$U+=A*V\f$, add the product of a \p ShellMatrix \p
//   void add_vector (const NumericVectorM& v,
// 		   const ShellMatrix<double>& a);
//
    /// \f$ U+=V \f$ where U and V are type
    virtual void add_vector (const DenseVectorM& V,
                             const std::vector< int>& dof_indices) = 0;
    /// \f$U+=A^T*V\f$, add the product of the transpose of a \p SparseMatrix \p A_trans
    /// and a \p NumericVector \p V to \p this, where \p this=U.
    virtual void add_vector_transpose (
        const NumericVectorM&,
        const SparseMatrixM&
    ) = 0;


    /// Scale each element
    virtual void scale (const double factor) = 0;
    /// v = abs(v)... that is, each entry in v is replaced
    virtual void abs() = 0;
    /// Computes the dot product, p = U.V
    virtual double dot(const NumericVectorM&) const = 0;

    /// Exchanges the values/sizes of two vectors.
    virtual void swap (NumericVectorM &v);

    // =====================================
    // PARALLEL OPERATIONS
    // =====================================

    /// Creates a copy of the global vector in the local vector \p v_local.
    virtual void localize (std::vector<double>& v_local) const = 0;
    /// Same, but fills a \p NumericVectorM instead
    virtual void localize (NumericVectorM& v_local) const = 0;
    /// Creates a local vector \p v_local
    virtual void localize (NumericVectorM& v_local,
                           const std::vector< int>& send_list) const = 0;
    /// Updates a local vector with selected values from neighboring
    virtual void localize (const  int first_local_idx,
                           const  int last_local_idx,
                           const std::vector< int>& send_list) = 0;
    /// Creates a local copy of the global vector
    virtual void localize_to_one (std::vector<double>& v_local,
                                  const  int proc_id=0) const = 0;
    /// @returns \p -1 when \p this is equivalent to \p other_vector,
    virtual int compare (const NumericVectorM &other_vector,
                         const double threshold = TOLERANCEM) const;
    /// @returns \p -1 when \p this is equivalent to \p other_vector,
    /// up to the given local relative \p threshold.  When differences
    /// occur, the return value contains the first index where
    /// the difference \p (a[i]-b[i])/max(a[i],b[i]) exceeded the
    /// threshold.  When no threshold is given, the \p libMesh\p TOLERANCEM is used.
    virtual int local_relative_compare (
        const NumericVectorM &other_vector,
        const double threshold = TOLERANCEM) const;

    /*
     * @returns \p -1 when \p this is equivalent to \p other_vector,
     * up to the given local relative \p threshold.  When differences
     * occur, the return value contains the first index where
     * the difference \p (a[i]-b[i])/max_j(a[j],b[j]) exceeded the
     * threshold.  When no threshold is given, the \p libMesh
     * \p TOLERANCEM is used.
     */
    virtual int global_relative_compare (
        const NumericVectorM &other_vector,
        const double threshold = TOLERANCEM) const;

    /// Computes the pointwise (i.e. component-wise) product of \p vec1
    virtual void pointwise_mult (const NumericVectorM& vec1,
                                 const NumericVectorM& vec2) = 0;


    // =====================================
    // PRINTING FUNCTIONS
    // =====================================
    /// Prints the local contents of the vector to the screen.
    virtual void print(std::ostream& os=std::cout) const;
    /// Prints the local contents of the vector to the screen.
    virtual void print_personal(std::ostream& os=std::cout) const=0;
    /// Prints the global contents of the vector to the screen.
    virtual void print_global(std::ostream& os=std::cout) const;
    /// Same as above but allows you to use stream syntax.
    friend std::ostream& operator << (std::ostream& os, const NumericVectorM& v) {
        v.print_global(os);
        return os;
    }

    /// Print the contents of the matrix in hdf5 sparse matrix format.
    virtual void print_hdf5(const std::string /*name="NULL"*/) const=0;


};


/*----------------------- Inline functions ----------------------------------*/



// ==============================================
inline NumericVectorM::NumericVectorM (
    const ParallelM::Communicator &comm_in,
    const ParallelTypeM type) :
    ParallelObjectM(comm_in),
    _is_closed(false),  _is_initialized(false),  _type(type) {}
// ==============================================
inline NumericVectorM::NumericVectorM (
    const ParallelM::Communicator &comm_in,
    const  int /*n*/,
    const ParallelTypeM type) :
    ParallelObjectM(comm_in),_is_closed(false),_is_initialized(false), _type(type) {
    std::cout<< "Abstract base class! ";
    exit(0); // Abstract base class!
    // init(n, n, false, type);
}

// ==============================================
inline NumericVectorM::NumericVectorM (
    const ParallelM::Communicator &comm_in,
    const int /*n*/,const int /*n_local*/,
    const ParallelTypeM type) :
    ParallelObjectM(comm_in),_is_closed(false),  _is_initialized(false),  _type(type) {
    std::cout<< "Abstract base class! ";
    exit(0); // Abstract base class!
    // init(n, n_local, false, type);
}

// ==============================================
inline NumericVectorM::NumericVectorM (
    const ParallelM::Communicator &comm_in,
    const int /*n*/,const int /*n_local*/,
    const std::vector<int>& /*ghost*/,
    const ParallelTypeM type) :
    ParallelObjectM(comm_in),_is_closed(false),  _is_initialized(false),  _type(type) {
    std::cout<< "Abstract base class! ";
    exit(0); // Abstract base class!
    // init(n, n_local, ghost, false, type);
}


// ==============================================
inline void NumericVectorM::get(const std::vector<int>& index,
                                std::vector<double>& values) const {
    const std::size_t num = index.size();
    values.resize(num);
    for( std::size_t i=0; i<num; i++) values[i] = (*this)(index[i]);
}

// ==============================================
inline void  NumericVectorM::swap (NumericVectorM &v) {
    std::swap(_is_closed, v._is_closed);
    std::swap(_is_initialized, v._is_initialized);
    std::swap(_type, v._type);
}

// ==============================================
inline void NumericVectorM::print(std::ostream& os) const {
    assert (this->initialized());
    os << "Size\tglobal =  " << this->size()
       << "\t\tlocal =  " << this->local_size() << std::endl;
    os << "#\tValue" << std::endl;
    for ( int i=this->first_local_index(); i<this->last_local_index(); i++)
        os << i << "\t" << (*this)(i) << std::endl;
}
// ==============================================
inline void NumericVectorM::print_global(std::ostream& os) const {
    assert (this->initialized());
    std::vector<double> v(this->size());
    this->localize(v);
    // Right now we only want one copy of the output
    if (this->processor_id()) return;
    os << "Size\tglobal =  " << this->size() << std::endl;
    os << "#\tValue" << std::endl;
    for (unsigned int i=0; i!=v.size(); i++) os << i << "\t" << v[i] << std::endl;
//   printf("Size\tglobal =  %d \n", this->size());
//   printf("#\tValue \n");
//   for (int i=0; i!=v.size(); i++) printf(" %d \t %g \n", i, v[i]);
}




#endif  // #ifdef __numeric_vector_h__
