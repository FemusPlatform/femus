#ifndef __petsc_matrixMn_h__
#define __petsc_matrixMn_h__

#include "Solverlib_conf.h"  // solver conf library
// ======================================
#ifdef HAVE_PETSCM
// ======================================

// This class
#include "sparse_matrixM.h" // parent
#include "petsc_macroM.h"
// Petsc include files.
EXTERN_C_FOR_PETSC_BEGIN
# include <petscmat.h>
EXTERN_C_FOR_PETSC_END

// C++ includes -------------------
#include <algorithm> // equal_range

// Local includes
#include "parallelM.h"        // parallel commands



// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMatrixM methods
#undef semiparallel_onlyM
#ifndef NDEBUG
#include <cstring>

#define semiparallel_onlyM() do { if (this->initialized()) { const char *mytype; \
    MatGetType(_mat,&mytype); \
    if (!strcmp(mytype, MATSEQAIJ)) \
      parallel_onlyM(); } } while (0)
#else
#define semiparallel_onlyM()
#endif




// =======================================================
// Petsc matrix
// =======================================================

class PetscMatrixM : public SparseMatrixM {

private:
  // data ------------------------------------
  Mat _mat;                 ///< Petsc matrix pointer
  bool _destroy_mat_on_exit;///< Boolean value (false)

public:
  // Constructor ---------------------------------------------------------
  /// Constructor I;  initialize the matrix before usage with \p init(...).
  PetscMatrixM(
    const ParallelM::Communicator &comm
  ): SparseMatrixM(comm), _destroy_mat_on_exit(true) {}
  /// Constructor II.
  PetscMatrixM(
    Mat m, 
    const ParallelM::Communicator &comm
  ): SparseMatrixM(comm), _destroy_mat_on_exit(false) {
  this->_mat = m;  this->_is_initialized = true;/* this->_communicator=comm;*/
}

  /// Initialize a Petsc matrix
  void init(const int m,const int n,const int m_l,const int n_l,
            const int nnz=0, const int noz=0);
  void init (const int m,  const int n){_m=m;_n=n;  }
  void init (){}; 
  // Destructors ----------------------------
  /// Destructor
  ~PetscMatrixM();

  void clear();    ///< Release all memory
  void zero();     ///< set to zero
  void zero_rows(  ///< set  rows to zero
    std::vector<int> & rows, 
    double diag_value = 0.0
  );
  void close() const;///< close
 

  // Returns -------------------------------------------
  /// PETSc matrix context pointer
  Mat mat() { assert(_mat != NULL); return _mat; }

  // matrix dimensions
  int m() const;          ///< row-dimension
  int n() const;          ///< column dimension
  int row_start() const;  ///< row-start
  int row_stop() const;   ///< row-stop

  ///    Return the value of the entry  (i,j). Do not use
  double operator()(const int i,const int j) const; 
  /// Return the row values
  int MatGetRowM(const int i_val,int cols[]=PETSC_NULL , double vals[]=PETSC_NULL ); 
  /// Petsc matrix has been closed and fully assembled
  bool closed() const;
  
  // Setting -------------------------------------
  // update pattern
  void update_sparsity_pattern(const Graph &sparsity_pattern); ///<   sparsity patter update (Graph)
  void update_sparsity_pattern(int m,int n,int m_l,int n_l,    ///<   sparsity patter update (petsc)
    const std::vector<int>  n_oz,const std::vector<int>  n_nz);
  // set values
  void set(const int i,const int j, const double value); ///< Set the value.
  void add(const int i,const int j,const double value);///<  add  value
  // add
  /// Add the full matrix to the Petsc matrix.
  void add_matrix(const DenseMatrixM &dm,
                  const std::vector<int> &rows,
                  const std::vector<int> &cols);
  /// Add the full matrix to the Petsc matrix.
  void add_matrix(const DenseMatrixM &dm,
                  const std::vector<int> &dof_indices);

  /// Add a Sparse matrix
  void add(const double a, SparseMatrixM &X);


  // functions ------------------------------------
  double l1_norm() const;/// Return the l1-norm of the matrix
  double linfty_norm() const; /// Return the linfty-norm of the matrix

 


  // print -----------------------------------------------------
  void print_personal(       ///< print personal
    std::ostream& os=std::cout
  ) const;   
  void print_personal(       ///< print
    const std::string name="NULL"
  ) const; 
  void print_hdf5(            ///< print hdf5
    const std::string name="NULL"
  ) const;     

  // functions ---------------------------------------------------
  
  virtual void get_diagonal(///< Copies the diagonal part of the matrix
    NumericVectorM& dest
  ) const;
  virtual void get_transpose(///< Transpose Matrix.
    SparseMatrixM& dest
  ) const;
  void swap(PetscMatrixM &); ///< Swaps the raw PETSc matrix context pointers.


// ========================================================
 

protected:
  ///  get "submatrix" from sparse matrix.
  virtual void _get_submatrix(
    SparseMatrixM& submatrix,
    const std::vector<int>& rows,
    const std::vector<int>& cols,
    const bool reuse_submatrix) const;
};

// ===============================================
// PetscMatrixM inline members
// ===============================================

// ===============================================
// inline PetscMatrixM::PetscMatrixM(const ParallelM::Communicator &comm)  : _destroy_mat_on_exit(true) {}

// =================================================================
// inline PetscMatrixM::PetscMatrixM(Mat m, const ParallelM::Communicator &comm): _destroy_mat_on_exit(false) {
//   this->_mat = m;  this->_is_initialized = true;/* this->_communicator=comm;*/
// }

// ===============================================
inline PetscMatrixM::~PetscMatrixM() {  this->clear();}

// ==========================================
/// This function checks if the matrix is closed
inline void PetscMatrixM::close() const {
  parallel_onlyM();
  int ierr=0;
  ierr = MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(this->comm().get(),ierr);
  ierr = MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(this->comm().get(),ierr);
}

// ==================================================
/// This function returns the row number
inline int PetscMatrixM::m() const {
  assert(this->initialized());
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize(_mat, &petsc_m, &petsc_n); CHKERRABORT(this->comm().get(),ierr);
  return static_cast<int>(petsc_m);
}

// ============================================
/// This function returns the column number
inline int PetscMatrixM::n() const {
  assert(this->initialized());
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize(_mat, &petsc_m, &petsc_n); CHKERRABORT(this->comm().get(),ierr);
  return static_cast<int>(petsc_n);
}

// ==========================================
/// This function returns the start row location
inline int PetscMatrixM::row_start() const {
  assert(this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop); CHKERRABORT(this->comm().get(),ierr);
  return static_cast<int>(start);
}

// =========================================
/// This function returns the stop row location
inline int PetscMatrixM::row_stop() const {
  assert(this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop); CHKERRABORT(this->comm().get(),ierr);

  return static_cast<int>(stop);
}

// ======================================================
/// This function sets the value in the (i,j) pos
inline void PetscMatrixM::set(const int i,
                              const int j,
                              const double value) {
  assert(this->initialized());
  int ierr=0, i_val=i, j_val=j;
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
                      &petsc_value, INSERT_VALUES);
  CHKERRABORT(this->comm().get(),ierr);
}

// =================================================
/// This function add the  value to the element (i,j).
inline void PetscMatrixM::add(const int i,    // index i
                              const int j, // index j
                              const double value      // value
                             ) {
  assert(this->initialized());  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
                      &petsc_value, ADD_VALUES);
  CHKERRABORT(this->comm().get(),ierr);
}

// =====================================================
/// This function adds a square dense matrix to the main sparse matrix
inline void PetscMatrixM::add_matrix(
  const DenseMatrixM& dm,                       // dense matrix
  const std::vector<int>& dof_indices  // index vector
) { this->add_matrix(dm, dof_indices, dof_indices);}

// =========================================================
/// This function performs \f$\texttt{this} = a*X + \texttt{this} \f$.
inline void PetscMatrixM::add(const double a_in,      // double constant
                              SparseMatrixM &X_in  // sparse matrix
                             ) {
  assert(this->initialized());

  // crash due to incompatible sparsity structure...
  assert(this->m() == X_in.m());  assert(this->n() == X_in.n());

  PetscScalar     a = static_cast<PetscScalar>(a_in);
  PetscMatrixM*   X = dynamic_cast<PetscMatrixM*>(&X_in);

  // the matrix from which we copy the values has to be assembled/closed
  // X->close ();
  assert(X != NULL); assert(X->closed());
  semiparallel_onlyM(); int ierr=0;

#if PETSC_VERSION_LESS_THAN(2,3,0)  // 2.2.x & earlier style  
  ierr = MatAXPY(&a,  X->_mat, _mat, SAME_NONZERO_PATTERN);
  CHKERRABORT(this->comm().get(),ierr);
#else    // 2.3.x & newer 
  ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(this->comm().get(),ierr);
#endif
}

// ========================================================
inline double PetscMatrixM::operator()(const int i, // index i
                                     const int j  // index j
                                    ) const {

  // do no use this function !!!!!!!!
  assert(this->initialized());
#if PETSC_VERSION_LESS_THAN(2,2,1)  // PETSc 2.2.0 & older
  PetscScalar *petsc_row;
  int* petsc_cols;
#else   // PETSc 2.2.1 & newer
  const PetscScalar *petsc_row;
  const PetscInt    *petsc_cols;
#endif

  // If the entry is not in the sparse matrix, it is 0.
  double value=0.;
  int  ierr=0, ncols=0,  i_val=static_cast<int>(i),    j_val=static_cast<int>(j);

  assert(this->closed());// the matrix needs to be closed for this to work
  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
  CHKERRABORT(this->comm().get(),ierr);
  // Perform a binary search to find the contiguous index in
  // petsc_cols (resp. petsc_row) corresponding to global index j_val
  std::pair<const int*, const int*> p =
    std::equal_range(&petsc_cols[0], &petsc_cols[0] + ncols, j_val);
  // Found an entry for j_val
  if (p.first != p.second)    {
    // The entry in the contiguous row corresponding to the j_val column of interest
    const int j = std::distance(const_cast<int*>(&petsc_cols[0]),const_cast<int*>(p.first));

    assert(j < ncols);    assert(petsc_cols[j] == j_val);
    value = static_cast<double>(petsc_row[j]);
  }
  ierr  = MatRestoreRow(_mat, i_val,
                        &ncols, &petsc_cols, &petsc_row);
  CHKERRABORT(this->comm().get(),ierr);
  return value;
}


// ========================================================
inline int PetscMatrixM::MatGetRowM(const int i_val, // index i
                                    int cols[],
                                    double vals[]
                                   ) {
  // Set up
  assert(this->initialized());
  const PetscScalar *petsc_row;
  const PetscInt    *petsc_cols;
  int  ierr=0, ncols=0;
  // Get row
  assert(this->closed());// the matrix needs to be closed for this to work
  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
  CHKERRABORT(this->comm().get(),ierr);
  // close row
  ierr  = MatRestoreRow(_mat, i_val,&ncols, &petsc_cols, &petsc_row);
  CHKERRABORT(this->comm().get(),ierr);
  // print row
  if(&cols[0] !=PETSC_NULL)   for (int j=0;j<ncols;j++)  {cols[j]=petsc_cols[j]; vals[j]=petsc_row[j];}
  // return # of columns
  return ncols;
}


// ================================================
inline bool PetscMatrixM::closed() const {
  assert(this->initialized());
  int ierr=0;
  PetscBool assembled;
  ierr = MatAssembled(_mat, &assembled);
  CHKERRABORT(this->comm().get(),ierr);
  return (assembled == PETSC_TRUE);
}

// =========================================
inline void PetscMatrixM::swap(
  PetscMatrixM &m // pointer to swap
) {// =========================================
  std::swap(_mat, m._mat);
  std::swap(_destroy_mat_on_exit, m._destroy_mat_on_exit);
}






#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
