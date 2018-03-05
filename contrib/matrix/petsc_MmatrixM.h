#ifndef __petsc_MmatrixMn_h__
#define __petsc_MmatrixMn_h__

#include "Solverlib_conf.h"
// ===================================
#ifdef HAVE_PETSCM
// ===================================

// C++ includes
#include <algorithm>    // equal_range

// Local includes -------------------------------
#include "sparse_MmatrixM.h"

// petsc ----------------------------------------------
#include "parallelM.h"// #include "parallel.h"
#include "petsc_macroM.h"

EXTERN_C_FOR_PETSC_BEGIN
# include <petscmat.h>
EXTERN_C_FOR_PETSC_END
// Macro to identify and debug functions which should be called in
// parallel on parallel matrices but which may be called in serial on
// serial matrices.  This macro will only be valid inside non-static
// PetscMMatrixM methods
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


// ==============================================
// Rectangular Petsc matrix
// ==============================================
 /**
   * Constructor; initializes the matrix to
   * be empty, without any structure, i.e.
   * the matrix is not usable at all. This
   * constructor is therefore only useful
   * for matrices which are members of a
   * class. All other matrices should be
   * created at a point in the data flow
   * where all necessary information is
   * available.
   *
   * You have to initialize
   * the matrix before usage with
   * \p init(...).
   */

class PetscMMatrixM : public SparseMMatrixM {
  
  private:
   // data ------------------------------------
  Mat _mat;                 ///< Petsc matrix pointer
  bool _destroy_mat_on_exit;///< Boolean value (false)

public:
  public:
  // Constructor ---------------------------------------------------------
  /// Constructor I;  initialize the matrix before usage with \p init(...).
  explicit
  PetscMMatrixM(
    const ParallelM::Communicator &comm CAN_DEFAULT_TO_COMMWORLD
  );   
//   PetscMMatrixM();
  /// Constructor II.
  explicit
  PetscMMatrixM(
    Mat m,                          ///<  matrix Mat m  already exists
    const ParallelM::Communicator &comm CAN_DEFAULT_TO_COMMWORLD
  );  
//   PetscMMatrixM(Mat m);

  // init
   void init(const int  m,const int  n,const int  m_l,const int  n_l, ///< Initialize a Petsc matrix
            const  int nnz=0, const  int noz=0);
  void init();                                            ///< Initialize using sparsity structure
  void init(const int  m,  const int  n) {  _m=m; _n=n; } ///< Initialize using sparsity structure
  // Destructors ----------------------------
  /// Destructor
  ~PetscMMatrixM();

  void clear(); /// Release all memory
  void zero();///< set to zero
  void zero_rows(std::vector<int> & rows, double diag_value = 0.0);///< set  rows to zero
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
  double operator()(const  int i,const  int j) const;
  
  int MatGetRowM(const  int i_val,int *cols, double *vals);

  // Setting -------------------------------------
  ///   sparsity patter update
  void update_sparsity_pattern(const Graph &sparsity_pattern);
  void update_sparsity_pattern(int m,int n,int m_l,int n_l,
     const std::vector< int>  n_oz, const std::vector< int>  n_nz);
  // set values
  void set(const  int i,const  int j, const double value); ///< Set the value.
  void add(const  int i,const  int j,const double value);///<  add  value
  // set flags
  bool closed() const; ///< set closed flag
  /// Add the full matrix to the Petsc matrix.
  void add_matrix(const DenseMatrixM &dm,
                  const std::vector< int> &rows,
                  const std::vector< int> &cols);
  /// Add the full matrix to the Petsc matrix.
  void add_matrix(const DenseMatrixM &dm,
                  const std::vector< int> &dof_indices);

  /// Add a Sparse matrix
  void add(const double a, SparseMMatrixM &X);


  // functions ------------------------------------
  /// Return the l1-norm of the matrix
  double l1_norm() const;
  /// Return the linfty-norm of the matrix
  double linfty_norm() const;

 
    // Setting -------------------------------------

 
  /// print 
    void print(const std::string&/* name*/)const{
      std::cerr << "No implemented";std::abort();     
    };
  /// Print the contents of the matrix to the screen
  virtual  void print(std::ostream& /*os=std::cout*/) const{
      std::cerr << "No implemented";std::abort();     
    };

  /**
   * Print the contents of the matrix to the screen
   * with the PETSc viewer.  This function only allows
   * printing to standard out, this is because we have
   * limited ourselves to one PETSc implementation for
   * writing.
   */
  void print_personal(std::ostream& os=std::cout) const;

  /**
   * Print the contents of the matrix in Matlab's
   * sparse matrix format. Optionally prints the
   * matrix to the file named \p name.  If \p name
   * is not specified it is dumped to the screen.
   */
  void print_matlab(const std::string name="NULL") const;
void print_hdf5(const std::string name) const;
// void read_hdf5(const std::string namefile,const int mode,const int ml_init);
  /**
   * Copies the diagonal part of the matrix into \p dest.
   */
  virtual void get_diagonal (NumericVectorM& dest) const;

  /**
   * Copies the transpose of the matrix into \p dest, which may be
   * *this.
   */
  virtual void get_transpose (SparseMMatrixM& dest) const;

  /**
   * Swaps the raw PETSc matrix context pointers.
   */
  void swap (PetscMMatrixM &);

 

protected:

  /**
   * This function either creates or re-initializes
   * a matrix called "submatrix" which is defined
   * by the row and column indices given in the "rows" and "cols" entries.
   * This function is implemented in terms of the MatGetSubMatrix()
   * routine of PETSc.  The boolean reuse_submatrix parameter determines
   * whether or not PETSc will treat "submatrix" as one which has already
   * been used (had memory allocated) or as a new matrix.
   */
  virtual void _get_submatrix(SparseMMatrixM& submatrix,
			      const std::vector< int>& rows,
			      const std::vector< int>& cols,
			      const bool reuse_submatrix) const;


};




//-----------------------------------------------------------------------
// PetscMMatrixM inline members
inline PetscMMatrixM::PetscMMatrixM(
  const ParallelM::Communicator &comm 
)  : SparseMMatrixM(comm), 
    _destroy_mat_on_exit(true){}
// =================================================================
inline PetscMMatrixM::PetscMMatrixM(
  Mat m,                                    ///< matrix already exists
  const ParallelM::Communicator &comm 
): SparseMMatrixM(comm), 
  _destroy_mat_on_exit(false){
  this->_mat = m;  this->_is_initialized = true;
}
// ===============================================
inline PetscMMatrixM::~PetscMMatrixM(){  this->clear();}
// ==========================================
inline void PetscMMatrixM::close () const{
#ifdef NDEBUG  
  parallel_onlyM();
#endif
  // BSK - 1/19/2004
  // strictly this check should be OK, but it seems to
  // fail on matrix-free matrices.  Do they falsely
  // state they are assembled?  Check with the developers...
//   if (this->closed())
//     return;
  
  int ierr=0;
 
  ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
         CHKERRABORT(this->comm().get(),ierr);
  ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
         CHKERRABORT(this->comm().get(),ierr);
}
// ==================================================
inline  int PetscMMatrixM::m () const{
  assert (this->initialized());
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize (_mat, &petsc_m, &petsc_n); CHKERRABORT(this->comm().get(),ierr);
  return static_cast< int>(petsc_m);
}
// ============================================
inline  int PetscMMatrixM::n () const{
  assert (this->initialized());  
  int petsc_m=0, petsc_n=0, ierr=0;
  ierr = MatGetSize (_mat, &petsc_m, &petsc_n); CHKERRABORT(this->comm().get(),ierr);
  return static_cast< int>(petsc_n);
}
// ==========================================
inline  int PetscMMatrixM::row_start () const{
  assert (this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop);CHKERRABORT(this->comm().get(),ierr);
  return static_cast< int>(start);
}
// =========================================
inline  int PetscMMatrixM::row_stop () const{
  assert (this->initialized());
  int start=0, stop=0, ierr=0;
  ierr = MatGetOwnershipRange(_mat, &start, &stop);
         CHKERRABORT(this->comm().get(),ierr);

  return static_cast< int>(stop);
}



inline void PetscMMatrixM::set (const  int i,
			  const  int j,
			  const double value){  
  assert (this->initialized());
  int ierr=0, i_val=i, j_val=j;
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
		      &petsc_value, INSERT_VALUES);
         CHKERRABORT(this->comm().get(),ierr);
}



// =================================================
inline
void PetscMMatrixM::add (const  int i,
			  const  int j,
			  const double value){
  assert (this->initialized());
  
  int ierr=0, i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
		      &petsc_value, ADD_VALUES);
         CHKERRABORT(this->comm().get(),ierr);
}



// =====================================================
inline void PetscMMatrixM::add_matrix(const DenseMatrixM& dm,
				const std::vector< int>& dof_indices){
  this->add_matrix (dm, dof_indices, dof_indices);
}
// =========================================================
inline void PetscMMatrixM::add (const double a_in, SparseMMatrixM &X_in){
  assert (this->initialized());

  // sanity check. but this cannot avoid 
  // crash due to incompatible sparsity structure...
  assert (this->m() == X_in.m());
  assert (this->n() == X_in.n());

  PetscScalar     a = static_cast<PetscScalar>      (a_in);
  PetscMMatrixM* X = dynamic_cast<PetscMMatrixM*> (&X_in);

  assert (X != NULL);
  
  int ierr=0;

  // the matrix from which we copy the values has to be assembled/closed
  // X->close ();
  assert(X->closed());

  semiparallel_onlyM();

// 2.2.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,3,0)  
  
  ierr = MatAXPY(&a,  X->_mat, _mat, SAME_NONZERO_PATTERN);
         CHKERRABORT(this->comm().get(),ierr);
	 
// 2.3.x & newer
#else
  
  ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
         CHKERRABORT(this->comm().get(),ierr);
	 
#endif
}




// ========================================================
inline double PetscMMatrixM::operator () (const  int i,
			       const  int j) const{
  assert (this->initialized());
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
  // the matrix needs to be closed for this to work
  // this->close();
  // but closing it is a semiparallel operation; we want operator()
  // to run on one processor.
  assert(this->closed());
  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
         CHKERRABORT(this->comm().get(),ierr);
  // Perform a binary search to find the contiguous index in
  // petsc_cols (resp. petsc_row) corresponding to global index j_val
  std::pair<const int*, const int*> p =
    std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);
  // Found an entry for j_val
  if (p.first != p.second)    {
      // The entry in the contiguous row corresponding
      // to the j_val column of interest
      const int j = std::distance (const_cast<int*>(&petsc_cols[0]),
				   const_cast<int*>(p.first));
      
      assert (j < ncols);
      assert (petsc_cols[j] == j_val);
      value = static_cast<double> (petsc_row[j]);
    } 
  ierr  = MatRestoreRow(_mat, i_val,
			&ncols, &petsc_cols, &petsc_row);
          CHKERRABORT(this->comm().get(),ierr);
  return value;
}
// ================================================
inline bool PetscMMatrixM::closed() const{
  assert (this->initialized());
  int ierr=0;
  PetscBool assembled;
  ierr = MatAssembled(_mat, &assembled);
         CHKERRABORT(this->comm().get(),ierr);
  return (assembled == PETSC_TRUE);
}
// =========================================
inline void PetscMMatrixM::swap(PetscMMatrixM &m){
  std::swap(_mat, m._mat);
  std::swap(_destroy_mat_on_exit, m._destroy_mat_on_exit);
}
// ===================================================
inline void PetscMMatrixM::print_personal(std::ostream& os) const{
  assert (this->initialized());
#ifndef NDEBUG
  if (os != std::cout)
    std::cerr << "Warning! PETSc can only print to std::cout!" << std::endl;
#endif
  int ierr=0;
  ierr = MatView(_mat, PETSC_VIEWER_STDOUT_SELF);
         CHKERRABORT(this->comm().get(),ierr);
}


typedef  PetscMMatrixM PetscMMatrixM;


#endif // #ifdef LIBMESH_HAVE_PETSC
#endif // #ifdef __petsc_matrix_h__
