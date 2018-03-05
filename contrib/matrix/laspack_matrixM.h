#ifndef __laspack_matrixM_h__
#define __laspack_matrixM_h__
#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM

#include "Typedefs_conf.h"
#include "MGGraph.h"

// C++ includes
#include <algorithm>

// Local includes
#include "sparse_matrixM.h"
#include "dense_matrixM.h"

#include <qmatrix.h>


// Forward declarations
class LaspackVectorM;
// template <typename T> class LaspackLinearSolver;

// ===========================================
// Generic laspack matrix. This class contains
// pure virtual members that must be overloaded
// in derived classes.  Using a derived class
// allows for uniform access to laspack matrices
// from various different solver packages in
// different formats.
// ================================================

class LaspackMatrixM : public SparseMatrixM {
public:

  private:
  // data ------------------------------------
  QMatrix _QMat;             ///< Petsc matrix pointer
  bool _destroy_mat_on_exit;///< Boolean value (false)
  int *_uncomp_row;         /// Aux vector to compress and decompress
  bool _closed; ///< Flag closure
  
public:
  // Constructor ---------------------------------------------------------
  /// Constructor I;  initialize the matrix before usage with \p init(...).
  LaspackMatrixM();
  /// Constructor II.
  LaspackMatrixM(QMatrix m);

  /// Initialize a Petsc matrix
  void init(const uint m,const uint n,const uint m_l,const uint n_l,
            const unsigned int nnz=0, const unsigned int noz=0);

  // Destructors ----------------------------
  /// Destructor
  ~LaspackMatrixM();
  void clear();                                        /// Release all memory

 
  // Returns -------------------------------------------
  // matrix dimensions
  unsigned int m() const;          ///< row-dimension
  unsigned int n() const;          ///< column dimension
  unsigned int row_start() const;  ///< row-start
  unsigned int row_stop() const;   ///< row-stop
  // flags
  bool closed() const {return _closed;}                 ///< Return the closed flag value
  bool need_full_sparsity_pattern() const {return true;}///< Return the  sparsity flag value
  // values
  int MatGetRowM(const unsigned int i_val,int *cols, double *vals){}; ///< Return the row values
  Real operator()(const uint i,const uint j) const;     ///< Return the value 
  
  // Setting -------------------------------------
  //  sparsity patter update
  void update_sparsity_pattern(const Graph &sparsity_pattern);///< Full sparsity update
  void update_sparsity_pattern(int m,int n,int m_l,int n_l,const uint ml_start, ///< Partial sparsity update
			       const std::vector<uint>  n_oz,const std::vector<uint>  n_nz){};
  // values
  void set(const unsigned int i,const unsigned int j, const Real value);///< Set the value.
  void add(const unsigned int i,const unsigned int j,const Real value);///<  add  value
  // set zeros
  void zero_rows(std::vector<int> & rows, Real diag_value);///< set a row to zero
  void zero();                                            ///< set to zero
  // adding
  /// Add the full matrix to the Laspack matrix.
  void add_matrix(const DenseMatrixM &dm,const std::vector<uint> &rows,const std::vector<uint> &cols);
  /// Add the full matrix to the Laspack matrix.
  void add_matrix(const DenseMatrixM &dm,const std::vector<uint> &dof_indices);
  /// Add a Sparse matrix
  void add(const Real a, SparseMatrixM &X);
  // flags
  void close() const {const_cast<LaspackMatrixM*>(this)->_closed = true;}///< close
    
  // functions ------------------------------------
  // norm
  Real l1_norm() const;     /// Return the l1-norm of the matrix
  Real linfty_norm() const; /// Return the linfty-norm of the matrix

  // print -----------------------------------------------------
  void print_personal(std::ostream& os=std::cout) const;    ///< print personal
  void print_personal(const std::string name="NULL") const; ///< print
  void print_hdf5(const std::string name="NULL") const;     ///< print hdf5
  friend std::ostream& operator << (std::ostream& os, const SparseMatrixM& m); /// Same as the print method above
 
  // function
  /// Copies the diagonal part of the matrix
  virtual void get_diagonal(NumericVectorM& dest) const;
  /// Transpose Matrix.
  virtual void get_transpose(SparseMatrixM& dest) const;
  /// Swaps the raw PETSc matrix context pointers.
  void swap(LaspackMatrixM &);

private:
  /// @returns the position in the compressed row. This is an expensive operation
  unsigned int pos(const uint i, const uint j) const;

  /// Make other Laspack datatypes friends
  friend class LaspackVectorM;
  friend class LaspackLinearSolverM;
};



// ========================================
// LaspackMatrixM class inline members
// ========================================
inline LaspackMatrixM::LaspackMatrixM():_closed(false) {}
// ========================================
inline LaspackMatrixM::~LaspackMatrixM() {this->clear();}
// ==========================================
inline void LaspackMatrixM::clear() {
  if (this->initialized())   Q_Destr(&_QMat);
  delete []_uncomp_row;
  _closed = false;
  this->_is_initialized = false;
}
// ========================================
inline void LaspackMatrixM::zero() {
  for (unsigned int row=0; row<this->m(); row++) {
    const unsigned int len =  Q__GetLen(&_QMat, row+1);
    for (unsigned int l=0; l<len; l++)  Q__SetEntry(&_QMat,row+1,l,Q__GetPos(&_QMat,row+1,l), 0.);
  }
  return;
}
// ==============================================
inline unsigned int LaspackMatrixM::m() const {
  assert(this->initialized());
  return static_cast<uint>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}
// ========================================================
inline unsigned int LaspackMatrixM::n() const {
  assert(this->initialized());
  return static_cast<unsigned int>(Q_GetDim(const_cast<QMatrix*>(&_QMat)));
}
// ========================================================
inline unsigned int LaspackMatrixM::row_start() const {  return 0;}
// =======================================================
inline unsigned int LaspackMatrixM::row_stop() const {  return this->m();}
// ===================================================
inline void LaspackMatrixM::set(const uint i, const uint j,const Real value) {
  // very expensive
  assert(this->initialized());  assert(i < this->m());
  assert(j < this->n());
  const unsigned int position = this->pos(i,j);
  assert((j+1) == Q_GetPos(&_QMat, i+1, position));
  Q_SetEntry(&_QMat, i+1, position, j+1, value);
}
// ================================================
inline void LaspackMatrixM::add(const uint i,const uint j, const Real value) {
  // very expensive
  assert(this->initialized());  assert(i < this->m());  assert(j < this->n());
  const unsigned int position = this->pos(i,j);
  Q_AddVal(&_QMat, i+1, position, value);
}
// ===========================================
inline void LaspackMatrixM::add_matrix(const DenseMatrixM& dm,
                                       const std::vector<uint>& dof_indices) {
  this->add_matrix(dm, dof_indices, dof_indices);
}
// ===========================================
inline void LaspackMatrixM::add_matrix(const DenseMatrixM& dm,
                                       const std::vector<uint>& rows,
                                       const std::vector<uint>& cols) {
  assert(this->initialized());  assert(dm.m() == rows.size());
  assert(dm.n() == cols.size());
  for (unsigned int i=0; i<rows.size(); i++) {
    int irow=rows[i]+1;
    for (unsigned int k=0; k<Q__GetLen(&_QMat,irow); k++)  _uncomp_row[Q__GetPos(&_QMat,irow,k)-1]=k;
    for (unsigned int j=0; j<cols.size(); j++) Q__AddVal(&_QMat,irow,_uncomp_row[cols[j]],dm(i,j));
  }
}
// ===========================================
inline void LaspackMatrixM::add(const Real a_in, SparseMatrixM &X_in) {
  // the matrices must have the same structure
  assert(this->initialized());  assert(this->m() == X_in.m());
  assert(this->n() == X_in.n());

  LaspackMatrixM* X = dynamic_cast<LaspackMatrixM*>(&X_in);
  assert(X != NULL);
  _LPNumber         a = static_cast<_LPNumber>(a_in);

  // loops taken from LaspackMatrixM::zero ()
  const unsigned int n_rows = this->m();
  for (unsigned int row=0; row<n_rows; row++)    {
    const unsigned int len = Q__GetLen(&_QMat,row+1);
    for (unsigned int l=0; l<len; l++) {
      const _LPNumber value = a * Q__GetVal(const_cast<QMatrix*>(&(X->_QMat)), row+1,l);
      Q__AddVal(&_QMat, row+1, l, value);
    }
  }
}
// ===========================================
inline Real LaspackMatrixM::operator()(const uint i,const uint j) const {
  // very expensive call
  assert(this->initialized());    assert(i < this->m());    assert(j < this->n());
  return Q_GetEl(const_cast<QMatrix*>(&_QMat), i+1, j+1);
}
// ===========================================
inline unsigned int LaspackMatrixM::pos(const uint i,const uint j) const {
  // very expensive call
  assert(i < this->m());  assert(j < this->n());
  const uint length=Q__GetLen(&_QMat,i+1);
  for (uint k=0;k<length;k++) if (Q__GetPos(&_QMat,i+1,k) ==j) return k;
  std::cout << " LaspackMatrixM::pos position not found \n";
  abort();  return 0;
}


#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifdef __laspack_matrix_h__
