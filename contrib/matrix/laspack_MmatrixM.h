#ifndef __laspack_MmatrixM_h__
#define __laspack_MmatrixM_h__
#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM

#include "Typedefs_conf.h"
#include "MGGraph.h"

// C++ includes
#include <algorithm>

// Local includes
#include "sparse_MmatrixM.h"
#include "dense_matrixM.h"

#include <matrix.h>


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

class LaspackMMatrixM : public SparseMMatrixM {
public:

  ///   The Laspack sparse matrix pointer.
  Matrix _Mat;
  int *_uncomp_row;
private:
  bool _closed;  /// Flag indicating if the matrix has been closed yet.

public:
  // Constructor ---------------------------------------------------------
  /// Constructor I;  initialize the matrix before usage with \p init(...).
  LaspackMMatrixM():_closed(false){};
  /// Constructor II.
  LaspackMMatrixM(Matrix m);

 
  void init(const uint m,const uint n,const uint m_l,const uint n_l, /// Initialize a Petsc matrix
            const unsigned int nnz=0, const unsigned int noz=0);
  void init(); /// Initialize using sparsity structure

  // Destructors ----------------------------
  /// Destructor
  ~LaspackMMatrixM(){this->clear();};

  void clear();                                        /// Release all memory

  // ==============================
  // DATA RETURN
  // ==============================
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
  // zeroes
  void zero_rows(std::vector<int> & rows, Real diag_value);///< set a row to zero
  void zero();                                            ///< set to zero
  // add
  /// Add the full matrix to the Laspack matrix.
  void add_matrix(const DenseMatrixM &dm,const std::vector<uint> &rows,const std::vector<uint> &cols);
  /// Add the full matrix to the Laspack matrix.
  void add_matrix(const DenseMatrixM &dm,const std::vector<uint> &dof_indices);
  /// Add a Sparse matrix
  void add(const Real a, SparseMMatrixM &X);
  // flags
  void close() const {  const_cast<LaspackMMatrixM*>(this)->_closed = true;;}///< close
    
  // functions ------------------------------------
  // norm
  Real l1_norm() const;     /// Return the l1-norm of the matrix
  Real linfty_norm() const; /// Return the linfty-norm of the matrix

  // print -----------------------------------------------------
  void print_personal(std::ostream& os=std::cout) const{};    ///< print personal
  void print_personal(const std::string name="NULL") const{}; ///< print
  void print_hdf5(const std::string name="NULL") const{};     ///< print hdf5
  friend std::ostream& operator << (std::ostream& os, const SparseMMatrixM& m); /// Same as the print method above
 
  // function
  /// Copies the diagonal part of the matrix
  virtual void get_diagonal(NumericVectorM& dest) const;
  /// Transpose Matrix.
  virtual void get_transpose(SparseMMatrixM& dest) const;
  /// Swaps the raw PETSc matrix context pointers.
  void swap(LaspackMMatrixM &);

private:
  /// This function returns the position in the compressed row. Very expensive
  unsigned int pos(const uint i, const uint j) const;

  /// Make other Laspack datatypes friends
  friend class LaspackVectorM;
  friend class LaspackLinearSolverM;
};
 



// ========================================
// LaspackMMatrixM class inline members
// ==========================================
inline void LaspackMMatrixM::clear() { 
  if (this->initialized())   M_Destr(&_Mat); delete []_uncomp_row;
  _closed = false;  this->_is_initialized = false;
}
// ========================================
inline void LaspackMMatrixM::zero() {
  for (unsigned int row=0; row<this->m(); row++)    {
    const unsigned int len =  M__GetLen(&_Mat, row+1);
    for (unsigned int l=0; l<len; l++)  M__SetEntry(&_Mat,row+1,l,M__GetPos(&_Mat,row+1,l), 0.);
  }
  return;
}
// ==============================================
inline unsigned int LaspackMMatrixM::m() const {
  assert(this->initialized());
  return static_cast<uint>(M_GetRowDim(const_cast<Matrix*>(&_Mat)));
}
// ========================================================
inline unsigned int LaspackMMatrixM::n() const {
  assert(this->initialized());
  return static_cast<unsigned int>(M_GetClmDim(const_cast<Matrix*>(&_Mat)));
}
// ========================================================
inline unsigned int LaspackMMatrixM::row_start() const {  return 0;}
// =======================================================
inline unsigned int LaspackMMatrixM::row_stop() const {  return this->m();}
// ===================================================
inline void LaspackMMatrixM::set(const unsigned int i, const unsigned int j,const Real value) { 
  // very expensive
  assert(this->initialized());  assert(i < this->m());  assert(j < this->n());
  const unsigned int position = this->pos(i,j);
  assert((j+1) == M_GetPos(&_Mat, i+1, position));
  M_SetEntry(&_Mat, i+1, position, j+1, value);
}
// ================================================
inline void LaspackMMatrixM::add(const unsigned int i,
                                const unsigned int j, const Real value) { 
  // very expensive
  assert(this->initialized());
  assert(i < this->m());
  assert(j < this->n());
  const unsigned int position = this->pos(i,j);
  M_AddVal(&_Mat, i+1, position, value);
}
// ===========================================
inline void LaspackMMatrixM::add_matrix(const DenseMatrixM& dm,
                                       const std::vector<unsigned int>& dof_indices) {
  this->add_matrix(dm, dof_indices, dof_indices);
}
// ===========================================
inline void LaspackMMatrixM::add_matrix(const DenseMatrixM& dm,
                                       const std::vector<unsigned int>& rows,
                                       const std::vector<unsigned int>& cols)

{
  assert(this->initialized());  assert(dm.m() == rows.size());  assert(dm.n() == cols.size());
  for (unsigned int i=0; i<rows.size(); i++) {
    int irow=rows[i]+1;
    for (unsigned int k=0; k<M__GetLen(&_Mat,irow); k++)  _uncomp_row[M__GetPos(&_Mat,irow,k)-1]=k;
    for (unsigned int j=0; j<cols.size(); j++)  M__AddVal(&_Mat,irow,_uncomp_row[cols[j]],dm(i,j));
  }
}
// ===========================================
inline void LaspackMMatrixM::add(const Real a_in, SparseMMatrixM &X_in) {
  // the matrices must have the same structure
  assert(this->initialized());  assert(this->m() == X_in.m());
  assert(this->n() == X_in.n());

  LaspackMMatrixM* X = static_cast<LaspackMMatrixM*>(&X_in);
  assert(X != NULL);
  _LPNumber a = static_cast<_LPNumber>(a_in);

  // loops taken from LaspackMMatrixM::zero ()
  const unsigned int n_rows = this->m();
  for (unsigned int row=0; row<n_rows; row++)    {
    const unsigned int len = M__GetLen(&_Mat,row+1);
    for (unsigned int l=0; l<len; l++) {
      const _LPNumber value = a * M__GetVal(const_cast<Matrix*>(&(X->_Mat)), row+1,l);
      M__AddVal(&_Mat, row+1, l, value);
    }
  }
}
// ===========================================
inline Real LaspackMMatrixM::operator()(const uint i,
                                       const uint j) const { 
  // very expensive call
  assert(this->initialized());  assert(i < this->m());  assert(j < this->n());
  return M_GetEl(const_cast<Matrix*>(&_Mat), i+1, j+1);
}
// ===========================================
inline unsigned int LaspackMMatrixM::pos(const uint i,
                                        const uint j) const {
  // very expensive call				  
  assert(i < this->m());  assert(j < this->n());
  const uint length=M__GetLen(&_Mat,i+1);
  for (uint k=0;k<length;k++) if(M__GetPos(&_Mat,i+1,k) ==j) return k;
  std::cout << " LaspackMMatrixM::pos position not found \n";abort();
  return 0;
}


#endif // #ifdef LIBMESH_HAVE_LASPACK
#endif // #ifdef __laspack_matrix_h__
