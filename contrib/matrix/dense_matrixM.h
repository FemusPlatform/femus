#ifndef __dense_matrixM_h__
#define __dense_matrixM_h__

// C++ includes
#include <vector>
#include <cmath>
#include <algorithm>

// Local Includes
#include "dense_matrix_baseM.h"

// Forward Declarations
class DenseVectorM;


// #define MATRIX_HAVE_PETSC 1

// ===========================================================
// Defines a dense matrix for use in Finite Element-type computations.
// Useful for storing element stiffness matrices before summation
// into a global matrix.
// ==============================================================

// ------------------------------------------------------------
// Dense Matrix class definition
class DenseMatrixM : public DenseMatrixBaseM
{
  
  // ===============================
  // Data
  // =============================
public:  
  /// Run-time selectable option to turn on/off blas support.
  bool use_blas;
private:
  /// The actual data values, stored as a 1D array.
  std::vector<double> _val;
  
   ///The decomposition schemes above change the entries of the matrix
  enum DecompositionType {LU=0, CHOLESKY=1, NONE};
  /// This flag keeps track of which type of decomposition has been
  DecompositionType _decomposition_type;

public:
  // ==================================
   // Construction/Destr
  // ===============================
  /// Constructor.  Creates a dense matrix of dimension \p m by \p n.
  DenseMatrixM(const int m=0,const int n=0);
  /// Copy-constructor.
  //DenseMatrixM (const DenseMatrixM& other_matrix);
  /// Destructor.  Empty.     
  virtual ~DenseMatrixM() {}
 
  /// Set every element in the matrix to 0.
  virtual void zero();
  /// Resize the matrix.  Will never free memory, but may allocate more.  Sets all elements to 0.
  void resize(const int m, const int n);
  
  // ===============================
  // Return functions
  // ===============================
  /// @returns the \p (i,j) element of the matrix.
  double operator() (const int i,
		const int j) const;
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  double & operator() (const int i,
		  const int j);
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual double el(const int i,
	       const int j) const { return (*this)(i,j); }
  /// @returns the \p (i,j) element of the matrix as a writeable reference.
  virtual double & el(const int i,const int j)     { return (*this)(i,j); } 
  
  /// Access to the values array
  std::vector <double>& get_values() { return _val; }

  /// Return a constant reference to the matrix values.
  const std::vector <double>& get_values() const { return _val; }

  // =============================== 
  // Setting function
  // ===============================
  ///Put the \p sub_m x \p sub_n principal submatrix into \p dest.
  void get_principal_submatrix (int sub_m, int sub_n, DenseMatrixM& dest) const;
  /// Put the \p sub_m x \p sub_m principal submatrix into \p dest.
  void get_principal_submatrix (int sub_m, DenseMatrixM& dest) const;

  
  // ===============================
  // Algebra
  // ===============================	
  /// Left multipliess by the matrix \p M2.
  virtual void left_multiply (const DenseMatrixBaseM& M2);
  /// Right multiplies by the matrix \p M3.
  virtual void right_multiply (const DenseMatrixBaseM& M3);

  /// Perform matrix vector multiplication.
  void vector_mult(DenseVectorM& dest, const DenseVectorM& arg) const;
  /// Perform matrix vector multiplication and add scaled result to \p dest.
  void vector_mult_add (DenseVectorM& dest,const double factor,const DenseVectorM& arg) const;

  /// Assignment operator
  DenseMatrixM& operator = (const DenseMatrixM& other_matrix);
  /// STL-like swap method
  void swap(DenseMatrixM& other_matrix);
 
  /// Multiplies every element in the matrix by \p factor.
  void scale (const double factor);
  /// Multiplies every element in the matrix by \p factor.
  DenseMatrixM& operator *= (const double factor);
  
  /// Adds \p mat to this matrix.
  DenseMatrixM& operator+= (const DenseMatrixM &mat);
  /// Adds \p factor times \p mat to this matrix.
  void add (const double factor, const DenseMatrixM& mat);
  
  /// Tests if \p mat is exactly equal to this matrix.
  bool operator== (const DenseMatrixM &mat) const;
  /// Tests if \p mat is not exactly equal to this matrix.
  bool operator!= (const DenseMatrixM &mat) const;


  /// Subtracts \p mat from this matrix.
  DenseMatrixM& operator-= (const DenseMatrixM &mat);

  /// this returns the minimum
  double min () const;
  /// @returns the maximum element in the matrix.
  double max () const;

  /// Return the l1-norm of the matrix, that is max. sum of columns
  double l1_norm () const;
  /// Return the linfty-norm of the max. sum of rows.
  double linfty_norm () const;

  /// Left multiplies by the transpose of the matrix \p A.
  void left_multiply_transpose (const DenseMatrixM& A);
  /// Right multiplies by the transpose of the matrix \p A
  void right_multiply_transpose (const DenseMatrixM& A);
  
  /// @returns the \p (i,j) element of the transposed matrix.
  double transpose (const int i,
	       const int j) const;
    /// Put the tranposed matrix into \p dest.
  void get_transpose(DenseMatrixM& dest) const;
 

  // =========================================
  //  SOLVE
  // ========================================
  
  /**
   * Condense-out the \p (i,j) entry of the matrix, forcing
   * it to take on the value \p val.  This is useful in numerical
   * simulations for applying boundary conditions.  Preserves the
   * symmetry of the matrix.
   !!!!!!!!!!!!!!!!!!!!!
   TODO
   DenseVectorBase -> DenseVector
   */
  void condense(const int i,
		const int j,
		const double val,
		DenseVectorBaseM& rhs)
  { DenseMatrixBaseM::condense (i, j, val, rhs); }

  /// Solve the system Ax=b given the input vector b.
  void lu_solve (DenseVectorM& b,DenseVectorM& x,const bool partial_pivot = false);
  /// A Cholesky factorizationof A such that A = L L^T   
  void cholesky_solve(DenseVectorM& b,DenseVectorM& x);

  ///   * @returns the determinant of the matrix.  
  double det();

  /// Computes the inverse of the dense matrix (assuming it is invertible)
  // void inverse();


private:

  /// Form the LU decomposition of the matrix.  
  void _lu_decompose (const bool partial_pivot = false);
  /// Solves the system Ax=b through back substitution.  
  void _lu_back_substitute (DenseVectorM& b,DenseVectorM& x,const bool partial_pivot = false) const;
  
  ///Decomposes a symmetric positive definite matrix 
  void _cholesky_decompose();
  /// Solves the equation Ax=b for the unknown value x and rhs b based on the Cholesky factorization 
  void _cholesky_back_substitute(DenseVectorM& b, DenseVectorM& x) const;

  // ======================
  // BLAS operations
  // ======================
  /// Computes A <- op(A) * op(B) using BLAS gemm function.
  // Used in the right_multiply(), left_multiply(), right_multiply_transpose(),
  // and left_multiply_tranpose() routines.
  enum _BLAS_Multiply_Flag {
    LEFT_MULTIPLY = 0,
    RIGHT_MULTIPLY,
    LEFT_MULTIPLY_TRANSPOSE,
    RIGHT_MULTIPLY_TRANSPOSE
  };	
   /// BLAS operations
  void _multiply_blas(const DenseMatrixBaseM& other,
		      _BLAS_Multiply_Flag flag);
		      

  
  
};












// ------------------------------------------------------------
// Dense Matrix member functions

inline
DenseMatrixM::DenseMatrixM(const int m,const int n)
  : DenseMatrixBaseM(m,n),
#ifdef MATRIX_HAVE_PETSC
    use_blas(true),
#else
    use_blas(false),
#endif
    _decomposition_type(NONE){
  this->resize(m,n);
}



// FIXME[JWP]: This copy ctor has not been maintained along with
// the rest of the class...
// Can we just use the compiler-generated copy ctor here?
// 
// inline
// DenseMatrixM::DenseMatrixM (const DenseMatrixM& other_matrix)
//   : DenseMatrixMBaseM <double>(other_matrix._m, other_matrix._n)
// {
//   _val = other_matrix._val;
// }



// =========================================
inline
void DenseMatrixM::swap(DenseMatrixM& other_matrix){
  std::swap(this->_m, other_matrix._m);
  std::swap(this->_n, other_matrix._n);
  _val.swap(other_matrix._val);
  DecompositionType _temp = _decomposition_type;
  _decomposition_type = other_matrix._decomposition_type;
  other_matrix._decomposition_type = _temp;
}


// ===============================================  
inline
void DenseMatrixM::resize(const int m,const int n)
{
  _val.resize(m*n);
  this->_m = m;  this->_n = n;
  _decomposition_type = NONE;
  this->zero();
}

// ==========================================
inline
void DenseMatrixM::zero(){
  _decomposition_type = NONE;
  std::fill (_val.begin(), _val.end(), 0.);
}

// =============================================
inline
DenseMatrixM& DenseMatrixM::operator = (const DenseMatrixM& other_matrix)
{
  this->_m = other_matrix._m; this->_n = other_matrix._n;
  _val                = other_matrix._val;
  _decomposition_type = other_matrix._decomposition_type;
  return *this;
}

// =========================================
inline
double DenseMatrixM::operator () (const int i,
			       const int j) const
{
  assert (i*j<(int)_val.size());assert (i < this->_m);assert (j < this->_n);
  //  return _val[(i) + (this->_m)*(j)]; // col-major
  return _val[(i)*(this->_n) + (j)]; // row-major
}

// ========================================================
inline
double & DenseMatrixM::operator () (const int i,const int j)
{
  assert (i*j< (int)_val.size());assert (i < this->_m); assert (j < this->_n);
  //return _val[(i) + (this->_m)*(j)]; // col-major
  return _val[(i)*(this->_n) + (j)]; // row-major
}

// =================================================
inline
void DenseMatrixM::scale (const double factor)
{
  for (int i=0; i<(int)_val.size(); i++) _val[i] *= factor;
}

// ===============================================
inline
DenseMatrixM& DenseMatrixM::operator *= (const double factor)
{
  this->scale(factor);
  return *this;
}

// ===================================================
inline void DenseMatrixM::add (const double factor, const DenseMatrixM& mat)
{
  for (int i=0; i<(int)_val.size(); i++)  _val[i] += factor * mat._val[i];
}

// =================================================
inline bool DenseMatrixM::operator == (const DenseMatrixM &mat) const
{
  for (int i=0; i<(int)_val.size(); i++)  if (_val[i] != mat._val[i])  return false;
  return true;
}

// ==================================================
inline bool DenseMatrixM::operator != (const DenseMatrixM &mat) const
{
  for (int i=0; i<(int)_val.size(); i++)  if (_val[i] != mat._val[i])  return true;
  return false;
}

// ========================================
inline DenseMatrixM& DenseMatrixM::operator += (const DenseMatrixM &mat)
{
  for (int i=0; i<(int)_val.size(); i++)  _val[i] += mat._val[i];
  return *this;
}

// =======================================================
inline
DenseMatrixM& DenseMatrixM::operator -= (const DenseMatrixM &mat)
{
  for (int i=0; i<(int)_val.size(); i++)   _val[i] -= mat._val[i];
  return *this;
}

// ================================================
inline double DenseMatrixM::min () const
{
  assert (this->_m);  assert (this->_n);
  double my_min = 0.;

  for (int i=0; i!=this->_m; i++)   {
      for (int j=0; j!=this->_n; j++)   {
          double current = (double)((*this)(i,j));
          my_min = (my_min < current? my_min : current);
        }
    }
  return my_min;
}

// ==============================================
inline  double DenseMatrixM::max () const 
{
  assert (this->_m); assert (this->_n);
  double my_max = (double)((*this)(0,0));
  for (int i=0; i!=this->_m; i++)   {
      for (int j=0; j!=this->_n; j++)  {
          double current = (double)((*this)(i,j));
          my_max = (my_max > current? my_max : current);
        }
    }
  return my_max;
}

// =============================================
inline double DenseMatrixM::l1_norm () const
{
  assert (this->_m);assert (this->_n);
  double columnsum = 0.;
  for (int i=0; i!=this->_m; i++)   columnsum += fabs((*this)(i,0));
  double my_max = columnsum;
  for (int j=1; j!=this->_n; j++){
    columnsum = 0.;
      for (int i=0; i!=this->_m; i++) columnsum += fabs((*this)(i,j));
      my_max = (my_max > columnsum? my_max : columnsum);
    }
  return my_max;
}

// =================================================
inline double DenseMatrixM::linfty_norm () const
{
  assert (this->_m);assert (this->_n);
  double rowsum = 0.;
  for (int j=0; j!=this->_n; j++)  rowsum += std::abs((*this)(0,j));
  double my_max = rowsum;
  for (int i=1; i!=this->_m; i++)    {
      rowsum = 0.;
      for (int j=0; j!=this->_n; j++)  rowsum += std::abs((*this)(i,j));
      my_max = (my_max > rowsum? my_max : rowsum);
    }
  return my_max;
}

// ================================================
inline double DenseMatrixM::transpose (const int i,const int j) const{
  // Implement in terms of operator()
  return (*this)(j,i);
}

// =======================================================
// inline void DenseMatrixM::condense(const int iv,
// 			      const int jv,
// 			      const double val,
// 			      DenseVector <double>& rhs)
// {
//   assert (this->_m == rhs.size()); assert (iv == jv);
// 
//   // move the known value into the RHS
//   // and zero the column
//   for (int i=0; i<this->m(); i++)  {
//       rhs(i) -= ((*this)(i,jv))*val;
//       (*this)(i,jv) = 0.;
//     }
// 
//   // zero the row
//   for (int j=0; j<this->n(); j++)  (*this)(iv,j) = 0.;
//   (*this)(iv,jv) = 1.;  rhs(iv) = val;
//   return;
// }


#endif // #ifndef __dense_matrix_h__

