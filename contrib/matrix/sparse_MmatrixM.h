#ifndef __sparse_MmatrixM_h__
#define __sparse_MmatrixM_h__

// C++ includes -------------------
#include <memory>
#include <iostream>

// configure files ----------------
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"  // solver library

// local includes -----------------------
#include "MGGraph.h"        // graph class
#include "SolverPackage_enum.h"// solver package
#include "parallelM.h"
#include "parallel_objectM.h"

// forward declarations ----------------------
class SparseMatrixM;     // sparse matrix
class DenseMatrixM;      // dense matrix
class NumericVectorM;    // vector

#define numeric_index_type int

// ==========================================
// Rectangular sparse matrix.
// =========================================
// ==========================================
//             Generic sparse matrix.
// ===========================================
///   Matrix Data parallel/non parall matrix dimensions
///    _m                        = number of global rows
///    _n                        = number of global columns
///    _m_l                      = number of local rows
///    _n_l                      = number of local columns
///    _ml_start                 = processor start dof
///    bool _is_initialized; the matrix has been initialized (y/n)
///
///   Build from the solver package:
///    ParallelM::Communicator &comm,            parallel communicator
///    SolverPackageM  solver_package            solver
///
///   Matrix init based on
///    int _m;                         number of global rows
///    int _n;                         number of global columns
///    int _m_l;                       number of local rows
///    int _n_l;                       number of local columns

class SparseMMatrixM:
#ifdef LM_REFCOUNT
  public ReferenceCountedObject<SparseMMatrixM >,
#endif
  public ParallelObjectM {
protected:
  // =================
  // Data
  // =================
  // matrix dimensions
  int _m;  ///< number of global rows
  int _n;  ///< number of global columns
  // parallel matrix dimensions
  int _m_l;       ///< number of local rows
  int _n_l;       ///< number of local columns
  int _ml_start;  ///< processor start dof
  // Flags
  bool _is_initialized; ///<  the matrix has been initialized (y/n)
public:
  // ===============================
  // Constructor - Destructor
  // ===============================
  /// Constructor;
  explicit
  SparseMMatrixM(
    const ParallelM::Communicator &comm CAN_DEFAULT_TO_COMMWORLD ///< matrix communicator
  );

  /// Builds a \p SparseMMatrixM using the linear solver package
  /// specified by \p solver_package
  static std::unique_ptr<SparseMMatrixM>  build(
    const ParallelM::Communicator &comm,           ///< matrix communicator
    const SolverPackageM solver_package = LSOLVER  ///< solver package
  );

  // Init functions
  virtual void init() = 0;///< Initialize  matrix 0
  virtual void init(      ///< Initialize  matrix 1
    const int  m,  const int  n           ///< global indices
  ) {    _m=m;    _n=n;  }
  virtual void init(      ///< Initialize  matrix 2
    const int  m,  const int  n,             ///< global indices
    const int  m_l,const int  n_l            ///< lobal indices
  ) {    _m=m;    _n=n;    _m=m_l;    _n=n_l;  }

//   virtual void init(      ///< Initialize  matrix 3
//      const int  m,  const int  n,         ///< global indices
//      const int  m_l,const int  n_l,       ///< local indices
//      const int  nnz=30, const int  noz=10, ///< diag-offset indices
//      const numeric_index_type blocksize=1
//    ) { _m=m; _n=n; _m=m_l; _n=n_l; }

  // Destructor
  virtual ~SparseMMatrixM() {}; ///< Destructor. Free all memory
  virtual void clear() = 0;  ///< Release all memory and return

  // Setting data -------------------------------------
  // setting single values
  virtual void set(
    const int i, const int j,  ///< (i,j)
    const double value         ///< M(i,j) =value
  ) = 0;   ///< Set a value
  virtual void add(
    const int i, const int j, ///< (i,j)
    const double value        ///< M(i,j) +=value
  ) = 0;   ///< Add a value
  // setting zero
  virtual void zero() = 0; /// Set all entries to 0.
  virtual void zero_rows( ///< Set all row entries to 0 and the diagonal entry
    std::vector<numeric_index_type> & rows,   ///< row entries
    double diag_value = 0.0                   ///< row diag values
  );
  // setting flag
  virtual void close() const = 0;  /// set close flag

// Return data -------------------------------------------------
  // matrix values
  virtual double operator()( /// Return the value
    const numeric_index_type i,const numeric_index_type j     ///< (i,j)
  ) const = 0;
//   virtual int MatGetRowM(const int i_val,int *cols,double *vals)=0;      ///< Return the values of a row

  // flag
  virtual bool initialized() const {
    return _is_initialized;    ///< Initialization flag
  }
  virtual bool need_full_sparsity_pattern() const {
    return false;    ///< Full sparsity flag (laspack)
  }
  virtual bool closed() const = 0;                                  ///< Close flag

  // Updates the matrix sparsity pattern
  virtual void update_sparsity_pattern(const Graph &) =0;    ///< Full sparsity update
  virtual void update_sparsity_pattern(                      ///< Partial sparsity update
    int m,int n,                ///< global indices
    int m_l,int n_l,            ///< lobal indices
    const std::vector<int >  n_oz,      ///<  offline dofs
    const std::vector<int >  n_nz       ///<  diag dofs
  ) =0;

  // dimensions
  virtual numeric_index_type m() const = 0; ///  Row-dimension  m return
  virtual numeric_index_type n() const = 0; ///  Column-dimension n return

  virtual numeric_index_type row_start() const = 0; /// return row_start, the index of the first matrix row
  virtual numeric_index_type row_stop() const = 0;  /// return row_stop, the index of the last matrix row (+1)

  // ==========================================
  // Operations 
  // ==========================================
  // Add
  /// Add the full matrix to the Sparse matrix.
  virtual void add_matrix(
    const DenseMatrixM &dm,            ///< dense matrix
    const std::vector<int> &rows,      ///< dof rows
    const std::vector<int> &cols       ///< dof columns
  ) = 0;
  /// Same, but assumes the row and column maps are the same.
  virtual void add_matrix(
    const DenseMatrixM &dm,              ///< dense matrix
    const std::vector<int> &dof_indices  ///< dof indices
  ) = 0;
  /// Add a Sparse matrix \p _X, scaled with \p _a, to \p  A += cB
  virtual void add(
    const double /*c*/, 
    SparseMMatrixM & /*B*/          ///< dense matrix
  ) = 0;

  // norm of the matrix
  virtual double l1_norm() const = 0;      ///< Return the l1-norm
  virtual double linfty_norm() const = 0;  ///< Return the linfty-norm

  /// Multiplies the matrix with  arg and stores the result in dest.
  void vector_mult(NumericVectorM& dest,const NumericVectorM& arg) const;
  /// Multiplies the matrix with  arg and adds the result to  dest.
  void vector_mult_add(NumericVectorM& dest,const NumericVectorM& arg) const;

  // form new matrices
  virtual void get_diagonal(NumericVectorM& dest) const = 0;  ///<  Diagonal
  virtual void get_transpose(SparseMMatrixM& dest) const = 0;  ///< Matrix transpose

  // Print - Read  ---------------------------------------------------
  // print
  void print(std::ostream& /*os=std::cout*/) const {};                          ///< Print inot an hdf5 file
  friend std::ostream& operator << (std::ostream& os, const SparseMMatrixM& m); ///< Print  to the screen
  virtual void print_personal(std::ostream& os=std::cout) const = 0;            ///< Print  to file
  // print into an hdf5 file
  virtual  void print_hdf5(const std::string name) const=0;                     ///< Print  to hdf5 file
  // read
  void read(const std::string& name);                                 ///< Read
  // Read from an hdf5 file
  virtual void read_len_hdf5(const std::string namefile,  const std::string mode,
                             int len_row[],int leng_off_row[]);      ///< Read lengths
  virtual void read_pos_hdf5(const std::string namefile,  const std::string mode    ,
                             int pos[],double val[]);                ///< Read pos
  virtual int read_dim_hdf5(const std::string namefile,const std::string mode, int dim[]);   ///< Read dimensions


// Submatrices ------------------------------------------------------
  /// SubMatrix creation
  virtual void create_submatrix(SparseMMatrixM& submatrix,
                                const std::vector<int>& rows,
                                const std::vector<int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,false); // false means DO NOT REUSE submatrix
  }

  /// SubMatrix reinitialization (quick init  mode same size)
  virtual void reinit_submatrix(SparseMMatrixM& submatrix,
                                const std::vector<int>& rows,
                                const std::vector<int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,true); // true means REUSE submatrix
  }

protected:

  /// Get submatrix from main matrix
  virtual void _get_submatrix(SparseMMatrixM& ,
                              const std::vector<int>& ,
                              const std::vector<int>& ,
                              const bool) const {
    std::cerr << " SparseMMatrixM::_get_submatrix not yet implemented \n";
  }
};


// ==============================
// SparseMMatrixM inline members
// ==============================

inline SparseMMatrixM::SparseMMatrixM(
  const ParallelM::Communicator &comm_in  ///< parallel communicator
) :
  ParallelObjectM(comm_in),
  _is_initialized(false) {}


inline std::ostream& operator << (std::ostream& os, const SparseMMatrixM& m) {
  m.print(os);
  return os;
}


#endif // #ifndef __sparse_matrix_h__
