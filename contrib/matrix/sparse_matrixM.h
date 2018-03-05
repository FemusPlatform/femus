#ifndef __sparse_matrixM_h__
#define __sparse_matrixM_h__

// mydef
#define numeric_index_type  int

// C++ includes -------------------
#include <memory>     // auto pointer
#include <iostream>   // std::cout

// configure files ----------------
#include "Printinfo_conf.h"
#include "Solverlib_conf.h"  // solver library

// local includes ----------------------------
#include "SolverPackage_enum.h"// solver package
#include "MGGraph.h"          // graph class
#include "parallelM.h"          // parallel class
#include "parallel_objectM.h"          // parallel class

// forward declarations ----------------------
class SparseMatrixM;     // sparse matrix
class DenseMatrixM;      // dense matrix
class NumericVectorM;    // vector

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



class SparseMatrixM :
#ifdef LM_REFCOUNT
 public ReferenceCountedObject<SparseMatrixM >,
#endif
public ParallelObjectM

{
// =================================
//  Data
// ================================
protected:
    // =================
    // Data
    // =================
    // matrix dimensions
    int _m;                        ///< number of global rows
    int _n;                        ///< number of global columns
    // parallel matrix dimensions
    int _m_l;                      ///< number of local rows
    int _n_l;                      ///< number of local columns
    int _ml_start;                 ///< processor start dof
    // Flags
    bool _is_initialized; ///<  the matrix has been initialized (y/n)   
 
// =================================
//  Constructor - Destructor 
// ================================    
    public:
      explicit
    SparseMatrixM(
        const ParallelM::Communicator &comm CAN_DEFAULT_TO_COMMWORLD ///< matrix communicator
    );

  /// Builds a \p SparseMatrixM using the linear solver package 
  /// specified by \p solver_package
  static std::unique_ptr<SparseMatrixM>  build(
    const ParallelM::Communicator &comm,            ///< parallel communicator ->
    const SolverPackageM solver_package = LSOLVER   ///< linear solver ->
  );
  
  //   SparseMatrixM(): _is_initialized(false){}
  // Init functions
  virtual void init() = 0;///< Initialize  matrix 0
  virtual void init(      ///< Initialize  matrix 1
    const int  m,  const int  n           ///< global indices
  ) {  _m=m; _n=n; }
  virtual void init(      ///< Initialize  matrix 2
    const int  m,  const int  n,          ///< global indices
    const int  m_l,const int  n_l         ///< local indices
  ) { _m=m; _n=n; _m=m_l; _n=n_l; }
//   virtual void init(      ///< Initialize  matrix 3
//      const int  m,  const int  n,         ///< global indices
//      const int  m_l,const int  n_l,       ///< local indices
//      const int  nnz=30, const int  noz=10, ///< diag-offset indices
//      const numeric_index_type blocksize=1
//    ) { _m=m; _n=n; _m=m_l; _n=n_l; }
                                                                
      
  // Destructors
  virtual ~SparseMatrixM (){} ///< Destructor. Free all memory
  virtual void clear () = 0;  ///< Release all memory and return
  
  // Setting data -------------------------------------
  // setting single values 
  virtual void set (const  int i, const  int j,const double value) = 0;  ///< Set a value 
  virtual void add (const  int i, const  int j,const double value) = 0;  ///< Add a value
  // setting zero
  virtual void zero () = 0;/// Set all entries to 0.
  virtual void zero_rows (std::vector<int> & rows, double diag_value = 0.0);///< Set all row entries to 0 and the diagonal entry
  
  // setting flag
  virtual void close () const = 0; /// set close flag
  
  // Return data -------------------------------------------------
  // matrix values
  virtual double operator () (const  int i,const  int j) const = 0; /// Return the value 
  virtual int MatGetRowM(const  int i_val,int *cols=NULL,double *vals=NULL)=0; ///< Return the row values
  
  // flag
  virtual bool initialized() const { return _is_initialized; }      ///< Initialization flag
  virtual bool need_full_sparsity_pattern() const { return false; } ///< Full sparsity flag (laspack) 
  virtual bool closed() const = 0;                                  ///< Close flag
  
  // Updates the matrix sparsity pattern
  virtual void update_sparsity_pattern (const Graph &) =0;                                  ///< Full sparsity update
  virtual void update_sparsity_pattern(int m,int n, int m_l,int n_l,		   
                        const std::vector<int >  n_oz, const std::vector<int >  n_nz  ) =0; ///< Partial sparsity update 

  // dimensions    
  virtual  int m () const = 0;///  Row-dimension  m return
  virtual  int n () const = 0;///  Column-dimension n return  

  virtual  int row_start () const = 0;/// return row_start, the index of the first matrix row 
  virtual  int row_stop () const = 0; /// return row_stop, the index of the last matrix row (+1) 

// =================================
//  Operations
// ================================   
  // Add
  /// Add the full matrix to the Sparse matrix.  
  virtual void add_matrix (const DenseMatrixM &dm,
			   const std::vector< int> &rows,
			   const std::vector< int> &cols) = 0;
  /// Same, but assumes the row and column maps are the same.
  virtual void add_matrix (const DenseMatrixM &dm,const std::vector< int> &dof_indices) = 0;
  /// Add a Sparse matrix \p _X, scaled with \p _a, to \p  A += cB
  virtual void add (const double /*c*/, SparseMatrixM & /*B*/) = 0;

  // norm of the matrix
  virtual double l1_norm () const = 0;     ///< Return the l1-norm
  virtual double linfty_norm () const = 0; ///< Return the linfty-norm
  
  /// Multiplies the matrix with \p arg and stores the result in \p dest.
  void vector_mult (NumericVectorM& dest,const NumericVectorM& arg) const;
  ///Multiplies the matrix with \p arg and adds the result to \p dest.
  void vector_mult_add (NumericVectorM& dest,const NumericVectorM& arg) const;

  // form new matrices
  virtual void get_diagonal (NumericVectorM& dest) const = 0;///<  Diagonal
  virtual void get_transpose (SparseMatrixM& dest) const = 0; ///< Matrix transpose 
  
// =================================
//  Read - Print 
// ================================ 
  // print
  virtual void print(const std::string& name)const;                           ///< Print  to file
  virtual  void print(std::ostream& os=std::cout) const{print_personal(os);}  ///< Print  to the screen
  friend std::ostream& operator << (std::ostream& os, const SparseMatrixM& m);///< Print  to the screen
  virtual void print_personal(std::ostream& os=std::cout) const = 0;          ///< Print  to the text file 
  // print to hdf5 files 
  virtual void print_hdf5(const std::string name="NULL") const=0;             ///< Print  an hdf5 file
  
   // Read 
  void read(const std::string& name);                                ///< Read
  // Read from hdf5 files 
  virtual void read_len_hdf5(const std::string namefile,std::string mode ,
			      int len_row[],int leng_off_row[]);      ///< Read lengths
  virtual void read_pos_hdf5(const std::string namefile,const int mode,
			      int pos[]);                             ///< Read pos
  virtual int read_dim_hdf5(const std::string namefile, // file name 				
				  const std::string mode,    // linear or quadratic
				  int ldim[]);   ///< Read dimensions 	      
// =================================
//  Submatrices
// ================================ 
  /// SubMatrix creation
  virtual void create_submatrix(SparseMatrixM& submatrix,
                                const std::vector< int>& rows,
                                const std::vector< int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,false); // false means DO NOT REUSE submatrix
  }

  /// SubMatrix reinitialization (quick init  mode same size)
  virtual void reinit_submatrix(SparseMatrixM& submatrix,
                                const std::vector< int>& rows,
                                const std::vector< int>& cols) const {
    this->_get_submatrix(submatrix,rows,cols,true); // true means REUSE submatrix
  }

// protected:

  /// Get submatrix from main matrix
  virtual void _get_submatrix(SparseMatrixM& ,
                              const std::vector< int>& ,
                              const std::vector< int>& ,
                              const bool) const {
    std::cerr << " SparseMatrixM::_get_submatrix not yet implemented \n";
  }
 
};


// ==============================
// SparseMatrixM inline members
// ==============================
inline std::ostream& operator << (std::ostream& os, const SparseMatrixM& m){m.print(os);return os;}


#endif // #ifndef __sparse_matrix_h__
