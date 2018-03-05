#include "Solverlib_conf.h"

#ifdef HAVE_PETSCM

// Local includes
#include "petsc_macroM.h"
#include "petsc_MmatrixM.h"
// #include "dof_map.h"
#include "dense_matrixM.h"
#include "petsc_vectorM.h"
#include <sstream>
// hdf5
#include "hdf5.h"


//-----------------------------------------------------------------------
// PetscMatrix members
void PetscMMatrixM::init (const  int m,
                         const  int n,
                         const  int m_l,
                         const  int n_l,
                         const  int /*nnz*/,
                         const  int /*noz*/) {
    _m=m; _n=n; _m_l=m_l; _n_l=n_l;
//  
//     // We allow 0x0 matrices now
//     //if ((m==0) || (n==0))
//     //  return;
// 
//     // Clear initialized matrices
//     if (this->initialized())  this->clear();
//     this->_is_initialized = true;
// 
// 
//     int ierr     = 0;
//     int m_global = static_cast<int>(m);
//     int n_global = static_cast<int>(n);
//     int m_local  = static_cast<int>(m_l);
//     int n_local  = static_cast<int>(n_l);
//     int n_nz     = static_cast<int>(nnz);
//     int n_oz     = static_cast<int>(noz);

    // create a sequential matrix on one processor
//     if (libMesh::n_processors() == 1) {
//         assert ((m_l == m) && (n_l == n));
// 
//         // Create matrix.  Revisit later to do preallocation and make more efficient
//         ierr = MatCreateSeqAIJ (this->comm().get(), m_global, n_global,
//                                 n_nz, PETSC_NULL, &_mat);
//         CHKERRABORT(this->comm().get(),ierr);
//         ierr = MatSetFromOptions (_mat);
//         CHKERRABORT(this->comm().get(),ierr);
//     }
// 
//     else{
//         parallel_only();
//         ierr = MatCreateMPIAIJ (this->comm().get(), m_local, n_local, m_global, n_global,
//                                 n_nz, PETSC_NULL, n_oz, PETSC_NULL, &_mat);
//         CHKERRABORT(this->comm().get(),ierr);
//         ierr = MatSetFromOptions (_mat);
//         CHKERRABORT(this->comm().get(),ierr);
//     }
//     this->zero ();
}

// =====================================0
void PetscMMatrixM::init () {
    std::cout << "  PetscMMatrixM::init: Our PetscMatrix  needs dimensions";
    abort();
}
// =====================================0

// =====================================0
void PetscMMatrixM::update_sparsity_pattern(
                        int m_global,                          // # global rows                           
                        int n_global,                          // # global columns 
                        int m_local,                           // # local rows (local proc)
                        int n_local,                           // # local columns (local proc)	   
                        const std::vector< int>  n_oz, // # diagoanl entries 
                        const std::vector< int>  n_nz  // # offset entries
                        ) {

  // Clear initialized matrices
  if (this->initialized())    this->clear();  this->_is_initialized = true;
  assert((int)n_nz.size()== m_local);assert ((int)n_oz.size()== m_local);
  
//   // processor info
//   int proc_id;    MPI_Comm_rank(this->comm().get(), &proc_id);
//   int numprocs; MPI_Comm_size(this->comm().get(),&numprocs);
   int ierr     = 0;
  
  ierr = MatCreate(this->comm().get(), &_mat);
  CHKERRABORT(this->comm().get(),ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
  CHKERRABORT(this->comm().get(),ierr);
  

  // create a sequential matrix on one processor -------------------
//   if (numprocs == 1) {
//     assert((m_local == m_global) && (n_local == n_global));
//     if (n_nz.empty())
//       ierr = MatCreateSeqAIJ(this->comm().get(), m_global, n_global,
//                              PETSC_NULL, (int*) PETSC_NULL, &_mat);
//     else
//       ierr = MatCreateSeqAIJ(this->comm().get(), m_global, n_global,
//                              PETSC_NULL, (int*) &n_nz[0], &_mat);
//     CHKERRABORT(this->comm().get(),ierr);
// 
//     ierr = MatSetFromOptions(_mat);     CHKERRABORT(this->comm().get(),ierr);
//   } else { // multi processors ----------------------------
//     parallel_onlyM();     //TODO
//     if (n_nz.empty()) {

//old piece of code with Petsc lib version 3.1       
//       ierr = MatCreateMPIAIJ(this->comm().get(),
//                              m_local, n_local,
//                              m_global, n_global,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              PETSC_NULL, (int*) PETSC_NULL,
//                              &_mat
//                             );
    
//     ierr = MatCreate(this->comm().get(), &_mat);
//     CHKERRABORT(this->comm().get(),ierr);
//      
//     ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
//     CHKERRABORT(this->comm().get(),ierr);
//      
//     ierr = MatSetType(_mat, MATMPIAIJ); 
//     CHKERRABORT(this->comm().get(),ierr);
// 
//     ierr = MatMPIAIJSetPreallocation(_mat, 0, 0, 0, 0);
//     CHKERRABORT(this->comm().get(),ierr);
//       
//     }   else {
      
  //old piece of code with Petsc lib version 3.1     
//       ierr = MatCreateMPIAIJ(this->comm().get(),
//                              m_local, n_local,
//                              m_global, n_global,
//                              PETSC_NULL, (int*) &n_nz[0],
//                              PETSC_NULL, (int*) &n_oz[0],
//                              &_mat
//                             );
    
    
//      ierr = MatCreate(this->comm().get(), &_mat);
//      CHKERRABORT(this->comm().get(),ierr);
//      
//      ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
//      CHKERRABORT(this->comm().get(),ierr);
//      
//      ierr = MatSetType(_mat, MATMPIAIJ); // Automatically chooses seqaij or mpiaij
//      CHKERRABORT(this->comm().get(),ierr);
// 
//      ierr = MatMPIAIJSetPreallocation(_mat, 0, (int*)&n_nz[0], 0, (int*)&n_oz[0]);
//      CHKERRABORT(this->comm().get(),ierr);
    
//      } 
     
      ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
      CHKERRABORT(this->comm().get(),ierr);

      ierr = MatSeqAIJSetPreallocation(_mat, 0, (PetscInt*)(n_nz.empty()?NULL:&n_nz[0]));
       CHKERRABORT(this->comm().get(),ierr);
      ierr = MatMPIAIJSetPreallocation(_mat,
				       0, (PetscInt*)(n_nz.empty()?NULL:&n_nz[0]),
				       0, (PetscInt*)(n_oz.empty()?NULL:&n_oz[0]));
       CHKERRABORT(this->comm().get(),ierr);
      
    
     
    ierr = MatSetOptionsPrefix(_mat, "");  CHKERRABORT(this->comm().get(),ierr);
    ierr = MatSetFromOptions(_mat);  CHKERRABORT(this->comm().get(),ierr);
    
//   } // ----------------------------------------------------

  this->zero();
  return;
}

// =====================================0
void PetscMMatrixM::update_sparsity_pattern(const Graph &sparsity_pattern) {

  // dimension
  const  int m   =  sparsity_pattern._m; // global rows;
  const  int n   =  sparsity_pattern._n; // global columns
  const  int n_l =sparsity_pattern._nl;  // local rows
  const  int m_l =sparsity_pattern._ml;  // local columns
  const  int ml_start =sparsity_pattern._ml_start;

  // vectors n_nz (diagonal) and n_oz (offset) -------------------------------
  std::vector< int> n_nz(m_l);   std::vector< int> n_oz(m_l);
  for (int i=0;i<m_l;i++) {
    int len=sparsity_pattern[ml_start+i].size()-1;
    n_oz[i]=sparsity_pattern[ml_start+i][len];
    n_nz[i]=len-n_oz[i];
  }
  update_sparsity_pattern(m,n,m_l,n_l,n_oz,n_nz);
  return;
}

// void PetscMMatrixM::update_sparsity_pattern (const Graph &sparsity_pattern){
//     
//     // Clear initialized matrices
//     if (this->initialized())    this->clear();
//     this->_is_initialized = true;
//     int proc_id = 0;   MPI_Comm_rank (this->comm().get(), &proc_id);
//     
//        // Clear initialized matrices
//     const  int m   =  sparsity_pattern._m;//this->_dof_map->n_dofs();
//     const  int n   =  sparsity_pattern._n;
//     const  int n_l =sparsity_pattern._nl;
//     const  int m_l =sparsity_pattern._ml;
//     const  int ml_start =sparsity_pattern._ml_start;
//     std::vector< int> n_nz(m_l);
//     std::vector< int> n_oz(m_l);
//    for(uint i=0;i<m_l;i++) {
//      uint len=sparsity_pattern[ml_start+i].size()-1;
//      n_oz[i]=sparsity_pattern[ml_start+i][len];
//       n_nz[i]=len-n_oz[i];
//    }
//     
// //    const  int n_l = this->_dof_map->n_dofs_on_processor(proc_id);
// //   const  int m_l = n_l;
// //   const std::vector< int>& n_nz = this->_dof_map->get_n_nz();
// //   const std::vector< int>& n_oz = this->_dof_map->get_n_oz();
// //
// //   // Make sure the sparsity pattern isn't empty unless the matrix is 0x0
// //   assert (n_nz.size() == n_l);  assert (n_oz.size() == n_l);
// //
// //   // We allow 0x0 matrices now
// //   //if (m==0)
// //   //  return;
// //
//   int ierr     = 0;
//   int m_global = static_cast<int>(m);
//   int n_global = static_cast<int>(n);
//   int m_local  = static_cast<int>(m_l);
//   int n_local  = static_cast<int>(n_l);
// 
//   int nprocs; MPI_Comm_size(this->comm().get(), &nprocs);
//   // create a sequential matrix on one processor
//   if ( nprocs == 1 )    {
// //       assert ((m_l == m) && (n_l == n));
//       if (n_nz.empty())
//         ierr = MatCreateSeqAIJ (this->comm().get(), m_global, n_global,
// 			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
//       else
//         ierr = MatCreateSeqAIJ (this->comm().get(), m_global, n_global,
// 			        PETSC_NULL, (int*) &n_nz[0], &_mat);
//              CHKERRABORT(this->comm().get(),ierr);
// 
//       ierr = MatSetFromOptions (_mat);
//              CHKERRABORT(this->comm().get(),ierr);
//     }
//   else    {
//       parallel_onlyM();
//       if (n_nz.empty())
//         ierr = MatCreateMPIAIJ (this->comm().get(),
// 			        m_local, n_local,
// 			        m_global, n_global,
// 			        PETSC_NULL, (int*) PETSC_NULL,
// 			        PETSC_NULL, (int*) PETSC_NULL, &_mat);
//       else
//         ierr = MatCreateMPIAIJ (this->comm().get(),
// 			        m_local, n_local,
// 			        m_global, n_global,
// 			        PETSC_NULL, (int*) &n_nz[0],
// 			        PETSC_NULL, (int*) &n_oz[0], &_mat);
//              CHKERRABORT(this->comm().get(),ierr);
//       ierr = MatSetFromOptions (_mat);
//              CHKERRABORT(this->comm().get(),ierr);
//     }
//    this->zero();
// }

// ===================================
void PetscMMatrixM::zero () {
    assert (this->initialized());
    semiparallel_onlyM();
    int ierr=0;
    ierr = MatZeroEntries(_mat);
    CHKERRABORT(this->comm().get(),ierr);
}

// =======================================================
/// This function sets to zero the desired rows:

void PetscMMatrixM::zero_rows (std::vector<int> & rows, double diag_value) 
{
/// MatZeroRows now takes two additional optional arguments.
///  The optional arguments (x,b) can be used to specify the
/// solutions for the zeroed rows (x) and right hand side (b) to update.
/// Could be useful for setting boundary conditions
// =======================================================

  assert(this->initialized());
  semiparallel_onlyM(); int ierr=0;

#if PETSC_VERSION_RELEASE && PETSC_VERSION_LESS_THAN(3,1,1)
  if(!rows.empty()) ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value);
  else  ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value);
#else

  if(!rows.empty())  ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value, PETSC_NULL, PETSC_NULL);
  else  ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value, PETSC_NULL, PETSC_NULL);
#endif
  CHKERRABORT(this->comm().get(),ierr);
//     assert (this->initialized());
//     semiparallel_onlyM();
//     int ierr=0;
//     if (!rows.empty())
//         ierr = MatZeroRows(_mat, rows.size(), &rows[0], diag_value);
//     else
//         ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value);
//     CHKERRABORT(this->comm().get(),ierr);
  return;
}


// =================================================
/// This function clears the matrix
void PetscMMatrixM::clear () 
{// ============================================
    int ierr=0;
    if ((this->initialized()) && (this->_destroy_mat_on_exit)) {
        semiparallel_onlyM();
        ierr = MatDestroy (&_mat);
        CHKERRABORT(this->comm().get(),ierr);
        this->_is_initialized = false;
    }
}

// ============================================
double PetscMMatrixM::l1_norm () const {
    assert (this->initialized());

    semiparallel_onlyM();
    int ierr=0;
    PetscReal petsc_value;
    double value;
    assert (this->closed());
    ierr = MatNorm(_mat, NORM_1, &petsc_value);
    CHKERRABORT(this->comm().get(),ierr);
    value = static_cast<double>(petsc_value);
    return value;
}

// ==========================================
double PetscMMatrixM::linfty_norm () const {
    assert (this->initialized());
    semiparallel_onlyM();
    int ierr=0;
    PetscReal petsc_value;
    double value;
    assert (this->closed());
    ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
    CHKERRABORT(this->comm().get(),ierr);
    value = static_cast<double>(petsc_value);
    return value;
}



// =====================================================
void PetscMMatrixM::print_matlab (const std::string name) const {
    assert (this->initialized());
    semiparallel_onlyM();
    // assert (this->closed());
    this->close();  int ierr=0;
    PetscViewer petsc_viewer;
    ierr = PetscViewerCreate (this->comm().get(), &petsc_viewer);
    CHKERRABORT(this->comm().get(),ierr);

    /**
     * Create an ASCII file containing the matrix
     * if a filename was provided.
     */
    if (name != "NULL")  {
        ierr = PetscViewerASCIIOpen( this->comm().get(),name.c_str(),&petsc_viewer);
        CHKERRABORT(this->comm().get(),ierr);
        ierr = PetscViewerSetFormat (petsc_viewer,PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT(this->comm().get(),ierr);
        ierr = MatView (_mat, petsc_viewer);
        CHKERRABORT(this->comm().get(),ierr);
    }

    // Otherwise the matrix will be dumped to the screen.
    else {
        ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT(this->comm().get(),ierr);
        ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRABORT(this->comm().get(),ierr);
    }
    // Destroy the viewer.
    ierr = PetscViewerDestroy (&petsc_viewer); CHKERRABORT(this->comm().get(),ierr);
}

// ============================================================
/// This Petsc adds a dense matrix to a sparse matrix
void PetscMMatrixM::add_matrix(
  const DenseMatrixM& dm,        // dense matrix
  const std::vector< int>& rows, // row vector
  const std::vector< int>& cols) // column vector 
{// ==================================================
    assert (this->initialized());  int ierr=0;
    const  int m = dm.m();  const  int n = dm.n();
    assert ((int)rows.size() == m);  assert ((int)cols.size() == n);
    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues(_mat,m, (int*) &rows[0],n, (int*) &cols[0],
                        (PetscScalar*) &dm.get_values()[0],ADD_VALUES);
    CHKERRABORT(this->comm().get(),ierr);
    return;
}

// ===========================================================
/// This Petsc function exctract submatrix from matrix
void PetscMMatrixM::_get_submatrix(
  SparseMMatrixM& submatrix,       // submatrix
  const std::vector< int> &rows,   // row vector
  const std::vector< int> &cols,   // column vector
  const bool reuse_submatrix       // flag submatrix
) const 
{ // =====================================================
    this->close();  // Can only extract submatrices from closed matrices

    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMMatrixM* petsc_submatrix = dynamic_cast<PetscMMatrixM*>(&submatrix);

    // If we're not reusing submatrix and submatrix is already initialized
    // then we need to clear it, otherwise we get a memory leak.
    if ( !reuse_submatrix && submatrix.initialized() )  submatrix.clear();
    // Construct row and column index sets.
    int ierr=0;    IS isrow, iscol;
    ierr = ISCreateGeneral(this->comm().get(),rows.size(),
                           (int*) &rows[0],PETSC_USE_POINTER,&isrow);
    CHKERRABORT(this->comm().get(),ierr);

    ierr = ISCreateGeneral(this->comm().get(),cols.size(),
                           (int*) &cols[0],PETSC_USE_POINTER,&iscol);
    CHKERRABORT(this->comm().get(),ierr);

    // Extract submatrix
#if !PETSC_VERSION_LESS_THAN(3,0,1) || !PETSC_VERSION_RELEASE
    ierr = MatGetSubMatrix(_mat, isrow, iscol,
//                            PETSC_DECIDE,   //PETSC version 3.1
                           (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                           &(petsc_submatrix->_mat));
    CHKERRABORT(this->comm().get(),ierr);
#else
   ierr = MatGetSubMatrix(_mat,isrow, iscol,
                            PETSC_DECIDE,
                           (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                           &(petsc_submatrix->_mat));
    CHKERRABORT(this->comm().get(),ierr);

#endif
    // Specify that the new submatrix is initialized and close it.
    petsc_submatrix->_is_initialized = true;
    petsc_submatrix->close();

    // Clean up PETSc data structures
    ierr = ISDestroy(&isrow);  CHKERRABORT(this->comm().get(),ierr);
    ierr = ISDestroy(&iscol);  CHKERRABORT(this->comm().get(),ierr);
    return;
}

// ======================================================
void PetscMMatrixM::get_diagonal (NumericVectorM& dest) const {
    // Make sure the NumericVector passed in is really a PetscVector
    PetscVectorM& petsc_dest = dynamic_cast<PetscVectorM&>(dest);
    // Call PETSc function.
#if PETSC_VERSION_LESS_THAN(2,3,1)

    std::cout << "This method has been developed with PETSc 2.3.1.  "
              << "No one has made it backwards compatible with older "
              << "versions of PETSc so far; however, it might work "
              << "without any change with some older version." << std::endl;
    libmesh_error();
#else
    // Needs a const_cast since PETSc does not work with const.
    int ierr =  MatGetDiagonal(const_cast<PetscMMatrixM*>(this)->mat(),petsc_dest.vec());
    CHKERRABORT(this->comm().get(),ierr);
#endif
}


// ===========================================================
void PetscMMatrixM::get_transpose (SparseMMatrixM& dest) const {
    // Make sure the SparseMatrix passed in is really a PetscMatrix
    PetscMMatrixM& petsc_dest = dynamic_cast<PetscMMatrixM&>(dest);

    // If we aren't reusing the matrix then need to clear dest,
    // otherwise we get a memory leak
    if (&petsc_dest != this)    dest.clear();

    int ierr;
#if PETSC_VERSION_LESS_THAN(3,0,0)
    if (&petsc_dest == this)
        ierr = MatTranspose(_mat,PETSC_NULL);
    else
        ierr = MatTranspose(_mat,&petsc_dest._mat);
    CHKERRABORT(this->comm().get(),ierr);
#else
    // FIXME - we can probably use MAT_REUSE_MATRIX in more situations
    if (&petsc_dest == this)
        ierr = MatTranspose(_mat,MAT_REUSE_MATRIX,&petsc_dest._mat);
    else
        ierr = MatTranspose(_mat,MAT_INITIAL_MATRIX,&petsc_dest._mat);
    CHKERRABORT(this->comm().get(),ierr);
#endif
    // Specify that the transposed matrix is initialized and close it.
    petsc_dest._is_initialized = true;
    petsc_dest.close();
}
// =====================================================
void PetscMMatrixM::print_hdf5(const std::string name) const {
  
  
  
  
//   hid_t status =0;
//  
// //   std::ostringstream name;
// //   name << femus_dir << "/" << input_dir << f_matrix  << Level1 << ext_h5;
//   // print quadratic  --------------------------------------------------
// //    hid_t fileM= H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//    hid_t fileM = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
// //   hid_t fileM = H5Fopen(name.c_str(),  H5P_DEFAULT,H5P_DEFAULT);
//   // matrix dimensions
//   hsize_t dimsf[2];     dimsf[0]=2;  dimsf[1] = 1;
//    int rowcln[2];   rowcln[0]=_n_l;    rowcln[1]=_m_l;
//    
//    
//   hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
//   hid_t dataset = H5Dcreate(fileM,"DIM",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,rowcln);
//   H5Sclose(dataspace);
//   H5Dclose(dataset); 
//    
//    
// //    status =_mgutils.print_Ihdf5(fileM,"DIM", dimsf,rowcln);
//   // print matrix position of P
// //    std::ostringstream name1;name1 << "POS" << qq;
// //   dimsf[0]=count_q;  status =_mgutils.print_Ihdf5(fileM,name1.str().c_str(), dimsf,Mat_q);
//   // print row length of R
// //   std::ostringstream name2;name2 << "LEN";
//   dimsf[0]=_m_l;//    status =_mgutils.print_Ihdf5(fileM,name2.str().c_str(), dimsf,_len_row[0]);
//    dataspace = H5Screate_simple(2,dimsf, NULL);
//    dataset = H5Dcreate(fileM,"LEN",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&_len_row[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   
//   // print row length of R
// //   std::ostringstream name3;name3 << "OFFLEN";
// //   status =_mgutils.print_Ihdf5(fileM,name3.str().c_str(), dimsf,_len_diag_row[0]);
//     dataspace = H5Screate_simple(2,dimsf, NULL);
//    dataset = H5Dcreate(fileM,"DIAGLEN",H5T_NATIVE_INT,
//                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&_len_diag_row[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   //   clean
//    H5Fclose(fileM);
//   
  
//   assert(this->initialized());  semiparallel_onlyM();
//   // assert (this->closed());
//   this->close();  int ierr=0;
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerCreate(this->comm().get(), &petsc_viewer);
//   CHKERRABORT(this->comm().get(),ierr);
// 
//   // Otherwise the matrix will be dumped to the screen.
//   ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
//                               PETSC_VIEWER_ASCII_MATLAB);
//   CHKERRABORT(this->comm().get(),ierr);
// 
//   ierr = MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
//   CHKERRABORT(this->comm().get(),ierr);
// 
//   // Destroy the viewer.
//   ierr = PetscViewerDestroy(petsc_viewer);
//   CHKERRABORT(this->comm().get(),ierr);
// #ifdef PRINT_INFO
std::cout << " PetscMatrixM::print_hdf5 in file   " <<  name   <<   "\n";
// #endif-
}


// // =====================================================
// void PetscMMatrixM::read_hdf5(const std::string namefile,const int mode,const int ml_init) {
//  
//   
//    hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
// 
//   // quadratic-quadratic  ------------------------------------------------
//   std::ostringstream name_dst;  name_dst << "DIM" <<mode;  
//   int rowcln_q[2];hid_t dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   hid_t status=H5Dread(dataset,
//                        H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln_q);
//   int n_row_q=rowcln_q[0];
//   // matrix offset (subdomainsand levelels)
// //   int *off_nd_q= _mgmesh._off_nd;
//   // row legth
//   name_dst.str(""); name_dst << "LEN" <<mode; 
//   uint *length_row=new uint[n_row_q+1];
//   dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_row);
//   // matrix off diagonal
//   name_dst.str(""); name_dst <<"OFFLEN" <<mode; 
//   uint *length_offrow=new uint[n_row_q+1];
//   dataset=H5Dopen(file_id,name_dst.str().c_str(), H5P_DEFAULT);
//   status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_offrow);
// //   // matrix off diagonal
// //   uint *pos_row=new uint[length_row[n_row_q]];dataset=H5Dopen(file_id,"POS2", H5P_DEFAULT);
// //   status=H5Dread(dataset,
// //                  H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row);
// 
//     for ( int i=ml_init; i<ml_init+_m_l; i++) {
//       // row index
//       _len_row[i-ml_init]= length_row[i+1]-length_row[i];
//       int  loffrow= length_offrow[i+1]-length_offrow[i];
//       _len_diag_row[i-ml_init]=_len_row[i-ml_init]-loffrow;
//     }
//     
//   delete   []length_row;  
//   delete   []length_offrow;
//   H5Fclose(file_id);
//   return;
// 
// }

#endif // #ifdef LIBMESH_HAVE_PETSC
