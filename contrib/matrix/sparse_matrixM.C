// ========================================
//            SparseMatrix class
// ========================================
// std lib includes -------------------------
#include <fstream>  // std::ofsream
#include <sstream>  // std::ostringstream
#include "hdf5.h"   // hdf5

// this class ----------------------------
#include "sparse_matrixM.h" // definition header 

// configuration files -------------------
#include "Solverlib_conf.h"  // in ./config

// local includes ------------------------
#include "numeric_vectorM.h" // numerical base vector class
#include "laspack_matrixM.h" // laspack matrix (son of this class)
#include "petsc_matrixM.h"   // petsc matrix   (son of this class)
// #include "trilinos_epetra_matrix.h"  // trilinos matrix

// ====================================================================
//            SparseMatrix Methods: Constructor/Destructor/init/build
// ===================================================================

// =====================================================================================
  /// Constructor;  before usage call init(...).
  SparseMatrixM::SparseMatrixM(
    const ParallelM::Communicator &comm   // parallel communicator <-
  ): ParallelObjectM(comm), _is_initialized(false){}

// =====================================================================================
/// This function builds a  SparseMatrixM using the linear solver 
/// package specified by  solver_package
std::unique_ptr<SparseMatrixM > SparseMatrixM::build( // -----
const ParallelM::Communicator &comm,
  const SolverPackageM solver_package //  solver_package
) { // =================================================================================
  // Build the appropriate vector
  switch (solver_package) {
#ifdef HAVE_LASPACKM // ----------------------------
  case LASPACK_SOLVERSM: {std::unique_ptr<SparseMatrixM> ap(new LaspackMatrixM(comm));return ap;}
#endif
#ifdef HAVE_PETSCM // ------------------------------
  case PETSC_SOLVERSM: {std::unique_ptr<SparseMatrixM > ap(new PetscMatrixM(comm));return ap;}
#endif
#ifdef HAVE_TRILINOSM // ----------------------------
  case TRILINOS_SOLVERSM: {std::unique_ptr<SparseMatrixM > ap(new EpetraMatrix<double>);return ap;}
#endif
  default: std::cerr << "SolverPackageM solver_package:  Unrecognized: " << solver_package;
    abort();
  }
  std::cerr << "SolverPackageM solver_package:  Unrecognized: " << solver_package;
//   std::unique_ptr<SparseMatrixM > ap(NULL);
//   return ap;
}

// =================================================
//            SparseMatrix Methods: Add/mult
// =================================================

// =========================================================
void SparseMatrixM::vector_mult(NumericVectorM& dest,
                                const NumericVectorM& arg) const {
  dest.zero();  this->vector_mult_add(dest,arg);
}

// ===========================================================
void SparseMatrixM::vector_mult_add(NumericVectorM& dest,
                                    const NumericVectorM& arg) const {
  /* This functionality is actually implemented in the \p NumericVector class.  */
  dest.add_vector(arg,*this);
}

// =================================================
//            SparseMatrix Methods: setting/return
// =================================================

// ==========================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
void SparseMatrixM::zero_rows(std::vector<int> &, double) {
  std::cerr << " SparseMatrixM::zero_rows: Not implememnted "; abort();
}


// ===============================================
//            SparseMatrix Methods: Print/Read
// ==============================================
// // ----------------------------------------------
// void SparseMatrixM::print(const std::string& name) const{
// std::cerr << "SparseMatrixM::print:not implemented";
// //    std::ofstream outfile(name.c_str());
// // //    outfile<<" Matrix Level "<<std::endl <<0<<"  "<<m()<< std::endl;
// //    // data
// //    for (uint i=1; i<= m(); i++) {
// //      std::vector<uint> row;
// //      for(uint j=0; j < n(); j++)
// //        if((*this)(i-1,j) > 0.)
// // 	 row.push_back(j);
// //      outfile << row.size() << " ";
// //      for (uint j=0; j< row.size(); j++) outfile  << row[j] << "  ";
// //      outfile << std::endl;
// //    }
// //    outfile.close();
// //   //  std::cout  << Level << " level matrix written in  " <<name.c_str()<< std::endl;  std::cout  << " dim  " << Q_GetDim(&matrix)<< std::endl;
// //   return;
// }


// void SparseMatrixM::print(std::ostream& os) const
//{
//   parallel_only();
//
//   libmesh_assert (this->initialized());
//
//   // We'll print the matrix from processor 0 to make sure
//   // it's serialized properly
//   if (libMesh::processor_id() == 0){
// //       libmesh_assert(this->_dof_map->first_dof() == 0);
//       for (unsigned int i=this->_dof_map->first_dof();
//            i!=this->_dof_map->end_dof(); ++i){
//           for (unsigned int j=0; j<this->n(); j++)   os << (*this)(i,j) << " ";
//           os << std::endl;
//         }
//
//       std::vector<unsigned int> ibuf, jbuf;
//       std::vector<double> cbuf;
//       unsigned int currenti = this->_dof_map->end_dof();
//       for (unsigned int p=1; p < libMesh::n_processors(); ++p)  {
//           Parallel::receive(p, ibuf);
//           Parallel::receive(p, jbuf);
//           Parallel::receive(p, cbuf);
//           libmesh_assert(ibuf.size() == jbuf.size());
//           libmesh_assert(ibuf.size() == cbuf.size());
//
//           if (ibuf.empty())
//             continue;
//           libmesh_assert(ibuf.front() >= currenti);
//           libmesh_assert(ibuf.back() >= ibuf.front());
//
//           unsigned int currentb = 0;
//           for (;currenti <= ibuf.back(); ++currenti)            {
//               for (unsigned int j=0; j<this->n(); j++)                {
//                   if (currentb < ibuf.size() && ibuf[currentb] == currenti && jbuf[currentb] == j) {
// 	              os << cbuf[currentb] << " ";
// 	              currentb++;
//                     }
//                   else	os << static_cast<double>(0.0) << " ";
//                 }
//               os << std::endl;
//             }
//         }
//       for (; currenti != this->m(); ++currenti)   {
//           for (unsigned int j=0; j<this->n(); j++)  os << static_cast<double>(0.0) << " ";
//           os << std::endl;
//         }
//     }
//   else
//     {
//       std::vector<unsigned int> ibuf, jbuf;
//       std::vector<double> cbuf;
//
//       // We'll assume each processor has access to entire
//       // matrix rows, so (*this)(i,j) is valid if i is a local index.
//       for (unsigned int i=this->_dof_map->first_dof();
//            i!=this->_dof_map->end_dof(); ++i)        {
//           for (unsigned int j=0; j<this->n(); j++)	    {
//               double c = (*this)(i,j);
//               if (c != static_cast<double>(0.0))                {
//                   ibuf.push_back(i); jbuf.push_back(j); cbuf.push_back(c);
//                 }
// 	    }
//         }
//       Parallel::send(0,ibuf); Parallel::send(0,jbuf); Parallel::send(0,cbuf);
//     }
//}
// =====================================================
 void SparseMatrixM::print(const std::string& name)const{
     std::ofstream out(name.c_str());   print(out);
 }
 
 
 #if HDF5_VERSION == 188
// =====================================================
/// This function reads len sparse matrix structures
void SparseMatrixM::read_len_hdf5(const std::string namefile,  // file name    
				  const std::string mode,              // linear or quadratic
                                  int len_row[],               // row lengths
				  int len_off_row[]            // row off entry lengths
                                  ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
  
  // row lengths -------------------------------------------------------------
  std::ostringstream name_dst;
  if (mode=="0" || mode=="1" || mode=="2") name_dst <<"LEN0"<< mode;
  else name_dst <<"LEN"<< mode;
  
  hid_t dataset=H5Dopen(file_id,name_dst.str().c_str());
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len_row);
  assert(status==0);
  // matrix row off entry lengths --------------------------------------------
  name_dst.str("");
  if (mode=="0" || mode=="1" || mode=="2") name_dst <<"OFFLEN0"<< mode;
  else name_dst <<"OFFLEN"<< mode;

  dataset=H5Dopen(file_id,name_dst.str().c_str());
  status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len_off_row);
  assert(status==0);
  H5Fclose(file_id);
  return;

}

// =====================================================
/// This function reads pos sparse matrix structure
void SparseMatrixM::read_pos_hdf5(const std::string namefile,  // file name 
				  const int mode,              // linear or quadratic
                                  int pos_row[]                // compressed row positions
				 ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  std::ostringstream name_dst;  name_dst.str(""); name_dst << "POS" <<mode;
  hid_t dataset=H5Dopen(file_id,name_dst.str().c_str());
  H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row);
  H5Fclose(file_id);
  return;
}
// =====================================================
/// This function reads quad-linear matrix dimensions
int SparseMatrixM::read_dim_hdf5(const std::string namefile, // file name 				
				  const std::string mode,    // linear or quadratic
				  int ldim[]                  // dimensions
				 ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  
  std::ostringstream name_label;
  if (mode=="0" || mode=="1" || mode=="2") name_label <<"DIM0"<< mode;
  else name_label <<"DIM"<< mode;
  
  hid_t dataset=H5Dopen(file_id,name_label.str().c_str());
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ldim);
  assert(status==0);
  H5Fclose(file_id);
  return( (int) status);
}
#else

// =====================================================
/// This function reads len sparse matrix structures
void SparseMatrixM::read_len_hdf5(const std::string namefile,  // file name    
				  const std::string mode,              // linear or quadratic
                                  int len_row[],               // row lengths
				  int len_off_row[]            // row off entry lengths
                                  ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
  
  // row lengths -------------------------------------------------------------
  std::ostringstream name_dst;
  if (mode=="0" || mode=="1" || mode=="2") name_dst <<"LEN0"<< mode;
  else name_dst <<"LEN"<< mode;
  
  hid_t dataset=H5Dopen(file_id,name_dst.str().c_str()
     #if HDF5_VERSIONM!=1808                        
                        , H5P_DEFAULT                    
    #endif
  );
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len_row);
  assert(status==0);
  // matrix row off entry lengths --------------------------------------------
  name_dst.str("");
  if (mode=="0" || mode=="1" || mode=="2") name_dst <<"OFFLEN0"<< mode;
  else name_dst <<"OFFLEN"<< mode;

  dataset=H5Dopen(file_id,name_dst.str().c_str()
     #if HDF5_VERSIONM!=1808                        
                        , H5P_DEFAULT                    
    #endif
  );
  status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len_off_row);
  assert(status==0);
  H5Fclose(file_id);
  return;

}

// =====================================================
/// This function reads pos sparse matrix structure
void SparseMatrixM::read_pos_hdf5(const std::string namefile,  // file name 
				  const int mode,              // linear or quadratic
                                  int pos_row[]                // compressed row positions
				 ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  std::ostringstream name_dst;  name_dst.str(""); name_dst << "POS" <<mode;
  hid_t dataset=H5Dopen(file_id,name_dst.str().c_str()
     #if HDF5_VERSIONM!=1808                        
                        , H5P_DEFAULT                    
    #endif
  );
  H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row);
  H5Fclose(file_id);
  return;
}
// =====================================================
/// This function reads quad-linear matrix dimensions
int SparseMatrixM::read_dim_hdf5(const std::string namefile, // file name 				
				  const std::string mode,    // linear or quadratic
				  int ldim[]                  // dimensions
				 ) {

  hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  
  std::ostringstream name_label;
  if (mode=="0" || mode=="1" || mode=="2") name_label <<"DIM0"<< mode;
  else name_label <<"DIM"<< mode;
  
  hid_t dataset=H5Dopen(file_id,name_label.str().c_str()
     #if HDF5_VERSIONM!=1808                        
                        , H5P_DEFAULT                    
    #endif
  );
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,ldim);
  assert(status==0);
  H5Fclose(file_id);
  return( (int) status);
}
#endif



//------------------------------------------------------------------
