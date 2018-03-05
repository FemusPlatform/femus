#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM
// =======================
// LaspackMatrix members
// =======================

#include "Typedefs_conf.h"

// Local includes
#include "laspack_matrixM.h"

// C++ includes
#include <iostream>
#include <fstream>

#include "hdf5.h"

// =======================
//   INITIALIZATION
// =======================


// =======================================================
/// This function do  initialize the matrix with dimension mxn
void LaspackMatrixM::init(const unsigned int m,
                          const uint /* n*/, const uint /*m_l*/,
                          const unsigned int /*n_l*/,
                          const unsigned int /*nnz*/,
                          const unsigned int) {
  // Ignore calls on initialized objects
  if (this->initialized()) {std::cout<<" LaspackMatrixM::init: already init"; return;}
  // initialize auxilary row
//   _uncomp_row=new int[n];
  // initialize matrix
  if (m==0)  return;
  Q_Constr(&_QMat, const_cast<char*>("Mat"),m,_LPFalse, Rowws, Normal,_LPTrue);
  this->_is_initialized = true; assert(m == this->m());
}

// ========================================================
/// This function initializes the laspack matrix structure
void LaspackMatrixM::update_sparsity_pattern(const Graph &sparsity_pattern) {
  assert(this->initialized());
  _uncomp_row=new int[this->m()];
  // The matrix structure.
  for (unsigned int i=0; i<this->m(); i++) {// row loop
    _uncomp_row[i]=0; // _uncomp_row=zero.
    const uint length =sparsity_pattern[i].size()-1; // get length
    Q_SetLen(&_QMat, i+1, length);  // set length
    for (unsigned int l=0; l<length; l++) { // Matrix=zero.
      Q_SetEntry(&_QMat, i+1, l, sparsity_pattern[i][l]+1, 0.);
    }
  }
  return;
}


// =========================================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
void LaspackMatrixM::zero_rows(std::vector<int> & rows, Real diag_value) {
  std::cout << "LaspackMatrixM::zero_rows: not implemented \n"; abort();
}

// =============================================
//                     COMPUTING
// =============================================
// ========================================================
/// This function computes  the l1-norm: norm=max over (sum over j-cln).
Real LaspackMatrixM::l1_norm() const {
  Real norm = 0;
  for (uint irow = 0; irow < m(); irow++) {
    Real sum = 0; // set to zero sum coln
    for (uint jcol=0; jcol < Q__GetLen(&_QMat,irow+1); jcol++) sum += fabs(Q__GetVal(&_QMat,irow+1,jcol));
    if (sum > norm) norm = sum; // norm=max (sum over j-cln)
  }
  return norm;
}
// =============================================================
 /// This function computes  the linfty-norm of the matrix
 Real LaspackMatrixM::linfty_norm() const{
    Real norm = 0;
  for (uint irow = 0; irow < m(); irow++) {
    Real sum = 0; // set to zero sum coln
    for (uint jcol=0; jcol < Q__GetLen(&_QMat,irow+1); jcol++) {
      sum = fabs(Q__GetVal(&_QMat,irow+1,jcol));
    if (sum > norm) norm = sum; // norm=max (i-row j-cln)
    }
  }
  return norm;
};

void LaspackMatrixM::get_diagonal(NumericVectorM& /*dest*/) const {
  std::cout << "not_implemented";abort();
}
// ===================================
void LaspackMatrixM::get_transpose(SparseMatrixM& /*dest*/) const {
  std::cout << "not_implemented";abort();
}

// =============================================
//                      PRINT
// =============================================
// ===================================
/// This function prints the contents of the matrix to the screen,
void LaspackMatrixM::print_personal(std::ostream& os) const {
  uint ndim=n();
  os  << " LaspackMatrixM   "  <<   "   \n ";
  for (uint i=0;i<ndim;i++) {
    uint Len=Q__GetLen(&_QMat,i+1);
    os  << i << " row; len="<< Len << "  \n ";
    for (uint j=0;j<Len;j++) {
      os  << "("<<Q__GetPos(&_QMat,i+1,j)-1 << "," << Q__GetVal(&_QMat,i+1,j) << ")";
    }
    os << "  \n ";
  }
  return;
}
// ==============================================================
#include <sstream>
/// This function prints the contents of the matrix to an hdf5 file,
void LaspackMatrixM::print_hdf5(std::string file) const {
  uint ndim=n();
  std::ostringstream name;
  name.str("");    name << file.c_str() << ".h5";
  std::cout << "print file"<< name.str().c_str() << "\n";
  hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  // Matrix dimensions
  uint *offset_pos=new uint[2];
  hsize_t dimsf[2];     dimsf[0]=2;    dimsf[1] = 1;
  offset_pos[0]=m(); offset_pos[1]=n();
  name.str("");    name << "DIM" ;
  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,offset_pos);
  H5Sclose(dataspace);  H5Dclose(dataset);

  // print offset vector -------------------------------------------
  dimsf[0]=ndim+1;    dimsf[1] = 1;
  delete []offset_pos; offset_pos=new uint[ndim+1];  
  offset_pos[0] =0;
  for (uint i=0;i<ndim;i++) { offset_pos[i+1] =offset_pos[i]+Q__GetLen(&_QMat,i+1); }
  name.str("");    name << "OFFSET" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,offset_pos);
  H5Sclose(dataspace);  H5Dclose(dataset);
  // print value vector -------------------------------------------
  int count=offset_pos[ndim];
  delete []offset_pos; offset_pos=new uint[count];  double *val=new double[count];
  int jcount=0;
  for (uint i=0;i<ndim;i++) {
    uint Len=Q__GetLen(&_QMat,i+1);
    for (uint j=0;j<Len;j++) {
      offset_pos[jcount]=Q__GetPos(&_QMat,i+1,j)-1;
      val[jcount]= Q__GetVal(&_QMat,i+1,j);
      jcount++;
    }
  }
  // print matrix pos ----------------------------------------------
  name.str(""); name << "POS" ; dimsf[0]=count;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,offset_pos);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  // print matrix values ------------------------------------------
  name.str(""); name << "VAL" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_DOUBLE,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,val);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Fclose(fileP);

  return;
}

#endif // #ifdef HAVE_LASPACK



// #include "libmeshM.h"
// 
// #ifdef HAVE_LASPACKM
// 
// // C++ includes
// 
// // Local includes
// #include "laspack_matrixM.h" 
// // =======================
// // LaspackMatrix members
// // =======================
// 
// // ========================================================
// void LaspackMatrixM::update_sparsity_pattern (const Graph &sparsity_pattern)
// {
//   assert (this->initialized());
//   // Tell the matrix about its structure.  Initialize it to zero.
//   for (unsigned int i=0; i<this->m(); i++)    {
//       const uint length =sparsity_pattern[i].size();
//       Q_SetLen (&_QMat, i+1, length);
// 
//       for (unsigned int l=0; l<length; l++){
// 	  Q_SetEntry (&_QMat, i+1, l, sparsity_pattern[i][l]+1, 0.);
// 	  // sanity check
// // 	  std::cout << "m()=" << m() << std::endl;
// // 	  std::cout << "(i,j,l) = (" << i << "," << j<< "," << l << ")" << std::endl;
// // 	  std::cout << "pos(i,j)=" << pos(i,j) << std::endl;
// 	}          
//     }
// }        
// 
// // ======================================================= 
// void LaspackMatrixM::init (const unsigned int m,
// 			     const uint  n, const uint /*m_l*/,
// 			     const unsigned int /*n_l*/,
// 			     const unsigned int /*nnz*/,
// 			     const unsigned int)
// {
//   // Ignore calls on initialized objects
//   if (this->initialized()) {
//     std::cout<<" LaspackMatrixM::init: already init"; return;
//   }
//   
//   // initialize auxilary row
//   _uncomp_row=new int[n];
//   // initialize matrix
//   if (m==0)  return;
//   Q_Constr(&_QMat, const_cast<char*>("Mat"),m,_LPFalse, Rowws, Normal,_LPTrue);
//   this->_is_initialized = true; assert (m == this->m());
// }
// // ======================================================== 
// Real LaspackMatrixM::l1_norm() const {
//   Real norm = 0;
//     for (uint irow = 0; irow < m(); irow++) {
//     Real sum = 0;
//     for(uint jcol=0; jcol < Q__GetLen(&_QMat,irow+1); jcol++) sum += fabs(Q__GetVal(&_QMat,irow+1,jcol));
//     if (sum > norm) norm = sum;
//   }
//   return norm;
// }
// 
// // ======================================================== 
// void LaspackMatrixM::init (){
//   std::cout << "  LaspackMatrixM::init: Our LaspackMatrix  needs dimensions"; abort();
// }
// // ===================================
// void LaspackMatrixM::get_diagonal (NumericVectorM& /*dest*/) const{
//   std::cout << "not_implemented";exit(0);
// }
// // ===================================
// void LaspackMatrixM::get_transpose (SparseMatrixM& /*dest*/) const{ 
//   std::cout << "not_implemented";exit(0);
// }
// 
// 
// 
// 
// #endif // #ifdef LIBMESH_HAVE_LASPACK
