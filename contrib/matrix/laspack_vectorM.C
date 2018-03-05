#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM

#include "Typedefs_conf.h"

#include "hdf5.h"
// --------------------------------
// C++ includes
#include <algorithm> // for std::min
#include <limits>
#include <sstream>  //for print hdf5

// --------------------------------
// Local Includes
#include "dense_subvectorM.h"
#include "dense_vectorM.h"
#include "laspack_vectorM.h"
#include "laspack_matrixM.h"
#include "laspack_MmatrixM.h"
// ---------------------------------



// ================================
Real LaspackVectorM::sum() const {
  assert(this->closed());
  Real sum = 0;
  const unsigned int n = this->size();
  for (unsigned int i=1; i<=n; ++i) sum += V__GetCmp(&_vec,i);
  return sum;
}
// ================================
Real LaspackVectorM::l1_norm() const {
  assert(this->closed());
  return static_cast<Real>(l1Norm_V(const_cast<QVector*>(&_vec)));
}
// ================================
Real LaspackVectorM::l2_norm() const {
  assert(this->closed());
  return static_cast<Real>(l2Norm_V(const_cast<QVector*>(&_vec)));
}
// ================================
Real LaspackVectorM::linfty_norm() const {
  assert(this->closed());
  return static_cast<Real>(MaxNorm_V(const_cast<QVector*>(&_vec)));
}
// ================================
// ================================
NumericVectorM& LaspackVectorM::operator += (const NumericVectorM& v) {
  assert(this->closed());   this->add(1., v);  return *this;
}
// ================================
NumericVectorM& LaspackVectorM::operator -= (const NumericVectorM& v) {
  assert(this->closed());  this->add(-1., v);  return *this;
}
// ================================
void LaspackVectorM::add(const Real v) {
  const unsigned int n = this->size();
  for (unsigned int i=0; i<n; i++)  this->add(i, v);
}
// ===================================
void LaspackVectorM::add(const NumericVectorM& v) {  this->add(1., v);}
// ===================================
void LaspackVectorM::add(const Real a, const NumericVectorM& v_in) {
  // Make sure the vector passed in is really a LaspackVectorM
  const LaspackVectorM* v = static_cast<const LaspackVectorM*>(&v_in);
  assert(v != NULL);  assert(this->size() == v->size());
  for (unsigned int i=0; i<v->size(); i++)    this->add(i,a*(*v)(i));
}
// ===================================
// ===================================
void LaspackVectorM::add_vector(const std::vector<Real>& v,
                                const std::vector<unsigned int>& dof_indices) {
  assert(!v.empty());  assert(v.size() == dof_indices.size());
  for (unsigned int i=0; i<v.size(); i++)    this->add(dof_indices[i], v[i]);
}
// ===================================
void LaspackVectorM::add_vector(const NumericVectorM& V,
                                const std::vector<unsigned int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)    this->add(dof_indices[i], V(i));
}
// ================================
void LaspackVectorM::add_vector(const DenseVectorM& V,
                                const std::vector<unsigned int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)    this->add(dof_indices[i], V(i));
}
// ==================================
void LaspackVectorM::add_vector(const NumericVectorM &vec_in,
                                const SparseMatrixM &mat_in) {
  // Make sure the data passed in are really in Laspack types
  const LaspackVectorM* vec = static_cast<const LaspackVectorM*>(&vec_in);
  const LaspackMatrixM* mat = static_cast<const LaspackMatrixM*>(&mat_in);
  assert(vec != NULL);  assert(mat != NULL);
  // += mat*vec
  AddAsgn_VV(&_vec, Mul_QV(const_cast<QMatrix*>(&mat->_QMat),
                           const_cast<QVector*>(&vec->_vec)));
  return;
}

// ==================================
void LaspackVectorM::add_vector(const NumericVectorM &vec_in,
                                const SparseMMatrixM &mat_in) {
  // Make sure the data passed in are really in Laspack types
  const LaspackVectorM* Lvec = static_cast<const LaspackVectorM*>(&vec_in);
  const LaspackMMatrixM* Lmat = static_cast<const LaspackMMatrixM*>(&mat_in);
  assert(Lvec != NULL);  assert(Lmat != NULL);
  // += mat*vec
//   const Matrix *mat=const_cast<Matrix*>(&Lmat->_Mat);
//   const QVector * vec=const_cast<QVector*>(&Lvec->_vec);
//    assert(mat != NULL);  assert(vec != NULL);

  AddAsgn_VV(&_vec, Mul_MV(const_cast<Matrix*>(&Lmat->_Mat),const_cast<QVector*>(&Lvec->_vec)));
//   uint n_row=this->size();
//   for (uint irow=1;irow<=n_row;irow++) {
//     double sum=0.;
//     for (uint j=0;j<M__GetLen(mat,irow);j++) sum +=M__GetVal(mat,irow,j)*V__GetCmp(vec,M__GetPos(mat,irow,j));
//   V__AddCmp(&this->_vec,irow,sum);
// }
  return;
}

// ==================================
// void LaspackVectorM::add_vector(const NumericVectorM &vec_in,
//                                 const SparseMMatrixM &mat_in) {
//   // Make sure the data passed in are really in Laspack types
//   const LaspackVectorM* Lvec = static_cast<const LaspackVectorM*>(&vec_in);
//   const LaspackMMatrixM* Lmat = static_cast<const LaspackMMatrixM*>(&mat_in);
//   assert(Lvec != NULL);  assert(Lmat != NULL);
//   // += mat*vec
//   const Matrix *mat=const_cast<Matrix*>(&Lmat->_Mat);
//   const QVector * vec=const_cast<QVector*>(&Lvec->_vec);
//    assert(mat != NULL);  assert(vec != NULL);
//    
//   //   AddAsgn_VV (&_vec, Mul_MV(const_cast<Matrix*>(&mat->_Mat),const_cast<QVector*>(&vec->_vec)));
//   uint n_row=this->size();
//   for (uint irow=1;irow<=n_row;irow++) {
//     double sum=0.;
//     for (uint j=0;j<M__GetLen(mat,irow);j++) sum +=M__GetVal(mat,irow,j)*V__GetCmp(vec,M__GetPos(mat,irow,j));
//   V__AddCmp(&this->_vec,irow,sum);
// }
// }

// =====================================================
void LaspackVectorM::resid (const NumericVectorM &rhs_in,const NumericVectorM &vec_in,
                                const SparseMatrixM &mat_in) {
  // Make sure the data passed in are really in Laspack types
  const LaspackVectorM* Lrhs = static_cast<const LaspackVectorM*>(&rhs_in);
  const LaspackVectorM* Lvec = static_cast<const LaspackVectorM*>(&vec_in);
  const LaspackMatrixM* Lmat = static_cast<const LaspackMatrixM*>(&mat_in);
  assert(Lvec != NULL);  assert(Lmat != NULL);
  // += mat*vec
  const QMatrix *mat=const_cast<QMatrix*>(&Lmat->_QMat);
   const QVector * vec=const_cast<QVector*>(&Lvec->_vec);
   const QVector * rhs=const_cast<QVector*>(&Lrhs->_vec);
   assert(mat != NULL);  assert(vec != NULL);
   
  //   AddAsgn_VV (&_vec, Mul_MV(const_cast<Matrix*>(&mat->_Mat),const_cast<QVector*>(&vec->_vec)));
  uint n_row=this->size();
  for (uint irow=1;irow<=n_row;irow++) {
    double sum=0.;
    for (uint j=0;j<Q__GetLen(mat,irow);j++) sum +=Q__GetVal(mat,irow,j)*V__GetCmp(vec,Q__GetPos(mat,irow,j));
  V__SetCmp(&this->_vec,irow,V__GetCmp(rhs,irow)-sum);
}
return;
}
// ============================================================
void LaspackVectorM::matrix_mult(const NumericVectorM &vec_in,
                                 const SparseMMatrixM &mat_in) {
  // Make sure the data passed in are really in Laspack types
  const LaspackVectorM* Lvec = static_cast<const LaspackVectorM*>(&vec_in);
  const LaspackMMatrixM* Lmat = static_cast<const LaspackMMatrixM*>(&mat_in);
  assert(Lvec != NULL);  assert(Lmat != NULL);
  // += mat*vec
  const Matrix *mat=const_cast<Matrix*>(&Lmat->_Mat);
  const QVector * vec=const_cast<QVector*>(&Lvec->_vec);
  assert(mat != NULL);  assert(vec != NULL);

  //   AddAsgn_VV (&_vec, Mul_MV(const_cast<Matrix*>(&mat->_Mat),const_cast<QVector*>(&vec->_vec)));
  uint n_row=this->size();
  for (uint irow=1;irow<=n_row;irow++) {
    double sum=0.;
    for (uint j=0;j<M__GetLen(mat,irow);j++) sum +=M__GetVal(mat,irow,j)*V__GetCmp(vec,M__GetPos(mat,irow,j));
    V__SetCmp(&this->_vec,irow,sum);
  }
  return;
}
// =====================================================
// =====================================================
void LaspackVectorM::insert(const std::vector<Real>& v,
                            const std::vector<unsigned int>& dof_indices) {
  assert(!v.empty());  assert(v.size() == dof_indices.size());
  for (unsigned int i=0; i<v.size(); i++)  V__SetCmp(&this->_vec,dof_indices[i]+1, v[i]);
  //this->set(dof_indices[i], v[i]);
}
// ====================================================
void LaspackVectorM::insert(const NumericVectorM& V,
                            const std::vector<unsigned int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)  V__SetCmp(&this->_vec,dof_indices[i]+1, V(i));
  // this->set(dof_indices[i], V(i));
}
// ===================================================
void LaspackVectorM::insert(const DenseVectorM& V,
                            const std::vector<unsigned int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)   V__SetCmp(&_vec,dof_indices[i]+1, V(i));
  // this->set(dof_indices[i], V(i));
}

// ===================================================
void LaspackVectorM::insert(const DenseSubVectorM& V,
                            const std::vector<unsigned int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for (unsigned int i=0; i<V.size(); i++)   V__SetCmp(&_vec,dof_indices[i]+1, V(i));
  //this->set(dof_indices[i], V(i));
}



// ==================================
void LaspackVectorM::scale(const Real factor) {
  assert(this->initialized());   //   Mul_SV (factor, &_vec);
  for (uint i=1;i<=this->size();i++) V__SetCmp(&_vec,i,factor*V__GetCmp(&_vec,i));

}
// ==================================
void LaspackVectorM::pointwise_mult(const NumericVectorM& /*vec1*/,
                                    const NumericVectorM& /*vec2*/) {
  std::cout << "  libmesh_not_implemented ";abort();
}

// ==================================
void LaspackVectorM::abs() {
  assert(this->initialized());  const unsigned int n = this->size();
  for (unsigned int i=1; i<=n; ++i)  V__SetCmp(&_vec,i,fabs(V__GetCmp(&_vec,i)));
  //this->set(i,std::abs((*this)(i)));
}
// ==================================
Real LaspackVectorM::dot(const NumericVectorM& V) const {
  assert(this->initialized());
  // Make sure the NumericVector passed in is really a LaspackVectorM
  const LaspackVectorM* v = static_cast<const LaspackVectorM*>(&V);
  assert(v != NULL);
  Real sum=0.;
  for (unsigned int i=1; i<=V.size(); i++)  sum +=V__GetCmp(&_vec,i)*
        V__GetCmp(const_cast<QVector*>(&(v->_vec)),i);
  return sum;
//   return Mul_VV(const_cast<QVector*>(&(this->_vec)),const_cast<QVector*>(&(v->_vec)));
}

// ==================================
NumericVectorM& LaspackVectorM::operator = (const Real s) {
  assert(this->initialized());  V_SetAllCmp(&_vec, s);  return *this;
}
// ======================================
NumericVectorM& LaspackVectorM::operator = (const NumericVectorM& v_in) {
  // Make sure the NumericVector passed in is really a LaspackVectorM
  const LaspackVectorM* v =static_cast<const LaspackVectorM*>(&v_in);
  assert(v != NULL);  *this = *v;  return *this;
}

// =======================================
LaspackVectorM& LaspackVectorM::operator = (const LaspackVectorM& v) {
  assert(this->initialized());  assert(this->size() == v.size());
  this->_is_closed = v._is_closed;
  if (v.size() != 0) for (unsigned int i=1; i<=V_GetDim(&_vec); ++i)
      V__SetCmp(&_vec,i,V__GetCmp(const_cast<QVector*>(&v._vec),i));
  //   if (v.size() != 0) Asgn_VV(const_cast<QVector*>(&_vec),const_cast<QVector*>(&v._vec));
  return *this;
}
// ==================================
NumericVectorM& LaspackVectorM::operator = (const std::vector<Real>& v) {
  // The vector is the same size of the global vector.  Only add the local components.
  if (this->size() == v.size()) for (uint i=0;i<v.size();i++) this->set(i, v[i]);
  else abort();
  return *this;
}

// ===============================================
void LaspackVectorM::localize(NumericVectorM& v_local_in) const {
  // Make sure the NumericVector passed in is really a LaspackVectorM
  LaspackVectorM* v_local =    static_cast<LaspackVectorM*>(&v_local_in);
  assert(v_local != NULL);
  *v_local = *this;
}
// =================================================
void LaspackVectorM::localize(NumericVectorM& v_local_in,
                              const std::vector<unsigned int>& send_list) const {
  // Make sure the NumericVector passed in is really a LaspackVectorM
  LaspackVectorM* v_local =    static_cast<LaspackVectorM*>(&v_local_in);
  assert(v_local != NULL);  assert(send_list.size() <= v_local->size());
  *v_local = *this;
}

// ================================================
void LaspackVectorM::localize(const unsigned int first_local_idx,
                              const unsigned int last_local_idx,
                              const std::vector<unsigned int>& send_list) {
  assert(first_local_idx  == 0);  assert(last_local_idx+1 == this->size());
  assert(send_list.size() <= this->size());
}
// ==================================
void LaspackVectorM::localize(std::vector<Real>& v_local) const {
  v_local.resize(this->size());
  for (unsigned int i=0; i<v_local.size(); i++)  v_local[i] = (*this)(i);
}
// ==================================
void LaspackVectorM::localize_to_one(std::vector<Real>& v_local,
                                     const unsigned int pid) const {
  assert(pid == 0);  this->localize(v_local);
}

// ==================================
Real LaspackVectorM::max() const {
  assert(this->initialized());
  if (!this->size())    return -std::numeric_limits<Real>::max();
  Real max = (Real)((*this)(0));
  for (unsigned int i=1; i<this->size(); i++)   max = std::max(max, (Real)((*this)(i)));
  return max;
}
// ==================================
Real LaspackVectorM::min() const {
  assert(this->initialized());
  if (!this->size())    return std::numeric_limits<Real>::max();
  Real min = (Real)((*this)(0));
  const unsigned int n = this->size();
  for (unsigned int i=1; i<n; i++)    min = std::min(min, (Real)((*this)(i)));
  return min;
}
// =====================================================
/// Print the contents of the vector to the screen,
void LaspackVectorM::print_personal(std::ostream& os) const {
  uint ndim=this->size();
  os  << " LaspackVectorM  dim= "  << ndim <<  "   \n ";
  for (uint i=0;i<ndim;i++) os  << "("<< i << "," << V__GetCmp(&_vec,i+1) << ")";
  os << "  \n ";
  return;
}

/// Print the contents of the vector to hdf5 file,
void LaspackVectorM::print_hdf5(std::string file) const {
  uint ndim=this->size();
  std::ostringstream name;
  name.str("");    name << file.c_str() << ".h5";
  hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  // Value storage
  double *val=new double[ndim];  for (uint i=0;i<ndim;i++)  val[i]= V__GetCmp(&_vec,i+1);

// print matrix values
  hsize_t dimsf[2]; dimsf[0]=1;    dimsf[1] = 1;
  int *nd=new int[1];nd[0]=ndim;
  name.str(""); name << "DIM" ;
  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,nd);
  H5Sclose(dataspace);  H5Dclose(dataset);

  dimsf[0]=ndim;
  name.str(""); name << "VAL" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_DOUBLE,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,val);
//    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,&((*this)._vec) /*val*/);
  H5Sclose(dataspace);  H5Dclose(dataset);

 
  H5Fclose(fileP);
  delete []val; delete []nd;
  return;
}

#endif // #ifdef LIBMESH_HAVE_LASPACK
