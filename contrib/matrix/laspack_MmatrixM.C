#include "Solverlib_conf.h"

#ifdef HAVE_LASPACKM
// =======================
// LaspackMatrix members
// =======================

#include "Typedefs_conf.h"

// Local includes
#include "laspack_MmatrixM.h" 


// ========================================================
void LaspackMMatrixM::update_sparsity_pattern (const Graph &sparsity_pattern)
{
  assert (this->initialized());
  // Tell the matrix about its structure.  Initialize it to zero.
  for (unsigned int i=0; i<this->m(); i++)    {
      const uint length =sparsity_pattern[i].size();
      M_SetLen (&_Mat, i+1, length);

      for (unsigned int l=0; l<length; l++){
	  M_SetEntry (&_Mat, i+1, l, sparsity_pattern[i][l]+1, 0.);
	  // sanity check
// 	  std::cout << "m()=" << m() << std::endl;
// 	  std::cout << "(i,j,l) = (" << i << "," << j<< "," << l << ")" << std::endl;
// 	  std::cout << "pos(i,j)=" << pos(i,j) << std::endl;
	}          
    }
}        

// ======================================================= 
void LaspackMMatrixM::init (const unsigned int m,
			     const uint  n, const uint /*m_l*/,
			     const unsigned int /*n_l*/,
			     const unsigned int /*nnz*/,
			     const unsigned int)
{
  // Ignore calls on initialized objects
  if (this->initialized()) {
    std::cout<<" LaspackMMatrixM::init: already init"; return;
  }
  
  // initialize auxilary row
  _uncomp_row=new int[n];
  // initialize matrix
  if (m==0)  return;
  M_Constr(&_Mat, const_cast<char*>("Mat"),m, n, Rowws, Normal,_LPTrue);

  this->_is_initialized = true; assert (m == this->m());
}
// ======================================================== 
/// This function computes  the l1-norm: norm=max over (sum over j-cln).
Real LaspackMMatrixM::l1_norm() const {
  Real norm = 0;
    for (uint irow = 0; irow < m(); irow++) {
    Real sum = 0;
    for(uint jcol=0; jcol < M__GetLen(&_Mat,irow+1); jcol++) sum += fabs(M__GetVal(&_Mat,irow+1,jcol));
    if (sum > norm) norm = sum;
  }
  return norm;
}
// ======================================================== 
/// This function computes  the l1-infinity: norm=max over (i-rw j-cln).
Real LaspackMMatrixM::linfty_norm() const {
  Real norm = 0;
    for (uint irow = 0; irow < m(); irow++) {
    Real sum = 0;
    for(uint jcol=0; jcol < M__GetLen(&_Mat,irow+1); jcol++) {
      sum = fabs(M__GetVal(&_Mat,irow+1,jcol));  if (sum > norm) norm = sum;
    }
  }
  return norm;
}
// ======================================================== 
void LaspackMMatrixM::init (){
  std::cout << "  LaspackMMatrixM::init: Our LaspackMMatrix  needs dimensions"; abort();
}
// ===================================
void LaspackMMatrixM::get_diagonal (NumericVectorM& /*dest*/) const{
  std::cout << "not_implemented";exit(0);
}
// ===================================
void LaspackMMatrixM::get_transpose (SparseMMatrixM& /*dest*/) const{ 
  std::cout << "not_implemented";exit(0);
}
// =========================================================================
/// This function sets all row entries to 0 then puts diag_value in the diagonal entry
void LaspackMMatrixM::zero_rows(std::vector<int> & rows, Real diag_value) {
  std::cout << "LaspackMatrixM::zero_rows: not implemented \n"; abort();
}
/*
void LaspackMMatrixM::vector_mult (NumericVectorM &vec_out_in,const NumericVectorM &vec_in) const {
  // Make sure the data passed in are really in Laspack types
  LaspackVectorM* Lvec2 = dynamic_cast<LaspackVectorM*>(&vec_out_in);
  const LaspackVectorM* Lvec = static_cast<const LaspackVectorM*>(&vec_in);
  assert(Lvec != NULL);  assert(Lvec2 != NULL);
  // += mat*vec
  const QVector * vec=const_cast<QVector*>(&Lvec->_vec);
  QVector * vec2=static_cast<QVector*>(&Lvec2->_vec);
   assert(mat != NULL);  assert(vec != NULL);
   
  //   AddAsgn_VV (&_vec, Mul_MV(const_cast<Matrix*>(&mat->_Mat),const_cast<QVector*>(&vec->_vec)));
  uint n_row=vec_out.size();
  for (uint irow=1;irow<=n_row;irow++) {
    double sum=0.;
    for (uint j=0;j<M__GetLen(&_Mat,irow);j++) sum +=M__GetVal(&_Mat,irow,j)*V__GetCmp(vec,Q__GetPos(&_Mat,irow,j));
  V__SetCmp(&vec2,irow,sum);
}*/


#endif // #ifdef LIBMESH_HAVE_LASPACK
