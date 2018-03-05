#include "Solverlib_conf.h"

#ifdef HAVE_PETSCM

#include <sstream>
#include "hdf5.h"

// C++ includes
// Local Includes
#include "petsc_vectorM.h"
#include "petsc_matrixM.h"
#include "petsc_MmatrixM.h"
#include "dense_subvectorM.h"
#include "dense_vectorM.h"
#include "sparse_MmatrixM.h"
#include "parallelM.h"
#include "petsc_macroM.h"

#include "MGCasts.h"  // TODO #include "utility.h"  

//-----------------------------------------------------------------------
// PetscVector members

// void PetscVectorM::init (const NumericVectorM& v, const bool fast)
// {
//   libmesh_error();
//   init (v.local_size(), (int)v.size(), fast);
//   vec = libmesh_cast_ref<const PetscVectorM&>(v).vec;
// }


// ============================================
///< This function returns the sum of values in a vector
double PetscVectorM::sum() const {
  this->_restore_array(); assert(this->closed());
  int ierr=0;  PetscScalar value=0.;
  ierr = VecSum(_vec, &value);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ====================================================
/// This function returns the \f$l_1\f$-norm of the vector
double PetscVectorM::l1_norm() const {
  this->_restore_array();  assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_1, &value); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// =============================================
/// This function returns the \f$l_2\f$-norm of the vector
double PetscVectorM::l2_norm() const {
  this->_restore_array();  assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_2, &value); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ============================================
/// This function returns the maximum absolute value of the elements of this vector
double PetscVectorM::linfty_norm() const {
  this->_restore_array(); assert(this->closed());
  int ierr=0;  PetscReal value=0.;
  ierr = VecNorm(_vec, NORM_INFINITY, &value); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// =======================================
NumericVectorM& PetscVectorM::operator += (const NumericVectorM& v) {
  this->_restore_array(); assert(this->closed());
  this->add(1., v);
  return *this;
}

// ============================================================
NumericVectorM& PetscVectorM::operator -= (const NumericVectorM& v) {
  this->_restore_array();  assert(this->closed());
  this->add(-1., v);
  return *this;
}

// =============================================================
void PetscVectorM::set(const  int i, const double value) {
  this->_restore_array();  assert(i<size());
  int ierr=0;  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = VecSetValues(_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}
// ===============================
void PetscVectorM::reciprocal(){
  PetscErrorCode ierr = 0;
  // VecReciprocal has been in PETSc since at least 2.3.3 days
  ierr = VecReciprocal(_vec);CHKERRABORT(MPI_COMM_WORLD,ierr);
}
// ============================================================
void PetscVectorM::add (const  int i, const double value) {
  this->_restore_array();  assert(i<size());
  int ierr=0;  int i_val = static_cast<int>(i);
  PetscScalar petsc_value = static_cast<PetscScalar>(value);

  ierr = VecSetValues(_vec, 1, &i_val, &petsc_value, ADD_VALUES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->_is_closed = false;
}

// ===============================================================
void PetscVectorM::add_vector(const std::vector<double>& v,
                              const std::vector< int>& dof_indices) {
  this->_restore_array(); assert(v.size() == dof_indices.size());
  for(int i=0; i<(int)v.size(); i++)    this->add(dof_indices[i], v[i]);
}

// ===========================================================
void PetscVectorM::add_vector(const NumericVectorM& V,
                              const std::vector< int>& dof_indices) {
  assert((int)V.size() == (int)dof_indices.size());
  for(int i=0; i<(int)V.size(); i++) this->add(dof_indices[i], V(i));
}

// =========================================================
void PetscVectorM::add_vector(const NumericVectorM& V_in,
                              const SparseMatrixM& A_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* V = dynamic_cast<const PetscVectorM*>(&V_in);
  const PetscMatrixM* A = dynamic_cast<const PetscMatrixM*>(&A_in);
  int ierr=0;
  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultAdd(const_cast<PetscMatrixM*>(A)->mat(), V->_vec, _vec, _vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}
// =========================================================
void PetscVectorM::add_vector(const NumericVectorM& V_in,
                              const SparseMMatrixM& A_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* V = dynamic_cast<const PetscVectorM*>(&V_in);
  const PetscMMatrixM* A = dynamic_cast<const PetscMMatrixM*>(&A_in);
  int ierr=0;
  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultAdd(const_cast<PetscMMatrixM*>(A)->mat(), V->_vec, _vec, _vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ====================================================
void PetscVectorM::add_vector(const DenseVectorM& V,
                              const std::vector< int>& dof_indices) {
  assert((int)V.size() == (int)dof_indices.size());
  for(int i=0; i<(int)V.size(); i++)  this->add(dof_indices[i], V(i));
}
// ==========================================
void PetscVectorM::add_vector_transpose (
  const NumericVectorM& V_in,
  const SparseMatrixM& A_in)
{
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* V = dynamic_cast<const PetscVectorM*>(&V_in);
  const PetscMatrixM* A = dynamic_cast<const PetscMatrixM*>(&A_in);

  PetscErrorCode ierr=0;

  A->close();

  // The const_cast<> is not elegant, but it is required since PETSc
  // is not const-correct.
  ierr = MatMultTransposeAdd(const_cast<PetscMatrixM*>(A)->mat(), V->_vec, _vec, _vec);
         LIBMESH_CHKERRABORT(ierr);
}

// ====================================================
void PetscVectorM::matrix_mult(const NumericVectorM &vec_in,const SparseMMatrixM &mat_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&vec_in);
  const PetscMMatrixM* A = static_cast<const PetscMMatrixM*>(&mat_in);
  int ierr=0;  A->close();
  ierr = MatMult(const_cast<PetscMMatrixM*>(A)->mat(),v->_vec,_vec);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}

/// This function computes the residual r=b-Ax
// =========================================================
void PetscVectorM::resid(const NumericVectorM& b_in,const NumericVectorM& x_in,
                         const SparseMatrixM& A_in) {
  this->_restore_array();
  // Make sure the data passed in are really of Petsc types
  const PetscVectorM* b = static_cast<const PetscVectorM*>(&b_in);
  const PetscVectorM* x = static_cast<const PetscVectorM*>(&x_in);
  const PetscMatrixM* A = static_cast<const PetscMatrixM*>(&A_in);
  int ierr=0; /* A->close();*/
  // residual computation r=b-Ax
  ierr = MatMult(const_cast<PetscMatrixM*>(A)->mat(), x->_vec, _vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecAYPX(_vec,-1,b->_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
}

// ====================================================
void PetscVectorM::add(const double v_in) {
  this->_restore_array();  int ierr=0;
  PetscScalar* values;  const PetscScalar v = static_cast<PetscScalar>(v_in);

  if(this->type() != GHOSTEDM)    {
    const int n   = static_cast<int>(this->local_size());
    const int fli = static_cast<int>(this->first_local_index());

    for(int i=0; i<n; i++) {
      ierr = VecGetArray(_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
      int ig = fli + i;  PetscScalar value = (values[i] + v);
      ierr = VecRestoreArray(_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(_vec, 1, &ig, &value, INSERT_VALUES); CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  } else    {
    /* Vectors that include ghost values require a special
    handling.  */
    Vec loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);

    int n=0;
    ierr = VecGetSize(loc_vec, &n);  CHKERRABORT(MPI_COMM_WORLD,ierr);

    for(int i=0; i<n; i++) {
      ierr = VecGetArray(loc_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
      PetscScalar value = (values[i] + v);
      ierr = VecRestoreArray(loc_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecSetValues(loc_vec, 1, &i, &value, INSERT_VALUES);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  this->_is_closed = false;
}

// ==============================================
void PetscVectorM::add(const NumericVectorM& v) {
  this->add(1., v);
}

// ====================================================
void PetscVectorM::add(const double a_in, const NumericVectorM& v_in) {
  this->_restore_array();
  int ierr = 0;  PetscScalar a = static_cast<PetscScalar>(a_in);

  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&v_in);
  v->_restore_array();  assert(this->size() == v->size());
  if(this->type() != GHOSTEDM)  {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecAXPY(&a, v->_vec, _vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else  // 2.3.x & later style
    ierr = VecAXPY(_vec, a, v->_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
  } else    {
    Vec loc_vec;     Vec v_loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);      CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(v->_vec,&v_loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecAXPY(&a, v_loc_vec, loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else    // 2.3.x & later style
    ierr = VecAXPY(loc_vec, a, v_loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    ierr = VecGhostRestoreLocalForm(v->_vec,&v_loc_vec);   CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}


// ========================================================
void PetscVectorM::insert(const std::vector<double>& v,
                          const std::vector< int>& dof_indices) {
  assert(v.size() == dof_indices.size());
  for(int i=0; i<(int)v.size(); i++)    this->set(dof_indices[i], v[i]);
}

// =======================================================
void PetscVectorM::insert(const NumericVectorM& V,
                          const std::vector< int>& dof_indices) {
  assert((int)V.size() == (int)dof_indices.size());
  for(int i=0; i<(int)V.size(); i++)   this->set(dof_indices[i], V(i));
}

// ========================================================
void PetscVectorM::insert(const DenseVectorM& V,
                          const std::vector< int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for(int i=0; i<(int)V.size(); i++)  this->set(dof_indices[i], V(i));
}

// =========================================================
void PetscVectorM::insert(const DenseSubVectorM& V,
                          const std::vector< int>& dof_indices) {
  assert(V.size() == dof_indices.size());
  for(int i=0; i<(int)V.size(); i++)  this->set(dof_indices[i], V(i));
}

// ================================================
void PetscVectorM::scale(const double factor_in) {
  this->_restore_array();  int ierr = 0;
  PetscScalar factor = static_cast<PetscScalar>(factor_in);

  if(this->type() != GHOSTEDM)    {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecScale(&factor, _vec);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
    // 2.3.x & later style
    ierr = VecScale(_vec, factor);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
  } else   {
    Vec loc_vec; ierr = VecGhostGetLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);

// #if PETSC_VERSION_LESS_THAN(2,3,0)
//       // 2.2.x & earlier style
//       ierr = VecScale(&factor, loc_vec);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
    // 2.3.x & later style
    ierr = VecScale(loc_vec, factor); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}

// ================================================
void PetscVectorM::abs() {
  this->_restore_array();
  int ierr = 0;
  if(this->type() != GHOSTEDM)   {
    ierr = VecAbs(_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else    {
    Vec loc_vec; ierr = VecGhostGetLocalForm(_vec,&loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecAbs(loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
}

// ===========================================
double PetscVectorM::dot(const NumericVectorM& V) const {
  this->_restore_array();
  // Error flag and Return value
  int ierr = 0; PetscScalar value=0.;
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&V);
  // 2.3.x (at least) style.  Untested for previous versions.
  ierr = VecDot(this->_vec, v->_vec, &value);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return static_cast<double>(value);
}

// ===============================================
NumericVectorM& PetscVectorM::operator = (const double s_in) {
  this->_restore_array();  assert(this->closed());

  int ierr = 0;  PetscScalar s = static_cast<PetscScalar>(s_in);

  if(this->size() != 0)    {
    if(this->type() != GHOSTEDM)  {
// #if PETSC_VERSION_LESS_THAN(2,3,0)
//    // 2.2.x & earlier style
//    ierr = VecSet(&s, _vec);
//    CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
      // 2.3.x & later style
      ierr = VecSet(_vec, s); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif
    } else {
      Vec loc_vec; ierr = VecGhostGetLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);

// #if PETSC_VERSION_LESS_THAN(2,3,0)
//    // 2.2.x & earlier style
//    ierr = VecSet(&s, loc_vec);
//    CHKERRABORT(MPI_COMM_WORLD,ierr);
// #else
      // 2.3.x & later style
      ierr = VecSet(loc_vec, s); CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif

      ierr = VecGhostRestoreLocalForm(_vec,&loc_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  }

  return *this;
}



// =========================================================
NumericVectorM& PetscVectorM::operator = (const NumericVectorM& v_in) {
  // Make sure the NumericVector passed in is really a PetscVector
  const PetscVectorM* v = static_cast<const PetscVectorM*>(&v_in);
  *this = *v;
  return *this;
}

// =====================================================
PetscVectorM& PetscVectorM::operator = (const PetscVectorM& v) {
  this->_restore_array();  v._restore_array();
  assert(this->_type == v._type);  assert(this->size() == (int)v.size());
  assert(this->local_size() == v.local_size());  assert(this->_global_to_local_map == v._global_to_local_map);
  
    int ierr = 0;
    if (((this->type()==PARALLELM) && (v.type()==GHOSTEDM)) ||
      ((this->type()==GHOSTEDM) && (v.type()==PARALLELM)) ||
      ((this->type()==GHOSTEDM) && (v.type()==SERIALM))   ||
      ((this->type()==SERIALM) && (v.type()==GHOSTEDM)))
    {
    /* Allow assignment of a ghosted to a parallel vector since this
       causes no difficulty.  See discussion in libmesh-devel of
       June 24, 2010.  */
    ierr = VecCopy (v._vec, this->_vec);CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  else
  {
    /* In all other cases, we assert that both vectors are of equal
       type.  */
    assert (this->_type== v._type);
    assert (this->_global_to_local_map == v._global_to_local_map);

    if((int)v.size() != 0)    {
    if(this->type() != GHOSTEDM)  
    {
      ierr = VecCopy(v._vec, this->_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    } else {
      Vec loc_vec;    Vec v_loc_vec;
      ierr = VecGhostGetLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostGetLocalForm(v._vec,&v_loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecCopy(v_loc_vec, loc_vec);    CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm(v._vec,&v_loc_vec);    CHKERRABORT(MPI_COMM_WORLD,ierr);
      ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
  }
  }
  
  close();
  return *this;
}

// =====================================================
NumericVectorM& PetscVectorM::operator = (const std::vector<double>& v) {
  this->_restore_array();
  const  int nl   = this->local_size();
  const  int ioff = this->first_local_index();
  int ierr=0;
  PetscScalar* values;

  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if(this->size() == (int)v.size())    {
    ierr = VecGetArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    for(int i=0; i<nl; i++)  values[i] =  static_cast<PetscScalar>(v[i+ioff]);
    ierr = VecRestoreArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else    {
    assert(this->local_size() == (int)v.size());
    ierr = VecGetArray(_vec, &values);   CHKERRABORT(MPI_COMM_WORLD,ierr);
    for(int i=0; i<nl; i++)  values[i] = static_cast<PetscScalar>(v[i]);
    ierr = VecRestoreArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // Make sure ghost dofs are up to date
  if(this->type() == GHOSTEDM)   this->close();

  return *this;
}

// ============================================================
/// This function localizes data from any processor to a local vector
void PetscVectorM::localize(
  NumericVectorM& v_local_in  //  local vector
) const { // ===================================================

  this->_restore_array();
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVectorM* v_local = dynamic_cast<PetscVectorM*>(&v_local_in);
  assert(v_local != NULL); assert(v_local->size() == this->size());
  int ierr = 0;  const int n = this->size();
  IS is;  VecScatter scatter;
  // Create idx, idx[i] = i;
  std::vector<int> idx(n); iota(idx.begin(), idx.end(), 0);
//   Utility::iota (idx.begin(), idx.end(), 0);

  // Create the index set & scatter object
  ierr = ISCreateGeneral(this->comm().get(), n, &idx[0], PETSC_USE_POINTER, &is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecScatterCreate(_vec,   is,v_local->_vec, is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,SCATTER_FORWARD, scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = VecScatterEnd(_vec, v_local->_vec, INSERT_VALUES,SCATTER_FORWARD, scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, _vec, v_local->_vec,
                       INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif

  // Clean up
  ierr = ISDestroy(&is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Make sure ghost dofs are up to date
  if(v_local->type() == GHOSTEDM)   v_local->close();
  return;
}

// ==========================================================
/// This function localizes data from any processor to a local vector
void PetscVectorM::localize(
  NumericVectorM& v_local_in,         //  local vector
  const std::vector< int>& send_list  //  trasfer data vactor
) const { // ===========================================
  // Workaround for a strange bug at large-scale.
  // If we have ghosting, PETSc lets us just copy the solution, and
  // doing so avoids a segfault?
  if(v_local_in.type() == GHOSTEDM &&  this->type() == PARALLELM) {  v_local_in = *this; return;}
  this->_restore_array();
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVectorM* v_local = dynamic_cast<PetscVectorM*>(&v_local_in);

  assert(v_local != NULL);  assert((int)v_local->size() == (int)this->size());
  assert((int)send_list.size() <= (int)v_local->size());

  int ierr=0;  const  int n_sl = send_list.size();
  IS is;  VecScatter scatter;

  std::vector<int> idx(n_sl + this->local_size());
  for(int i=0; i<n_sl; i++)   idx[i] = static_cast<int>(send_list[i]);
  for(int i = 0; i != this->local_size(); ++i)   idx[n_sl+i] = i + this->first_local_index();

  // Create the index set & scatter object
  if(idx.empty())  ierr = ISCreateGeneral(this->comm().get(),n_sl+this->local_size(),
                                            PETSC_NULL, PETSC_USE_POINTER, &is);
  else  ierr = ISCreateGeneral(this->comm().get(),n_sl+this->local_size(),
                                 &idx[0],  PETSC_USE_POINTER,&is);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterCreate(_vec,is,v_local->_vec, is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
  ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
                         SCATTER_FORWARD, scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(_vec, v_local->_vec, INSERT_VALUES,
                       SCATTER_FORWARD, scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#else  // API argument order change in PETSc 2.3.3
  ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, _vec, v_local->_vec,
                       INSERT_VALUES, SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
  // Clean up
  ierr = ISDestroy(&is); CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Make sure ghost dofs are up to date
  if(v_local->type() == GHOSTEDM)  v_local->close();
  return;

}

// ================================================================
void PetscVectorM::localize(
  const  int first_local_idx,  // first index
  const  int last_local_idx,   // last index
  const std::vector< int>& send_list // index list
) { // ==================================

  assert((int)send_list.size() <= this->size());
  assert((int)last_local_idx+1 <= this->size());
  // -------------------------------------------
  const unsigned int size       = this->size();
  const unsigned int local_size = (last_local_idx - first_local_idx + 1);
  int ierr=0; //int npi=1;

  // Don't bother for serial cases
//  if ((first_local_idx == 0) &&
//      (my_local_size == my_size))
  // But we do need to stay in sync for degenerate cases
  if (this->n_processors() == 1)    return;


  // Build a parallel vector, initialize it with the local
  // parts of (*this)
  PetscVectorM parallel_vec(this->comm(), PARALLELM);

  parallel_vec.init (size, local_size, true, PARALLELM);
  
  
  // But we do need to stay in sync for degenerate cases
//    MPI_Comm_size(MPI_COMM_WORLD, &npi);
//    if(npi == 1)    return;
  // Build a parallel vector, initialize it with the local parts of (*this)
//   PetscVectorM parallel_vec;
//   parallel_vec.init(size, local_size, true, PARALLELM);

  {
    // Copy part of *this into the parallel_vec -------------
    IS is;   VecScatter scatter;
    // Create idx, idx[i] = i+first_local_idx;
    std::vector<int> idx(local_size);
    iota(idx.begin(), idx.end(), first_local_idx);
//     Utility::iota (idx.begin(), idx.end(), first_local_idx);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(this->comm().get(), local_size,
                           local_size ? &idx[0] : NULL, PETSC_USE_POINTER, &is);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterCreate(_vec,is, parallel_vec._vec, is, &scatter);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Perform the scatter
#if PETSC_VERSION_LESS_THAN(2,3,3)
    ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
                           SCATTER_FORWARD, scatter);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterEnd(_vec, parallel_vec._vec, INSERT_VALUES,
                         SCATTER_FORWARD, scatter);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
#else   // API argument order change in PETSc 2.3.3
    ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,
                           INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterEnd(scatter, _vec, parallel_vec._vec,
                         INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
#endif
    // Clean up
    ierr = ISDestroy(&is);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  } // -------------------------------------------------------

  // localize like normal
  parallel_vec.close(); parallel_vec.localize(*this, send_list);
  this->close();
  return;
}



// =============================================================
void PetscVectorM::localize(std::vector<double>& v_local) const {
  this->_restore_array();

  // This function must be run on all processors at once
  parallel_object_onlyM();

  int ierr=0; PetscScalar *values;
  const int n = this->size();  const int nl = this->local_size();
  v_local.clear();  v_local.resize(n, 0.);

  ierr = VecGetArray(_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);

  int ioff = first_local_index();
  for(int i=0; i<nl; i++)   v_local[i+ioff] = static_cast<double>(values[i]);
  ierr = VecRestoreArray(_vec, &values); CHKERRABORT(MPI_COMM_WORLD,ierr);
  this->comm().sum(v_local);
  return;
}

//  ===========================================================
void PetscVectorM::localize_to_one(
  std::vector<double>& v_local,
  const  int pid
) const {
  this->_restore_array();
  int ierr=0; PetscScalar *values;
  const int n  = size(); const int nl = local_size();

  v_local.resize(n);
  // only one processor
  if(n == nl)    {
    ierr = VecGetArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    for(int i=0; i<n; i++) v_local[i] = static_cast<double>(values[i]);
    ierr = VecRestoreArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  // otherwise multiple processors
  else    {
    int ioff = this->first_local_index();
    std::vector<double> local_values(n, 0.);
    {
      ierr = VecGetArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
      for(int i=0; i<nl; i++) local_values[i+ioff] = static_cast<double>(values[i]);
      ierr = VecRestoreArray(_vec, &values);  CHKERRABORT(MPI_COMM_WORLD,ierr);
    }
    MPI_Reduce(&local_values[0], &v_local[0],n,MPI_REAL,MPI_SUM,pid,this->comm().get());
  }
  return;
}

// =======================================================
void PetscVectorM::pointwise_mult(const NumericVectorM& vec1,
                                  const NumericVectorM& vec2) {
  this->_restore_array();

  int ierr = 0;

  // Convert arguments to PetscVector*.
  const PetscVectorM* vec1_petsc = static_cast<const PetscVectorM*>(&vec1);
  const PetscVectorM* vec2_petsc = static_cast<const PetscVectorM*>(&vec2);


  if(this->type() != GHOSTEDM) {
    ierr = VecPointwiseMult(this->vec(),
                            const_cast<PetscVectorM*>(vec1_petsc)->vec(),
                            const_cast<PetscVectorM*>(vec2_petsc)->vec());
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    Vec loc_vec;
    Vec v1_loc_vec;
    Vec v2_loc_vec;
    ierr = VecGhostGetLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(const_cast<PetscVectorM*>(vec1_petsc)->vec(),&v1_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostGetLocalForm(const_cast<PetscVectorM*>(vec2_petsc)->vec(),&v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecPointwiseMult(loc_vec,v1_loc_vec,v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecGhostRestoreLocalForm(const_cast<PetscVectorM*>(vec1_petsc)->vec(),&v1_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(const_cast<PetscVectorM*>(vec2_petsc)->vec(),&v2_loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostRestoreLocalForm(_vec,&loc_vec);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }


}

#if HDF5_VERSION==188
// ======================================================
/// Print the contents of the vector to hdf5 file,
void PetscVectorM::print_hdf5(std::string file) const {
  this->_restore_array(); assert(this->closed());
  int  ndim=this->size();
  std::ostringstream name;
  name.str("");    name << file.c_str() << ".h5";
  hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  // Value storage
//   double *val=new double[ndim];  for (int  i=0;i<ndim;i++)  val[i]= V__GetCmp(&_vec,i+1);

// print matrix values
  hsize_t dimsf[2]; dimsf[0]=1;    dimsf[1] = 1;
  int *nd=new int[1]; nd[0]=ndim;
  name.str(""); name << "DIM" ;
  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                            dataspace, H5P_DEFAULT);
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,nd);
  assert(status==0);
  H5Sclose(dataspace);  H5Dclose(dataset);

  dimsf[0]=ndim;
  name.str(""); name << "VAL" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_DOUBLE,
                      dataspace, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, &(_vec));
  H5Sclose(dataspace);  H5Dclose(dataset);


  H5Fclose(fileP);
  /* delete []val; */delete []nd;
  return;
}
#else
// ======================================================
/// Print the contents of the vector to hdf5 file,
void PetscVectorM::print_hdf5(std::string file) const {
  this->_restore_array(); assert(this->closed());
  int  ndim=this->size();
  std::ostringstream name;
  name.str("");    name << file.c_str() << ".h5";
  hid_t fileP = H5Fcreate(name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  // Value storage
//   double *val=new double[ndim];  for (int  i=0;i<ndim;i++)  val[i]= V__GetCmp(&_vec,i+1);

// print matrix values
  hsize_t dimsf[2]; dimsf[0]=1;    dimsf[1] = 1;
  int *nd=new int[1]; nd[0]=ndim;
  name.str(""); name << "DIM" ;
  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_INT,
                            dataspace , H5P_DEFAULT
   #if HDF5_VERSIONM!=1808                        
                           , H5P_DEFAULT, H5P_DEFAULT                    
    #endif
  );
  hid_t status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,nd);
  assert(status==0);
  H5Sclose(dataspace);  H5Dclose(dataset);

  dimsf[0]=ndim;
  name.str(""); name << "VAL" ;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset = H5Dcreate(fileP,name.str().c_str(),H5T_NATIVE_DOUBLE,
                      dataspace, H5P_DEFAULT
  #if HDF5_VERSIONM!=1808                        
                           , H5P_DEFAULT, H5P_DEFAULT                    
    #endif                      
                     );
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, &(_vec));
  H5Sclose(dataspace);  H5Dclose(dataset);


  H5Fclose(fileP);
  /* delete []val; */delete []nd;
  return;
}
#endif
// void PetscVectorM::print_hdf5(const std::string name) const{
//   this->_restore_array(); assert (this->closed());
//
//   int ierr=0;
//   PetscViewer petsc_viewer;
//   ierr = PetscViewerCreate (MPI_COMM_WORLD, &petsc_viewer);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);
//
//   // Create an ASCII file containing the matrix if a filename was provided.
//   if (name != "NULL"){
//       ierr = PetscViewerASCIIOpen( MPI_COMM_WORLD,
//           name.c_str(),&petsc_viewer);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = PetscViewerSetFormat (petsc_viewer,
//           PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = VecView (_vec, petsc_viewer);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
//
//   // Otherwise the matrix will be dumped to the screen.
//   else {
//       ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
//           PETSC_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//       ierr = VecView (_vec, PETSC_VIEWER_STDOUT_WORLD);
//              CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
//   // Destroy the viewer.
//   ierr = PetscViewerDestroy (petsc_viewer);
//          CHKERRABORT(MPI_COMM_WORLD,ierr);


// ========================================================
void PetscVectorM::print_personal(std::ostream&/* name*/) const {
  this->_restore_array();
  assert(this->closed());

  int ierr=0;
  PetscViewer petsc_viewer;
  ierr = PetscViewerCreate(MPI_COMM_WORLD,&petsc_viewer);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
                              PETSC_VIEWER_ASCII_MATLAB);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecView(_vec, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Destroy the viewer.
  ierr = PetscViewerDestroy(&petsc_viewer);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}



// ==============================================================
/// This function creates a subvector
void PetscVectorM::create_subvector(
  NumericVectorM& subvector,     //  subvector                             
  const std::vector< int>& rows  //  row index vector
) const {// ======================
  this->_restore_array();

  // PETSc data structures
  IS parent_is, subvector_is;
  VecScatter scatter;
  int ierr = 0;

  // Make sure the passed in subvector is really a PetscVector
  PetscVectorM* petsc_subvector = static_cast<PetscVectorM*>(&subvector);

  // If the petsc_subvector is already initialized, we assume that the
  // user has already allocated the *correct* amount of space for it.
  // If not, we use the appropriate PETSc routines to initialize it.
  if(!petsc_subvector->initialized())  {
    // Initialize the petsc_subvector to have enough space to hold
    // the entries which will be scattered into it.  Note: such an
    // init() function (where we let PETSc decide the number of local
    // entries) is not currently offered by the PetscVector
    // class.  Should we differentiate here between sequential and
    // parallel vector creation based on libMesh::n_processors() ?
    ierr = VecCreateMPI(this->comm().get(),
                        PETSC_DECIDE,          // n_local
                        rows.size(),           // n_global
                        &(petsc_subvector->_vec)); CHKERRABORT(MPI_COMM_WORLD,ierr);

    ierr = VecSetFromOptions(petsc_subvector->_vec); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Mark the subvector as initialized
    petsc_subvector->_is_initialized = true;
  }
  else { petsc_subvector->_restore_array(); }

  // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
  std::vector<int> idx(rows.size());
  iota(idx.begin(), idx.end(), 0);
//   Utility::iota (idx.begin(), idx.end(), 0);

  // Construct index sets
  ierr = ISCreateGeneral(this->comm().get(),rows.size(),(int*) &rows[0],
			 PETSC_USE_POINTER,&parent_is); 
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(this->comm().get(),rows.size(),(int*) &idx[0], 
			 PETSC_USE_POINTER,&subvector_is); 
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Construct the scatter object
  ierr = VecScatterCreate(this->_vec,parent_is,petsc_subvector->_vec,
                          subvector_is,&scatter);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Actually perform the scatter
  ierr = VecScatterBegin(scatter, this->_vec, petsc_subvector->_vec,
                         INSERT_VALUES, SCATTER_FORWARD);               
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterEnd(scatter, this->_vec, petsc_subvector->_vec,
                       INSERT_VALUES,SCATTER_FORWARD);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Clean up
  ierr = ISDestroy(&parent_is);       CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = ISDestroy(&subvector_is);    CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRABORT(MPI_COMM_WORLD,ierr);
  return;
}




//------------------------------------------------------------------
// Explicit instantiations
// template class PetscVector<Number>;



#endif // #ifdef LIBMESH_HAVE_PETSC
