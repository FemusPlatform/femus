#ifndef TYPE_VECTOR_H
#define TYPE_VECTOR_H

// Local includes
// #include "libmesh/libmesh_common.h"
// #include "libmesh/compare_types.h"
// #include "libmesh/tensor_tools.h"
#include <assert.h> 
// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
// #include <complex>
#include <iostream>
#include "Domain_conf.h"
// namespace libMesh
// {

// Forward declarations
//   class TypeTensor;
//   class VectorValue;
//   class TensorValue;
// #define DIMENSION 2

class VectorND{
// template <typename T2>
// friend class VectorND;

// friend class TypeTensor<T>;

protected:
  /// Empty constructor.  Gives the vector 0 in \p DIMENSION
  VectorND  ();
  /// Constructor-from-T.  By default sets higher dimensional
  VectorND (const double x,const double y=0,const double z=0);

//   ///
//    * Constructor-from-scalars.  By default sets higher dimensional
//    * entries to 0.
//    */
//   template <typename Scalar>
//   VectorND (const Scalar x,
// 	      const Scalar y=0,
//               typename
//                 boostcopy::enable_if_c<ScalarTraits<Scalar>::value,
//                                        const Scalar>::type z=0);

public:

  /// Copy-constructor.
  VectorND (const VectorND & p);

  /// Destructor.
  ~VectorND ();

  /// Assign to a vector without creating a temporary.
//   template <typename T2>
  void assign (const VectorND &);

  /**
   * Assignment-from-scalar operator.  Used only to zero out vectors.
   */
//   template <typename Scalar>
//   typename boostcopy::enable_if_c<
//     ScalarTraits<Scalar>::value,
//     VectorND&>::type
//   operator = (const Scalar& p)
//   { libmesh_assert_equal_to (p, Scalar(0)); this->zero(); return *this; }

  /**
   * Return the \f$ i^{th} \f$ element of the vector.
   */
  const double & operator () (const unsigned int i) const;
  const double & slice (const unsigned int i) const { return (*this)(i); }

  /// Return a writeable reference to the \f$ i^{th} \f$ element of the vector.
  double & operator () (const unsigned int i);

  double & slice (const unsigned int i) { return (*this)(i); }

  // Add two vectors.
  VectorND  operator + (const VectorND &) const;

  /// Add to this vector.
//   template <typename T2>
  const VectorND & operator += (const VectorND &);

  /// Add to this vector without creating a temporary.
//   template <typename T2>
  void add (const VectorND &);

  /**
   * Add a scaled value to this vector without
   * creating a temporary.
   */
//   template <typename T2>
  void add_scaled (const VectorND &, const double);

  /**
   * Subtract two vectors.
   */
//   template <typename T2>
  VectorND  operator - (const VectorND &) const;

  /// Subtract from this vector.
//   template <typename T2>
  const VectorND & operator -= (const VectorND &);

  /// Subtract from this vector without creating a temporary.
  void subtract (const VectorND &);
  /// Subtract a scaled value from this vector
  void subtract_scaled (const VectorND &, const double);
  /// Return the opposite of a vector
  VectorND operator - () const;

  /// Multiply a vector by a number, i.e. scale.
//   template <typename Scalar>
//   typename boostcopy::enable_if_c<
//     ScalarTraits<Scalar>::value,
//     VectorND<typename CompareTypes<T, Scalar>::supertype> >::type
//   operator * (const Scalar) const;

  /// Multiply this vector by a number, i.e. scale.
  const VectorND & operator *= (const double);

//   /// Divide a vector by a number, i.e. scale.
//   template <typename Scalar>
//   typename boostcopy::enable_if_c<
//     ScalarTraits<Scalar>::value,
//     VectorND<typename CompareTypes<T, Scalar>::supertype> >::type
//   operator / (const Scalar) const;

  /**
   * Divide this vector by a number, i.e. scale.
   */
  const VectorND & operator /= (const double);

//   /// Multiply 2 vectors together, i.e. dot-product.
//   /// The vectors may be of different types.
//    
//   template <typename T2>
//   typename CompareTypes<T, T2>::supertype
//   operator * (const VectorND<T2> &) const;

  /**
   * Multiply 2 vectors together, i.e. dot-product.
   * The vectors may be of different types.
   */
//   template <typename T2>
//   typename CompareTypes<T, T2>::supertype  contract (const VectorND<T2> &) const;

  /**
   * Cross 2 vectors together, i.e. cross-product.
   */
//   template <typename T2>
  VectorND
//   <typename CompareTypes<T, T2>::supertype> 
  cross(const VectorND &) const;

  /// Think of a vector as a \p dim dimensional vector.  
  VectorND unit() const;

  /// Returns the magnitude of the vector, i.e. the square-root 
  double size() const;
  /// Returns the magnitude of the vector squared, i.e. the sum of the element magnitudes squared.
  double size_sq() const;
  /// Zero the vector in any dimension.
  void zero();

//   /**
//    * @returns \p true iff two vectors occupy approximately the same
//    * physical location in space, to within a relative tolerance of \p tol.
//    */
//   bool relative_fuzzy_equals(const VectorND& rhs, Real tol = TOLERANCE) const;
// 
//   /**
//    * @returns \p true iff two vectors occupy approximately the same
//    * physical location in space, to within an absolute tolerance of \p tol.
//    */
//   bool absolute_fuzzy_equals(const VectorND& rhs, Real tol = TOLERANCE) const;

  /// @returns \p true iff two vectors occupy approximately the same
  /// physical location in space, to within an absolute tolerance of \p TOLERANCE.  
  bool operator == (const VectorND& rhs) const;
  bool operator != (const VectorND& rhs) const;

  /**
   * @returns \p true if this vector is "less"
   * than another.  Useful for sorting.
   * Also used for choosing some arbitrary basis function
   * orientations
   */
  bool operator < (const VectorND& rhs) const;

  /**
   * @returns \p true if this vector is "less"
   * than or equal to another.  Useful for sorting.
   * Also used for choosing some arbitrary constraint
   * equation directions
   */
  bool operator <= (const VectorND& rhs) const;

  /**
   * @returns \p true if this vector is "greater"
   * than another.  Useful for sorting.
   * Also used for choosing some arbitrary basis function
   * orientations
   */
  bool operator > (const VectorND& rhs) const;

  /**
   * @returns \p true if this vector is "greater"
   * than or equal to another.  Useful for sorting.
   * Also used for choosing some arbitrary constraint
   * equation directions
   */
  bool operator >= (const VectorND& rhs) const;

  /**
   * Formatted print, by default to \p libMesh::out.
   */
  void print(std::ostream& os = std::cout) const;

  /// Formatted print as above but allows you to do
  /// Point p(1,2,3); std::cout << p << std::endl;
  friend std::ostream& operator << (std::ostream& os, const VectorND& t){
    t.print(os);  return os;
  }

  /**
   * Unformatted print to the stream \p out.  Simply prints the elements
   * of the vector separated by spaces.  Optionally prints a newline,
   * which it does by default.
   */
  void write_unformatted (std::ostream &out, const bool newline = true) const;

//  protected:

  /// The coordinates of the \p VectorND
  double _coords[DIMENSION];
};





//------------------------------------------------------
// Inline functions
inline VectorND::VectorND (){
  _coords[0] = 0;
#if DIMENSION > 1
  _coords[1] = 0;
#endif
#if DIMENSION > 2
  _coords[2] = 0;
#endif
}
// ==========================================
inline VectorND::VectorND (
  const double x,const double y,const double z){
  _coords[0] = x;
#if DIMENSION > 1
  _coords[1] = y;
#else
  assert(y==0);
#endif
#if DIMENSION > 2
  _coords[2] = z;
#else
  assert(z==0);
#endif
}
// ========================================================
inline VectorND::VectorND (const VectorND &p){
  // copy the nodes from vector p to me
  for (unsigned int i=0; i< DIMENSION; i++) _coords[i] = p._coords[i];
}
// ========================================================
inline VectorND::~VectorND (){}

// ========================================================
inline void VectorND::assign (const VectorND &p){
  for (unsigned int i=0; i<DIMENSION; i++) _coords[i] = p._coords[i];
}



 // =========================================================
inline const double & VectorND::operator () (const unsigned int i) const{
  assert(i< DIMENSION);
  return _coords[i];
}

// ===========================================================
inline
double & VectorND::operator () (const unsigned int i){
  assert (i< DIMENSION);
  return _coords[i];
}

// // ===========================================================
// 
// inline
// VectorND<typename CompareTypes<T, T2>::supertype>
// VectorND::operator + (const VectorND<T2> &p) const
// {
//   typedef typename CompareTypes<T, T2>::supertype TS;
// #if DIMENSION == 1
//   return VectorND<TS> (_coords[0] + p._coords[0]);
// #endif
// 
// #if DIMENSION == 2
//   return VectorND<TS> (_coords[0] + p._coords[0],
//                          _coords[1] + p._coords[1]);
// #endif
// 
// #if DIMENSION == 3
//   return VectorND<TS> (_coords[0] + p._coords[0],
//                          _coords[1] + p._coords[1],
//                          _coords[2] + p._coords[2]);
// #endif
// 
// }



 
// template <typename T2>
inline const VectorND & VectorND::operator += (const VectorND &p){
  this->add (p);
  return *this;
}

// ======================================================================
inline void VectorND::add (const VectorND &p){
#if DIMENSION == 1
  _coords[0] += p._coords[0];
#endif

#if DIMENSION == 2
  _coords[0] += p._coords[0];
  _coords[1] += p._coords[1];
#endif

#if DIMENSION == 3
  _coords[0] += p._coords[0];
  _coords[1] += p._coords[1];
  _coords[2] += p._coords[2];
#endif

}

// ========================================================
inline
void VectorND::add_scaled (const VectorND &p, const double factor){
#if DIMENSION == 1
  _coords[0] += factor*p(0);
#endif
#if DIMENSION == 2
  _coords[0] += factor*p(0);
  _coords[1] += factor*p(1);
#endif
#if DIMENSION == 3
  _coords[0] += factor*p(0);
  _coords[1] += factor*p(1);
  _coords[2] += factor*p(2);
#endif
}


/*
 
template <typename T2>
inline
VectorND<typename CompareTypes<T, T2>::supertype>
VectorND::operator - (const VectorND<T2> &p) const
{
  typedef typename CompareTypes<T, T2>::supertype TS;

#if DIMENSION == 1
  return VectorND<TS>(_coords[0] - p._coords[0]);
#endif

#if DIMENSION == 2
  return VectorND<TS>(_coords[0] - p._coords[0],
		        _coords[1] - p._coords[1]);
#endif

#if DIMENSION == 3
  return VectorND<TS>(_coords[0] - p._coords[0],
		        _coords[1] - p._coords[1],
		        _coords[2] - p._coords[2]);
#endif

}*/



// ============================================================
inline
const VectorND & VectorND::operator -= (const VectorND &p){
  this->subtract (p); return *this;
}
// ==========================================================
inline void VectorND::subtract (const VectorND& p){
  for (unsigned int i=0; i<DIMENSION; i++)  _coords[i] -= p._coords[i];
}
// ==========================================================
inline void VectorND::subtract_scaled (const VectorND &p, const double factor){
  for (unsigned int i=0; i<DIMENSION; i++)   _coords[i] -= factor*p(i);
}
//=========================================================
inline VectorND VectorND::operator - () const{
#if DIMENSION == 1
  return VectorND(-_coords[0]);
#endif
#if DIMENSION == 2
  return VectorND(-_coords[0],
		    -_coords[1]);
#endif
#if DIMENSION == 3
  return VectorND(-_coords[0],
		    -_coords[1],
		    -_coords[2]);
#endif
}
// ==================================================


 
// template <typename Scalar>
// inline
// typename boostcopy::enable_if_c<
//   ScalarTraits<Scalar>::value,
//   VectorND<typename CompareTypes<T, Scalar>::supertype> >::type
// VectorND::operator * (const Scalar factor) const
// {
//   typedef typename CompareTypes<T, Scalar>::supertype SuperType;
// 
// #if DIMENSION == 1
//   return VectorND<SuperType>(_coords[0]*factor);
// #endif
// 
// #if DIMENSION == 2
//   return VectorND<SuperType>(_coords[0]*factor,
// 		        _coords[1]*factor);
// #endif
// 
// #if DIMENSION == 3
//   return VectorND<SuperType>(_coords[0]*factor,
// 		        _coords[1]*factor,
// 		        _coords[2]*factor);
// #endif
// }



// template <typename T, typename Scalar>
// inline
// typename boostcopy::enable_if_c<
//   ScalarTraits<Scalar>::value,
//   VectorND<typename CompareTypes<T, Scalar>::supertype> >::type
// operator * (const Scalar factor,
//             const VectorND &v)
// {
//   return v * factor;
// }



 
inline const VectorND & VectorND::operator *= (const double factor){
#if DIMENSION == 1
  _coords[0] *= factor;
#endif
#if DIMENSION == 2
  _coords[0] *= factor;
  _coords[1] *= factor;
#endif
#if DIMENSION == 3
  _coords[0] *= factor;
  _coords[1] *= factor;
  _coords[2] *= factor;
#endif
  return *this;
}



 
// template <typename Scalar>
// inline
// typename boostcopy::enable_if_c<
//   ScalarTraits<Scalar>::value,
//   VectorND<typename CompareTypes<T, Scalar>::supertype> >::type
// VectorND::operator / (const Scalar factor) const
// {
//   libmesh_assert_not_equal_to (factor, static_cast<T>(0.));
// 
//   typedef typename CompareTypes<T, Scalar>::supertype TS;
// 
// #if DIMENSION == 1
//   return VectorND<TS>(_coords[0]/factor);
// #endif
// 
// #if DIMENSION == 2
//   return VectorND<TS>(_coords[0]/factor,
// 		        _coords[1]/factor);
// #endif
// 
// #if DIMENSION == 3
//   return VectorND<TS>(_coords[0]/factor,
// 		        _coords[1]/factor,
// 		        _coords[2]/factor);
// #endif
// 
// }




 
inline const VectorND & VectorND::operator /= (const double factor){
  assert (factor != static_cast<double>(0.));
  for (unsigned int i=0; i<DIMENSION; i++)  _coords[i] /= factor;
  return *this;
}
// ========================================================



 
// template <typename T2>
// inline
// typename CompareTypes<T, T2>::supertype
// VectorND::operator * (const VectorND<T2> &p) const
// {
// #if DIMENSION == 1
//   return _coords[0]*p._coords[0];
// #endif
// 
// #if DIMENSION == 2
//   return (_coords[0]*p._coords[0] +
// 	  _coords[1]*p._coords[1]);
// #endif
// 
// #if DIMENSION == 3
//   return (_coords[0]*p(0) +
// 	  _coords[1]*p(1) +
// 	  _coords[2]*p(2));
// #endif
// }

 
// template <typename T2>
// inline
// typename CompareTypes<T, T2>::supertype
// VectorND::contract(const VectorND<T2> &p) const
// {
//   return (*this)*(p);
// }



 
// template <typename T2>
// VectorND<typename CompareTypes<double, double>::supertype>
// VectorND::cross(const VectorND<T2>& p) const
// {
//   typedef typename CompareTypes<T, T2>::supertype TS;
//   libmesh_assert_equal_to (DIMENSION, 3);
// 
//   // |     i          j          k    |
//   // |(*this)(0) (*this)(1) (*this)(2)|
//   // |   p(0)       p(1)       p(2)   |
// 
//   return VectorND<TS>( _coords[1]*p._coords[2] - _coords[2]*p._coords[1],
//                         -_coords[0]*p._coords[2] + _coords[2]*p._coords[0],
//                          _coords[0]*p._coords[1] - _coords[1]*p._coords[0]);
// }



// =============================================== 
inline double VectorND::size() const{
  return std::sqrt(this->size_sq());
}
// =============================================== 
inline void VectorND::zero(){
  for (unsigned int i=0; i<DIMENSION; i++)  _coords[i] = 0.;
}
// ===============================================
inline double VectorND::size_sq() const{
#if DIMENSION == 1
  return (_coords[0]*_coords[0]);
#endif
#if DIMENSION == 2
  return (_coords[0]*_coords[0] +_coords[1]*_coords[1]);
#endif
#if DIMENSION == 3
  return (_coords[2]*_coords[2] +_coords[0]*_coords[0] +_coords[1]*_coords[1]);
#endif
}



 




// ===================================================== 
inline bool VectorND::operator == (const VectorND& rhs) const{
#if DIMENSION == 1
  return (_coords[0] == rhs._coords[0]);
#endif
#if DIMENSION == 2
  return (_coords[0] == rhs._coords[0] &&
	  _coords[1] == rhs._coords[1]);
#endif
#if DIMENSION == 3
  return (_coords[0] == rhs._coords[0] &&
	  _coords[1] == rhs._coords[1] &&
	  _coords[2] == rhs._coords[2]);
#endif
}
// ========================================================
inline bool VectorND::operator != (const VectorND& rhs) const{
  return (!(*this == rhs));
}





// ================================================================
//  ==================   Point =============================
// ===============================================================

class Point : public VectorND{
 public:

  /// Constructor
  Point  (const double x=0., const double y=0., const double z=0.);
  /// Copy-constructor.
  Point (const Point& p);
  /// Copy-constructor.
  Point (const VectorND& p);
  /// Empty.
  ~Point() {}

//   /**
//    * @returns a key associated with this point.  Useful for sorting.
//    */
//   dof_id_type key() const;
 protected:
  /// Make the derived class a friend
  friend class Node;
};

//------------------------------------------------------
// Inline functions
// ======================================================
inline Point::Point (const double x, const double y, const double z) :
  VectorND (x,y,z)
{}
// ====================================================
inline Point::Point (const Point& p) :  VectorND (p)
{}
// =======================================================
inline Point::Point (const VectorND& p) :  VectorND (p)
{}






#endif // TYPE_VECTOR_H
