

// C++ includes
#include <iostream>
#include <iomanip> // for std::setw, std::setiosflags

// Local includes
#include "type_vector.h"

// ------------------------------------------------------------
// VectorND class member funcions
VectorND VectorND::unit() const{

  const double length = size();
  assert(length != static_cast<double>(0.));
#if DIMENSION == 1
  return VectorND(_coords[0]/length);
#endif
#if DIMENSION == 2
  return VectorND(_coords[0]/length,
		       _coords[1]/length);
#endif
#if DIMENSION == 3
  return VectorND(_coords[0]/length,
		       _coords[1]/length,
		       _coords[2]/length);
#endif

}
// ================================================
void VectorND::print(std::ostream& os) const
{
#if DIMENSION == 1

  os << "x=" << (*this)(0);

#endif
#if DIMENSION == 2

  os << "(x,y)=("
     << std::setw(8) << (*this)(0) << ", "
     << std::setw(8) << (*this)(1) << ")";

#endif
#if DIMENSION == 3

  os <<  "(x,y,z)=("
     << std::setw(8) << (*this)(0) << ", "
     << std::setw(8) << (*this)(1) << ", "
     << std::setw(8) << (*this)(2) << ")";
#endif
}
// ===================================================
void VectorND::write_unformatted (
  std::ostream &os,const bool newline) const{
  assert(os);
  os << std::setiosflags(std::ios::showpoint) << (*this)(0) 
       << " " << (*this)(1) << " "<< (*this)(2) << " ";
  if (newline)  os << '\n';
  }
// ======================================================
bool VectorND::operator < (const VectorND& rhs) const{
  for (unsigned int i=0; i<DIMENSION; i++)    {
      if ((*this)(i) < rhs(i))   return true;
      if ((*this)(i) > rhs(i))  return false;
    }
  return false;
}
// ====================================================== 
bool VectorND::operator <= (const VectorND& rhs) const{
  for (unsigned int i=0; i<DIMENSION; i++)   {
      if ((*this)(i) < rhs(i))        return true;
      if ((*this)(i) > rhs(i))        return false;
    }
  return true;
}
// =====================================================
bool VectorND::operator > (const VectorND& rhs) const{
  for (unsigned int i=0; i<DIMENSION; i++)    {
      if ((*this)(i) > rhs(i))        return true;
      if ((*this)(i) < rhs(i))        return false;
    }
  return false;
}
// ====================================================
bool VectorND::operator >= (const VectorND& rhs) const{
  for (unsigned int i=0; i<DIMENSION; i++)    {
      if ((*this)(i) > rhs(i))        return true;
      if ((*this)(i) < rhs(i))        return false;
    }
  return true;
}




