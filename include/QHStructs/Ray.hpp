#ifndef QuickHull_Ray_hpp
#define QuickHull_Ray_hpp

#include "nemosys_export.h"

#include "Vector3.hpp"

namespace NEM {

namespace GEO {

namespace quickhull {

template <typename T>
struct NEMOSYS_EXPORT  Ray
{
  const Vector3<T> m_S;
  const Vector3<T> m_V;
  const T m_VInvLengthSquared;
  
  Ray(const Vector3<T>& S,const Vector3<T>& V) : m_S(S), m_V(V), 
    m_VInvLengthSquared(1/m_V.getLengthSquared()) 
    {}
};
	
} // namespace quickhull

} // namespace GEO

} // namespace NEM

#endif
