#ifndef VertexDataSource_h
#define VertexDataSource_h

#include "nemosys_export.h"

#include "Vector3.hpp"

namespace NEM {

namespace GEO {

namespace quickhull {
	
template<typename T>
class NEMOSYS_EXPORT VertexDataSource 
{
  const Vector3<T>* m_ptr;
  size_t m_count;
  
  public:
  VertexDataSource(const Vector3<T>* ptr, size_t count) : m_ptr(ptr), m_count(count) 
    {}
  
  VertexDataSource(const std::vector<Vector3<T>>& vec) : m_ptr(&vec[0]), m_count(vec.size()) 
    {}
  
  VertexDataSource() : m_ptr(nullptr), m_count(0) {}
  
  VertexDataSource& operator=(const VertexDataSource& other) = default;
  
  size_t size() const 
  {
  	return m_count;
  }
  
  const Vector3<T>& operator[](size_t index) const 
  {
  	return m_ptr[index];
  }
  
  const Vector3<T>* begin() const 
  {
  	return m_ptr;
  }
  
  const Vector3<T>* end() const 
  {
  	return m_ptr + m_count;
  }
};
	
} // namespace quickhull

} // namespace GEO

} // namespace NEM


#endif /* VertexDataSource_h */
