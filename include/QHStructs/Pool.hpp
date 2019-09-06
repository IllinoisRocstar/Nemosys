#ifndef Pool_h
#define Pool_h

#include "nemosys_export.h"

#include <vector>
#include <memory>

namespace NEM {

namespace GEO {

namespace quickhull {
	
template<typename T>
class NEMOSYS_EXPORT Pool 
{
  public:
  void clear() 
  {
  	m_data.clear();
  }
  
  void reclaim(std::shared_ptr<T>& ptr) 
  {
  	m_data.push_back(ptr);
  }
  
  std::shared_ptr<T> get() 
  {
  	if (m_data.size()==0) 
    {
  		return std::unique_ptr<T>(new T());
  	}
  	auto it = m_data.end()-1;
  	std::shared_ptr<T> r = std::move(*it);
  	m_data.erase(it);
  	return r;
  }

  private:
  std::vector<std::shared_ptr<T>> m_data;

};
	
} // namespace quickhull

} // namespace GEO

} // namespace NEM


#endif /* Pool_h */
