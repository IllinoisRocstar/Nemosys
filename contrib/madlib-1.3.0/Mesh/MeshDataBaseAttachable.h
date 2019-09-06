// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _MATTACHABLEDATACONTAINER_H_
#define _MATTACHABLEDATACONTAINER_H_

#include <algorithm>
#include <vector>

namespace MAd {

  /**
     Base class for attachable data's
  */
  class mAttachableData
  {
  public :
    virtual ~mAttachableData(){};
  };

  class mAttachableVoid : public mAttachableData
  {
  public :
    ~mAttachableVoid(){}
    void *veryNastyPointer;
  };

  /**
     Int as attachable data
  */
  class mAttachableInt : public mAttachableData
  {
  public :
    ~mAttachableInt(){}
    int i;
  };

  /**
     Double as attachable data
  */
  class mAttachableDouble : public mAttachableData
  {
  public :  
    ~mAttachableDouble(){}
    double d;
  };

  /**
     Container for attachable data's Internal, mEntity and mMesh provides interfaces.
  */

  class mAttachableDataContainer  
  {
  public :
    typedef std::pair <unsigned int, mAttachableData *> info;
    typedef std::vector< info > container;
    typedef container::iterator iter_attachdata;
    typedef container::const_iterator citer_attachdata;
  private :
    container *tab;
  public:
    mAttachableDataContainer():tab(0){}
    ~mAttachableDataContainer();
    inline void attachData(unsigned int, mAttachableData *);
    inline void deleteData(unsigned int);
    inline mAttachableData * getData(unsigned int) ;
    inline citer_attachdata begin_attachdata() const {return tab->begin();};
    inline citer_attachdata end_attachdata() const {return tab->end();};	
    /// specific data types for int
    inline void attachInt(unsigned int, int);
    /// specific data types for int
    inline int getAttachedInt(unsigned int);
    /// specific data types for double
    inline void attachDouble(unsigned int, double);
    /// specific data types for double
    inline double getAttachedDouble(unsigned int);
  };

  class equalInfoPred
  {
    const unsigned int c;
  public:
    equalInfoPred ( unsigned int i ) :c(i) {}
    inline bool operator () (const mAttachableDataContainer::info &i) const
    {
      return (i.first == c);
    }
  };

  inline mAttachableDataContainer::~mAttachableDataContainer()
  {
    if(!tab)return;
    citer_attachdata it = begin_attachdata();
    citer_attachdata itEnd = end_attachdata();
    for(;it!=itEnd;++it)
      {
        mAttachableData *a = (*it).second;
        delete a;
      }
    delete tab;
  }

  inline mAttachableData *mAttachableDataContainer::getData(unsigned int c) 
  {
    if (!tab)return 0;
    iter_attachdata it = std::find_if (tab->begin(),tab->end(),equalInfoPred(c));
    if(it == tab->end())return 0;
    return (*it).second;
  }

  inline void mAttachableDataContainer::attachData(unsigned int c, mAttachableData *v)
  {
    if (!tab) tab = new container;
    tab->push_back(mAttachableDataContainer::info(c,v));
  }

  inline void mAttachableDataContainer::deleteData(unsigned int c)
  {
    mAttachableData *data = getData (c);
    if (data)
      {
        delete data;
        tab->erase ( std::remove_if (tab->begin(),tab->end(),equalInfoPred(c)) , 
                     tab->end () );
      }
  }

  inline void mAttachableDataContainer::attachInt(unsigned int c, int i)
  {
    mAttachableInt *ai = (mAttachableInt *)getData(c);
    if(!ai)
      {
        ai = new mAttachableInt;
        attachData(c,ai);
      }
    ai->i = i;
  }

  inline int mAttachableDataContainer::getAttachedInt (unsigned int c)
  {
    mAttachableInt *ai = (mAttachableInt *)getData(c);
    if(!ai)return 0;
    return ai->i;
  }

  inline void mAttachableDataContainer::attachDouble(unsigned int c, double d)
  {
    mAttachableDouble *ai = (mAttachableDouble *)getData(c);
    if(!ai)
      {
        ai = new mAttachableDouble;
        attachData(c,ai);
      }
    ai->d = d;
  }

  inline double mAttachableDataContainer::getAttachedDouble (unsigned int c)
  {
    mAttachableDouble *ai = (mAttachableDouble *)getData(c);
    if(!ai)return 0;
    return ai->d;
  }

} // End of namespace MAd

#endif
