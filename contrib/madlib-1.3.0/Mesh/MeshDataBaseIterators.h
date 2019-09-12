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

#ifndef H_MESHDATABASEITERATORS
#define H_MESHDATABASEITERATORS

#include <stdio.h>

namespace MAd {

  template <class CONTAINER, class DATA, class CLASSIF>
  class MDB_Iter
  {
  public:
    CONTAINER *l;
    CLASSIF   g;
    int closure;
  
    typename CONTAINER::iterator current;
    typename CONTAINER::iterator end;
  
    MDB_Iter (CONTAINER *_l):l(_l),g(NULL),closure(0)
    {
      reset();
    }
    MDB_Iter (CONTAINER *_l,CLASSIF _g):l(_l),g((CLASSIF) _g),closure(0) {
      reset();
    }
  
    MDB_Iter (CONTAINER *_l,CLASSIF _g,int _c):l(_l),g((CLASSIF) _g),closure(_c) {

      /*     if ((DATA::getDim() == 3) && (_c == 1)) { */
      /*       printf("Closure iterator not allowed for regions\n"); */
      /*       exit(1); */
      /*     } */
    
      if (_c == 1) {
        printf("Closure iterator not yet implemented \n");
        exit(1);
      }
    
      reset();
    }
    
    virtual ~MDB_Iter() {}
  
    inline void reset()
    {
      end = l->end();
      current = l->begin();

      if (!g) {
        while(current != end && (*current)->deleted)
          {
            typename CONTAINER::iterator currentb = current;
            current++;
            delete *currentb;
            l->erase(currentb);       
          }
      }
      else {  
        while(current != end && ((*current)->deleted || (*current)->g != g))
          {
            typename CONTAINER::iterator currentb = current;
            current++;
            if ((*currentb)->deleted) {
              delete *currentb;
              l->erase(currentb);
            }
          }
      }
      //     if (!g  && current!= end && (*current)->deleted)next();
      //     if (g   && current!= end && ((*current)->deleted || (*current)->g != g)) next();
    
    }
    inline void cleanup()
    {
      current = l->begin();
      while(next()){}  
    }
  
    inline virtual DATA * next ()
    {
      if (current == end) return 0;
    
      DATA *r = *current;
      ++ current;

      if (!g) {
        while(current != end && (*current)->deleted)
          {
            typename CONTAINER::iterator currentb = current;
            current++;
            delete *currentb;
            l->erase(currentb);       
          }
      }
      else {  
        while(current != end && ((*current)->deleted || (*current)->g != g))
          {
            typename CONTAINER::iterator currentb = current;
            current++;
            if ((*currentb)->deleted) {
              delete *currentb;
              l->erase(currentb);
            }
          }
      }

      // if the result data was marked for deletion
      // after previous call to next() look for the next one
      if (r->deleted) {

        if (current == end) return 0;

        r = *current;
        ++ current;
      
        if (!g) {
          while(current != end && (*current)->deleted)
            {
              typename CONTAINER::iterator currentb = current;
              current++;
              delete *currentb;
              l->erase(currentb);       
            }
        }
        else {  
          while(current != end && ((*current)->deleted || (*current)->g != g))
            {
              typename CONTAINER::iterator currentb = current;
              current++;
              if ((*currentb)->deleted) {
                delete *currentb;
                l->erase(currentb);
              }
            }
        }

      }

      return r;
    }
  
  };

  template <class CONTAINER, class DATA, class CLASSIF>
  int countClassifiedElements (CONTAINER* c,CLASSIF g) {
  
    MDB_Iter<CONTAINER,DATA,CLASSIF>* iter = new MDB_Iter<CONTAINER,DATA,CLASSIF>(c,g);
  
    DATA* el = NULL;
    int count = 0;
    while ( (el = iter->next()) ) count++;
    delete iter;

    return count;
  
  }

}

  
/* //report errors to koen.hillewaert@cenaero.be */
/* template <class CONTAINER, class DATA> */
/* class MDB_CIter: public  MDB_Iter <CONTAINER,DATA> */
/* { */
/* public: */
/*   pGEntity ge; */
/*   MDB_CIter (CONTAINER *_l,pGEntity _g):MDB_Iter<CONTAINER,DATA>(_l),ge(_g) */
/*   { */
/*     reset(); */
/*   } */
/*   inline DATA * next () */
/*   { */
/*     if (current == end)return 0; */
/*     DATA *r = *current; */
/*     ++ current; */
/*     while(current != end && */
/*           (*current)->deleted && */
/*           EN_whatIn(current) != ge) */
/*       { */
/*         typename CONTAINER::iterator currentb = current; */
/*         current++; */
/*         delete *currentb; */
/*         l->erase(currentb);        */
/*       } */
/*     return r; */
/*   } */
/* }; */

#endif
