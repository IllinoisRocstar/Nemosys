/*
  Copyright (C) 2016-2017 The University of Edinburgh

  The file is part of libsupermesh
    
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation;
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 * The following code is derived from femtools/Data_structures_C.c in Fluidity
 * git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
 */
 
// Fluidity copyright information (note that AUTHORS mentioned in the following
// has been renamed to Fluidity_AUTHORS):

/*
  Copyright (C) 2006 Imperial College London and others.
  
  Please see the AUTHORS file in the main source directory for a full list
  of copyright holders.

  Prof. C Pain
  Applied Modelling and Computation Group
  Department of Earth Science and Engineering
  Imperial College London

  amcgsoftware@imperial.ac.uk
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation,
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  USA
*/

#include "Data_Structures_C.h"

#ifndef LIBSUPERMESH_ENABLE_JUDY

#include <iostream>
#include <iterator>
#include <map>
#include <set>

#include "libsupermesh_debug_C.h"

using namespace std;

namespace libsupermesh {

class integer_set {
  public:
    inline integer_set(void) : index(0) {};    
    inline ~integer_set(void) {}
    
    inline void insert(const int &value, int &changed) {
      if(this->index != 0) {
        this->index = 0;
      }
      pair<set<int>::iterator, bool> position = this->value_set.insert(value);
      changed = position.second ? 1 : 0;
    }
    
    inline void size(int &size) const {
      size = this->value_set.size();
    };
    
    inline void fetch(const int &index, int &value) {
      if(index <= 0 || index > this->value_set.size()) {
        std::cerr << "Failed to fetch integer set element with index " << index << endl;
        libsupermesh_abort("Failed to fetch integer set element");
      }
      if(this->index != index) {
        this->iter = this->value_set.begin();
        advance(this->iter, index - 1);
        this->index = index;
      }
      value = *this->iter;
      this->iter++;
      this->index++;
    };
    
    inline void remove(const int &value) {
      if(this->index != 0) {
        this->index = 0;
      }
      set<int>::size_type count = this->value_set.erase(value);
      if(count == 0) {
        cerr << "Failed to remove integer set element with value " << value << endl;
        libsupermesh_abort("Failed to remove integer set element");
      }
    }
    
    inline void has_value(const int &value, int &present) const {
      present = (this->value_set.find(value) == this->value_set.end()) ? 0 : 1;
    }
    
  private:
    set<int> value_set;
    set<int>::iterator iter;
    set<int>::size_type index;
};

class integer_map {
  public:
    inline integer_map(void) : index(0) {}    
    inline ~integer_map(void) {}
    
    inline void insert(const int &key, const int &value) {
      if(this->index != 0) {
        this->index = 0;
      }
      this->value_map[key] = value;
    }
    
    inline void size(int &size) const {
      size = this->value_map.size();
    }
    
    inline void fetch(const int &key, int &value) {
      if(this->value_map.count(key) == 0) {
        cerr << "Failed to fetch integer map element with key " << key << endl;
        libsupermesh_abort("Failed to fetch integer map element");
      }
      value = this->value_map[key];
    }
    
    inline void remove(const int &key) {
      if(this->index != 0) {
        this->index = 0;
      }
      map<int, int>::size_type count = this->value_map.erase(key);
      if(count == 0) {
        cerr << "Failed to remove integer map element with key " << key << endl;
        libsupermesh_abort("Failed to remove integer map element");
      }
    }
    
    inline void has_key(const int &key, int &present) const {
      present =  (this->value_map.find(key) == this->value_map.end()) ? 0 : 1;
    }
    
    inline void fetch_pair(const int &index, int &key, int &value) {
      if(index <= 0 || index > this->value_map.size()) {
        cerr << "Failed to fetch integer map element with index " << index << endl;
        libsupermesh_abort("Failed to fetch integer map element");
      }
      if(this->index != index) {
        this->iter = this->value_map.begin();
        advance(this->iter, index - 1);
        this->index = index;
      }
      key = this->iter->first;
      value = this->iter->second;
      this->iter++;
      this->index++;
    }
    
  private:
    map<int, int> value_map;
    map<int, int>::iterator iter;
    map<int, int>::size_type index;
};

}

extern "C" {
  void libsupermesh_integer_set_new(void **i) {
    *i = static_cast<void*>(new libsupermesh::integer_set());
  }
  
  void libsupermesh_integer_set_delete(void **i) {
    delete (static_cast<libsupermesh::integer_set*>(*i));
  }
  
  void libsupermesh_integer_set_insert(void **i, int value, int *changed) {
    (static_cast<libsupermesh::integer_set*>(*i))->insert(value, *changed);
  }
  
  void libsupermesh_integer_set_size(void **i, int *size) {
    (static_cast<libsupermesh::integer_set*>(*i))->size(*size);
  }
  
  void libsupermesh_integer_set_fetch(void **i, int index, int *value) {
    (static_cast<libsupermesh::integer_set*>(*i))->fetch(index, *value);
  }
  
  void libsupermesh_integer_set_remove(void **i, int value) {
    (static_cast<libsupermesh::integer_set*>(*i))->remove(value);
  }

  void libsupermesh_integer_set_has_value(void **i, int value, int *present) {
    (static_cast<libsupermesh::integer_set*>(*i))->has_value(value, *present);
  }
  
  void libsupermesh_integer_map_new(void **i) {
    *i = static_cast<void*>(new libsupermesh::integer_map());
  }
  
  void libsupermesh_integer_map_delete(void **i) {
    delete (static_cast<libsupermesh::integer_map*>(*i));
  }
  
  void libsupermesh_integer_map_insert(void **i, int key, int value) {
    (static_cast<libsupermesh::integer_map*>(*i))->insert(key, value);
  }
  
  void libsupermesh_integer_map_size(void **i, int *size) {    
    (static_cast<libsupermesh::integer_map*>(*i))->size(*size);
  }
  
  void libsupermesh_integer_map_fetch(void **i, int key, int *value) {
    (static_cast<libsupermesh::integer_map*>(*i))->fetch(key, *value);
  }
  
  void libsupermesh_integer_map_remove(void **i, int key) {
    (static_cast<libsupermesh::integer_map*>(*i))->remove(key);
  }
  
  void libsupermesh_integer_map_has_key(void **i, int key, int *present) {
    (static_cast<libsupermesh::integer_map*>(*i))->has_key(key, *present);
  }

  void libsupermesh_integer_map_fetch_pair(void **i, int index, int *key, int *value) {
    (static_cast<libsupermesh::integer_map*>(*i))->fetch_pair(index, *key, *value);
  }
}

#endif
