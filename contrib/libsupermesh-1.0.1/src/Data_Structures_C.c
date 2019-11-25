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
 
/*
 * Fluidity copyright information (note that AUTHORS mentioned in the following
 * has been renamed to Fluidity_AUTHORS):
 */

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

#ifdef LIBSUPERMESH_ENABLE_JUDY

#include <Judy.h>
#include <stdio.h>

#include "libsupermesh_debug_C.h"

/* To understand these, read
   http://judy.sourceforge.net/doc/Judy1_3x.htm and
   http://judy.sourceforge.net/doc/JudyL_3x.htm */

void libsupermesh_integer_set_new(Pvoid_t *i) {
  *i = NULL;
}

void libsupermesh_integer_set_delete(Pvoid_t *i) {
  Word_t mem_freed;
  J1FA(mem_freed, *i);
}

void libsupermesh_integer_set_insert(Pvoid_t *i, int value, int *changed) {
  Word_t lvalue = value;
  J1S(*changed, *i, lvalue);
}

void libsupermesh_integer_set_size(Pvoid_t *i, int *size) {
  Word_t lsize;
  J1C(lsize, *i, 0, -1);
  *size = lsize;
}

void libsupermesh_integer_set_fetch(Pvoid_t *i, int index, int *value) {
  Word_t lindex = index, lvalue;
  int worked;
  J1BC(worked, *i, lindex, lvalue);
  if(!worked) {
    fprintf(stderr, "Failed to fetch integer set element with index %i\n", index);
    libsupermesh_abort("Failed to fetch integer set element");
  }
  *value = lvalue;
}

void libsupermesh_integer_set_remove(Pvoid_t *i, int value) {
  Word_t lvalue = value;
  int worked;
  J1U(worked, *i, lvalue); 
  if(!worked) {
    fprintf(stderr, "Failed to remove integer set element with value %i\n", value);
    libsupermesh_abort("Failed to remove integer set element");
  }
}

void libsupermesh_integer_set_has_value(Pvoid_t *i, int value, int *present) {
  Word_t lvalue = value;
  J1T(*present, *i, lvalue); 
}

void libsupermesh_integer_map_new(Pvoid_t *i) {
  *i = NULL;
}

void libsupermesh_integer_map_delete(Pvoid_t *i) {
  Word_t mem_freed;
  JLFA(mem_freed, *i);
}

void libsupermesh_integer_map_insert(Pvoid_t *i, int key, int value) {
  Word_t lkey = key;
  PWord_t lvalue;
  JLI(lvalue, *i, lkey);
  *lvalue = value;
}

void libsupermesh_integer_map_size(Pvoid_t *i, int *size) {
  Word_t lsize;
  JLC(lsize, *i, 0, -1);
  *size = lsize;
}

void libsupermesh_integer_map_fetch(Pvoid_t *i, int key, int *value) {
  Word_t lkey = key;
  PWord_t lvalue;
  JLG(lvalue, *i, lkey); 
  if(!lvalue) {
    fprintf(stderr, "Failed to fetch integer map element with key %i\n", key);
    libsupermesh_abort("Failed to fetch integer map element");
  }
  *value = *lvalue;
}

void libsupermesh_integer_map_remove(Pvoid_t *i, int key) {
  Word_t lkey = key;
  int worked;
  JLD(worked, *i, lkey); 
  if(!worked) {
    fprintf(stderr, "Failed to remove integer map element with key %i\n", key);
    libsupermesh_abort("Failed to remove integer map element");
  }
}

void libsupermesh_integer_map_has_key(Pvoid_t *i, int key, int *present) {
  Word_t lkey = key;
  PWord_t lvalue;
  JLG(lvalue, *i, lkey); 
  *present = (lvalue != NULL);
}

void libsupermesh_integer_map_fetch_pair(Pvoid_t *i, int index, int *key, int *value) {
  Word_t lindex = index;
  Word_t lkey;
  PWord_t lvalue;
  JLBC(lvalue, *i, lindex, lkey); 
  if(!lvalue) {
    fprintf(stderr, "Failed to fetch integer map element with index %i\n", index);
    libsupermesh_abort("Failed to fetch integer map element");
  }
  *key = lkey;
  *value = *lvalue;
}

#endif
