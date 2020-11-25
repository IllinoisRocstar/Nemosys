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

#ifndef LIBSUPERMESH_DATA_STRUCTURES_C_H
#define LIBSUPERMESH_DATA_STRUCTURES_C_H

#include "libsupermesh_configuration.h"

#ifdef __cplusplus
extern "C" {
#endif
  void libsupermesh_integer_set_new(void **i);
  void libsupermesh_integer_set_delete(void **i);
  void libsupermesh_integer_set_insert(void **i, int value, int *changed);
  void libsupermesh_integer_set_size(void **i, int *size);
  void libsupermesh_integer_set_fetch(void **i, int index, int *value);
  void libsupermesh_integer_set_remove(void **i, int value);
  void libsupermesh_integer_set_has_value(void **i, int value, int *present);
  
  void libsupermesh_integer_map_new(void **i);
  void libsupermesh_integer_map_delete(void **i);
  void libsupermesh_integer_map_insert(void **i, int key, int value);
  void libsupermesh_integer_map_size(void **i, int *size);
  void libsupermesh_integer_map_fetch(void **i, int key, int *value);
  void libsupermesh_integer_map_remove(void **i, int key);
  void libsupermesh_integer_map_has_key(void **i, int key, int *present);
  void libsupermesh_integer_map_fetch_pair(void **i, int index, int *key, int *value);
#ifdef __cplusplus
}
#endif

#endif
