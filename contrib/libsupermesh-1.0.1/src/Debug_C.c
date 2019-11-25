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
 * The following code is derived from debug/C++_Debug.cpp in Fluidity git
 * revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
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

#include "libsupermesh_configuration.h"

#ifdef Backtrace_FOUND
#include Backtrace_HEADER
#include <stdio.h>
#include <stdlib.h>

#include "libsupermesh_debug_C.h"

void libsupermesh_print_backtrace(int max_size) {
  char **symbols = NULL;
  int i, size;
  void **buffer = NULL;
  
  buffer = (void**)malloc(max_size * sizeof(void*));
  if(!buffer) goto cleanup;
  
  size = backtrace(buffer, max_size);
  symbols = backtrace_symbols(buffer, size);
  if(!symbols) goto cleanup;
  
  for(i = 0;i < size;i++) {
    fprintf(stderr, "%s\n", symbols[i]);
  }
  
cleanup:
  if(buffer) free(buffer);
  if(symbols) free(symbols);
}
#else
void libsupermesh_print_backtrace(int max_size) {}
#endif
