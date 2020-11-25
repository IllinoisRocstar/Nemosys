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
 * The following code is derived from include/Element_Intersection.h, 
 * include/MeshDataStream.h, and femtools/Element_Intersection.cpp in
 * Fluidity git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated
 * 2015-02-25). This uses modified versions of code from Rtree 0.4.3, added
 * 2016-02-24. This uses modified code from libspatialindex 1.8.5, added
 * 2016-02-24.
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

// Rtree 0.4.3 copyright information:

/*
# =============================================================================
# Rtree spatial index. Copyright (C) 2007 Sean C. Gillies
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License 
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact email: sgillies@frii.com
# =============================================================================
*/

// libspatialindex 1.8.5 copyright information:

/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include "libsupermesh_configuration.h"
#include "spatialindex/SpatialIndex.h"
#include "rtree/RTree.h"
#include "rtree/BulkLoader.h"

#include <cmath>

#include "R-Tree_Intersection_Finder_C++.h"

using namespace libsupermesh;
using namespace libsupermesh::SpatialIndex;

// Modified version of code from rtree/gispyspatialindex.h
// (GISPySpatialIndex), rtree/gispyspatialindex.cc (GISPySpatialIndex), and
// rtree/wrapper.cc (RtreeIndex_intersects) in Rtree 0.4.3. Added 2016-02-24.
// Original Fluidity ElementIntersectionFinder had comment:
  // Interface to spatialindex to calculate element intersection lists between
  // meshes using bulk storage
  // Uses code from gispatialindex.{cc,h} in Rtree 0.4.1
libsupermesh::RTree::RTree(const int &dim, const double *positions,
  const int &loc, const int &nelements, const int *enlist)
  : dim(dim), visitor(nelements) {
  this->memory = StorageManager::createNewMemoryStorageManager();
  this->buffer = StorageManager::createNewRandomEvictionsBuffer(*this->memory, capacity, bWriteThrough);
    
  // Modified version of code from createAndBulkLoadNewRTree in
  // src/rtree/RTree.cc in libspatialindex 1.8.5. Added 2016-02-24.     
  id_type indexIdentifier = 0;  // ??
  this->tree = SpatialIndex::RTree::createNewRTree(*this->buffer, fillFactor,
    indexCapacity, leafCapacity, this->dim, SpatialIndex::RTree::RV_RSTAR,
    indexIdentifier);
  uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor)));
  uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));
  MeshDataStream stream(dim, positions, loc, nelements, enlist);
  uint32_t pageSize = std::numeric_limits<uint32_t>::max(), numberOfPages = 1;  // Never cache on disk
  SpatialIndex::RTree::BulkLoader bl;
  bl.bulkLoadUsingSTR(static_cast<SpatialIndex::RTree::RTree*>(this->tree), stream, bindex, bleaf, pageSize, numberOfPages); 
  // End of modified code from createAndBulkLoadNewRTree in
  // src/rtree/RTree.cc in libspatialindex 1.8.5
}
// End of modified code from rtree/gispyspatialindex.h,
// rtree/gispyspatialindex.cc, and rtree/wrapper.cc

extern "C" {  
  void libsupermesh_build_rtree(void **rtree, int dim, int nnodes,
    const double *positions, int loc, int nelements, const int *enlist) {
    *rtree = static_cast<void*>(new libsupermesh::RTree(dim, positions, loc, nelements, enlist));
  }
  
  void libsupermesh_query_rtree(void **rtree, int dim, int loc_a,
    const double *element_a, int *neles_b) {
    assert(dim == (static_cast<libsupermesh::RTree*>(*rtree))->dim);
    *neles_b = (static_cast<libsupermesh::RTree*>(*rtree))->query(loc_a, element_a);
  }
  
  void libsupermesh_query_rtree_intersections(void **rtree, int *eles_b) {
    (static_cast<libsupermesh::RTree*>(*rtree))->query_intersections(eles_b);
  }
  
  void libsupermesh_deallocate_rtree(void **rtree) {
    delete (static_cast<libsupermesh::RTree*>(*rtree));
  }
}
