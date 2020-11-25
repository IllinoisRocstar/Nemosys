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

#ifndef LIBSUPERMESH_R_TREE_INTERSECTION_FINDER_CPP_H
#define LIBSUPERMESH_R_TREE_INTERSECTION_FINDER_CPP_H

#include "libsupermesh_configuration.h"
#include "spatialindex/SpatialIndex.h"

#include <cstring>
#include <limits>
#include <vector>

#include "libsupermesh_debug_C.h"

namespace libsupermesh {
  // Modified version of code from rtree/gispyspatialindex.h (GISPySpatialIndex)
  // and rtree/gispyspatialindex.cc (GISPySpatialIndex) in Rtree 0.4.3. Added
  // 2016-02-24.
  const uint32_t capacity = 10;
  const bool bWriteThrough = false;

  // R-Tree parameters
  const double fillFactor = 0.7;
  const uint32_t indexCapacity = 10;
  const uint32_t leafCapacity = 10;
  // End of modified code from rtree/gispyspatialindex.h and
  // rtree/gispyspatialindex.cc in Rtree 0.4.3
  
  // Modified version of MyDataStream class from test/rtree/RTreeBulkLoad.cc in
  // libspatialindex 1.8.5. Added 2016-02-24. Original Fluidity code had
  // comment:
    // Customised version of MyDataStream class in
    // regressiontest/rtree/RTreeBulkLoad.cc in spatialindex 1.2.0
  class MeshDataStream : public SpatialIndex::IDataStream {
    public:
      inline MeshDataStream(const int &dim, const double *positions,
        const int &loc, const int &nelements, const int *enlist)
     : dim(dim), loc(loc), nelements(nelements), positions(positions),
       enlist(enlist), ele(0) {};
      inline virtual ~MeshDataStream(void) {}
      inline virtual SpatialIndex::IData* getNext(void) {
        if(this->ele >= this->nelements) return NULL;

        // Indexing nodes from one
        double low[this->dim], high[this->dim];
        for(uint32_t i = 0;i < this->dim;i++) {
          low[i] = this->positions[this->dim * (this->enlist[this->ele * this->loc] - 1) + i];
          high[i] = low[i];
        }
        for(uint32_t j = 1;j < this->loc;j++) {
          for(uint32_t i = 0;i < this->dim;i++) {
            low [i] = std::min(low [i], this->positions[this->dim * (this->enlist[this->ele * this->loc + j] - 1) + i]);
            high[i] = std::max(high[i], this->positions[this->dim * (this->enlist[this->ele * this->loc + j] - 1) + i]);
          }
        }        
        SpatialIndex::Region bbox(low, high, this->dim);
        SpatialIndex::id_type id = this->ele + 1;  // Indexing elements from one
                                                                           // No data
        SpatialIndex::RTree::Data *m_pNext = new SpatialIndex::RTree::Data((uint32_t)0, NULL,
                                                                           bbox, id);
        
        this->ele++;
        return m_pNext;
      }
      inline virtual bool hasNext(void) {return this->ele < this->nelements;}
      inline virtual uint32_t size(void) {return this->nelements;}
      inline virtual void rewind(void) {this->ele = 0;}
      
    private:
      const uint32_t dim, loc, nelements;
      const double *positions;
      const int *enlist;
      uint32_t ele;
  };
  // End of modified code from test/rtree/RTreeBulkLoad.cc in libspatialindex
  // 1.8.5
  
  // Modified version of PyListVisitor class from rtree/wrapper.cc in Rtree
  // 0.4.3. Added 2016-02-24. Original Fluidity code had comment:
    // Customised version of PyListVisitor class in
    // wrapper.cc in Rtree 0.4.1
  class ElementListVisitor : public SpatialIndex::IVisitor {
    public:
        inline ElementListVisitor(const int &nelements) {
          this->eles = new int[nelements];
          this->neles = 0;
        }
        inline virtual ~ElementListVisitor(void) {delete[] this->eles;}
        inline virtual void visitNode(const SpatialIndex::INode &n) {}
        inline virtual void visitData(const SpatialIndex::IData &d) {this->eles[this->neles++] = d.getIdentifier();}
        inline virtual void visitData(std::vector<const SpatialIndex::IData*> &d) {
          libsupermesh_abort("void ElementListVisitor::visitData(std::vector<const SpatialIndex::IData*>& not implemented");
        }
        inline virtual void clear(void) {this->neles = 0;}
        inline virtual const int size(void) const {return this->neles;}
        inline virtual void data(int *eles) const {memcpy(eles, this->eles, this->neles * sizeof(int));}
      
      private:
        int *eles, neles;
  };  
  // End of modified code from rtree/wrapper.cc in Rtree 0.4.3

  // Modified version of code from rtree/gispyspatialindex.h
  // (GISPySpatialIndex), rtree/gispyspatialindex.cc (GISPySpatialIndex), and
  // rtree/wrapper.cc (RtreeIndex_intersects) in Rtree 0.4.3. Added 2016-02-24.
  // Original Fluidity ElementIntersectionFinder had comment:
    // Interface to spatialindex to calculate element intersection lists between
    // meshes using bulk storage
    // Uses code from gispatialindex.{cc,h} in Rtree 0.4.1
  class RTree {
    public:
      RTree(const int &dim, const double *positions, const int &loc,
        const int &nelements, const int *enlist);
      
      inline ~RTree(void) {
        delete this->tree;
        delete this->buffer;
        delete this->memory;
      }      
      inline const uint32_t query(const int &loc_a, const double *element_a) {
        // See MeshDataStream::getNext above
        double low[this->dim], high[this->dim];
        for(uint32_t i = 0;i < this->dim;i++) {
          low[i] = element_a[i];
          high[i] = low[i];
        }
        for(int j = 1;j < loc_a;j++) {
          for(uint32_t i = 0;i < this->dim;i++) {
            low [i] = std::min(low [i], element_a[j * this->dim + i]);
            high[i] = std::max(high[i], element_a[j * this->dim + i]);
          }
        }
        SpatialIndex::Region bbox(low, high, this->dim);
        this->visitor.clear();
        this->tree->intersectsWithQuery(bbox, this->visitor);
        
        return this->visitor.size();
      }         
      inline void query_intersections(int *eles_b) const {this->visitor.data(eles_b);}
    
      const uint32_t dim;
    private:
      SpatialIndex::IStorageManager *memory;
      SpatialIndex::StorageManager::IBuffer *buffer;
      SpatialIndex::ISpatialIndex *tree;
      ElementListVisitor visitor;
  };
  // End of modified code from rtree/gispyspatialindex.h,
  // rtree/gispyspatialindex.cc, and rtree/wrapper.cc
}

#endif
