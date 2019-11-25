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

/*
This is a modified version of libspatialindex for use with libsupermesh.
Code first added 2016-03-01.
*/

#include <spatialindex/SpatialIndex.h>

#include "../rtree/RTree.h"
#include "../mvrtree/MVRTree.h"
#include "../tprtree/TPRTree.h"

libsupermesh::SpatialIndex::InvalidPageException::InvalidPageException(id_type id)
{
	std::ostringstream s;
	s << "Unknown page id " << id;
	m_error = s.str();
}

std::string libsupermesh::SpatialIndex::InvalidPageException::what()
{
	return "InvalidPageException: " + m_error;
}

std::ostream& libsupermesh::SpatialIndex::operator<<(std::ostream& os, const ISpatialIndex& i)
{
	const libsupermesh::SpatialIndex::RTree::RTree* pRTree = dynamic_cast<const libsupermesh::SpatialIndex::RTree::RTree*>(&i);
	if (pRTree != 0)
	{
		os << *pRTree;
		return os;
	}

	const libsupermesh::SpatialIndex::MVRTree::MVRTree* pMVRTree = dynamic_cast<const libsupermesh::SpatialIndex::MVRTree::MVRTree*>(&i);
	if (pMVRTree != 0)
	{
		os << *pMVRTree;
		return os;
	}

	const libsupermesh::SpatialIndex::TPRTree::TPRTree* pTPRTree = dynamic_cast<const libsupermesh::SpatialIndex::TPRTree::TPRTree*>(&i);
	if (pTPRTree != 0)
	{
		os << *pTPRTree;
		return os;
	}

	std::cerr << "ISpatialIndex operator<<: Not implemented yet for this index type." << std::endl;
	return os;
}

std::ostream& libsupermesh::SpatialIndex::operator<<(std::ostream& os, const IStatistics& s)
{
	const libsupermesh::SpatialIndex::RTree::Statistics* pRTreeStats = dynamic_cast<const libsupermesh::SpatialIndex::RTree::Statistics*>(&s);
	if (pRTreeStats != 0)
	{
		os << *pRTreeStats;
		return os;
	}

	const libsupermesh::SpatialIndex::MVRTree::Statistics* pMVRTreeStats = dynamic_cast<const libsupermesh::SpatialIndex::MVRTree::Statistics*>(&s);
	if (pMVRTreeStats != 0)
	{
		os << * pMVRTreeStats;
		return os;
	}

	const libsupermesh::SpatialIndex::TPRTree::Statistics* pTPRTreeStats = dynamic_cast<const libsupermesh::SpatialIndex::TPRTree::Statistics*>(&s);
	if (pTPRTreeStats != 0)
	{
		os << * pTPRTreeStats;
		return os;
	}

	std::cerr << "IStatistics operator<<: Not implemented yet for this index type." << std::endl;
	return os;
}

