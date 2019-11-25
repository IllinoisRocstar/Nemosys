/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2004, Marios Hadjieleftheriou
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

#pragma once

namespace libsupermesh { namespace SpatialIndex
{
	namespace MVRTree
	{
		SIDX_DLL enum MVRTreeVariant
		{
			RV_LINEAR = 0x0,
			RV_QUADRATIC,
			RV_RSTAR
		};

		SIDX_DLL enum PersistenObjectIdentifier
		{
			PersistentIndex = 0x1,
			PersistentLeaf = 0x2
		};

		SIDX_DLL enum RangeQueryType
		{
			ContainmentQuery = 0x1,
			IntersectionQuery = 0x2
		};

		class SIDX_DLL Data : public IData, public libsupermesh::Tools::ISerializable
		{
		public:
			Data(uint32_t len, byte* pData, TimeRegion& r, id_type id);
			virtual ~Data();

			virtual Data* clone();
			virtual id_type getIdentifier() const;
			virtual void getShape(IShape** out) const;
			virtual void getData(uint32_t& len, byte** data) const;
			virtual uint32_t getByteArraySize();
			virtual void loadFromByteArray(const byte* data);
			virtual void storeToByteArray(byte** data, uint32_t& len);

			id_type m_id;
			TimeRegion m_region;
			byte* m_pData;
			uint32_t m_dataLength;
		}; // Data

		SIDX_DLL ISpatialIndex* returnMVRTree(IStorageManager& ind, libsupermesh::Tools::PropertySet& in);
		SIDX_DLL ISpatialIndex* createNewMVRTree(
			IStorageManager& in,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			MVRTreeVariant rv,
			id_type& out_indexIdentifier
		);
		SIDX_DLL ISpatialIndex* loadMVRTree(
			IStorageManager& in,
			id_type indexIdentifier
		);
	}
} }

