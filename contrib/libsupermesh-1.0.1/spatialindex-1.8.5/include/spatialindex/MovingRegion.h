/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2003, Marios Hadjieleftheriou
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
	class SIDX_DLL MovingRegion : public TimeRegion, public IEvolvingShape
	{
        using Region::getLow;
        using Region::getHigh;
        using TimeRegion::intersectsRegionInTime;
        using TimeRegion::containsRegionInTime;
        using TimeRegion::combineRegionInTime;
        using TimeRegion::getCombinedRegionInTime;
        using TimeRegion::containsPointInTime;
        
	public:
		MovingRegion();
		MovingRegion(
			const double* pLow, const double* pHigh,
			const double* pVLow, const double* pVHigh,
			const libsupermesh::Tools::IInterval& ti, uint32_t dimension);
		MovingRegion(
			const double* pLow, const double* pHigh,
			const double* pVLow, const double* pVHigh,
			double tStart, double tEnd, uint32_t dimension);
		MovingRegion(
			const Point& low, const Point& high,
			const Point& vlow, const Point& vhigh,
			const libsupermesh::Tools::IInterval& ti);
		MovingRegion(
			const Point& low, const Point& high,
			const Point& vlow, const Point& vhigh,
			double tStart, double tEnd);
		MovingRegion(const Region& mbr, const Region& vbr, const libsupermesh::Tools::IInterval& ivI);
		MovingRegion(const Region& mbr, const Region& vbr, double tStart, double tEnd);
		MovingRegion(const MovingPoint& low, const MovingPoint& high);
		MovingRegion(const MovingRegion& in);
		virtual ~MovingRegion();

		virtual MovingRegion& operator=(const MovingRegion& r);
		virtual bool operator==(const MovingRegion&) const;

		bool isShrinking() const;

		virtual double getLow(uint32_t index, double t) const;
		virtual double getHigh(uint32_t index, double t) const;
		virtual double getExtrapolatedLow(uint32_t index, double t) const;
		virtual double getExtrapolatedHigh(uint32_t index, double t) const;
		virtual double getVLow(uint32_t index) const;
		virtual double getVHigh(uint32_t index) const;

		virtual bool intersectsRegionInTime(const MovingRegion& r) const;
		virtual bool intersectsRegionInTime(const MovingRegion& r, libsupermesh::Tools::IInterval& out) const;
		virtual bool intersectsRegionInTime(const libsupermesh::Tools::IInterval& ivI, const MovingRegion& r, libsupermesh::Tools::IInterval& ret) const;
		virtual bool containsRegionInTime(const MovingRegion& r) const;
		virtual bool containsRegionInTime(const libsupermesh::Tools::IInterval& ivI, const MovingRegion& r) const;
		virtual bool containsRegionAfterTime(double t, const MovingRegion& r) const;

		virtual double getProjectedSurfaceAreaInTime() const;
		virtual double getProjectedSurfaceAreaInTime(const libsupermesh::Tools::IInterval& ivI) const;

		virtual double getCenterDistanceInTime(const MovingRegion& r) const;
		virtual double getCenterDistanceInTime(const libsupermesh::Tools::IInterval& ivI, const MovingRegion& r) const;

		virtual bool intersectsRegionAtTime(double t, const MovingRegion& r) const;
		virtual bool containsRegionAtTime(double t, const MovingRegion& r) const;

		virtual bool intersectsPointInTime(const MovingPoint& p) const;
		virtual bool intersectsPointInTime(const MovingPoint& p, libsupermesh::Tools::IInterval& out) const;
		virtual bool intersectsPointInTime(const libsupermesh::Tools::IInterval& ivI, const MovingPoint& p, libsupermesh::Tools::IInterval& out) const;
		virtual bool containsPointInTime(const MovingPoint& p) const;
		virtual bool containsPointInTime(const libsupermesh::Tools::IInterval& ivI, const MovingPoint& p) const;

		//virtual bool intersectsPointAtTime(double t, const MovingRegion& in) const;
		//virtual bool containsPointAtTime(double t, const MovingRegion& in) const;

		virtual void combineRegionInTime(const MovingRegion& r);
		virtual void combineRegionAfterTime(double t, const MovingRegion& r);
		virtual void getCombinedRegionInTime(MovingRegion& out, const MovingRegion& in) const;
		virtual void getCombinedRegionAfterTime(double t, MovingRegion& out, const MovingRegion& in) const;

		virtual double getIntersectingAreaInTime(const MovingRegion& r) const;
		virtual double getIntersectingAreaInTime(const libsupermesh::Tools::IInterval& ivI, const MovingRegion& r) const;

		//
		// IObject interface
		//
		virtual MovingRegion* clone();

		//
		// ISerializable interface
		//
		virtual uint32_t getByteArraySize();
		virtual void loadFromByteArray(const byte* data);
		virtual void storeToByteArray(byte** data, uint32_t& len);

		//
		// IEvolvingShape interface
		//
		virtual void getVMBR(Region& out) const;
		virtual void getMBRAtTime(double t, Region& out) const;

		//
		// ITimeShape interface
		//
		virtual double getAreaInTime() const;
		virtual double getAreaInTime(const libsupermesh::Tools::IInterval& ivI) const;
		virtual double getIntersectingAreaInTime(const ITimeShape& r) const;
		virtual double getIntersectingAreaInTime(const libsupermesh::Tools::IInterval& ivI, const ITimeShape& r) const;

		virtual void makeInfinite(uint32_t dimension);
		virtual void makeDimension(uint32_t dimension);

	private:
		void initialize(
			const double* pLow, const double* pHigh,
			const double* pVLow, const double* pVHigh,
			double tStart, double tEnd, uint32_t dimension);

	public:
		class CrossPoint
		{
		public:
			double m_t;
			uint32_t m_dimension;
			uint32_t m_boundary;
			const MovingRegion* m_to;

			struct ascending: public std::binary_function<CrossPoint&, CrossPoint&, bool>
			{
				bool operator()(const CrossPoint& __x, const CrossPoint& __y) const { return __x.m_t > __y.m_t; }
			};
		}; // CrossPoint

	public:
		double* m_pVLow;
		double* m_pVHigh;

		friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const MovingRegion& r);
	}; // MovingRegion

	typedef libsupermesh::Tools::PoolPointer<MovingRegion> MovingRegionPtr;
	SIDX_DLL std::ostream& operator<<(std::ostream& os, const MovingRegion& r);
} }
