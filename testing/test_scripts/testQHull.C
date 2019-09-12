#include <gtest/gtest.h>
#include "QHQuickHull.H"
#include "QHMathUtils.H"
#include <iostream>
#include <random>
#include <chrono>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

using namespace NEM::GEO::quickhull;

using FloatType = double;
using vec3 = Vector3<FloatType>;

// setup random number generator
static std::mt19937 rng;
static std::uniform_real_distribution<> dist(0,1);

class TestQHull : public testing::Test 
{
    protected:
    FloatType rnd(FloatType from, FloatType to) 
    {
        return from + (FloatType)dist(rng)*(to-from);
    }

    //void assertSameValue(FloatType a, FloatType b) 
    //{
    //    assert(std::abs(a-b)<0.0001f);
    //}
    
    template <typename T>
    std::vector<Vector3<T>> createSphere(T radius, size_t M, 
            Vector3<T> offset = Vector3<T>(0,0,0)) 
    {
        std::vector<Vector3<T>> pc;
        const T pi = M_PI;
        for (int i=0;i<=M;i++) 
        {
            FloatType y = std::sin(pi/2 + (FloatType)i/(M)*pi);
            FloatType r = std::cos(pi/2 + (FloatType)i/(M)*pi);
            FloatType K = FloatType(1) -
                std::abs((FloatType)((FloatType)i-M/2.0))/(FloatType)(M/2.0);
            const size_t pcount = (size_t)(1 + K*M + FloatType(1)/2);
            for (size_t j=0;j<pcount;j++) 
            {
                FloatType x = pcount>1 ? r*std::cos((FloatType)j/pcount*pi*2) : 0;
                FloatType z = pcount>1 ? r*std::sin((FloatType)j/pcount*pi*2) : 0;
                pc.emplace_back(x+offset.x,y+offset.y,z+offset.z);
            }
        }
        return pc;
    }


    virtual void SetUp() 
    {
        // Setup test env
        _N = 200;
        // Seed RNG using Unix time
        std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        auto seed = std::chrono::duration_cast<std::chrono::seconds>(
                now.time_since_epoch()).count();
        rng.seed((unsigned int)seed);

    }

    public:
    size_t _N;
    std::vector<vec3> _pc;
    QuickHull<FloatType> _qh;
    ConvexHull<FloatType> _hull;
};

TEST_F(TestQHull, Vector3) 
{
    vec3 a(1,0,0);
    vec3 b(1,0,0);

    vec3 c = a.projection(b);
    EXPECT_EQ( (c-a).getLength(), 0.0);
}

TEST_F(TestQHull, UnitCube) 
{
    // Put N points inside unit cube. Result mesh must have 
    // exactly 8 vertices because the convex hull is the unit cube.
    _pc.clear();
    for (int i=0;i<8;i++) 
      _pc.emplace_back(i&1 ? -1 : 1,i&2 ? -1 : 1,i&4 ? -1 : 1);
    for (size_t i=0; i<_N; i++)
      _pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
    _hull = _qh.getConvexHull(_pc,true,false);
    EXPECT_EQ(_hull.getVertexBuffer().size(), 8);
    // 6 cube faces, 2 triangles per face, 3 indices per triangle
    EXPECT_EQ(_hull.getIndexBuffer().size(), 3*2*6); 
    EXPECT_NE(&(_hull.getVertexBuffer()[0]), &(_pc[0]));
    auto hull2 = _hull;
    EXPECT_EQ(hull2.getVertexBuffer().size(), _hull.getVertexBuffer().size());
    EXPECT_EQ(hull2.getVertexBuffer()[0].x, _hull.getVertexBuffer()[0].x);
    EXPECT_EQ(hull2.getIndexBuffer().size(), _hull.getIndexBuffer().size());
    auto hull3 = std::move(_hull);
    EXPECT_EQ(_hull.getIndexBuffer().size(), 0);
}

TEST_F(TestQHull, UnitCubeOriginal) 
{
	// same test as unit cube, but using the original indices.
    // Put N points inside unit cube. Result mesh must have 
    // exactly 8 vertices because the convex hull is the unit cube.
    _pc.clear();
    for (int i=0;i<8;i++) {
    	_pc.emplace_back(i&1 ? -1 : 1,i&2 ? -1 : 1,i&4 ? -1 : 1);
    }
    for (size_t i=0; i<_N; i++)
    {
    	_pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
    }
    _hull = _qh.getConvexHull(_pc,true,true);
    EXPECT_EQ(_hull.getIndexBuffer().size(), 3*2*6);
    EXPECT_EQ(_hull.getVertexBuffer().size(), _pc.size());
    EXPECT_EQ(&(_hull.getVertexBuffer()[0]), &(_pc[0]));
}

TEST_F(TestQHull, UnitSphere) 
{
	// random N points from the boundary of unit sphere. 
    // Result mesh must have exactly N points.
	_pc = createSphere<FloatType>(1, 50);
	_hull = _qh.getConvexHull(_pc,true,false);
	EXPECT_EQ(_pc.size(), _hull.getVertexBuffer().size());
	_hull = _qh.getConvexHull(_pc,true,true);
	// Add every vertex twice. This should not affect final mesh
	auto s = _pc.size();
	for (size_t i=0;i<s;i++) {
		const auto& p = _pc[i];
		_pc.push_back(p);
	}
	_hull = _qh.getConvexHull(_pc,true,false);
	EXPECT_EQ(_pc.size()/2, _hull.getVertexBuffer().size());
}

//TEST_F(TestQHull, 0D) 
//{
//	const FloatType mul = 2*2*2;
//	while (true) {
//		for (auto& p : _pc) {
//			p.x *= mul;
//		}
//		_hull = _qh.getConvexHull(_pc,true,false);
//		if (_hull.getVertexBuffer().size() == 4) {
//			break;
//		}
//	}
//
//	_pc.clear();
//	vec3 centerPoint(2,2,2);
//	_pc.push_back(centerPoint);
//	for (size_t i=0;i<100;i++) {
//		auto newp = centerPoint + 
//            vec3(   rnd(-0.000001f,0.000001f),
//                    rnd(-0.000001f,0.000001f),
//                    rnd(-0.000001f,0.000001f) );
//		_pc.push_back(newp);
//	}
//	_hull = _qh.getConvexHull(_pc,true,false);
//	EXPECT_EQ(_hull.getIndexBuffer().size(), 12);
//}

TEST_F(TestQHull, PlanarCase) 
{
	QuickHull<FloatType> qh;
	std::vector<vec3> pc;
	pc.emplace_back(-3.000000f, -0.250000f, -0.800000f);
	pc.emplace_back(-3.000000f, 0.250000f, -0.800000f);
	pc.emplace_back(-3.125000f, 0.250000f, -0.750000);
	pc.emplace_back(-3.125000f, -0.250000f, -0.750000);
	auto hull = qh.getConvexHull(pc,true,false);
	EXPECT_EQ(hull.getIndexBuffer().size(), 12);
	EXPECT_EQ(hull.getVertexBuffer().size(), 4);
}

TEST_F(TestQHull, Cylinder) 
{
	// first a planar circle, then make a cylinder out of it
	_pc.clear();
	for (size_t i=0;i<_N;i++) {
		const FloatType alpha = (FloatType)i/_N*2*3.14159f;
		_pc.emplace_back(std::cos(alpha),0,std::sin(alpha));
	}
	_hull = _qh.getConvexHull(_pc,true,false);
	
	EXPECT_EQ(_hull.getVertexBuffer().size(), _pc.size());
	for (size_t i=0;i<_N;i++) {
		_pc.push_back(_pc[i] + vec3(0,1,0));
	}
	_hull = _qh.getConvexHull(_pc,true,false);
	EXPECT_EQ(_hull.getVertexBuffer().size(), _pc.size());
	EXPECT_EQ(_hull.getIndexBuffer().size()/3, _pc.size()*2-4);
}


TEST_F(TestQHull, Planes) 
{
	Vector3<FloatType> N(1,0,0);
	Vector3<FloatType> p(2,0,0);
	Plane<FloatType> P(N,p);
	auto dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(3,0,0), P);
	EXPECT_EQ(dist,1);
	dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
	EXPECT_EQ(dist,-1);
	dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
	EXPECT_EQ(dist,-1);
	N = Vector3<FloatType>(2,0,0);
	P = Plane<FloatType>(N,p);
	dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(6,0,0), P);
	EXPECT_EQ(dist,8);
}

TEST_F(TestQHull, HalfEdgeOutput) 
{
	QuickHull<FloatType> qh;
	
	// 8 corner vertices of a cube + tons of vertices inside. 
    // Output should be a half edge mesh with 12 faces 
    // (6 cube faces with 2 triangles per face) and 36 
    // half edges (3 half edges per face).
	std::vector<vec3> pc;
	for (int h=0;h<1000;h++) {
		pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
	}
	for (int h=0;h<8;h++) {
		pc.emplace_back(h&1?-2:2, h&2?-2:2, h&4?-2:2);
	}
	HalfEdgeMesh<FloatType, size_t> mesh = qh.getConvexHullAsMesh(&pc[0].x, pc.size(), true);
	EXPECT_EQ(mesh.m_faces.size(), 12);
	EXPECT_EQ(mesh.m_halfEdges.size(), 36);
	EXPECT_EQ(mesh.m_vertices.size(), 8);
}

TEST_F(TestQHull, Sphere) 
{
    QuickHull<FloatType> qh;
    FloatType y = 1;
    for (;;) {
        auto pc = createSphere<FloatType>(1, 100, Vector3<FloatType>(0,y,0));
        auto hull = qh.getConvexHull(pc,true,false);
        y *= 15;
        y /= 10;
        if (hull.getVertexBuffer().size()==4) {
            break;
        }
    }

    // Test worst case scenario: more and more points on the unit sphere. 
    // All points should be part of the convex hull, as long as we can make 
    // epsilon smaller without running out of numerical accuracy.
    bool failed = false;
    size_t i =  1;
    FloatType eps = 0.002;
    for (;;) {
        auto pc = createSphere<FloatType>(1, i, Vector3<FloatType>(0,0,0));
        auto hull = qh.getConvexHull(pc,true,false,eps);
        std::cout << i 
            << ":" << pc.size() 
            << " : " << hull.getVertexBuffer().size() 
            << " at eps=" << eps << std::endl;
        if (qh.getDiagnostics().m_failedHorizonEdges) {
            // This should not happen
            failed = true;
            break;
        }
        if (pc.size() == hull.getVertexBuffer().size()) {
            // Fine, all the points on unit sphere do belong to the convex mesh.
            i += 1;
        }
        else {
            eps *= 0.5f;
            std::cout << "Epsilon to " << eps << std::endl;
        }
        
        if (i == 140) {
            break;
        }
    }

    EXPECT_FALSE( failed );
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

