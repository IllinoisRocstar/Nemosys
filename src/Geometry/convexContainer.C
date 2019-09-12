#include"convexContainer.H"
#include<sstream>
#include<string.h>

using namespace std;

namespace NEM {

namespace GEO {

convexContainer::convexContainer():
    _isReady(false)
{}


convexContainer::convexContainer(std::vector<std::vector<double> >& inVrts)
{
    // deep copying
    setVertex(inVrts);
}

convexContainer::convexContainer(std::vector<quickhull::Vector3<double> >& inVrts)
{
    // deep copying
    for (auto itr=inVrts.begin(); itr!=inVrts.end(); itr++)
        vrts.emplace_back(*itr);
}

void convexContainer::setVertex(std::vector<std::vector<double> >& inVrts)
{
    // deep copying
    for (auto itr=inVrts.begin(); itr!=inVrts.end(); itr++)
        vrts.emplace_back((*itr)[0], (*itr)[1], (*itr)[2]);
}

void convexContainer::computeConvexHull()
{
    // computes/re-computes convex hull of the verices
    std::vector<size_t> indxBuf; 
    NEM::GEO::quickhull::VertexDataSource<double> vrtBuf; 
    qHull = new NEM::GEO::quickhull::QuickHull<double>();
    auto hull = qHull->getConvexHull(vrts, false, false, 1.e-14);
    indxBuf = hull.getIndexBuffer();
    vrtBuf = hull.getVertexBuffer();

    // converting to face-vector format
    //std::cout << "Size of index buffer = " 
    //    << indxBuf.size() << std::endl;

    // saving to vf
    for (int it=0; it<indxBuf.size(); it+=3)
    {
        Face f;
        for (int idx=0; idx<3; idx++)
        {
            NEM::MTH::Vector pnt;
            pnt.x = (vrtBuf[indxBuf[it+idx]]).x; 
            pnt.y = (vrtBuf[indxBuf[it+idx]]).y; 
            pnt.z = (vrtBuf[indxBuf[it+idx]]).z; 
            f.v.push_back(pnt);
        }
        fv.push_back(f);
    }

    _isReady = true;
    
}

bool convexContainer::isInConvexPoly(NEM::MTH::Vector const& p) {
    // check readiness
    if (!_isReady)
        computeConvexHull();
    
    // looping through faces
    for (Face const& f : fv) {
        NEM::MTH::Vector p2f = f.v[0] - p;         // f.v[0] is an arbitrary point on f
        double d = p2f.dot(f.normal());
        d /= p2f.norm();                 // for numeric stability

        constexpr double bound = -1e-15; // use 1e15 to exclude boundaries
        if (d < bound)
          return false;
    }
    return true;
}

bool convexContainer::isInConvexPoly(const std::vector<double>& p) {
    if (p.size() != 3)
        throw;
    NEM::MTH::Vector pnt;
    pnt.x = p[0]; 
    pnt.y = p[1]; 
    pnt.z = p[2]; 
    return(isInConvexPoly(pnt));
}

void convexContainer::toSTL(std::string file_name) const  {

    //ascii file
    std::ofstream myfile;
    myfile.open((file_name + ".stl").c_str(),  std::ios::out | std::ios::binary);
    myfile << "solid "
           << file_name
           << "\n";

    //write down every triangle
    for (auto it = fv.begin(); it!=fv.end(); it++){

        //normal vector coordinates
        NEM::MTH::Vector n = it->normal();
        myfile << "facet normal ";
        myfile << n.x << " " << n.y << " " << n.z << "\n";

        myfile << "\touter loop\n";
        //vertex coordinates
        for (int iv=0; iv<3; iv++)
        {
            myfile << "\t\tvertex "
                   << (it->v[iv]).x << " "
                   << (it->v[iv]).y << " "
                   << (it->v[iv]).z << "\n";
        }
        myfile << "\tendloop\n";
        myfile << "endfacet\n";
    }
    myfile << "endsolid " << file_name << "\n";
    myfile.close();

}

} // namespace GEO

} // namespace NEM
