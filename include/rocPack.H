#ifndef ROCPACK_H
#define ROCPACK_H

#include "nemosys_export.h"
#include "meshBase.H"

#include <string>

namespace NEM {
    
namespace GEO {

/**
    @brief A struct to respresent quaternion.
    @param Components of a quaternion (q = w + xi + yj + zk)
**/
struct rocQuaternion
{
    double w, x, y, z;
};

/** @brief This class converts Rockpack output file into periodic geometry and 
           writes into STL, VTK, .MSH, and ExodusII file format. This class can also
           be used for generating 3D periodic tetrahedral mesh using periodic 
           geometry. Currently the shapes supported are as follows.

           - Shapes
              -# Sphere
              -# Ellipsoid
              -# Cylinder
              -# HMX (Octogen crystal)
              -# PETN (Pentaerythritol tetranitrate crystal)
              -# ICOSIDODECAHEDRON

           Methods for creating spheres, cylinders and ellisoids are contained
           in this class. For other shapes, methods wrapped around abstract 
           class "rocPackShape" are used.
**/

// TODOs 
// 1. Boolean(OCC) is not working with ellipsoid volumes created by GEO 
//    methods. Find a way to make them periodic using some other boolean methods.
// 2. createCohesiveElement method needs attention. Finding proper set of interface
//    node is important. One idea is to take all surrounding nodes or take geometry
//    mesh nodes and substract periodic nodes from that.
// 3. Cohesive elements needs to be extended to finite thickness. This should be 
//    fairly very easy once we have zero thickness method working flawlessly. Also,
//    this method only works with just 2 physical groups enabled.
// 4. Affine transfrom is not set up for anything other than box of 1 X 1 X 1.

class NEMOSYS_EXPORT rocPack
{
  // Constructors, Deconstructor, Access
  public:
  /** @brief rocPack default constructor. Options for periodic geometry
             and removing boundary packs are off by default. To enable these
             options, user can specify following keywords anywhere in RocPack
             output file. (Default for both keyword is false.)
             - "SetPeriodicity = true"
             - "RemoveBoundaryPacks = true"

             These can also be switched on using newly added public methods.
             - removeBoundaryVolumes()
             - setPeriodicGeometry()
             - setPeriodicMesh()

      @param fname Rocpack output file name
      @param outName Output STL, VTK, .MSH, .exo file names
  **/
  rocPack(const std::string &fname, const std::string &outName);

  /** @brief rocPack standard destructor
  **/
  ~rocPack()
  { std::cout << "rockPack class destroyed" << std::endl; }

  /** @brief rocPack2Surf method converts Rocpack output file into STL/VTK 
             surface mesh and writes into current folder.
  **/
  void rocPack2Surf();

  /** @brief This method creates cohesive elements using existing vtu file
      @param filename Input vtu file name
      @param output file name
  **/
  void createCohesiveElements(
      const std::string & filename,const std::string & outname);

  /** @brief rocPack2Periodic3D method converts Rocpack output file into
             3D periodic mesh. Writes the output file in current path
  **/
  void rocPack2Periodic3D();

  /** @brief Enables an option to remove volumes intersecting boundary
  **/
  void removeBoundaryVolumes();

  /** @brief Enables an option to enforce periodicity on generated geometry.
  **/
  void setPeriodicGeometry();

  /** @brief Enables an option to enforce periodicity on generated mesh.
  **/
  void setPeriodicMesh();

  /** @brief Shrinks the volumes by defined percentage
      @param percntg Percentage between 0 and 1 by which the pack should be 
             shrinked
  **/
  void shrinkVolumes(const double percntg);

  /** @brief Performs surface smoothing over every volume in pack geometry
      @param smoothingParam Smoothing Iterations.
  **/
  void smoothSurfaces(const int smoothingParam);

  /** @brief Set mesh size for periodic mesh
      @param size Mesh Size
  **/
  void setMeshSize(const double size);

  /** @brief Enables physical grouping in final mesh
  **/
  void enablePhysicalGrps();

  /** @brief Enables two physical groups in final mesh
  **/
  void enableTwoPhysGrps();

  /** @brief Enables 6 surface patches for all sides of box
  **/
  void enableSurfacePatches();

  /** @brief Translates the whole geometry along user-defined parameters
      @param X Coordinate
      @param Y Coordinate
      @param Z Coordinate
  **/
  void translateAll(const double& X, const double& Y, const double& Z);

  /** @brief Sets custom domain size defined by users
  **/
  void setCustomDomain(const std::vector<double>& domainBounds);

  /** @brief Sets meshing algorithm of user's choice
  **/
  void setMeshingAlgorithm(const int& mshAlg);

  /** @brief For internal testing purpose.
  **/
  void setNodeLocations(const int &x, const int &y, const int &z);

  /** @brief For internal testing purpose.
  **/
  void setRandomSurface(const int &surf);

  /** @brief This method is for internal testing purpose only. Compares spatial 
             location of several periodic nodes in periodic mesh.
      @return Boolean for test
  **/
  bool getTestResult();

  // rocPack2STL internal use only
  private:
  /** @brief This method parses Rocpack output file and collects data for 
             creation of pack geometries
  **/
  void rocParser();

  /** @brief This method obtains parsed data from Rocpack output and creates
             periodic pack geometries.
  **/
  void rocToGeom();

  /** @brief This method creates surface mesh for final geometry and writes 
             into surface files (.stl, .vtk, .msh)
  **/
  void geomToSurf();

  /** @brief Performs operations on generated geometry
  **/
  void geomToVolFrac();


  // rocParser Internal methods
  private:
  /** @brief Finds particular word in file stream
      @param word Word that user wants to find
      @return line containing that word in file
  **/
  std::string findWord(const std::string &word);

  /** @brief Tokenizes string into words and stores them into a vector.
      @param String which needs to be tokenized
      @param Delimiter
      @return A vector of token strings
  **/
  std::vector<std::string> Tokenize(
     const std::string &lineIn, const char &delim);

  /** @brief Tokenizes string into words and stores them into a vector.
      @param String which needs to be tokenized
      @param Vector of delimiters
      @return A vector of token strings
  **/
  std::vector<std::string> Tokenize(
      const std::string &lineIn, const std::string &delims);

  /** @brief Reads the current line and returns the translate, rotate, and
             scale data for the pack shape
      @param iter Line number to read in file
      @param a string of delimiters for parsing
      @param L Line to read
      @return Vector of tokens for input line L
  **/
  std::vector<std::string> getShapeData(
      const int &iter, const std::string &a, const std::vector<std::string> &L);

  /** @brief Converts parsed string to double
      @param strng String to convert
      @return Char array
  **/
  char* strToChar(const std::string &strng);


  // rocToGeom and geomToSurface Internal methods
  private:
  /** @brief Initializes GMSH for workflow
  **/
  void initialize();

  /** @brief Makes geometry periodic
      @param rmbPacks boolean for removing boundary packs
  **/
  void makePeriodic(const bool rmbPacks);

  /** @brief Normalizes shape vertices
  **/
  void normalizeVerts();

  /** @brief Scales all volumes by shrinking percentage
  **/
  void scaleVols(const int &vol, const int &index);

  /** @brief Performs laplacian smoothing on all packs if requested
  **/
  void performSmoothing();

  /** @brief This method writes periodic pack geometries into STL format
      @param writeFile Output File Name
  **/
  void geomToSTL(const std::string &writeFile);

  /** @brief This method writes periodic pack geometries into VTK format
      @param writeFile Output File Name
  **/
  void geomToVTK(const std::string &writeFile);

  /** @brief This method writes periodic pack geometries into .MSH format
      @param writeFile Output File Name
  **/
  void geomToMsh(const std::string &writeFile);

  /** @brief Generates a Quaternion using Euler 3D rotation
             parameters
      @param r Rotation parameters in form (x,y,z,angle)
      @return Quaternion in form \f$ q = w + x i + y j + z k \f$.
  **/
  rocQuaternion toQuaternion(const std::vector<double> &r);

  /** @brief Rotates a 3D Vector by unit Quaternion. First, converts
             input vector into quaternion form by making "w" component
             zero (i.e \f$ v = 0 + x i + y j + z k \f$). Then performs
             transformation using \f$ v_{'} = q.v.q\textsuperscript{-1} \f$.
             if we write quaternion in this form \f$ q = w + p(x,y,z) \f$,
             then \f$ v_{'} = 2\(p.v\)p + \(w_2 - p.p\)v + 2w\(p\timesv\) \f$.
             This method is computationally efficient compared to performing
             direct quaternion products (i.e 
             \f$ v_{'} = q.v.q\textsuperscript{-1} \f$)
      @param q A Quaternion in form \f$ q = w + x i + y j + z k \f$.
      @param v Vector that needs to be rotated
      @return Rotated vector
  **/
  std::vector<double> rotateByQuaternion(
      const rocQuaternion &q, const std::vector<double> &v);

  /** @brief This method removes the pack shapes intersecting boundary
      @param n Index for pack number
  **/
  void tagBoundaryPacks();


  // Internal methods for 3D periodic mesh module
  private:
  /** @brief This method maps surfaces in X,Y, and Z direction which are 
             periodic and then enforces that into gmsh model
      @param prevTags Currently existing model entities
  **/
  void mapPeriodicSurfaces(const std::vector<std::pair<int,int>>& prevTags);

  /** @brief This method writes 3D periodic mesh into .msh and vtu format
  **/
  void geomToPeriodic3D();

  /** @brief This method gets points that belongs to surfaces
      @param surfaces A vector pair of surface tags
      @return Vector of vertices.
  **/
  std::vector<std::vector<std::pair<int,int>>> getAllPoints(
                      std::vector<std::pair<int,int>> surfaces);

  /** @brief This method find which surfaces are periodic on opposite sides
      @param vertsOneSide Vertices on one side
      @param vertsOtherSide Vertices on other side
      @param indexTranslate Index for which direction to translate
      @param amountTranslate Distance for translation
      @return Vector pair of periodic surfaces
  **/
  std::vector<std::pair<int,int>> getPeriodicSurfs(
      const std::vector<std::vector<std::pair<int,int>>> & vertsOneSide,
      const std::vector<std::vector<std::pair<int,int>>> & vertsOtherSide,
      const int & indexTranslate,
      const double & amountTranslate);

  /** @brief This method allows mapping of nodes on periodic boundaries. Writes
             periodic nodes surface wise in CSV file.
  **/
  void writePeriodicNodes();


  // Basic shape methods
  private:
    /** @brief Creates sphere pack shapes
        @param n Index for pack number
  **/
  void makeSphere(const int &n);

  /** @brief Creates Ellipsoid pack shapes
      @param n Index for pack number
  **/
  void makeEllipsoid(const int &n);

  /** @brief Creates cylinder pack shapes
      @param n Index for pack number
  **/
  void makeCylinder(const int &n);

  /** @brief Creates a crystal shape based on Rocpack output file
      @param n Index for pack number
      @param index Index for pack data access
  **/
  void makeCrystalShape(const int &n, const int &index);

  /** @brief creates VTK cell
  **/
  void createVtkCell(
    vtkSmartPointer<vtkUnstructuredGrid> dataSet,
    const int cellType, std::vector<int>& vrtIds);


  // Pointers and/or parameters
  private:
  /** @brief rocPack output file name
  **/
  std::string InFile;

  /** @brief Output STL/VTK file name
  **/
  std::string OutFile;

  /** @brief Vector for storing base shape names
  **/
  std::vector<std::string> shapeNames;

  /** @brief Vector for storing shapes
  **/
  std::vector<std::string> uniqueNames;

  /** @brief Vector for storing crystal names
  **/
  std::vector<std::string> crystalNames;

  // Box Parameters
  /** @brief X dimension of box geometry
  **/    
  double Xdim;

  /** @brief Y dimension of box geometry
  **/  
  double Ydim;

  /** @brief Z dimension of box geometry
  **/  
  double Zdim;

  /** @brief Vector of box starting coordinates
  **/
  std::vector<double> boxPt;

  // Pack shapes parameters
  /** @brief Vector of Ellipsoid radii
  **/
  std::vector<double> ellipsoidRad;

  /** @brief Vector of cylinder parameters
  **/
  std::vector<double> cylParams;

  /** @brief Vector of crystal shape vertices
  **/
  std::vector<std::vector<std::vector<double>>> verts;

  /** @brief Vector of crystal shape faces
  **/
  std::vector<std::vector<std::vector<int>>> faces;

  /** @brief Vector of translate coordinates for all packs
  **/
  std::vector<std::vector<double>> translateParams;

  /** @brief Vector of rotate coordinates for all packs
  **/
  std::vector<std::vector<double>> rotateParams;

  /** @brief Vector of scales for all packs
  **/
  std::vector<double> scaleOfPack;

  /** @brief Vector of pack shape names
  **/
  std::vector<std::string> nameOfPcks;

  /** @brief Name of pack shape
  **/
  std::string packName;

  /** @brief Ellipsoid packs cannot be made periodic as of now.
             This boolean will be true when ellipsoid packs are present.
  **/
  bool ellipsoidPresent = false;

  /** @brief Boolean to opt for non-periodic geometry (as generated by rocPack)
  **/
  bool enablePeriodicity = false;

  /** @brief Boolean to opt for removing packs on boundaries. (Resultant  
             geometry will automatically be periodic)
  **/
  bool removeBoundaryPacks = false;

  /** @brief Stores volume index of packs interseting boundary.
  **/
  std::vector<std::pair<int,int>> bndryPackTags;

  /** @brief Boolean for 3D periodic mesh generation
  **/
  bool periodic3D = false;

  /** @brief Surfaces linked to master periodic surface in X direction
  **/
  std::vector<int> slaveX;

  /** @brief Surfaces linked to master periodic surface in Y direction
  **/
  std::vector<int> slaveY;

  /** @brief Surfaces linked to master periodic surface in Z direction
  **/
  std::vector<int> slaveZ;

  /** @brief Master periodic surfaces in X direction
  **/
  std::vector<int> MasterX;

  /** @brief Master periodic surfaces in Y direction
  **/
  std::vector<int> MasterY;

  /** @brief Master periodic surfaces in Z direction
  **/
  std::vector<int> MasterZ;

  /** @brief For testing
  **/
  int randomNodeX = 0;

  /** @brief For testing
  **/
  int randomNodeY = 0;

  /** @brief For testing
  **/
  int randomNodeZ = 0;

  /** @brief Nodes for testing
  **/
  bool matchingCoords = false;

  /** @brief Shrink Scale
  **/
  double shrinkScale = 1.0;

  /** @brief some smoothing parameter
  **/
  int smoothingIter = 0;

  /** @brief Enables smoothing
  **/
  bool enableSmoothing = false;

  /** @brief Enables Shrinking
  **/
  bool enableScaling = false;

  /** @brief Mesh Size
  **/
  double meshSz = 0.1;

  /** @brief Enable Physical Grouping in Final Mesh
  **/
  bool enablePhysGrp = false;

  /** @brief Enable Two Physical Groups in Final Mesh
  **/
  bool just2Physgrps = false;

  /** @brief Stores node tags for surronding physical group
  **/
  std::vector<std::size_t> surrNodeTags;

  /** @brief Stores node coordinates for surronding physical group
  **/
  std::vector<double> surrCoords;

  /** @brief Stores node tags for geometries physical group
  **/
  std::vector<std::size_t> geomsNodeTags;

  /** @brief Stores node coordinates for geometries physical group
  **/
  std::vector<double> geomsCoords;

  /** @brief Tag number
  **/
  int surroundingGrp;

  /** @brief Tag number
  **/
  int packGrp;

  /** @brief Vector for stroing nodes at conformal interface
  **/
  std::vector<std::size_t> interfaceNodes;

  /** @brief Vector for storing phys group tags
  **/
  std::vector<int> surrondingEntities;

  /** @brief Vector for storing phys group tags
  **/
  std::vector<int> packEntities;

  /** @brief A boolean for user to enable assignement of surface patches
  **/
  bool assignSidePatches = false;

  /** @brief X coordinate for user-defined translate.
  **/
  double xUDF = 0;

  /** @brief Y coordinate for user-defined translate
  **/
  double yUDF = 0;
  
  /** @brief Z coordinate for user-defined translate
  **/
  double zUDF = 0;

  /** @brief Custom Domain Boolean
  **/
  bool cstmDomain = false;

  /** @brief Random location for a surface
  **/
  int randomSurfTest = 0;

  /** @brief Meshing Algorithm
  **/
  int meshingAlgorithm = 1;


};

}

}

#endif