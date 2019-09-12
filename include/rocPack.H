#ifndef ROCPACK_H
#define ROCPACK_H

#include "nemosys_export.h"

#include <string>

namespace NEM {
    
namespace GEO {

/**
    @brief A struct to respresent quaternion.
**/
struct rocQuaternion
{
    double w, x, y, z; /**< Components of a quaternion **/
};

/** @brief This class converts Rockpack output file into periodic geometry and 
           writes into STL, VTK, and .MSH file format. Currently the shapes supported
           are as follows.

           ## Shapes
           > * Sphere
           > * Ellipsoid
           > * Cylinder
           > * HMX (Octogen crystal)
           > * PETN (Pentaerythritol tetranitrate crystal)
           > * ICOSIDODECAHEDRON

           Methods for creating spheres, cylinders and ellisoids are contained
           in this class. For other shapes, methods wrapped around abstract 
           class "rocPackShape" are used.
**/

// TODOs Boolean(OCC) is not working with ellipsoid volumes created by GEO 
// methods. Find a way to make them periodic using some other boolean methods.

class NEMOSYS_EXPORT rocPack
{
  // Constructors, Deconstructor, Access
  public:
  /** @brief rocPack default constructor. Options for non-periodic geometry
             and removing boundary packs are off by default. To enable these
             options, user can specify following keywords anywhere in RocPack
             output file. (Default for both keyword is false.)
             > * "NoPeriodicity = true"
             > * "RemoveBoundaryPacks = true"
      @param fname Rocpack output file name
      @param outName Output STL, VTK, and .MSH file names
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
  std::vector<std::string> Tokenize(const std::string &lineIn, 
                                    const char &delim);

  /** @brief Tokenizes string into words and stores them into a vector.
      @param String which needs to be tokenized
      @param Vector of delimiters
      @return A vector of token strings
  **/
  std::vector<std::string> Tokenize(const std::string &lineIn, 
                                    const std::string &delims);

  /** @brief Reads the current line and returns the translate, rotate, and
             scale data for the pack shape
      @param iter Line number to read in file
      @param a string of delimiters for parsing
      @param L Line to read
      @return Vector of tokens for input line L
  **/
  std::vector<std::string> getShapeData(const int &iter, 
                                        const std::string &a,
                                        const std::vector<std::string> &L);

  /** @brief Converts parsed string to double
      @param strng String to convert
      @return Char array
  **/
  char* strToChar(const std::string &strng);


  // rocToGeom Internal methods
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
  std::vector<double> rotateByQuaternion(const rocQuaternion &q,
                                         const std::vector<double> &v);

  /** @brief This method removes the pack shapes intersecting boundary
      @param n Index for pack number
  **/
  void tagBoundaryPacks();


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
  bool noPeriodicity = false;

  /** @brief Boolean to opt for removing packs on boundaries. (Resultant  
             geometry will automatically be periodic)
  **/
  bool removeBoundaryPacks = false;

  /** @brief Stores volume index of packs interseting boundary.
  **/
  std::vector<std::pair<int,int>> bndryPackTags;

};

}

}

#endif