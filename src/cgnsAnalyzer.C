#include <cgnsAnalyzer.H>
#include <GmshEntities.h>
#include <string.h>
#include <iostream>
#include <meshBase.H>
#include <algorithm>


/********************************************
    solutionData class implementation
*********************************************/

void solutionData::appendData(const vecSlnType& inBuff, int inNData, int inNDim)
{
  // sanity check
  if (!slnData.empty())
  {
    if (nDim != inNDim || inNDim == 0)
    {
      std::cerr << "Incompatible data size can not be appended."
                << "Current nDim = " << nDim << " Requested nDim = "
                << inNDim << std::endl;
      return;
    }
  }
  else
  {
    nDim = inNDim;
  }

  // deep copy data
  slnData.insert(slnData.end(), inBuff.begin(), inBuff.end());
  nData += inNData;
}

void solutionData::getData(vecSlnType& outBuff, int& outNData, int& outNDim)
{
  // deep copy data
  outBuff.insert(outBuff.begin(), slnData.begin(), slnData.end());
  outNData = nData;
  outNDim = nDim;
}

// only copies data that are defined in the given mask
void solutionData::getData(vecSlnType& outBuff, int& outNData, int& outNDim,
                           std::vector<bool> mask)
{
  // deep copy data
  int cntr1 = 0;
  int cntr2 = 0;
  for (auto &&id : mask)
  {
    if (id)
    {
      outBuff.insert(outBuff.end(), slnData[cntr1]);
      cntr2++;
    }
    cntr1++;
  }
  outNData = cntr2;
  outNDim = nDim;
}

void solutionData::rmvDataIdx(const std::vector<int> rmvIdx)
{
  vecSlnType nsdv;
  int nnd = 0;
  for (int id=0; id<nData; id++)
  {
    auto it = std::find(rmvIdx.begin(), rmvIdx.end(), id);
    if (it != rmvIdx.end())
      continue;
    else
    {
      nnd++;
      for (int iDim=0; iDim<nDim; iDim++)
        nsdv.push_back(slnData[id*nDim + iDim]); 
    }
  }
  slnData = nsdv;
  nData = nnd; 
}


/********************************************
    cgnsAnalyzer class implementation
*********************************************/

void cgnsAnalyzer::loadGrid(int verb)
{
  // cgns related variables
  int i, j, k;
  char basename[33], zonename[33];

  // open CGNS file for read
  if (cg_open(cgFileName.c_str(), CG_MODE_MODIFY, &indexFile)) cg_error_exit();

  // reading base information
  cg_nbases(indexFile, &nBase);
  if (nBase > 1)
    std::cerr << "There are " << nBase << " bases in the file (not supported).\n";
  else
    indexBase = 1;
  if (cg_base_read(indexFile, indexBase, basename, &cellDim, &physDim)) cg_error_exit();
  baseName = basename;

  // reading units
  if (cg_goto(indexFile, indexBase, "end")) cg_error_exit();
  if (cg_units_read(&massU, &lengthU, &timeU, &tempU, &angleU)) cg_error_exit();
  if (verb > 0)
    std::cout << " ------------------------------------------------\n"
              << "Base name = " << baseName
              << "\ncellDim = " << cellDim
              << "\nphysDim = " << physDim
              << "\nUnits " << massU << " " << lengthU << " " << timeU << " " << tempU << " " << angleU
              << std::endl;

  // reading base iterative data
  char bitername[33];
  if (cg_biter_read(indexFile, indexBase, bitername, &nTStep)) cg_error_exit();
  baseItrName = bitername;
  if (nTStep != 1)
    std::cerr << "More than one time step is not supported.\n";
  int nArrays;
  std::cout << baseItrName << std::endl;
  if (cg_goto(indexFile, indexBase, bitername, 0, "end")) cg_error_exit();
  if (cg_array_read_as(1, CG_RealDouble, &timeLabel)) cg_error_exit();
  if (verb > 0) std::cout << "Time label = " << timeLabel << std::endl;

  // reading zone information
  cg_nzones(indexFile, indexBase, &nZone);
  if (nZone > 1)
  {
    //std::cout << "There are " << nZone << " zones in the file (not supported).\n";
    isMltZone = true;
  }
  loadZone(1, verb);
}

void cgnsAnalyzer::loadZone(int zIdx, int verb)
{
  char zonename[33];
  indexZone = zIdx;
  // reading zone type and name
  cg_zone_read(indexFile, indexBase, indexZone, zonename, cgCoreSize);
  cg_zone_type(indexFile, indexBase, indexZone, &zoneType);
  zoneName = zonename;
  zoneNames.push_back(zoneName);
  if (verb) std::cout << "Zone name = " << zoneName << std::endl;

  // type of zone
  if (zoneType == CG_ZoneTypeNull)
  {
    if (verb) std::cout << "Zone type = CG_ZoneTypeNull\n";
  }
  else if (zoneType == CG_Structured)
  {
    isUnstructured = false;
    if (verb) std::cout << "Zone type = CG_Structured\n";
    std::cerr << "Warning: Structured meshes only supported experimentally.\n";
  }
  else if (zoneType == CG_Unstructured)
  {
    isUnstructured = true;
    if (verb) std::cout << "Zone type = CG_Unstructured\n";
  }
  else if (zoneType == CG_ZoneTypeUserDefined)
  {
    if (verb) std::cout << "Zone type = CG_ZoneTypeUserDefined\n";
  }

  // reading zone iterative data, in the case of Rocstar outputs
  // there should be two additional nodes bellow this
  // containing grid coordinate pointer and flow solution pointer
  char zitername[33];
  if (cg_ziter_read(indexFile, indexBase, indexZone, zitername)) cg_error_exit();
  zoneItrName = zitername;
  if (verb) std::cout << "Zone iterative name = " << zoneItrName << std::endl;
  char gridcoorpntr[33], flowslnpntr[33];
  if (cg_goto(indexFile, indexBase, zonename, 0, zitername, 0, "end")) cg_error_exit();
  // readin flow solution and grid coordiate pointers if they exist at all
  //if (!cg_goto(indexFile, indexBase, zonename, 0, zitername, 0, "GridCoordiatesPointers", 0, "end"))
  if (!cg_goto(indexFile, indexBase, "Zone_t", indexZone, "ZoneIterativeData_t", 1, "end"))
  {
    if (!cg_array_read_as(1,  CG_Character, gridcoorpntr))
      gridCrdPntr = gridcoorpntr;
    else
      gridCrdPntr = "NA";
    if (verb) std::cout << "Grid coordinate pointer = " << gridcoorpntr << std::endl;
    if (!cg_array_read_as(2,  CG_Character, flowslnpntr))
      flowSlnPntr = flowslnpntr;
    else
      flowSlnPntr = "NA";
    if (verb) std::cout << "Flow solution pointer = " << flowslnpntr << std::endl;
  }
  else
  {
    std::cout << "This zone does not containing grid/flow pointers.\n";
  }

  // finding number of vertices and elements
  if (zoneType == CG_Unstructured)
  {
    nVertex = cgCoreSize[0];
    nElem = cgCoreSize[1];
    rmax[0] = nVertex;
    rmax[1] = nVertex;
    rmax[2] = nVertex;
  }
  else
  {
    // for strucutred meshes we have to 
    // take into account number of rind nodes in each
    // direction. In rocstar case they are equal in each direction
    int rdata[2];
    if (cg_goto(indexFile, indexBase, "Zone_t", zIdx,
              "GridCoordinates", 0, "end"))
        cg_error_exit();
    if (cg_rind_read(rdata))
        cg_error_exit();
    nRindNdeStr = rdata[1]; 
    // assuming same number in each direction and same number
    // of rind cells layers as the number of rind node layers
    // (true condition for Rocflo)
    std::cout << "nRindNdeStr = " << nRindNdeStr << "\n";
    // number of real vertices + nRindNode in each direction
    cgCoreSize[0] += 2*nRindNdeStr;
    cgCoreSize[1] += 2*nRindNdeStr;
    cgCoreSize[2] += 2*nRindNdeStr;
    // each rind node adds a row of elements as well
    cgCoreSize[3] += 2*nRindNdeStr;
    cgCoreSize[4] += 2*nRindNdeStr;
    cgCoreSize[5] += 2*nRindNdeStr;
    nVertex = cgCoreSize[0] * cgCoreSize[1] * cgCoreSize[2] ; 
    nElem = cgCoreSize[3] * cgCoreSize[4] * cgCoreSize[5];
    rmin[0] = 1;
    rmin[1] = 1;
    rmin[2] = 1;
    rmax[0] = cgCoreSize[0];
    rmax[1] = cgCoreSize[1];
    rmax[2] = cgCoreSize[2];
    std::cout << " Number of rind node for structured zone = " << nRindNdeStr << std::endl;
  }
  if (verb)
    std::cout << "Number of vertex = " << nVertex
              << "\nNumber of elements = " << nElem
              << std::endl;

  // capturing rinde node ids
  int ndId = 0;
  if (nRindNdeStr != 0)
    for (int i=rmin[2]; i<=rmax[2]; i++)
      for (int j=rmin[1]; j<=rmax[1]; j++)
        for (int k=rmin[0]; k<=rmax[0]; k++)
        {
          ndId++;
          if ( (i<=nRindNdeStr || (rmax[2]-i)<nRindNdeStr) ||
               (j<=nRindNdeStr || (rmax[1]-j)<nRindNdeStr) ||
               (k<=nRindNdeStr || (rmax[0]-k)<nRindNdeStr) )
              cgRindNodeIds.push_back(ndId);
        }
  std::cout << "Total number of rind nodes = " << cgRindNodeIds.size() << "\n";

  // reading coordinates
  int one = 1;
  if (zoneType == CG_Unstructured)
  {
    // reading coordinates X
    xCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateX", CG_RealDouble, &one, &rmax[0], &xCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;

    // reading coordinates Y
    yCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateY", CG_RealDouble, &one, &rmax[1], &yCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;

    // reading coordinates Z
    zCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateZ", CG_RealDouble, &one, &rmax[2], &zCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
  }
  else if (zoneType == CG_Structured)
  {
    // rind information in case any (testing)

    // reading coordinates X
    xCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateX", CG_RealDouble, &rmin[0], &rmax[0], &xCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;

    // reading coordinates Y
    yCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateY", CG_RealDouble, &rmin[0], &rmax[0], &yCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;

    // reading coordinates Z
    zCrd.resize(nVertex, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateZ", CG_RealDouble, &rmin[0], &rmax[0], &zCrd[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
  }
  else
  {
    std::cerr << "Error in load, only CG_Structured and CG_Unstructured girds are supported.\n ";
    exit(0);
  }

  // reading connectivity only if CG_Unstructured
  if (zoneType == CG_Unstructured)
  {
    int indexSection, eBeg, eEnd, nBdry, parentFlag;
    char sectionname[33];
    if (cg_nsections(indexFile, indexBase, indexZone, &nSection) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error();
    if (nSection > 1)
      if (verb)
        std::cout << nSection << " element connectivity sections were found (real+virtual in case of rocstar)"
                  << ", only first section will be read for now." << std::endl;
    indexSection = 1;
    if (cg_section_read(indexFile, indexBase, indexZone, indexSection,
                        sectionname, &sectionType, &eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    sectionName = sectionname;
    if (verb)
      std::cout << "Section " << sectionname
                << " eBeg = " << eBeg
                << " eEnd = " << eEnd
                << " nBdry = " << nBdry
                << std::endl;
    switch (sectionType)
    {
      case CG_TETRA_4:
        nVrtxElem = 4;
        break;
      case CG_HEXA_8:
        nVrtxElem = 8;
        break;
      case CG_TRI_3:
        nVrtxElem = 3;
        break;
      case CG_QUAD_4:
        nVrtxElem = 4;
        break;
      case CG_TETRA_10:
        nVrtxElem = 10;
        break;
      default:
        std::cerr << "Unknown element type " << sectionType << std::endl;
        break;
    }
    elemConn.resize(nVrtxElem * nElem, -1);
    if (cg_elements_read(indexFile, indexBase, indexZone, indexSection, &elemConn[0], nullptr) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    // reduce by 1 to base the node index to zero
    // using lambda function
    //std::for_each(elemConn.begin(), elemConn.end(), [](int& d) { d -= 1; });
    if (verb) std::cout << "Size of connectivity vector = " << elemConn.size() << std::endl;
  }
  else if (zoneType == CG_Structured)
  {
    // Generating connectivity for structured meshes since file does not contain data
    std::cerr << "Warning: converting implicit structured connectivity to 8-Node hexahedrals.\n";
    sectionType = CG_HEXA_8;
    elemConn.resize(8*nElem, -1);
    int iElm;
    int bs,d1,d2,d3,d4,d5;
    iElm = 0;
    bs = 0;
    d4 = cgCoreSize[1]*cgCoreSize[0];
    for (int k=1; k < cgCoreSize[2]; k++)
    {
      d1 = d4*(k-1);      
      for (int j=1; j < cgCoreSize[1]; j++)
      {
        d2 = d1 + (cgCoreSize[0])*(j-1);
        for (int i=1; i < cgCoreSize[0]; i++)
        {
          // rind detector
          d3 = d2 + i;
          elemConn[bs] = d3;
          elemConn[bs+1] = d3 + 1;
          elemConn[bs+2] = d3 + 1 + cgCoreSize[0];
          elemConn[bs+3] = d3 + cgCoreSize[0];
          elemConn[bs+4] = d3 + d4;
          elemConn[bs+5] = d3 + 1 + d4;
          elemConn[bs+6] = d3 + 1 + cgCoreSize[0] + d4;
          elemConn[bs+7] = d3 + cgCoreSize[0] + d4;
          //std::cout << "Elem conn = "
          //    << elemConn[bs] << " "
          //    << elemConn[bs+1] << " "
          //    << elemConn[bs+2] << " "
          //    << elemConn[bs+3] << " "
          //    << elemConn[bs+4] << " "
          //    << elemConn[bs+5] << " "
          //    << elemConn[bs+6] << " "
          //    << elemConn[bs+7] << "\n";
          bs += 8;
          iElm +=1;
          if ( (k <= nRindNdeStr || (cgCoreSize[2]-k) <= nRindNdeStr) ||
               (j <= nRindNdeStr || (cgCoreSize[1]-j) <= nRindNdeStr) ||
               (i <= nRindNdeStr || (cgCoreSize[0]-i) <= nRindNdeStr) )
             cgRindCellIds.push_back(iElm);
        }
      }
    }
    //auto it1 = max_element(std::begin(elemConn), std::end(elemConn));
    //auto it2 = min_element(std::begin(elemConn), std::end(elemConn));
    //std::cout << "Max Conn = " << *it1 << std::endl;
    //std::cout << "Min Conn = " << *it2 << std::endl;
    //std::cerr << "nElm = " << iElm << "\n";
  }
}

int cgnsAnalyzer::getIndexFile()
{
  return indexFile;
}

int cgnsAnalyzer::getIndexBase()
{
  return indexBase;
}

int cgnsAnalyzer::getCellDim()
{
  return cellDim;
}

std::string cgnsAnalyzer::getFileName()
{
  return cgFileName;
}

std::string cgnsAnalyzer::getBaseName()
{
  return baseName;
}

std::string cgnsAnalyzer::getZoneName()
{
  return zoneName;
}

std::string cgnsAnalyzer::getZoneName(int iCgFile)
{
  return zoneNames[iCgFile];
}

std::string cgnsAnalyzer::getSectionName()
{
  return sectionName;
}

std::string cgnsAnalyzer::getBaseItrName()
{
  return baseItrName;
}

int cgnsAnalyzer::getNZone()
{
  return nZone;
}

int cgnsAnalyzer::getNTStep()
{
  return nTStep;
}

std::string cgnsAnalyzer::getZoneItrName()
{
  return zoneItrName;
}

std::string cgnsAnalyzer::getGridCrdPntr()
{
  return gridCrdPntr;
}

std::string cgnsAnalyzer::getSolutionPntr()
{
  return flowSlnPntr;
}

int cgnsAnalyzer::getNVertex()
{
  return nVertex;
}

int cgnsAnalyzer::getNElement()
{
  return nElem;
}

int cgnsAnalyzer::getPhysDim()
{
  return physDim;
}

int cgnsAnalyzer::getElementType()
{
  return sectionType;
}

int cgnsAnalyzer::getNVrtxElem()
{
  return nVrtxElem;
}

double cgnsAnalyzer::getTimeStep()
{
  return timeLabel;
}

CG_MassUnits_t cgnsAnalyzer::getMassUnit()
{
  return massU;
}

CG_LengthUnits_t cgnsAnalyzer::getLengthUnit()
{
  return lengthU;
}

CG_TimeUnits_t cgnsAnalyzer::getTimeUnit()
{
  return timeU;
}

CG_TemperatureUnits_t cgnsAnalyzer::getTemperatureUnit()
{
  return tempU;
}

CG_AngleUnits_t cgnsAnalyzer::getAngleUnit()
{
  return angleU;
}

CG_ZoneType_t cgnsAnalyzer::getZoneType()
{
  return zoneType;
}

bool cgnsAnalyzer::isStructured()
{
  return !isUnstructured;
}

vtkSmartPointer<vtkDataSet> cgnsAnalyzer::getVTKMesh()
{
  if (vtkMesh)
    return vtkMesh;
  else
  {
    exportToVTKMesh();
    return vtkMesh;
  }
}

std::vector<double> cgnsAnalyzer::getVertexCoords()
{
  std::vector<double> crd;
  for (int iVrt = 0; iVrt < nVertex; ++iVrt)
  {
    crd.push_back(xCrd[iVrt]);
    crd.push_back(yCrd[iVrt]);
    crd.push_back(zCrd[iVrt]);
  }
  return crd;
}

std::vector<double> cgnsAnalyzer::getVertexCoords(int vrtxId)
{
  std::vector<double> crd;
  if (vrtxId > nVertex || vrtxId < 0)
  {
    std::cerr << "Requested index is out of bounds." << std::endl;
    return crd;
  }
  crd.push_back(xCrd[vrtxId]);
  crd.push_back(yCrd[vrtxId]);
  crd.push_back(zCrd[vrtxId]);
  return crd;
}

double cgnsAnalyzer::getVrtXCrd(int vrtxId)
{
  return xCrd[vrtxId];
}

std::vector<double> cgnsAnalyzer::getVrtXCrd()
{
  return xCrd;
}

double cgnsAnalyzer::getVrtYCrd(int vrtxId)
{
  return yCrd[vrtxId];
}

std::vector<double> cgnsAnalyzer::getVrtYCrd()
{
  return yCrd;
}

double cgnsAnalyzer::getVrtZCrd(int vrtxId)
{
  return zCrd[vrtxId];
}

std::vector<double> cgnsAnalyzer::getVrtZCrd()
{
  return zCrd;
}

/*
  Returns element connectivity for given element index.
  Input
    elemId : element index (-1 all elements)

  Output
    connectivity of the element.
  
*/
std::vector<int> cgnsAnalyzer::getElementConnectivity(int elemId)
{
  std::vector<int> elmConn;
  int nNdeElm = 0;
  if (elemId > nElem)
  {
    std::cerr << "Element index is out of bounds.\n";
    return elmConn;
  }
  // return the whole connectivity if requested
  if (elemId == -1)
  {
    return elemConn;
  }
  // returning individual element connectivity
  switch (sectionType)
  {
    case CG_TETRA_4:
      nNdeElm = 4;
      break;
    case CG_HEXA_8:
      nNdeElm = 8;
      break;
    case CG_TRI_3:
      nNdeElm = 3;
      elemConn.resize(3 * nElem, -1);
      break;
    case CG_QUAD_4:
      nNdeElm = 4;
      break;
    case CG_TETRA_10:
      nNdeElm = 10;
      break;
    default:
      std::cerr << "Unknown element type " << sectionType << std::endl;
      break;
  }
  for (int iNde = elemId * nNdeElm; iNde < (elemId + 1) * nNdeElm; ++iNde)
    elmConn.push_back(elemConn[iNde]);

  return elmConn;
}


void cgnsAnalyzer::getSectionNames(std::vector<std::string> &names)
{
  // reading section names 
  int eBeg, eEnd, nBdry, parentFlag;
  CG_ElementType_t secTyp;
  char sectionname[33];
  if (cg_nsections(indexFile, indexBase, indexZone, &nSection) != CG_OK)
    std::cerr << "Error in reading sections, " << cg_get_error();
  for (int iSec=1; iSec<= nSection; iSec++)
  {
    if (cg_section_read(indexFile, indexBase, indexZone, iSec, sectionname, 
                &secTyp, &eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
        std::cerr << "Error in load, " << cg_get_error() << std::endl;
    names.push_back(sectionname);
    if (_verb)
      std::cout << "Section " << names[iSec]
              << " Type = " << secTyp
              << " eBeg = " << eBeg
              << " eEnd = " << eEnd
              << " nBdry = " << nBdry
              << std::endl;
  }
}

CG_ElementType_t cgnsAnalyzer::getSectionType(std::string secName)
{
  std::vector<std::string> secNames;
  getSectionNames(secNames);
  auto it = std::find(secNames.begin(), secNames.end(), secName);
  if (it == secNames.end())
  {
      std::cerr << "No section found named " << secName << "\n";
      throw;
  }
  int iSec = it-secNames.begin();

  // reading section info 
  int eBeg, eEnd, nBdry, parentFlag;
  CG_ElementType_t secTyp;
  char sectionname[33];
  if (cg_section_read(indexFile, indexBase, indexZone, iSec, sectionname, 
                &secTyp, &eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
        std::cerr << "Error in load, " << cg_get_error() << std::endl;
  return secTyp;
}

void cgnsAnalyzer::getSectionConn(std::string secName, std::vector<int>& conn, int& nElm)
{
    std::vector<std::string> secNames;
    getSectionNames(secNames);
    auto it = std::find(secNames.begin(), secNames.end(), secName);
    if (it == secNames.end())
    {
        std::cerr << "No section found named " << secName << "\n";
        nElm = 0;
        return;
    }
    // reading connectivities
    int eBeg, eEnd, nBdry, parentFlag, nVrtxElem;
    char sectionname[33];
    CG_ElementType_t secTyp;
    int iSec = it-secNames.begin()+1;
    if (cg_section_read(indexFile, indexBase, indexZone, iSec,
                        sectionname, &secTyp, &eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    switch (secTyp)
    {
      case CG_TETRA_4:
        nVrtxElem = 4;
        break;
      case CG_HEXA_8:
        nVrtxElem = 8;
        break;
      case CG_TRI_3:
        nVrtxElem = 3;
        break;
      case CG_QUAD_4:
        nVrtxElem = 4;
        break;
      case CG_TETRA_10:
        nVrtxElem = 10;
        break;
      default:
        std::cerr << "Unknown element type " << secTyp << std::endl;
        break;
    }
    nElm = eEnd - eBeg + 1;
    conn.resize(nVrtxElem * nElm, -1);
    if (cg_elements_read(indexFile, indexBase, indexZone, iSec, &conn[0], nullptr) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    if (_verb) std::cout << "Size of connectivity vector = " << conn.size() << std::endl;
}


vtkSmartPointer<vtkDataSet> cgnsAnalyzer::getSectionMesh(std::string secName)
{
    // rind information test
    int rdata[2];
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone,
              "GridCoordinates", 0, "end"))
        cg_error_exit();
    if (cg_rind_read(rdata))
        cg_error_exit();
    int nRindNde = rdata[1];
    int one=1;
    int rmax = nVertex + nRindNde;
    // reading all coordinates
    std::vector<double> xCrdR, yCrdR, zCrdR;
    xCrdR.resize(rmax, 0);
    yCrdR.resize(rmax, 0);
    zCrdR.resize(rmax, 0);
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateX", CG_RealDouble, &one, &rmax, &xCrdR[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateY", CG_RealDouble, &one, &rmax, &yCrdR[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    if (cg_coord_read(indexFile, indexBase, indexZone,
                      "CoordinateZ", CG_RealDouble, &one, &rmax, &zCrdR[0]) != CG_OK)
      std::cerr << "Error in load, " << cg_get_error() << std::endl;
    int nElmSec;
    std::vector<int> connSec;
    getSectionConn(secName, connSec, nElmSec);
    CG_ElementType_t secTyp = getSectionType(secName);

    // points to be pushed into dataSet
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp
            = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // allocate size for vtk point container 
    points->SetNumberOfPoints(rmax);
    for (int i = 0; i < rmax; ++i)
      points->SetPoint(i, xCrdR[i], yCrdR[i], zCrdR[i]);

    // add points to vtk mesh data structure 
    dataSet_tmp->SetPoints(points);

    // allocate space for elements
    dataSet_tmp->Allocate(nElmSec);
    // add the elements
    int nNdeElm = connSec.size()/nElmSec;
    for (int i = 0; i < nElmSec; i++)
    {
        vtkSmartPointer<vtkIdList> vtkElmIds = vtkSmartPointer<vtkIdList>::New();
        vtkElmIds->SetNumberOfIds(nNdeElm);
        for (int j = 0; j < nNdeElm; ++j)
          vtkElmIds->SetId(j, connSec[i*nNdeElm+j] - 1);
        switch (secTyp)
        {
          case CG_TETRA_4:
            dataSet_tmp->InsertNextCell(VTK_TETRA, vtkElmIds);
            break;
          case CG_HEXA_8:
            dataSet_tmp->InsertNextCell(VTK_HEXAHEDRON, vtkElmIds);
            break;
          case CG_TRI_3:
            dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkElmIds);
            break;
          case CG_QUAD_4:
            dataSet_tmp->InsertNextCell(VTK_QUAD, vtkElmIds);
            break;
          case CG_TETRA_10:
            dataSet_tmp->InsertNextCell(VTK_QUADRATIC_TETRA, vtkElmIds);
            break;
          default:
            std::cerr << "Unknown element type " << sectionType << std::endl;
            break;
        }
    }
    return(dataSet_tmp);
}


void cgnsAnalyzer::writeSampleStructured()
{
  /*
     dimension statements (note that tri-dimensional arrays
     x,y,z must be dimensioned exactly as [N][17][21] (N>=9)
     for this particular case or else they will be written to
     the CGNS file incorrectly!  Other options are to use 1-D
     arrays, use dynamic memory, or pass index values to a
     subroutine and dimension exactly there):
  */
  double x[9][17][21], y[9][17][21], z[9][17][21];
  cgsize_t isize[3][3];
  int ni, nj, nk, i, j, k;
  int indexFile, icelldim, iphysdim, indexBase;
  int indexZone, indexCoord;
  char basename[33], zonename[33];

  /* create grid points for simple example: */
  ni = 21;
  nj = 17;
  nk = 9;
  for (k = 0; k < nk; ++k)
  {
    for (j = 0; j < nj; ++j)
    {
      for (i = 0; i < ni; ++i)
      {
        x[k][j][i] = i;
        y[k][j][i] = j;
        z[k][j][i] = k;
      }
    }
  }
  printf("\ncreated simple 3-D grid points");

  /* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
  /* open CGNS file for write */
  if (cg_open(cgFileName.c_str(), CG_MODE_WRITE, &indexFile)) cg_error_exit();
  /* create base (user can give any name) */
  strcpy(basename, "Base");
  icelldim = 3;
  iphysdim = 3;
  cg_base_write(indexFile, basename, icelldim, iphysdim, &indexBase);
  /* define zone name (user can give any name) */
  strcpy(zonename, "Zone  1");
  /* vertex size */
  isize[0][0] = 21;
  isize[0][1] = 17;
  isize[0][2] = 9;
  /* cell size */
  isize[1][0] = isize[0][0] - 1;
  isize[1][1] = isize[0][1] - 1;
  isize[1][2] = isize[0][2] - 1;
  /* boundary vertex size (always zero for structured grids) */
  isize[2][0] = 0;
  isize[2][1] = 0;
  isize[2][2] = 0;
  /* create zone */
  cg_zone_write(indexFile, indexBase, zonename, *isize, CG_Structured, &indexZone);
  /* write grid coordinates (user must use SIDS-standard names here) */
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateX",
                 x, &indexCoord);
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateY",
                 y, &indexCoord);
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateZ",
                 z, &indexCoord);
}

void cgnsAnalyzer::writeSampleUnstructured()
{
  /*
     dimension statements (note that tri-dimensional arrays
     x,y,z must be dimensioned exactly as [N][17][21] (N>=9)
     for this particular case or else they will be written to
     the CGNS file incorrectly!  Other options are to use 1-D
     arrays, use dynamic memory, or pass index values to a
     subroutine and dimension exactly there):
  */
  double x[9][17][21], y[9][17][21], z[9][17][21];
  cgsize_t isize[3][3];
  int ni, nj, nk, i, j, k;
  int indexFile, icelldim, iphysdim, indexBase;
  int indexZone, indexCoord, indexSection;
  int index_flow, index_field;
  char basename[33], zonename[33];

  /* create gridpoints for simple example: */
  ni = 21;
  nj = 17;
  nk = 9;
  for (k = 0; k < nk; ++k)
  {
    for (j = 0; j < nj; ++j)
    {
      for (i = 0; i < ni; ++i)
      {
        x[k][j][i] = i;
        y[k][j][i] = j;
        z[k][j][i] = k;
      }
    }
  }
  printf("\nCreated simple 3-D grid points");

  /* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
  /* open CGNS file for write */
  if (cg_open(cgFileName.c_str(), CG_MODE_WRITE, &indexFile)) cg_error_exit();
  /* create base (user can give any name) */
  strcpy(basename, "Base");
  icelldim = 3;
  iphysdim = 3;
  cg_base_write(indexFile, basename, icelldim, iphysdim, &indexBase);
  /* define zone name (user can give any name) */
  strcpy(zonename, "Zone  1");
  /* vertex size */
  isize[0][0] = 21 * 17 * 9;
  isize[0][1] = 20 * 16 * 8;
  isize[0][2] = 0;
  /* create zone */
  cg_zone_write(indexFile, indexBase, zonename, *isize, CG_Unstructured, &indexZone);
  /* write grid coordinates (user must use SIDS-standard names here) */
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateX",
                 x, &indexCoord);
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateY",
                 y, &indexCoord);
  cg_coord_write(indexFile, indexBase, indexZone, CG_RealDouble, "CoordinateZ",
                 z, &indexCoord);
  /* write connectivities */
  cgsize_t elemConn[8][20 * 16 * 8];
  int iElem = -1;
  /* index of first element */
  int nElemStart = 1;
  for (k = 1; k < nk; ++k)
  {
    for (j = 1; j < nj; ++j)
    {
      for (i = 1; i < ni; ++i)
      {
        iElem++;
        int ifirstNode = i + (j - 1) * ni + (k - 1) * ni * nj;
        elemConn[0][iElem] = ifirstNode;
        elemConn[1][iElem] = ifirstNode + 1;
        elemConn[2][iElem] = ifirstNode + 1 + ni;
        elemConn[3][iElem] = ifirstNode + ni;
        elemConn[4][iElem] = ifirstNode + ni * nj;
        elemConn[5][iElem] = ifirstNode + ni * nj + 1;
        elemConn[6][iElem] = ifirstNode + ni * nj + 1 + ni;
        elemConn[7][iElem] = ifirstNode + ni * nj + ni;
      }
    }
  }
  /* index of last element */
  int nElemEnd = iElem;
  /* unsorted boundary elements */
  int nBdyElem = 0;
  /* write HEX_8 element connectivity (user can give any name) */
  cg_section_write(indexFile, indexBase, indexZone,
                   "Elem", CG_HEXA_8, nElemStart, nElemEnd,
                   nBdyElem, elemConn[0], &indexSection);
  /* write vertex based field data */
  std::string solName = "solution";
  cg_sol_write(indexFile, indexBase, indexZone,
               solName.c_str(), CG_Vertex, &index_flow);
  std::vector<double> T;
  T.resize(isize[0][0], 10);
  cg_field_write(indexFile, indexBase, indexZone,
                 index_flow, CG_RealDouble, "T", &T[0], &index_field);
  /* close CGNS file */
  cg_close(indexFile);
}

void cgnsAnalyzer::clearAllSolutionData()
{
  // clearing all data related maps and fields
  nSolution = 0;
  nField = 0;
  solutionNameLocMap.clear();
  solutionName.clear();
  solutionGridLocation.clear();
  solutionMap.clear();
  appendedSolutionName.clear();
  // clearing all solution data objects
  for (auto &it : slnDataCont)
    delete it;
  slnDataCont.clear();
  solutionDataPopulated = false;
}

// lists all solution based data
void cgnsAnalyzer::populateSolutionDataNames()
{
  // only one time needed for each zone
  if (solutionDataPopulated)
    return;
  // populating solution data
  char fieldName[33];
  char solName[33];
  CG_GridLocation_t gloc;
  int nFlds;
  int cntr = 0;
  CG_DataType_t dt;
  // number of solution fields
  cg_nsols(indexFile, indexBase, indexZone, &nSolution);
  for (int iSol = 1; iSol <= nSolution; ++iSol)
  {
    std::pair<int, keyValueList> slnPair;
    keyValueList fldIndxSln;
    cg_sol_info(indexFile, indexBase, indexZone, iSol,
                solName, &gloc);
    cg_nfields(indexFile, indexBase, indexZone, iSol, &nFlds);
    solutionNameLocMap[solName] = gloc;
    for (int iFld = 1; iFld <= nFlds; ++iFld)
    {
      // name and location
      solutionName.push_back(solName);
      solutionGridLocation.push_back(gloc);
      // field information
      cg_field_info(indexFile, indexBase, indexZone,
                    iSol, iFld, &dt, fieldName);
      fldIndxSln[iFld] = fieldName;
    }
    slnPair.first = iSol;
    slnPair.second = fldIndxSln;
    solutionMap[cntr++] = slnPair;
  }
  solutionDataPopulated = true;
}

void cgnsAnalyzer::getSolutionDataNames(std::vector<std::string>& list)
{
  populateSolutionDataNames();
  for (auto &it : solutionMap)
  {
    std::pair<int, keyValueList> slnIndxListPair = it.second;
    keyValueList fldIndxSln = slnIndxListPair.second;
    for (auto &it2 : fldIndxSln)
    {
      list.push_back(it2.second);
    }
  }
}

/*
   Given solution name provides solution vector by reading from the CGNS file
   and processing it. Also passes the type of solution in the return. The 
   solution type (2) for CG_Vertex based and (3) for the element based. Types
   are in agreement with CG_GridLocation_t definition provided by CGNS API.
   NOTE: It is assumed that vector-valued nodal and element data are 
   decomposed into scalar fields.
*/
solution_type_t cgnsAnalyzer::getSolutionData(std::string sName, std::vector<double>& slnData)
{
  CG_DataType_t dt;
  char fieldName[33];
  int slnCntr = -1;
  int slnIndx = -1;
  int fldIndx = -1;
  int dataType = -1;
  // find the solution index
  for (auto &it : solutionMap)
  {
    std::pair<int, keyValueList> slnIndxListPair = it.second;
    keyValueList fldIndxSln = slnIndxListPair.second;
    for (auto &it2 : fldIndxSln)
    {
      slnCntr++;
      if (!((it2.second).compare(sName)))
      {
        fldIndx = it2.first;
        slnIndx = slnIndxListPair.first;
        break;
      }
    }
    if (slnIndx != -1)
      break;
  }
  // fail check
  if (slnIndx == -1)
  {
    std::cerr << "The solution name "
              << sName << " does not exist.\n";
    return UNKNOWN;
  }
  // check if data exists
  if (cg_field_info(indexFile, indexBase, indexZone, slnIndx, fldIndx,
                    &dt, fieldName) != CG_OK)
    std::cerr << "Error in reading solution, " << cg_get_error() << std::endl;
  // reading actual data from the file
  int one = 1;
  dataType = solutionGridLocation[slnCntr];
  if (isUnstructured)
  {
    if (dataType == CG_Vertex)
    {
      slnData.resize(nVertex, -1.0);
      rmax[0] = nVertex;
    }
    else if (dataType == CG_CellCenter)
    {
      slnData.resize(nElem, -1.0);
      rmax[0] = nElem;
    }
    else
    {
      std::cerr << "Unknown data gird location " << solutionGridLocation[slnCntr]
                << std::endl;
      return UNKNOWN;
    }
  }
  else
  {
    if (dataType == CG_Vertex)
    {
      slnData.resize(nVertex, -1.0);
      rmax[0] = cgCoreSize[0]; 
      rmax[1] = cgCoreSize[1]; 
      rmax[2] = cgCoreSize[2]; 
    }
    else if (dataType == CG_CellCenter)
    {
      slnData.resize(nElem, -1.0);
      rmax[0] = cgCoreSize[3]; 
      rmax[1] = cgCoreSize[4]; 
      rmax[2] = cgCoreSize[5]; 
    }
    else
    {
      std::cerr << "Unknown data gird location " << solutionGridLocation[slnCntr]
                << std::endl;
      return UNKNOWN;
    }
  
  }

  if (isUnstructured)
  {
    if (dt == CG_RealDouble)
    {
      // for double solution data
      if (cg_field_read(indexFile, indexBase, indexZone, slnIndx, fieldName,
                        dt, &one, &rmax[0], &slnData[0]) != CG_OK)
        std::cerr << "Error in reading solution data, " << cg_get_error() << std::endl;
    }
    else if (dt == CG_Integer)
    {
      // for integer solution data
      std::vector<int> tmpSlnData;
      tmpSlnData.resize(rmax[0], -1);
      if (cg_field_read(indexFile, indexBase, indexZone, slnIndx, fieldName,
                        dt, &one, &rmax[0], &tmpSlnData[0]) != CG_OK)
        std::cerr << "Error in reading solution data, " << cg_get_error() << std::endl;
      slnData.clear();
      for (int &it : tmpSlnData)
        slnData.push_back(it);
    }
  }
  else
  {
    if (cg_field_read(indexFile, indexBase, indexZone, slnIndx, fieldName,
                      dt, &rmin[0], &rmax[0], &slnData[0]) != CG_OK)
      std::cerr << "Error in reading solution data, " << cg_get_error() << std::endl;
  }

  // returns the type of the data
  return (dataType == 2 ? NODAL : ELEMENTAL);
}

/*
   Returns pointer to the solutionData class containing solution
   information with the given name.
*/
solutionData *cgnsAnalyzer::getSolutionDataObj(std::string sName)
{
  // return pointer if already loaded
  if (!slnDataCont.empty())
  {
    for (auto &is : slnDataCont)
    {
      if (!strcmp((is->getDataName()).c_str(), sName.c_str()))
      {
        return is;
      }
    }
  }
  // load if needed
  std::vector<double> slnData;
  solution_type_t dataType = getSolutionData(sName, slnData);
  if (dataType == UNKNOWN)
  {
    std::cerr << "Unknown data type is not supported.!\n";
    return nullptr;
  }
  // allocating new data object
  solutionData *slnDataObjPtr = new solutionData(sName, dataType);
  if (dataType == NODAL)
  {
    slnDataObjPtr->appendData(slnData, slnData.size(), 1);
  }
  else
  {
    slnDataObjPtr->appendData(slnData, slnData.size(), 1);
  }
  return slnDataObjPtr;
}

int cgnsAnalyzer::getNVertexSolution()
{
  int nVrtData = 0;

  // loading data if not yet
  if (slnDataCont.empty())
    loadSolutionDataContainer();

  // number of vertex solution data
  for (auto &is : slnDataCont)
    if (is->getDataType() == NODAL)
      nVrtData++;

  return nVrtData;
}

int cgnsAnalyzer::getNCellSolution()
{
  int nCellData = 0;

  // loading data if not yet
  if (slnDataCont.empty())
    loadSolutionDataContainer();

  // number of vertex solution data
  for (auto &is : slnDataCont)
    if (is->getDataType() == ELEMENTAL)
      nCellData++;

  return nCellData;
}

solution_type_t cgnsAnalyzer::getSolutionDataStitched(std::string sName, std::vector<double>& slnData,
                                                      int& outNData, int& outNDim)
{
  // load data if needed
  if (slnDataCont.empty())
    loadSolutionDataContainer();

  for (auto &is : slnDataCont)
  {
    if (!strcmp((is->getDataName()).c_str(), sName.c_str()))
    {
      is->getData(slnData, outNData, outNDim);
      return (is->getDataType());
    }
  }
  return UNKNOWN;
}

void cgnsAnalyzer::appendSolutionData(std::string sName, std::vector<double>& slnData,
                                      solution_type_t dt, int inNData, int inNDim)
{
  // load data if needed
  if (slnDataCont.empty())
    loadSolutionDataContainer();

  solutionData* nwSlnPtr = new solutionData(sName, dt);
  nwSlnPtr->appendData(slnData, inNData, inNDim);
  slnDataCont.push_back(nwSlnPtr);
  appendedSolutionName.push_back(sName);
}

void cgnsAnalyzer::appendSolutionData(std::string sName, double slnData,
                                      solution_type_t dt, int inNData, int inNDim)
{
  // load data if needed
  if (slnDataCont.empty())
    loadSolutionDataContainer();
  solutionData* nwSlnPtr = new solutionData(sName, dt);
  // convert to vector
  std::vector<double> slnDataVec;
  slnDataVec.resize(inNData, slnData);
  nwSlnPtr->appendData(slnDataVec, inNData, inNDim);
  // register
  slnDataCont.push_back(nwSlnPtr);
  appendedSolutionName.push_back(sName);
}

bool cgnsAnalyzer::delAppSlnData(std::string sName)
{
  for (auto is = appendedSolutionName.begin(); is != appendedSolutionName.end(); ++is)
    if (strcmp(sName.c_str(), (*is).c_str()) == 0)
    {
      delete getSolutionDataObj(sName);
      appendedSolutionName.erase(is);
      return true;
    }
  return false;
}

void cgnsAnalyzer::getAppendedSolutionDataName(std::vector<std::string>& appSName)
{
  appSName.insert(appSName.end(), appendedSolutionName.begin(), appendedSolutionName.end());
}

/* provides node names in which solution data are stored in */
std::vector<std::string> cgnsAnalyzer::getSolutionNodeNames()
{
  return solutionName;
}

/* provides the list of solution types */
std::vector<CG_GridLocation_t> cgnsAnalyzer::getSolutionGridLocations()
{
  return solutionGridLocation;
}

/* provides solution map */
std::map<int, std::pair<int, keyValueList> > cgnsAnalyzer::getSolutionMap()
{
  return solutionMap;
}

/* gets a map between solution node name and solution location type */
std::map<std::string, CG_GridLocation_t> cgnsAnalyzer::getSolutionNameLocMap()
{
  return solutionNameLocMap;
}

void cgnsAnalyzer::exportToVTKMesh()
{
  if (!vtkMesh)
  {
    // remove rind data if any
    cleanRind();
    // points to be pushed into dataSet
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp
            = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // allocate size for vtk point container 
    points->SetNumberOfPoints(nVertex);
    for (int i = 0; i < nVertex; ++i)
    {
      points->SetPoint(i, xCrd[i], yCrd[i], zCrd[i]);
    }
    // add points to vtk mesh data structure 
    dataSet_tmp->SetPoints(points);
    // allocate space for elements
    dataSet_tmp->Allocate(nElem);
    // add the elements
    for (int i = 0; i < nElem; ++i)
    {
      vtkSmartPointer<vtkIdList> vtkElmIds = vtkSmartPointer<vtkIdList>::New();
      std::vector<int> cgnsElmIds(getElementConnectivity(i));
      vtkElmIds->SetNumberOfIds(cgnsElmIds.size());
      for (int j = 0; j < cgnsElmIds.size(); ++j)
      {
        vtkElmIds->SetId(j, cgnsElmIds[j] - 1);
      }
      switch (sectionType)
      {
        case CG_TETRA_4:
          dataSet_tmp->InsertNextCell(VTK_TETRA, vtkElmIds);
          break;
        case CG_HEXA_8:
          dataSet_tmp->InsertNextCell(VTK_HEXAHEDRON, vtkElmIds);
          break;
        case CG_TRI_3:
          dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkElmIds);
          break;
        case CG_QUAD_4:
          dataSet_tmp->InsertNextCell(VTK_QUAD, vtkElmIds);
          break;
        case CG_TETRA_10:
          dataSet_tmp->InsertNextCell(VTK_QUADRATIC_TETRA, vtkElmIds);
          break;
        default:
          std::cerr << "Unknown element type " << sectionType << std::endl;
          break;
      }
    }
    vtkMesh = dataSet_tmp;
  }
}

void cgnsAnalyzer::overwriteSolData(meshBase* mbObj)
{
  // loading data if not yet
  if (slnDataCont.empty())
    loadSolutionDataContainer();
  // write individual data fields
  int iSol = -1;
  for (auto &is : solutionMap)
  {
    std::pair<int, keyValueList> slnPair = is.second;
    int slnIdx = slnPair.first;
    keyValueList fldLst = slnPair.second;
    for (auto &ifl : fldLst)
    {
      iSol++;
      std::vector<double> newData;
      if (solutionGridLocation[iSol] == CG_Vertex)
      {
        mbObj->getPointDataArray(ifl.second, newData);
      }
      else
      {
        //gs field is weird in irocstar files we don't write it back
        if (/*!(ifl->second).compare("gs") ||*/ !(ifl.second).compare("mdot_old"))
        {
          continue;
        }
        mbObj->getCellDataArray(ifl.second, newData);
      }
      std::cout << "Writing "
                << newData.size()
                << " to "
                << ifl.second
                << " located in "
                << solutionName[iSol]
                << std::endl;
      // write to file
      //if (!(ifl.second).compare("bflag")) // cg_io complains if this isn't CG_Integer type
      if (slnDataCont[slnIdx]->getDataType() == CG_Integer || !(ifl.second).compare("bflag"))
      {
        // need to cast to int before passing as void*
        std::vector<int> newDataToInt(newData.begin(), newData.end());
        overwriteSolData(ifl.second, solutionName[iSol], slnIdx, CG_Integer, &newDataToInt[0]);
      }
      else
        overwriteSolData(ifl.second, solutionName[iSol], slnIdx, CG_RealDouble, &newData[0]);
    }
  }
}

void cgnsAnalyzer::overwriteSolData(const std::string& fname,
                                    const std::string& ndeName,
                                    int slnIdx, CG_DataType_t dt, void* data)
{
  // write solution to file
  CG_GridLocation_t gloc(solutionNameLocMap[fname]);
  int fldIdx;
  if (cg_field_write(indexFile, indexBase, indexZone, slnIdx,
                     dt, fname.c_str(), data, &fldIdx))
    cg_error_exit();
  // finding range of data
  double* tmpData = (double*) data;
  double min = tmpData[0];
  double max = tmpData[0];
  int nItr = 0;
  if (gloc == CG_Vertex)
  {
    nItr = nVertex;
  }
  else
  {
    nItr = nElem;
  }
  for (int it = 0; it < nItr; ++it)
  {
    min = std::min(tmpData[it], min);
    max = std::max(tmpData[it], max);
  }
  // writing range descriptor
  std::ostringstream os;
  os << min << ", " << max;
  std::string range = os.str();
  if (cg_goto(indexFile, indexBase, "Zone_t", indexZone,
              "FlowSolution_t", slnIdx, "DataArray_t", fldIdx, "end"))
    cg_error_exit();
  if (cg_descriptor_write("Range", range.c_str())) cg_error_exit();
  // write DimensionalExponents and units for cell data
  if (gloc == CG_CellCenter)
  {
    if (cg_goto(indexFile, indexBase, "Zone_t", indexZone, "FlowSolution_t", slnIdx,
                "DataArray_t", fldIdx, "end"))
      cg_error_exit();
    // dummy exponents and units
    float exponents[5] = {0, 0, 0, 0, 0};
    if (cg_exponents_write(CG_RealSingle, exponents)) cg_error_exit();
    if (cg_descriptor_write("Units", "dmy")) cg_error_exit();
  }
}

void cgnsAnalyzer::exportToMAdMesh(const MAd::pMesh MAdMesh)
{
  // --- Build the vertices ---
  MAdToCgnsIds.clear();
  cgnsToMAdIds.clear();
  for (int iV = 0; iV < nVertex; ++iV)
  {
    MAdMesh->add_point(iV + 1, xCrd[iV], yCrd[iV], zCrd[iV]);
    cgnsToMAdIds[iV] = iV + 1;
    MAdToCgnsIds[iV + 1] = iV;
  }

  // --- Build the elements ---
  if (physDim == 3)
    switch (sectionType)
    {
      case CG_TETRA_4:
      {
        for (int iC = 0; iC < nElem; ++iC)
        {
          int iN = iC * 4;
          // the last argument for the next function is just a dummy now
          MAd::pGEntity geom = (MAd::pGEntity) MAd::GM_regionByTag(MAdMesh->model, 0);
          //geom->setPhysical(3,1);
          MAd::MDB_Tet *tet = MAdMesh->add_tet(elemConn[iN], elemConn[iN + 1],
                                               elemConn[iN + 2], elemConn[iN + 3], geom);
          // changing tet ID
          tet->iD = iC + 1;
        }
      }
        break;
      case CG_TRI_3:
      {
        for (int iC = 0; iC < nElem; ++iC)
        {
          int iN = iC * 3;
          // the last argument for the next function is just a dummy now
          MAd::pGEntity geom = (MAd::pGEntity) MAd::GM_faceByTag(MAdMesh->model, 0);
          MAdMesh->add_triangle(elemConn[iN], elemConn[iN + 1],
                                elemConn[iN + 2], geom);
        }
      }
        break;
      default:
        std::cerr << "Current version only works for TRI and TET elements. Element type "
                  << sectionType << " is not supported." << std::endl;
        break;
    }
}

void cgnsAnalyzer::classifyMAdMeshOpt(const MAd::pMesh MAdMesh)
{
  /*
    Here, the entities of the MAdLib mesh should be classified
    on their corresponding geometrical entities, like for boundary 
    faces in 3D for instance. The implementation of this step 
    is highly dependent on the implementation of Solver_mesh and 
    Solver_model so it is up to the reader to add the right 
    instructions here. 

    Note that the geometrical entities have been created in the 
    execution of 'exportToMAdModel'. Any mesh entity can be 
    associated to a geometrical entity using the EN_setWhatIn(...) 
    function of the MAdLib mesh interface.

    Note that all the steps involving geometrical entities can be
    replaced by appropriate constraints on boundary mesh entities
    (see AdaptInterface.h) but no mesh modification will therefore
    be applied on the boundaries, which can be problematic for some
    computations.
  */
  // finding boundary faces and classifying them as 2 dimensional
  //MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(MAdMesh->model, 0);
  //MAdMesh->classify_grid_boundaries(bnd);
  // classifying the rest of the domain (faces, edges, vertices)
  // the element geometric region properties will be used by them.
  MAdMesh->classify_unclassified_entities();
  MAdMesh->destroyStandAloneEntities();
}

void cgnsAnalyzer::classifyMAdMeshBnd(const MAd::pMesh MAdMesh)
{
  // finding boundary faces and classifying them as 2 dimensional
  MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(MAdMesh->model, 0);
  MAdMesh->classify_grid_boundaries(bnd);
}

void cgnsAnalyzer::unclassifyMAdMeshBnd(const MAd::pMesh MAdMesh)
{
  // finding boundary faces and classifying them as 2 dimensional
  MAdMesh->unclassify_grid_boundaries();
}

void cgnsAnalyzer::buildVertexKDTree()
{
  //ANNpointArray vrtxCrd;
  if (vrtxCrd)
    annDeallocPts(vrtxCrd);
  // clearing onld instance
  vrtxCrd = annAllocPts(nVertex, physDim);
  if (kdTree)
  {
    delete kdTree;
  }

  // filling up vertex coordinate array for the current mesh
  for (int iVrt = 0; iVrt < nVertex; ++iVrt)
  {
    vrtxCrd[iVrt][0] = xCrd[iVrt];
    vrtxCrd[iVrt][1] = yCrd[iVrt];
    vrtxCrd[iVrt][2] = zCrd[iVrt];
  }
  // building kdTree
  kdTree = new ANNkd_tree(vrtxCrd, nVertex, physDim);
}

void cgnsAnalyzer::buildElementKDTree()
{
  //ANNpointArray vrtxIdx;
  if (vrtxIdx)
    annDeallocPts(vrtxIdx);
  // clearing onld instance
  vrtxIdx = annAllocPts(nElem, nVrtxElem);
  if (kdTreeElem)
    delete kdTreeElem;

  // filling up vertex coordinate array for the current mesh
  for (int iElem = 0; iElem < nElem; ++iElem)
  {
    for (int iVrtx = 1; iVrtx <= nVrtxElem; ++iVrtx)
    {
      vrtxIdx[iElem][iVrtx] = elemConn[(iElem - 1) * nVrtxElem + iVrtx];
    }
  }
  // building kdTree
  kdTreeElem = new ANNkd_tree(vrtxIdx, nElem, nVrtxElem);
}

/*
   Check for duplicated vertices in the grid.
*/
void cgnsAnalyzer::checkVertex()
{
  // (re)building the kdTree
  buildVertexKDTree();
  // loop through vertices to find duplicated ones
  int nDupVer = 0;
  for (int iVrt = 0; iVrt < getNVertex(); ++iVrt)
  {
    ANNpoint qryVrtx;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = getVrtXCrd(iVrt);
    qryVrtx[1] = getVrtYCrd(iVrt);
    qryVrtx[2] = getVrtZCrd(iVrt);
    nnIdx = new ANNidx[1];
    dists = new ANNdist[1];
    kdTree->annkSearch(qryVrtx, 2, nnIdx, dists);
    if (dists[1] < 1e-10)
    {
      nDupVer++;
      std::cout << "Vertex " << iVrt << " is duplicated."
                << " Distances = " << dists[0]
                << " " << dists[1] << std::endl;
    }
    delete[] nnIdx;
    delete[] dists;
    annDeallocPt(qryVrtx);
  }
  std::cout << "Found " << nDupVer << " duplicate vertex.\n";
}

/*
   Check element to element connectivity making sure there is at least
   given number of nodes shared between neighboring elements.
*/
bool cgnsAnalyzer::checkElmConn(int nSharedNde)
{
  /*
  MatrixInt eConn(nElem, nVrtxElem);
  MatrixInt dummy(nElem, nElem);
  VectorInt elmIdx(nElem);
  for (int iElm = 0; iElm < nElem; ++iElm)
  {
    elmIdx(iElm) = iElm;
    for (int iNde = 0; iNde < nVrtxElem; ++iNde)
      eConn(iElm, iNde) = elemConn[iElm * nVrtxElem + iNde];
  }
  clock_t startTime = clock();
  dummy = eConn * eConn.transpose();
  std::cout << double(clock() - startTime) / (double) CLOCKS_PER_SEC << " seconds." << std::endl;
  */

  // building node to element map
  std::map<int, std::set<int> > nde2Elm;
  for (int iElm = 0; iElm < nElem; ++iElm)
    for (int iNde = 0; iNde < nVrtxElem; ++iNde)
      nde2Elm[elemConn[iElm * nVrtxElem + iNde]].insert(iElm);

  /*
  std::cout << "nde2Elm = " << std::endl;
  for (auto it = nde2Elm.begin(); it != nde2Elm.end(); ++it)
  {
    std::set<int> tmp = it->second;
    for (auto it2 = tmp.begin(); it2 != tmp.end(); ++it2)
      std::cout << *it2 << " ";
    std::cout << std::endl;
  }
  */

  // going through the list of elements and finding 
  // those with less than nSharedNde 
  std::vector<int> hangingElmIdx;
  for (int iElm = 0; iElm < nElem; ++iElm)
  {
    int nShrdNde = 0;
    for (int iNde = 0; iNde < nVrtxElem; ++iNde)
      nShrdNde += nde2Elm[iNde].size() > 1;
    if (nShrdNde < nSharedNde)
      hangingElmIdx.push_back(iElm);
  }

  std::cout << "Number of elements with less than " << nSharedNde
            << " shared node is " << hangingElmIdx.size() << std::endl;
  return !hangingElmIdx.empty();
}

/* 
   Gets element center coordinates.
*/
std::vector<double> cgnsAnalyzer::getElmCntCoords(MAd::pMesh msh)
{
  std::vector<double> elmCntCrds;
  MAd::RIter ri = M_regionIter(msh);
  int rCnt = 0;
  while (MAd::pRegion pr = RIter_next(ri))
  {
    double xc[3];
    MAd::R_center(pr, xc);
    elmCntCrds.push_back(xc[0]);
    elmCntCrds.push_back(xc[1]);
    elmCntCrds.push_back(xc[2]);
  }
  return elmCntCrds;
}

/*
   Stitches the mesh of the given instance to the current one.
*/
void cgnsAnalyzer::stitchMesh(cgnsAnalyzer* inCg, bool withFields)
{
  // sanity check
  if (physDim != inCg->getPhysDim())
  {
    std::cerr << "Error : Stitching mesh should have same dimensions."
              << std::endl;
    return;
  }

  // removing rind data (if any)
  cleanRind();
  inCg->cleanRind();

  // (re)building the kdTree
  buildVertexKDTree();

  // clear old masks 
  vrtDataMask.clear();
  elmDataMask.clear();

  // adding new mesh non-repeating vertices
  std::vector<int> newVrtIdx;
  std::vector<int> rptVrtIdx;
  std::map<int, int> rptVrtMap; // <newMeshIdx, currentMeshIdx>
  std::vector<double> newXCrd;
  std::vector<double> newYCrd;
  std::vector<double> newZCrd;
  int nNewVrt = 0;
  for (int iVrt = 0; iVrt < inCg->getNVertex(); ++iVrt)
  {
    ANNpoint qryVrtx;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    qryVrtx = annAllocPt(physDim);
    qryVrtx[0] = inCg->getVrtXCrd(iVrt);
    qryVrtx[1] = inCg->getVrtYCrd(iVrt);
    qryVrtx[2] = inCg->getVrtZCrd(iVrt);
    nnIdx = new ANNidx[1];
    dists = new ANNdist[1];
    kdTree->annkSearch(qryVrtx, 1, nnIdx, dists);
    if (dists[0] > searchEps)
    {
      nNewVrt++;
      vrtDataMask.push_back(true);
      newVrtIdx.push_back(nVertex + nNewVrt);
      newXCrd.push_back(qryVrtx[0]);
      newYCrd.push_back(qryVrtx[1]);
      newZCrd.push_back(qryVrtx[2]);
    }
    else
    {
      vrtDataMask.push_back(false);
      newVrtIdx.push_back(nnIdx[0] + 1);
      rptVrtIdx.push_back(iVrt);
      rptVrtMap[iVrt] = nnIdx[0] + 1;
    }
    delete[] nnIdx;
    delete[] dists;
    annDeallocPt(qryVrtx);
  }
  std::cout << "Found " << nNewVrt << " new vertices.\n";
  std::cout << "Number of repeating index " << rptVrtIdx.size()
            << std::endl;

  // currently implemented to add all new elements
  std::vector<int> newElemConn;
  int nNewElem = 0;
  for (int iElem = 0; iElem < inCg->getNElement(); ++iElem)
  {
    std::vector<int> rmtElemConn = inCg->getElementConnectivity(iElem);

    // just adding all elements
    elmDataMask.push_back(true);
    nNewElem++;
    newElemConn.insert(newElemConn.end(),
                       rmtElemConn.begin(), rmtElemConn.end());
  }
  std::cout << "Found " << nNewElem << " new elements.\n";

  // switching connectivity table to global
  for (int iIdx = 0; iIdx < newElemConn.size(); ++iIdx)
  {
    newElemConn[iIdx] = newVrtIdx[newElemConn[iIdx] - 1];
  }

  // stitching field values if requested
  if (withFields)
    stitchFields(inCg);

  // updating internal data structure
  nVertex += nNewVrt;
  nElem += nNewElem;
  xCrd.insert(xCrd.end(), newXCrd.begin(), newXCrd.end());
  yCrd.insert(yCrd.end(), newYCrd.begin(), newYCrd.end());
  zCrd.insert(zCrd.end(), newZCrd.begin(), newZCrd.end());
  elemConn.insert(elemConn.end(), newElemConn.begin(), newElemConn.end());
  zoneNames.push_back(inCg->getZoneName());
}

void cgnsAnalyzer::loadSolutionDataContainer(int verb)
{
  // look at the solution names of current grid
  std::vector<std::string> crntCgList;
  getSolutionDataNames(crntCgList);
  if (verb > 0)
    for (auto &it : crntCgList)
      std::cout << it << std::endl;
  // load solution data container if empty
  if (slnDataCont.empty())
  {
    for (auto &is : crntCgList)
    {
      solutionData* slnDataPtr = getSolutionDataObj(is);
      if (verb > 0)
        std::cout << is
                  << " number of data read "
                  << slnDataPtr->getNData()
                  << " "
                  << slnDataPtr->getNDim()
                  << std::endl;
      slnDataCont.push_back(slnDataPtr);
    }
  }
  // information
  std::cout << "  Number of Solution Field = " << slnDataCont.size()
            << std::endl;
}

void cgnsAnalyzer::stitchFields(cgnsAnalyzer* inCg)
{
  // look at the solution names of current grid
  std::vector<std::string> crntCgList;
  getSolutionDataNames(crntCgList);

  // load solution data container if empty
  if (slnDataCont.empty())
    loadSolutionDataContainer();

  // now going through the list of inCg and stitch only repeating data
  std::vector<std::string> inCgList;
  inCg->getSolutionDataNames(inCgList);
  for (auto &id : inCgList)
  {
    bool isRptData = false;
    for (auto &is : crntCgList)
      if (!strcmp(id.c_str(), is.c_str()))
      {
        isRptData = true;
        std::cout << id << " is an existing field.\n";
        break;
      }

    if (isRptData)
    {
      std::vector<double> inCgSlnData;
      int outNData, outNDim;
      solutionData* inCgSlnDataPtr = inCg->getSolutionDataObj(id);
      if (inCgSlnDataPtr->getDataType() == NODAL)
      {
        inCgSlnDataPtr->getData(inCgSlnData, outNData, outNDim, vrtDataMask);
        //testing
        //std::cout << " will append " << outNData << " data points.\n";
      }
      else
      {
        inCgSlnDataPtr->getData(inCgSlnData, outNData, outNDim, elmDataMask);
        //testing
        //std::cout << " will append " << outNData << " data points.\n";
      }
      // attaching data
      int nData = slnDataCont.size();
      for (int icd=0; icd<nData; icd++)
      {
        if (!strcmp((slnDataCont[icd]->getDataName()).c_str(), id.c_str()))
        {
          std::cout << "To -> " << slnDataCont[icd]->getDataName() << std::endl;
          //std::cout << "Current size = " << icnt->getNData() << std::endl;
          slnDataCont[icd]->appendData(inCgSlnData, inCgSlnData.size(), 1);
          //std::cout << "Size after = " << icnt->getNData() << std::endl;
        }
      }
    }
  }

  // now go through appended data to the current grid and see if they
  // also existing in inCg and thus stitch them as well.
  std::vector<std::string> inCgAppLst;
  inCg->getAppendedSolutionDataName(inCgAppLst);
  if (appendedSolutionName.empty() || inCgAppLst.empty())
    return;
  for (auto &id : inCgAppLst)
  {
    bool isRptData = false;
    //std::cout << "Remote CG field name = " << id.c_str() << std::endl;
    for (auto &is : appendedSolutionName)
    {
      //std::cout << "Current CG field name = " << is.c_str() << std::endl;
      //std::cout << " strcmp = " << strcmp(id.c_str(), is.c_str()) << std::endl;
      if (strcmp(id.c_str(), is.c_str()) == 0)
      {
        isRptData = true;
        //std::cout << id << " is an existing field.\n";
        break;
      }
    }

    if (isRptData)
    {
      std::vector<double> inCgSlnData;
      int outNData, outNDim;
      solutionData* inCgSlnDataPtr = inCg->getSolutionDataObj(id);
      if (inCgSlnDataPtr->getDataType() == NODAL)
      {
        inCgSlnDataPtr->getData(inCgSlnData, outNData, outNDim, vrtDataMask);
        //testing
        //std::cout << " will append " << outNData << " data points.\n";
      }
      else
      {
        inCgSlnDataPtr->getData(inCgSlnData, outNData, outNDim, elmDataMask);
        //testing
        std::cout << " will append " << outNData << " data points.\n";
      }
      // attaching data
      for (auto &icnt : slnDataCont)
      {
        if (!strcmp((icnt->getDataName()).c_str(), id.c_str()))
        {
          //std::cout << "To -> " << icnt->getDataName() << std::endl;
          //std::cout << "Current size = " << icnt->getNData() << std::endl;
          icnt->appendData(inCgSlnData, inCgSlnData.size(), 1);
          //std::cout << "Size after = " << icnt->getNData() << std::endl;
        }
      }
    }
  }
}

bool cgnsAnalyzer::isCgRindNode(int cgNodeId)
{
  // for structured zones we should avoid writing rind cells
  if (!cgRindNodeIds.empty())
  {
      auto it = std::find(cgRindNodeIds.begin(), cgRindNodeIds.end(), cgNodeId);
      return ( it != cgRindNodeIds.end() );
  } 
  else
      return false;
}

bool cgnsAnalyzer::isCgRindCell(int cgCellId)
{
  // for structured zones we should avoid writing rind cells
  if (!cgRindCellIds.empty())
  {
      auto it = std::find(cgRindCellIds.begin(), cgRindCellIds.end(), cgCellId);
      return ( it != cgRindCellIds.end() );
  } 
  else
      return false;
}

void cgnsAnalyzer::cleanRind()
{
   // sanity check
   // only supports structured meshes for now
   if (!isStructured())
      return;
   if (_rindOff)
      return;
   std::cout << "Cleaning up rind data from the grid.\n";
   // create map btw real and rind node ids
   std::map<int, int> old2NewNdeIds;
   int nNewNde = 1;
   int nRindNde = 0;
   for (int iNde=1; iNde<=nVertex; iNde++)
     if (isCgRindNode(iNde))
     {
       old2NewNdeIds[iNde] = -1;
       nRindNde++;
     }
     else 
       old2NewNdeIds[iNde] = nNewNde++;
   // sanity check
   if (nRindNde != cgRindNodeIds.size())
   {
       std::cerr << "Problem occured during rind node removal.\n";
       throw;
   }
   // remove rind node coords
   nRindNde = 0;
   std::vector<double> nxCrd, nyCrd, nzCrd;
   for (int iVrt=0; iVrt<nVertex; iVrt++)
   {
     if (old2NewNdeIds[iVrt+1] < 0)
     {
       nRindNde++;
       continue;
     }
     else
     {
       nxCrd.push_back(xCrd[iVrt]);
       nyCrd.push_back(yCrd[iVrt]);
       nzCrd.push_back(zCrd[iVrt]);
     }
   }
   xCrd = nxCrd;
   yCrd = nyCrd;
   zCrd = nzCrd;
   // sanity check
   if (nRindNde != cgRindNodeIds.size())
   {
     std::cerr << "Problem occured during rind coord removal.\n";
     throw;
   }
   // remove rind elements
   std::vector<int> newElemConn;
   for (int iElm=0; iElm <nElem; iElm++)
   {
     if (isCgRindCell(iElm+1))
         continue;
     for (int id=0; id<8; id++)
       newElemConn.push_back(old2NewNdeIds[ elemConn[iElm*8+id] ] ); 
   }
   elemConn = newElemConn;
   // update connectivity
   for (auto i : elemConn)
       i = old2NewNdeIds[i];
   // remove rind nodal data
   // remove rind elemental data
   if (!slnDataCont.empty())
   {
      // internal data index start from zero
      std::vector<int> intRindCellId, intRindNodeId;
      for (auto i : cgRindNodeIds)
          intRindNodeId.push_back(i-1);
      for (auto i : cgRindCellIds)
          intRindCellId.push_back(i-1);
      // removing rind data
      std::cout << "Cleaning up rind solution data ";
      int cntr = 29;
      for (auto sd : slnDataCont)
      {
        std::cout << "..";
        cntr+=2;
        if (cntr > 70)
        {
            cntr = 0;
            std::cout << "\n";
        }
        //std::cerr << sd->getDataName() << std::endl;
        if (sd->getDataType() == NODAL)
          sd->rmvDataIdx(intRindNodeId);
        else if (sd->getDataType() == ELEMENTAL)
          sd->rmvDataIdx(intRindCellId);
      }
      std::cout << "\n";
   }
   // fix number of nodes
   nVertex -= cgRindNodeIds.size();
   // fix number of elements
   nElem -= cgRindCellIds.size();
   // setting flag
   _rindOff = true;
}
