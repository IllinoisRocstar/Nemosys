#include <algorithm>
#include <map>

// Nemosys
#include <ConversionDriver.H>
#include "meshSrch.H"
#include <AuxiliaryFunctions.H>

// vtk
#include <vtkIdList.h>
#include <vtkCellData.h>

//----------------------- Conversion Driver -----------------------------------------//
ConversionDriver::ConversionDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname,
                               json inputjson)
    : source(NULL)
{
  std::cout << "ConversionDriver created" << std::endl;
  Timer T;
  T.start();
  // convert to pnt mesh
  if (!method.compare("VTK->PNT"))
  {
    source = meshBase::Create(srcmsh);
    // number of dimensions
    int dim = inputjson["Conversion Options"]
                       ["Dimension"].as<int>();
    // looping through blocks
    int nBlk = inputjson["Conversion Options"]["Block Data"].size();
    std::cout << "Number of Blocks : " << nBlk << std::endl;
    PNTMesh::BlockMap elmBlkMap;
    elmBlkMap.resize(nBlk);
    int iBlk = 0;
    for (auto it = inputjson["Conversion Options"]
                            ["Block Data"].array_range().begin();
              it!= inputjson["Conversion Options"]
                            ["Block Data"].array_range().end();
              ++it)
    {
      // element ids for block
      std::cout << (*it)["Name"].as<std::string>() << std::endl;
      std::vector<int> ids;
      if (it->has_key("Element ID Range"))
      {
        int rs,re;
        rs = (*it)["Element ID Range"][0].as<int>();
        re = (*it)["Element ID Range"][1].as<int>();
        std::cout << "Range " << rs << " to " << re << std::endl;
        for (int eid= rs;
                 eid <= re;
                 eid++)
          ids.push_back(eid);
      }
      else
      {
        for (auto it2 = (*it)["Element ID"].array_range().begin();
                  it2 != (*it)["Element ID"].array_range().end();
                  it2++)
        {
          ids.push_back( (*it2).as<int>() );
          std::cout << (*it2).as<int>() << " ";
        }
        std::cout << std::endl;
      }
      // registering block information
      elmBlkMap[iBlk].ordIntrp = (*it)["Element Order"].as<int>();
      elmBlkMap[iBlk].ordEquat = (*it)["Equation Order"].as<int>();
      elmBlkMap[iBlk].eTpe = PNTMesh::elmTypeNum( (*it)["Element Type"].as<std::string>() );
      elmBlkMap[iBlk].regionName = (*it)["Name"].as<std::string>();
      std::string bcTagName = (*it)["BC Tag"].as<std::string>();
      elmBlkMap[iBlk].srfBCTag.push_back( PNTMesh::bcTagNum( bcTagName ) );
      elmBlkMap[iBlk].elmIds = ids;
      iBlk++;
    }
    PNTMesh::pntMesh* pm = new PNTMesh::pntMesh(source, dim, nBlk, elmBlkMap);
    pm->write(trgmsh);
    delete pm;
    pm = NULL;
  }
  else if (!method.compare("EXODUSII"))
  {
#ifdef HAVE_EXODUSII // NEMOSYS is compiled with exodus
    // reading vitals
    std::cout << "Converting to EXODUSII...\n";
    json opts = inputjson["Conversion Options"];
    int nMsh = opts.get_with_default("Number of Mesh", 0);
    bool needsPP = opts.get_with_default("Post Processing", false);
    
    // sanity check
    if (nMsh == 0)
    {
        std::cerr << "Error: At least one mesh should be provided!\n";
        exit(-1);
    }

    // starting conversion operation
    EXOMesh::exoMesh* em = new EXOMesh::exoMesh(ofname);
    
    // reading meshes
    int ndeIdOffset = 0;
    int ins = 0;
    int ieb = 0;
    for (int iMsh=0; iMsh < nMsh; iMsh++)
    {
        json mshOpts = opts["Mesh Data"];
        std::string mshFName = mshOpts[iMsh].get_with_default("File", "");
        std::string mshName = mshOpts[iMsh].get_with_default("Name", "default");
        bool usePhys = mshOpts[iMsh].get_with_default("Use Physical Groups", false);

        // reading input mesh
        meshBase* mb = meshBase::Create(mshFName);
        //mb->write("exo_inp_mesh_"+std::to_string(iMsh)+".vtu");

        // adding information to exodusII object
        if (!usePhys)
        {
            // node coordinate to one nodeSet
            EXOMesh::ndeSetType ns;
            ns.id = ++ins;
            ns.name = mshName;
            ns.usrNdeIds = false; // automatically node ids
            ns.nNde = mb->getNumberOfPoints();
            for (int iNde=0; iNde<ns.nNde; iNde++)
                ns.crds.push_back(mb->getPoint(iNde));
            em->addNdeSet(ns);

            // all elements into one element block
            // assuming input mesh has only one element type
            EXOMesh::elmBlockType eb;
            eb.id = ++ieb;
            eb.name = mshName;
            eb.ndeIdOffset = ndeIdOffset;
            int iTet = 0; 
            for (int iElm=0; iElm<mb->getNumberOfCells(); iElm++)
            {
                EXOMesh::elementType eTpe = EXOMesh::v2eEMap((VTKCellType) (mb->getDataSet()->GetCellType(iElm)));
                if (eTpe != EXOMesh::elementType::TETRA)
                    continue;
                iTet++;
                vtkIdList* nids = vtkIdList::New();
                mb->getDataSet()->GetCellPoints(iElm, nids);
                for (int in=0; in< nids->GetNumberOfIds(); in++)
                    eb.conn.push_back(nids->GetId(in)+1);  // offset node ids by 1
            }
            std::cout << "Number of tetrahedral elements = " << iTet << std::endl;
            std::cout << "Min nde indx = " << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
                      << "Max nde indx = " << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
            std::cout << "Size connectivity check = " << eb.conn.size()/4 << std::endl;
            std::cout << "Starting node offset = " << ndeIdOffset << std::endl;
            eb.nElm = iTet;
            eb.eTpe = EXOMesh::elementType::TETRA;
            eb.ndePerElm = 4;
            em->addElmBlk(eb);
            // offseting starting node id for next file
            ndeIdOffset += ns.nNde;
        }
        else
        {
            // get number of physical groups
            // loop through cell data and identify physical groups
            // we also filter only TETRA cells
            int nc = mb->getNumberOfCells();
            // check physical group exist and obtain id
            vtkCellData *cd = mb->getDataSet()->GetCellData();
            int physGrpArrIdx = -1;
            if (cd)
            {
                for (int ida=0; ida < cd->GetNumberOfArrays(); ida++)
                {
                    std::cout << cd->GetArrayName(ida) << std::endl;
                    if (!strcmp( cd->GetArrayName(ida), "PhysGrpId"))
                    {
                        physGrpArrIdx = ida;
                        break;
                    }
                }
            }
            if (physGrpArrIdx < 0)
            {
                std::cerr << "Error : Input dataset does not have PhyGrpId cell data! Aborting.\n";
                exit(-1);
            }
            vtkDataArray* physGrpIds = cd->GetArray(physGrpArrIdx);
            // loop trhough elements and obtain physical group ids
            std::vector<int> elmIds;
            std::vector<int> elmPhysGrp;
            for (int ic=0; ic<nc; ic++)
            {
                EXOMesh::elementType eTpe = EXOMesh::v2eEMap((VTKCellType) (mb->getDataSet()->GetCellType(ic)));
                if (eTpe != EXOMesh::elementType::TETRA)
                    continue;
                elmIds.push_back(ic);
                double* tmp = physGrpIds->GetTuple(ic);
                elmPhysGrp.push_back(int(*tmp));
            }
            // number of unique physical groups
            std::set<int> unqPhysGrpIds;
            int nPhysGrp;
            for (auto it = elmPhysGrp.begin(); it != elmPhysGrp.end(); it++)
                unqPhysGrpIds.insert(*it);
            nPhysGrp = unqPhysGrpIds.size();
            std::cout << "Number of physical groups : " << nPhysGrp << std::endl;

            // one nodeset for all groups
            // node coordinate to nodeSet
            EXOMesh::ndeSetType ns;
            ns.id = ++ins;
            ns.name = mshName;
            ns.usrNdeIds = false; // automatically node ids
            ns.nNde = mb->getNumberOfPoints();
            for (int iNde=0; iNde<ns.nNde; iNde++)
                ns.crds.push_back(mb->getPoint(iNde));
            em->addNdeSet(ns);
            
            // for each physical group one element block
            for (auto it1=unqPhysGrpIds.begin(); it1!=unqPhysGrpIds.end(); it1++)
            {
                // element to elementBlock
                EXOMesh::elmBlockType eb;
                eb.id = ++ieb;
                eb.name = mshName+"_PhysGrp_"+std::to_string(ieb);
                eb.ndeIdOffset = ndeIdOffset;
                int iTet = 0; 
                for (auto it2=elmIds.begin(); it2!=elmIds.end(); it2++)
                {
                    double* tmp2 = physGrpIds->GetTuple(*it2);
                    if (int(*tmp2) != *it1)
                        continue;
                    iTet++;
                    vtkIdList* nids = vtkIdList::New();
                    mb->getDataSet()->GetCellPoints(*it2, nids);
                    for (int in=0; in< nids->GetNumberOfIds(); in++)
                        eb.conn.push_back(nids->GetId(in)+1);  // offset node ids by 1
                }
                std::cout << "Number of group tetrahedral elements = " << iTet << std::endl;
                std::cout << "Min nde indx = " << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
                          << "Max nde indx = " << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
                std::cout << "Size connectivity check = " << eb.conn.size()/4 << std::endl;
                std::cout << "Starting node offset = " << ndeIdOffset << std::endl;
                eb.nElm = iTet;
                eb.eTpe = EXOMesh::elementType::TETRA;
                eb.ndePerElm = 4;
                em->addElmBlk(eb);
            }
            // offset starting node id for next file
            ndeIdOffset += ns.nNde;
        }

        // clean up
        delete mb;
        mb = NULL;
    }

    // writing the file
    em->write();
    em->report();

    // performing post-processing tasks
    if (needsPP)
    {   
        int nTsk = opts.get_with_default("Number of Tasks", 0);
        json ppTsk = opts["Tasks"];
        for (int iTsk=0; iTsk<nTsk; iTsk++)
        {
            std::string ppFName= ppTsk[iTsk].get_with_default("File", "");
            std::cout << "Reading Post Processing JSON file "<< iTsk << std::endl;
            std::ifstream inputStream(ppFName);
            if (!inputStream.good() || find_ext(ppFName) != ".json")
            {
              std::cout << "Error opening file " << ppFName << std::endl;
              exit(1);
            }
            if (find_ext(ppFName) != ".json")
            {
              std::cout << "Input File must be in .json format" << std::endl;
              exit(1);
            }
            json ppJson;
            inputStream >> ppJson;
            procExo(ppJson, ofname, em);
        }

        // writing augmented exo file
        em->write();
        em->report();
    }
    
    // clean up
    delete em;
    em = NULL;
#else
    std::cerr << "Error: Compile NEMOSYS with ENABLE_EXODUS to use this option.\n";
    exit(-1);
#endif
  }
  T.stop();
}


ConversionDriver::~ConversionDriver()
{
  if (source)
  {
    delete source;
    source = 0;
  }
  std::cout << "ConversionDriver destroyed" << std::endl;
}

ConversionDriver* ConversionDriver::readJSON(json inputjson)
{

  std::cout << "Reading JSON object\n";
  std::string srcmsh; 
  std::string trgmsh; 
  std::string outmsh; 
  std::string method;
 
  srcmsh = "";
  trgmsh = "";
  outmsh = "";

  if (inputjson.has_key("Mesh File Options"))
  {
      if (inputjson["Mesh File Options"].has_key("Input Mesh Files"))
      {
          if (inputjson["Mesh File Options"]["Input Mesh Files"].has_key("Source Mesh"))
            srcmsh = inputjson["Mesh File Options"]
                                ["Input Mesh Files"]
                                ["Source Mesh"].as<std::string>();

          if (inputjson["Mesh File Options"]["Input Mesh Files"].has_key("Target Mesh"))
            trgmsh = inputjson["Mesh File Options"]
                                ["Input Mesh Files"]
                                ["Target Mesh"].as<std::string>();
          
          if (inputjson["Mesh File Options"].has_key("Output Mesh File"))
            outmsh = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();
      }
  }

  // minimal json input is "Conversion Options"->"Method"
  method = inputjson["Conversion Options"]
                    ["Method"].as<std::string>(); 

  // starting proper conversion driver, right now just one
  ConversionDriver* convdrvobj;
  convdrvobj = new ConversionDriver(srcmsh, trgmsh, method, outmsh, inputjson);
  return convdrvobj;

}

ConversionDriver* ConversionDriver::readJSON(std::string ifname)
{
  std::cout << "Reading JSON file\n";
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
  
  // checking if array
  if (inputjson.is_array())
  {
    std::cout << "Warning: Input is an array. Only first element will be processed\n";
    return ConversionDriver::readJSON(inputjson[0]);
  } 
  else
  {
    return ConversionDriver::readJSON(inputjson);
  }
    
}


#ifdef HAVE_EXODUSII
void ConversionDriver::procExo(json ppJson, std::string fname, EXOMesh::exoMesh* em)
{
  // converting to mesh base for geometric inquiry
  meshBase* mb = meshBase::Create(fname); 
  meshSrch* ms = meshSrch::Create(mb);

  // performing requested operation
  std::string opr = ppJson.get_with_default("Operation", "");
  if (!opr.compare("Material Assignment"))
  {
      // gathering information about all zones
      std::map<std::string, std::set<int> > zoneGeom;
      json zones = ppJson["Zones"];
      int nZn = zones.size();
      for (int iZn=0; iZn<nZn; iZn++)
      {
          // assuming first element is zone infomration
          // keyed by zone name that we do not care about
          // yet
          json znInfo = zones[iZn][0];
          std::string matName = znInfo.get_with_default("Material Name","N/A");
          std::string shape = znInfo.get_with_default("Shape","N/A");
          std::cout <<"Processing zone "<<iZn<<" Material "<<matName<<" Shape "<<shape<<std::endl;

          if (!shape.compare("Box"))
          {
              std::vector<double> bb;
              bb.push_back( znInfo["Params"]["Min"][0].as<double>() ); 
              bb.push_back( znInfo["Params"]["Max"][0].as<double>() ); 
              bb.push_back( znInfo["Params"]["Min"][1].as<double>() ); 
              bb.push_back( znInfo["Params"]["Max"][1].as<double>() ); 
              bb.push_back( znInfo["Params"]["Min"][2].as<double>() ); 
              bb.push_back( znInfo["Params"]["Max"][2].as<double>() ); 
              
              std::vector<int> lst;
              ms->FindCellsWithinBounds(bb, lst, true);
              zoneGeom[matName].insert(lst.begin(), lst.end());
              
          }
          else
              std::cout << "WARNNING: Skipiing unkown zone shape: " << shape << std::endl;     
      }

      // adjusting exodus database accordingly
      for (auto it1=zoneGeom.begin(); it1!=zoneGeom.end(); it1++)
      {
          std::vector<int> elmLst;
          elmLst.insert(elmLst.end(), (it1->second).begin(), (it1->second).end());
          em->addElmBlkByElmIdLst(it1->first, elmLst);
      }
  }
  else if (!opr.compare("Check Duplicate Elements"))
  {
      std::cout << "Checking for existance of duplicate elements ... ";
      bool ret = ms->chkDuplElm(); 
      if (ret)
      {
          std::cerr << " The exodus database contains duplicate elements.\n";
          exit(-1);
      } 
      else 
        std::cout << "False\n";
  }
  else if (!opr.compare("Remove Block"))
  {
      std::string blkName = ppJson.get_with_default("Block Name", "");
      std::cout << "Removing Block " << blkName << std::endl;
      em->removeElmBlkByName(blkName);
  }
  else if (!opr.compare("Snap Node Coords To Zero"))
  {
      double tol = ppJson.get_with_default("Tolerance", 0.0);
      std::cout << "Snapping nodal coordinates to zero using tolerance " << tol << std::endl;
      em->snapNdeCrdsZero(tol);
  }
  else
  {
      std::cout << "Unknown operation requested : " << opr << std::endl;
  }
}
#endif

