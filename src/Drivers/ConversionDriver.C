#include <ConversionDriver.H>

//----------------------- Conversion Driver -----------------------------------------//
ConversionDriver::ConversionDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname,
                               json inputjson)
{
  source = meshBase::Create(srcmsh);
  std::cout << "ConversionDriver created" << std::endl;
  Timer T;
  T.start();
  // convert to pnt mesh
  if (!method.compare("VTK->PNT"))
  {
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
  if (target)
  {
    delete target;
    target = 0;
  }
  std::cout << "ConversionDriver destroyed" << std::endl;
}

ConversionDriver* ConversionDriver::readJSON(json inputjson)
{
  std::string srcmsh; 
  std::string trgmsh; 
  std::string outmsh; 
  std::string method;
  
  srcmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Source Mesh"].as<std::string>();
  trgmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Target Mesh"].as<std::string>();
  outmsh = inputjson["Mesh File Options"]
                    ["Output Mesh File"].as<std::string>();
  method = inputjson["Conversion Options"]
                    ["Method"].as<std::string>(); 

  ConversionDriver* convdrvobj;
  convdrvobj = new ConversionDriver(srcmsh, trgmsh, method, outmsh, inputjson);
  return convdrvobj;

}

ConversionDriver* ConversionDriver::readJSON(std::string ifname)
{
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
  return ConversionDriver::readJSON(inputjson);
    
}
