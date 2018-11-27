#include <TransferDriver.H>
#include <AuxiliaryFunctions.H>

//----------------------- Transfer Driver -----------------------------------------//
TransferDriver::TransferDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname,
                               bool checkQuality)
{
  source = meshBase::Create(srcmsh);
  target = meshBase::Create(trgmsh);
  std::cout << "TransferDriver created" << std::endl;
  Timer T;
  T.start();
  source->setCheckQuality(checkQuality);
  source->transfer(target, method);
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname); 

}


TransferDriver::TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                               std::vector<std::string> arrayNames, std::string ofname,
                               bool checkQuality)
{
  source = meshBase::Create(srcmsh);
  target = meshBase::Create(trgmsh);
  Timer T;
  T.start();
  source->setCheckQuality(checkQuality);
  source->transfer(target, method, arrayNames);
  //source->write("new.vtu");
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname); 
  std::cout << "TransferDriver created" << std::endl;
}

TransferDriver::~TransferDriver()
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
  std::cout << "TransferDriver destroyed" << std::endl;
}

TransferDriver* TransferDriver::readJSON(json inputjson)
{
  std::string srcmsh; 
  std::string trgmsh; 
  std::string outmsh; 
  std::string method;
  std::string transferAll;
  std::string checkQual;
  bool transferall = 1;
  bool checkQuality = 0;
  std::vector<std::string> arrayNames;

  srcmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Source Mesh"].as<std::string>();
  trgmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Target Mesh"].as<std::string>();
  outmsh = inputjson["Mesh File Options"]
                    ["Output Mesh File"].as<std::string>();
  method = inputjson["Transfer Options"]
                    ["Method"].as<std::string>(); 

  transferAll = inputjson["Transfer Options"]
                         ["Transfer All Arrays"].as<std::string>();
  
  if (!transferAll.compare("False") || !transferAll.compare("false"))
  {
    arrayNames = inputjson["Transfer Options"]
                          ["Array Names"].as<std::vector<std::string>>();
    transferall = 0;
  }
  
  checkQual = inputjson["Transfer Options"]
                       ["Check Transfer Quality"].as<std::string>();
  if (!checkQual.compare("True") || !checkQual.compare("true"))
  {
    checkQuality = 1;
  } 

  TransferDriver* trnsdrvobj;
  if (transferall)
  {
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method, outmsh, checkQuality);
  } 
  else
  {
    std::cout << "Transferring selected arrays:" << std::endl;
    for (int i = 0; i < arrayNames.size(); ++i)
    {
      std::cout << "\t" << arrayNames[i] << std::endl;
    }
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method, arrayNames, outmsh, checkQuality); 
  }
  
  return trnsdrvobj;

}

TransferDriver* TransferDriver::readJSON(std::string ifname)
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
    
  // checking if array
  if (inputjson.is_array())
  {
    std::cout << "Warning: Input is an array. Only first element will be processed\n";
    return TransferDriver::readJSON(inputjson[0]);
  } 
  else
  {
    return TransferDriver::readJSON(inputjson);
  }
}
