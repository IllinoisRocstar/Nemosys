#ifdef HAVE_EPIC

#include "AuxiliaryFunctions.H"
#include "inputGen.H"
#include "ep16Prep.H"
 

ep16Prep* ep16Prep::readJSON(std::string ifname)
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
    inputjson = inputjson[0];
  } 
  
  ep16Prep* ep = new ep16Prep(inputjson);
  ep->readJSON();
  return ep;
}

ep16Prep* ep16Prep::readJSON(json inputjson)
{
  ep16Prep* ep = new ep16Prep(inputjson);
  ep->readJSON();
  return ep;
}


void ep16Prep::readJSON()
{

    std::cout << "Reading Epic 2016 Input Generation JSON.\n";
    std::string fname;
    std::string type;

    // reading mandatory fields
    fname = _jstrm["File Name"].as<std::string>();
    if (_jstrm.has_key("Type"))
    {
        type = _jstrm["Type"].as<std::string>();
        type = toLower(type);
    }
   
    // sanity checking
    if (type.compare("short_form"))
    {
        std::cerr << "Error: Only short form is supported.\n";
        exit(-1);
    }

    _set_key(fname); 
    setNameType(fname, INPGEN_TXT);
    std::vector<std::string> order = {
        "prep.description", 
        "prep.run", 
        "prep.exodus", 
        "prep.array_size",
        "ndeset.projectile_node_set",
        "ndeset.target_node_set",
        "elmset.projectile_element_set",
        "elmset.target_element_set",
        "misc.velocity",
        "misc.detonation"          
    };
    setOrder(order);
    std::cout << findToStr("sdfsdfsd.555555",".") << "\n";
    std::cout << findFromStr("sdfsdfsd.555555",".") << "\n";

}


#endif
