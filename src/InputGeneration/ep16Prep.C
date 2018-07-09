#ifdef HAVE_EPIC

#include "AuxiliaryFunctions.H"
#include "inputGen.H"
#include "ep16Prep.H"
#include "exoMesh.H"
 

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
  ep->process();
  ep->write();
  return ep;
}


void ep16Prep::readJSON()
{

    std::cout << "Reading Epic 2016 Input Generation JSON.\n";
    std::string fname;
    std::string type = "short_form";

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

    // order definition
    _set_key(fname); 
    setNameType(fname, INPGEN_TXT);
    std::vector<std::string> order = {
        "prep.case", 
        "prep.run", 
        "prep.exodus", 
        "prep.array_size",
        "mesh.ndeset.projectile_node_set",
        "mesh.ndeset.target_node_set",
        "mesh.elmset.projectile_element_set",
        "mesh.elmset.target_element_set",
        "misc.velocity",
        "misc.detonation"  
    };
    setOrder(order);

    // other preps
    setCmntStr("$");
}


void ep16Prep::process()
{
    // Top level information
    wrtCmnt("EPIC 2016 INPUT FILE");
    wrtCmnt("Generated on " + getTimeStr() + " by NEMoSys");
    wrtCmnt("Short Form Description Card for ExodusII/CUBIT Data");

    // reading exodus databse
    std::string exo_fname = _jstrm["Mesh"]["File"].as<std::string>();
    EXOMesh::exoMesh* ex = new EXOMesh::exoMesh(exo_fname);
    ex->read();

    // begin processing based on order specified
    std::vector<std::string> ord = getOrder();
    for (auto oi=ord.begin(); oi!=ord.end(); oi++)
    {
        std::string tsk = findToStr(*oi, ".");
        std::string _tsk = findFromStr(*oi, ".");
        std::cout << tsk << " -- " << _tsk << std::endl;
        if (!tsk.compare("prep"))
            wrtPre(_tsk);

    }
}


void ep16Prep::wrtCmnt(std::string cmnt)
{
    std::string cmntStr = getCmntStr() + cmnt;
    _write(cmntStr);
}


void ep16Prep::wrtPre(std::string tsk)
{
    std::stringstream _tcmnt;
    std::stringstream _tstr;
    if (!tsk.compare("case"))
    {
        _tstr.clear();
        _tcmnt.clear();
        int _ctype = _jstrm["Case"]["Type"].as<int>();
        int _cid =  _jstrm["Case"]["Id"].as<int>();
        std::string _des = " " + _jstrm["Case"]["Description"].as<std::string>();
        _tcmnt << "CASE DESCRIPTION";
        _tstr << std::setw(5) << _ctype 
              << std::setw(5) << _cid
              << std::setw(70) << std::left << _des;
    }
    else if (!tsk.compare("run"))
    {
        _tstr.clear();
        _tcmnt.clear();
        int _rmde, _unt, _pat, _dplt, _ssld;
        double _tmx, _cpmx;
        _pat = 3;
        _dplt = 3;
        _cpmx = 0.0;
        _rmde = _jstrm["Run"]["Mode"].as<int>();
        _unt = _jstrm["Run"]["Unit"].as<int>();
        _tmx = _jstrm["Run"]["Tmax"].as<double>();
        _tcmnt << "Run Card";
        _tstr.precision(4);
        _tstr << std::setw(5) << _rmde
              << std::setw(5) << _unt
              << std::setw(5) << _pat
              << std::setw(5) << _dplt
              << std::setw(10) << std::scientific << _tmx
              << std::setw(10) << _cpmx
              << std::setw(5) << _ssld;
    }

    // write to stream
    if (!_tcmnt.str().empty())
        wrtCmnt(_tcmnt.str());
    if (!_tstr.str().empty())
        _write(_tstr.str());
}


#endif
