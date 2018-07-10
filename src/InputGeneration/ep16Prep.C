#ifdef HAVE_EPIC

#include "AuxiliaryFunctions.H"
#include "inputGen.H"
#include "ep16Prep.H"
#include "exoMesh.H"

ep16Prep::~ep16Prep()
{
    if (_mdb)
        delete _mdb;
}

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
        "prep.arraysize",
        "mesh.nodeset.projectile",
        "mesh.nodeset.target",
        "mesh.elementset.projectile",
        "mesh.elementset.target",
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

    // reading mesh databse
    std::string exo_fname = _jstrm["Mesh"]["File"].as<std::string>();
    _mdb = new EXOMesh::exoMesh(exo_fname);
    _mdb->read();

    // material information
    json _matj = _jstrm["Material"];
    for (const auto& km: _matj.object_range())
    {
        std::cout << toLower(km.key()) << std::endl;
        _mat[toLower(km.key())] = km.value().as<int>();
    }

    // begin processing based on order specified
    std::vector<std::string> ord = getOrder();
    for (auto oi=ord.begin(); oi!=ord.end(); oi++)
    {
        std::string tsk = findToStr(*oi, ".");
        std::string _tsk = findFromStr(*oi, ".");
        std::string __tsk = findFromStr(_tsk, ".");
        _tsk = findToStr(_tsk,".");
        std::cout << tsk << " -- " << _tsk << " -- " << __tsk <<std::endl;
        if (!tsk.compare("prep"))
            wrtPre(_tsk,__tsk);
        else if (!tsk.compare("mesh"))
            wrtMsh(_tsk,__tsk);
        else if (!tsk.compare("misc"))
            wrtMisc(_tsk,__tsk);

    }

}


void ep16Prep::wrtCmnt(std::string cmnt)
{
    std::stringstream ss;
    ss << getCmntStr() << cmnt;
    _write(ss);
}


void ep16Prep::wrtPre(std::string _tsk, std::string __tsk)
{
    std::stringstream _tcmnt;
    std::stringstream _tstr;
    if (!_tsk.compare("case"))
    {
        _tstr.clear();
        _tstr.str(std::string());
        _tcmnt.clear();
        _tcmnt.str(std::string());    
        int _ctype = _jstrm["Case"]["Type"].as<int>();
        int _cid =  _jstrm["Case"]["Id"].as<int>();
        std::string _des = " " + _jstrm["Case"]["Description"].as<std::string>();
        _tcmnt << "CASE DESCRIPTION";
        _tstr << std::setw(5) << _ctype 
              << std::setw(5) << _cid
              << std::setw(70) << std::left << _des;
    }
    else if (!_tsk.compare("run"))
    {
        _tstr.clear();
        _tstr.str(std::string());
        _tcmnt.clear();
        _tcmnt.str(std::string());        
        int _rmde, _unt, _pat, _dplt, _ssld;
        double _tmx, _cpmx;
        _pat = 3;
        _dplt = 3;
        _cpmx = 0.0;
        _rmde = _jstrm["Run"]["Mode"].as<int>();
        _unt = _jstrm["Run"]["Unit"].as<int>();
        _tmx = _jstrm["Run"]["Tmax"].as<double>();
        _tcmnt << "Run Card\n";
        _tcmnt << getCmntStr() << " RUN UNIT  PAT DPLT      TMAX     CPMAX";
        _tstr.precision(4);
        _tstr << std::setw(5) << _rmde
              << std::setw(5) << _unt
              << std::setw(5) << _pat
              << std::setw(5) << _dplt
              << std::setw(10) << std::scientific << _tmx
              << std::setw(10) << _cpmx;
    }
    else if (!_tsk.compare("exodus"))
    {
        // data gathering
        int _npns,_ntns,_npes,_ntes,_ctact,_conv,_sym,_tet;
        double _gap;
        std::vector<std::string> ns_names,eb_names;
        ns_names = _mdb->getNdeSetNames();
        eb_names = _mdb->getElmBlkNames();
        _npns = 0;
        _ntns = 0;
        _npes = 0;
        _ntes = 0;
        for (auto ni=ns_names.begin(); ni!=ns_names.end(); ni++)
        {
            std::string pre = toLower( findToStr(*ni, "_") );
            if (!pre.compare("prj"))
                _npns++;
            else if (!pre.compare("trg"))
                _ntns++;
        }
        for (auto ei=eb_names.begin(); ei!=eb_names.end(); ei++)
        {
            std::string pre = toLower( findToStr(*ei, "_") );
            if (!pre.compare("prj"))
                _npes++;
            else if (!pre.compare("trg"))
                _ntes++;
        }
        _ctact = 1;
        _conv = 2;
        _sym = 0;
        _tet = 1;
        _gap = 0.1;

        _tstr.clear();
        _tstr.str(std::string());
        _tcmnt.clear();
        _tcmnt.str(std::string());        
        _tcmnt << "Exodus II Description Card\n";
        _tcmnt << getCmntStr() << "NPNS NTNS NPES NTES CTCT CONV  SYM  TET       GAP";
        _tstr.precision(4);
        _tstr<< std::setw(5) << _npns
             << std::setw(5) << _ntns
             << std::setw(5) << _npes
             << std::setw(5) << _ntes
             << std::setw(5) << _ctact
             << std::setw(5) << _conv
             << std::setw(5) << _sym
             << std::setw(5) << _tet
             << std::setw(10) << std::scientific << _gap;
    }
    else if (!_tsk.compare("arraysize"))
    {
        // data gathering
        int _mxn,_mxl,_mxmn,_mxsn;
        _mxn = _mdb->getNumberOfNode();
        _mxl = _mdb->getNumberOfElement();
        _mxmn = _mxn;
        _mxsn = _mxn;

        _tstr.clear();
        _tstr.str(std::string());
        _tcmnt.clear();
        _tcmnt.str(std::string());        
        _tcmnt << "Array Size/Dimension Card\n";
        _tcmnt << getCmntStr() << "      MXN       MXL      MXMN      MXSN";
        _tstr << std::setw(10) << _mxn
              << std::setw(10) << _mxl
              << std::setw(10) << _mxmn
              << std::setw(10) << _mxsn;
    }
    // write to stream
    if (!_tcmnt.str().empty())
        wrtCmnt(_tcmnt.str());
    if (!_tstr.str().empty())
        _write(_tstr);
}

void ep16Prep::wrtMsh(std::string _tsk, std::string __tsk)
{
    std::stringstream _tcmnt;
    std::stringstream _tstr;
    if (!_tsk.compare("nodeset"))
    {
        _tcmnt.clear();
        _tcmnt.str(std::string());        
        _tstr.clear();
        _tstr.str(std::string());
        std::string preStr;
        if (!__tsk.compare("projectile"))
        {
            _tcmnt << "Projectile Node Set Cards\n";
            preStr = "prj";
        }
        else if (!__tsk.compare("target"))
        {
            _tcmnt << "Target Node Set Cards\n";
            preStr = "trg";
        }
        else
        {
            std::cerr << "Error: Unknown node set type " << __tsk << "\n";
            throw;
        }
        _tcmnt << getCmntStr() << "NSET     XYZ";
        wrtCmnt(_tcmnt.str());
        for (int ns=0; ns<_mdb->getNumberOfNodeSet(); ns++)
        {
            std::string pre = findToStr(toLower(_mdb->getNdeSetName(ns)),"_");
            if (!pre.compare(preStr))
            {
                int _nset = _mdb->getNdeSetId(ns);
                _tstr << std::setw(5) << _nset << "     000";
                _write(_tstr);    
            }
        }
    }
    else if (!_tsk.compare("elementset"))
    {
        _tcmnt.clear();
        _tcmnt.str(std::string());
        _tstr.clear();
        _tstr.str(std::string());
        std::string preStr;
        if (!__tsk.compare("projectile"))
        {
            _tcmnt << "Projectile Element/Particle Set Cards\n";
            preStr = "prj";
        }
        else if (!__tsk.compare("target"))
        {
            _tcmnt << "Target Element/Particle Set Cards\n";
            preStr = "trg";
        }
        else
        {
            std::cerr << "Error: Unknown node set type " << __tsk << "\n";
            throw;
        }
        _tcmnt << getCmntStr() << "LSET MATL TYPE SSET                                                  THICK/AREA";
        wrtCmnt(_tcmnt.str());
        int _lset,_matl,_type,_sset;
        double _ta;
        std::string _blnk;
        _blnk.insert(_blnk.end(),50,' ');
        for (int es=0; es<(_mdb->getNumberOfElementBlock()); es++)
        {
            std::string es_name = toLower(_mdb->getBlockName(es));
            std::string pre = findToStr(es_name,"_");
            if (!pre.compare(preStr))
            {
                _lset = _mdb->getElmBlkId(es);
                _matl = _mat[es_name];
                if (_matl == 0)
                {
                    std::cerr << "Error: Keyword " << es_name << " is not defined in input file.\n";
                    throw;
                }
                _type = 0;
                _sset = 0;
                _ta = 0.; 
                _tstr << std::setw(5) << _lset
                      << std::setw(5) << _matl
                      << std::setw(5) << _type
                      << std::setw(5) << _sset
                      << _blnk
                      << std::setw(10) << std::scientific << _ta;
                _write(_tstr); 
            }
        }
    }
}


void ep16Prep::wrtMisc(std::string _tsk, std::string __tsk)
{
    std::stringstream _tcmnt;
    std::stringstream _tstr;
    if (!_tsk.compare("velocity"))
    {
        _tcmnt.clear();
        _tcmnt.str(std::string());
        _tcmnt << "Projectile and Target Velocities Card\n";
        _tcmnt << getCmntStr() << "    PXVEL     PYVEL     PZVEL     TXVEL     TYVEL     TZVEL";
        wrtCmnt(_tcmnt.str());

        // velocity information
        std::vector<double> _prjv,_trgv;
        json _vj = _jstrm["Velocity"];
        for (int i=0; i<3; i++)
        {
            _prjv.push_back( _vj["Projectile"][i].as<double>() );
            _trgv.push_back( _vj["Target"][i].as<double>() );
        }

        _tstr.clear();
        _tstr.str(std::string());
        _tstr.precision(8);
        _tstr << std::setw(10) << _prjv[0]
              << std::setw(10) << _prjv[1]
              << std::setw(10) << _prjv[2]
              << std::setw(10) << _trgv[0]
              << std::setw(10) << _trgv[1]
              << std::setw(10) << _trgv[2];
        _write(_tstr);
    }
    else if (!_tsk.compare("detonation"))
    {
        _tcmnt.clear();
        _tcmnt.str(std::string());
        _tcmnt << "Detonation Point Card\n";
        _tcmnt << getCmntStr() << "     XDET      YDET      ZDET";
        wrtCmnt(_tcmnt.str());

        // detonation information
        std::vector<double> _det;
        json _vj = _jstrm["Detonation"];
        for (int i=0; i<3; i++)
        {
            _det.push_back( _vj["Coordinate"][i].as<double>() );
        }

        _tstr.clear();
        _tstr.str(std::string());
        _tstr.precision(8);
        _tstr << std::setw(10) << _det[0]
              << std::setw(10) << _det[1]
              << std::setw(10) << _det[2];
        _write(_tstr);
    }

}
#endif
