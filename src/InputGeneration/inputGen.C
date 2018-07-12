// NEMoSys headers
#include "inputGen.H"

// other headers
#include <iostream>
#include <fstream>
#include <algorithm>

void inputGen::setNameType(std::string fname, inpFileType ftyp, std::string key) 
{ 
    if (!key.empty())
        _key = key;
    _fn[_key] = fname; 
    _tpe[_key] = ftyp; 
}

void inputGen::setOrder(std::vector<std::string>& ord, std::string key)
{
    if (!key.empty())
        _key = key;
    for (auto it=ord.begin(); it!=ord.end(); it++)
        std::transform(it->begin(), it->end(), it->begin(), ::tolower);
    _ord[_key].insert(_ord[_key].end(), ord.begin(), ord.end());
}

std::vector<std::string> inputGen::getOrder(std::string key)
{
    if (!key.empty())
        _key = key;
    return(_ord[_key]);
}


void inputGen::pushOrder(std::string ord, std::string key)
{
    if (!key.empty())
        _key = key;
    std::transform(ord.begin(), ord.end(), ord.begin(), ::tolower);
    _ord[_key].push_back(ord);
}

void inputGen::setMsh(meshBase* mb, std::string key)
{
    if (!key.empty())
        _key = key;
    _mb[_key].push_back(mb);
}


void inputGen::setCmntStr(std::string cmstr, std::string key)
{
    if (!key.empty())
        _key = key;
    _cmnt[_key] = cmstr;
}

std::string inputGen::getCmntStr(std::string key)
{
    if (!key.empty())
        _key = key;
    return(_cmnt[_key]);
}

void inputGen::write(std::string key)
{
    bool onlyKey = false;
    if (!key.empty())
        onlyKey = true;
    for (auto it=_inp.begin(); it!=_inp.end(); it++)
    {
        if (onlyKey && (it->first).compare(key))
            continue;

        std::string fname = _fn[(it->first)]; 
        std::ofstream ofile;
        ofile.open (fname);
        if (!ofile.good())
        {
            std::cerr << "Error opening file "<< fname << std::endl;
            throw;
        }
        ofile << (it->second).str();
        ofile.close();
    }
}

