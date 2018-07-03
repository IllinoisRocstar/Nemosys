// NEMoSys headers
#include "inputGen.H"

// other headers
#include <iostream>
#include <algorithm>

void inputGen::setNameType(std::string fname, inpFileType ftyp, std::string key) 
{ 
    if (!key.empty())
        _key = key;
    _fn[_key].push_back(fname); 
    _tpe[_key].push_back(ftyp); 
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

