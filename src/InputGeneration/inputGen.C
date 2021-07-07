// Nemosys headers
#include "InputGeneration/inputGen.H"

// other headers
#include <iostream>
#include <fstream>

#include "AuxiliaryFunctions.H"


void inputGen::setNameType(const std::string &fname,
                           inpFileType ftyp,
                           const std::string &key)
{
  if (!key.empty())
    _key = key;
  _fn[_key] = fname;
  _tpe[_key] = ftyp;
}


void inputGen::setOrder(const std::vector<std::string> &__ord,
                        const std::string &key)
{
  if (!key.empty())
    _key = key;
  std::vector<std::string> ord = __ord;
  for (auto &&it : ord)
    nemAux::toLower(it);
  _ord[_key].insert(_ord[_key].end(), ord.begin(), ord.end());
}


std::vector<std::string> inputGen::getOrder(const std::string &key)
{
  if (!key.empty())
    _key = key;
  return _ord[_key];
}


void inputGen::pushOrder(const std::string &__ord, const std::string &key)
{
    if (!key.empty())
        _key = key;
    std::string ord = __ord;
    nemAux::toLower(ord);
    _ord[_key].push_back(ord);
}


void inputGen::setMsh(meshBase *mb, const std::string &key)
{
  if (!key.empty())
    _key = key;
  _mb[_key].push_back(mb);
}


void inputGen::setCmntStr(const std::string &cmstr, const std::string &key)
{
  if (!key.empty())
    _key = key;
  _cmnt[_key] = cmstr;
}


std::string inputGen::getCmntStr(const std::string &key)
{
  if (!key.empty())
    _key = key;
  return _cmnt[_key];
}

void inputGen::write(const std::string &key) const
{
  bool onlyKey = false;
  if (!key.empty())
    onlyKey = true;
  for (const auto &it : _inp)
  {
    if (onlyKey && it.first != key)
      continue;

    std::string fname = _fn.at(it.first);
    std::ofstream ofile;
    ofile.open(fname);
    if (!ofile.good())
    {
      std::cerr << "Error opening file " << fname << std::endl;
      throw;
    }
    ofile << (it.second)->str();
    ofile.close();
  }
}

