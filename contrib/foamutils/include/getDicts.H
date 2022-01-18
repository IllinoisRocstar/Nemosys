#ifndef GETDICTS_H
#define GETDICTS_H

#include "IOdictionary.H"
#include "OFstream.H"
#include "Time.H"
#include "argList.H"

class getDicts {
 public:
  /** @brief Constructor
   **/
  getDicts(){};

  /** @brief Destructor
   **/
  ~getDicts(){};

  /** @brief Create control dict
      @param write Write the dictionary if needed
      @param pth Path where it should be written
      @param startT Start time
      @param dt Time step
      @param endT End time
      @return Dictionary pointer for in memory operations
  **/
  std::unique_ptr<Foam::dictionary> createControlDict(
      const bool write, const std::string &pth = "system/",
      const double &startT = 0.0, const double &dt = 1.0,
      const double &endT = 1.0);

  /** @brief Create fvSchemes
      @param write write Write the dictionary if needed
      @param pth Path where it should be written
      @return Dictionary pointer for in memory operations
  **/
  std::unique_ptr<Foam::dictionary> createFvSchemes(
      const bool write, const std::string &pth = "system/");

  /** @brief Create fvSolution
      @param write write Write the dictionary if needed
      @param pth Path where it should be written
      @return Dictionary pointer for in memory operations
  **/
  std::unique_ptr<Foam::dictionary> createFvSolution(
      const bool write, const std::string &pth = "system/");

  /** @brief Create dynamicMeshDict for AMR only
      @param write write Write the dictionary if needed
      @param pth Path where it should be written
      @return Dictionary pointer for in memory operations
  **/
  std::unique_ptr<Foam::dictionary> createDynamicMeshDict(
      const bool write, const std::string &pth = "constant/");

  /** @brief Write all basic dictionaries
      @param pth Path where it should be written
   **/
  void writeBasicDicts(const std::string &pth = "system/",
                       const double &startT = 0.0, const double &dt = 1.0,
                       const double &endT = 1.0);

 private:
  /** @brief Writes the foam dictionary
      @param pth Path where it should be written
      @param dictToWrite Pointer for dictionary to write
   **/
  void writeDict(const std::string &path,
                 std::unique_ptr<Foam::dictionary> &dictToWrite);
};

#endif