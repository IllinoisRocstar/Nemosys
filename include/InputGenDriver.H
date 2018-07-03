#ifndef INPUTGENDRIVER_H
#define INPUTGENDRIVER_H

#include <NemDriver.H>

class InputGenDriver : public NemDriver
{
  public:
    
    InputGenDriver(std::string srvName, json inputjson);
    ~InputGenDriver();

    static InputGenDriver* readJSON(json inputjson);
    static InputGenDriver* readJSON(std::string ifname);

};


#endif
