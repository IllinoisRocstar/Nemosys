// -*- C++ -*-

#ifndef MAD_H_MADSINGLETON
#define MAD_H_MADSINGLETON

template<typename T> class MAdSingleton
{
 public:

  static T& instance()
    {
      static T theSingleInstance; // suppose T has a default constructor
      return theSingleInstance;
    }

 private:
  MAdSingleton();
};

#endif
