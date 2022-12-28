#ifndef _ODEFUNCS_T_
#define _ODEFUNCS_T_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
using namespace std;
class ODEFuncs
{
    private:
      double func1(double t);
    public:
      double getFunc(int funcType, double t, double y);
};

double ODEFuncs::func1(double t){
  return t;
}

double ODEFuncs::getFunc(int funcType, double t, double y){
  switch (funcType)
  {
  case 1:
    return ODEFuncs::func1(t);
    break;
  
  default:
    return ODEFuncs::func1(t);
    break;
  }
}

#endif