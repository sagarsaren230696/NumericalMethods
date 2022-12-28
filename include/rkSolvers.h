// y[n+1] = y[n] + h*sum(b[i]*k[i]){i=1 to s}
// k[i] = f(t[n]+c[i]*h,y[n]+ h*sum(a[i][j]*k[j]){j=1 to i-1})
// Butcher tableau
// c1 | a11 a12 ... a1s
// c2 | a21 a22 ... a2s
// .
// .
// cs | as1 as2 ... ass
// ---------------------
//    | b1  b2  ... bs
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <functional>

#include "odeFuncs.h"

using namespace std;

class ExplicitMethods
{
  string typeOfSolver;
  double initVal;
  int funcType;
  double h;

  public:
    ExplicitMethods(string,double,double,int);
    vector<double> ForwardEuler(vector<double> &tArr); 
    vector<double> Midpoint(vector<double> &tArr);
};

ExplicitMethods::ExplicitMethods(string typeOfMethod,double initVal,double h,int funcType){
  this->typeOfSolver = typeOfMethod;
  this->initVal = initVal;
  this->funcType = funcType;
  this->h = h;
}

vector<double> ExplicitMethods::ForwardEuler(vector<double> &tArr){
  // Butcher tableau
  // 0 | 0
  // ------
  //   | 1
  // y[n+1] = y[n] + h*k[0]
  // k[0] = f(t[n],y[n])
  vector<double> y;
  ODEFuncs odes;
  y.push_back(this->initVal);
  cout << y[0] << endl;
  for (int i = 1; i<tArr.size(); i++){
    y.push_back(y[i-1]+this->h*odes.getFunc(funcType,tArr[i-1],y[i-1]));
  }
  return y;
}

vector<double> ExplicitMethods::Midpoint(vector<double> &tArr){
  // Butcher tableau
  // 0   | 0    0
  // 0.5 | 0.5  0
  // --------------
  //     | 0    1
  // y[n+1] = y[n] + h*sum(b[i]*k[i]){i=1 to s} sum2
  // k[i] = f(t[n]+c[i]*h,y[n]+ h*sum(a[i][j]*k[j]){j=1 to i-1}) sum1
  double c[2] = {0.0,0.5};
  double b[2] = {0.0,1.0};
  double a[2][2] = {{0,0},{0.5,0}};
  vector<double> y;
  ODEFuncs odes;
  y.push_back(this->initVal);
  cout << y[0] << endl;
  
  for (int n = 0; n<tArr.size()-1; n++){
    vector<double> k;
    k.push_back(odes.getFunc(funcType,tArr[n],y[n]));
    double sum2 = 0;
    for (int i = 1; i < 2; i++){
      double sum1 = 0;
      for (int j = 0; j < i; j++){
        sum1 = sum1 + a[i][j]*k[j];  
      }
      k.push_back(odes.getFunc(1,tArr[n]+c[i]*this->h,y[n]+this->h*sum1));
      sum2 = sum2 + b[i]*k[i];
    }
    y.push_back(y[n]+this->h*sum2);
  }
  return y;
}