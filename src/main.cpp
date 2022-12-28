#include "rkSolvers.h"
#include "odeFuncs.h"

// https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
// y[n+1] = y[n] + h*sum(b[i]*k[i]){i=1 to s}
// k[i] = f(t[n]+c[i]*h,y[n]+ h*sum(a[i][j]*k[j]){j=1 to i-1})
// Butcher tableau
// c1 | a11 a12 ... a1s
// c2 | a21 a22 ... a2s
// .
// .
// cs | as1 as2 ... ass
// ---------------------
//      b1  b2  ... bs

#include <iostream>
#include <vector>
using namespace std;
template<typename T>
vector<double> linspace(T start_in, T end_in, int num_in)
{

  vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);
  
  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
      // cout << start + delta*i << endl;
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  // cout << "size: " << linspaced.size() << endl;
  return linspaced;
}

void print_vector(vector<double> &vec)
{
  cout << "size: " << vec.size() << endl;
  for (double d : vec)
    cout << d << " ";
  cout << endl;
}

int main()
{
  vector<double> tArr = linspace(0.0,1.0,21);
  ExplicitMethods solver = ExplicitMethods("FE",0.0,tArr[1]-tArr[0],1);
  vector<double> yArr1 = solver.ForwardEuler(tArr);
  vector<double> yArr2 = solver.Midpoint(tArr);
  cout << "Forward Euler = "<< yArr1[yArr1.size()-1] << " Explicit Midpoint = " << yArr2[yArr2.size()-1] << endl;
  return 0;
}