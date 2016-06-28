#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

double powerlaw(double xin, int nin);

int main ()
{
  string mystr;
  int n;
  double x, xton;

  cout << "Enter real number x for x^n: ";
  getline (cin,mystr);
  stringstream(mystr) >> x;
  cout << "The value you entered is " << x << endl;

  cout << "Enter integer n for x^n: ";
  getline (cin,mystr);
  stringstream(mystr) >> n;
  cout << "The value you entered is " << n << endl;

  xton = powerlaw(x, n);

  cout << "x^n for x = " << x << " and " << " n = " << n << " is " << xton << endl;

  return(0);
}

double powerlaw(double xin, int nin)
{

int j;
double xout;

xout = 1;

if(nin > 0) 
  for(j=0; j<nin; j++)
    xout = xout * xin;
else
  for(j=0; j>nin; j--)
    xout = xout / xin;

if(xin == 0) xout = 0;

return(xout);

}

