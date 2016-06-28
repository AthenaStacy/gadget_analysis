#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d));
 
double probfunc(double r, double y, double rho, double b, double c, double d);



double probfunc(double epsilon, double y, double sci, double fdis, double nu, double d)
{

  double dens, g, result;
   
  // if(epsilon > sci)
  //   printf("oops, sci = %g, epsilon = %g\n", sci, epsilon);

  //result=4.*pi*fdis/nu*pow(2*(sci - epsilon), 0.5);
  result = fdis*pow((sci - epsilon), 0.5); 

  if(result != result)
    result = 0;

  return(result);
}

double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d))
{
double k1, k2, k3, k4;
k1 = step*f(x,y, a, b, c, d);
k2 = step*f(x+ step/2.0, y+k1/2.0, a, b, c, d);
k3 = step*f(x + step/2.0, y+k2/2.0, a, b, c, d);
k4 = step*f(x + step, y + k3, a, b, c, d);
return(y+(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
}
