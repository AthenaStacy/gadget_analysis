#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

main()
{

double ti_base, ti_step, ti_begin, ti_end;
double z1, z2, t1, t2;

ti_base = 4.88263e-15;
ti_step = exp(1.e9*ti_base);
//ti_step = 0.0001473;
printf("ti_step = %lg\n", ti_step);

z1=0.046;
z2=0.046+3.6e-15;

t1=5.4e8*pow(10.*z1,1.5);
t2=5.4e8*pow(10.*z2,1.5);

printf(" t1 = %20.15g\n t2 = %20.15g\n", t1, t2);

}
