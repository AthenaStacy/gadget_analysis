#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv)
{
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  i, j, k, n, imax, jmax, kmax, n_length, random, ncount, ncounthalo1, ncount2, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax;
  double NHmax, NCmax, NOmax, NSimax, NFemax, dis, theta, phi, nh, temp, tmax, h2, h2max, gam, gammin;
  double NH, NC, NO, NSi, NFe, vrotx, vroty, vrotz, vrot, mh_mass;
  double CII, OII, SiII, Metallicity, Cmax; 
  FILE *outfile;

  sprintf(output_fname, "/nobackupp1/astacy/grb/NZB_pisn.dat");

  outfile=fopen(output_fname, "r");

  ncount = 0;
  n_length = 6502800;

  NHmax = NCmax = NOmax = NSimax = NFemax = 0.0;
  tmax=0;
  h2max=0;
  gammin=0;
  Metallicity=0;
  typemax = 5;
  Cmax = 1.75e14;

  for(n=1; n<=n_length; n++) { 

    fscanf(outfile, "%d %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg \n",
                        &i, &j, &k, &theta, &phi, &dis, &nh, &NH, &NO, &NC, &NSi, &NFe);  

    if(NH > NHmax)
     {
     NHmax = NH;
    // imax = i;
    // jmax = j;
   //  kmax = k;
     }
   
    if(NO > NOmax && NO < Cmax && k > 190)
      {
      printf("NO = %lg yes! \n", NO);
      NOmax = NO;
      imax = i;
      jmax = j;
      kmax = k;
      }

 
    if(NC > NCmax && NO < Cmax && k > 190)
      {
      NCmax = NC;
      //imax = i;
      //jmax = j;
      //kmax = k;
      }

    if(NSi > NSimax && NO < Cmax && k > 190)
      NSimax = NSi;

    if(NFe > NFemax && NO < Cmax && k > 190)
      {
      NFemax = NFe;
      //imax = i;
      //jmax = j;
      //kmax = k;
      }

    if(k == 190 && NO < 3e13)
      printf("NO = %lg\n", NO);

    }


   printf("NHmax= %lg\n", NHmax);
   printf("SiII = %lg\n", NSimax);
   printf("OII = %lg\n", NOmax);
   printf("CII = %lg\n", NCmax);
   printf("Metallicity = %g\n", Metallicity);
   printf("imax = %d\n", imax);
   printf("jmax = %d\n", jmax);
   printf("kmax = %d\n", kmax);

  fclose(outfile);

}









  











