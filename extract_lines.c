#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv)
{
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  i, j, n, type, snapshot_number, files, random, ncount, ncounthalo1, NumCurrentTimestep, ID1;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, ztime;
  double sinkposx, sinkposy, sinkposz, sinkmass, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  FILE *infile, *outfile;


  sprintf(path, "/global/scratch/minerva");

  sprintf(basename, "sinkmasses_hiacc");
  sprintf(basenameout, "sinkmasses_cut");


  sprintf(input_fname, "%s/%s", path, basename);
  sprintf(output_fname, "%s/%s", path, basenameout);

  infile=fopen(input_fname,"r");
  outfile=fopen(output_fname,"a");

  for(i=0; i<=55971; i++)
    {
    fscanf(infile, "%lg %d %lg %d %lg %lg %lg", &ztime, &ID1, &sinkmass, &NumCurrentTimestep, &sinkposx, &sinkposy, &sinkposz);

    fprintf(outfile,"%17.13g %8d %15.6g %15.6d %15.11g %15.11g %15.11g\n", ztime, ID1, sinkmass, NumCurrentTimestep, sinkposx, sinkposy, sinkposz);

    if(i%100 == 0)
      printf("%17.13g %8d %15.6g %15.6d %15.11g %15.11g %15.11g\n", ztime, ID1, sinkmass, NumCurrentTimestep, sinkposx, sinkposy, sinkposz);
    }

  fclose(infile);
  fclose(outfile);

}






