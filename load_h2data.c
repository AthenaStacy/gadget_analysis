#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv)
{
  char path[200], input_fname[200], output_fname1[200], output_fname2[200], output_fname3[200]; 
  char basename[200], out1[200], out2[200], out3[200];
  int  i, j, n, type, snapshot_number, files, random, ncount, ncounthalo1, NumCurrentTimestep, ID1;
  int nh2op = 101;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, ztime;
  double h2_in1, h2_in2, h2_in3;
  double h2_opac_temp[nh2op], h2_opac_column[nh2op], h2_opac[nh2op][nh2op], vrotx, vroty, vrotz, vrot, mh_mass;
  FILE *infile, *outfile1, *outfile2, *outfile3;


  sprintf(path, "/home1/00863/minerva/gadget_runs/");

  sprintf(basename, "H2-cooling-ratios.dat");
  sprintf(out1, "h2_opac_temp.dat");
  sprintf(out2, "h2_opac_column.dat");
  sprintf(out3, "h2_opac.dat");

  sprintf(input_fname, "%s/%s", path, basename);
  sprintf(output_fname1, "%s/%s", path, out1);
  sprintf(output_fname2, "%s/%s", path, out2);
  sprintf(output_fname3, "%s/%s", path, out3);

  infile=fopen(input_fname,"r");
  outfile1=fopen(output_fname1,"a");
  outfile2=fopen(output_fname2,"a");
  outfile3=fopen(output_fname3,"a");

  for(i=0; i<nh2op; i++)
    for(j=0; j<nh2op; j++)
     {
     fscanf(infile, "%lg %lg %lg", &h2_in1, &h2_in2, &h2_in3);

     h2_opac_temp[i] = h2_in1;
     h2_opac_column[j] = h2_in2;
     h2_opac[i][j] = h2_in3;
     }

  for(i=0; i<nh2op; i++)
    {
    fprintf(outfile1,"%10.5g, ", h2_opac_temp[i]);
    fprintf(outfile2,"%10.5g, ", h2_opac_column[i]);
    }

  for(i=0; i<nh2op; i++)
    for(j=0; j<nh2op; j++)
    {
    if(j==0) fprintf(outfile3,"{ ");

    if(j==nh2op-1) fprintf(outfile3,"%10.5g },\n ",  h2_opac[i][j]);
    else fprintf(outfile3,"%10.5g, ", h2_opac[i][j]);
    }


  fclose(infile);
  fclose(outfile1);
  fclose(outfile2);
  fclose(outfile3);

}






