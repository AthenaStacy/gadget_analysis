//This program makes a file containing the redshifts at which each of the snapshots were taken in a simulation run

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double load_snapshot(char *fname, int files);

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


int     NumPart, Ngas;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, nh, Ne, Density;
  float  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink;
} *P;

int *Id;

double  Time, zred;

main()
{
  char path[200], input_fname[200], basename[200];
  int  n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2, j;
  float x,y,z,x1,y1, z1, delr;
  FILE *outfile;

  sprintf(path, "/work/utexas/ao/minerva");
  sprintf(basename, "tgrid");

  outfile=fopen("tgrid_z","a"); 

  for(j=101;j<=140;j++){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  zred=load_snapshot(input_fname, files);
 
  fprintf(outfile," %6.2f \n",zred);

}

fclose(outfile);
 return 0;
}

double load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.4e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  return(zred);
}


