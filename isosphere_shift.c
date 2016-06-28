//This file reads in float point values for r^-2 minihalo particles, and writes out new file with double precision values.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 6
#define MAXREF2 20

int write_snapshot(char *fname,int files, char *outname, double delx, double dely, double delz, double delx_dm, double dely_dm, double delz_dm, double  boxsize);
int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(void);
int allocate_memory_dm(void);

struct io_header_1
{
  int      npart[MAXREF];
  double   mass[MAXREF];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[MAXREF];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


struct io_header_2
{
  int      npart[MAXREF2];
  double   mass[MAXREF2];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[MAXREF2];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header2;


int     NumPart, Ngas;

struct particle_data
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  U;
} *P;

int *Id;

double  Time, zred;



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char path[200],  input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2;
  float x,y,z,x1,y1, z1, delr;
  double delx, dely, delz, delx_dm, dely_dm, delz_dm, boxsize;
  FILE *outfile;

  files = 1;
  boxsize = 100.0;
  delx = 50.339110402;  //NOTE: In this code the delx, delx_dm, etc. values are NOT RELEVANT
  dely = 50.122129376;
  delz = 49.486652655;

  delx_dm = 3.4331373996e-06;
  dely_dm = 3.4336503546e-06;
  delz_dm = 3.432893655e-06;
 

  sprintf(input_fname, "/work/00863/minerva/minihalo_128");  //"nfw" denotes an r^-1 profile!
  sprintf(output_fname, "/work/00863/minerva/iso_map_000");

  //sprintf(output_fname, "/work/00863/minerva/ds10_hr_000");

  Ngas = write_snapshot(input_fname, files, output_fname, delx, dely, delz, delx_dm, dely_dm, delz_dm, boxsize);

  //unit_conversion();  

  ncount = 0;
  ncount2 = 0;

  printf("ncount = %d.\n", ncount);
  printf("ncount2 = %d.\n", ncount2);

  do_what_you_want();
}


/* here the particle data is at your disposal 
 */
int do_what_you_want(void)
{

}





/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */




/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double delx_dm, double dely_dm, double delz_dm, double  boxsize)
{
  FILE *fd;
  FILE *fd_dm;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  int NumPart_new = 0, Ngas_new = 0;
  double dum1[3], dum2, dis_cm;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define WSKIP fwrite(&dummy, sizeof(dummy),1, outfile);

  outfile=fopen(outname,"w");

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


      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npart[k];
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      
////////////////////reset header for the new output file
      for(k=0; k<MAXREF; k++)
        {
        header2.npart[k] = header1.npart[k];
        header2.npartTotal[k] = header1.npartTotal[k];
        header2.mass[k] = header1.npartTotal[k];
        }
      for(k=MAXREF; k<MAXREF2; k++)
        {
        header2.npart[k] = 0;
        header2.npartTotal[k] = 0;
        header2.mass[k] = 0;
        }

      header2.time = header1.time;
      header2.redshift = header1.redshift;
      header2.flag_sfr = header1.flag_sfr;
      header2.flag_feedback = header1.flag_feedback;
      header2.flag_cooling = header1.flag_cooling;
      header2.num_files = header1.num_files;
      header2.BoxSize = header1.BoxSize;
      header2.Omega0 = header1.Omega0;
      header2.OmegaLambda = header1.OmegaLambda;
      header2.HubbleParam = header1.HubbleParam;
/////////////////////////////

      printf("Ngas= %6d \n",Ngas); 
      printf("N_DM %6d\n", header1.npart[1]);
      printf("Omega0 = %lg\n", header2.Omega0);
      printf("HubbleParam = %lg\n", header2.HubbleParam);

      fwrite(&dummy, sizeof(dummy), 1, outfile);
      fwrite(&header2, sizeof(header2), 1, outfile);
      fwrite(&dummy, sizeof(dummy), 1, outfile);


      if(i==0)
        {
	allocate_memory();
        }
  
      SKIP;
      WSKIP;


 printf("line 324\n");
      for(k=0,pc_new=pc;k<1;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
              dum1[0] = (double)P[pc_new].Pos[0];
              dum1[1] = (double)P[pc_new].Pos[1];
              dum1[2] = (double)P[pc_new].Pos[2];
              fwrite(&dum1, sizeof(double), 3, outfile);
              pc_new++;
            }
        }
      SKIP;
      WSKIP;

      printf("line 309\n");
      printf("x = %g, y = %g, z = %g\n", P[2].Pos[0], P[2].Pos[1], P[2].Pos[2]);
      printf("x = %lg, y = %lg, z = %lg\n", P[pc_new-1].Pos[0], P[pc_new-1].Pos[1], P[pc_new-1].Pos[2]);

      SKIP;
      WSKIP;
      for(k=0,pc_new=pc;k<1;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
            fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
            dum1[0] = (double)P[pc_new].Vel[0];
            dum1[1] = (double)P[pc_new].Vel[1];
            dum1[2] = (double)P[pc_new].Vel[2];
            fwrite(&dum1, sizeof(double), 3, outfile);
            pc_new++;
            }
        }
      SKIP;
      WSKIP;

      printf("vx = %lg, vy = %lg, vz = %lg\n", P[2].Vel[0], P[2].Vel[1], P[2].Vel[2]);
      printf("vx = %lg, vy = %lg, vz = %lg\n", P[pc_new-2].Vel[0], P[pc_new-2].Vel[1], P[pc_new-2].Vel[2]);

      SKIP;
      WSKIP;
      for(k=0,pc_new=pc;k<1;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&Id[pc_new], sizeof(int), 1, fd);
              fwrite(&Id[pc_new], sizeof(int), 1, outfile);
              pc_new++;
            }
        }
      SKIP;
      WSKIP;
     
      printf("line 360\n");
      printf("Id = %d\n", Id[0]);

      if(ntot_withmasses>0)
      {
	SKIP;
        WSKIP;
      }

      for(k=0, pc_new=pc; k<1; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      if(header1.mass[k]==0)
                {
       	        fread(&P[pc_new].Mass, sizeof(double), 1, fd);
                fwrite(&P[pc_new].Mass, sizeof(double), 1, outfile);
                }
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      for(k=0, pc_new=pc; k<1; k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              P[pc_new].Type=k;
              if(header1.mass[0]==0)
                {
                fread(&P[pc_new].Mass, sizeof(float), 1, fd);
                dum2 = (double)P[pc_new].Mass; 
                fwrite(&dum2, sizeof(double), 1, outfile);
                }
              else
                {
                P[pc_new].Mass= header1.mass[0];
                dum2 = (double)P[pc_new].Mass;
                fwrite(&dum2, sizeof(double), 1, outfile);
                }
              pc_new++;
            }
        }
      if(ntot_withmasses>0)
      {
	SKIP;
        WSKIP;
      }
     
      printf("Mass = %lg\n", P[1].Mass);
      printf("Mass = %lg\n", P[pc_new-1].Mass); 
 
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      

      fclose(fd);
      fclose(outfile);
       
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  return(Ngas);
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=(struct particle_data *) malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}



/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
}






  











