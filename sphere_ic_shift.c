
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int write_snapshot(char *fname, int files, char *outname, float delx, float dely, float delz, float boxsize);
int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(void);

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
  char path[200], input_fname[200], basename[200], output_fname[200];
  int  n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2;
  float x, y, z, x1, y1, z1, delr;
  float delx, dely, delz, boxsize; 



  sprintf(path, "/work/utexas/ao/minerva");
  sprintf(basename, "sphere");
  snapshot_number= 40;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */ 

  boxsize=100.0;
  delx=81.5814 - (boxsize/2.e0);
  dely=91.4107 - (boxsize/2.e0);
  delz=43.2321 - (boxsize/2.e0);

  //delx= 6.05978 - (boxsize/2.e0);
  //dely= 10.4559 - (boxsize/2.e0);
  //delz=  39.694 - (boxsize/2.e0);

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(output_fname, "/work/utexas/ao/minerva/ICs/sphere_shift_040");
  Ngas = write_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  unit_conversion();  

  ncount = 0;
  ncount2 = 0;
  for(n=1;n<=NumPart;n++) { 
    x=P[n].Pos[0];
    y=P[n].Pos[1];
    z=P[n].Pos[2];
  }


  printf("ncount = %d.\n", ncount);
  printf("ncount2 = %d.\n", ncount2);

  n=0;

}




/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(void)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;

  /* physical constants in cgs units */
  GRAVITY   = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
  UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s= 1.0e5;

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
  UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

  G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);


  Xh= 0.76e0;  /* mass fraction of hydrogen */
  HubbleParam= 0.7e0;


  for(i=1; i<=NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
	  MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	}
    }
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int write_snapshot(char *fname, int files, char *outname, float delx, float dely, float delz, float boxsize)
{
  FILE *fd;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  float x, y, z, *pos, mass;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define WSKIP fwrite(&dummy, sizeof(dummy),1, outfile);

  outfile=fopen(outname,"w");
  pos= (float*)malloc(sizeof(float) * 3);

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
            header1.mass[k]=header1.mass[k]/pow(8.0,3.0);
	}

      fwrite(&dummy, sizeof(dummy), 1, outfile);
      fwrite(&header1, sizeof(header1), 1, outfile);
      fwrite(&dummy, sizeof(dummy), 1, outfile);


      if(i==0)
	allocate_memory();

      SKIP;
      WSKIP;

      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {

  	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
              
              x=P[pc_new].Pos[0];
              y=P[pc_new].Pos[1];
              z=P[pc_new].Pos[2];

              x=x/8.0;
              y=y/8.0;
              z=z/8.0;

              pos[0]=x;
              pos[1]=y;
              pos[2]=z;

              fwrite(pos, sizeof(float), 3, outfile);
	      pc_new++;
	    }
	}
      SKIP;
      WSKIP;
      printf("Ngas= %6d \n",Ngas); 
/*      printf("N_DM= %6d \n",NumPart-Ngas); */

      SKIP;
      WSKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      fwrite(&P[pc_new].Vel[0], sizeof(float), 3, outfile);
              pc_new++;
	    }
	}
      SKIP;
      WSKIP;

      SKIP;
      WSKIP;
      for(k=0,pc_new=pc;k<6;k++)
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

      if(ntot_withmasses>0)
	SKIP;
        WSKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      
	      if(header1.mass[k]==0)
		{
	      	fread(&P[pc_new].Mass, sizeof(float), 1, fd);
                //mass=&P[pc_new].Mass/pow(8.0,3.0);
                //fwrite(mass, sizeof(float), 1, outfile);
		}
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;
        WSKIP;
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      

      if(header1.npart[0]>0)
	{
	  SKIP;
          WSKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      //fread(&P[pc_sph].U, sizeof(float), 1, fd);
              P[pc_sph].U=0;
              fwrite(&P[pc_sph].U, sizeof(float), 1, outfile);
	      pc_sph++;
	    }
	  SKIP;
          WSKIP;
	}

      fclose(fd);
      fclose(outfile);
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.4e \n",Time);
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
