//Adds gamma values of the particles to a new snapshot file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define PI 3.14159265359
#define BFF 1

int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(void);

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


struct io_header_old
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
} header_old;


int     NumPart, Ngas;

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  double dis;
  double disx, disy, disz;
  double  Rho, U, Temp, nh, Density, hsm;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
   double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy;
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
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  j, n, type, snapshot_number, output_number,files, Ngas, random, ncount, ncounthalo1, ncount2;
  float x,y,z,x1,y1, z1, delr;
  double delx, dely, delz, boxsize;
  FILE *outfile;

  sprintf(path, "/work/00863/minerva/");
  sprintf(basename, "bin_HR10_map");

output_number = 0;
for(j=0; j<=0; j=j+1)
  {

  snapshot_number=j;
  
  files=1;     

  boxsize = 140.0;

  delx = 512.65800353;  //del's used in making bin_zoom10_cut_??? (bin_zoom10_???)
  dely = 505.21836888;
  delz = 497.39349349;

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  if(snapshot_number > 999)
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  if(snapshot_number > 9999)
    sprintf(input_fname, "%s/%s_%05d", path, basename, snapshot_number);
  
  sprintf(basenameout, "ot");
  sprintf(output_fname, "%s/%s_%04d", path, basenameout, output_number);
  if(snapshot_number > 999)
    sprintf(output_fname, "%s/%s_%04d", path, basenameout, output_number);
  if(snapshot_number > 9999)
    sprintf(output_fname, "%s/%s_%05d", path, basenameout, output_number);
  output_number++;

  Ngas = write_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  unit_conversion();  

  do_what_you_want();
  free(P);
}
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


  for(i=0; i<NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
	  MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  //gamma= 5.0/3.0;
          gamma = P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs;
	  P[i].nh= P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
	  /*  printf("zred = %g", zred);*/
	}
    }
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize)
{
  FILE *fd;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,l,dummy,ntot_withmasses; 
  int    t,n,off,pc=0,pc_new,pc_sph;
  int NumPart_new = 0, Ngas_new = 0;
  int Idnew, nnew=1;
  double *pos, *binit, massnew, hsmnew, Unew, x, y, z, dis, nnew_doub=1.0;
  double randomx, randomy, randomz, nh, nhmax, rho, rhomax; 
  double ngrid_xleft=1, ngrid_yleft=1, ngrid_zleft=1, ngrid_x=1600, ngrid_y=1600, ngrid_z=1;
  double ngrid_tot = 0;
  double frac_x=0, frac_y=0, frac_z=0;  
  double num_left=0, num=0, num_tot;
  double dmax = 1.e10;  //size of 'cut-out' square in pc
  double mass_max = 1.e12;

  pos= (double*)malloc(sizeof(double) * 3);

  binit= (double*)malloc(sizeof(double) * 3);
  binit[0] = 1.e-20; binit[1] = 0, binit[2] = 0;
	
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define WSKIP fwrite(&dummy, sizeof(dummy),1, outfile);

  outfile=fopen(outname,"w");

  for(i=0, pc=0; i<files; i++, pc=pc_new)
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
            {
	    NumPart+= header1.npart[k];
            printf("NumPart[%d] = %d\n", k, header1.npart[k]);
            }
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

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
              pc_new++;
            }
        }
      SKIP;


      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&Id[pc_new], sizeof(int), 1, fd);
              pc_new++;
            }
        }
      SKIP;

      if(ntot_withmasses>0)
      {
        SKIP;
     }

      for(k=0, pc_new=pc; k<6; k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              P[pc_new].Type=k;
              if(header1.mass[k]==0)
                {
                fread(&P[pc_new].Mass, sizeof(double), 1, fd);
                }
              else
                P[pc_new].Mass= header1.mass[k];
              pc_new++;
            }
        }
      if(ntot_withmasses>0)
      {
        SKIP;
      }


		if(header1.npart[0]>0)
		{
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].U, sizeof(double), 1, fd);
			    pc_sph++;
			}
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].Rho, sizeof(double), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].hsm, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
				fread(&P[pc_sph].HII, sizeof(double), 1, fd);
				fread(&P[pc_sph].DII, sizeof(double), 1, fd);
				fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
				fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
				fread(&P[pc_sph].HeIII, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].gam, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;
		
            printf("gam = %lg\n",P[100].gam);
	
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].sink, sizeof(double), 1, fd); 
				pc_sph++;
            }
			SKIP;
			
        }		

 nhmax = 0;
 for(n=0;n<Ngas;n++) 
          {
          nh = P[n].Rho;
          if(nh > nhmax)
            {
            nhmax = nh;
            delx = P[n].Pos[0];
            dely = P[n].Pos[1];
            delz = P[n].Pos[2];
            }
          }

      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
              {
              P[pc_new].dis = 1.e3*header1.time/(0.7)*pow(((P[pc_new].Pos[0]-delx)*(P[pc_new].Pos[0]-delx) + (P[pc_new].Pos[1]-dely)*(P[pc_new].Pos[1]-dely) + (P[pc_new].Pos[2]-delz)*(P[pc_new].Pos[2]-delz)), 0.5);
              P[pc_new].disx =fabs(1.e3*header1.time/(0.7)*(P[pc_new].Pos[0] - delx));
              P[pc_new].disy =fabs(1.e3*header1.time/(0.7)*(P[pc_new].Pos[1] - dely));
              P[pc_new].disz =fabs(1.e3*header1.time/(0.7)*(P[pc_new].Pos[2] - delz));
              pc_new++;
              }
          }

      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              //if(P[pc_new].dis < 10.0)
              if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax)  
                NumPart_new++;
              //if(P[pc_new].dis < 10.0 && k ==0)
              if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
                Ngas_new++;
              pc_new++;
            }
        }

      for(k=0;k<5;k++)
        header_old.npart[k] = header1.npart[k];

//////////////////////////////////////////////////////////////////////////////////////////
////Reset header values!
/////////////////////////////////////////////////////////////////////////////////////////

      NumPart_new = int(ngrid_x) * int(ngrid_y) * int(ngrid_z);
      Ngas_new = int(ngrid_x) * int(ngrid_y) * int(ngrid_z);
      header1.npartTotal[0] = Ngas_new;
      header1.npartTotal[1] = 0;
      header1.npart[0] = Ngas_new;
      header1.npart[1] = 0;

      header1.time = 0.01;
      header1.redshift = 0;
      header1.BoxSize = 1.0;
      header1.HubbleParam = 1.0;

      printf("Ngas= %6d \n",Ngas); 
      for(k=0;k<6;k++)
        printf("Ngas_new %6d\n",header1.npart[k]);
     for(k=0;k<6;k++)
       printf("Ngas_new_tot %6d\n",header1.npartTotal[k]);
     for(k=0;k<6;k++)
       printf("Ngas_old_tot %6d\n",header1.npartTotal[k]);
      printf("NumPart_new %6d\n",NumPart_new);
      //printf("N_DM= %6d \n",NumPart-Ngas); 


//////////write new file!!!!!!!!!!!
      sprintf(buf,"%s",outname);
      printf("writing `%s' ...\n",buf);
			
      fwrite(&dummy, sizeof(dummy), 1, outfile);
      fwrite(&header1, sizeof(header1), 1, outfile);
      fwrite(&dummy, sizeof(dummy), 1, outfile);

      printf("header finished\n");

      num_tot = double(NumPart_new);
      ngrid_tot = pow(num_tot, 1./2.);

      double frac_x, frac_y, frac_z;
      int i, j, ngrid_x_int, ngrid_y_int, ngrid_z_int;
      ngrid_x_int = int(ngrid_x);
      ngrid_y_int = int(ngrid_y);
      ngrid_z_int = int(ngrid_z);

      WSKIP;
      pc_new = 0;
      for(i=0; i<ngrid_x_int; i++)
         for(j=0; j<ngrid_y_int; j++)
             for(k=0; k<ngrid_z_int; k++)
                 {
                 x = double(i+1) / ngrid_x / header1.BoxSize /*- header1.BoxSize/2.*/;
                 y = double(j+1) / ngrid_y / header1.BoxSize /*- header1.BoxSize/2.*/; 
                 z = (double(k + 0.5) - ngrid_z/2.) / ngrid_x  / header1.BoxSize;  
                 //z = (double(k) - ngrid_z/2.) / ngrid_x  / header1.BoxSize;


                 P[pc_new].Pos[0]=x;
                 P[pc_new].Pos[1]=y;
                 P[pc_new].Pos[2]=z;
			
	         pos[0]=x;
		 pos[1]=y;
		 pos[2]=z;
						
		 fwrite(pos, sizeof(double), 3, outfile);	
                 pc_new++;
                 }
      WSKIP;

      printf("positions finished, num = %lg, num_left = %lg\n", num, num_left);

      WSKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {

              if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
                {
                P[pc_new].Vel[0] = -sin(2 * PI * P[pc_new].Pos[1]);
                P[pc_new].Vel[1] =  sin(2 * PI * P[pc_new].Pos[0]);
                P[pc_new].Vel[2] = 0; 
                fwrite(&P[pc_new].Vel[0], sizeof(double), 3, outfile);
                }
	      pc_new++;
	    }
	}
      WSKIP;

      printf("velocities finished\n");

      WSKIP;
      Idnew = 0;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
              if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
                {
                fwrite(&Idnew, sizeof(int), 1, outfile);
                Idnew++;
		}
	      pc_new++;
	    }
	}
      WSKIP;


      if(ntot_withmasses>0)
      {
        WSKIP;
     }

     double rho_init, mass_tot, volume;
     rho_init = 25. / (36. * PI);
     volume = header1.BoxSize * header1.BoxSize * header1.BoxSize * 0.015625;
     mass_tot = rho_init * volume;     
     massnew = mass_tot / num_tot; 

      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      if(header1.mass[k]==0)
                {
                if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
                  {  
                  fwrite(&massnew, sizeof(double), 1, outfile);
                  }
                }
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
      {
        WSKIP;
      }

      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      

   double pres_init = 5.0 / (12.*PI);
   Unew = pres_init / ((5.0/3.0) - 1.0) / rho_init;
   printf("Unew = %lg\n", Unew);

   if(header_old.npart[0]>0)
	{
          WSKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max) 
                {
                if(n%100000 == 0)
                     printf("U = %lg\n",P[pc_sph].U);
	         fwrite(&Unew, sizeof(double), 1, outfile);
                }
	      pc_sph++;
	    }
          WSKIP;

          WSKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                  fwrite(&rho_init, sizeof(double), 1, outfile);  
                }
	      pc_sph++;
	    }
          WSKIP;


//divergence between gadget1 and gadget2 begins here!

          WSKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
           {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                if(n%100000 == 0)
                     printf("hsm = %lg\n",P[pc_sph].hsm);
		   //hsmnew = P[pc_sph].hsm*pow(nnew_doub,-.3333);
		   if(P[pc_sph].Pos[0] < 0.5*header1.BoxSize)
                     hsmnew = 0.01;
                   if(P[pc_sph].Pos[0] >= 0.5*header1.BoxSize)
                     hsmnew = 0.01;
                   fwrite(&hsmnew, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;

          WSKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                fwrite(&P[pc_sph].H2I, sizeof(double), 1, outfile);
                fwrite(&P[pc_sph].HII, sizeof(double), 1, outfile);
                fwrite(&P[pc_sph].DII, sizeof(double), 1, outfile);
                fwrite(&P[pc_sph].HDI, sizeof(double), 1, outfile);
                fwrite(&P[pc_sph].HeII, sizeof(double), 1, outfile);
                fwrite(&P[pc_sph].HeIII, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;


          WSKIP;
          double gam_init = 5./3.;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                 if(n%100000 == 0)
                   printf("gam = %lg\n",P[pc_sph].gam);
                fwrite(&gam_init, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;

          WSKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max) 
                {
                fwrite(&P[pc_sph].sink, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;

#if(BFF)
          double bwrite;
          printf("writing bfields!\n");
          WSKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                for(l=0;l<nnew;l++)
                  {
                  rho = P[pc_sph].Rho;

                  bwrite = binit[0];
                  fwrite(&bwrite, sizeof(double), 1, outfile);

                  bwrite = binit[1];
                  fwrite(&bwrite, sizeof(double), 1, outfile);

                  bwrite = binit[2];
                  fwrite(&bwrite, sizeof(double), 1, outfile);
                  }
                }
              pc_sph++;
            }
          WSKIP;
#endif

 
	}

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

  if(!(P=(struct particle_data *) malloc(3*NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
//  P--;   /* start with offset 1 */

  
  if(!(Id=(int *) malloc(3*NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
//  Id--;   /* start with offset 1 */

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






  











