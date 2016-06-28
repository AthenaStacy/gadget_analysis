#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20

int load_snapshot(char *fname, int files);
//int reordering(void);
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



int     NumPart, Ngas;

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  int    Id;

  double  U, Temp, nh, Rho, hsm;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
   double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy;
} *P;


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
  int  i, j, n, type, snapshot_number, files, random, ncount, ncounthalo1, ncount2, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double if_mass, if_rad, if_rmax, hii_frac, if_temp, if_dens;
  FILE *outfile;

  sprintf(path, "/nobackupp7/astacy");

  sprintf(basename, "ifront_hiacc");
  sprintf(basenameout, "snapifront_hiacc");

  //for(j=3699;j<=3700;j=j+1){
  for(j=3770;j<=3770;j=j+10){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/if_evol_hiacc", path);

  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */

  printf("hi 1, Ngas = %d, NumPart = %d\n", Ngas, NumPart);

  unit_conversion();  

  ncount = 0;
  ncount2 = 0;

  nhmax = if_mass = if_rad = if_rmax = if_temp = if_dens = hii_frac = 0;
  mmax = mh_mass = 0;
  xmax = ymax = zmax = 0;
  sl = tmax = h2max = gammin=0;
  masstot = dismax = 0;
  typemax = 5;

  printf("hi 4\n");

  for(n=0;n<Ngas;n++) { 

          nh = P[n].nh;
          mass=P[n].Mass;
          temp = P[n].Temp;
          h2 = P[n].H2I;
          gam=P[n].gam;

          //if(nh > nhmax)
          if(mass > mmax)
            {
              //printf("Found the sink!\n");
              nhmax = nh;
              mmax=mass;
              xmax=P[n].Pos[0];
              ymax=P[n].Pos[1];
              zmax=P[n].Pos[2];
              sl = P[n].hsm;
              tmax=P[n].Temp;
              h2max=P[n].H2I;
              gammin=P[n].gam;
              typemax = P[n].Type;
              idmax = P[n].Id;
              vx = P[n].Vel[0];
              vy = P[n].Vel[1];
              vz = P[n].Vel[2];
            }
    }

         sinkposx=xmax;
         sinkposy=ymax;
         sinkposz=zmax;
         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g, nhmax = %15.11g, mmax = %15.11g\n", sinkposx, sinkposy, sinkposz, nhmax, mmax);
 

         if(nhmax < 1.e1)
           {
           sinkposx=header1.BoxSize/2.0;
           sinkposy=header1.BoxSize/2.0;
           sinkposz=header1.BoxSize/2.0; 
           }
         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g\n", sinkposx, sinkposy, sinkposz);
         printf("vx = %lg, vy = %lg, vz = %lg\n", vx, vy, vz);


        for(i = 0; i < NumPart; i++)
             {
             dis = pow(((P[i].Pos[0]-sinkposx)*(P[i].Pos[0]-sinkposx) + (P[i].Pos[1]-sinkposy)*(P[i].Pos[1]-sinkposy) + (P[i].Pos[2]-sinkposz)*(P[i].Pos[2]-sinkposz)), 0.5);
             dis=dis*1.e3*Time/(0.7);
             disAU=dis*206264.806;

             if(disAU < 10.0)
               {
               mh_mass = mh_mass + P[i].Mass*1.e10/0.7;
               }
             }
         printf("mh_mass = %lg\n", mh_mass);
         printf("vx = %lg, vy = %lg, vz = %lg\n", vx, vy, vz);
 

         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);

             if(dis > dismax)  dismax = dis;

             dis=dis*1.e3*Time/(0.7);
             disAU=dis*206264.806;

             if(P[n].HII > 0.9 && P[n].sink < 0.5)
              {
              if_mass = if_mass + P[n].Mass * 1.e10 / 0.7;
              if_rad  = if_rad  + dis * (P[n].Mass * 1.e10 / 0.7);
              if_temp = if_temp + P[n].Temp * (P[n].Mass * 1.e10 / 0.7);
              if_dens = if_dens + P[n].nh * (P[n].Mass * 1.e10 / 0.7);
              hii_frac = hii_frac + P[n].HII * (P[n].Mass * 1.e10 / 0.7);
              if(dis > if_rmax) if_rmax = dis;              
              }


             vel = P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1] + P[n].Vel[2]*P[n].Vel[2];
             vel = pow(vel,0.5)*pow(Time, 0.5); 

             vrad =  (P[n].Vel[0]*(P[n].Pos[0]-sinkposx) + P[n].Vel[1]*(P[n].Pos[1]-sinkposy) + P[n].Vel[2]*(P[n].Pos[2]-sinkposz))/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             vrad = vrad*pow(Time, 0.5); 

             vrotx = (P[n].Pos[1]-sinkposy)*P[n].Vel[2]  - (P[n].Pos[2]-sinkposz)*P[n].Vel[1];
             vroty = (P[n].Pos[2]-sinkposz)*P[n].Vel[0] - (P[n].Pos[0]-sinkposx)*P[n].Vel[2];
             vrotz = (P[n].Pos[0]-sinkposx)*P[n].Vel[1]  -  (P[n].Pos[1]-sinkposy)*P[n].Vel[0];


             vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5); 
             vrot = vrot*pow(Time, 0.5)/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
;  //convert to km/s

             if(P[n].nh > 1.e0 && P[n].Temp > 2.e4 && P[n].sink == 0)
               printf("ID = %d, nh = %lg, temp = %lg, elec %lg, H2 = %lg HeII = %lg mass = %lg\n", P[n].Id, P[n].nh, P[n].Temp, P[n].HII, P[n].H2I, P[n].HeII, P[n].Mass);

             if(P[n].nh > 1.e10 && P[n].sink > -4)
               masstot = masstot + P[n].Mass/1.e-10/.7;

             if(P[n].nh > 1.0e3 && P[n].Temp < 0.1*2.73/Time)
               {
               printf("Mass = %g sink = %g ID = %d nh = %lg Dens = %lg Temp = %lg dis = %lg, vel = %lg, vrad = %lg, vrot = %lg\n", P[n].Mass, P[n].sink, P[n].Id, P[n].nh, P[n].Rho, P[n].Temp, disAU, vel, vrad, vrot);
                }

  }

  
   if(if_mass > 0)
     {
     if_rad = if_rad/if_mass;  
     if_temp = if_temp / if_mass;
     if_dens = if_dens / if_mass;
     hii_frac = hii_frac / if_mass;
     }
   else
    {
    if_mass = 1.e-30;
    if_rad = 1.e-30;
    if_temp = 1.e-30;
    if_dens = 1.e-30;
    hii_frac = 1.e-30;
    if_rmax = 1.e-30;
    }

   outfile=fopen(output_fname, "a");
   fprintf(outfile, "%d %15.10g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n", snapshot_number, Time, if_mass, if_rad, if_rmax, hii_frac, if_temp, if_dens); 
   fclose(outfile);

   printf("ncount = %d.\n", ncount);
   printf("ncount2 = %d.\n", ncount2);

   printf("nhmax= %g\n", nhmax);
   printf("mmax= %g\n", mmax);
   printf("xmax= %15.11g\n", xmax);
   printf("ymax= %15.11g\n", ymax);
   printf("zmax= %15.11g\n", zmax);
   printf("sl = %g\n", sl);
   printf("tmax = %g\n", tmax);
   printf("h2max = %g\n", h2max);
   printf("gammin = %g\n", gammin);
   printf("typemax = %lg\n", typemax);
   printf("idmax = %d\n", idmax);
   printf("sink = %g\n", P[n-1].sink);
   printf("masstot = %g\n", masstot);
   printf("dismax = %g\n", dismax);

   free(P);
}

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
int unit_conversion(void)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;
  double h2frac, muh2, muh2in;
 
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
 
       MeanWeight=1.2195;
       h2frac=2.0*P[i].H2I;

       muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
       muh2=pow(muh2in, -1.0);

       if(muh2 >= 1.22)
         {
          MeanWeight=muh2;
         }

          MeanWeight=MeanWeight*PROTONMASS;

	  //MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  //gamma= 5.0/3.0;
          gamma=P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs;
	  P[i].nh= P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
	  /*  printf("zred = %g", zred);*/
	}
    }
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

//  P--;   /* start with offset 1 */

 /* 
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  */
//  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc=0,pc_new=0,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

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
      printf("Ngas= %6d \n",Ngas); 
      printf("NumPart= %6d \n",NumPart); 

     for(k=0;k<6;k++)
       printf("npartTotal %6d\n",header1.npartTotal[k]);
     for(k=0;k<6;k++)
       printf("npart %6d\n",header1.npart[k]);


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
	      fread(&P[pc_new].Id, sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;



      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      
	      if(header1.mass[k]==0)
	      	fread(&P[pc_new].Mass, sizeof(double), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;

      printf("Mass = %lg\n", P[1000].Mass);
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      

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
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);              
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
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

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].sink, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;


	}

      fclose(fd);
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %17.13e \n",Time);
  printf("L= %15.11g \n",header1.BoxSize);
  return(Ngas);
}




/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
/*
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
*/





  











