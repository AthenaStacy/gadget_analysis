#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 6

int load_snapshot(char *fname, int files);
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



int     NumPart, Ngas;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, nh, Density, hsm;
  float H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  float CII;
  float SiII;
  float SiIII;
  float OII;
  float Metallicity;
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
  int  i, j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax;
  double nh, nhmax, mass, metal_max, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double CII, OII, SiII, Metallicity; 
  FILE *outfile;


  sprintf(path, "/nobackupp1/astacy");

  sprintf(basename, "sng1");
  sprintf(basenameout, "snapsng1");

  for(j=526;j<=526;j=j+5){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%03d", path, basenameout, snapshot_number);
  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */

  unit_conversion();  

  outfile=fopen(output_fname, "w");

  ncount = 0;
  ncount2 = 0;

  nhmax = metal_max = 0;
  mmax=mh_mass=0;
  xmax=0;
  ymax=0;
  zmax=0;
  sl = 0;
  tmax=0;
  h2max=0;
  gammin=0;
  masstot=Metallicity=0;
  typemax = 5;

  for(n=1;n<=NumPart;n++) { 

          nh = P[n].nh;
          mass=P[n].Mass;
          temp = P[n].Temp;
          h2 = P[n].H2I;
          gam = P[n].gam;
          Metallicity = Metallicity + (P[n].Metallicity*P[n].Mass/1.e-10/.7);

          //if(nh > nhmax && Id[n] > 0)
          if(P[n].Metallicity > metal_max && 
            P[n].Pos[0] > 348.0 && P[n].Pos[0] < 352.0 && P[n].Pos[1] > 348.0 && P[n].Pos[1] < 352.0 && P[n].Pos[2] > 348.0 && P[n].Pos[2] < 352.0)
          //if(mass > mmax)
          //if(Id[n] > idmax)
          //if(Id[n] == 16035519)
          //if(Id[n] == 18338822)
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
              idmax = Id[n];
              SiII = P[n].SiII;
              CII = P[n].SiII;
              OII = P[n].OII;
              metal_max = P[n].Metallicity;
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

		for(i = 1; i < NumPart; i++)
		     {
		     dis = pow(((P[i].Pos[0]-sinkposx)*(P[i].Pos[0]-sinkposx) + (P[i].Pos[1]-sinkposy)*(P[i].Pos[1]-sinkposy) + (P[i].Pos[2]-sinkposz)*(P[i].Pos[2]-sinkposz)), 0.5);
		     dis=dis*1.e3*Time/(0.7);
		     disAU=dis*206264.806;

		     if(disAU < 10.0)
		       {
		       //vx = vx + P[i].Vel[0];
		       //vy = vy + P[i].Vel[1];
		       //vz = vz + P[i].Vel[2];
		       mh_mass = mh_mass + P[i].Mass*1.e10/0.7;
		       }
		     }
		 printf("mh_mass = %lg\n", mh_mass);
		 //vx = vx/(mh_mass)*(7.23601e-13*1.e10/0.7);
		 //vy = vy/(mh_mass)*(7.23601e-13*1.e10/0.7);
		 //vz = vz/(mh_mass)*(7.23601e-13*1.e10/0.7); 
		 printf("vx = %lg, vy = %lg, vz = %lg\n", vx, vy, vz);
	 
		 //for(n = 1; n <= NumPart; n++)
		 for(n = 1; n <= Ngas; n++)
		     {

		     x=P[n].Pos[0];
		     y=P[n].Pos[1];
		     z=P[n].Pos[2];

		     P[n].Vel[0] = P[n].Vel[0] - vx;
		     P[n].Vel[1] = P[n].Vel[1] - vy;
		     P[n].Vel[2] = P[n].Vel[2] - vz;

		     dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
		     //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
		     dis=dis*1.e3*Time/(0.7);
		     disAU=dis*206264.806;

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

		     //if(P[n].nh > 1.e5)
                     //if(P[n].nh > 1.e2 && P[n].Temp > 1.e4)
                     //if(P[n].nh > 1.e0 && P[n].Temp > 1.e4) 
                     //if(P[n].Metallicity > 1.e-4)
		       //printf("nh = %lg, gam = %lg temp = %lg, dis = %lg elec %lg, H2 = %lg HeII = %lg HDI = %lg CII = %lg SiII = %lg, Metals = %lg\n", P[n].nh, P[n].gam, P[n].Temp, dis, P[n].HII, P[n].H2I, P[n].HeII, P[n].HDI, P[n].CII, P[n].SiII, P[n].Metallicity);


		     if(P[n].nh > 1.e5 && P[n].sink > -4)
		     //if(disAU <= 1.e1 /*&& P[n].sink < -4*/)
		       masstot = masstot + P[n].Mass/1.e-10/.7;

		     if(P[n].sink > 0)
		       {
		       printf("Mass = %g sink = %g ID = %d nh = %lg Dens = %lg Temp = %lg dis = %lg, vel = %lg, vrad = %lg, vrot = %lg\n", P[n].Mass, P[n].sink, Id[n], P[n].nh, P[n].Rho, P[n].Temp, disAU, vel, vrad, vrot);
			//printf("Mass = %g sink = %g ID = %d Temp = %lg dis = %lg, vel = %lg, vrad = %lg, vrot = %lg\n", P[n].Mass, P[n].sink, Id[n], P[n].Temp, disAU, P[n].Vel[0], P[n].Vel[1], P[n].Vel[2]);
			}



		       //if(P[n].Mass > 1.5e-10 && dis < 1.e1)
		//	 printf("Problem! Mass = %g sink = %g ID = %d nh = %lg Temp = %lg dis = %lg, vel = %lg, vrad = %lg, vrot = %lg\n", P[n].Mass, P[n].sink, Id[n], P[n].nh, P[n].Temp, disAU, vel, vrad, vrot);


	    random = rand();

	    //if((random*(1.e0/RAND_MAX)< 0.05e0 && P[n].nh > 3.e-2) || P[n].nh > 2.e1 /*&& P[n].sink < 0*/ /* ||  P[n].sink > 0.5*/ ) {
	    if(P[n].nh > 1.e5 || P[n].sink > 0.5 || (random*(1.e0/RAND_MAX))< 0.02e0) 
	    //if(P[n].nh > 1.e6 || P[n].sink > 0.5 /*|| (random*(1.e0/RAND_MAX))< 0.02e0*/)
	      {
	      fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n", P[n].sink, Id[n],x,y,z,P[n].Temp,P[n].nh,P[n].HII,P[n].H2I,P[n].HDI, P[n].gam, disAU, vrad, vrot,P[n].Vel[0], P[n].Mass);
	      //fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %lg %lg %15.6g %15.6g\n", P[n].sink, Id[n],x,y,z,P[n].Temp,P[n].nh,P[n].HII,P[n].H2I,P[n].HDI, P[n].gam, disAU, P[n].Vel[0]*pow(Time, 0.5), P[n].Vel[2]*pow(Time, 0.5), P[n].hsm, P[n].Mass);
	      ncount = ncount + 1;
	      }
	  }
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
   printf("SiII = %g\n", SiII);
   printf("OII = %g\n", OII);
   printf("CII = %g\n", CII);
   printf("metal_max = %g\n", metal_max);
   printf("Metallicity = %g\n", Metallicity);
   printf("gammin = %g\n", gammin);
   printf("typemax = %lg\n", typemax);
   printf("idmax = %d\n", idmax);
   printf("sink = %g\n", P[n-1].sink);
   printf("masstot = %g\n", masstot);

  fclose(outfile);
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


  for(i=1; i<=NumPart; i++)
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






/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;
  int metal_flag = 1, element_flag=1, mass_flag=0;

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

      for(k=0, NumPart=0, ntot_withmasses=0; k<MAXREF; k++)
        printf("header1.npartTotal[%d] = %d\n", k, header1.npartTotal[k]);
      for(k=0, NumPart=0, ntot_withmasses=0; k<MAXREF; k++)
        printf("header1.npart[%d] = %d\n", k, header1.npart[k]);

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
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
      printf("Ngas= %6d \n",Ngas); 
/*      printf("N_DM= %6d \n",NumPart-Ngas); */

     for(k=0;k<6;k++)
       printf("npartTotal %6d\n",header1.npartTotal[k]);
     for(k=0;k<6;k++)
       printf("npart %6d\n",header1.npart[k]);


      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
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

      printf("ntot_withmasses = %d\n", ntot_withmasses);

   for(k = 0; k < 6; k++)
      if(header1.mass[k] == 0 && header1.npart[k] > 0)
        mass_flag = 1;

      if(mass_flag==1)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      
	      if(header1.mass[k]==0)
	      	fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(mass_flag==1)
	SKIP;

      printf("Mass = %lg\n", P[1000].Mass);
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      

      if(header1.npart[0]>1)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;
printf("hi, line 499 \n");

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

         SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].hsm, sizeof(float), 1, fd);
              pc_sph++;
            }
          SKIP;
printf("hi, line 516 \n");

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].H2I, sizeof(float), 1, fd);
              fread(&P[pc_sph].HII, sizeof(float), 1, fd);

             if(element_flag == 1)
                 {
                   fread(&P[pc_sph].CII, sizeof(float), 1, fd);
                   fread(&P[pc_sph].SiII, sizeof(float), 1, fd);
                   fread(&P[pc_sph].SiIII, sizeof(float), 1, fd);
                   fread(&P[pc_sph].OII, sizeof(float), 1, fd);
                 }

              fread(&P[pc_sph].DII, sizeof(float), 1, fd);
              fread(&P[pc_sph].HDI, sizeof(float), 1, fd);
              fread(&P[pc_sph].HeII, sizeof(float), 1, fd);
              fread(&P[pc_sph].HeIII, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);              
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              pc_sph++;
            }
          SKIP;
printf("hi, line 534 \n");

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].gam, sizeof(float), 1, fd);
              pc_sph++;
            }
          SKIP;

printf("hi, line 556 \n");
    if(metal_flag == 1)
      {
        SKIP;
        for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
          {
            fread(&P[pc_sph].Metallicity, sizeof(float), 1, fd);
            pc_sph++;
          }
        SKIP;
      }
printf("hi, line 563 \n");
/*
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].sink, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;
*/

	}

      fclose(fd);
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);
  printf("L= %15.11g \n",header1.BoxSize);
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






  











