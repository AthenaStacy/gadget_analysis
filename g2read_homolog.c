#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define hubble_param 1.0
#define PI 3.14159265359

#define MAXREF 20

#define readB 0
#define vpot 0
#define radprof 0
#define BFF 1
#define def_fac 0
#define exp_fac 1.0

#define width_small 0.1

#define width_outer 0.1

#define ref_lev 32

//#define snapbegin 100
//#define snapend   1700
#define snapbegin 6
#define snapend   12

#define snapcheck  snapend
#define snapcenter snapend

#define snapnum (snapend - snapbegin + 1)

#define hfac_out 1.0
#define neighb_num 100

int load_snapshot(char *fname, int files);
//int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(void);
double calc_det(double matrix11, double matrix12, double matrix13, double matrix21, double matrix22, double matrix23, double matrix31, double matrix32, double matrix33);
double calc_kernel_spline(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double cosmo_fac);
double del_kernel(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double del_rho(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double curl(int dir, double x1, double y1, double z1, double x2, double y2, double z2);

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
//#ifndef BFF
  double  Vel[3];
//#endif
  double  Mass;
  int    Type;
  int    Id;

  //double  U; 
  double Temp, nh, Rho, hsm; 
  //double gam, sink;
  //double gam, H2I, HII, DII, HDI, HeII, HeIII;
  double dummy;
  double Bfieldx, Bfieldy, Bfieldz;
#if(readB)
  double Density, nh_test, eror;
  double Afield[3], fcor;
#endif
#if(radprof)
  double disx, disy, disz, hsm_pc;
#endif
  int to_print;
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
  char path[200], pathout[200], path2[200], basename2[200], input_fname[200], input_fname2[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  int  i, j, k, m, n, type, snapshot_number, files, random, ncount, ncounthalo1, ncount2, idmax;
  int    snaparr[snapnum];
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vxCOM=0, vyCOM=0, vzCOM=0, vel, typemax, dismax; 
  double sq_error=0, rms_error=0, tot_error=0;
  double delfac=0.9, n0, b0, a0, a0_x, a0_y, n_anal, b_anal;
  double nh, nhmax, mass, mmax, dis, xmax=0, ymax=0, zmax=0, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double disx, disy, disz, time_fac, cosmo_fac, cosmo_fac0;
  FILE *outfile, *outfile2, *infile;

  //sprintf(basename, "homolog3");
  //sprintf(basenameout, "snaphomolog3");

  //sprintf(basename, "homolog4");
  //sprintf(basenameout, "snaphomolog4");

  sprintf(basename, "homologR3");
  sprintf(basenameout, "snaphomologR3");

  //sprintf(basename, "homologR4");
  //sprintf(basenameout, "snaphomologR4");

  sprintf(path, "/scratch/00863/minerva");
  sprintf(pathout, "/work/00863/minerva");

  for(n=0;n<snapnum;n++)
    snaparr[n] = snapbegin+(1*n);

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

  int nskip = 1;
  for(j=0;j<snapnum;j=j+nskip){

  snapshot_number= snaparr[j];                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */
  double arrnum = 0;

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);

  sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  

  Ngas = load_snapshot(input_fname, files);
 
  printf("t_Hubble = %lg \n", 5.4e8/pow((1.e0+header1.redshift)/10.e0, 1.5e0));
  printf("hi 1, Ngas = %d, NumPart = %d\n", Ngas, NumPart);

  unit_conversion();  

  time_fac = Time;

  //For NON-cosmological runs
  time_fac = 1.0;
  cosmo_fac = time_fac/(hubble_param);

  printf("hi 2\n");

  outfile=fopen(output_fname, "w");

  ncount = 0;
  ncount2 = 0;

  printf("hi 3\n");

//////////////////////////////////////////////////////////////////////////////////

  nhmax = mmax = mh_mass = sl = tmax = h2max = gammin = masstot = dismax = 0;
  typemax = 5;

  printf("cosmo_fac = %lg \n", cosmo_fac);

  for(n=0;n<Ngas;n++) { 

          nh = P[n].nh;
          mass=P[n].Mass;

          P[n].Pos[0] = P[n].Pos[0]*cosmo_fac*exp_fac;
          P[n].Pos[1] = P[n].Pos[1]*cosmo_fac*exp_fac;
          P[n].Pos[2] = P[n].Pos[2]*cosmo_fac*exp_fac;
          P[n].hsm = P[n].hsm*cosmo_fac*hfac_out*exp_fac;

          P[n].to_print = 0;

          if(nh > nhmax /*&& P[n].sink > -1*/)
          //if(mass > mmax)
            {
              //printf("Found the sink!\n");
              nhmax = nh;
              mmax=mass;

              xmax=P[n].Pos[0];
              ymax=P[n].Pos[1];
              zmax=P[n].Pos[2];

              sl = P[n].hsm;
              idmax = P[n].Id;
//#ifndef BFF
              vx = P[n].Vel[0];
              vy = P[n].Vel[1];
              vz = P[n].Vel[2];
//#endif
            }
    }
  printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);
 
  double ref_lev_doub, DeltaX, hsm_pc, num[ref_lev], dis_arr[ref_lev]; 
  double Bx[ref_lev], By[ref_lev], Bz[ref_lev], Bmag[ref_lev], rho[ref_lev];
  double velx[ref_lev], vely[ref_lev], velz[ref_lev], temp[ref_lev];
  double xgridL, ygridL, zgridL, xgridR, ygridR, zgridR;
  int take_profile = 1; 
 
  xgridL = ygridL = zgridL = xgridR = ygridR = zgridR = -width_small/2.0;

  ref_lev_doub = (double) ref_lev;
  DeltaX = width_small / ref_lev_doub;
  printf("DeltaX = %lg\n", DeltaX);

  for(i=0;i<ref_lev;i++) {
  num[i] = Bx[i] = By[i] = Bz[i] = Bmag[i] = rho[i] = velx[i] = vely[i] = velz[i] = temp[i] = 0; 
  if(take_profile != 1) dis_arr[i] = xgridL + DeltaX/2;
  if(take_profile == 1) dis_arr[i] = (i+0)*width_small/2/ref_lev; 
  xgridL = xgridL + DeltaX;
  printf("dis_arr[%d] = %lg\n", i, dis_arr[i]); 
  }


  double nh_typ=0, b_typ=0, u_typ=0, num_typ=0, bfield, ufield, u_norm, bx_typ=0, by_typ=0, bz_typ=0, ax_typ=0, ay_typ=0, az_typ=0, a_anal=0, a_typ;
  double ax_diff, ay_diff, az_diff, rcheck, dAxdy, dAxdz, dAydx, dAydz, dAzdx, dAzdy;
  double nmin = 1.e-2;
  double dwdu[3][neighb_num]; 
  double pos[3][neighb_num], Afield[3][neighb_num];
  double jacob[3][3], jacob_inv[3][3], det_jacob, jacob2[3][3], jacob_fin[3][3];

  int dest, p, q;
  struct particle_data psave, psource;

  printf("arrnum = %d\n", arrnum);

         sinkposx=header1.BoxSize/2.0;
         sinkposy=header1.BoxSize/2.0;
         sinkposz=header1.BoxSize/2.0; 

         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

//#ifndef BFF
             vel = P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1] + P[n].Vel[2]*P[n].Vel[2];
             vel = pow(vel,0.5)*pow(Time, 0.5); 

             vrad =  (P[n].Vel[0]*(P[n].Pos[0]-sinkposx) + P[n].Vel[1]*(P[n].Pos[1]-sinkposy) + P[n].Vel[2]*(P[n].Pos[2]-sinkposz))/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             vrad = vrad*pow(Time, 0.5); 

             vrotx = (P[n].Pos[1]-sinkposy)*P[n].Vel[2]  - (P[n].Pos[2]-sinkposz)*P[n].Vel[1];
             vroty = (P[n].Pos[2]-sinkposz)*P[n].Vel[0] - (P[n].Pos[0]-sinkposx)*P[n].Vel[2];
             vrotz = (P[n].Pos[0]-sinkposx)*P[n].Vel[1]  -  (P[n].Pos[1]-sinkposy)*P[n].Vel[0];
//#endif

             vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5); 
             vrot = vrot*pow(Time, 0.5)/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
;  //convert to km/s

              vx = P[n].Vel[0]*pow(Time, 0.5); 
              vy = P[n].Vel[1]*pow(Time, 0.5);  
              vz = P[n].Vel[2]*pow(Time, 0.5);

    x = P[n].Pos[0];  y = P[n].Pos[1];  z = P[n].Pos[2]; 

    hsm_pc = P[n].hsm * 1.e3;

    random = rand();
    int to_print = 0; double bmin = 1.e-19;
    if(random*(1.e0/RAND_MAX)< 0.1e0) to_print++;


    if(to_print > 0)
      {
       //printf("hello line 891\n");

       fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
        P[n].dummy, P[n].Id,x,y,z,P[n].Temp,P[n].nh,vx,vy,vz, P[n].dummy, disAU, vrad, vrot, hsm_pc, P[n].Mass, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].dummy);

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
   printf("gammin = %g\n", gammin);
   printf("typemax = %lg\n", typemax);
   printf("idmax = %d\n", idmax);
   printf("masstot = %g\n", masstot);
   printf("dismax = %g\n", dismax);

   sq_error = sq_error/tot_error;
   rms_error = pow(sq_error, 0.5);

  printf("tot_error = %lg, sq_error = %lg, rms_error = %lg\n", tot_error, sq_error, rms_error);

  fclose(outfile);
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

/*
       h2frac=2.0*P[i].H2I;

       muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
       muh2=pow(muh2in, -1.0);

       if(muh2 >= 1.22)
         {
          MeanWeight=muh2;
         }

          MeanWeight=MeanWeight*PROTONMASS;
*/
	  MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].Temp * UnitEnergy_in_cgs/ UnitMass_in_g;

	  gamma= 5.0/3.0;
	  //gamma = 1.1;
          //gamma=P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs;
	  P[i].nh=  P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
          P[i].Rho= P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
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
//#ifndef BFF
	      fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
/*
#else
              fread(&P[pc_new].dummy, sizeof(double), 1, fd);
              fread(&P[pc_new].dummy, sizeof(double), 1, fd);
              fread(&P[pc_new].dummy, sizeof(double), 1, fd);
#endif
*/
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
      printf("M_B= %15.6e\n",header1.mass[0]);
      printf("M_DM= %15.6e %lg\n",header1.mass[1], P[NumPart-100].Mass);
      

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Temp, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

         printf("U = %lg\n", P[1000].Temp);

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

         printf("Density = %lg\n", P[1000].Rho);

         SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].hsm, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

          printf("hsm = %lg\n", P[1000].hsm);

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
/*
              fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
              fread(&P[pc_sph].HII, sizeof(double), 1, fd);
              fread(&P[pc_sph].DII, sizeof(double), 1, fd);
              fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeIII, sizeof(double), 1, fd);
*/
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);

              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);              
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              //fread(&P[pc_sph].dummy, sizeof(float), 1, fd);
              pc_sph++;
            }
          SKIP;

          printf("H2I = %lg\n", P[1000].dummy);

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

         printf("gam = %lg\n", P[1000].dummy);

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

         printf("sink = %lg\n", P[1000].dummy);

#if(BFF)
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].Bfieldx, sizeof(double), 1, fd);
              fread(&P[pc_sph].Bfieldy, sizeof(double), 1, fd);
              fread(&P[pc_sph].Bfieldz, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

          printf("Bfieldx = %lg\n", P[1000].Bfieldx);

#endif

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

double calc_kernel_spline(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
  double hsm_phys, fac=1.0;
  double radx, rady, radz;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  radx = fabs(x_part - x);
  rady = fabs(y_part - y);
  radz = fabs(z_part - z);

  ratio = rad/hsm;

  if(ratio >= 1. && radx < grid_size_half && rady < grid_size_half && radz < grid_size_half)
     {
     ratio = rad/(grid_size_half);
     if(hsm < grid_size_half)
       hsm = grid_size_half;
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, hsm*1.e3);
     }

  if(ratio <= 0.5)
    kernel = fac*(8./PI/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./PI/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;

  return(kernel);
}


double del_kernel(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp)
{
  double dkernel;
  double rad, del_kernel, ratio, drdu;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  ratio = rad/hsm;
  drdu = 0.5 * pow((x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z), -0.5 ) * 2. * rad_comp;

  dkernel = 0;

  if(ratio <= 0.5)
    dkernel = (8./3.14159/pow(hsm,3)) * ( -12.*pow(ratio,1)*(1./hsm)*drdu + 18.*pow(ratio,2)*(1./hsm)*drdu );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = (8./3.14159/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(-1./hsm)*drdu;
  if(ratio > 1.)
    dkernel = 0.;

  return(dkernel);
}

double del_rho(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp)
{
  double dkernel;
  double rad, del_kernel, ratio, drdu;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  ratio = rad/hsm;
  drdu = 0.5 * pow((x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z), -0.5 ) * 2. * rad_comp;

  dkernel = 0;

  if(ratio <= 0.5)
    dkernel = -3.*(8./3.14159/pow(hsm,4)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3))
             + (8./3.14159/pow(hsm,3)) * ( - 12.*pow(ratio,1)*(-rad/hsm/hsm) +  18*pow(ratio,2)*(-rad/hsm/hsm) );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = -3.*(8./3.14159/pow(hsm,4)) * 2.*pow(1. - ratio, 3) + (8./3.14159/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(rad/hsm/hsm);
  if(ratio > 1.)
    dkernel = 0.;

  return(dkernel);
}

double curl(int dir, double x1, double y1, double z1, double x2, double y2, double z2)
{

double xnew, ynew, znew;

xnew = y1*z2 - z1*y2;
ynew = z1*x2 - x1*z2;
znew = x1*y2 - y1*x2;

if(dir == 0) return(xnew);
else if(dir == 1) return(ynew);
else if(dir == 2) return(znew);
else return(0);
}


double calc_det(double matrix00, double matrix01, double matrix02, double matrix10, double matrix11, double matrix12,
                double matrix20, double matrix21, double matrix22)
{
int p, q, r;
double determinant, trA1, trA2, trA3, A2[3][3], A3[3][3], matrix[3][3];

matrix[0][0] = matrix00 + def_fac;
matrix[0][1] = matrix01;
matrix[0][2] = matrix02;
matrix[1][0] = matrix10;
matrix[1][1] = matrix11 + def_fac;
matrix[1][2] = matrix12;
matrix[2][0] = matrix20;
matrix[2][1] = matrix21;
matrix[2][2] = matrix22 + def_fac;

determinant =       matrix[0][0]*matrix[1][1]*matrix[2][2]
                  + matrix[0][1]*matrix[1][2]*matrix[2][0]
                  + matrix[0][2]*matrix[1][0]*matrix[2][1]
                  - matrix[0][2]*matrix[1][1]*matrix[2][0]
                  - matrix[0][1]*matrix[1][0]*matrix[2][2]
                  - matrix[0][0]*matrix[1][2]*matrix[2][1];

return(determinant);
}

