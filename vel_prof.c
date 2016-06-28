#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20

//#define width 1.0
#define width 0.2

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
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;

  double  Rho, U, Temp, nh;
  double H2I, HII, HDI, HeII, gam, sink, c_s;
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
  char path[200], input_fname[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  int  i, j, jnum, k, arrnum=200, halo, n, nsink, nsink_bin, type, snapshot_number, files, Ngas, random, ncount, pair_count, nsinkmax=20, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax, ke, pe, vesc, omega0, omega_crit;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, disCOM, disCOM1, disCOM2, vCOM, vCOM1, ncount_doub;
  double nh, nhmax, nthresh, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, dis_pc, disAU, vrad, vradx, vrady, vradz, vrotx, vroty, vrotz, vrot, mh_mass;
  double posx[nsinkmax], posy[nsinkmax], posz[nsinkmax], velx, vely, velz, sinkmass[nsinkmax], vrot_sink[nsinkmax], a[nsinkmax],mass1[nsinkmax],dis_cur[nsinkmax], mass2[nsinkmax], q[nsinkmax];
  double m_enc, m_enc_arr[arrnum], d_arr[arrnum], num_arr[arrnum], turb, turb_arr[arrnum], c_s_arr[arrnum], vrad_arr[arrnum], vrot_arr[arrnum];
  double vel_arr[arrnum], velx_arr[arrnum], vely_arr[arrnum], velz_arr[arrnum], vradx_arr[arrnum], vrady_arr[arrnum], vradz_arr[arrnum]; 
  double vrotx_arr[arrnum], vroty_arr[arrnum], vrotz_arr[arrnum], amom_arr[arrnum];
  double rmin, rmax, disx, disy, disz, hubble;
  FILE *outfile;


  sprintf(path, "/work/00863/minerva/bin_map");

  halo = 10;

  if(halo == 1)
     {
     sprintf(basename, "hires_test10");
     sprintf(basenameout, "snaphires_test10");
     jnum = 1;
     }
  if(halo == 2)
     {
     sprintf(path, "/nobackupp1/astacy/midhires");
     sprintf(basename, "midhires");
     sprintf(basenameout, "snapmidhires");
     jnum = 800;
     }
  if(halo == 3)
     {
     sprintf(basename, "hires3");
     sprintf(basenameout, "hires3");
     jnum = 1;
     }
  if(halo == 4)
     {
     sprintf(basename, "dma_alt");
     sprintf(basenameout, "dma_alt");
     jnum = 37;
     }
  if(halo == 5)
     {
     sprintf(basename, "dmaH100");
     sprintf(basenameout, "dmaH100");
     jnum = 307;
     }
  if(halo == 6)
     {
     sprintf(basename, "dmaH60");
     sprintf(basenameout, "dmaH60");
     jnum = 50;
     }
  if(halo == 7)
     {
     sprintf(basename, "dma_prof1_ng");
     sprintf(basenameout, "dma_prof1_ng");
     jnum = 1;
     }
  if(halo == 8)
     {
     sprintf(basename, "dmaH101");
     sprintf(basenameout, "dmaH101");
     jnum = 496;
     }
  if(halo == 9)
     {
     sprintf(basename, "bin_HR10_map");
     sprintf(basenameout, "bin_HR10");
     jnum = 0;
     //jnum = 10;
     //jnum = 19;  //(2e11 seconds)
     }
  if(halo == 10)
     {
     sprintf(path, "/work/00863/minerva/");
     sprintf(basename, "bin_HR10");
     sprintf(basenameout, "bin_HR10");
     jnum = 7;
     }

  printf("snapshot_number = %d \n", jnum);

  for(j=jnum;j<=jnum;j=j+10){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  printf("snapshot_number = %d \n", snapshot_number);

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s_vel_%04d.dat",basename, snapshot_number);
  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */

  unit_conversion();  

  nhmax = nsink = 0;
  mmax = mh_mass = 0;
  xCOM = yCOM = zCOM = vxCOM = vyCOM = vzCOM = ncount_doub = 0.0;
  xmax = ymax = zmax = 0;
  sl = tmax = h2max = gammin=masstot = dismax = omega0 = 0;
  typemax = 5;
  nthresh = 2.e11;

   rmin = 100;
   rmax = 3e5;
   for(n=0; n<arrnum; n++)
     {
     m_enc_arr[n] = num_arr[n] = turb_arr[n] = c_s_arr[n] = vrad_arr[n] = vrot_arr[n] = velx_arr[n] = vely_arr[n] = velz_arr[n] = 0.0;
     //d_arr[n] = pow(5.e8,double(n)/double(arrnum));
     d_arr[n] = rmin + (rmax - rmin)*double(n)/double(arrnum);
     //d_arr[n] = rmin * pow(rmax/rmin, double(n)/double(arrnum));
     }


  for(n=0;n<=Ngas;n++) { 

          nh = P[n].nh;

          if(nh > nhmax)
          //if(mass > mmax)
          //if(Id[n] > idmax)
          //  if(Id[n] == 994680) //binHR5
            {
              //printf("Found the sink!\n");
              nhmax = nh;
              mmax=P[n].Mass;
              xmax=P[n].Pos[0];
              ymax=P[n].Pos[1];
              zmax=P[n].Pos[2];
              typemax = P[n].Type;
              idmax = Id[n];
              vx = P[n].Vel[0];
              vy = P[n].Vel[1];
              vz = P[n].Vel[2];
            }
    }

  for(n=0;n<=Ngas;n++) {
            nh = P[n].nh;
            dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
            dis=dis*1.e3*Time/(0.7);
            disx = fabs((P[n].Pos[0]-xmax))*1.e3*Time/(0.7);
            disy = fabs((P[n].Pos[1]-ymax))*1.e3*Time/(0.7);
            disz = fabs((P[n].Pos[2]-zmax))*1.e3*Time/(0.7);
            disAU=dis*206264.806;
            //if(disAU < 400.)
            //if(dis < width)
            if(fabs(disx) < width/2.0 && fabs(disy) < width/2.0 && fabs(disz) < width/2.0)
            //if(P[n].nh > nhmax/10.0)
            {
            vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
            vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
            vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
            xCOM = xCOM + P[n].Pos[0]*P[n].Mass;
            yCOM = yCOM + P[n].Pos[1]*P[n].Mass;
            zCOM = zCOM + P[n].Pos[2]*P[n].Mass;
            ncount_doub = ncount_doub + P[n].Mass;
            }
/*
     if(disAU < 10000.)
       {
       P[n].Vel[0] = 2./pow(Time, 0.5);
       P[n].Vel[1] = 3./pow(Time, 0.5);
       P[n].Vel[2] = 4./pow(Time, 0.5);
       }
*/
    }


         vx = vxCOM/ncount_doub;
         vy = vyCOM/ncount_doub;
         vz = vzCOM/ncount_doub;
         //vx = vy = vz = 0;

         sinkposx=xmax;
         sinkposy=ymax;
         sinkposz=zmax;

         //sinkposx=xCOM/ncount_doub;
         //sinkposy=yCOM/ncount_doub;
         //sinkposz=zCOM/ncount_doub;

         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g, nhmax = %15.11g, mmax = %15.11g\n", sinkposx, sinkposy, sinkposz, nhmax, mmax);
 
         if(nhmax < 1.e1)
           {
           sinkposx=header1.BoxSize/2.0;
           sinkposy=header1.BoxSize/2.0;
           sinkposz=header1.BoxSize/2.0; 
           }
         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g\n", sinkposx, sinkposy, sinkposz);
         printf("vx = %lg, vy = %lg, vz = %lg\n", vx, vy, vz);

        for(i = 0; i <= Ngas; i++)
             {
             dis = pow(((P[i].Pos[0]-sinkposx)*(P[i].Pos[0]-sinkposx) + (P[i].Pos[1]-sinkposy)*(P[i].Pos[1]-sinkposy) + (P[i].Pos[2]-sinkposz)*(P[i].Pos[2]-sinkposz)), 0.5);
             dis=dis*1.e3*Time/(0.7);
             disAU=dis*206264.806;

             if(disAU < 5.e8)
               {
               for(n=0; n<arrnum; n++)
                  {
                  if(disAU < d_arr[n] /*|| P[i].sink > 0.5*/)
                     m_enc_arr[n] = m_enc_arr[n] + P[i].Mass*1.e10/0.7;
                  }
               }

             if(disAU < 10.0)
               {
               mh_mass = mh_mass + P[i].Mass*1.e10/0.7;
               }
             }
         printf("mh_mass = %lg\n", mh_mass);

         //for(k=0; k<arrnum; k++)
         //   printf("m_enc_arr = %lg\n", m_enc_arr[k]);
 
         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n <= Ngas; n++)
             {

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
             dis_pc=dis*1.e3*Time/(0.7);
             disAU=dis_pc*206264.806;

             vel = P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1] + P[n].Vel[2]*P[n].Vel[2];
             vel = pow(vel,0.5)*pow(Time, 0.5); 

             velx = P[n].Vel[0]*pow(Time, 0.5);
             vely = P[n].Vel[1]*pow(Time, 0.5);
             velz = P[n].Vel[2]*pow(Time, 0.5);

             vrad =  (P[n].Vel[0]*(P[n].Pos[0]-sinkposx) + P[n].Vel[1]*(P[n].Pos[1]-sinkposy) + P[n].Vel[2]*(P[n].Pos[2]-sinkposz))/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             vrad = vrad*pow(Time, 0.5); 

             vradx = vrad*(P[n].Pos[0]-sinkposx)/dis;
             vrady = vrad*(P[n].Pos[1]-sinkposy)/dis;
             vradz = vrad*(P[n].Pos[2]-sinkposz)/dis;

             vrotx = (P[n].Pos[1]-sinkposy)*P[n].Vel[2]  - (P[n].Pos[2]-sinkposz)*P[n].Vel[1];
             vroty = (P[n].Pos[2]-sinkposz)*P[n].Vel[0] - (P[n].Pos[0]-sinkposx)*P[n].Vel[2];
             vrotz = (P[n].Pos[0]-sinkposx)*P[n].Vel[1]  -  (P[n].Pos[1]-sinkposy)*P[n].Vel[0];

             vrotx = vrotx*pow(Time, 0.5)/dis;   //convert to km/s
             vroty = vroty*pow(Time, 0.5)/dis;
             vrotz = vrotz*pow(Time, 0.5)/dis;  

             vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5); 

             if(P[n].nh > nthresh && disAU > 0)
             //if(disAU <= 1.e1 /*&& P[n].sink < -4*/)
               {
               masstot = masstot + P[n].Mass/1.e-10/.7;
               omega0 = omega0 + (P[n].Mass/1.e-10/.7)*vrot/(disAU*1.5e8);
               }

             if(P[n].nh > 1.e-2 && disAU > 0)
             //if(disAU <= 1.e1 /*&& P[n].sink < -4*/)
               {
               for(k=0; k<arrnum; k++)
                  if(disAU < d_arr[k])
                     m_enc = m_enc_arr[k];
               for(k=0; k<arrnum; k++)
                  if(disAU < d_arr[k+1] && disAU > d_arr[k])
                  //if(disAU < d_arr[k] && disAU > d_arr[k-1])
                     {
                     num_arr[k] = num_arr[k] + 1.0; 
                     c_s_arr[k] = c_s_arr[k] + P[n].c_s/1.e5;
                     vel_arr[k] = vel_arr[k] + vel;
                     velx_arr[k] = velx_arr[k] + velx;
                     vely_arr[k] = vely_arr[k] + vely;
                     velz_arr[k] = velz_arr[k] + velz;
                     vrot_arr[k] = vrot_arr[k] + vrot;
                     vrotx_arr[k] = vrotx_arr[k] + vrotx;
                     vroty_arr[k] = vroty_arr[k] + vroty;
                     vrotz_arr[k] = vrotz_arr[k] + vrotz;
                     vrad_arr[k] = vrad_arr[k] + vrad;
                     vradx_arr[k] = vradx_arr[k] + vradx;
                     vrady_arr[k] = vrady_arr[k] + vrady;
                     vradz_arr[k] = vradz_arr[k] + vradz;
                     }
               }


            }

		  	  
   outfile=fopen(output_fname, "a");
   for(n=0; n<=arrnum; n++)
      {
      c_s_arr[n] = c_s_arr[n]/num_arr[n];
      vel_arr[n] = vel_arr[n]/num_arr[n];
      velx_arr[n] = velx_arr[n]/num_arr[n];
      vely_arr[n] = vely_arr[n]/num_arr[n];
      velz_arr[n] = velz_arr[n]/num_arr[n];
      vrot_arr[n] = vrot_arr[n]/num_arr[n];
      vrotx_arr[n] = vrotx_arr[n]/num_arr[n];
      vroty_arr[n] = vroty_arr[n]/num_arr[n];
      vrotz_arr[n] = vrotz_arr[n]/num_arr[n];
      vrad_arr[n] = vrad_arr[n]/num_arr[n];
      vradx_arr[n] = vradx_arr[n]/num_arr[n];
      vrady_arr[n] = vrady_arr[n]/num_arr[n];
      vradz_arr[n] = vradz_arr[n]/num_arr[n];

      turb_arr[n] = pow(velx_arr[n] - vradx_arr[n] - vrotx_arr[n],2) + pow(vely_arr[n] - vrady_arr[n] - vroty_arr[n],2) + pow(velz_arr[n] - vradz_arr[n] - vrotz_arr[n],2);
      turb_arr[n] = pow(turb_arr[n], 0.5);

      vrot_arr[n] = pow(vrotx_arr[n]*vrotx_arr[n] + vroty_arr[n]*vroty_arr[n] + vrotz_arr[n]*vrotz_arr[n], 0.5);
      if(n > 0)
        amom_arr[n] = (m_enc_arr[k] - m_enc_arr[k-1])*vrot_arr[n]*d_arr[n]*1.5e13*1.e5*2.e33; 
      if(n == 0)
        amom_arr[n] = m_enc_arr[0]*vrot_arr[0]*d_arr[0]*1.5e13*1.e5*2.e33;

      //turb_arr[n] = turb_arr[n]/num_arr[n];
      fprintf(outfile, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n",  d_arr[n], c_s_arr[n], turb_arr[n], vrad_arr[n], vrot_arr[n], velx_arr[n], vely_arr[n], velz_arr[n], vradx_arr[n], vrady_arr[n], vradz_arr[n], vrotx_arr[n], vroty_arr[n], vrotz_arr[n], amom_arr[n], m_enc_arr[n]);
      }
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
  double h2frac, muh2, muh2in, pressure;
 
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

          //MeanWeight=MeanWeight*PROTONMASS;

	  MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  gamma= 5.0/3.0;
          //gamma=P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs* HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  P[i].nh= P[i].Rho / MeanWeight;
          pressure = P[i].nh * BOLTZMANN * P[i].Temp;
          //P[i].c_s = pow(1.38e-16*P[i].Temp*gamma/MeanWeight,0.5);
          P[i].c_s = pow(gamma*pressure/P[i].Rho, 0.5);
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
/*
      for(k=0, NumPart=0, ntot_withmasses=0; k<MAXREF; k++)
        printf("header1.npartTotal[%d] = %d\n", k, header1.npartTotal[k]);
      for(k=0, NumPart=0, ntot_withmasses=0; k<MAXREF; k++)
        printf("header1.npart[%d] = %d\n", k, header1.npart[k]);
*/
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
	      fread(&Id[pc_new], sizeof(int), 1, fd);
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
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd); //hsm
              pc_sph++;
            }
          SKIP;


          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
              fread(&P[pc_sph].HII, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);   
              fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
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

  
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
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






  











