#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20

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
  char path[200], input_fname[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  int  i, j, jnum, k, arrnum=200, halo, n, nsink, nsink_bin, type, snapshot_number, files, Ngas, random, ncount, pair_count, nsinkmax=20, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax, ke, pe, vesc, omega0, omega_crit;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, disCOM, disCOM1, disCOM2, vCOM, vCOM1, ncount_doub;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double posx[nsinkmax], posy[nsinkmax], posz[nsinkmax], velx[nsinkmax], vely[nsinkmax], velz[nsinkmax], sinkmass[nsinkmax], vrot_sink[nsinkmax], a[nsinkmax],mass1[nsinkmax],dis_cur[nsinkmax], mass2[nsinkmax], q[nsinkmax];
  double m_enc, m_enc_arr[arrnum], d_arr[arrnum];
  int sinkid1[nsinkmax], sinkid2[nsinkmax], sinkid3[nsinkmax], id1[nsinkmax], id2[nsinkmax], id3[nsinkmax], id4[nsinkmax];
  FILE *outfile;


  sprintf(path, "/nobackupp1/astacy/bin_zoom");

  halo = 6;

  if(halo == 1)
     {
     sprintf(basename, "bin_zoom1");
     sprintf(basenameout, "snapbin_zoom1");
     jnum = 80;
     }
  if(halo == 2)
     {
     sprintf(basename, "bin_zoom2");
     sprintf(basenameout, "snapbin_zoom2");
     jnum = 75;
     }
  if(halo == 3)
     {
     sprintf(basename, "bin_zoom3");
     sprintf(basenameout, "snapbin_zoom3");
     jnum = 80;
     }
  if(halo == 4)
     {
     sprintf(basename, "bin_zoom4");
     sprintf(basenameout, "snapbin_zoom4");
     jnum = 125;
     }
  if(halo == 5)
     {
     sprintf(basename, "bin_zoom5");
     sprintf(basenameout, "snapbin_zoom5");
     jnum = 110;
     }
  if(halo == 6)
     {
     sprintf(basename, "bin_zoom6");
     sprintf(basenameout, "snapbin_zoomb");
     jnum = 138;
     }
  if(halo == 7)
     {
     sprintf(basename, "bin_zoom7");
     sprintf(basenameout, "snapbin_zoom7");
     jnum = 88;
     }
  if(halo == 8)
     {
     sprintf(basename, "bin_zoom8");
     sprintf(basenameout, "snapbin_zoom8");
     jnum = 86;
     }
  if(halo == 9)
     {
     sprintf(basename, "bin_zoom9");
     sprintf(basenameout, "snapbin_zoom9");
     jnum = 109;
     }
  if(halo == 10)
     {
     sprintf(basename, "bin_zoom10");
     sprintf(basenameout, "snapbin_zoom10");
     jnum = 87;
     }
  
  //jnum = 8;

  printf("snapshot_number = %d \n", jnum);

  for(j=jnum;j<=jnum;j=j+10){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  printf("snapshot_number = %d \n", snapshot_number);

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(output_fname, "%s_pairs.dat",basename);
  sprintf(output_fname2, "%s_vel.dat",basename);
  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */

  unit_conversion();  

  nhmax = nsink = 0;
  mmax = mh_mass = 0;
  xCOM = yCOM = zCOM = vxCOM = vyCOM = vzCOM = ncount_doub = 0.0;
  xmax = ymax = zmax = 0;
  sl = tmax = h2max = gammin=masstot = dismax = omega0 = 0;
  typemax = 5;

  for(n=0;n<=nsinkmax-1;n++) 
    {
    posx[n] = 0;
    posy[n] = 0;
    posz[n] = 0;
    velx[n] = 0;
    vely[n] = 0;
    velz[n] = 0;
    sinkmass[n] = 0;
    mass1[n] = mass2[n] = vrot_sink[n] = a[n] = q[n] =  0.0;
    sinkid1[n] = sinkid2[n] = sinkid3[n] = id1[n] = id2[n] = id3[n] = id4[n] = 0;
    }

   for(n=0; n<arrnum; n++)
     {
     m_enc_arr[n] = 0.0;
     d_arr[n] = pow(5.e8,double(n)/double(arrnum));
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
          if(nh > nhmax/2.0)
            {
            vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
            vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
            vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
            xCOM = xCOM + P[n].Pos[0]*P[n].Mass;
            yCOM = yCOM + P[n].Pos[1]*P[n].Mass;
            zCOM = zCOM + P[n].Pos[2]*P[n].Mass;
            ncount_doub = ncount_doub + P[n].Mass;
            }
    }


         vx = vxCOM/ncount_doub;
         vy = vyCOM/ncount_doub;
         vz = vzCOM/ncount_doub;

         //sinkposx=xmax;
         //sinkposy=ymax;
         //sinkposz=zmax;

         sinkposx=xCOM/ncount_doub;
         sinkposy=yCOM/ncount_doub;
         sinkposz=zCOM/ncount_doub;

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

         for(k=0; k<arrnum; k++)
            printf("m_enc_arr = %lg\n", m_enc_arr[k]);
 
         //for(n = 1; n <= NumPart; n++)
         for(n = 1; n <= Ngas; n++)
             {

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
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


             if(P[n].nh > nhmax/2.0 && disAU > 0)
             //if(disAU <= 1.e1 /*&& P[n].sink < -4*/)
               {
               masstot = masstot + P[n].Mass/1.e-10/.7;
               omega0 = omega0 + (P[n].Mass/1.e-10/.7)*vrot/(disAU*1.5e8);
               }

             P[n].Vel[0] = P[n].Vel[0] + vx;
             P[n].Vel[1] = P[n].Vel[1] + vy;
             P[n].Vel[2] = P[n].Vel[2] + vz;



  }

   omega_crit = 4e-17*pow(nhmax/1.e3,.66666666);
   omega0 = omega0/masstot;
   printf("masstot = %lg, omega0 = %lg, omega_crit = %lg\n", masstot, omega0, omega_crit);
   ncount=nsink_bin=0;
   pair_count=0;
   for(i=0; i<nsinkmax; i++)
      {
      for(ncount=i+1; ncount<nsinkmax; ncount++)
          {
          if(sinkid1[i] == sinkid1[ncount])
            continue;
          if(sinkid1[i] == sinkid2[ncount])
            continue;
          if(sinkid2[i] == sinkid1[ncount])
            continue;
          if(sinkid2[i] == sinkid2[ncount] && sinkid2[i] != 0)
            continue;


          if(sinkid1[i] == 0)
            continue;
          if(sinkid1[ncount] == 0)
            continue;


          vx = velx[i] - velx[ncount];
          vy = vely[i] - vely[ncount];
          vz = velz[i] - velz[ncount]; 

          vel = pow(vx*vx + vy*vy + vz*vz,0.5);
          vel = vel*pow(Time, 0.5);
 
          dis = pow(((posx[i]-posx[ncount])*(posx[i]-posx[ncount]) + (posy[i]-posy[ncount])*(posy[i]-posy[ncount]) + (posz[i]-posz[ncount])*(posz[i]-posz[ncount])), 0.5);
          dis=dis*1.e3*Time/(0.7);
          disAU=dis*206264.806;

          vrad =  (vx*(posx[i]-posx[ncount]) + vy*(posy[i]-posy[ncount]) + vz*(posz[i]-posz[ncount]))/pow(((posx[i]-posx[ncount])*(posx[i]-posx[ncount]) + (posy[i]-posy[ncount])*(posy[i]-posy[ncount]) + (posz[i]-posz[ncount])*(posz[i]-posz[ncount])), 0.5);
          vrad = vrad*pow(Time, 0.5);

          vrotx = (posy[i]-posy[ncount])*vz  - (posz[i]-posz[ncount])*vy;
          vroty = (posz[i]-posz[ncount])*vx - (posx[i]-posx[ncount])*vz;
          vrotz = (posx[i]-posx[ncount])*vy  -  (posy[i]-posy[ncount])*vx;

          vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5);
          vrot = vrot*pow(Time, 0.5)/pow(((posx[i]-posx[ncount])*(posx[i]-posx[ncount]) + (posy[i]-posy[ncount])*(posy[i]-posy[ncount]) + (posz[i]-posz[ncount])*(posz[i]-posz[ncount])), 0.5);    ;  //convert to km/s

          pe = 6.67e-8*(sinkmass[ncount]+sinkmass[i])/(dis*3.0856802e18); 
          pe = pe*(1.9891e33);

          ke = 0.5*pow(1.e5*vel,2);

          xCOM = (posx[i]*sinkmass[i] + posx[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          yCOM = (posy[i]*sinkmass[i] + posy[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          zCOM = (posz[i]*sinkmass[i] + posz[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          vxCOM = (velx[i]*sinkmass[i] + velx[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          vyCOM = (vely[i]*sinkmass[i] + vely[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          vzCOM = (velz[i]*sinkmass[i] + velz[ncount]*sinkmass[ncount]) / (sinkmass[i] + sinkmass[ncount]);
          vCOM = pow(vx*vx + vy*vy + vz*vz,0.5);	
          vCOM1 = pow(vx*vx + vy*vy + vz*vz,0.5);
	  disCOM =  pow(((posx[i]-xCOM)*(posx[i]-xCOM) + (posy[i]-yCOM)*(posy[i]-yCOM) + (posz[i]-zCOM)*(posz[i]-zCOM)), 0.5); 
	  disCOM1 =  pow(((posx[ncount]-xCOM)*(posx[ncount]-xCOM) + (posy[ncount]-yCOM)*(posy[ncount]-yCOM) + (posz[ncount]-zCOM)*(posz[ncount]-zCOM)), 0.5);
			  
          disCOM=disCOM*1.e3*Time/(0.7);
          disCOM=disCOM*206264.806;
	  disCOM1=disCOM1*1.e3*Time/(0.7);
	  disCOM1=disCOM1*206264.806;
		  	  


          if(pe > ke)
             {
             printf("pair_count = %d\n", pair_count);
             id1[pair_count] = sinkid1[i]; 
             id2[pair_count] = sinkid2[i];
             id3[pair_count] = sinkid1[ncount];
             id4[pair_count] = sinkid2[ncount]; 
             //a[pair_count] = disAU;
	     a[pair_count] = 6.67e-8*2e33*(sinkmass[ncount]+sinkmass[i])/(2.*(pe-ke))/3.0856802e18*206264.806;
             dis_cur[pair_count] = disAU;
             mass1[pair_count] = sinkmass[i]; 
             mass2[pair_count] = sinkmass[ncount];
             vrot_sink[pair_count] = vrot; 
             pair_count++;

             nsink_bin++;
             if(sinkid2[i] == 0 && sinkid2[ncount] == 0)
               {
               posx[nsink+nsink_bin] = xCOM;
               posy[nsink+nsink_bin] = yCOM;
               posz[nsink+nsink_bin] = zCOM;
               velx[nsink+nsink_bin] = vxCOM;
               vely[nsink+nsink_bin] = vyCOM;
               velz[nsink+nsink_bin] = vzCOM;
               sinkmass[nsink+nsink_bin] = sinkmass[i]+sinkmass[ncount];
               sinkid1[nsink+nsink_bin] = sinkid1[i];
               sinkid2[nsink+nsink_bin] = sinkid1[ncount];
               }
             }
          }
      }



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






  











