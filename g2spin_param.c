#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int load_snapshot(char *fname, int files);
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
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  int Id;

  double  Rho, Temp, nh;
  double  U, HI, HII, HeI, HeII,  HeIII, H2II, HM, DII, DM;
  double elec, H2I, hsm, HDI, FosHII, sink, gam;
  double dummy;
} *P;

//int *Id;

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
  int  i, k=0, j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2;
  double x,y,z,x1,y1,z1, delr;
  double nh, nhmax, mass, mmax, dis, disAU, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  int IDmax;
  double rho_b, rho_DM[200], rho_vir, r_vir, v_vir;
  double kB, G;
  double massHIT, diskmass, massin=0.0, massout=0.0;
  double mass_sink1, mass_sink2, mass_sink3, mass_sink4, mass_tot;
  double dis12, dis13, x_com, y_com, z_com, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2, sinkposx3, sinkposy3, sinkposz3, sinkposx4, sinkposy4, sinkposz4;
  double num, numz, tempz, vz, rhoz, amom, amomtot, amomtotx, amomtoty, amomtotz, amomx, amomy, amomz, amomx1, amomx2, amomy1, amomy2, amomz1, amomz2;
  double xfac, yfac, zfac, kfloat, nh_avg, num_tot;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, jzcent;
  double m_enc_arr[200], d_arr[200], mass_enc;
  double KE, KE_tot, PE, PE_tot, vel, spin_param;
  int ID_sink1, ID_sink2, ID_sink3, ID_sink4, sinknumk;
  FILE *outfile;


  sprintf(path, "/work/utexas/ao/minerva/bf");
  sprintf(basename, "bf10");

  for(j=0;j<=58;j++){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  Ngas = load_snapshot(input_fname, files);

          /*    reordering();*/ /* call this routine only if your ID's are set properly */

  rho_b = 2.94e-30; //background density of the universe

  ID_sink1 = 3755078;
  ID_sink2 = 4333881;
  ID_sink3 = 3668992;
  ID_sink3 = 4653020;

          unit_conversion();

          nh_avg=0;
          num_tot=0;
          h2max=0;
          nhmax=0;
          mmax=0;
          kfloat=0.0;
          xtimesm = 0.0;
          ytimesm = 0.0;
          ztimesm = 0.0;
          vxtimesm = 0.0;
          vytimesm = 0.0;
          vztimesm = 0.0;
          diskmass = 0;

          amomtot = 0.0;
          PE_tot = 0.0;
          KE_tot = 0.0;


          for(i = 0; i < Ngas; i++)
              {
                if(P[i].nh > nhmax && P[i].sink < 0.5)
                  {
                    nhmax = P[i].nh;
                    IDmax = P[i].Id;
                    xmax = P[i].Pos[0];
                    ymax = P[i].Pos[1];
                    zmax = P[i].Pos[2];
                  }
              }
           printf("xmax = %lg ymax = %lg zmax = %lg\n", xmax, ymax, zmax);

           for(i = 0; i < Ngas; i++)
              {
               dis = pow(((P[i].Pos[0]-xmax)*(P[i].Pos[0]-xmax) + (P[i].Pos[1]-ymax)*(P[i].Pos[1]-ymax) + (P[i].Pos[2]-zmax)*(P[i].Pos[2]-zmax)), 0.5);
               dis=dis*1.e3*Time/(0.7);                   //dis is in pc
              //if(P[i].Id == ID_sink1)
               if(P[i].Id == IDmax)
                 {
                 sinkposx1 = P[i].Pos[0];
                 sinkposy1 = P[i].Pos[1];
                 sinkposz1 = P[i].Pos[2];
                 vx1 = P[i].Vel[0];
                 vy1 = P[i].Vel[1];
                 vz1 = P[i].Vel[2];
                 mass_sink1 = P[i].Mass*1.e10/0.7;
                 printf("Found the sink! mass= %lg\n", mass_sink1);
                 }
  
                if(P[i].Id == ID_sink2)
                  {
                  mass_sink2 = P[i].Mass*1.e10/0.7;
                  sinkposx2 = P[i].Pos[0];
                  sinkposy2 = P[i].Pos[1];
                  sinkposz2 = P[i].Pos[2];
                  }
                             
                 if(/*P[i].nh > 1.e3 &&*/ dis < 200.0)
                  {
                  diskmass = diskmass + (P[i].Mass*1.e10/0.7);
                  xtimesm = xtimesm + P[i].Pos[0]*(P[i].Mass*1.e10/0.7);
                  ytimesm = ytimesm + P[i].Pos[1]*(P[i].Mass*1.e10/0.7);
                  ztimesm = ztimesm + P[i].Pos[2]*(P[i].Mass*1.e10/0.7);
                  vxtimesm = vxtimesm + P[i].Vel[0]*(P[i].Mass*1.e10/0.7);
                  vytimesm = vytimesm + P[i].Vel[1]*(P[i].Mass*1.e10/0.7);
                  vztimesm = vztimesm + P[i].Vel[2]*(P[i].Mass*1.e10/0.7);
                  }
                
  /*     
                  if(P[i].sink > 0.5 )
                  //if(P[i].Id == ID_sink1 || P[i].Id == ID_sink2)
                  {
                  xtimesm = xtimesm + P[i].Pos[0]*(P[i].Mass*1.e10/0.7);
                  ytimesm = ytimesm + P[i].Pos[1]*(P[i].Mass*1.e10/0.7);
                  ztimesm = ztimesm + P[i].Pos[2]*(P[i].Mass*1.e10/0.7);
                  vxtimesm = vxtimesm + P[i].Vel[0]*(P[i].Mass*1.e10/0.7);
                  vytimesm = vytimesm + P[i].Vel[1]*(P[i].Mass*1.e10/0.7);
                  vztimesm = vztimesm + P[i].Vel[2]*(P[i].Mass*1.e10/0.7);
                  diskmass = diskmass + (P[i].Mass*1.e10/0.7);
                  }
    */    
             }
        
    
             xCOM = xtimesm/diskmass;
             yCOM = ytimesm/diskmass;
             zCOM = ztimesm/diskmass;
             vxCOM = vxtimesm/diskmass;
             vyCOM = vytimesm/diskmass;
             vzCOM = vztimesm/diskmass;
             printf("%15.11g %15.11g %15.11g\n", xCOM, yCOM, zCOM); 
             printf("%15.11g %15.11g %15.11g\n", vxCOM, vyCOM, vzCOM);            

	
          //  if(j>=4)
          //  {
            sinkposx1 = xCOM;
            sinkposy1 = yCOM;
            sinkposz1 = zCOM;
            mass_sink1 = diskmass;
            vx1 = vxCOM;
            vy1 = vyCOM;
            vz1 = vzCOM;
         //   }
	 

          for(i = 0; i < Ngas; i++)
              {

              }

         diskmass=0.0;
         ncount=0;  

        for(n=0; n<200; n++)
               {
               m_enc_arr[n] = 0.0;
               rho_DM[n] = 0.0;
               d_arr[n] = (double(n)/199.0)*5.e8;
               //d_arr[n] = pow(1.e12,(double(n)/200.0));
               printf("d_arr[%d] = %lg\n", n, d_arr[n]);
               }

        for(i = Ngas; i < NumPart; i++)
               {
                   dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                   dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                   disAU = dis*206264.8060;
                  if(disAU < 5.e8)
                   {
                   for(n=0; n<200; n++)
                     {
                     if(disAU < d_arr[n] /*|| P[i].sink > 0.5*/)
                       m_enc_arr[n] = m_enc_arr[n] + P[i].Mass*1.e10/0.7;
                     }
                   }
                  }

        for(n=0; n<200; n++)
               {
               rho_DM[n] = m_enc_arr[n]*1.98892e33/((4.0/3.0)*3.14159*pow(d_arr[n]*1.5e13, 3.0));  //DM density in cgs
               printf("m_enc_arr[%d] = %lg\n", n, m_enc_arr[n]);
               printf("rho_DM[%d] = %lg\n", n, rho_DM[n]/(rho_b*pow(Time, -3.0)));
               if(rho_DM[n] > 200.0*rho_b*pow(Time, -3.0))
                 {
                 rho_vir = rho_DM[n];
                 r_vir = d_arr[n];
                 }
               }


         for(i = Ngas; i < NumPart; i++)
                {

                 P[i].Vel[0] = P[i].Vel[0] - vx1;
                 P[i].Vel[1] = P[i].Vel[1] - vy1;
                 P[i].Vel[2] = P[i].Vel[2] - vz1;
          
                 amomx = (P[i].Pos[1]-sinkposy1)*P[i].Vel[2]  - (P[i].Pos[2]-sinkposz1)*P[i].Vel[1];
                 amomy = (P[i].Pos[2]-sinkposz1)*P[i].Vel[0] - (P[i].Pos[0]-sinkposx1)*P[i].Vel[2];
                 amomz = (P[i].Pos[0]-sinkposx1)*P[i].Vel[1]  -  (P[i].Pos[1]-sinkposy1)*P[i].Vel[0];

                 amomx =  amomx*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);  //convert to cgs units
                 amomy =  amomy*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);
                 amomz =  amomz*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);

                 amom = (P[i].Mass*1.e10*1.98892e33/.7)*pow(amomx*amomx + amomy*amomy + amomz*amomz, 0.5);  //ang. mom. of particle in cgs

                 vel = pow(P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2], 0.5);
                 vel = vel*pow(Time, 0.5)*1.e5;  //convert vel to cgs 

                 KE = 0.5*(P[i].Mass*1.e10*1.98892e33/.7)*pow(vel,2.0);   //kinetic energy in cgs

                 dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5); 
                 //dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                 dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                 disAU = dis*206264.8060;
                 dis=dis*3.086e18;   //dis in now in cm                
 
                 if(disAU < r_vir && disAU > 0.0)
                 {
                   mass_enc = m_enc_arr[0];
                   for(n=0; n<200; n++)
                     { 
                     if(disAU > d_arr[n])
                       mass_enc = m_enc_arr[n];
                     }
                 jzcent = pow(6.67e-8 * mass_enc * 1.98892e33*dis, 0.5);
                 PE = -6.67e-8*(P[i].Mass*1.e10*1.98892e33/.7)*mass_enc*1.98892e33/dis;  //binding energy of particle in cgs
                 }

                 if(disAU < r_vir /*&& KE + PE < 0.0*/)
                  {
                  KE_tot = KE_tot + KE;
                  PE_tot = PE_tot + PE;
                  amomtot = amomtot + amom;
                  num_tot = num_tot + P[i].Mass;

                    ncount++;  
                  }
                }


            spin_param = amomtot*pow(double(fabs(KE_tot+PE_tot)),0.5)/(6.67e-8*pow(num_tot*1.e10*1.98892e33/.7, 2.5));
            v_vir = 6.67e-8*(num_tot*1.e10/.7)*2.e33/(r_vir*1.5e13); 
            v_vir = pow(v_vir,0.5);

             outfile=fopen("bf10_spin", "a");             
             fprintf(outfile, "%15.11g %8d %15.11g %15.11g %15.11g %15.11g  %15.11g %15.11g  %15.11g %15.11g %15.11g\n", Time, nhmax, IDmax, KE_tot, PE_tot, amomtot, num_tot*1.e10/.7,  rho_vir/(rho_b*pow(Time, -3.0)), r_vir, spin_param, v_vir);
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

	  u  = P[i].Temp * UnitEnergy_in_cgs/ UnitMass_in_g;

	  //gamma= 5.0/3.0;
          gamma=P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs;
	  P[i].nh= P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
          P[i].Rho = P[i].Rho * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
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
/*      printf("N_DM= %6d \n",NumPart-Ngas); */

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
  printf("Time= %12.10e \n",Time);
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
  
  //P--;   /* start with offset 1 */

  /*
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  //Id--;   // start with offset 1 
*/
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

  for(i=0; i<NumPart; i++)
    {
      if(P[i].Id != i)
	{
	  psource= P[i];
	  idsource=P[i].Id;
	  dest=P[i].Id;

	  do
	    {
	      psave= P[dest];
	      idsave=P[dest].Id;

	      P[dest]= psource;
	      P[dest].Id= idsource;
	      
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

  //Id++;   
  //free(Id);

  printf("space for particle ID freed\n");
}






  











