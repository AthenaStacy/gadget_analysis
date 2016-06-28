#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define Hubble 3.2407789e-18

int load_snapshot(char *fname, int files);
int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(int N_part);
int allocate_memory2(int N_part);
int allocate_memory3(int N_part);

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


int     NumPart, Ngas, N_gas, N_gas3;

struct particle_data 
{
    double Pos[3];
    double Vel[3];
    double Mass;
    int Type;
    double EgySpec;
    int Id;
    double sink;
    double nh;
    double Temp;
    double gam;
    double dummy;

} *P;

struct particle_data2 
{
    double Pos[3];
    double Vel[3];
    double Mass;
    int Type;
    double EgySpec;
    double t0;
    int Id;
    double pres_x;
    double pres_y;
    double pres_z;
    double visc_x;
    double visc_y;
    double visc_z;
    double grav_x;
    double grav_y;
    double grav_z;  
    double Hsml;
    double sink;
    double nh;
    double Temp;
    double gam;
    double dummy;
    double H2I;
    double HII;
} *P2;

struct particle_data3
{
    double Pos[3];
    double nHcgs;
    double lamH2;
    double lamHl;
    double Temp;
    double dummy;
    int Id;
} *P3;


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
  int  i, k, j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2;
  double x,y,z,x1,y1,z1, delr;
  double nh, nhmax, mass, mmax, dis, disAU, disxy, disxz, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double kB, G;
  double tgrav, tpres, tvisc, tgravx, tgravy, tgravz, tpresx, tpresy, tpresz, tviscx, tviscy, tviscz;
  double mass_sink1, mass_sink2, mass_sink3, mass_sink4, mass_tot, mass_part;
  double dis12, dis13, x_com, y_com, z_com, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double ax1, ay1, az1;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2, sinkposx3, sinkposy3, sinkposz3, sinkposx4, sinkposy4, sinkposz4;
  double numz, tempz, vz, rhoz, amomx, amomy, amomz, amom, amomx2, amomy1, amomy2, amomz1, amomz2;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double axtimesm, aytimesm, aztimesm;
  double kin, pot;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, axCOM, ayCOM, azCOM, jzcent, diskmass;
  double tgrav_IN, tgrav_MID, tgrav_OUT, tpres_IN, tpres_OUT, tvisc_IN, tvisc_OUT, num_IN, num_MID, num_OUT;
  double  omega, cs_squared_S1, h, pres, talpha, ttot;
  int ID_sink1, ID_sink2, ID_sink3, ID_sink4;
  int weight;
  FILE *infile, *infile3;
  FILE *outfile;
  FILE *outfile2, *outfile3;
  double m_enc_arr[200], d_arr[200], angmom_arr[199], n_arr[199], mass_enc, mass_enc_cent;
  double hubble_a, AllHub, OmegaLambda, dtconv, tcool;
  double UnitEnergy_in_cgs, UnitTime_in_s, UnitMass_in_g, UnitVelocity_in_cm_per_s, UnitLength_in_cm, UnitDensity_in_cgs;
  double Omega0=0.3;
  double HubbleParam=0.7;

  UnitLength_in_cm= 3.085678e21;
  UnitVelocity_in_cm_per_s= 1.0e5;
  UnitMass_in_g= 1.989e43;
  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
  OmegaLambda = 0.7;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);

  sprintf(path, "/work/utexas/ao/minerva");
  sprintf(basename, "tgridsink_am4");

  for(j=120;j<=120;j++){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  kB = 1.38e-16;
  G = 6.67e-8;
  diskmass = 0.0;
  tgrav = 0.0;
  tpres = 0.0;
  tvisc = 0.0;

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  //Ngas = load_snapshot(input_fname, files);

  ID_sink1 = 3755078;
  ID_sink2 = 4025630;
  ID_sink3 = 4001858;
  ID_sink3 = 4653020;

	  /*    reordering();*/ /* call this routine only if your ID's are set properly */  
	 
	  nhmax=0; 
	  tmax=0;
	  h2max=0;
	  masstot=0;
	  mmax=0;

          xtimesm=0; 
          ytimesm=0; 
          ztimesm=0; 
          vxtimesm=0; 
          vytimesm=0; 
          vztimesm=0;
          axtimesm=0; 
          aytimesm=0; 
          aztimesm=0;

          tgrav_IN=0.0;
          tpres_IN=0.0;
          tvisc_IN=0.0;
          tgrav_OUT=0.0;
          tpres_OUT=0.0;
          tvisc_OUT=0.0;
          num_IN=0.0;
          num_OUT=0.0;

          diskmass=0.0;
          free(P);

          N_gas = 26833;  //_120
          //N_gas = 29302;  //004   

	  N_gas3 = 429336;  //_120
          //N_gas3 = 468832;  //004

          infile = fopen("/work/utexas/ao/minerva/torque_120.dat", "r");

          for(i = 0; i < N_gas; i++)
             {
             if(i==0)
              allocate_memory2(N_gas);
 
             fscanf(infile, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",   
                        &P2[i].t0, &P2[i].dummy, &P2[i].Mass,
                        &P2[i].Pos[0], &P2[i].Pos[1], &P2[i].Pos[2],
                        &P2[i].nh,
                        &P2[i].pres_x,
                        &P2[i].pres_y,
                        &P2[i].pres_z,
                        &P2[i].visc_x,
                        &P2[i].visc_y,
                        &P2[i].visc_z,
                        &P2[i].grav_x,
                        &P2[i].grav_y,
                        &P2[i].grav_z,
                        &P2[i].dummy,
                        &P2[i].Id,
                        &P2[i].Hsml,
                        &P2[i].sink,
                        &P2[i].gam,
                        &P2[i].Temp,
                        &P2[i].Vel[0],
                        &P2[i].Vel[1],
                        &P2[i].Vel[2],
                        &P2[i].H2I,
                        &P2[i].HII); 
              }
           fclose(infile);

              Time = P2[0].t0;

              printf("ID_sink1 = %d, Time = %lg \n",ID_sink1, Time);
              printf("ID[25000] = %d\n", P2[0].Id);

               for(i = 0; i < N_gas; i++)
                {
                if(P2[i].Id == ID_sink2)
                  {
                  sinkposx1 = P2[i].Pos[0];
                  sinkposy1 = P2[i].Pos[1];
                  sinkposz1 = P2[i].Pos[2];
                  vx1 = P2[i].Vel[0];
                  vy1 = P2[i].Vel[1];
                  vz1 = P2[i].Vel[2];
                  mass_sink1 = P2[i].Mass*1.e10/0.7;
                  printf("Found the sink! mass= %lg\n", mass_sink1);
                  }
                }

        for(n=0; n<200; n++)
               {
               m_enc_arr[n] = 0.0;
               d_arr[n] = 2.e4*(double(n+1)/199.0);
               printf("d_arr[%d] = %lg\n", n, d_arr[n]);
               }


        for(n=0; n<199; n++)
               {
               angmom_arr[n] = 0.0;
               n_arr[n]=0.0;
               }

        for(i = 0; i < N_gas; i++)
               {
                if(P2[i].nh > 1.e1)
                  {
                   dis = pow(((P2[i].Pos[0]-sinkposx1)*(P2[i].Pos[0]-sinkposx1) + (P2[i].Pos[1]-sinkposy1)*(P2[i].Pos[1]-sinkposy1) + (P2[i].Pos[2]-sinkposz1)*(P2[i].Pos[2]-sinkposz1)), 0.5);
                   dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                   disAU = dis*206264.8060;
                   for(n=0; n<200; n++)
                     {
                     if(disAU < d_arr[n] /*|| P2[i].sink > 0.5*/ /*&& P2[i].Id != ID_sink1*/ /*|| P2[i].ID == ID_sink2*/)
                       m_enc_arr[n] = m_enc_arr[n] + P2[i].Mass*1.e10/0.7;
                     }
                    }
                  }

        for(n=0; n<200; n++)
               {
               printf("m_enc_arr[%d] = %lg\n", n, m_enc_arr[n]);
               }


           printf("point 1\n");  

           AllHub = Hubble*UnitTime_in_s;
           hubble_a = AllHub * sqrt(Omega0 / (Time * Time * Time)
                                   + (1 - Omega0 - OmegaLambda) / (Time * Time) +
                                   OmegaLambda);
           dtconv = (1.0/(hubble_a*Time))*UnitTime_in_s/HubbleParam;

           printf("point 2\n");

          outfile2=fopen("/work/utexas/ao/minerva/torque_am2_120", "a");

          printf("point 3\n");

          for(i = 0; i < N_gas; i++)
               {

                P2[i].Vel[0] = P2[i].Vel[0] - vx1;
                P2[i].Vel[1] = P2[i].Vel[1] - vy1;
                P2[i].Vel[2] = P2[i].Vel[2] - vz1;
 
                amomx = (P2[i].Pos[1]-sinkposy1)*P2[i].Vel[2]  - (P2[i].Pos[2]-sinkposz1)*P2[i].Vel[1];
                amomy = (P2[i].Pos[2]-sinkposz1)*P2[i].Vel[0] - (P2[i].Pos[0]-sinkposx1)*P2[i].Vel[2];
                amomz = (P2[i].Pos[0]-sinkposx1)*P2[i].Vel[1]  -  (P2[i].Pos[1]-sinkposy1)*P2[i].Vel[0];

                amomx =  amomx*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);  //convert to cgs units
                amomy =  amomy*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);
                amomz =  amomz*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);
 
                amom = pow(amomx*amomx + amomy*amomy + amomz*amomz, 0.5);

		 //printf("Time again = %lg", Time);
                P2[i].pres_x = (P2[i].pres_x - ax1)*pow(Time, 0.5)*1.e5/dtconv;  //convert from code units to cm/s per time
                P2[i].pres_y = (P2[i].pres_y - ay1)*pow(Time, 0.5)*1.e5/dtconv;  //also convert from cm/s per time to cm/s^2
                P2[i].pres_z = (P2[i].pres_z - az1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].visc_x = (P2[i].visc_x - ax1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].visc_y = (P2[i].visc_y - ay1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].visc_z = (P2[i].visc_z - az1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].grav_x = (P2[i].grav_x - ax1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].grav_y = (P2[i].grav_y - ay1)*pow(Time, 0.5)*1.e5/dtconv;
                P2[i].grav_z = (P2[i].grav_z - az1)*pow(Time, 0.5)*1.e5/dtconv;

                dis = pow(((P2[i].Pos[0]-sinkposx1)*(P2[i].Pos[0]-sinkposx1) + (P2[i].Pos[1]-sinkposy1)*(P2[i].Pos[1]-sinkposy1) + (P2[i].Pos[2]-sinkposz1)*(P2[i].Pos[2]-sinkposz1)), 0.5); 
                dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                disAU = dis*206264.8060;
                dis=dis*3.086e18;   //dis in now in cm

                disxz = pow(((P2[i].Pos[0]-sinkposx1)*(P2[i].Pos[0]-sinkposx1) + (P2[i].Pos[2]-sinkposz1)*(P2[i].Pos[2]-sinkposz1)), 0.5);
                disxz=disxz*1.e3*Time/(0.7);                   //dis is in pc
                disxz=disxz*3.086e18;   //dis in now in cm

                disxy = pow(((P2[i].Pos[0]-sinkposx1)*(P2[i].Pos[0]-sinkposx1) + (P2[i].Pos[1]-sinkposy1)*(P2[i].Pos[1]-sinkposy1)), 0.5);
                disxy=disxy*1.e3*Time/(0.7);                   //dis is in pc
                disxy=disxy*3.086e18;   //dis in now in cm

                mass_part = P2[i].Mass*1.e10/0.7;

                x1 = P2[i].Pos[0];
                y1 = P2[i].Pos[1];
                z1 = P2[i].Pos[2];

                tpresx = (y1-sinkposy1)*P2[i].pres_z  - (z1-sinkposz1)*P2[i].pres_y;
                tpresx = tpresx*3.086e18*1.e3*Time/(0.7);   //convert to cm*(cm/s^2) 
                tpresy = (z1-sinkposz1)*P2[i].pres_x - (x1-sinkposx1)*P2[i].pres_z;
                tpresy = tpresy*3.086e18*1.e3*Time/(0.7);
                tpresz = (x1-sinkposx1)*P2[i].pres_y  -  (y1-sinkposy1)*P2[i].pres_x;
                tpresz = tpresz*3.086e18*1.e3*Time/(0.7);
                tpres = pow((tpresx*tpresx + tpresy*tpresy + tpresz*tpresz),0.5);
   
                tviscx = (y1-sinkposy1)*P2[i].visc_z  - (z1-sinkposz1)*P2[i].visc_y;
                tviscx = tviscx*3.086e18*1.e3*Time/(0.7);
                tviscy = (z1-sinkposz1)*P2[i].visc_x - (x1-sinkposx1)*P2[i].visc_z;
                tviscy = tviscy*3.086e18*1.e3*Time/(0.7);
                tviscz = (x1-sinkposx1)*P2[i].visc_y  -  (y1-sinkposy1)*P2[i].visc_x;
                tviscz = tviscz*3.086e18*1.e3*Time/(0.7);
                tvisc = pow((tviscx*tviscx + tviscy*tviscy + tviscz*tviscz),0.5);

                tgravx = (y1-sinkposy1)*P2[i].grav_z  - (z1-sinkposz1)*P2[i].grav_y;
                tgravx = tgravx*3.086e18*1.e3*Time/(0.7);
                tgravy = (z1-sinkposz1)*P2[i].grav_x - (x1-sinkposx1)*P2[i].grav_z;
                tgravy = tgravy*3.086e18*1.e3*Time/(0.7);
                tgravz = (x1-sinkposx1)*P2[i].grav_y  -  (y1-sinkposy1)*P2[i].grav_x;
                tgravz = tgravz*3.086e18*1.e3*Time/(0.7);
                tgrav = pow((tgravx*tgravx + tgravy*tgravy + tgravz*tgravz),0.5);
 
                if(P2[i].sink > 0.5)
                {
                tpresx = 0;
                tpresy = 0;
                tpresz = 0;
                tviscx = 0;
                tviscy = 0;
                tviscz = 0;
                }
 
               mass_enc = m_enc_arr[0];
               for(n=0; n<200; n++)
                 {
                 if(disAU > d_arr[n])
                   {
                   mass_enc_cent = m_enc_arr[n] + mass_sink1;
                   mass_enc = m_enc_arr[n];
                   }
                 }

                if(tvisc != tvisc)
                  tvisc = 0;

                if(tpres != tpres)
                  tpres = 0;

                omega=amom/(pow(dis,2));
                cs_squared_S1 = 850.0*1.38e-16/1.67e-24;
                h = pow(cs_squared_S1,0.5)/omega;
                pres = 0.1*P2[i].nh*1.38e-16*850.0;
                talpha = pres*2.0*3.14159*(pow(dis,2))*2.0*h;
                ttot = tgrav+tpres+tvisc;

               jzcent = pow(6.67e-8 * mass_enc_cent * 1.98892e33*dis, 0.5);
               kin = pow(P2[i].Vel[0]*P2[i].Vel[0] + P2[i].Vel[1]*P2[i].Vel[1] + P2[i].Vel[2]*P2[i].Vel[2], 0.5);
               pot =  pow(6.67e-8 * mass_enc * 1.98892e33/dis, 0.5);

               fprintf(outfile2, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", P2[i].Temp, tgrav, tpres, tvisc, P2[i].nh, disAU, amom, amom/ttot, amom/tgrav, amom*mass_enc* 1.98892e33/talpha);
                       
           }

           fclose(outfile2);

           infile3 = fopen("/work/utexas/ao/minerva/lam_120.dat", "r");
           outfile3= fopen("/work/utexas/ao/minerva/tcool_am2_120", "a");
           for(i = 0; i < N_gas3; i++)
              {
              if(i==0)
               allocate_memory3(N_gas3);

              fscanf(infile, " %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n",
                        &P3[i].dummy, &P3[i].nHcgs, &P3[i].Temp, &P3[i].lamHl, &P3[i].lamH2, &P3[i].dummy,
                        &P3[i].Pos[0], &P3[i].Pos[1], &P3[i].Pos[2], &P3[i].dummy, &P3[i].Id);
               }
            fclose(infile3);

            for(n = 0; n < N_gas3; n++)
               {
                if(n>0 && P3[n].Id == P3[n-1].Id)
                {
                P3[n-1].lamHl = 0;  //each particle appears in file twice, so only take every alternate entry
                P3[n-1].lamH2 = 0;
                P3[n-1].nHcgs = 0 ;
                }
                //f(P3[n].Temp < 500.0) P3[n].Temp = 500.0;
                if (n%100 ==0) printf("%lg %lg\n", P3[n].nHcgs, P3[n].lamHl);
               }

            for(i = 0; i < N_gas3; i++) 
               {
               dis = pow(((P3[i].Pos[0]-sinkposx1)*(P3[i].Pos[0]-sinkposx1) + (P3[i].Pos[1]-sinkposy1)*(P3[i].Pos[1]-sinkposy1) + (P3[i].Pos[2]-sinkposz1)*(P3[i].Pos[2]-sinkposz1)), 0.5);
               dis=dis*1.e3*Time/(0.7);                   //dis is in pc
               disAU = dis*206264.8060;
               dis=dis*3.086e18;   //dis in now in cm

               tcool = 0;
           
               if(P3[i].nHcgs > 0)
                 {
                 tcool = P3[i].nHcgs*kB*P3[i].Temp/(P3[i].lamHl + P3[i].lamH2);
                 fprintf(outfile3, "%15.11g %15.11g %15.11g %15.11g\n", P3[i].nHcgs, disAU, P3[i].Temp, tcool);
                 }
               } 
	    fclose(outfile3);

  free(P2);
}

  do_what_you_want();
}





/* here the particle data is at your disposal 
 */
int do_what_you_want(void)
{

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
	allocate_memory(NumPart);

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
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

/*	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(double), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	      {
		P[pc_sph].Ne= 1.0;
		pc_sph++;
	      }
*/

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      //fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
/*	       printf("h_smooth= %15.6g \n",P[pc_sph].hsm); */
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      //fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
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

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;


	}

      fclose(fd);
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
int allocate_memory(int N_part)
{
  printf("allocating memory...\n");

  if(!(P=(struct particle_data *) malloc((N_part+1)*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  printf("allocating memory...done\n");
}


int allocate_memory2(int N_part)
{
  printf("allocating memory...\n");

  if(!(P2=(struct particle_data2 *) malloc(N_part*sizeof(struct particle_data2))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
   
  printf("allocating memory...done\n");
}

int allocate_memory3(int N_part)
{
  printf("allocating memory...\n");

  if(!(P3=(struct particle_data3 *) malloc(N_part*sizeof(struct particle_data3))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

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






  











