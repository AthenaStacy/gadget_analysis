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
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  i, k=0, j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, ncount2, halo;
  double x,y,z,x1,y1,z1, delr;
  double nh, nhmax, mass, mmax, dis, disAU, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double kB, G;
  double massHIT, diskmass, massin=0.0, massout=0.0;
  double mass_sink1, mass_sink2, mass_sink3, mass_sink4, mass_tot;
  double dis12, dis13, x_com, y_com, z_com, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2, sinkposx3, sinkposy3, sinkposz3, sinkposx4, sinkposy4, sinkposz4;
  double num, numz, tempz, vz, rhoz, amom, amomtot, amomtotx, amomtoty, amomtotz, amomx, amomy, amomz, amomx1, amomx2, amomy1, amomy2, amomz1, amomz2;
  double xfac, vrad, vrot, vrotx, vroty, vrotz, omega, nh_avg, temp_avg, vrad_avg, vrot_avg, omega_avg, num_tot, num_tot2;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, jzcent, jnum, jnum2;
  double m_enc_arr[200], d_arr[200], mass_enc, xshield, b5, fshield, fshield_tot, fac;
  int ID_sink1, ID_sink2, ID_sink3, ID_sink4, sinknum;
  FILE *outfile;
  FILE *outfile2;


  //sprintf(path, "/nobackupp1/astacy/hires");
  sprintf(path, "/nobackupp1/astacy/");
  //sprintf(path, "/nobackupp1/astacy/binHR");

  halo = 0;

    if(halo == 0)
     {
     //sprintf(basename, "hires_test10");
     //sprintf(basename, "hires_test10b");
     sprintf(basename, "midhires_testLW");
     //sprintf(basename, "hires3");
     //jnum = 228;
     jnum = 200;
     jnum2 = 209;
     }
  if(halo == 1)
     {
     sprintf(basename, "bin_HR1");
     //jnum = 228;
     jnum = 6;
     }
  if(halo == 2)
     {
     sprintf(basename, "bin_HR2");
     //jnum = 207;
     jnum = 5;
     }
  if(halo == 3)
     {
     sprintf(basename, "bin_HR3");
     //jnum = 206;
     jnum = 7;
     }
  if(halo == 4)
     {
     sprintf(basename, "bin_HR4");
     //jnum = 204;
     jnum = 7;
     }
  if(halo == 5)
     {
     sprintf(basename, "bin_HR5");
     //jnum = 209;
     jnum = 7;
     }
  if(halo == 6)
     {
     sprintf(basename, "bin_HR6");
     //jnum = 224;
     jnum = 7;
     }
  if(halo == 7)
     {
     sprintf(basename, "bin_HR7b");
     //jnum = 655;
     jnum = 9;
     //jnum=38;
     }
  if(halo == 8)
     {
     sprintf(basename, "bin_HR8");
     //jnum = 219;
     jnum = 7;
     }
  if(halo == 9)
     {
     sprintf(basename, "bin_HR9");
     //jnum = 207;
     jnum = 5;
     }
  if(halo == 10)
     {
     //sprintf(basename, "bin_HR10");
     sprintf(basename, "bin_HR10_alt");
     //jnum = 241;
     //jnum = 7;
     jnum=8;
     //jnum = 24;
     }


  //for(j=2;j<=50;j=j+1){
  for(j=jnum;j<=jnum2;j=j+1){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  kB = 1.38e-16;
  G = 6.67e-8;
  diskmass = 0.0;

  //sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s_fshield.dat",basename);
  Ngas = load_snapshot(input_fname, files);

  ID_sink1 = 4352905;
  ID_sink2 = 4333881;
  sinknum = 0;
          /*    reordering();*/ /* call this routine only if your ID's are set properly */

printf("hello line 187\n");

          unit_conversion();

printf("hello line 191\n");

          nh_avg= vrad_avg = vrot_avg = omega_avg = 0;

printf("hello line 195\n");

          temp_avg=0;
          num_tot=num_tot2=0;
          h2max=nhmax=0;
          masstot=0;
          mmax=0;
          xtimesm = 0.0;
          ytimesm = 0.0;
          ztimesm = 0.0;
          vxtimesm = 0.0;
          vytimesm = 0.0;
          vztimesm = 0.0;


             for(i = 0; i < Ngas; i++)
                {
                if(Id[i] == ID_sink1)
                  {
                  sinkposx1 = P[i].Pos[0];
                  sinkposy1 = P[i].Pos[1];
                  sinkposz1 = P[i].Pos[2];
                  vx1 = P[i].Vel[0];
                  vy1 = P[i].Vel[1];
                  vz1 = P[i].Vel[2];
                  mass_sink1 = P[i].Mass*1.e10/0.7;
                  printf("Found the sink! mass= %lg, nh = %lg\n", mass_sink1, P[i].nh);
                  }

                 if(P[i].nh > 1e9)
                  {
                  diskmass = diskmass + (P[i].Mass*1.e10/0.7);
                  xtimesm = xtimesm + P[i].Pos[0]*(P[i].Mass*1.e10/0.7);
                  ytimesm = ytimesm + P[i].Pos[1]*(P[i].Mass*1.e10/0.7);
                  ztimesm = ztimesm + P[i].Pos[2]*(P[i].Mass*1.e10/0.7);
                  vxtimesm = vxtimesm + P[i].Vel[0]*(P[i].Mass*1.e10/0.7);
                  vytimesm = vytimesm + P[i].Vel[1]*(P[i].Mass*1.e10/0.7);
                  vztimesm = vztimesm + P[i].Vel[2]*(P[i].Mass*1.e10/0.7);
                  }
          
                 if(P[i].nh > nhmax)
                   nhmax = P[i].nh;
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

         diskmass=fshield_tot=0.0;
         ncount=ncount2=0;

        for(n=0; n<200; n++)
               {
               m_enc_arr[n] = 0.0;
               d_arr[n] = (double(n)/199.0)*10000.0;
               printf("d_arr[%d] = %lg\n", n, d_arr[n]);
               }

       for(i = 0; i < Ngas; i++)
               {
                if(P[i].nh > 1.e8)
                  {
                   dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                   dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                   disAU = dis*206264.8060;
                   for(n=0; n<200; n++)
                     {
                     if(disAU < d_arr[n] /*|| P[i].sink > 0.5*/)
                       m_enc_arr[n] = m_enc_arr[n] + P[i].Mass*1.e10/0.7;
                     }
                    }
                  }

        for(n=0; n<200; n++)
               {
               printf("m_enc_arr[%d] = %lg\n", n, m_enc_arr[n]);
               }


         for(i = 0; i < Ngas; i++)
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

                 amom = pow(amomx*amomx + amomy*amomy + amomz*amomz, 0.5);

                 dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                 //dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                 dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                 disAU = dis*206264.8060;
                 dis=dis*3.086e18;   //dis in now in cm                

                 vrad =  (P[n].Vel[0]*(P[n].Pos[0]-sinkposx1) + P[n].Vel[1]*(P[n].Pos[1]-sinkposy1) + P[n].Vel[2]*(P[n].Pos[2]-sinkposz1))/pow(((P[n].Pos[0]-sinkposx1)*(P[n].Pos[0]-sinkposx1) + (P[n].Pos[1]-sinkposy1)*(P[n].Pos[1]-sinkposy1) + (P[n].Pos[2]-sinkposz1)*(P[n].Pos[2]-sinkposz1)), 0.5);
                 vrad = vrad*pow(Time, 0.5);


                 vrotx = (P[n].Pos[1]-sinkposy1)*P[n].Vel[2]  - (P[n].Pos[2]-sinkposz1)*P[n].Vel[1];
                 vroty = (P[n].Pos[2]-sinkposz1)*P[n].Vel[0] - (P[n].Pos[0]-sinkposx1)*P[n].Vel[2];
                 vrotz = (P[n].Pos[0]-sinkposx1)*P[n].Vel[1]  -  (P[n].Pos[1]-sinkposy1)*P[n].Vel[0];

                 vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5);
                 vrot = vrot*pow(Time, 0.5)/pow(((P[n].Pos[0]-sinkposx1)*(P[n].Pos[0]-sinkposx1) + (P[n].Pos[1]-sinkposy1)*(P[n].Pos[1]-sinkposy1) + (P[n].Pos[2]-sinkposz1)*(P[n].Pos[2]-sinkposz1)), 0.5);

                 omega = vrot/(dis/1.e5);  // pi factors CANCEL OUT -- 2 pi r to get circumference CANCELS with the 2 pi radians per circumference

                 if(P[i].nh > 1.e9 && P[i].H2I > 1.e-3 || P[i].sink > 0.5)
                  {
                  diskmass = diskmass + (P[i].Mass*1.e10/0.7);
                  nh_avg = nh_avg +(1- 2.*P[i].H2I);

                  if(P[i].sink < 0.5)
                   temp_avg = temp_avg + P[i].Temp;
                  
                  vrad_avg = vrad_avg + fabs(vrad);
                  vrot_avg = vrot_avg + vrot;
                  omega_avg = omega_avg + omega;
                  num_tot = num_tot + 1.0;

                  if(P[i].sink > 0.5)
                   sinknum++;

                  ncount++;
                  }

                if(P[i].nh > nhmax/10.0)
                  {
                  xshield = 4.*(P[i].H2I/1e-6)*pow(P[i].nh/10.,0.5)*pow(P[i].Temp/1.e3,0.5);
                  b5 = 2.9*pow(P[i].Temp/1e3,0.5);
                  fac = pow(1. + xshield, 0.5);
                  fshield = 0.965*pow(1 + xshield/b5,-1.1) + (0.035/fac)*exp(-8.5e-4*fac);
                  fshield_tot = fshield_tot + fshield;
                  num_tot2 = num_tot2 + 1.0;
                  ncount2++;
                  }

               }
             printf("temp = %lg, num_tot = %lg, temp_avg = %lg\n", temp_avg, num_tot, temp_avg/num_tot);

             outfile=fopen(output_fname, "a");
             fprintf(outfile, "%15.11g %15.11g %15.11g \n", Time, nhmax, fshield_tot/num_tot2);
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


  for(i=1; i<=NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
        {
/*        MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */

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
                                                                  




  











