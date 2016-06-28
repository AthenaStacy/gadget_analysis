#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define width 0.01

int load_snapshot(char *fname, int files);
int reordering(void);
int unit_conversion(void);
int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1);
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
  double  Pos_rot[3];
  double Pos_new[3];

  double  Vel[3];
  double Vel_rot[3];  
  double Vel_new[3];

  double  Mass;
  int    Type;
  int Id;

  double  Rho, Temp, nh;
  double  U, HI, HII; 
  //double HeII,  HeIII, H2II, DII, DM;
  double H2I, hsm, HDI, sink, gam;
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
  char path[200], input_fname[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  int  i, k=0, j, n, type, snapshot_number, files, Ngas, random, ncount, ncounthalo1, arrnum=200;
  double x,y,z,x1,y1,z1, delr;
  double nh, nhmax, mass, mmax, dis, dis_sim, disAU, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  int IDmax;
  double rho_b, rho_DM[arrnum], rho_vir, r_vir, v_vir;
  double kB, G;
  double bind, diskmass, massin=0.0, massout=0.0;
  double vrad, vrot, disAUxy, disAUz, disAUxyz, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2; 
  double sinkposx3, sinkposy3, sinkposz3, sinkposx4, sinkposy4, sinkposz4;
  double num, amom, amom_arr[arrnum], amomtot, amomx_tot, amomy_tot; 
  double amomz_tot, amomx, amomy, amomz, amomx1, amomx2, amomy1, amomy2, amomz1, amomz2;
  double xfac, yfac, zfac, kfloat, nh_avg, num_tot, omega, omega0, nthresh;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, jzcent;
  double temp_arr[arrnum], mbe_arr[arrnum], omega_arr[arrnum];
  double m_enc_arr[arrnum], d_arr[arrnum], mass_enc, n_arr[arrnum], n_arr2[arrnum], surf_dens_arr[arrnum];
  double rmin, rmax, nmin, nmax, numdens_arr[arrnum]; 
  double vradx, vrady, vradz, vrotx, vroty, vrotz, vrad_arr[arrnum], vrot_arr[arrnum], dens_arr[arrnum], c_s_arr[arrnum]; 
  double vradx_arr[arrnum], vrady_arr[arrnum], vradz_arr[arrnum], vrotx_arr[arrnum], vroty_arr[arrnum], vrotz_arr[arrnum];
  double velx, vely, velz, velx_arr[arrnum], vely_arr[arrnum], velz_arr[arrnum];
  double mdm_enc_arr[arrnum], turb_arr[arrnum], disx, disy, disz;
  double KE, KE_tot, PE, PE_tot, vel, spin_param, ncount_doub;
  int ID_sink1, ID_sink2, ID_sink3, ID_sink4, sinknumk;
  FILE *outfile, *outfile2;

  sprintf(path, "/work/00863/minerva/");
  //sprintf(basename, "ifront_hiacc");

  //sprintf(path, "/nobackupp1/astacy/binHR");
  //sprintf(basename, "bin_HR1");
  //sprintf(basename, "bin_HR7");
  //sprintf(basename, "bin_HR10");

  //sprintf(path, "/nobackupp1/astacy/bin_zoom");
  //sprintf(basename, "bin_zoom1");
  //sprintf(basename, "bin_zoom2");
  //sprintf(basename, "bin_zoom3");
  //sprintf(basename, "bin_zoom4");
  //sprintf(basename, "bin_zoom5");
  //sprintf(basename, "bin_zoom6");
  sprintf(basename, "bin_zoom7");
  //sprintf(basename, "bin_zoom10");


  //for(j=453;j<=453;j++)
  //for(j=6;j<=6;j++)	
  //for(j=9;j<=9;j++)
  //for(j=1;j<=1;j++)
  //for(j=90;j<=90;j++)
  //for(j=77;j<=77;j++)
  //for(j=86;j<=86;j++)
  //for(j=127;j<=127;j++)
  //for(j=112;j<=112;j++)
  //for(j=139;j<=139;j++)
  //for(j=91;j<=91;j++)
  for(j=101;j<=101;j++)

  //for(j=80;j<=80;j++)
  {
  
  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  //sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  Ngas = load_snapshot(input_fname, files);

  sprintf(output_fname, "%s_nhprof_%03d.dat",basename, snapshot_number);

          /*    reordering();*/ /* call this routine only if your ID's are set properly */

  rho_b = 2.94e-30; //background density of the universe

  unit_conversion();

  nh_avg=ncount_doub=0;
  num_tot=0;
  h2max=nhmax=mmax=kfloat=0.0;
  xtimesm = ytimesm = ztimesm = 0.0;
  vxtimesm = vytimesm = vztimesm = 0.0;
  vxCOM=vyCOM=vzCOM=xCOM=yCOM=zCOM=0;
  amomx_tot = amomy_tot = amomz_tot = 0;
  amomtot = omega0 = masstot = 0.0;
  
  nthresh = 1.e4;	  	
	
       for(i = 0; i < Ngas; i++)
        {
       if(P[i].nh > nhmax /*&& P[i].sink < 0.5*/)
   	   {
           nhmax = P[i].nh;
           IDmax = P[i].Id;
           xmax = P[i].Pos[0];
           ymax = P[i].Pos[1];
           zmax = P[i].Pos[2];
           mmax = P[i].Mass;
           }
         }
	printf("xmax = %lg ymax = %lg zmax = %lg mmax = %lg nhmax = %lg\n", xmax, ymax, zmax, mmax, nhmax);

        for(n=0;n<=Ngas;n++) 
         {
           nh = P[n].nh;
           dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
           dis=dis*1.e3*Time/(0.7);
           disx = fabs((P[n].Pos[0]-xmax))*1.e3*Time/(0.7);
           disy = fabs((P[n].Pos[1]-ymax))*1.e3*Time/(0.7);
           disz = fabs((P[n].Pos[2]-zmax))*1.e3*Time/(0.7);
           disAU=dis*206264.806;
           if(disx < width && disy < width && disz < width)
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

         vx1 = vxCOM/ncount_doub;
         vy1 = vyCOM/ncount_doub;
         vz1 = vzCOM/ncount_doub;
         //vx = vy = vz = 0;

	sinkposx1 = xmax;
	sinkposy1 = ymax;
	sinkposz1 = zmax;

        printf("vx1 = %lg, vy1 = %lg, vz1 = %lg, sinkposx1 = %lg, sinkposy1 = %lg, sinkposz1 = %lg\n", vx1, vy1, vz1, sinkposx1, sinkposy1, sinkposz1);

	diskmass=0.0;
	ncount=0;  
	ncount_doub=0.;
	  

	rmin = 100.;
        rmax = 3e8;
        nmin = 1.e-2;
        nmax = 1.e12;

        for(n=0; n<arrnum; n++)
          {
          m_enc_arr[n] = amom_arr[n] = n_arr[n] = n_arr2[n]  = surf_dens_arr[n] = c_s_arr[n] = 0.0;
          temp_arr[n] = mbe_arr[n] = omega_arr[n] = 0.0; 
          vrot_arr[n] = vrotx_arr[n] = vroty_arr[n] = vrotz_arr[n] = 0;
          vrad_arr[n] = vradx_arr[n] = vrady_arr[n] = vradz_arr[n] = 0;
          rho_DM[n] = 0.0;
          //d_arr[n] = (double(n)/199.0)*5.e8;
          //d_arr[n] = pow(5.e8,double(n)/double(arrnum));
          d_arr[n] = rmin * pow(rmax/rmin, double(n)/double(arrnum));
          numdens_arr[n] = nmin * pow(nmax/nmin, double(n)/double(arrnum));
          //printf("d_arr[%d] = %lg\n", n, d_arr[n]);
          }
 
         //for(i = Ngas; i < NumPart; i++)
         for(i = 0; i < NumPart; i++)
               {
	       dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
	       dis=dis*1.e3*Time/(0.7);                   //dis is in pc
	       disAU = dis*206264.8060;
	       if(disAU < 5.e8)
                    {
                    for(n=0; n<arrnum; n++)
                      {
                      if(disAU < d_arr[n] && P[i].Type == 0)
                        m_enc_arr[n] = m_enc_arr[n] + P[i].Mass*1.e10/0.7;
                      if(disAU < d_arr[n] && P[i].Type == 1)
                        mdm_enc_arr[n] = mdm_enc_arr[n] + P[i].Mass*1.e10/0.7;
                      } 
                     }
                  }

          for(n=0; n<arrnum; n++)
               {
               rho_DM[n] = mdm_enc_arr[n]*1.98892e33/((4.0/3.0)*3.14159*pow(d_arr[n]*1.5e13, 3.0));  //DM density in cgs
               printf("mdm_enc_arr[%d] = %lg\n", n, mdm_enc_arr[n]);
               //printf("rho_DM[%d] = %lg\n", n, rho_DM[n]/(rho_b*pow(Time, -3.0)));
               if(rho_DM[n] > 200.0*rho_b*pow(Time, -3.0))
                 {
                 rho_vir = rho_DM[n];
                 r_vir = d_arr[n];
                 }
               }
	  
	
             //rotate(sinkposx1,sinkposy1,sinkposz1, vx1, vy1, vz1);
 
         //for(i = Ngas; i < NumPart; i++)
         for(i = 0; i < Ngas; i++)
                {

                 P[i].Vel[0] = P[i].Vel[0] - vx1;
                 P[i].Vel[1] = P[i].Vel[1] - vy1;
                 P[i].Vel[2] = P[i].Vel[2] - vz1;
          
                 amomx = (P[i].Pos[1]-sinkposy1)*P[i].Vel[2]  - (P[i].Pos[2]-sinkposz1)*P[i].Vel[1];
                 amomy = (P[i].Pos[2]-sinkposz1)*P[i].Vel[0] - (P[i].Pos[0]-sinkposx1)*P[i].Vel[2];
                 amomz = (P[i].Pos[0]-sinkposx1)*P[i].Vel[1]  -  (P[i].Pos[1]-sinkposy1)*P[i].Vel[0];

                 amomx = (P[i].Mass*1.e10*1.98892e33/.7)*amomx*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);  //convert to cgs units
                 amomy = (P[i].Mass*1.e10*1.98892e33/.7)*amomy*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);
                 amomz = (P[i].Mass*1.e10*1.98892e33/.7)*amomz*pow(Time, 0.5)*1.e5*3.086e18*1.e3*Time/(0.7);

                 amom = pow(amomx*amomx + amomy*amomy + amomz*amomz, 0.5);  //ang. mom. of particle in cgs

                 vel = pow(P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2], 0.5);
                 vel = vel*pow(Time, 0.5)*1.e5;  //convert vel to cgs 

                 KE = 0.5*(P[i].Mass*1.e10*1.98892e33/.7)*pow(vel,2.0);   //kinetic energy in cgs

                 velx = P[i].Vel[0]*pow(Time, 0.5);
                 vely = P[i].Vel[1]*pow(Time, 0.5);
                 velz = P[i].Vel[2]*pow(Time, 0.5);

                 dis_sim = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5);
                 dis = pow(((P[i].Pos[0]-sinkposx1)*(P[i].Pos[0]-sinkposx1) + (P[i].Pos[1]-sinkposy1)*(P[i].Pos[1]-sinkposy1) + (P[i].Pos[2]-sinkposz1)*(P[i].Pos[2]-sinkposz1)), 0.5); 
                 dis=dis*1.e3*Time/(0.7);                   //dis is in pc
                 disAU = dis*206264.8060;
                 dis=dis*3.086e18;   //dis in now in cm                

                 vrot  = amom/(P[i].Mass*1.e10*1.98892e33/.7)/dis/1e5;
                 vrotx = amomx/(P[i].Mass*1.e10*1.98892e33/.7)/dis/1e5;              
                 vroty = amomy/(P[i].Mass*1.e10*1.98892e33/.7)/dis/1e5;
                 vrotz = amomz/(P[i].Mass*1.e10*1.98892e33/.7)/dis/1e5; //vrot in km/s 
                 omega = (P[i].Mass/1.e-10/.7)*vrot/(disAU*1.5e8); 

             if(P[i].nh > nthresh && disAU > 0)
                   {
                   masstot = masstot + P[i].Mass/1.e-10/.7;
                   omega0 = omega0 + (P[i].Mass/1.e-10/.7)*vrot/(disAU*1.5e8);  //divide vrot by distance in km
                   }

                 disAUxyz = pow((P[i].Pos_new[0])*(P[i].Pos_new[0]) + (P[i].Pos_new[1])*(P[i].Pos_new[1]) + (P[i].Pos_new[2])*(P[i].Pos_new[2]), 0.5);
                 disAUxyz = disAUxyz*1.e3*Time/(0.7)*206264.8060;

                 disAUxy = pow(((P[i].Pos_new[0])*(P[i].Pos_new[0]) + (P[i].Pos_new[1])*(P[i].Pos_new[1])), 0.5);
                 disAUxy = disAUxy*1.e3*Time/(0.7)*206264.8060;
 
                 disAUz = fabs(P[i].Pos_new[2]*1.e3*Time/(0.7)*206264.8060);  

                 vrad =  (P[i].Vel[0]*(P[i].Pos[0]-sinkposx1) + P[i].Vel[1]*(P[i].Pos[1]-sinkposy1) + P[i].Vel[2]*(P[i].Pos[2]-sinkposz1))/dis_sim;
                 vrad = vrad*pow(Time, 0.5);
                 vradx = (P[i].Vel[0]*(P[i].Pos[0]-sinkposx1))/dis_sim*pow(Time, 0.5);
                 vrady = (P[i].Vel[1]*(P[i].Pos[1]-sinkposy1))/dis_sim*pow(Time, 0.5);
                 vradz = (P[i].Vel[2]*(P[i].Pos[2]-sinkposz1))/dis_sim*pow(Time, 0.5); 

                 if(disAU < r_vir /*&& KE + PE < 0.0*/)
                   {
                   KE_tot = KE_tot + KE;
                   PE_tot = PE_tot + PE;
                   amomtot = amomtot + amom;
                   amomx_tot = amomx_tot + amomx;          
                   amomy_tot = amomy_tot + amomy;
                   amomz_tot = amomz_tot + amomz;
                   num_tot = num_tot + P[i].Mass;
		   ncount++;  
                   }

                 if(P[i].nh > 1.e-2)
                   for(n=0; n<arrnum; n++)
                    {
                    if(P[i].nh > numdens_arr[n] && P[i].nh < numdens_arr[n+1] && P[i].sink < 0.5)
                      {
                      amom_arr[n] = amom_arr[n] + amom/(P[i].Mass*1.e10*1.98892e33/.7);
                      dens_arr[n] = dens_arr[n] + P[i].nh;
                      c_s_arr[n] = c_s_arr[n] + pow(1.38e-16*P[i].Temp/1.e-24,0.5)/1.e5;
                      temp_arr[n] = temp_arr[n] + P[i].Temp;
                      mbe_arr[n] = mbe_arr[n] + 1050.*pow(P[i].Temp/200.,1.5)*pow(P[i].nh/1.e4,-0.5);

                      velx_arr[k] = velx_arr[k] + velx;
                      vely_arr[k] = vely_arr[k] + vely;
                      velz_arr[k] = velz_arr[k] + velz;
                      vrad_arr[n] = vrad_arr[n] + vrad;
                      //vrot_arr[n] = vrot_arr[n] + vrot;
                      vradx_arr[n] = vradx_arr[n] + vradx;
                      vrady_arr[n] = vrady_arr[n] + vrady;
                      vradz_arr[n] = vradz_arr[n] + vradz;
                      vrotx_arr[n] = vrotx_arr[n] + vrotx;
                      vroty_arr[n] = vroty_arr[n] + vroty;
                      vrotz_arr[n] = vrotz_arr[n] + vrotz;
                      n_arr[n]=n_arr[n] + 1.0;
                      }
                     if(disAU < d_arr[n+1] && disAU > 0)
                       omega_arr[n] = omega_arr[n] + omega;
                     }

                 if(P[i].nh > 1.e-2)
                   for(n=0; n<arrnum; n++) 
                    {   
                    if(disAUxy > d_arr[n] && disAUxy < d_arr[n+1] && disAUz < 10.  && P[i].sink < 0.5)
                    //if(disAUxyz > d_arr[n] && disAUxyz < d_arr[n+1] && P[i].sink < 0.5)
                     {
                     surf_dens_arr[n] = surf_dens_arr[n] + P[i].Mass*1.e10*1.98892e33/.7;
                     n_arr2[n]=n_arr2[n] + 1.0;                     
                     }
                    }
                }

             outfile=fopen(output_fname, "a");
             for(n=0; n<arrnum; n++)
                   {
                   printf("n_arr2 = %lg , c_s = %lg, vrad = %lg, surf_dens = %lg omega = %lg\n", n_arr2[n], c_s_arr[n], vrad_arr[n], surf_dens_arr[n], omega_arr[n]);
                   dens_arr[n] = dens_arr[n]/n_arr[n];

                   velx_arr[n] = velx_arr[n]/n_arr[n];
                   vely_arr[n] = vely_arr[n]/n_arr[n];
                   velz_arr[n] = velz_arr[n]/n_arr[n];
                   vradx_arr[n] = vradx_arr[n]/n_arr[n];
                   vrady_arr[n] = vrady_arr[n]/n_arr[n];
                   vradz_arr[n] = vradz_arr[n]/n_arr[n];
                   vrotx_arr[n] = vrotx_arr[n]/n_arr[n];
                   vroty_arr[n] = vroty_arr[n]/n_arr[n];
                   vrotz_arr[n] = vrotz_arr[n]/n_arr[n];
                   //vrad_arr[n] = pow(vradx_arr[n]*vradx_arr[n] + vrady_arr[n]*vrady_arr[n] + vradz_arr[n]*vradz_arr[n], 0.5)/n_arr[n];
                   vrad_arr[n] = vrad_arr[n]/n_arr[n];
                   vrot_arr[n] = pow(vrotx_arr[n]*vrotx_arr[n] + vroty_arr[n]*vroty_arr[n] + vrotz_arr[n]*vrotz_arr[n], 0.5)/n_arr[n];
                   c_s_arr[n] = c_s_arr[n]/n_arr[n];
                   temp_arr[n] = temp_arr[n]/n_arr[n];
                   mbe_arr[n] = mbe_arr[n]/n_arr[n];
                   surf_dens_arr[n] = surf_dens_arr[n]/3.14159/(pow(1.4958e13*d_arr[n+1],2) - pow(1.4958e13*d_arr[n],2));

	           turb_arr[n] = pow(velx_arr[n] - vradx_arr[n] - vrotx_arr[n],2) + pow(vely_arr[n] - vrady_arr[n] - vroty_arr[n],2) + pow(velz_arr[n] - vradz_arr[n] - vrotz_arr[n],2);
      		   turb_arr[n] = pow(turb_arr[n], 0.5);

      		   vrot_arr[n] = pow(vrotx_arr[n]*vrotx_arr[n] + vroty_arr[n]*vroty_arr[n] + vrotz_arr[n]*vrotz_arr[n], 0.5);
      		   if(n > 0)
        		amom_arr[n] = (m_enc_arr[k] - m_enc_arr[k-1])*vrot_arr[n]*d_arr[n]*1.5e13*1.e5*2.e33;
      		   if(n == 0)
        		amom_arr[n] = m_enc_arr[0]*vrot_arr[0]*d_arr[0]*1.5e13*1.e5*2.e33;


                   fprintf(outfile, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", numdens_arr[n], m_enc_arr[n], mdm_enc_arr[n], dens_arr[n], surf_dens_arr[n], c_s_arr[n], temp_arr[n], mbe_arr[n], omega_arr[n]/m_enc_arr[n], vrad_arr[n], vrot_arr[n], turb_arr[n], amom_arr[n]);
                   }
             fclose(outfile);
 
             free(P);
  }

}





/* here the particle data is at your disposal 
 */
int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1)
{
double dis, disAU, vec[3], r, alpha, beta, d1[3][3], d2[3][3], d12[3][3];
double dmax = 200.0;
int i, j,k, n;

vec[0] = vec[1] = vec[2] = 0; 

for(i=0; i<NumPart; i++)
   {
   P[i].Pos_new[0] = 0;
   P[i].Pos_new[1] = 0;
   P[i].Pos_new[2] = 0;
   P[i].Vel_new[0] = 0;
   P[i].Vel_new[1] = 0;
   P[i].Vel_new[2] = 0;
   P[i].Pos_rot[0] = P[i].Pos[0] - x1;
   P[i].Pos_rot[1] = P[i].Pos[1] - y1;
   P[i].Pos_rot[2] = P[i].Pos[2] - z1;
   P[i].Vel_rot[0] = P[i].Vel[0] - vx1;
   P[i].Vel_rot[1] = P[i].Vel[1] - vy1;
   P[i].Vel_rot[2] = P[i].Vel[2] - vz1;
   }

for(i=0; i<NumPart; i++)
  {
  dis = pow(((P[i].Pos[0]-x1)*(P[i].Pos[0]-x1) + (P[i].Pos[1]-y1)*(P[i].Pos[1]-y1) + (P[i].Pos[2]-z1)*(P[i].Pos[2]-z1)), 0.5);
  dis=dis*1.e3*Time/(0.7);                   //dis is in pc
  disAU = dis*206264.8060;

  if(disAU < dmax)
    {
    vec[0] = vec[0] + P[i].Mass * (P[i].Pos_rot[1] *  P[i].Vel_rot[2] - P[i].Pos_rot[2] *  P[i].Vel_rot[1]);
    vec[1] = vec[1] + P[i].Mass * (P[i].Pos_rot[2] *  P[i].Vel_rot[0] - P[i].Pos_rot[0] *  P[i].Vel_rot[2]);
    vec[2] = vec[2] + P[i].Mass * (P[i].Pos_rot[0] *  P[i].Vel_rot[1] - P[i].Pos_rot[1] *  P[i].Vel_rot[0]);
    }
  }

r = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
r = pow(r, 0.5);

for(i=0; i<=2; i++)
  {
  vec[i] = vec[i]/r;
  vec[i] <= 1;
  vec[i] >= -1;
  }

alpha = acos(vec[2] / pow(pow(vec[1],2) + pow(vec[2],2), 0.5));
beta = asin(vec[0]);

printf("alpha = %lg, beta = %lg \n", alpha, beta);

for(i = 0; i<=2; i++)
  for(j = 0; j<=2; j++)
    { 
    d1[i][j] = 0;
    d2[i][j] = 0;
    d12[i][j] = 0;
    }


d1[0][0] = 1;
d1[1][1] = cos(alpha);
d1[1][2] = sin(alpha);
d1[2][1] = -sin(alpha);
d1[2][2] = cos(alpha);

d2[0][0] = cos(beta);
d2[0][2] = -sin(beta);
d2[1][1] = 1;
d2[2][0] = sin(beta);
d2[2][2] = cos(beta);

for(i = 0; i<=2; i++)  
  for(j = 0; j<=2; j++)  
    for(k = 0; k<= 2; k++)  
      d12[i][j] += d2[i][k] * d1[k][j];

printf("d1 = %lg, d1 = %lg, d1 = %lg, d2 = %lg, d2 = %lg, d2 = %lg\n", d1[1][1] , d1[1][2], d1[2][1], d2[0][0], d2[0][2], d2[2][0]);
printf("d2 = %lg, d2 = %lg d2 = %lg d2 = %lg\n", d12[0][2],  d12[0][0],  d12[1][1], d12[2][2]);

for(n=0; n<NumPart; n++)
  {

  
  dis = pow(((P[n].Pos[0]-x1)*(P[n].Pos[0]-x1) + (P[n].Pos[1]-y1)*(P[n].Pos[1]-y1) + (P[n].Pos[2]-z1)*(P[n].Pos[2]-z1)), 0.5);
  dis=dis*1.e3*Time/(0.7);                   //dis is in pc
  disAU = dis*206264.8060;


  if(P[n].nh > 1.e1)
    {
    for(i = 0; i<=2; i++)  
      for(j = 0; j<=2; j++)  
        P[n].Pos_new[i] += d12[i][j] * P[n].Pos_rot[j];

    for(i = 0; i<=2; i++) 
      for(j = 0; j<= 2; j++) 
        P[n].Vel_new[i] += d12[i][j] * P[n].Vel[j];

    if(n%100000 == 0)
      {
      printf("P[n].Pos_rot[0] = %lg, P[n].Pos_rot[1] = %lg, P[n].Pos_rot[2] = %lg\n", P[n].Pos_rot[0], P[n].Pos_rot[1], P[n].Pos_rot[2]);
      printf("P[n].Pos_new[0] = %lg, P[n].Pos_new[1] = %lg, P[n].Pos_new[2] = %lg\n", P[n].Pos_new[0], P[n].Pos_new[1], P[n].Pos_new[2]); 
      }
    }
  }
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
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
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






  











