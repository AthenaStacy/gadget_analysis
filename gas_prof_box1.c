#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define Msun 1.98892e33

#define width 12.8
//#define width 6.4
//#define width 3.2
//#define width 1.6
//#define width 0.8
//#define width 0.2

#define wpot 1

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

  double  Rho, U, Pres, Temp, nh;
  double H2I, HII, HDI, HeII, gam, sink;
  double pot;
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
  char path[200], input_fname[200], output_fname[200], output_fname2[200], output_fname3[200], output_fname4[200], basename[200], basenameout[200];
  int  i, j, jnum, k, arrnum=128, halo, n, nsink, nsink_bin, type, snapshot_number, files, Ngas, random, ncount, pair_count, idmax;
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, dismax, ke, pe, vesc, omega0, omega_crit;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, disCOM, disCOM1, disCOM2, vCOM, vCOM1, ncount_doub;
  double nh, nhmax, nthresh, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, dis_pc, disAU, vrad, vradx, vrady, vradz, vrotx, vroty, vrotz, vrot, mh_mass;
  double m_enc, m_enc_arr[arrnum], m_shell_arr[arrnum], d_arr[arrnum], n_arr[arrnum], num_arr[arrnum], num2_arr[arrnum]; 
  double vol_arr[arrnum], c_s_arr[arrnum], vrad_arr[arrnum], vrot_arr[arrnum];
  double velx, vely, velz, disx, disy, disz, mass_conv;
  double vel_arr[arrnum], velx_arr[arrnum], vely_arr[arrnum], velz_arr[arrnum], vradx_arr[arrnum]; 
  double hd_arr[arrnum], nh_arr[arrnum], h2_arr[arrnum], elec_arr[arrnum], temp_arr[arrnum], temp_nh_arr[arrnum];
  double eint_arr[arrnum], gpot_arr[arrnum], xmom_arr[arrnum], ymom_arr[arrnum], zmom_arr[arrnum];
  double gpot_xL[arrnum], gpot_xR[arrnum], gpot_yL[arrnum], gpot_yR[arrnum], gpot_zL[arrnum], gpot_zR[arrnum];
  double gacc_xL[arrnum], gacc_xR[arrnum], gacc_yL[arrnum], gacc_yR[arrnum], gacc_zL[arrnum], gacc_zR[arrnum], DeltaX, dfac;
  double num_xL[arrnum], num_xR[arrnum], num_yL[arrnum], num_yR[arrnum], num_zL[arrnum], num_zR[arrnum];
  double rmax, rmin, nmin, nmax, volume, hubble;
  double au_cm = 1.49597871e13, pc_cm = 3.08567758e18;
  double pot_min=1.e30, pot_max = -1.e30;
  FILE *outfile, *outfile2, *outfile3, *outfile4;
  
  sprintf(path, "/work/00863/minerva");

  halo = 8;

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
     sprintf(basename, "iso_map");
     sprintf(basenameout, "iso_map");
     jnum = 1;
     }
  if(halo == 7)
     {
     sprintf(basename, "dma_prof1_ng");
     sprintf(basenameout, "dma_prof1_ng");
     jnum = 1;
     }
  if(halo == 8)
     {
     sprintf(basename, "/work/00863/minerva/");
     sprintf(basenameout, "bin_HR10_ref3_wpot");
     sprintf(basename, "bin_HR10_ref3_wpot");
     jnum = 7131;
     }
  if(halo == 9)
     {
     sprintf(path, "/work/00863/minerva/bin_map");
     sprintf(basename, "bin_HR10_map");
     sprintf(basenameout, "bin_HR10_map");
     jnum = 0;
     //jnum = 10;
     //jnum = 19;
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

  printf("hello line 139\n");

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s_gas_%04d.dat",basename, snapshot_number);
  sprintf(output_fname2, "%s_gas_nh_%04d.dat",basename, snapshot_number);
  sprintf(output_fname3, "%s_energy_%04d.dat",basename, snapshot_number);
  sprintf(output_fname4, "%s_gpot_%04d.dat",basename, snapshot_number);
  Ngas = load_snapshot(input_fname, files);

  hubble = 0.7;

  //for NON-comoving integrations
  //Time = hubble = 1.0;

  printf("hello line 145\n");

  /*    reordering();*/ /* call this routine only if your ID's are set properly */

  unit_conversion(); 

  printf("hello line 151\n");

  nhmax = nsink = mmax = mh_mass = 0;
  xCOM = yCOM = zCOM = vxCOM = vyCOM = vzCOM = ncount_doub = 0.0;
  xmax = ymax = zmax = 0;
  sl = tmax = h2max = gammin=masstot = dismax = omega0 = 0;
  nthresh = 2.e11;
  mass_conv = (1.e10/hubble)*Msun; 

  DeltaX = width * pc_cm / 2.0 / double(arrnum);
  dfac = 1. - .2*(1./double(arrnum));  
 
   //rmin = 100.;
   //rmax = 3e5;
   rmin = 0;
   rmax = width * pc_cm / au_cm / 2.0;
   nmin = 1.e-2;
   nmax = 1.e12;
   for(n=0; n<arrnum; n++)
     {
     m_enc_arr[n] = m_shell_arr[n] = num_arr[n] = num2_arr[n] = nh_arr[n] = temp_arr[n] = temp_nh_arr[n] = h2_arr[n] = elec_arr[n] = hd_arr[n] = 0.0;
     xmom_arr[n] = ymom_arr[n] = zmom_arr[n] = eint_arr[n] = gpot_arr[n] = 0;
     gpot_xL[n] = gpot_xR[n] = gpot_yL[n] = gpot_yR[n] = gpot_zL[n] = gpot_zR[n] = 0;
     num_xL[n] = num_xR[n] = num_yL[n] = num_yR[n] = num_zL[n] = num_zR[n] = 0;
     //d_arr[n] = pow(5.e8,double(n)/double(arrnum));
     d_arr[n] = rmin + (rmax - rmin)*double(n+1)/double(arrnum);
     //d_arr[n] = rmin * pow(rmax/rmin, double(n)/double(arrnum));
     n_arr[n] = nmin * pow(nmax/nmin, double(n)/double(arrnum));
     }

  //for(n=1; n<arrnum; n++)  d_arr[n] = (d_arr[n-1] + d_arr[n]) / 2.; 

  printf("hello line 167\n");

  for(n=0;n<Ngas;n++) 
     { 
          nh = P[n].nh;
          if(nh > nhmax && P[n].sink > -1)
            {
              //printf("Found the sink!\n");
              nhmax = nh;
              mmax=P[n].Mass;
              xmax=P[n].Pos[0];
              ymax=P[n].Pos[1];
              zmax=P[n].Pos[2];
              idmax = Id[n];
            }
       }

  printf("hello line 189\n");

  for(n=0;n<Ngas;n++) {
            nh = P[n].nh;
            disx = fabs((P[n].Pos[0]-xmax))*1.e3*Time/(0.7);
            disy = fabs((P[n].Pos[1]-ymax))*1.e3*Time/(0.7);
            disz = fabs((P[n].Pos[2]-zmax))*1.e3*Time/(0.7);
            if(fabs(disx) < width/2.0 && fabs(disy) < width/2.0 && fabs(disz) < width/2.0)
            //if(nh > nhmax/10.0)
            {
            vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
            vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
            vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
            xCOM = xCOM + P[n].Pos[0]*P[n].Mass;
            yCOM = yCOM + P[n].Pos[1]*P[n].Mass;
            zCOM = zCOM + P[n].Pos[2]*P[n].Mass;
            ncount_doub = ncount_doub + P[n].Mass;
 
           if(P[n].pot < pot_min)
             pot_min = P[n].pot;
           if(P[n].pot > pot_max)
             pot_max = P[n].pot;

           }
    }

         printf("pot_min = %lg, pot_max = %lg\n", pot_min, pot_max);
         for(n = 0; n < Ngas; n++)
            P[n].pot = P[n].pot - pot_min;

         vx = vxCOM/ncount_doub;
         vy = vyCOM/ncount_doub;
         vz = vzCOM/ncount_doub;

         sinkposx=xmax;
         sinkposy=ymax;
         sinkposz=zmax;

         //sinkposx=xCOM/ncount_doub;
         //sinkposy=yCOM/ncount_doub;
         //sinkposz=zCOM/ncount_doub;

         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g, nhmax = %15.11g, mmax = %15.11g\n", sinkposx, sinkposy, sinkposz, nhmax, mmax);
         printf("Ngas = %d, idmax = %d\n", Ngas, idmax);
 
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
             dis=dis*1.e3*Time/(hubble);
             disAU=dis*206264.806;


             if(disAU < 10.0)
               {
               mh_mass = mh_mass + P[i].Mass*1.e10/hubble;
               }
             }
         printf("mh_mass = %lg\n", mh_mass);

         //for(k=0; k<arrnum; k++)
         //   printf("m_enc_arr = %lg\n", m_enc_arr[k]);
 
         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

             if(n%100000 == 0) printf("n =%d, nh =%lg\n", n, P[n].nh);              

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
             dis_pc=dis*1.e3*Time/(hubble);
             disAU=dis_pc*206264.806;

             disx = (P[n].Pos[0]-sinkposx) * 1.e3*Time/(hubble) * 206264.806 - 0.5*rmax/double(arrnum);
             disy = (P[n].Pos[1]-sinkposy) * 1.e3*Time/(hubble) * 206264.806 - 0.5*rmax/double(arrnum);
             disz = (P[n].Pos[2]-sinkposz) * 1.e3*Time/(hubble) * 206264.806 - 0.5*rmax/double(arrnum);

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
               masstot = masstot + P[n].Mass/1.e-10/hubble;
               omega0 = omega0 + (P[n].Mass/1.e-10/hubble)*vrot/(disAU*1.5e8);
               }

             volume = pow(au_cm*d_arr[arrnum-1],3);
             for(k=0; k<arrnum; k++)
                //if(disAU < d_arr[k])
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k])
                   {
                   if(k == 0) vol_arr[k] = pow(au_cm*d_arr[k],3);
                   if(k > 0)  vol_arr[k] = pow(au_cm*d_arr[k],3);
                   m_enc_arr[k] = m_enc_arr[k] + P[n].Mass*1.e10/hubble;
                   //eint_arr[k] = eint_arr[k] + 1.5*P[n].Pres/P[n].Rho;   //sum the internal energy
                   eint_arr[k] = eint_arr[k] + P[n].U*P[n].Mass*mass_conv; 
                   xmom_arr[k] = xmom_arr[k] + P[n].Mass * (P[n].Vel[0])*1.e5*pow(Time, 0.5)*mass_conv;
                   ymom_arr[k] = ymom_arr[k] + P[n].Mass * (P[n].Vel[1])*1.e5*pow(Time, 0.5)*mass_conv;
                   zmom_arr[k] = zmom_arr[k] + P[n].Mass * (P[n].Vel[2])*1.e5*pow(Time, 0.5)*mass_conv;
                   }            
 

             //if(disAU < (d_arr[0]+d_arr[1])/2.)
             //if(disAU < d_arr[0])
             if(fabs(disx) < d_arr[0] && fabs(disy) < d_arr[0] && fabs(disz) < d_arr[0])
                   {
                   volume = (4*3.14159/3)*(pow(1.5e13*d_arr[k],3) - pow(1.5e13*d_arr[k-1],3));
                   num_arr[0] = num_arr[0] + 1.0;
                   nh_arr[0]  = nh_arr[0] + P[n].nh;
                   //nh_arr[0] = Msun*(m_enc_arr[k] - m_enc_arr[k-1])/volume/(1.22*1.67e-24);
                   temp_arr[0] = temp_arr[0] + P[n].Temp;
                   h2_arr[0]   = h2_arr[0] + P[n].H2I;
                   elec_arr[0] = elec_arr[0] + P[n].HII + P[n].HeII;
                   hd_arr[0]   = hd_arr[0] + P[n].HDI;
                   gpot_arr[0] = gpot_arr[0] + P[n].pot;

                   gpot_xR[0] = gpot_xR[0] + P[n].pot;
                   num_xR[0] = num_xR[0] + 1.0;
                   gpot_xL[0] = gpot_xL[0] + P[n].pot;
                   num_xL[0] = num_xL[0] + 1.0;
                   gpot_yR[0] = gpot_yR[0] + P[n].pot;
                   num_yR[0] = num_yR[0] + 1.0;
                   gpot_yL[0] = gpot_yL[0] + P[n].pot;
                   num_yL[0] = num_yL[0] + 1.0;
                   gpot_zR[0] = gpot_zR[0] + P[n].pot;
                   num_zR[0] = num_zR[0] + 1.0;
                   gpot_zL[0] = gpot_zL[0] + P[n].pot;
                   num_zL[0] = num_zL[0] + 1.0;
                   }

             for(k=1; k<arrnum; k++)
                //if(disAU < d_arr[k+1] && disAU > d_arr[k])
                //if(disAU < d_arr[k] && disAU > d_arr[k-1])
                //if(disAU < (d_arr[k+1]+d_arr[k])/2. && disAU > (d_arr[k-1]+d_arr[k])/2.)
                //if(disAU < d_arr[k+1] && disAU > d_arr[k-1])
                if(   fabs(disx) < d_arr[k]   && fabs(disy) < d_arr[k]   && fabs(disz) < d_arr[k] 
                  && (fabs(disx) > d_arr[k-1] || fabs(disy) > d_arr[k-1] || fabs(disz) > d_arr[k-1]))
                   {
                   volume = (4.*3.14159/3)*(pow(1.5e13*d_arr[k+1],3) - pow(1.5e13*d_arr[k],3)); 
                   num_arr[k] = num_arr[k] + 1.0; 
                   //nh_arr[k] = nh_arr[k] + P[n].nh;
                   nh_arr[k] = Msun*(m_enc_arr[k+1] - m_enc_arr[k])/volume/(1.22*1.67e-24);
                   temp_arr[k] = temp_arr[k] + P[n].Temp;
                   h2_arr[k] = h2_arr[k] + P[n].H2I;
                   elec_arr[k] = elec_arr[k] + P[n].HII /*+ P[n].HeII*/;
                   hd_arr[k] = hd_arr[k] + P[n].HDI;
                   gpot_arr[k] = gpot_arr[k] + P[n].pot;
                   }

              for(k=1; k<arrnum; k++)
                {
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disx > d_arr[k-1]  && fabs(disx-0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm )
                   {
                   gpot_xR[k] = gpot_xR[k] + P[n].pot;
                   num_xR[k] = num_xR[k] + 1.0; 
                   }                          
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disx < -d_arr[k-1] && fabs(disx+0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm)
                   {
                   gpot_xL[k] = gpot_xL[k] + P[n].pot;
                   num_xL[k] = num_xL[k] + 1.0;
                   }
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disy > d_arr[k-1]  && fabs(disy-0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm)
                   {
                   gpot_yR[k] = gpot_yR[k] + P[n].pot;
                   num_yR[k] = num_yR[k] + 1.0;
                   }
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disy < -d_arr[k-1] && fabs(disy+0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm)
                   {
                   gpot_yL[k] = gpot_yL[k] + P[n].pot;
                   num_yL[k] = num_yL[k] + 1.0;
                   } 
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disz > d_arr[k-1]  && fabs(disz-0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm)
                   {
                   gpot_zR[k] = gpot_zR[k] + P[n].pot;
                   num_zR[k] = num_zR[k] + 1.0;
                   }
                if(fabs(disx) < d_arr[k] && fabs(disy) < d_arr[k] && fabs(disz) < d_arr[k] && disz < -d_arr[k-1] && fabs(disz+0.5*(d_arr[k-1]+d_arr[k])) < 0.2*DeltaX/au_cm)
                   {
                   gpot_zL[k] = gpot_zL[k] + P[n].pot;
                   num_zL[k] = num_zL[k] + 1.0;
                   }
                 }

             if(P[n].nh > nmin && disAU > 0 && P[n].sink > -4)
               {
                for(k=1; k<arrnum; k++)
                  if(P[n].nh < n_arr[k] && P[n].nh > n_arr[k-1] && disAU < 1.3e5)
                     {
                     num2_arr[k] = num2_arr[k] + P[n].Mass;
                     temp_nh_arr[k] = temp_nh_arr[k] + P[n].Temp*P[n].Mass;
                     }
               }


            }

		  	  
   outfile=fopen(output_fname, "a");
   outfile2=fopen(output_fname2, "a");
   outfile3=fopen(output_fname3, "a");
   outfile4=fopen(output_fname4, "a");
   for(n=0; n<arrnum; n++)
      {
      //nh_arr[n] = nh_arr[n]/num_arr[n];
      temp_arr[n] = temp_arr[n]/num_arr[n];
      h2_arr[n] = h2_arr[n]/num_arr[n];
      elec_arr[n] = elec_arr[n]/num_arr[n];
      hd_arr[n] = hd_arr[n]/num_arr[n];
      gpot_arr[n] = gpot_arr[n]/num_arr[n];
      temp_nh_arr[n] = temp_nh_arr[n]/num2_arr[n];
      gpot_xR[n] = gpot_xR[n]/num_xR[n];
      gpot_xL[n] = gpot_xL[n]/num_xL[n];
      gpot_yR[n] = gpot_yR[n]/num_yR[n];
      gpot_yL[n] = gpot_yL[n]/num_yL[n];
      gpot_zR[n] = gpot_zR[n]/num_zR[n];
      gpot_zL[n] = gpot_zL[n]/num_zL[n];    

      printf("%lg, %lg, %lg, %lg, %lg, %lg \n", num_xR[n], num_xL[n], num_yR[n], num_yL[n], num_zR[n], num_zL[n]);
      }

   for(n=0; n<arrnum; n++)
      {

      if(n == 0)
        {
        gacc_xR[n] = fabs(gpot_xL[1] - gpot_xR[0]) / DeltaX;
        gacc_xL[n] = fabs(gpot_xR[1] - gpot_xL[0]) / DeltaX;
        gacc_yR[n] = fabs(gpot_yL[1] - gpot_yR[0]) / DeltaX;
        gacc_yL[n] = fabs(gpot_yR[1] - gpot_yL[0]) / DeltaX;
        gacc_zR[n] = fabs(gpot_zL[1] - gpot_zR[0]) / DeltaX;
        gacc_zL[n] = fabs(gpot_zR[1] - gpot_zL[0]) / DeltaX;
        }
      else
       {
        gacc_xR[n] = fabs(gpot_xR[n] - gpot_xR[n-1]) / DeltaX;
        gacc_xL[n] = fabs(gpot_xL[n] - gpot_xL[n-1]) / DeltaX;
        gacc_yR[n] = fabs(gpot_yR[n] - gpot_yR[n-1]) / DeltaX;
        gacc_yL[n] = fabs(gpot_yL[n] - gpot_yL[n-1]) / DeltaX;
        gacc_zR[n] = fabs(gpot_zR[n] - gpot_zR[n-1]) / DeltaX;
        gacc_zL[n] = fabs(gpot_zL[n] - gpot_zL[n-1]) / DeltaX;
       }


      fprintf(outfile, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n",  d_arr[n], nh_arr[n], temp_arr[n], h2_arr[n], elec_arr[n], hd_arr[n], m_enc_arr[n]);
      fprintf(outfile2, "%15.11g %15.11g \n",  n_arr[n], temp_nh_arr[n]); 
      fprintf(outfile3, "%15.11g %15.11g %15.11g %15.11g %15.11g\n",  eint_arr[n], xmom_arr[n], ymom_arr[n], zmom_arr[n], gpot_arr[n]);
      fprintf(outfile4, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", d_arr[n]*au_cm,  gpot_xL[n], gpot_xR[n], gpot_yL[n], gpot_yR[n], gpot_zL[n], gpot_zR[n],  gacc_xL[n], gacc_xR[n], gacc_yL[n], gacc_yR[n], gacc_zL[n], gacc_zR[n]);
     }
   fclose(outfile);
   fclose(outfile2);
   fclose(outfile3);
   fclose(outfile4);

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
  double h2frac, muh2, muh2in, pot_fac;
 
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

  pot_fac = (1. / G) * GRAVITY * (UnitMass_in_g / header1.HubbleParam) /  (UnitLength_in_cm * pow(header1.time/header1.HubbleParam,1)) ;

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
	  P[i].Rho= P[i].Rho * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  P[i].nh= P[i].Rho / MeanWeight;
          P[i].Pres = (gamma-1)*u*P[i].Rho;
          P[i].U = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

          //for ISOTHERMAL run
          //gamma = 1.0001; 
          //P[i].Temp = PROTONMASS/BOLTZMANN * u;  

	  /*  printf("zred = %g", zred);*/
	}

      if(wpot == 1)
        P[i].pot = P[i].pot * pot_fac;
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


         if(wpot == 1)
         {
         printf("Reading in potentials!\n");
         SKIP;
         for(k=0,pc_new=pc;k<6;k++)
           {
             for(n=0;n<header1.npart[k];n++)
               {
                 fread(&P[pc_new].pot, sizeof(double), 1, fd);  //gravitational potential
                 pc_new++;
               }
           }
         SKIP;
         printf("P[100].potential = %lg\n", P[100].dummy);
         }

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
              fread(&P[pc_sph].HII, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);    //DII
              fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);  //HeIII
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






  











