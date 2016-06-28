#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define hubble_param 0.7
#define PI 3.14159265359

#define readB 0
#define vpot 0
#define radprof 0
#define BFF 1
#define def_fac 0
#define exp_fac 1.0

#define width_small 0.1  //Orion2 box distance in pc
//#define ref_lev 128
#define ref_lev 32

#define const_gamma 1
#define which_sim 6


//#define snapbegin 41  //(which_sim==0)
//#define snapbegin 6601 //which_sim == 3
//#define snapbegin 6900 //which_sim == 3
//#define snapbegin 1 //which_sim == 4
//#define snapbegin 33
//#define snapbegin 1000
//#define snapbegin 7132 //which_sim == 3
//#define snapbegin 5979
//#define snapbegin 679 //which_sim == 4
//#define snapbegin 282 //which_sim == 5
#define snapbegin 445 //which_sim == 6

//#define snapend 5980  //(which_sim==0)
//#define snapend 7133 //which_sim == 3
//#define snapend 680 //which_sim == 4
//#define snapend 283 //which_sim == 5
#define snapend 445  //which_sim == 6

#define snapcheck  snapend
#define snapcenter snapend

#define snapnum (snapend - snapbegin + 1)

int load_snapshot(char *fname, int files);
int unit_conversion(void);
int do_what_you_want(void);
int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1);
int allocate_memory(void);
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

struct plist_data
{
  int Id;
  double disx, disy, disz;
}*plist;

struct particle_data 
{
  double  Pos[3], Pos_new[3], Pos_rot[3];
  double  Vel[3], Vel_new[3], Vel_rot[3];
  double  Mass;
  int    Type;
  int    Id;
  int shell;

  double  U, Temp, nh, Density, Rho, hsm;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
#if(const_gamma != 1)
  double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
#endif
  double dummy;
  double turb, soundspeed;
  double disx, disy, disz, hsm_pc;
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
  int  i, j, k, n, type, snapshot_number, files, random, ncount, ncounthalo1, ncount2, idmax;
  int    snaparr[snapnum];
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vxCOM=0, vyCOM=0, vzCOM=0, vel, vel_new, typemax, dismax; 
  double sq_error=0, rms_error=0, tot_error=0;
  double delfac=0.9, n0, b0, a0, a0_x, a0_y, n_anal, b_anal;
  double sl, masstot, temp, tmax, h2, h2max, gam, gammin; 
  double nh, nhmax, mass, mmax, dis, dis_pc, dis_new, xmax=0, ymax=0, zmax=0; 
  double sinkposx, sinkposy, sinkposz, disAU, vradx, vrady, vradz, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double disx, disy, disz, time_fac, cosmo_fac, cosmo_fac0;
  FILE *outfile, *outfile2, *infile;

  sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
  sprintf(basenameout, "turbbin_zoom10_new_cut");

  sprintf(path2, "/work/00863/minerva/orion/");
  sprintf(basename2, "bin_zoom10");

  int arrnum = 0, idnum = 0;
  int snapmin = 6, snapmax = 40, jinc = 5;

  if(which_sim == 0)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
   sprintf(basenameout, "turbbin_zoom10_new_cut");
   }
  if(which_sim == 1)
   {
   sprintf(path, "/scratch/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_zoom10_new_cut_ref");
   sprintf(basenameout, "turbbin_zoom10_new_cut_ref");
   }
  if(which_sim == 3)
   {
   sprintf(path, "/scratch/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_zoom10_new_cut_ref3");
   sprintf(basenameout, "turbbin_zoom10_new_cut_ref3");
   }
  if(which_sim == 4)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_HR10_ideal");
   sprintf(basename2,"bin_HR10_ideal");
   sprintf(basenameout, "turbbin_HR10_ideal");
   }
  if(which_sim == 5)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_MR10_ideal");
   sprintf(basename2,"bin_MR10_ideal");
   sprintf(basenameout, "turbbin_MR10_ideal");
   }
  if(which_sim == 6)
   {
   sprintf(path, "/scratch/00863/minerva");
   sprintf(path2, "/work/00863/minerva/orion");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_zoom10_ref4_nosmooth");
   sprintf(basenameout, "turbbin_zoom10_ref4");
   }


  for(n=0;n<snapnum;n++)
    snaparr[n] = snapbegin+(1*n);

/*
  for(n=0;n<snapnum;n++)
    printf("snapnum = %d, snaparr[%d] = %d\n", snapnum, n, snaparr[n]);
*/
///////////////////////////////get ID of particles whose evolution was followed////////////////
  for(j=snapcheck;j<=snapcheck;j=j+5){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);

  if(which_sim == 4 || which_sim == 5)
  {
    sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  }

  Ngas = load_snapshot(input_fname, files);

  unit_conversion();

  time_fac = Time;

  nhmax=0;
  for(i = 0; i < Ngas; i++)
      {
      if(P[i].nh > nhmax /*&& P[i].sink < 0.5*/)
           {
           nhmax = P[i].nh;
           xmax=P[i].Pos[0];
           ymax=P[i].Pos[1];
           zmax=P[i].Pos[2];
           }
       }
   printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);

   int iskip = 10;
 
  for(i = 0; i < Ngas; i++)
     {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%iskip == 0)
        idnum++;
     }

  printf("idnum = %d\n", idnum);

  if(!(plist=(struct plist_data *) malloc(idnum*sizeof(struct plist_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }

  n=0;
  for(i = 0; i < Ngas; i++)
      {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%iskip == 0)
        {
        P[i].to_print = 10;
        plist[n].Id = P[i].Id;
        n++;
        }
      }

free(P);
}
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////Use this section if we want central position to be taken from an earlier snapshot
  for(j=snapcenter;j<=snapcenter;j=j+5){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  if(which_sim == 4 || which_sim == 5)
    sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);

  Ngas = load_snapshot(input_fname, files);

  unit_conversion();

  time_fac = Time;
  cosmo_fac = time_fac/(hubble_param);

  nhmax=0;
  for(i = 0; i < Ngas; i++)
      {
      P[i].Pos[0] = P[i].Pos[0]*cosmo_fac;
      P[i].Pos[1] = P[i].Pos[1]*cosmo_fac;
      P[i].Pos[2] = P[i].Pos[2]*cosmo_fac;

      if(P[i].nh > nhmax /*&& P[i].sink < 0.5*/)
           {
           nhmax = P[i].nh;
           xmax=P[i].Pos[0];
           ymax=P[i].Pos[1];
           zmax=P[i].Pos[2];
           }
       }
   printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);

  free(P);
  }
////////////////////////////////////////////////////////////////////////////////////////////

  int nskip = 1000;
  for(j=0;j<snapnum;j=j+nskip){

  snapshot_number= snaparr[j];                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */
  arrnum = 0;

  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  if(which_sim == 4 || which_sim == 5)
    sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);

  sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);

  sprintf(output_fname2, "turb_comp_%04d", snapshot_number);
  

  Ngas = load_snapshot(input_fname, files);

 
  printf("t_Hubble = %lg \n", 5.4e8/pow((1.e0+header1.redshift)/10.e0, 1.5e0));
  printf("hi 1, Ngas = %d, NumPart = %d\n", Ngas, NumPart);

  unit_conversion();  

  time_fac = Time;
  //For NON-cosmological runs
  //time_fac = 1.0;
  cosmo_fac = time_fac/(hubble_param);

  printf("hi 2\n");

  outfile=fopen(output_fname, "w");

  ncount = 0;
  ncount2 = 0;

  printf("hi 3\n");


  nhmax = mmax = mh_mass = sl = tmax = h2max = gammin = masstot = dismax = vxCOM = vyCOM = vzCOM = 0;
  typemax = 5;

  printf("hi 4\n");

  for(n=0;n<Ngas;n++) { 

          nh = P[n].nh;
          mass=P[n].Mass;

          P[n].Pos[0] = P[n].Pos[0]*cosmo_fac;
          P[n].Pos[1] = P[n].Pos[1]*cosmo_fac;
          P[n].Pos[2] = P[n].Pos[2]*cosmo_fac;
          P[n].hsm = P[n].hsm*cosmo_fac;

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
              vx = P[n].Vel[0];
              vy = P[n].Vel[1];
              vz = P[n].Vel[2];
            }
    }
  printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);

  int iskip = 10;
  for(i = 0; i < Ngas; i++)
      {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%iskip == 0)
        {
        P[i].to_print = 10;
        }
      }


//////////////////////////////////////////////////////////////////////////////////////////////////////
//calculate disx, disy, and disz! 
  double ref_lev_doub, DeltaX, num[ref_lev], dis_arr[ref_lev]; 
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

  for(n=0;n<Ngas;n++) {
    P[n].disx = (P[n].Pos[0] - xmax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].disy = (P[n].Pos[1] - ymax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].disz = (P[n].Pos[2] - zmax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].hsm_pc = P[n].hsm * 1.e3;

   if(take_profile == 1)
    {
    P[n].disx = P[n].disx + DeltaX/2;
    P[n].disy = P[n].disy + DeltaX/2;
    P[n].disz = P[n].disz + DeltaX/2;
    }
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  double nh_typ=0, b_typ=0, u_typ=0, num_typ=0, bfield, ufield, u_norm, bx_typ=0, by_typ=0, bz_typ=0, ax_typ=0, ay_typ=0, az_typ=0, a_anal=0, a_typ;
  double dwdx, dwdy, dwdz, ax, ay, az, rcheck;
  double nmin = 1.e-1;

/*
  int dest;
  struct particle_data psave, psource;
  printf("start reordering\n");
  for(n=0;n<Ngas;n++)
      if(P[n].nh > nmin)
      for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
             {
             psource = P[n];
             psave   = P[k];
             P[n] = psave;
             P[k] = psource;  
             }
  printf("finished reordering\n");

  for(n=0;n<idnum;n++)
        if(P[n].nh > nmin)
          for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
            {
             arrnum++;
             P[n].Density = 0;
             P[n].to_print = 10;
            } //end if(P[n].Id == plist[k].Id) and Ngas for loop
*/

  printf("arrnum = %d\n", arrnum);

////////////////////////////////////////////////////////////////////////////////////////////////////////
//subtract out bulk velocities and find alternate particle coordinates!!
//
  for(n=0;n<Ngas;n++) {
    if(P[n].nh > nhmax/2)
      {
      masstot = masstot + P[n].Mass * 1.e10 / hubble_param;
      vxCOM = vxCOM + P[n].Vel[0] * P[n].Mass * 1.e10 / hubble_param;
      vyCOM = vyCOM + P[n].Vel[1] * P[n].Mass * 1.e10 / hubble_param;
      vzCOM = vzCOM + P[n].Vel[2] * P[n].Mass * 1.e10 / hubble_param;
      }
    }

  vxCOM = vxCOM/masstot; vyCOM = vyCOM/masstot; vzCOM = vzCOM/masstot;

  if(masstot <= 0)  {vxCOM = vx; vyCOM = vy; vzCOM = vz;}

  rotate(xmax, ymax, zmax, vxCOM, vyCOM, vzCOM);

  for(n=0;n<Ngas;n++) {
    P[n].Vel[0] = P[n].Vel[0] - vxCOM;
    P[n].Vel[1] = P[n].Vel[1] - vyCOM;
    P[n].Vel[2] = P[n].Vel[2] - vzCOM;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kfac, dens_conv, y_orion = -1*DeltaX/2, z_orion = -1*DeltaX/2;
  dens_conv = 1.e10 * 1.989e33 / hubble_param * pow(3.08567758e18,-3) / 1.2195 /  1.6726e-24;


         sinkposx=xmax;
         sinkposy=ymax;
         sinkposz=zmax;
         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g, nhmax = %15.11g, idmax = %d, mmax = %15.11g\n", sinkposx, sinkposy, sinkposz, nhmax, idmax, mmax);
 

         if(nhmax < 1.e0)
           {
           sinkposx=header1.BoxSize/2.0;
           sinkposy=header1.BoxSize/2.0;
           sinkposz=header1.BoxSize/2.0; 
           }
         printf("sinkposx = %15.11g, sinkposy = %15.11g, sinkposz = %15.11g\n", sinkposx, sinkposy, sinkposz);
         printf("vxCOM = %lg, vyCOM = %lg, vzCOM = %lg, vx = %lg, vy = %lg, vz = %lg, Time = %lg\n", vxCOM, vyCOM, vzCOM, vx, vy, vz, Time);

         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
             dis_pc = dis*1.e3*Time/(0.7);
             disAU  = dis_pc*206264.806;

             dis_new = pow(((P[n].Pos_new[0])*(P[n].Pos_new[0]) + (P[n].Pos_new[1])*(P[n].Pos_new[1]) + (P[n].Pos_new[2])*(P[n].Pos_new[2])), 0.5);

             vel = P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1] + P[n].Vel[2]*P[n].Vel[2];
             vel = pow(vel,0.5)*pow(Time, 0.5); 

             vel_new = P[n].Vel_new[0]*P[n].Vel_new[0] + P[n].Vel_new[1]*P[n].Vel_new[1] + P[n].Vel_new[2]*P[n].Vel_new[2];
             vel_new = pow(vel_new,0.5)*pow(Time, 0.5);

             if(n % 100000 == 0)
               printf("dis = %lg, dis_new = %lg, vel = %lg, vel_new = %lg, vx = %lg, vy = %lg, vz = %lg, vx_new = %lg, vy_new = %lg, vz_new = %lg\n", dis, dis_new, vel, vel_new, P[n].Vel[0], P[n].Vel[1], P[n].Vel[2], P[n].Vel_new[0], P[n].Vel_new[1], P[n].Vel_new[2]);

             vx = P[n].Vel[0]*pow(Time, 0.5);
             vy = P[n].Vel[1]*pow(Time, 0.5);
             vz = P[n].Vel[2]*pow(Time, 0.5);

             vrad =  (P[n].Vel[0]*(P[n].Pos[0]-sinkposx) + P[n].Vel[1]*(P[n].Pos[1]-sinkposy) + P[n].Vel[2]*(P[n].Pos[2]-sinkposz))/pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             vrad = vrad*pow(Time, 0.5); 

             vradx = vrad*(P[n].Pos[0]-sinkposx)/dis;
             vrady = vrad*(P[n].Pos[1]-sinkposy)/dis;
             vradz = vrad*(P[n].Pos[2]-sinkposz)/dis;

             vrotx = (P[n].Pos[1]-sinkposy)*P[n].Vel[2] - (P[n].Pos[2]-sinkposz)*P[n].Vel[1];
             vroty = (P[n].Pos[2]-sinkposz)*P[n].Vel[0] - (P[n].Pos[0]-sinkposx)*P[n].Vel[2];
             vrotz = (P[n].Pos[0]-sinkposx)*P[n].Vel[1] - (P[n].Pos[1]-sinkposy)*P[n].Vel[0];

             vrotx = vrotx*pow(Time, 0.5)/dis;   //convert to km/s
             vroty = vroty*pow(Time, 0.5)/dis;
             vrotz = vrotz*pow(Time, 0.5)/dis;
             
             vrot = pow(vrotx*vrotx + vroty*vroty + vrotz*vrotz, 0.5); 

             ////////////////////alternate vrot formulation/////////////
             vrot = P[n].Pos_new[0]*P[n].Vel_new[1] - P[n].Pos_new[1]*P[n].Vel_new[0];
             vrot = vrot / pow(P[n].Pos_new[0]*P[n].Pos_new[0] + P[n].Pos_new[1]*P[n].Pos_new[1], 0.5);
             vrot = vrot*pow(Time, 0.5);

             vrotx = vrot*(P[n].Pos[0]-sinkposx)/dis;
             vroty = vrot*(P[n].Pos[1]-sinkposy)/dis;
             vrotz = vrot*(P[n].Pos[2]-sinkposz)/dis;
             ///////////////////////////////////////////////////////////

             P[n].turb = pow(vx - vradx - vrotx,2) + pow(vy - vrady - vroty,2) + pow(vz - vradz - vrotz,2);
             P[n].turb = pow(P[n].turb, 0.5);

             P[n].soundspeed = pow(1.38e-16*P[n].Temp/1.e-24,0.5)/1.e5;

             random = rand();
             int to_print = 0;
             if(random*(1.e0/RAND_MAX)< 0.002e0) to_print++;

            if(to_print > 0)
              {
             #if(const_gamma == 0)
              fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
              P[n].soundspeed, P[n].Id,x,y,z,P[n].Temp,P[n].nh,P[n].HII,P[n].H2I,P[n].HDI, P[n].gam, disAU, vrad, vrot, P[n].hsm_pc, P[n].Mass, P[n].Vel[0]*pow(Time, 0.5), P[n].Vel[1]*pow(Time, 0.5), P[n].Vel[2]*pow(Time, 0.5), P[n].turb);
              #else
               fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
              P[n].soundspeed, P[n].Id,x,y,z,P[n].Temp,P[n].nh,P[n].dummy,P[n].dummy,P[n].dummy, P[n].dummy, disAU, vrad, vrot, P[n].hsm_pc, P[n].Mass, P[n].Vel[0]*pow(Time, 0.5), P[n].Vel[1]*pow(Time, 0.5), P[n].Vel[2]*pow(Time, 0.5), P[n].turb);
               #endif
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
#if(const_gamma==0)
       h2frac=2.0*P[i].H2I;
#else
       h2frac = 1.e-3;
#endif
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

#if(const_gamma)
	  gamma= 5.0/3.0;
	  //gamma = 1.1;
#else
          gamma=P[i].gam;	 
#endif

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Density * UnitDensity_in_cgs;
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
	      fread(&P[pc_sph].Density, sizeof(double), 1, fd);
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

#if(const_gamma == 0)
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {

              fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
              fread(&P[pc_sph].HII, sizeof(double), 1, fd);
              fread(&P[pc_sph].DII, sizeof(double), 1, fd);
              fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
              fread(&P[pc_sph].HeIII, sizeof(double), 1, fd);

/*
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
*/
              pc_sph++;
            }
          SKIP;


          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].gamma, sizeof(double), 1, fd);
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
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].hsm_pc);
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
    dkernel = (8./3.14159/pow(hsm,3)) * ( -12*pow(ratio,1)*(1./hsm)*drdu + 18*pow(ratio,2)*(1./hsm)*drdu );
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



int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1)
{
int shell_num = 1000;
double shell_num_doub;
double dis, disAU, vec[3][shell_num], r, alpha, beta, d1[3][3], d2[3][3], d12[3][3][shell_num];
double dmax=0, i_doub, rshell[shell_num+1];
int i, j, k, n, ks;

for(ks=0; ks<shell_num; ks++)  
   vec[0][ks] = vec[1][ks] = vec[2][ks] = 0;

for(i=0; i<NumPart; i++)
   {
   P[i].shell = 0;

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
   if(P[i].to_print > 0)
      {
      dis = pow(P[i].Pos_rot[0]*P[i].Pos_rot[0] + P[i].Pos_rot[1]*P[i].Pos_rot[1] + P[i].Pos_rot[2]*P[i].Pos_rot[2], 0.5);
      dis=dis*1.e3*Time/(0.7);                   //dis is in pc
      disAU = dis*206264.8060;
      if(disAU > dmax) dmax = disAU;
      }
   }

shell_num_doub = (double) shell_num;
for(ks=0; ks<shell_num+1; ks++)
  {
  i_doub = (double) ks;
  rshell[ks] = dmax * i_doub / shell_num_doub;
  printf("rshell[%d] = %lg \n", ks, rshell[ks]);
  }


for(i=0; i<NumPart; i++)
  {
  dis = pow(P[i].Pos_rot[0]*P[i].Pos_rot[0] + P[i].Pos_rot[1]*P[i].Pos_rot[1] + P[i].Pos_rot[2]*P[i].Pos_rot[2], 0.5);
  dis=dis*1.e3*Time/(0.7);                   //dis is in pc
  disAU = dis*206264.8060;

  for(ks=0; ks<shell_num; ks++)
  if(disAU < rshell[ks+1] && disAU > rshell[ks])
    {
    P[i].shell = ks;
    vec[0][ks] = vec[0][ks] + P[i].Mass * (P[i].Pos_rot[1] *  P[i].Vel_rot[2] - P[i].Pos_rot[2] *  P[i].Vel_rot[1]);
    vec[1][ks] = vec[1][ks] + P[i].Mass * (P[i].Pos_rot[2] *  P[i].Vel_rot[0] - P[i].Pos_rot[0] *  P[i].Vel_rot[2]);
    vec[2][ks] = vec[2][ks] + P[i].Mass * (P[i].Pos_rot[0] *  P[i].Vel_rot[1] - P[i].Pos_rot[1] *  P[i].Vel_rot[0]);
    }
  }

for(ks=0; ks<shell_num; ks++)
{
r = pow(vec[0][ks],2) + pow(vec[1][ks],2) + pow(vec[2][ks],2);
r = pow(r, 0.5);

for(i=0; i<=2; i++)
  {
  vec[i][ks] = vec[i][ks]/r;
  vec[i][ks] <= 1;
  vec[i][ks] >= -1;
  }

alpha = acos(vec[2][ks] / pow(pow(vec[1][ks],2) + pow(vec[2][ks],2), 0.5));
beta = asin(vec[0][ks]);

printf("alpha = %lg, beta = %lg \n", alpha, beta);

for(i = 0; i<=2; i++)
  for(j = 0; j<=2; j++)
    {
    d1[i][j] = 0;
    d2[i][j] = 0;
    d12[i][j][ks] = 0;
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
      d12[i][j][ks] += d2[i][k] * d1[k][j];


//printf("d1 = %lg, d1 = %lg, d1 = %lg, d2 = %lg, d2 = %lg, d2 = %lg\n", d1[1][1] , d1[1][2], d1[2][1], d2[0][0], d2[0][2], d2[2][0]);
//printf("d2 = %lg, d2 = %lg d2 = %lg d2 = %lg\n", d12[0][2][ks],  d12[0][0][ks],  d12[1][1][ks], d12[2][2][ks]);
}

for(n=0; n<NumPart; n++)
  {
    for(i = 0; i<=2; i++)
      for(j = 0; j<=2; j++)
        P[n].Pos_new[i] += d12[i][j][P[n].shell] * P[n].Pos_rot[j];

    for(i = 0; i<=2; i++)
      for(j = 0; j<= 2; j++)
        P[n].Vel_new[i] += d12[i][j][P[n].shell] * P[n].Vel_rot[j];

    if(n%100000 == 0)
      {
      printf("P[n].Pos_rot[0] = %lg, P[n].Pos_rot[1] = %lg, P[n].Pos_rot[2] = %lg\n", P[n].Pos_rot[0], P[n].Pos_rot[1], P[n].Pos_rot[2]);
      printf("P[n].Pos_new[0] = %lg, P[n].Pos_new[1] = %lg, P[n].Pos_new[2] = %lg\n", P[n].Pos_new[0], P[n].Pos_new[1], P[n].Pos_new[2]);
      }
  }
}


  











