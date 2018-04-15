#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define readB 1
#define vpot 0
#define hubble_param 0.7
#define PI 3.14159265359

#define width_small 1.0  //Orion2 box distance in pc
#define ref_lev 128

#define which_sim 3

//#define snapbegin 40  //(which_sim==0)
#define snapbegin 6601 //which_sim == 3
//#define snapbegin 6900 //which_sim == 3
//#define snapbegin 7130 //which_sim == 3

//#define snapend 5980  //(which_sim==0)
//#define snapend 7100 //which_sim == 3
#define snapend 7133

#define snapnum (snapend - snapbegin)

int load_snapshot(char *fname, int files);
//int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
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
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  int    Id;

  double  U, Temp, nh, Density, Rho, hsm;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
   double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy;
  double nh_test, Bfieldx, Bfieldy, Bfieldz, error;
  double Afieldx, Afieldy, Afieldz, fcor;
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
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax, sq_error=0, rms_error=0, tot_error=0;
  double delfac=0.9, n0, b0, a0, a0_x, a0_y, n_anal, b_anal;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double disx, disy, disz, time_fac, cosmo_fac, cosmo_fac0;
  FILE *outfile, *outfile2, *infile;

  sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
  sprintf(basenameout, "snapbin_zoom10_new_cut");

  sprintf(path2, "/work/00863/minerva/orion/");
  if(vpot == 1 && which_sim == 3) sprintf(path2, "/work/00863/minerva/orion/bfield_comp_vpot");
  sprintf(basename2, "bin_zoom10");

  int arrnum = 0, idnum = 0;
  int snapmin = 6, snapmax = 40, jinc = 5;

  if(which_sim == 0)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
   sprintf(basenameout, "snapbin_zoom10_new_cut");
   }
  if(which_sim == 1)
   {
   sprintf(path, "/scratch/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_zoom10_new_cut_ref");
   sprintf(basenameout, "snapbin_zoom10_new_cut_ref");
   }
  if(which_sim == 3)
   {
   sprintf(path, "/scratch/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_zoom10_new_cut_ref3");
   sprintf(basenameout, "snapbin_zoom10_new_cut_ref3");
   }

  for(n=0;n<snapnum;n++)
    snaparr[n] = snapbegin+(1*n);

///////////////////////////////get ID of particles whose evolution was followed////////////////
  int snapcheck = snapend;
  for(j=snapcheck;j<=snapcheck;j=j+5){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  if(j<=999)
  {
  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  if(j>999)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  if(j>9999 || which_sim == 3)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
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
        plist[n].Id = P[i].Id;
        n++;
        }
      }

free(P);
}
///////////////////////////////////////////////////////////////////////////////////////////////
  int nskip = 50;
  for(j=0;j<=snapnum;j=j+nskip){

  snapshot_number= snaparr[j];                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */
  arrnum = 0;

  if(j<=999)
  {
  sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  if(j>999)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  if(j>9999 || which_sim == 3)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  if(which_sim == 3 && vpot == 1)
    sprintf(output_fname, "%s/%s_vpot_%04d", pathout, basenameout, snapshot_number);

  sprintf(output_fname2, "bfield_comp_%04d", snapshot_number);

  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */
 
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

//Read in B-field info!///////////////////////////////////////////////////////////
    sprintf(input_fname2, "%s/%s_bfield_%04d", path2, basename2, snapshot_number);
    if(which_sim == 1)
      sprintf(input_fname2, "%s/%s_bfield_ref_%04d", path2, basename2, snapshot_number);
    if(which_sim == 3)
      sprintf(input_fname2, "%s/%s_bfield_ref3_%04d", path2, basename2, snapshot_number);
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }
    
    printf("reading bfield file!\n");

if(vpot != 1)
{
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldx, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldy, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldz, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile);
}
if(vpot == 1)
{
    for(n=0;n<Ngas;n++)
        P[n].Bfieldx = P[n].Bfieldy = P[n].Bfieldz = 0;
    for(n=0;n<Ngas;n++)
        P[n].fcor = 1.0;
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldx, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldy, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldz, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile);
}
    fclose(infile);

//////////////////////////////////////////////////////////////////////////////////

  nhmax = 0;
  mmax = mh_mass = 0;
  xmax = ymax = zmax = 0;
  sl = tmax = h2max = gammin=0;
  masstot = dismax = 0;
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

          if(nh > nhmax && P[n].sink > -1)
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
 
  double ref_lev_doub, DeltaX, num[ref_lev], dis_arr[ref_lev], Bx[ref_lev], By[ref_lev], Bz[ref_lev], Bmag[ref_lev], rho[ref_lev];
  double xgridL, ygridL, zgridL, xgridR, ygridR, zgridR;
  
  xgridL = ygridL = zgridL = xgridR = ygridR = zgridR = -width_small/2.0;

  ref_lev_doub = (double) ref_lev;
  DeltaX = width_small / ref_lev_doub;
  printf("DeltaX = %lg\n", DeltaX);

  for(i=0;i<ref_lev;i++) {
  num[i] = Bx[i] = By[i] = Bz[i] = Bmag[i] = rho[i] = 0; 
  dis_arr[i] = xgridL + DeltaX/2;
  xgridL = xgridL + DeltaX;
  printf("dis_arr[%d] = %lg\n", i, dis_arr[i]); 
  }

  for(n=0;n<Ngas;n++) {
    P[n].disx = (P[n].Pos[0] - xmax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].disy = (P[n].Pos[1] - ymax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].disz = (P[n].Pos[2] - zmax) * 1.e3 - DeltaX/2; //distance from densest point in parsecs
    P[n].hsm_pc = P[n].hsm * 1.e3;
    }

  double nh_typ=0, b_typ=0, u_typ=0, num_typ=0, bfield, ufield, u_norm, bx_typ=0, by_typ=0, bz_typ=0, ax_typ=0, ay_typ=0, az_typ=0, a_anal=0, a_typ;
  double dwdx, dwdy, dwdz, ax, ay, az, rcheck;

  for(n=0;n<Ngas;n++)
          //if(2. * (P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test) > 0.01)
          for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
            {
             arrnum++;
             P[n].Density = 0;
             P[n].to_print = 10;

             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {

               if(i == n) continue;
               rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue;
               P[n].Density = P[n].Density + P[i].Mass * calc_kernel_spline(n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
                }


             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {
               if(i == n) continue;
               rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue;
               P[n].fcor = P[n].fcor + (P[n].hsm /  3. / P[n].Density) * P[i].Mass * del_rho(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
               dwdx = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[0] - P[i].Pos[0]);
               dwdy = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[1] - P[i].Pos[1]);
               dwdz = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[2] - P[i].Pos[2]);
               ax = P[n].Afieldx - P[i].Afieldx;
               ay = P[n].Afieldy - P[i].Afieldy;
               az = P[n].Afieldz - P[i].Afieldz;
               P[n].Bfieldx = P[n].Bfieldx + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(0, ax, ay, az, dwdx, dwdy, dwdz);
               P[n].Bfieldy = P[n].Bfieldy + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(1, ax, ay, az, dwdx, dwdy, dwdz);
               P[n].Bfieldz = P[n].Bfieldz + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(2, ax, ay, az, dwdx, dwdy, dwdz);
               }

             nh_typ = nh_typ + P[n].nh;
             bfield = pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);
             ufield = pow(bfield,2) / (8.*3.14159);
             b_typ  = b_typ + bfield;
             u_typ = u_typ + ufield; 

             bx_typ = bx_typ + P[n].Bfieldx;
             by_typ = by_typ + P[n].Bfieldy;
             bz_typ = bz_typ + P[n].Bfieldz;
             ax_typ = ax_typ + P[n].Afieldx;
             ay_typ = ay_typ + P[n].Afieldy;
             az_typ = az_typ + P[n].Afieldz;

             num_typ = num_typ + 1.0;

             printf("n = %d x = %lg, y = %lg, z = %lg, bfieldx = %lg, bfieldy = %lg, bfieldz = %lg  ax = %lg, ay = %lg az = %lg, fcor = %lg\n", n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].Afieldx, P[n].Afieldy, P[n].Afieldz, P[n].fcor);
            }

  printf("arrnum = %d\n", arrnum);


  nh_typ = nh_typ / num_typ;
  b_typ  = bz_typ / num_typ;
  u_typ  = u_typ / num_typ;
  u_norm = u_typ / pow(nh_typ*1.67e-24, 1.3333);

  bx_typ  = bx_typ / num_typ;
  by_typ  = by_typ / num_typ;
  bz_typ  = bz_typ / num_typ;
  ax_typ  = ax_typ / num_typ;
  ay_typ  = ay_typ / num_typ;
  az_typ  = az_typ / num_typ;

  //b_typ = pow(bx_typ*bx_typ + by_typ*by_typ + bz_typ*bz_typ, 0.5);
  a_typ = pow(ax_typ*ax_typ + ay_typ*ay_typ + az_typ*az_typ, 0.5);

  if(j == nskip)
    {
    cosmo_fac0 = cosmo_fac;
    n0 = nh_typ;
    b0 = bz_typ;
    a0 = a_typ;
    a0_x = ax_typ;
    a0_y = ay_typ;
    }

  b_anal = b0 * pow(cosmo_fac0/cosmo_fac,2);
  a_anal = a0 * pow(cosmo_fac0/cosmo_fac,1); 

  printf("nh_typ = %lg, nh_anal = %lg, b_typ = %lg, b_anal = %lg, u_typ = %lg, u_norm = %lg a_typ = %lg, a_anal = %lg\n", nh_typ, n_anal, b_typ, b_anal, u_typ, u_norm, a_typ, a_anal);

  for(n=0;n<Ngas;n++)
     {
     if(P[n].disy-P[n].hsm_pc < -DeltaX/2 && P[n].disy+P[n].hsm_pc > -DeltaX/2 && P[n].disz-P[n].hsm_pc < -DeltaX/2 && P[n].disz+P[n].hsm_pc > -DeltaX/2)
     for(i=0;i<ref_lev;i++)
       {
       if(P[n].disx-P[n].hsm_pc < dis_arr[i] && P[n].disx+P[n].hsm_pc > dis_arr[i])
          {
          num[i]  = num[i] + 1.;
          rho[i]  = rho[i] + P[n].nh;
          Bx[i]   = Bx[i] + P[n].Bfieldx;
          By[i]   = By[i] + P[n].Bfieldy;
          Bz[i]   = Bz[i] + P[n].Bfieldz;  
          Bmag[i] = Bmag[i] + pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);
          }
       }
     }

  outfile2=fopen(output_fname2, "a");
  for(i=0;i<ref_lev;i++)
    {
    rho[i] = rho[i] / num[i];
    Bx[i] = Bx[i] / num[i];
    By[i] = By[i] / num[i];
    Bz[i] = Bz[i] / num[i];
    Bmag[i] = Bmag[i] / num[i];
    fprintf(outfile2, "%g %g %g %g %g %g\n", dis_arr[i], rho[i], Bx[i], By[i], Bz[i], Bmag[i]);
    printf("i = %d, num = %lg, rho = %lg, Bx = %lg, By = %lg, Bz = %lg, Bmag = %lg\n", i, num[i], rho[i], Bx[i], By[i], Bz[i], Bmag[i]);
    }
  //fprintf(outfile2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
  //zred, nh_typ, n_anal, b_typ, b_anal, a_typ, a_anal, u_typ, u_norm, bx_typ, by_typ, bz_typ, ax_typ, ay_typ, az_typ);
  fclose(outfile2);

  printf("hi 5\n");

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
         printf("vx = %lg, vy = %lg, vz = %lg\n", vx, vy, vz);

         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             P[n].error = 2. * (P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test); 

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
             dis=dis*1.e3*Time/(0.7);
             disAU=dis*206264.806;

             if(dis < 0.5 && fabs(P[n].error > 0.05))
               {
               sq_error = sq_error + P[n].error*P[n].error;
               tot_error++;
               }

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

             if(P[n].nh > 1.e0 && P[n].Temp > 1.e4 && P[n].sink == 0)
               printf("ID = %d, nh = %lg, temp = %lg, elec %lg, H2 = %lg HeII = %lg mass = %lg\n", P[n].Id, P[n].nh, P[n].Temp, P[n].HII, P[n].H2I, P[n].HeII, P[n].Mass);


    random = rand();
    int to_print = 0;
    //if(random*(1.e0/RAND_MAX)< 0.1e0) to_print++;
    if(fabs(P[n].Bfieldx) > 0 || fabs(P[n].Bfieldy) > 0 || fabs(P[n].Bfieldz) > 0) to_print++;


    if(P[n].to_print > 0)
    //if(P[n].nh > 1.e0 /*&& random*(1.e0/RAND_MAX)< 0.1e0*/ /*|| P[n].sink > 0.5 || (random*(1.e0/RAND_MAX))< 0.02e0*/) 
      {
      fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
      P[n].error, P[n].Id,x,y,z,P[n].Temp,P[n].nh,P[n].HII,P[n].H2I,P[n].HDI, P[n].gam, disAU, vrad, vrot, P[n].hsm*Time/0.7, P[n].Mass, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].nh_test);
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

double calc_kernel_spline(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double cosmo_fac)
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




  











