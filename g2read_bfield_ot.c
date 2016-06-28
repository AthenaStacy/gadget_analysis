#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define readB 1
#define vpot 0
#define hubble_param 0.7
#define PI 3.14159265359

#define def_fac 0
#define exp_fac 1.0

#define width_small 1.0
//#define width_small 0.1
//#define width_small 0.05  //Orion2 box distance in pc
//#define width_small (1.3/2/2)

#define width_outer 1.0

//#define ref_lev 128
#define ref_lev 32

#define which_sim 1

#define snapbegin 500  //(which_sim==1)
//#define snapbegin 1100  //(which_sim==1)
//#define snapbegin 2266  //(which_sim==1)

//#define snapbegin 1100  //(which_sim==2)
//#define snapbegin 2300  //(which_sim==2)
//#define snapbegin 4712  //(which_sim==2)

//#define snapbegin 1700  //(which_sim==3)
//#define snapbegin 3300  //(which_sim==3)
//#define snapbegin 7097  //(which_sim==3)

#define snapend snapbegin+1  //(which_sim==0)

#define snapcheck  snapbegin
#define snapcenter snapbegin

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

struct plist_data
{
  int Id;
  double disx, disy, disz;
}*plist;

struct plist_outer_data
{
  int Id;
  double disx, disy, disz;
}*plist_outer;

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
  double Afield[3], fcor;
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
  int  i, j, k, m, n, type, snapshot_number, files, random, ncount, ncounthalo1, ncount2, idmax;
  int    snaparr[snapnum];
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vxCOM=0, vyCOM=0, vzCOM=0, vel, typemax, dismax; 
  double sq_error=0, rms_error=0, tot_error=0;
  double delfac=0.9, n0, b0, a0, a0_x, a0_y, n_anal, b_anal;
  double nh, nhmax, mass, mmax, dis, xmax=0, ymax=0, zmax=0, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double disx, disy, disz, time_fac, cosmo_fac, cosmo_fac0;
  FILE *outfile, *outfile2, *infile;

  sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
  sprintf(basenameout, "snapbin_zoom10_new_cut");

  sprintf(path2, "/work/00863/minerva/orion/");
  sprintf(basename2, "ot");

  int arrnum = 0, idnum = 0, idnum_outer = 0;
  int snapmin = 6, snapmax = 40, jinc = 5;

  if(which_sim == 1)
   {
   sprintf(path, "/scratch/00863/minerva/ot_twodim");
   sprintf(path2, "/work/00863/minerva/orion/ot_twodim");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "ot");
   sprintf(basenameout, "snapot");
   }
  if(which_sim == 2)
   {
   sprintf(path, "/scratch/00863/minerva/ot_twodim_midres");
   sprintf(path2, "/work/00863/minerva/orion/ot_twodim_midres");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "ot");
   sprintf(basenameout, "snapot");
   }
  if(which_sim == 3)
   {
   sprintf(path, "/scratch/00863/minerva/ot_twodim_highres");
   sprintf(path2, "/work/00863/minerva/orion/ot_twodim_highres");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "ot");
   sprintf(basenameout, "snapot");
   }
  if(which_sim == 5)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_MR10_ideal");
   sprintf(basename2,"bin_MR10_ideal");
   sprintf(basenameout, "snapbin_MR10_ideal");
   }

  if(vpot == 1) sprintf(path2, "/work/00863/minerva/orion/bfield_comp_vpot_bwards");
  //if(vpot == 1) sprintf(path2, "/work/00863/minerva/orion/bfield_comp_vpot");

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

   int iskip;
   if(which_sim == 3) iskip = 5;
   if(which_sim == 0) iskip = 1; 

  for(i = 0; i < Ngas; i++)
     {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%iskip == 0)
        idnum++;
     if(disx < width_outer/2.0 && disy < width_outer/2.0 && disz < width_outer/2.0)
        idnum_outer++;
     }

  printf("idnum = %d\n", idnum);
  printf("idnum_outer = %d\n", idnum_outer);

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

/////////////////////////////////////////////////////////////////////////////////////////////////
//find outer circumference particles to include in calculation for inner circumference particles

  if(!(plist_outer=(struct plist_outer_data *) malloc(idnum_outer*sizeof(struct plist_outer_data))))
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
     if(disx < width_outer/2.0 && disy < width_outer/2.0 && disz < width_outer/2.0)
        {
        plist_outer[n].Id = P[i].Id;
        n++;
        }
      }
///////////////////////////////////////////////////////////////////////////////////////////////

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
  if(vpot == 1)  sprintf(output_fname, "%s/%s_vpot_%04d", pathout, basenameout, snapshot_number);

  sprintf(output_fname2, "bfield_comp_%04d", snapshot_number);
  if(vpot == 1)  sprintf(output_fname2, "bfield_comp_vpot_%04d", snapshot_number);
  

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

//Read in B-field info!///////////////////////////////////////////////////////////
if(readB == 1)
{
    sprintf(input_fname2, "%s/%s_bfield_%04d", path2, basename2, snapshot_number);
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }
    
    printf("reading bfield file %s!\n", input_fname2);

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
        fread(&P[n].Afield[0], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afield[1], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afield[2], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile);
}
    fclose(infile);
}
//////////////////////////////////////////////////////////////////////////////////

  nhmax = mmax = mh_mass = sl = tmax = h2max = gammin = masstot = dismax = 0;
  typemax = 5;

  printf("hi 4\n");

  for(n=0;n<Ngas;n++) { 

          nh = P[n].nh;
          mass=P[n].Mass;

          P[n].Pos[0] = P[n].Pos[0]*cosmo_fac*exp_fac;
          P[n].Pos[1] = P[n].Pos[1]*cosmo_fac*exp_fac;
          P[n].Pos[2] = P[n].Pos[2]*cosmo_fac*exp_fac;
          P[n].hsm = P[n].hsm*cosmo_fac*hfac_out*exp_fac;

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
  printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);
 
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

    dis = pow(P[n].disx*P[n].disx + P[n].disy*P[n].disy + P[n].disz*P[n].disz ,0.5);

    if(dis < width_small)
      {
      masstot = masstot + P[n].Mass * 1.e10 / hubble_param;
      vxCOM = vxCOM + P[n].Vel[0] * P[n].Mass * 1.e10 / hubble_param;
      vyCOM = vyCOM + P[n].Vel[1] * P[n].Mass * 1.e10 / hubble_param;
      vzCOM = vzCOM + P[n].Vel[2] * P[n].Mass * 1.e10 / hubble_param;
      }
    }

  vxCOM = vxCOM/masstot; vyCOM = vyCOM/masstot; vzCOM = vzCOM/masstot;

  for(n=0;n<Ngas;n++) {
    P[n].Vel[0] = P[n].Vel[0] - vxCOM;
    P[n].Vel[1] = P[n].Vel[1] - vyCOM;
    P[n].Vel[2] = P[n].Vel[2] - vzCOM;
    }

  double nh_typ=0, b_typ=0, u_typ=0, num_typ=0, bfield, ufield, u_norm, bx_typ=0, by_typ=0, bz_typ=0, ax_typ=0, ay_typ=0, az_typ=0, a_anal=0, a_typ;
  double ax_diff, ay_diff, az_diff, rcheck, dAxdy, dAxdz, dAydx, dAydz, dAzdx, dAzdy;
  double nmin = 1.e-2;
  double dwdu[3][neighb_num]; 
  double pos[3][neighb_num], Afield[3][neighb_num];
  double jacob[3][3], jacob_inv[3][3], det_jacob, jacob2[3][3], jacob_fin[3][3];

  int dest, p, q;
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

  double hsm_max=0.0;
  for(n=0;n<idnum;n++)
    {
    if(P[n].hsm > hsm_max)
      {
      hsm_max = P[n].hsm;
      nmin = P[n].nh;
      }
    }
  hsm_max = hsm_max * 1.e3 * Time / 0.7;
  nmin = 0.5*nmin; 
  printf("hsm_max = %lg pc, nmin = %lg\n", hsm_max, nmin);

  for(n=0;n<idnum;n++)
        if(P[n].nh > nmin)
          for(m=0;m<idnum;m++)
          if(P[n].Id == plist[m].Id)
            {
             arrnum++;
             P[n].Density = 0;
             P[n].fcor = 1.0;
             P[n].to_print = 10;

             for(k=0;k<neighb_num;k++) 
               {
               for(p=0;p<3;p++)
                 {
                 dwdu[p][k]= 0; 
                 pos[p][k] = 0;
                 Afield[p][k] = 0;
                 }
               }

             k=0;

             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {
               if(i == n) continue;
               if(P[i].nh < nmin) continue;
               rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue;
               if(k >= neighb_num) continue;              

               for(p=0;p<3;p++) 
                 {
                 dwdu[p][k] = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[p] - P[i].Pos[p]);
                 pos[p][k] = P[i].Pos[p];
                 Afield[p][k] = P[i].Afield[p];
                 }      

               P[n].Density = P[n].Density + P[i].Mass * calc_kernel_spline(n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
                k++;
                }

             if(vpot == 1)
             for(k=0;k<neighb_num;k++)
               {
               if(pos[0][k] == 0) continue;
               P[n].fcor = P[n].fcor + (P[n].hsm /  3. / P[n].Density) * P[1].Mass * del_rho(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], pos[0][k],  pos[1][k], pos[2][k], P[n].hsm, 0);
               }

             for(p=0; p<3; p++)
               for(q=0; q<3; q++)
                  jacob[p][q] = jacob_inv[p][q] = jacob_fin[p][q] = jacob2[p][q] = 0;

             for(k=0;k<neighb_num;k++)
               {
               for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                    {
                    jacob[p][q] = jacob[p][q]   + P[1].Mass*(pos[p][k]-P[n].Pos[p])*dwdu[q][k];
                    jacob2[p][q] = jacob2[p][q] - P[1].Mass*(Afield[p][k]-P[n].Afield[p])*dwdu[q][k];
                    }
               }

              det_jacob = calc_det(jacob[0][0], jacob[0][1], jacob[0][2], jacob[1][0], jacob[1][1], jacob[1][2], jacob[2][0], jacob[2][1], jacob[2][2]);

              jacob_inv[0][0] = jacob[1][1]*jacob[2][2] - jacob[1][2]*jacob[2][1];
              jacob_inv[0][1] = jacob[0][2]*jacob[2][1] - jacob[0][1]*jacob[2][2];
              jacob_inv[0][2] = jacob[0][1]*jacob[1][2] - jacob[0][2]*jacob[1][1];

              jacob_inv[1][0] = jacob[1][2]*jacob[2][0] - jacob[1][0]*jacob[2][2];
              jacob_inv[1][1] = jacob[0][0]*jacob[2][2] - jacob[0][2]*jacob[2][0];
              jacob_inv[1][2] = jacob[0][2]*jacob[1][0] - jacob[0][0]*jacob[1][2];

              jacob_inv[2][0] = jacob[1][0]*jacob[2][1] - jacob[1][1]*jacob[2][0];
              jacob_inv[2][1] = jacob[0][1]*jacob[2][0] - jacob[0][0]*jacob[2][1];
              jacob_inv[2][2] = jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0];

              for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                   jacob_inv[p][q] = jacob_inv[p][q] / det_jacob;

           dAxdy = jacob_inv[1][0]*jacob2[0][0] + jacob_inv[1][1]*jacob2[0][1] + jacob_inv[1][2]*jacob2[0][2];
           dAxdz = jacob_inv[2][0]*jacob2[0][0] + jacob_inv[2][1]*jacob2[0][1] + jacob_inv[2][2]*jacob2[0][2];

           dAydx = jacob_inv[0][0]*jacob2[1][0] + jacob_inv[0][1]*jacob2[1][1] + jacob_inv[0][2]*jacob2[1][2];
           dAydz = jacob_inv[2][0]*jacob2[1][0] + jacob_inv[2][1]*jacob2[1][1] + jacob_inv[2][2]*jacob2[1][2];

           dAzdx = jacob_inv[0][0]*jacob2[2][0] + jacob_inv[0][1]*jacob2[2][1] + jacob_inv[0][2]*jacob2[2][2];
           dAzdy = jacob_inv[1][0]*jacob2[2][0] + jacob_inv[1][1]*jacob2[2][1] + jacob_inv[1][2]*jacob2[2][2];

             if(vpot == 1)
             for(k=0;k<neighb_num;k++)
               {
               if(pos[0][k] == 0) continue;
               ax_diff = P[n].Afield[0] - Afield[0][k];
               ay_diff = P[n].Afield[1] - Afield[1][k];
               az_diff = P[n].Afield[2] - Afield[2][k];
               P[n].Bfieldx = P[n].Bfieldx + (1./P[n].Density/P[n].fcor)*P[1].Mass*curl(0, ax_diff, ay_diff, az_diff, dwdu[0][k], dwdu[1][k], dwdu[2][k]);
               P[n].Bfieldy = P[n].Bfieldy + (1./P[n].Density/P[n].fcor)*P[1].Mass*curl(1, ax_diff, ay_diff, az_diff, dwdu[0][k], dwdu[1][k], dwdu[2][k]);
               P[n].Bfieldz = P[n].Bfieldz + (1./P[n].Density/P[n].fcor)*P[1].Mass*curl(2, ax_diff, ay_diff, az_diff, dwdu[0][k], dwdu[1][k], dwdu[2][k]);
               if(n == 1) printf("n = %d x = %lg, y = %lg, z = %lg, bfieldx = %lg, bfieldy = %lg, bfieldz = %lg  ax = %lg, ay = %lg az = %lg, fcor = %lg, dwdu = %lg\n", n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, Afield[0][k], Afield[1][k], Afield[2][k], P[n].fcor, dwdu[0][k]);
               }

             if(vpot == 1)
               {
               P[n].Bfieldx = dAzdy - dAydz;
               P[n].Bfieldy = dAxdz - dAzdx;
               P[n].Bfieldz = dAydx - dAxdy;
               }

             nh_typ = nh_typ + P[n].nh;
             bfield = pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);
             ufield = pow(bfield,2) / (8.*3.14159);
             b_typ  = b_typ + bfield;
             u_typ = u_typ + ufield; 

             bx_typ = bx_typ + P[n].Bfieldx;
             by_typ = by_typ + P[n].Bfieldy;
             bz_typ = bz_typ + P[n].Bfieldz;
             ax_typ = ax_typ + P[n].Afield[0];
             ay_typ = ay_typ + P[n].Afield[1];
             az_typ = az_typ + P[n].Afield[2];

             num_typ = num_typ + 1.0;

             if(n % 100 == 0)
               printf("n = %d x = %lg, y = %lg, z = %lg, bfieldx = %lg, bfieldy = %lg, bfieldz = %lg  ax = %lg, ay = %lg az = %lg, fcor = %lg\n", n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].Afield[0], P[n].Afield[1], P[n].Afield[2], P[n].fcor);
            } //end if(P[n].Id == plist[k].Id) and Ngas for loop

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

  double kfac, dens_conv, y_orion = -1*DeltaX/2, z_orion = -1*DeltaX/2;
  dens_conv = 1.e10 * 1.989e33 / hubble_param * pow(3.08567758e18,-3) / 1.2195 /  1.6726e-24; 

  for(n=0;n<Ngas;n++)
     {
     dis = pow(P[n].disx*P[n].disx + P[n].disy*P[n].disy + P[n].disz*P[n].disz ,0.5);

///////ray-through-box option
/*
     if(P[n].nh > nmin)
     if((P[n].disy-P[n].hsm_pc < y_orion && P[n].disy+P[n].hsm_pc > y_orion 
         && P[n].disz-P[n].hsm_pc < z_orion && P[n].disz+P[n].hsm_pc > z_orion) 
         || (P[n].hsm_pc < DeltaX && fabs(P[n].disy) < DeltaX && fabs(P[n].disz) < DeltaX))
     for(i=0;i<ref_lev;i++)
       {
       if((P[n].disx-P[n].hsm_pc < dis_arr[i] && P[n].disx+P[n].hsm_pc > dis_arr[i]) 
          || (P[n].hsm_pc < DeltaX && fabs(P[n].disx) < dis_arr[i]))
          {
          kfac = calc_kernel_spline(n, P[n].disx, P[n].disy, P[n].disz, dis_arr[i],  y_orion, z_orion, P[n].hsm_pc, DeltaX/2);
*/
////////////////////////////

///////////radial profile option 
     for(i=0;i<ref_lev-1;i++)
       {
       if(dis-P[n].hsm_pc < dis_arr[i] && dis+P[n].hsm_pc > dis_arr[i])
       //if(dis < dis_arr[i+1] && dis > dis_arr[i])
          {
          kfac    = pow(P[n].Mass,1);
          //kfac = calc_kernel_spline(n, P[n].disx, P[n].disy, P[n].disz, P[n].disx + fabs(dis-dis_arr[i]), P[n].disy, P[n].disz, P[n].hsm_pc, DeltaX/2);
//////////////////////////

          num[i]  = num[i] + kfac;
          if(take_profile!=1 ) rho[i]  = rho[i] + P[n].Mass*kfac*dens_conv;
          if(take_profile == 1) rho[i] = rho[i] + P[n].nh*kfac;

          Bx[i]   = Bx[i] + P[n].Bfieldx*kfac;
          By[i]   = By[i] + P[n].Bfieldy*kfac;
          Bz[i]   = Bz[i] + P[n].Bfieldz*kfac;  
          Bmag[i] = Bmag[i] + pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5)*kfac;

          velx[i] = velx[i] + P[n].Vel[0]*pow(Time, 0.5)*kfac;
          vely[i] = vely[i] + P[n].Vel[1]*pow(Time, 0.5)*kfac;
          velz[i] = velz[i] + P[n].Vel[2]*pow(Time, 0.5)*kfac;
          temp[i] = temp[i] + P[n].Temp*kfac;
          }
       }
     }   //end Ngas for loop

  outfile2=fopen(output_fname2, "a");
  for(i=0;i<ref_lev;i++)
    {
    if(take_profile == 1) rho[i] = rho[i] / num[i];

    Bx[i] = Bx[i] / num[i];
    By[i] = By[i] / num[i];
    Bz[i] = Bz[i] / num[i];
    Bmag[i] = Bmag[i] / num[i];
 
    velx[i] = velx[i] / num[i];
    vely[i] = vely[i] / num[i];
    velz[i] = velz[i] / num[i]; 
    temp[i] = temp[i] / num[i];

    fprintf(outfile2, "%g %g %g %g %g %g %g %g %g %g\n", dis_arr[i], rho[i], Bx[i], By[i], Bz[i], Bmag[i],
            velx[i], vely[i], velz[i], temp[i]);
    printf("i = %d, num = %lg, rho = %lg, Bx = %lg, By = %lg, Bz = %lg, B = %lg, vx = %lg, vy = %lg, vz = %lg, temp = %lg\n", 
           i, num[i], rho[i], Bx[i], By[i], Bz[i], Bmag[i], velx[i], vely[i], velz[i], temp[i]);
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

              x = P[n].disx; y = P[n].disy, z = P[n].disz;

    random = rand();
    int to_print = 0;
    //if(random*(1.e0/RAND_MAX)< 0.1e0) to_print++;
    if(fabs(P[n].Bfieldx) > 0 || fabs(P[n].Bfieldy) > 0 || fabs(P[n].Bfieldz) > 0) to_print++;


    if(P[n].to_print > 0)
      {
      if(vpot == 0)
        fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
        P[n].error, P[n].Id,x,y,z,P[n].Temp,P[n].nh,P[n].HII,P[n].H2I,P[n].HDI, P[n].gam, disAU, vrad, vrot, P[n].hsm*Time/0.7, P[n].Mass, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].nh_test);
      if(vpot == 1)
        fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
        P[n].error, P[n].Id,x,y,z,P[n].Temp,P[n].nh,P[n].Afield[0],P[n].Afield[1],P[n].Afield[2], P[n].gam, disAU, vrad, vrot, P[n].hsm*Time/0.7, P[n].Mass, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].nh_test);
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
	  gamma = 1.1;
          //gamma=P[i].gam;	 

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

/*
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
*/

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

