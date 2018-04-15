#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define hubble_param 0.7
#define HUBBLE 3.2407789e-18
#define MAXREF 20
#define def_fac 1

#define which_sim 0

#define width 50.0
#define width_small 1.0  //full width in pc of cubic box of particles whose evolution will be tracked

#define restart 0
#define write_freq 100

#define snapbegin 40  //(which_sim==0)
//#define snapbegin 1003  //ref (which_sim==1)
//#define snapbegin 1003  //ref2 (which_sim==2)?  //note - ref2 has same mass resolution as ref3!
//#define snapbegin 32  //ref3  (which_sim==3)
//#define snapbegin 6600  //ref3, nhmax = 1.e6
//#define snapbegin 7100
//#define snapbegin 0

#define snapend 5980  //(which_sim==0)
//#define snapend 8227   //ref  (which_sim==1)
//#define snapend 15781  //ref2 (which_sim==2)? //note - ref2 has same mass resolution as ref3!
//#define snapend 7133  //zoom10_ref3 (which_sim==3)
//#define snapend 680  //(which_sim==4)
//#define snapend 283  //which_sim==5

#define snapcheck snapend

#define snapnum (snapend - snapbegin)
#define neighbnum 250 

int load_snapshot(char *fname, int files, int myrank);
int reordering(void);
int unit_conversion(void);
int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1);
int allocate_memory(int myrank);
double calc_det(double matrix11, double matrix12, double matrix13, double matrix21, double matrix22, double matrix23, 
                double matrix31, double matrix32, double matrix33);
void reassign_P(int numpart, int arraynum, int n1, int n2);
void reassign_plist(void);
double calc_kernel_spline(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double cosmo_fac);

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



int     NumPart, Ngas, Ngas0;

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];

  double  Mass;
  int    Type;
  int Id;

  double  Rho, Rho_phys, Temp, nh;
  double hsm; 
  double gam; 
  double dummy;
  int marker;
} *P;


struct particle_data0
{
  double  Pos[3];
  int Id;
  double Vel[3];
  double  Rho;
} *P0;

struct plist_data0
{
  double vel[3], veltot; 
  double pos[3];
  double hsm; 
  double dens_conv, dens_alt;
  double dis_avg_all; 
  double NumCalc;
  double dens;
  double B_anal[3], nh_anal, nh_anal_alt;
} *plist0;

struct plist_data1
{
  double vel[3], veltot;
  double pos[3];
  double hsm;
  int    NumNeighb;
  double detD;
  double dens_conv, dens_alt;
  //double dis_avg
  double dis_avg_all; 
  //double delx_avg, dxdx_avg; 
  //double delnh_avg;
  //double vel_avg0, vel_avg1, vel_disp0, vel_disp1;
  double NumCalc;
  double dens;
  double temp;
  double B_anal[3], nh_anal, nh_anal_alt;
} *plist1;


struct neighb_data
{
  double vel0[3][neighbnum], vel1[3][neighbnum];  
  double pos0[3][neighbnum], pos1[3][neighbnum];
  int id[neighbnum];
  double kernel0[neighbnum], kernel1[neighbnum];
  //double nh1[neighbnum], nh0[neighbnum];
  //double deform[3][3][neighbnum];
} *nb;
 
 

//int *Id;

double  Time, zred;
double Time_new, Time_old, delta_t, t_Hubble, t_Hubble_old;


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char path[200], pathout[200], input_fname[200], input_fname2[200], plist_fname[200], output_fname[200], output_fname2[300]; 
  char basename[200], basename2[200], basenameout[200];
  int  ID, arrnum, i, k=0, ksort, nskip, write_file, j, n, p, q, r, s, type, snapshot_number, files, Ngas, Ngas0, random, counter;
  double x,y,z,x1,y1,z1, delr, neighbnum_doub, nthresh = 1.e2, nthresh_cur, delfac = 1.e-3; 
  double hfac = 1.0, hfac_in = 0.0,  hfac_out = 1.5, hfac_avg; 
  double time_fac, time_fac0, cosmo_fac, cosmo_fac0, init_pow=0., fix_pow=2./3.;
  double nh, nhmax, nhmin_old=0.0, nhmin, mass, mmax, dis_min, dis_sim, disAU, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double vrad, vrot, disAUxy, disAUz, disAUxyz, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2; 
  double num, kernel0, kernel1, kernel_fac0, kernel_fac1; 
  double xfac, yfac, zfac, kfloat, nh_avg[snapnum], nh_avg2, num_tot;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, jzcent;
  double B0_100=1.e-14, B0, B0x, B0y, B0z, rho_baryon, Bfield[snapnum], deform_avg[3][3]; 
  double detD_temp, detD_neighb[neighbnum];
  int    snaparr[snapnum+1], *plist_id, ncount, ncount_tot;
  double vradx, vrady, vradz, vrotx, vroty, vrotz; 
  double velx, vely, velz, vel_conv;
  double rad, dmin, disx, disy, disz, del[3][3], del0[3][3];
  double vel[neighbnum], vel0[neighbnum], dis[neighbnum], dis0[neighbnum], ncount_doub;
  double davgHI0[3][3], davgLO0[3][3], davgHI1[3][3], davgLO1[3][3], davgHI2[3][3], numHI0[3][3], numLO0[3][3], numHI1[3][3], numLO1[3][3], numHI2[3][3];
  double ncheckHI[3], ncheckLO[3], ncheck_max = 100;

  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs, Hubble, hubble_a, hubble_a_old;
  UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
  UnitVelocity_in_cm_per_s= 1.0e5;
  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  Hubble = HUBBLE * UnitTime_in_s;

  FILE *pfile, *infile, *outfile, *outfile2;
 
  int npes, myrank, ierr, tot_proc=1, tot_proc_sum, nmin, nmax, arrnum_per_proc, send_size;
  double *Ptot, *PBfieldx, *PBfieldy, *PBfieldz, *Pnh_anal;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  hfac_avg = 0.0;
 
  sprintf(basename, "bin_zoom10");

  if(which_sim == 0)
    {
    sprintf(path, "/work/00863/minerva/bin_zoom_cut");
    sprintf(basename2, "bin_zoom10_new_cut");
    }
  if(which_sim == 1)
    {
    sprintf(path, "/scratch/00863/minerva");
    sprintf(basename2, "bin_zoom10_new_cut_ref");
    }
  if(which_sim == 2)
    {
    sprintf(path, "/scratch/00863/minerva");
    sprintf(basename2, "bin_zoom10_new_cut_ref2");
    }
  if(which_sim == 3)
    {
    sprintf(path, "/scratch/00863/minerva"); 
    sprintf(basename2, "bin_zoom10_new_cut_ref3");
    }
  if(which_sim == 4)
    {
    sprintf(path, "/work/00863/minerva/");
    sprintf(basename2, "bin_HR10_ideal");
    sprintf(basename, "bin_HR10_ideal");
    }
  if(which_sim == 5)
    {
    sprintf(path, "/work/00863/minerva/");
    sprintf(basename2, "bin_MR10_ideal");
    sprintf(basename, "bin_MR10_ideal");
    }

  sprintf(pathout, "/work/00863/minerva/orion/"); 


  for(n=0;n<=snapnum;n++)
    {
    snaparr[n] = snapbegin+(1*n);
//    if(which_sim == 3 && snapbegin + 1*n == 6636) snaparr[n] = 6635;
//    if(which_sim == 3 && snapbegin + 1*n == 6771) snaparr[n] = 6770;
//    if(which_sim == 3 && snapbegin + 1*n == 6846) snaparr[n] = 6845;
//    if(which_sim == 3 && snapbegin + 1*n == 6924) snaparr[n] = 6923;
    }


  for(n=0;n<snapnum;n++)
    {
    nh_avg[n] = 0;
    if(myrank == 0)
      printf("snaparr = %d \n", snaparr[n]);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////
//// identify particles whose evolution shold be analytically followed
  snapshot_number = snapcheck;
  files = 1;
  sprintf(input_fname, "%s/%s_%04d", path, basename2, snapshot_number);
  if(which_sim == 4 || which_sim == 5)  sprintf(input_fname, "%s/%s_%03d", path, basename2, snapshot_number);
  if(snapshot_number > 9999)
    sprintf(input_fname, "%s/%s_%05d", path, basename2, snapshot_number);  
  Ngas = load_snapshot(input_fname, files, myrank);

  unit_conversion();

  time_fac = Time;
  //For NON-cosmological runs
  //time_fac = 1.0;

  cosmo_fac = time_fac/(hubble_param);

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
  printf("myrank = %d, nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", myrank, nhmax, xmax, ymax, zmax);

  nskip = 1;

  arrnum=0;
  for(i = 0; i < Ngas; i++)
     {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%nskip == 0)
        {
        arrnum++;
        if(myrank == 0 && i % 10000 == 0)
          printf("disx = %lg, disy = %lg, disz = %lg\n", disx, disy, disz);
        }
     }

  plist_id =  (int *) malloc(arrnum * sizeof(int));

  n=0;
  for(i = 0; i < Ngas; i++)
      {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
     disy = fabs((P[i].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
     disz = fabs((P[i].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
     if(disx < width_small/2.0 && disy < width_small/2.0 && disz < width_small/2.0 && i%nskip == 0)
        {
        plist_id[n] = P[i].Id;
        n++;
        }
      }

  printf("arrnum = %d, Time = %lg, hubble_param = %lg\n", arrnum, Time, hubble_param);
  MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  arrnum_per_proc = ((int) (arrnum/tot_proc_sum)) + 1;
  nmin = arrnum_per_proc * (myrank);
  nmax = arrnum_per_proc * (myrank+1) -1;
  if(myrank == tot_proc_sum -1)
    nmax = arrnum-1;

  if(nmin >= arrnum-1)
    nmin = nmax = arrnum-1;
  if(nmax >= arrnum-1)
    nmax = arrnum-1;

  printf("myrank = %d, arrnum_per_proc = %d, nmin = %d, nmax = %d\n", myrank, arrnum_per_proc, nmin, nmax);
  free(P);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  int n_snapskip = 1, nbad, nbad2, nbad3, nbad_max;
  int Ngas_sum;

  for(j=0;j<=snapnum-1;j++)
  {
  //MPI_Barrier(MPI_COMM_WORLD);

  snapshot_number= snaparr[j];                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */
  Ngas_sum = 0;

  if(j > 0 && snapshot_number%n_snapskip != 0) continue;

  sprintf(input_fname, "%s/%s_%04d", path, basename2, snapshot_number);
  if(which_sim == 4 || which_sim == 5)  sprintf(input_fname, "%s/%s_%03d", path, basename2, snapshot_number);
  if(snapshot_number > 9999)
    sprintf(input_fname, "%s/%s_%05d", path, basename2, snapshot_number);

  cosmo_fac0 = cosmo_fac;
  Time_old = Time_new;
  t_Hubble_old = t_Hubble;
  hubble_a_old = hubble_a;

  if(which_sim == 3 && (snapshot_number == 5123 || snapshot_number == 5363 || snapshot_number == 6489))
    {
    printf("myank = %d, SKIP this snapshot!\n", myrank);
    continue;
    }

  Ngas = load_snapshot(input_fname, files, myrank);

  MPI_Allreduce(&Ngas, &Ngas_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(Ngas_sum < Ngas * (tot_proc_sum-5)) 
    {
    printf("myank = %d, Ngas_sum = %d, SKIP this snapshot!\n", myrank, Ngas_sum);
    continue;
    }
  if(j==1) printf("tot_proc_sum = %d, Ngas_sum = %d\n", tot_proc_sum, Ngas_sum);

  sprintf(output_fname, "%s_plist.dat",basename);

  hubble_a = header1.Omega0 / (header1.time *header1.time * header1.time)
        + (1 - header1.Omega0 - header1.OmegaLambda) / (header1.time * header1.time) + header1.OmegaLambda;
  hubble_a = Hubble * hubble_param * sqrt(hubble_a);
  Time_new = pow(header1.time, 1.5) /1.5 / HUBBLE / hubble_param /sqrt(header1.Omega0);
  delta_t = Time_new - Time_old;

  t_Hubble = 5.4e8/pow((1.e0+header1.redshift)/10.e0, 1.5e0) * 3.14e7;
  if(myrank == 0) printf("Time = %lg, Time_old = %lg, Time_new = %lg, t_Hub_est = %lg, hubble_a = %lg, delta_t = %lg\n", Time, Time_old, Time_new, t_Hubble - t_Hubble_old, hubble_a, delta_t);

  unit_conversion();

  time_fac = Time;
  //For NON-cosmological runs
  //time_fac = 1.0;
  
   cosmo_fac = time_fac/(hubble_param);
   vel_conv = pow(Time, 0.5) * 3.24077929e-17; //convert from sim. units to pc / s
   B0 = B0_100 * pow((1.+zred) / 100., 2);
   rho_baryon = 3.8e-31*pow(1.+zred, 3);  //average baryonic density of universe
   if(myrank == 0)
     printf("B0 = %lg, zred = %lg rho_baryon=%lg\n", B0, zred, rho_baryon);

///////////////////////////////////////////////////////////////////////////////////////////////////////
  if(restart == 1 && j == 0)
    {
    printf("begin file read!\n");
    PBfieldx = (double *) malloc(Ngas * sizeof(double));
    PBfieldy = (double *) malloc(Ngas * sizeof(double));
    PBfieldz = (double *) malloc(Ngas * sizeof(double));
    Pnh_anal = (double *) malloc(Ngas * sizeof(double));
    sprintf(input_fname2, "%s/%s_bfield_%04d", pathout, basename, snapbegin);
    if(which_sim == 1)
      {
      sprintf(input_fname2, "%s/%s_bfield_ref_%04d", pathout, basename, snapbegin);
      if(snapbegin > 9999) sprintf(input_fname2, "%s/%s_bfield_ref_%05d", pathout, basename, snapbegin);
      }
    if(which_sim == 2)
      {
      sprintf(input_fname2, "%s/%s_bfield_ref2_%04d", pathout, basename, snapbegin);
      if(snapbegin > 9999) sprintf(input_fname2, "%s/%s_bfield_ref2_%05d", pathout, basename, snapbegin);
      }
    if(which_sim == 3)
      sprintf(input_fname2, "%s/%s_bfield_ref3_%04d", pathout, basename, snapbegin);
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }
    for(n=0;n<Ngas;n++)
        fread(&PBfieldx[n], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&PBfieldy[n], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&PBfieldz[n], sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&Pnh_anal[n], sizeof(double), 1, infile);
    fclose(infile);
    printf("end file read!\n");
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////

  write_file = 0;
  if(j+n_snapskip == snapnum ||(snapshot_number%write_freq==0 && j >= 1) || j == 1)
    {
    write_file = 1;
    Ptot =  (double *) malloc(Ngas * sizeof(double));
    PBfieldx = (double *) malloc(Ngas * sizeof(double));
    PBfieldy = (double *) malloc(Ngas * sizeof(double));
    PBfieldz = (double *) malloc(Ngas * sizeof(double));
    Pnh_anal = (double *) malloc(Ngas * sizeof(double));
    }

  if(j > 0) nhmin_old = nhmin;
  nhmin = 1.e50;
  for(i = 0; i < Ngas; i++)
     {
     P[i].marker = 0;
     if(P[i].nh > 0.6*nhmin_old)
     for(n=0; n<arrnum; n++)
       if(P[i].Id == plist_id[n])
         {
         P[i].marker = 1;
         if(P[i].nh < nhmin) nhmin = P[i].nh;
         }
     }
     

  num_tot=ncount=0;
  ncount_doub = h2max = nhmax = mmax = kfloat = 0.0;
  xtimesm = ytimesm = ztimesm = 0.0;
  vxtimesm = vytimesm = vztimesm = 0.0;
  vxCOM=vyCOM=vzCOM=xCOM=yCOM=zCOM=nthresh_cur=0;

  if(!(plist1=(struct plist_data1 *) malloc(arrnum*sizeof(struct plist_data1))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }
  if(!(nb=(struct neighb_data *) malloc((nmax-nmin+1)*sizeof(struct neighb_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }

  for(n=0; n<arrnum; n++)
      {
      plist1[n].NumNeighb = 1;
      plist1[n].NumCalc = 0.0;
      plist1[n].dis_avg_all = 0.;
      plist1[n].dens_alt = 0;
      plist1[n].dens = plist1[n].dens_conv = plist1[n].hsm = 0;
      plist1[n].pos[0] = plist1[n].pos[1] = plist1[n].pos[2] = plist1[n].vel[0] = plist1[n].vel[1] = plist1[n].vel[2] = 0;
      plist1[n].temp = plist1[n].veltot = plist1[n].nh_anal = 0;
      }

   if(j == 0)
     for(n=0; n<arrnum; n++)
          for(p=0; p<3; p++)  plist1[n].B_anal[p] = 0;

       if(write_file == 1)
         for(i = 0; i < Ngas; i++)
           {     
           Ptot[i] = PBfieldx[i] = PBfieldy[i] = PBfieldz[i] = Pnh_anal[i] = 0;
           }	

       for(i = 0; i < Ngas; i++)
        {
        if(P[i].nh > nhmax /*&& P[i].sink < 0.5*/)
           {
           nhmax = P[i].nh;
           xmax = P[i].Pos[0];
           ymax = P[i].Pos[1];
           zmax = P[i].Pos[2];
           mmax = P[i].Mass;
           vxCOM = P[i].Vel[0];
           vyCOM = P[i].Vel[1];
           vzCOM = P[i].Vel[2];
           }
         }
   
        //nthresh_cur = 1.e-1;
        nthresh_cur = 0.4*nhmin; 

        if(myrank == 0)
          printf("nhmax = %lg, nthresh = %lg, nthresh_cur = %lg, vxCOM = %lg, vyCOM = %lg, vzCOM = %lg\n", nhmax, nthresh, nthresh_cur, vxCOM, vyCOM, vzCOM);
   
	sinkposx1 = xmax;
	sinkposy1 = ymax;
	sinkposz1 = zmax;


	ncount=ncount_tot=counter=0;  
	ncount_doub=0.;

        //printf("line 313\n");       
 
        if(myrank == 0)
          outfile=fopen(output_fname, "a");

         B0x = 0.;
         B0y = 0.;
         B0z = 1.e-20;

          for(i = 0; i < Ngas; i++)
             {
 
             P[i].Pos[0] = P[i].Pos[0]*cosmo_fac;
             P[i].Pos[1] = P[i].Pos[1]*cosmo_fac;
             P[i].Pos[2] = P[i].Pos[2]*cosmo_fac;
             P[i].hsm = P[i].hsm*cosmo_fac*hfac_out;     

             if(write_file == 1)
                { 
                 if(P[i].marker > 0)
                   {
                   Pnh_anal[i] = 0;
                   PBfieldx[i] = 0;
                   PBfieldy[i] = 0;
                   PBfieldz[i] = 0;
                   }
                 else
                   {
                   Pnh_anal[i] = P[i].Rho / (double) tot_proc_sum;
                   PBfieldx[i] = B0x * pow(P[i].Rho_phys/rho_baryon, init_pow) / (double) tot_proc_sum;
                   PBfieldy[i] = B0y * pow(P[i].Rho_phys/rho_baryon, init_pow) / (double) tot_proc_sum;
                   PBfieldz[i] = B0z * pow(P[i].Rho_phys/rho_baryon, init_pow) / (double) tot_proc_sum;
                   }
                 }

             if(P[i].marker > 0)
             //if(P[i].nh > nthresh_cur) 
               for(n=nmin; n<=nmax; n++)
                      {
                      if(P[i].Id == plist_id[n])
                         {
                         plist1[n].dens = P[i].Rho_phys;
                         plist1[n].dens_conv = P[i].nh/P[i].Rho_phys;
                         plist1[n].hsm = P[i].hsm;
                         plist1[n].pos[0] = P[i].Pos[0];
                         plist1[n].pos[1] = P[i].Pos[1];
                         plist1[n].pos[2] = P[i].Pos[2];
                         plist1[n].vel[0] = P[i].Vel[0];
                         plist1[n].vel[1] = P[i].Vel[1];
                         plist1[n].vel[2] = P[i].Vel[2];
                         plist1[n].temp = P[i].Temp;
                         plist1[n].veltot = pow(P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2], 0.5);
                         if(n < 5)
                           printf("myrank = %d, j = %d, n = %d, nh_anal1 = %lg\n", myrank, j, n, plist1[n].nh_anal);
                         if(plist1[n].B_anal[0] == 0 && plist1[n].B_anal[1] == 0 && plist1[n].B_anal[2] == 0)
                           {
                           plist1[n].nh_anal = P[i].Rho;
                           plist1[n].nh_anal_alt = P[i].Rho;
                           plist1[n].B_anal[0] = B0x * pow(P[i].Rho_phys/rho_baryon, init_pow);
                           plist1[n].B_anal[1] = B0y * pow(P[i].Rho_phys/rho_baryon, init_pow);
                           plist1[n].B_anal[2] = B0z * pow(P[i].Rho_phys/rho_baryon, init_pow);
                           }
                          if(j == 0 && restart == 1)
                           {
                           plist1[n].nh_anal = Pnh_anal[i];
                           plist1[n].nh_anal_alt = Pnh_anal[i];
                           plist1[n].B_anal[0] = PBfieldx[i];
                           plist1[n].B_anal[1] = PBfieldy[i];
                           plist1[n].B_anal[2] = PBfieldz[i];
                           }

                         if(write_file == 1)
                           {
                           Pnh_anal[i] = plist0[n].nh_anal;
                           PBfieldx[i] = plist0[n].B_anal[0];
                           PBfieldy[i] = plist0[n].B_anal[1];
                           PBfieldz[i] = plist0[n].B_anal[2];
                           //on processor that was actually following *_anal evolution, reset the P* values to match the calculated *_anal values
                           if(counter < 2)
                             {
                             printf("RESET PNH_ANAL Pnh_anal = %lg, PBfieldx = %lg, nh_anal = %lg, Bfieldx = %lg\n", 
                                   Pnh_anal[i], PBfieldx[i], plist0[n].nh_anal, plist0[n].B_anal[0]);
                             counter++;
                             }
                           }
                         ncount++;
                         }
                      }
              }

        //printf("line 383, ncount_tot = %d\n", ncount_tot);

        if(j==0)
          {
          Ngas0 = Ngas;
          reassign_P(Ngas0, arrnum, nmin, nmax);
          }

        //printf("line 444, ncount_tot = %d\n", ncount_tot);

        if(j > 0)
         {
          for(i = 0; i < Ngas; i++)
            {
             dis_min = 100.;
             if(P[i].nh > nthresh_cur)
             for(n=nmin; n<=nmax; n++)
                {
                 rad = pow(((P[i].Pos[0]-plist1[n].pos[0])*(P[i].Pos[0]-plist1[n].pos[0]) + (P[i].Pos[1]-plist1[n].pos[1])*(P[i].Pos[1]-plist1[n].pos[1]) + (P[i].Pos[2]-plist1[n].pos[2])*(P[i].Pos[2]-plist1[n].pos[2])), 0.5);
                 disx =  fabs(P[i].Pos[0] - plist1[n].pos[0]);
                 disy =  fabs(P[i].Pos[1] - plist1[n].pos[1]);
                 disz =  fabs(P[i].Pos[2] - plist1[n].pos[2]);
              
                 if(rad < hfac*plist1[n].hsm && rad > hfac_in*plist1[n].hsm && plist1[n].NumNeighb < neighbnum-1)
                 //if(disx < hfac*plist1[n].hsm && disy < hfac*plist1[n].hsm && disz < hfac*plist1[n].hsm && rad > hfac_in*plist1[n].hsm && plist1[n].NumNeighb < neighbnum-1)
                     {
                     dis_min = rad;
                     k = plist1[n].NumNeighb -1;
                     plist1[n].dens_alt = plist1[n].dens_alt + P[i].Mass;
                     nb[n-nmin].kernel1[k] = calc_kernel_spline(i, plist1[n].pos[0], plist1[n].pos[1], plist1[n].pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], plist1[n].hsm, cosmo_fac);
                     kernel1 = 1;
                     //kernel1 = nb[n].kernel1[k];
                     nb[n-nmin].pos1[0][k] = P[i].Pos[0] * pow(kernel1,1./3.);
                     nb[n-nmin].pos1[1][k] = P[i].Pos[1] * pow(kernel1,1./3.);
                     nb[n-nmin].pos1[2][k] = P[i].Pos[2] * pow(kernel1,1./3.);
                     nb[n-nmin].vel1[0][k] = P[i].Vel[0];
                     nb[n-nmin].vel1[1][k] = P[i].Vel[1];
                     nb[n-nmin].vel1[2][k] = P[i].Vel[2];
                     nb[n-nmin].id[k] = P[i].Id;
                     plist1[n].NumNeighb++;
                     }
                }
             }

          for(i = 0; i < Ngas0; i++)
             {
             if(P0[i].Rho > 0.5*1.67e-24*nthresh_cur)
             for(n=nmin; n<=nmax; n++)
               {
               rad = pow(((P0[i].Pos[0]-plist0[n].pos[0])*(P0[i].Pos[0]-plist0[n].pos[0]) + (P0[i].Pos[1]-plist0[n].pos[1])*(P0[i].Pos[1]-plist0[n].pos[1]) + (P0[i].Pos[2]-plist0[n].pos[2])*(P0[i].Pos[2]-plist0[n].pos[2])), 0.5);
               if(rad < 1.5*hfac*plist1[n].hsm)
               for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  if(P0[i].Id == nb[n-nmin].id[k])
                     {
                     nb[n-nmin].kernel0[k] = calc_kernel_spline(i, plist0[n].pos[0], plist0[n].pos[1], plist0[n].pos[2], P0[i].Pos[0], P0[i].Pos[1], P0[i].Pos[2], plist0[n].hsm, cosmo_fac);                          
                     kernel0 = 1;
                     //kernel0 = nb[n].kernel0[k];
                     nb[n-nmin].pos0[0][k] = P0[i].Pos[0] * pow(kernel0,1./3.);
                     nb[n-nmin].pos0[1][k] = P0[i].Pos[1] * pow(kernel0,1./3.);
                     nb[n-nmin].pos0[2][k] = P0[i].Pos[2] * pow(kernel0,1./3.);
                     nb[n-nmin].vel0[0][k] = P0[i].Vel[0];
                     nb[n-nmin].vel0[1][k] = P0[i].Vel[1];
                     nb[n-nmin].vel0[2][k] = P0[i].Vel[2];

/////////////////////////////////////////////////////////
//particles shifted by a factor 'delfac'
/*  
                     nb[n-nmin].pos1[0][k] = plist1[n].pos[0] + (1. - delfac)*(nb[n-nmin].pos0[0][k] - plist0[n].pos[0]);
                     nb[n-nmin].pos1[1][k] = plist1[n].pos[1] + (1. - delfac)*(nb[n-nmin].pos0[1][k] - plist0[n].pos[1]);
                     nb[n-nmin].pos1[2][k] = plist1[n].pos[2] + (1. - delfac)*(nb[n-nmin].pos0[2][k] - plist0[n].pos[2]);
*/ 
//particles shifted solely by cosmological expansion factor
/*
                     nb[n-nmin].pos1[0][k] = plist1[n].pos[0] + (cosmo_fac/cosmo_fac0)*(nb[n-nmin].pos0[0][k] - plist0[n].pos[0]);
                     nb[n-nmin].pos1[1][k] = plist1[n].pos[1] + (cosmo_fac/cosmo_fac0)*(nb[n-nmin].pos0[1][k] - plist0[n].pos[1]);
                     nb[n-nmin].pos1[2][k] = plist1[n].pos[2] + (cosmo_fac/cosmo_fac0)*(nb[n-nmin].pos0[2][k] - plist0[n].pos[2]);
*/
///////////////////////////////////////////////////////////

                     }
                   }
                 }
              }

            nbad = nbad2 = nbad3 = 0;
            nbad_max = 3;
            for(n=nmin; n<=nmax; n++)
                {
                plist1[n].detD = 0.0;
                for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                  davgHI0[p][q]=davgLO0[p][q]=davgHI1[p][q]=davgLO1[p][q]=davgHI2[p][q]=numHI0[p][q]=numLO0[p][q]=numHI1[p][q]=numLO1[p][q]=numHI2[p][q]=0;

                dmin = 0.05*plist1[n].hsm;

                for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                    deform_avg[p][q] = 0;

                for(k=0; k<neighbnum; k++)
                  detD_neighb[k] = 0.;

                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {                
                  dis0[k] = pow(((nb[n-nmin].pos0[0][k]-plist0[n].pos[0])*(nb[n-nmin].pos0[0][k]-plist0[n].pos[0])
                       + (nb[n-nmin].pos0[1][k]-plist0[n].pos[1])*(nb[n-nmin].pos0[1][k]-plist0[n].pos[1])
                       + (nb[n-nmin].pos0[2][k]-plist0[n].pos[2])*(nb[n-nmin].pos0[2][k]-plist0[n].pos[2])), 0.5);

                  dis[k] = pow(((nb[n-nmin].pos1[0][k]-plist1[n].pos[0])*(nb[n-nmin].pos1[0][k]-plist1[n].pos[0])
                       + (nb[n-nmin].pos1[1][k]-plist1[n].pos[1])*(nb[n-nmin].pos1[1][k]-plist1[n].pos[1])
                       + (nb[n-nmin].pos1[2][k]-plist1[n].pos[2])*(nb[n-nmin].pos1[2][k]-plist1[n].pos[2])), 0.5);

                  vel0[k] = pow(nb[n-nmin].vel0[0][k]*nb[n].vel0[0][k] + nb[n-nmin].vel0[1][k]*nb[n-nmin].vel0[1][k]
                       + nb[n-nmin].vel0[2][k]*nb[n-nmin].vel0[2][k], 0.5);

                  vel[k] = pow(nb[n-nmin].vel1[0][k]*nb[n-nmin].vel1[0][k] + nb[n-nmin].vel1[1][k]*nb[n-nmin].vel1[1][k]
                       + nb[n-nmin].vel1[2][k]*nb[n-nmin].vel1[2][k], 0.5);

                  plist1[n].dis_avg_all = plist1[n].dis_avg_all + dis[k];
                  }
 
                plist1[n].dis_avg_all = plist1[n].dis_avg_all/((double)plist1[n].NumNeighb-1);


                for(q=0; q<3; q++) ncheckLO[q] = 0;
                for(q=0; q<3; q++) ncheckHI[q] = 0;
                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  if(dis0[k] < hfac*plist0[n].hsm)
                     {
                     for(q=0; q<3; q++)
                       if((plist0[n].pos[q] - nb[n-nmin].pos0[q][k]) < hfac_avg*plist1[n].hsm)
                         ncheckLO[q]++;           
                     for(q=0; q<3; q++)
                       if((plist0[n].pos[q] - nb[n-nmin].pos0[q][k]) > hfac_avg*plist1[n].hsm)
                         ncheckHI[q]++;
                      }
                   }
                 for(q=0; q<3; q++)
                   { 
                   if(ncheckLO[q] <= ncheckHI[q]) ncheckHI[q] = ncheckLO[q];
                   if(ncheckHI[q] <= ncheckLO[q]) ncheckLO[q] = ncheckHI[q];
                   if(ncheckLO[q] > ncheck_max) ncheckLO[q] = ncheck_max;
                   if(ncheckHI[q] > ncheck_max) ncheckHI[q] = ncheck_max;
                   if(which_sim == 0) ncheckHI[q] = ncheckLO[q] = ncheck_max;
                   }

                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  kernel_fac0 = nb[n-nmin].kernel0[k];
                  kernel_fac1 = nb[n-nmin].kernel1[k];
                  kernel_fac0 = kernel_fac1 = 1;
                  if(dis0[k] < hfac*plist0[n].hsm)
                     {
                     for(q=0; q<3; q++)
                       if((plist0[n].pos[q] - nb[n-nmin].pos0[q][k]) < hfac_avg*plist1[n].hsm && numLO0[0][q] < ncheckLO[q])  
                         {
                         for(p=0; p<3; p++) davgLO0[p][q] = davgLO0[p][q] + (plist0[n].pos[p] - nb[n-nmin].pos0[p][k])*kernel_fac0;
                         for(p=0; p<3; p++) numLO0[p][q] = numLO0[p][q] + 1.0*kernel_fac0; 
                         for(p=0; p<3; p++) davgLO1[p][q] = davgLO1[p][q] + (plist1[n].pos[p] - nb[n-nmin].pos1[p][k])*kernel_fac1;
                         for(p=0; p<3; p++) numLO1[p][q] = numLO1[p][q] + 1.0*kernel_fac1;
                         }
                     for(q=0; q<3; q++) 
                       if((plist0[n].pos[q] - nb[n-nmin].pos0[q][k]) > hfac_avg*plist1[n].hsm && numHI0[0][q] < ncheckHI[q])
                         {
                         for(p=0; p<3; p++) davgHI0[p][q] = davgHI0[p][q] + (plist0[n].pos[p] - nb[n-nmin].pos0[p][k])*kernel_fac0;
                         for(p=0; p<3; p++) numHI0[p][q] = numHI0[p][q] + 1.0*kernel_fac0;
                         for(p=0; p<3; p++) davgHI1[p][q] = davgHI1[p][q] + (plist1[n].pos[p] - nb[n-nmin].pos1[p][k])*kernel_fac1;
                         for(p=0; p<3; p++) numHI1[p][q] = numHI1[p][q] + 1.0*kernel_fac1;
                         }
                      }
                   }                

                for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                  {
                  if(numLO0[p][q] <= 0) numLO0[p][q] = 1.0;
                  if(numLO1[p][q] <= 0) numLO1[p][q] = 1.0;
                  if(numHI0[p][q] <= 0) numHI0[p][q] = 1.0;
                  if(numHI1[p][q] <= 0) numHI1[p][q] = 1.0;
                  if(numHI2[p][q] <= 0) numHI2[p][q] = 1.0;
                  }

             
                for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                  {
                  if(n < 2)
                    printf("[%d][%d] davgLO0= %lg, davgHI0=%lg, davgLO1 = %lg davgHI1 = %lg, numLO0 = %lg, numHI0 = %lg, numLO1 = %lg, numHI1 = %lg\n",
                          p, q, davgLO0[p][q], davgHI0[p][q], davgLO1[p][q], davgHI1[p][q], numLO0[p][q], numHI0[p][q], numLO1[p][q], numHI1[p][q]);

                  davgLO0[p][q] = davgLO0[p][q]/numLO0[p][q];

                  davgLO1[p][q] = davgLO1[p][q]/numLO1[p][q];

                  davgHI0[p][q] = davgHI0[p][q]/numHI0[p][q];

                  davgHI1[p][q] = davgHI1[p][q]/numHI1[p][q];
                  }

 
                for(p=0; p<3; p++)
                   for(q=0; q<3; q++)
                     {
                     deform_avg[p][q] =  ((davgHI1[p][q] - davgLO1[p][q]) - (davgHI0[p][q] - davgLO0[p][q])) / (davgHI0[q][q] - davgLO0[q][q]);
                     //deform_avg[p][q] =  ((davgHI1[p][q] + fabs(davgLO1[p][q])) - (davgHI0[p][q] + fabs(davgLO0[p][q]))) / (davgHI0[q][q] + fabs(davgLO0[q][q]));
                     //deform_avg[p][q] = (davgHI1[p][q] - davgLO1[p][q]) / (davgHI0[q][q] - davgLO0[q][q]);
                     //if(n < 1) printf("n = %d, deform[%d][%d] = %lg\n", n, p, q, deform_avg[p][q]);
                     }

                 plist1[n].detD = calc_det(deform_avg[0][0], deform_avg[0][1], deform_avg[0][2], deform_avg[1][0],
                                    deform_avg[1][1], deform_avg[1][2], deform_avg[2][0], deform_avg[2][1], deform_avg[2][2]);


                 int neighb_min = 10;                 
                 if(which_sim == 0) neighb_min = 8;

                 for(p=0; p<3; p++)
                   {
                 for(q=0; q<3; q++)
                   {
                   if(numHI1[p][q] < neighb_min || numHI0[p][q] < neighb_min)
                     {
                     plist1[n].detD = 1;
                     deform_avg[p][q] = 0;
                     if(p == 0 && q == 0)
                        { 
                        if(nbad < nbad_max) printf("FIXING IT, n = %d, dens = %lg, nh_anal = %lg, num_neighb = %d, hsm0 = %lg, hsm1 = %lg\n", n, plist1[n].dens, plist0[n].nh_anal, plist1[n].NumNeighb, plist0[n].hsm, plist1[n].hsm);
                        //printf("x0 = %15.10g, y0 = %15.10g, z0 = %15.10g \n", plist0[n].pos[0], plist0[n].pos[1], plist0[n].pos[2]);
                        //for(k=0; k<plist1[n].NumNeighb-1; k++)
                        //   printf("x_nb = %15.10g, y_nb = %15.10g, z_nb = %15.10g \n", nb[n-nmin].pos0[0][k], nb[n-nmin].pos0[1][k], nb[n-nmin].pos0[2][k]);       
                        nbad++;
                        }
                     if(p == 0 && nbad < nbad_max) printf("ncheckLO[%d]=%lg, ncheckHI[%d]=%lg\n", q, ncheckLO[q], q, ncheckHI[q]);
                     printf("numLO1[%d][%d]=%lg, numLO0[%d][%d]=%lg, numHI[%d][%d]=%lg, numHI0[%d][%d]=%lg\n", 
                     p, q, numHI1[p][q], p, q, numHI1[p][q], p, q, numLO0[p][q], p, q, numLO0[p][q]);
       
                     }
                   if((numHI1[p][q] < neighb_min || numHI0[p][q] < neighb_min) && plist0[n].dens > 0 && plist1[n].dens > 0)
                     {
                     plist1[n].detD = plist0[n].dens / plist1[n].dens;
                     if(p == 0 && q == 0 && nbad2 < nbad_max) printf("FIXING IT2, n = %d, dens = %lg, nh_anal = %lg\n", n, plist1[n].dens, plist0[n].nh_anal);
                     nbad2++;
                     }
                   }
                   }

                 if(plist0[n].nh_anal <= 0 || plist0[n].nh_anal != plist0[n].nh_anal)
                   {
                   plist0[n].nh_anal = plist0[n].dens;
                   plist0[n].B_anal[0] = B0x * pow(plist0[n].dens/rho_baryon, fix_pow);
                   plist0[n].B_anal[1] = B0y * pow(plist0[n].dens/rho_baryon, fix_pow);
                   plist0[n].B_anal[2] = B0z * pow(plist0[n].dens/rho_baryon, fix_pow);
                   if(nbad3 < nbad_max) printf("FIXING IT3, n = %d, dens = %lg, plist0[n].B_anal[2] = %lg\n", n, plist1[n].dens, plist0[n].B_anal[2]);
                   nbad3++;
                   }

                 plist1[n].nh_anal  = plist0[n].nh_anal / plist1[n].detD;
                 if(plist0[n].dis_avg_all > 0)
                   plist1[n].nh_anal_alt = plist0[n].nh_anal_alt * pow(plist0[n].dis_avg_all/plist1[n].dis_avg_all,3);
                 plist1[n].dens_alt = plist1[n].dens_alt / (4.*3.14159 /3) / pow(hfac*plist1[n].hsm, 3); 
                 plist1[n].dens_alt = plist1[n].dens_alt*plist1[n].dens_conv;

                 //for(p=0; p<3; p++)
                 //  plist1[n].B_anal[p] = plist0[n].B_anal[p] * (1. + deform_avg[p][p]) / plist1[n].detD;
                
                 plist1[n].B_anal[0] = plist0[n].B_anal[0] * (def_fac + deform_avg[0][0]) + plist0[n].B_anal[1] * deform_avg[0][1] + plist0[n].B_anal[2] * deform_avg[0][2]; 
                 plist1[n].B_anal[1] = plist0[n].B_anal[1] * (def_fac + deform_avg[1][1]) + plist0[n].B_anal[0] * deform_avg[1][0] + plist0[n].B_anal[2] * deform_avg[1][2];
                 plist1[n].B_anal[2] = plist0[n].B_anal[2] * (def_fac + deform_avg[2][2]) + plist0[n].B_anal[0] * deform_avg[2][0] + plist0[n].B_anal[1] * deform_avg[2][1];

                for(p=0; p<3; p++)
                  plist1[n].B_anal[p] = plist1[n].B_anal[p] / plist1[n].detD;

                 double B_phys = plist0[n].B_anal[2] * pow(cosmo_fac0 /cosmo_fac, 2);
                 double error = (plist1[n].B_anal[2] - B_phys) / B_phys;

                 if(n < 5) 
                 printf("n = %d, n_neighb = %d, dens0 = %lg, nh_anal0 = %lg, dens1 = %lg, nh_anal1 = %lg, dens_alt1 = %lg, hsm0 = %lg, hsm1 = %lg, deform = %lg, detD = %lg, detD_anal = %lg, ERROR = %lg, B_anal[0] = %lg, B_anal[1] = %lg, B_anal[2] = %lg\n",
                    n, plist1[n].NumNeighb, plist0[n].dens, plist0[n].nh_anal, plist1[n].dens,  plist1[n].nh_anal, plist1[n].dens_alt, plist0[n].hsm, plist1[n].hsm, deform_avg[0][0],
                    plist0[n].dens/plist1[n].dens, plist1[n].detD, error, plist1[n].B_anal[0], plist1[n].B_anal[1], plist1[n].B_anal[2]);

                 if(j == snapnum-2 && myrank == 0)
                    fprintf(outfile, "%lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", Time, n, plist1[n].dens, 
                   plist1[n].nh_anal, Bfield[j], plist1[n].B_anal[0], plist1[n].B_anal[1], plist1[n].B_anal[2], 
                   pow(plist1[n].B_anal[0]*plist1[n].B_anal[0] + plist1[n].B_anal[1]*plist1[n].B_anal[1] + plist1[n].B_anal[2]*plist1[n].B_anal[2], 0.5), 
                   plist1[n].dis_avg_all, plist1[n].dis_avg_all, plist1[n].hsm, plist1[n].NumCalc, plist1[n].detD, plist0[n].dens/plist1[n].dens, 
                   plist1[n].dens_alt, plist0[n].dens_alt/plist1[n].dens_alt, plist1[n].dis_avg_all, plist1[n].dis_avg_all, plist1[n].veltot, plist1[n].temp, 
                   plist1[n].dens - plist0[n].dens, plist1[n].pos[0],  plist1[n].pos[1],  plist1[n].pos[2], plist1[n].nh_anal_alt);
                }

                free(plist0);
                free(P0);

                Ngas0 = Ngas;
                reassign_P(Ngas0, arrnum, nmin, nmax);
          }
          if(myrank == 0)
            fclose(outfile); 

          send_size = Ngas;
          if(write_file == 1)
            {

            printf("myrank = %d, nbad = %d, line 722\n", myrank, nbad);

            MPI_Barrier(MPI_COMM_WORLD);

            if(myrank == 0)
              {         
              sprintf(output_fname2, "%s/%s_bfield_%04d", pathout, basename, snapshot_number);

              if(which_sim == 1)
                {
                sprintf(output_fname2, "%s/%s_bfield_ref_%04d", pathout, basename, snapshot_number);
                if(snapshot_number > 9999) sprintf(output_fname2, "%s/%s_bfield_ref_%05d", pathout, basename, snapshot_number);
                }
              if(which_sim == 2)
                {
                sprintf(output_fname2, "%s/%s_bfield_ref2_%04d", pathout, basename, snapshot_number);
                if(snapshot_number > 9999) sprintf(output_fname2, "%s/%s_bfield_ref2_%05d", pathout, basename, snapshot_number);
                }
              if(which_sim == 3)
                sprintf(output_fname2, "%s/%s_bfield_ref3_%04d", pathout, basename, snapshot_number); 

              outfile2=fopen(output_fname2,"w");
              fclose(outfile2);
              outfile2=fopen(output_fname2,"a");
              }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&PBfieldx[0], &Ptot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myrank == 0)
              for(i = 0; i < Ngas; i++)
                 {
                 fwrite(&Ptot[i], sizeof(double), 1, outfile2);
                 if(i%1000 == 0)
		     printf("Ptot = %lg, Bfieldx = %lg\n", Ptot[i], PBfieldx[i]);
                 }

            MPI_Barrier(MPI_COMM_WORLD);         
            MPI_Reduce(&PBfieldy[0], &Ptot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myrank == 0)
              for(i = 0; i < Ngas; i++)
                 {
                 fwrite(&Ptot[i], sizeof(double), 1, outfile2);
                 if(i%1000 == 0)
                     printf("Ptot = %lg, Bfieldy = %lg\n", Ptot[i], PBfieldy[i]);
                 }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&PBfieldz[0], &Ptot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myrank == 0)
              for(i = 0; i < Ngas; i++)
                 {
                 fwrite(&Ptot[i], sizeof(double), 1, outfile2);
                 if(P[i].marker > 0)
                    printf("Ptot = %lg, Bfieldz = %lg\n", Ptot[i], PBfieldz[i]);
                 }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(&Pnh_anal[0], &Ptot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myrank == 0)
              for(i = 0; i < Ngas; i++)
                {
                fwrite(&Ptot[i], sizeof(double), 1, outfile2);
                 if(P[i].marker > 0)
                    printf("Ptot = %lg, nh = %lg, nh_anal = %lg\n", Ptot[i], P[i].Rho, Pnh_anal[i]);
                 if(Ptot[i] < 1.e-30)
                    printf("UH OH Ptot = %lg, nh = %lg, nh_anal = %lg\n", Ptot[i], P[i].Rho, Pnh_anal[i]);
                }

            if(myrank == 0)
              fclose(outfile2);
            }

          //printf("line 685\n");

          free(nb);
          free(plist1);
          free(P);

          if(write_file == 1)
          {
          free(Ptot); 
          free(PBfieldx); 
          free(PBfieldy); 
          free(PBfieldz); 
          free(Pnh_anal);
          }

          //printf("line 763\n");
  } 

  ierr=MPI_Finalize();

}





/* here the particle data is at your disposal 
 */



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
  HubbleParam= hubble_param;


  for(i=0; i<NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
 
       MeanWeight=1.2195;
       //h2frac=2.0*P[i].H2I;
       h2frac = 2.0*1.e-3;

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

	  gamma= 5.0/3.0;
          //gamma=P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].nh  = P[i].Rho * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
          P[i].Rho_phys = P[i].Rho * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
          P[i].Rho = P[i].Rho * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  /*  printf("zred = %g", zred);*/
	}
    }

return(0);

}






/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files, int myrank)
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

      if(myrank == 0)
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
	allocate_memory(myrank);

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
      if(myrank == 0)
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
 
      if(myrank == 0)
      {
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
      }

      if(header1.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Temp, sizeof(double), 1, fd);
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

/*
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
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
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;
*/
	}


      fclose(fd);
    }

  Time= header1.time;
  zred= header1.redshift;
  if(myrank == 0)
  {
  printf("z= %6.2f \n",zred);
  printf("Time= %12.10e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  }
  return(Ngas);
}




/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(int myrank)
{
  if(myrank == 0)
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
   if(myrank == 0)
    printf("allocating memory...done\n");

return(0);

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
  
  return(0);
}

void reassign_P(int numpart, int arrnum, int nmin, int nmax)
{
int i, n, p;

            if(!(P0=(struct particle_data0 *) malloc(numpart*sizeof(struct particle_data0))))
              {
              fprintf(stderr,"failed to allocate memory.\n");
              exit(0);
              }
            for(i = 0; i < numpart; i++)
               {
               P0[i].Pos[0] = P[i].Pos[0];
               P0[i].Pos[1] = P[i].Pos[1];
               P0[i].Pos[2] = P[i].Pos[2];
               P0[i].Vel[0] = P[i].Vel[0];
               P0[i].Vel[1] = P[i].Vel[1];
               P0[i].Vel[2] = P[i].Vel[2];
               P0[i].Rho = P[i].Rho;
               P0[i].Id = P[i].Id;
               }
           if(!(plist0=(struct plist_data0 *) malloc(arrnum*sizeof(struct plist_data0))))
              {
              fprintf(stderr,"failed to allocate memory.\n");
              exit(0);
              }
           for(n=0; n<arrnum; n++)
              {
                plist0[n].dens = plist1[n].dens;
                plist0[n].dens_alt = plist1[n].dens_alt;
                plist0[n].hsm = plist1[n].hsm;
                plist0[n].pos[0] = plist1[n].pos[0];
                plist0[n].pos[1] = plist1[n].pos[1];
                plist0[n].pos[2] = plist1[n].pos[2];
                plist0[n].vel[0] = plist1[n].vel[0];
                plist0[n].vel[1] = plist1[n].vel[1];
                plist0[n].vel[2] = plist1[n].vel[2];
                plist0[n].veltot = plist1[n].veltot;
                plist0[n].dis_avg_all = plist1[n].dis_avg_all;
                plist0[n].nh_anal = plist1[n].nh_anal;
                plist0[n].nh_anal_alt = plist1[n].nh_anal_alt;
                for(p=0; p<3; p++)
                   plist0[n].B_anal[p] = plist1[n].B_anal[p];
              }

}

void reassign_plist()
{
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


/*
 determinant =  matrix[0][0]*(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])
              - matrix[0][1]*(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0])
              + matrix[0][2]*(matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]);
*/

determinant =       matrix[0][0]*matrix[1][1]*matrix[2][2]
                  + matrix[0][1]*matrix[1][2]*matrix[2][0]
                  + matrix[0][2]*matrix[1][0]*matrix[2][1]
                  - matrix[0][2]*matrix[1][1]*matrix[2][0]
                  - matrix[0][1]*matrix[1][0]*matrix[2][2]
                  - matrix[0][0]*matrix[1][2]*matrix[2][1];


//Alternative calculation for determinant of the deformation tensor
//Tested this alternative, and it is correct!  (Especially for small density changes.)  
//Also, it does NOT yield negative densities like the original formulation did.  


for(p=0; p<3; p++)
  for(q=0; q<3; q++)
     {
     A2[p][q] = 0;
     A3[p][q] = 0;
     }

for(p=0; p<3; p++)
  for(q=0; q<3; q++)
     {
     for(r=0; r<3; r++)
         A2[p][q] = A2[p][q] + matrix[p][r]*matrix[r][q];
     }

for(p=0; p<3; p++)
   for(q=0; q<3; q++)
        {
         for(r=0; r<3; r++)
             A3[p][q] = A3[p][q] + A2[p][r]*matrix[r][q];
        }

trA1 = matrix[0][0] + matrix[1][1] + matrix[2][2];
trA2 = A2[0][0] + A2[1][1] + A2[2][2];
trA3 = A3[0][0] + A3[1][1] + A3[2][2];
/*
determinant = 1. + trA1 + (1./2.)*pow(trA1, 2) - (1./2.)*trA2 + (1./6.)*pow(trA1, 3) + (1./3.)*trA3 - (1./2.)*trA1*trA2
                 + (1./8.)*pow(trA2,2) + (1./3.)*trA1*trA3 - (1./6.)*pow(trA1,2)*trA2 - (1./12.)*pow(trA1,2)*trA2;
               //+ (1./8.)*pow(trA2,2) + (1./3.)*trA1*trA3 + (1./6.)*pow(trA1,2)*trA2 - (1./12.)*pow(trA1,2)*trA2;
*/
return(determinant);
}


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
    kernel = fac*(8./3.14159/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./3.14159/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;

  return(kernel);
}
  










