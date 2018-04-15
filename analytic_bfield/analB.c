#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define hubble_param 0.7
#define MAXREF 20
#define width 200.0
#define snapnum 15000
//#define snapnum 15000
#define arrnum 5
#define neighbnum 200 

int load_snapshot(char *fname, int files);
int reordering(void);
int unit_conversion(void);
int rotate(double x1, double y1, double z1, double vx1, double vy1, double vz1);
int allocate_memory(void);
double calc_det(double matrix11, double matrix12, double matrix13, double matrix21, double matrix22, double matrix23, 
                double matrix31, double matrix32, double matrix33);
void reassign_P(int numpart);
void reassign_P1(int numpart);

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



int     NumPart, Ngas, Ngas0, Ngas1;

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

struct particle_data0
{
  double  Pos[3];
  int Id;
  double Vel[3];
  double  Mass;;
  double  Rho;
} *P0;

struct particle_data1
{
  double  Pos[3];
  int Id;
  double Vel[3];
  double  Mass;
  double  Rho;
} *P1;

struct plist_data0
{
  double vel[3], veltot; 
  double pos[3];
  double hsm; 
  double dens_conv, dens_alt;
  double dis_avg_all; 
  double NumCalc;
  double dens;
  double temp;
  double B_anal[3], nh_anal, nh_anal_alt;
} *plist0;

struct plist_data2
{
  double vel[3], veltot;
  double pos[3];
  double hsm;
  double detD;
  double dens_conv, dens_alt;
  double dis_avg_all;
  double NumCalc;
  double dens;
  double temp;
  double B_anal[3], nh_anal, nh_anal_alt;
} *plist2;

struct plist_data1
{
  double vel[3], veltot;
  double pos[3];
  double hsm;
  int    NumNeighb, NumMatch;
  double detD;
  double dens_conv, dens_alt;
  double dis_avg, dis_avg_all, delx_avg, dxdx_avg, delnh_avg;
  double vel_avg0, vel_avg1, vel_disp0, vel_disp1;
  double NumCalc;
  double dens;
  double temp;
  double B_anal[3], nh_anal, nh_anal_alt;
} *plist1;


struct neighb_data
{
  double vel0[3][neighbnum], vel1[3][neighbnum], vel2[3][neighbnum];  
  double pos0[3][neighbnum], pos1[3][neighbnum], pos2[3][neighbnum];
  int id[neighbnum];
  double nh1[neighbnum], nh0[neighbnum], nh2[neighbnum];
  double deform[3][3][neighbnum];
} *nb;
 
 

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
  char path[200], input_fname[200], plist_fname[200], output_fname[200], output_fname2[300], basename[200], basename2[200], basenameout[200];
  int  ID, i, k=0, ksort, j, n, p, q, r, s, type, snapshot_number, files, Ngas, Ngas0, Ngas1, random, ncount, ncounthalo1;
  double x,y,z,x1,y1,z1, delr, arrnum_doub = 10.0, neighbnum_doub, hfac = 0.75, hfac_in = 0.01, hfac_avg, mfac = 0.0, time_fac, cosmo_fac;
  double nh, nhmax, mass, mmax, dis_min, dis_sim, disAU, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double vrad, vrot, disAUxy, disAUz, disAUxyz, vx_com, vy_com, vz_com, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4;
  double sinkposx1, sinkposy1, sinkposz1, sinkposx2, sinkposy2, sinkposz2; 
  double num, amom, amomtot, amomx_tot, amomy_tot; 
  double amomz_tot, amomx, amomy, amomz, amomx1, amomx2, amomy1, amomy2, amomz1, amomz2;
  double xfac, yfac, zfac, kfloat, nh_avg[snapnum], nh_avg2, num_tot, omega, omega0, nthresh=1.e-20;
  double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, jzcent;
  double rmin, rmax, nmin, nmax; 
  double B0=1.e-20, Bfield[snapnum], deform_avg[3][3]; 
  double detD_temp, detD_neighb[neighbnum];
  int    snaparr[snapnum], *plist_id;
  double vradx, vrady, vradz, vrotx, vroty, vrotz; 
  double velx, vely, velz;
  double rad, dmin, disx, disy, disz, del[3][3], del0[3][3];
  double vel[neighbnum], vel0[neighbnum], vel2[neighbnum], dis[neighbnum], dis0[neighbnum], dis2[neighbnum], ncount_doub;
  double davgHI0[3], davgLO0[3], davgHI1[3], davgLO1[3], davgHI2[3], davgLO2[3], davgHI3[3], davgLO3[3]; 
  double numHI0[3], numLO0[3], numHI1[3], numLO1[3], numHI2[3], numLO2[3], numHI3[3], numLO3[3];
  FILE *pfile, *outfile, *outfile2;

  hfac_avg = (hfac_in + hfac)*mfac;

  plist_id =  (int *) malloc(arrnum * sizeof(int));

  sprintf(path, "/work/00863/minerva/bin_zoom_cut");
  //sprintf(path, "/scratch/00863/minerva");

  sprintf(basename, "bin_zoom10");
  sprintf(basename2, "bin_zoom10_new_cut");

  //sprintf(basename, "bin_zoom10_ref2");
  //sprintf(basename2, "bin_zoom10_new_cut_ref2");

  sprintf(plist_fname, "%s_plist",basename);
  pfile=fopen(plist_fname, "r");
  for(n=0;n<arrnum;n++)
    {
    fscanf(pfile, "%d\n", &ID);
    plist_id[n] = ID;
    printf("ID = %d \n", plist_id[n]);
    }
  fclose(pfile);

 
  for(n=0;n<snapnum;n++)
    snaparr[n] = 40+(1*n);
    //snaparr[n] = 1003+(1*n); 
    //snaparr[n] = 3000+(1*n); 

  for(n=0;n<snapnum;n++)
    {
    nh_avg[n] = 0;
    printf("snaparr = %d \n", snaparr[n]);
    }

  for(j=0;j<snapnum-1;j++)
  {

  if(!(plist2=(struct plist_data2 *) malloc(arrnum*sizeof(struct plist_data2))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }
  if(!(nb=(struct neighb_data *) malloc(arrnum*sizeof(struct neighb_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }

  snapshot_number= snaparr[j];                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%04d", path, basename2, snapshot_number);
  if(snapshot_number > 999)
    sprintf(input_fname, "%s/%s_%04d", path, basename2, snapshot_number);
  if(snapshot_number > 9999)
    sprintf(input_fname, "%s/%s_%05d", path, basename2, snapshot_number);

  Ngas = load_snapshot(input_fname, files);

  //sprintf(output_fname, "%s_plist_%03d.dat",basename, snapshot_number);
  sprintf(output_fname, "%s_plist.dat_err_HI",basename);

          /*    reordering();*/ /* call this routine only if your ID's are set properly */

  printf("Line 237\n");

  unit_conversion();

  num_tot=ncount=0;
  ncount_doub = h2max = nhmax = mmax = kfloat = 0.0;
  xtimesm = ytimesm = ztimesm = 0.0;
  vxtimesm = vytimesm = vztimesm = 0.0;
  vxCOM=vyCOM=vzCOM=xCOM=yCOM=zCOM=0;
  amomx_tot = amomy_tot = amomz_tot = 0;
  amomtot = omega0 = masstot = 0.0;

  printf("Line 261\n");

  time_fac = Time;
  //For NON-cosmological runs
  //time_fac = 1.0;
  
  cosmo_fac = time_fac/(hubble_param);
	
       for(i = 0; i < Ngas; i++)
        {
       if(P[i].nh > nhmax /*&& P[i].sink < 0.5*/)
   	   {
           nhmax = P[i].nh;
           xmax = P[i].Pos[0];
           ymax = P[i].Pos[1];
           zmax = P[i].Pos[2];
           mmax = P[i].Mass;
           }
         }

  printf("Line 276\n");

        for(n=0;n<Ngas;n++) 
         {
           nh = P[n].nh;
           rad = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
           rad=rad*1.e3*time_fac/(hubble_param);
           disx = fabs((P[n].Pos[0]-xmax))*1.e3*time_fac/(hubble_param);
           disy = fabs((P[n].Pos[1]-ymax))*1.e3*time_fac/(hubble_param);
           disz = fabs((P[n].Pos[2]-zmax))*1.e3*time_fac/(hubble_param);
           disAU=rad*206264.806;
           if(disx < width && disy < width && disz < width)
            {
            vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
            vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
            vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
            ncount_doub = ncount_doub + P[n].Mass;
            }
          if(nh > nhmax/20.0)
            {
            //printf("Id = %d, nh = %lg\n", P[n].Id, nh);
            //printf("%d \n", P[n].Id);
            ncount++;
            }
        }

        vxCOM = vxCOM/ncount_doub;
        vyCOM = vyCOM/ncount_doub;
        vzCOM = vzCOM/ncount_doub;
   
        printf("vxCOM = %lg, vyCOM = %lg, vzCOM = %lg\n", vxCOM, vyCOM, vzCOM);
   
	sinkposx1 = xmax;
	sinkposy1 = ymax;
	sinkposz1 = zmax;


	ncount=0;  
	ncount_doub=0.;
	
             //rotate(sinkposx1,sinkposy1,sinkposz1, vx1, vy1, vz1);

        if(j == 0)
           for(n=0; n<arrnum; n++)
              plist2[n].dens = 0;

        outfile=fopen(output_fname, "a");

          for(i = 0; i < Ngas; i++)
             {
 
             P[i].Pos[0] = P[i].Pos[0]*cosmo_fac;
             P[i].Pos[1] = P[i].Pos[1]*cosmo_fac;
             P[i].Pos[2] = P[i].Pos[2]*cosmo_fac;
             P[i].hsm = P[i].hsm*cosmo_fac;     
 
             P[i].Vel[0] = P[i].Vel[0] - vxCOM;
             P[i].Vel[1] = P[i].Vel[1] - vyCOM;
             P[i].Vel[2] = P[i].Vel[2] - vzCOM;
     
             if(P[i].nh > nthresh) 
               for(n=0; n<arrnum; n++)
                      {
                      if(P[i].Id == plist_id[n])
                         {
                         plist2[n].dens = P[i].Rho;
                         plist2[n].dens_conv = P[i].nh/P[i].Rho / hubble_param / hubble_param / pow(1.e0+zred,3.e0);
                         plist2[n].hsm = P[i].hsm;
                         plist2[n].pos[0] = P[i].Pos[0];
                         plist2[n].pos[1] = P[i].Pos[1];
                         plist2[n].pos[2] = P[i].Pos[2];
                         plist2[n].vel[0] = P[i].Vel[0];
                         plist2[n].vel[1] = P[i].Vel[1];
                         plist2[n].vel[2] = P[i].Vel[2];
                         plist2[n].temp = P[i].Temp;
                         plist2[n].veltot = pow(P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2], 0.5);
                         if(j == 0)
                           {
                           plist2[n].nh_anal = P[i].Rho;
                           plist2[n].nh_anal_alt = P[i].Rho;
                           printf("j = %d, n = %d, nh_anal1 = %lg\n", j, n, plist2[n].nh_anal);
                           for(p=0; p<3; p++)
                             plist2[n].B_anal[p] = 1.e-20;
                           }
                         }
                      }
              }

        printf("line 281\n");

        if(j==0)
          {
          Ngas1 = Ngas;
          reassign_P(Ngas);
          }

        if(j==1)
          {
          Ngas0 = Ngas1;
          reassign_P1(Ngas1);

          Ngas1 = Ngas;
          reassign_P(Ngas);
          }


        if(j > 1)
         {

           printf("Ngas = %d, Ngas1 = %d, Ngas0 = %d\n", Ngas, Ngas1, Ngas0);
           for(n=0; n<arrnum; n++)
            {
            plist1[n].NumNeighb = 1;
            plist1[n].NumMatch = 1;
            plist1[n].NumCalc = 0.0;
            plist1[n].dis_avg = plist1[n].dis_avg_all = 0.;
            plist1[n].delx_avg = 0.;
            plist1[n].dxdx_avg = 0.;
            plist1[n].delnh_avg = 0.;
            plist1[n].vel_avg0 = plist1[n].vel_avg1 = plist1[n].vel_disp0 = plist1[n].vel_disp1 = 0;
            plist1[n].dens_alt = 0;
            }

          for(n=0; n<arrnum; n++)
            {
             dis_min = 100.;
             for(i = 0; i < Ngas1; i++)
                {
                 rad = pow(((P1[i].Pos[0]-plist1[n].pos[0])*(P1[i].Pos[0]-plist1[n].pos[0]) + (P1[i].Pos[1]-plist1[n].pos[1])*(P1[i].Pos[1]-plist1[n].pos[1]) + (P1[i].Pos[2]-plist1[n].pos[2])*(P1[i].Pos[2]-plist1[n].pos[2])), 0.5);
                 disx =  P1[i].Pos[0] - plist1[n].pos[0];
                 disy =  P1[i].Pos[1] - plist1[n].pos[1];
                 disz =  P1[i].Pos[2] - plist1[n].pos[2];
              
                 //if((dis < dis_min) && (dis > 0) && (disx0 > disx) && (disy0 > disy) && (disz0 > disz))
                 //if(dis < dis_min && dis > 0)
                 if(rad < 1.2*hfac*plist1[n].hsm && rad > hfac_in*plist1[n].hsm && plist1[n].NumNeighb < neighbnum-1)
                 //if(fabs(disx) < hfac*plist_hsm[j][n] && fabs(disy) < hfac*plist_hsm[j][n] && fabs(disz) < hfac*plist_hsm[j][n] && dis > 0)
                     {
                     dis_min = rad;
                     k = plist1[n].NumNeighb -1;
                     plist1[n].dens_alt = plist1[n].dens_alt + P1[i].Mass;
                     nb[n].pos1[0][k] = P1[i].Pos[0];
                     nb[n].pos1[1][k] = P1[i].Pos[1];
                     nb[n].pos1[2][k] = P1[i].Pos[2];
                     nb[n].vel1[0][k] = P1[i].Vel[0];
                     nb[n].vel1[1][k] = P1[i].Vel[1];
                     nb[n].vel1[2][k] = P1[i].Vel[2];
                     nb[n].nh1[k] = P1[i].Rho;
                     nb[n].id[k] = P1[i].Id;
                     plist1[n].NumNeighb++;
                     }
                }
             }


          for(i = 0; i < Ngas0; i++)
             {
             if(P0[i].Rho > nthresh*1.67e-24*.1)
             for(n=0; n<arrnum; n++)
               for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  if(P0[i].Id == nb[n].id[k])
                     {                          
                     nb[n].pos0[0][k] = P0[i].Pos[0];
                     nb[n].pos0[1][k] = P0[i].Pos[1];
                     nb[n].pos0[2][k] = P0[i].Pos[2];
                     nb[n].vel0[0][k] = P0[i].Vel[0];
                     nb[n].vel0[1][k] = P0[i].Vel[1];
                     nb[n].vel0[2][k] = P0[i].Vel[2];
                     nb[n].nh0[k] = P0[i].Rho;
                     plist1[n].NumMatch++;
                     }
                    }
              }



            for(n=0; n<arrnum; n++)
                {
                printf("NUMMATCH = %d\n", plist1[n].NumMatch);
                plist1[n].detD = plist2[n].detD = 0.0;
                for(p=0; p<3; p++)
                  {
                  davgHI0[p]=davgLO0[p]=davgHI1[p]=davgLO1[p]=davgHI2[p]=davgLO2[p]=davgHI3[p]=davgLO3[p]=0;
                  numHI0[p]=numLO0[p]=numHI1[p]=numLO1[p]=numHI2[p]=numLO2[p]=numHI3[p]=numLO3[p]=0;
                  }
                dmin = 0.1*plist1[n].hsm;
                for(p=0; p<3; p++)
                 for(q=0; q<3; q++)
                    deform_avg[p][q] = 0;
                for(k=0; k<neighbnum; k++)
                  detD_neighb[k] = 0.;

                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {                
                  dis0[k] = pow(((nb[n].pos0[0][k]-plist0[n].pos[0])*(nb[n].pos0[0][k]-plist0[n].pos[0])
                       + (nb[n].pos0[1][k]-plist0[n].pos[1])*(nb[n].pos0[1][k]-plist0[n].pos[1])
                       + (nb[n].pos0[2][k]-plist0[n].pos[2])*(nb[n].pos0[2][k]-plist0[n].pos[2])), 0.5);

                  dis[k] = pow(((nb[n].pos1[0][k]-plist1[n].pos[0])*(nb[n].pos1[0][k]-plist1[n].pos[0])
                       + (nb[n].pos1[1][k]-plist1[n].pos[1])*(nb[n].pos1[1][k]-plist1[n].pos[1])
                       + (nb[n].pos1[2][k]-plist1[n].pos[2])*(nb[n].pos1[2][k]-plist1[n].pos[2])), 0.5);


                  vel0[k] = pow(nb[n].vel0[0][k]*nb[n].vel0[0][k] + nb[n].vel0[1][k]*nb[n].vel0[1][k]
                       + nb[n].vel0[2][k]*nb[n].vel0[2][k], 0.5);

                  vel[k] = pow(nb[n].vel1[0][k]*nb[n].vel1[0][k] + nb[n].vel1[1][k]*nb[n].vel1[1][k]
                       + nb[n].vel1[2][k]*nb[n].vel1[2][k], 0.5);

                  plist1[n].dis_avg_all = plist1[n].dis_avg_all + dis[k];
                  plist1[n].dxdx_avg  = plist1[n].dxdx_avg + fabs(dis0[k]-dis[k])/dis0[k];
                  plist1[n].delx_avg  = plist1[n].delx_avg + fabs(dis0[k]-dis[k]);
                  plist1[n].delnh_avg = plist1[n].delnh_avg + fabs(nb[n].nh1[k]-nb[n].nh0[k]);
                  plist1[n].vel_avg1   = plist1[n].vel_avg1 + vel[k];
                  plist1[n].vel_avg0   = plist1[n].vel_avg0 + vel0[k]; 
                  }
 
                plist1[n].dis_avg_all = plist1[n].dis_avg_all/(double(plist1[n].NumNeighb-1));
                plist1[n].dxdx_avg    = plist1[n].dxdx_avg/(double(plist1[n].NumNeighb-1));
                plist1[n].delx_avg    = plist1[n].delx_avg/(double(plist1[n].NumNeighb-1));
                plist1[n].vel_avg1    = plist1[n].vel_avg1/(double(plist1[n].NumNeighb-1));
                plist1[n].vel_avg0    = plist1[n].vel_avg0/(double(plist1[n].NumNeighb-1));
                plist1[n].delnh_avg   = plist1[n].delnh_avg/(double(plist1[n].NumNeighb-1));

                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  plist1[n].vel_disp0 = plist1[n].vel_disp0 + fabs(vel[k] - plist1[n].vel_avg0);
                  plist1[n].vel_disp1 = plist1[n].vel_disp1 + fabs(vel[k] - plist1[n].vel_avg1);
                  }
                plist1[n].vel_disp0  = plist1[n].vel_disp0/(double(plist1[n].NumNeighb-1));
                plist1[n].vel_disp1  = plist1[n].vel_disp1/(double(plist1[n].NumNeighb-1));

                for(k=0; k<plist1[n].NumNeighb-1; k++)
                  {
                  for(p=0; p<3; p++)
                     for(q=0; q<3; q++)
                       {
                       nb[n].deform[p][q][k] =  (-(plist0[n].pos[p] - nb[n].pos0[p][k]) + (plist1[n].pos[p] - nb[n].pos1[p][k]))
                                         / (plist0[n].pos[q]  - nb[n].pos0[q][k]);
                       
                       del0[p][q] = plist0[n].pos[p] - nb[n].pos0[q][k];
                       del[p][q] = plist1[n].pos[p] - nb[n].pos1[q][k];
                       }

                   detD_neighb[k] = calc_det(nb[n].deform[0][0][k], nb[n].deform[0][1][k], nb[n].deform[0][2][k], nb[n].deform[1][0][k], 
                                   nb[n].deform[1][1][k], nb[n].deform[1][2][k], nb[n].deform[2][0][k], nb[n].deform[2][1][k], nb[n].deform[2][2][k]);
  
                   if(dis0[k] < hfac*plist0[n].hsm  && dis[k] < hfac*plist1[n].hsm)
                    { 
                    plist1[n].detD = plist1[n].detD + detD_neighb[k];
                    plist1[n].dis_avg = plist1[n].dis_avg + dis[k];
                    //plist1[n].NumCalc = plist1[n].NumCalc + 1.0;
                    }


                  if(dis0[k] < hfac*plist0[n].hsm
                    //&& fabs(del0[0][0])  > dmin && fabs(del0[0][1])  > dmin && fabs(del0[0][2])  > dmin  //alt8b and alt6b
                    //&& fabs(del0[1][0])  > dmin && fabs(del0[1][1])  > dmin && fabs(del0[1][2])  > dmin
                    //&& fabs(del0[2][0])  > dmin && fabs(del0[2][1])  > dmin && fabs(del0[2][2])  > dmin
                    //&& fabs(vel[k] - plist1[n].vel_avg1) < 1.0*plist1[n].vel_disp1 && fabs(vel0[k] - plist1[n].vel_avg0) < 1.0*plist1[n].vel_disp0
                          //^vdisp6; no good!
                    //&& fabs(del0[0][0])  > dmin && fabs(del0[1][1])  > dmin  && fabs(del0[2][2]) > dmin  //alt6_dmin
                              //^not much help!
                  
/* 
                    && fabs(nb[n].vel0[0][k]/plist0[n].vel[0]) < 2.0 && fabs(nb[n].vel0[0][k]/plist0[n].vel[0]) > 0.5
                    && fabs(nb[n].vel0[1][k]/plist0[n].vel[1]) < 2.0 && fabs(nb[n].vel0[1][k]/plist0[n].vel[1]) > 0.5
                    && fabs(nb[n].vel0[2][k]/plist0[n].vel[2]) < 2.0 && fabs(nb[n].vel0[2][k]/plist0[n].vel[2]) > 0.5
                    && fabs(nb[n].vel1[0][k]/plist1[n].vel[0]) < 2.0 && fabs(nb[n].vel1[0][k]/plist1[n].vel[0]) > 0.5
                    && fabs(nb[n].vel1[1][k]/plist1[n].vel[1]) < 2.0 && fabs(nb[n].vel1[1][k]/plist1[n].vel[1]) > 0.5
                    && fabs(nb[n].vel1[2][k]/plist1[n].vel[2]) < 2.0 && fabs(nb[n].vel1[2][k]/plist1[n].vel[2]) > 0.5
*/
                    )
                     {
                     for(p=0; p<3; p++)
                      { 
                       if(fabs(plist0[n].pos[p] - nb[n].pos0[p][k]) < hfac_avg*plist0[n].hsm)  
                                                      //alt1,mfac=0.5; alt3,mfac=0; alt4,mfac=0,no "fabs"; alt5,mfac=0.5,hfac=1.0 
                                                      //alt6,mfac=0,hfac=1; alt7,mfac=0,hfac=1.5,hfac_in=0.5; 
                                                      //alt8,hfac=0.75,hfac_in=.01,mfac=0
                                                      //if(dis0[k] < 0.5*hfac*plist0[n].hsm)             //alt2
                                                      //alt9, hfac=2.0,mfac=0
                         {
                         davgLO0[p] = davgLO0[p] + fabs(plist0[n].pos[p] - nb[n].pos0[p][k]);
                         numLO0[p] = numLO0[p] + 1.0; 
                         }
                        if(fabs(plist1[n].pos[p] - nb[n].pos1[p][k]) < hfac_avg*plist1[n].hsm)
                         {
                         davgLO1[p] = davgLO1[p] + fabs(plist1[n].pos[p] - nb[n].pos1[p][k]);
                         numLO1[p] = numLO1[p] + 1.0;
                         }
                       if(fabs(plist0[n].pos[p] - nb[n].pos0[p][k]) > hfac_avg*plist0[n].hsm)
                         {
                         davgHI0[p] = davgHI0[p] + fabs(plist0[n].pos[p] - nb[n].pos0[p][k]);
                         numHI0[p] = numHI0[p] + 1.0;
                         }
                        if(fabs(plist1[n].pos[p] - nb[n].pos1[p][k]) > hfac_avg*plist1[n].hsm)
                         {
                         davgHI1[p] = davgHI1[p] + fabs(plist1[n].pos[p] - nb[n].pos1[p][k]);
                         numHI1[p] = numHI1[p] + 1.0;
                         } 
                      }
                   plist1[n].NumCalc = plist1[n].NumCalc + 1.0;
                   }
                  if(dis[k] < hfac*plist1[n].hsm)
                   {
                     for(p=0; p<3; p++)
                      {
                        if(fabs(plist0[n].pos[p] - nb[n].pos0[p][k]) < hfac_avg*plist0[n].hsm)
                         {
                         davgLO2[p] = davgLO2[p] + fabs(plist0[n].pos[p] - nb[n].pos0[p][k]);
                         numLO2[p] = numLO2[p] + 1.0;
                         }
                        if(fabs(plist1[n].pos[p] - nb[n].pos1[p][k]) < hfac_avg*plist1[n].hsm)
                         {
                         davgLO3[p] = davgLO3[p] + fabs(plist1[n].pos[p] - nb[n].pos1[p][k]);
                         numLO3[p] = numLO3[p] + 1.0;
                         }
                        if(fabs(plist0[n].pos[p] - nb[n].pos0[p][k]) > hfac_avg*plist0[n].hsm)
                         {
                         davgHI2[p] = davgHI2[p] + fabs(plist0[n].pos[p] - nb[n].pos0[p][k]);
                         numHI2[p] = numHI2[p] + 1.0;
                         }
                        if(fabs(plist1[n].pos[p] - nb[n].pos1[p][k]) > hfac_avg*plist1[n].hsm)
                         {
                         davgHI3[p] = davgHI3[p] + fabs(plist1[n].pos[p] - nb[n].pos1[p][k]);
                         numHI3[p] = numHI3[p] + 1.0;
                         }
                      }
                   }
                 }                

                for(p=0; p<3; p++)
                  {
                  if(numLO0[p] <= 0) numLO0[p] = 1.0;
                  if(numLO1[p] <= 0) numLO1[p] = 1.0;
                  if(numLO2[p] <= 0) numLO2[p] = 1.0;
                  if(numLO3[p] <= 0) numLO3[p] = 1.0;
                  if(numHI0[p] <= 0) numHI0[p] = 1.0;
                  if(numHI1[p] <= 0) numHI1[p] = 1.0;
                  if(numHI2[p] <= 0) numHI2[p] = 1.0;
                  if(numHI3[p] <= 0) numHI3[p] = 1.0;
                  }

                for(p=0; p<3; p++)
                  {
                  printf("davgLO0= %lg, davgLO1 = %lg, davgLO2 = %lg, davgLO3 = %lg, davgHI0=%lg, davgHI1 = %lg, davgHI2 = %lg, davgHI3 = %lg, numLO0 = %lg, numLO1 = %lg, numLO2 = %lg, numLO3 = %lgnumHI0 = %lg, numHI1 = %lg, numHI2 = %lg, numHI3 = %lg\n",
                          davgLO0[p], davgLO1[p], davgLO2[p], davgLO3[p], davgHI0[p], davgHI1[p], davgHI2[p], davgHI3[p], numLO0[p], numLO1[p], numLO2[p], numLO3[p], numHI0[p], numHI1[p], numHI2[p], numHI3[p]);
                  davgLO0[p] = davgLO0[p]/numLO0[p];
                  davgLO1[p] = davgLO1[p]/numLO1[p];
                  davgLO2[p] = davgLO2[p]/numLO2[p];
                  davgLO3[p] = davgLO3[p]/numLO3[p];
                  davgHI0[p] = davgHI0[p]/numHI0[p];
                  davgHI1[p] = davgHI1[p]/numHI1[p];
                  davgHI2[p] = davgHI2[p]/numHI2[p];
                  davgHI3[p] = davgHI3[p]/numHI3[p];
                  }
 
                for(p=0; p<3; p++)
                   for(q=0; q<3; q++)
                     deform_avg[p][q] =  (-(davgHI0[p] - davgLO0[p]) + (davgHI1[p] - davgLO1[p]))/(davgHI0[q] - davgLO0[q]);
 
                 plist1[n].detD = calc_det(deform_avg[0][0], deform_avg[0][1], deform_avg[0][2], deform_avg[1][0],
                                    deform_avg[1][1], deform_avg[1][2], deform_avg[2][0], deform_avg[2][1], deform_avg[2][2]);

                for(p=0; p<3; p++)
                   for(q=0; q<3; q++)
                     //deform_avg[p][q] =  (-(davgHI2[p] - davgLO2[p]) + (davgHI3[p] - davgLO3[p]))/(davgHI2[q] - davgLO2[q]);
                     deform_avg[p][q] =  (-(davgHI3[p] - davgLO3[p]) + (davgHI2[p] - davgLO2[p]))/(davgHI3[q] - davgLO3[q]);

                 plist2[n].detD = calc_det(deform_avg[0][0], deform_avg[0][1], deform_avg[0][2], deform_avg[1][0],
                                    deform_avg[1][1], deform_avg[1][2], deform_avg[2][0], deform_avg[2][1], deform_avg[2][2]);

                 plist2[n].detD = 1./plist2[n].detD; 
                 //plist2[n].detD = plist1[n].detD;

                 if(plist1[n].NumCalc <= 1.)
                   {
                   printf("NO NEARBY PARTICLES, n = %d\n", n);
                   plist1[n].detD = 1.;
                   plist2[n].detD = 1.;
                   plist1[n].NumCalc = 1.;
                   }

                 for(p=0; p<3; p++)
                   if(numHI1[p] <= 1 || numHI0[p] <= 1)
                     {
                     plist1[n].detD = plist2[n].detD = 1;
                     deform_avg[p][p] = 0;
                     }

                 if(plist0[n].nh_anal != plist0[n].nh_anal || plist0[n].nh_anal < 1.e-40)
                   plist0[n].nh_anal = plist0[n].dens;

                 plist1[n].nh_anal  = plist0[n].nh_anal / ((plist1[n].detD + plist2[n].detD) / 2.);
                 if(plist0[n].dis_avg_all > 0)
                   plist1[n].nh_anal_alt = plist0[n].nh_anal_alt * pow(plist0[n].dis_avg_all/plist1[n].dis_avg_all,3);
                 plist1[n].dis_avg  = plist1[n].dis_avg/plist1[n].NumCalc;
                 plist1[n].dens_alt = plist1[n].dens_alt / (4.*3.14159 /3) / pow(hfac*plist1[n].hsm, 3); 
                 plist1[n].dens_alt = plist1[n].dens_alt*plist1[n].dens_conv;

                 for(p=0; p<3; p++)
                   plist1[n].B_anal[p] = plist0[n].B_anal[p] * (1. + deform_avg[p][p]) / plist1[n].detD;
                 
                 printf("n = %d, dens0 = %lg, nh_anal0 = %lg, dens1 = %lg, nh_anal1 = %lg, hsm0 = %lg, hsm1 = %lg, deform = %lg, detD = %lg, detD_anal = %lg,detd_anal2 = %lg,  detD_alt = %lg\n",
                    n, plist0[n].dens, plist0[n].nh_anal, plist1[n].dens,  plist1[n].nh_anal, plist0[n].hsm, plist1[n].hsm, nb[n].deform[0][0][0],
                    plist0[n].dens/plist1[n].dens, plist1[n].detD, plist2[n].detD, plist0[n].dens_alt/plist1[n].dens_alt);

                 fprintf(outfile, "%lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", Time, n, plist1[n].dens, 
                   plist1[n].nh_anal, Bfield[j], plist1[n].B_anal[0], plist1[n].B_anal[1], plist1[n].B_anal[2], 
                   pow(plist1[n].B_anal[0]*plist1[n].B_anal[0] + plist1[n].B_anal[1]*plist1[n].B_anal[1] + plist1[n].B_anal[2]*plist1[n].B_anal[2], 0.5), 
                   plist1[n].dis_avg, plist1[n].delx_avg, plist1[n].hsm, plist1[n].NumCalc, plist1[n].detD, plist0[n].dens/plist1[n].dens, 
                   plist1[n].dens_alt, plist0[n].dens_alt/plist1[n].dens_alt, plist1[n].vel_avg1, plist1[n].vel_disp1, plist1[n].veltot, plist1[n].temp, 
                   plist1[n].dens - plist0[n].dens, plist1[n].pos[0],  plist1[n].pos[1],  plist1[n].pos[2], plist1[n].nh_anal_alt);
                }

                free(plist0);
                free(P0);
                Ngas0 = Ngas1;
                reassign_P1(Ngas1);

                free(plist1);
                free(P1);
                Ngas1 = Ngas;
                reassign_P(Ngas);

          }
          fclose(outfile); 

          free(nb);
          free(plist2);
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
  dis=dis*1.e3*Time/(hubble_param);                   //dis is in pc
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
  dis=dis*1.e3*Time/(hubble_param);                   //dis is in pc
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
  HubbleParam= hubble_param;


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
	  P[i].nh  = P[i].Rho * UnitDensity_in_cgs *HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
          P[i].Rho = P[i].Rho *UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
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
  int    t,n,off,pc=0,pc_new,pc_sph;

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

void reassign_P(int numpart)
{
int i, n, p, num;


            if(!(P1=(struct particle_data1 *) malloc(numpart*sizeof(struct particle_data1))))
              {
              fprintf(stderr,"failed to allocate memory.\n");
              exit(0);
              }
            for(i = 0; i < numpart; i++)
               {
               P1[i].Pos[0] = P[i].Pos[0];
               P1[i].Pos[1] = P[i].Pos[1];
               P1[i].Pos[2] = P[i].Pos[2];
               P1[i].Vel[0] = P[i].Vel[0];
               P1[i].Vel[1] = P[i].Vel[1];
               P1[i].Vel[2] = P[i].Vel[2];
               P1[i].Rho = P[i].Rho;
               P1[i].Id = P[i].Id;
               P1[i].Mass = P[i].Mass;
               }
           if(!(plist1=(struct plist_data1 *) malloc(arrnum*sizeof(struct plist_data1))))
              {
              fprintf(stderr,"failed to allocate memory.\n");
              exit(0);
              }
           for(n=0; n<arrnum; n++)
              {
                plist1[n].dens = plist2[n].dens;
                plist1[n].dens_alt = plist2[n].dens_alt;
                plist1[n].dens_conv = plist2[n].dens_conv;
                plist1[n].hsm = plist2[n].hsm;
                plist1[n].temp = plist2[n].temp;
                plist1[n].pos[0] = plist2[n].pos[0];
                plist1[n].pos[1] = plist2[n].pos[1];
                plist1[n].pos[2] = plist2[n].pos[2];
                plist1[n].vel[0] = plist2[n].vel[0];
                plist1[n].vel[1] = plist2[n].vel[1];
                plist1[n].vel[2] = plist2[n].vel[2];
                plist1[n].veltot = plist2[n].veltot;
                plist1[n].dis_avg_all = plist2[n].dis_avg_all;
                plist1[n].nh_anal = plist2[n].nh_anal;
                plist1[n].nh_anal_alt = plist2[n].nh_anal_alt;
                for(p=0; p<3; p++)
                   plist1[n].B_anal[p] = plist2[n].B_anal[p];
              }

}


void reassign_P1(int numpart)
{
int i, n, p, num;

            if(!(P0=(struct particle_data0 *) malloc(numpart*sizeof(struct particle_data0))))
              {
              fprintf(stderr,"failed to allocate memory.\n");
              exit(0);
              }
            for(i = 0; i < numpart; i++)
               {
               P0[i].Pos[0] = P1[i].Pos[0];
               P0[i].Pos[1] = P1[i].Pos[1];
               P0[i].Pos[2] = P1[i].Pos[2];
               P0[i].Vel[0] = P1[i].Vel[0];
               P0[i].Vel[1] = P1[i].Vel[1];
               P0[i].Vel[2] = P1[i].Vel[2];
               P0[i].Rho = P1[i].Rho;
               P0[i].Id = P1[i].Id;
               P0[i].Mass = P1[i].Mass;
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
                plist0[n].dens_conv = plist1[n].dens_conv;
                plist0[n].hsm = plist1[n].hsm;
                plist0[n].temp = plist1[n].temp;
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


double calc_det(double matrix00, double matrix01, double matrix02, double matrix10, double matrix11, double matrix12,
                double matrix20, double matrix21, double matrix22)
{
int p, q, r;
double determinant, trA1, trA2, trA3, A2[3][3], A3[3][3], matrix[3][3];

matrix[0][0] = matrix00;
matrix[0][1] = matrix01;
matrix[0][2] = matrix02;
matrix[1][0] = matrix10;
matrix[1][1] = matrix11;
matrix[1][2] = matrix12;
matrix[2][0] = matrix20;
matrix[2][1] = matrix21;
matrix[2][2] = matrix22;

determinant = (1+matrix[0][0])*(1+matrix[1][1])*(1+matrix[2][2])
                  + matrix[0][1]*matrix[1][2]*matrix[2][0]
                  + matrix[0][2]*matrix[1][0]*matrix[2][1]
                  - matrix[0][2]*(1+matrix[1][1])*matrix[2][0]
                  - matrix[0][1]*matrix[1][0]*(1+matrix[2][2])
                  - (1+matrix[0][0])*matrix[1][2]*matrix[2][1];


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




  











