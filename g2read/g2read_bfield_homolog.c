#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define readB 1
#define vpot 0
#define hubble_param 1.0
#define PI 3.14159265359

#define which_sim 0  //homolog
//#define which_sim 1  //sod

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
  //double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy;
  double nh_test, Bfieldx, Bfieldy, Bfieldz, error;
  double Afieldx, Afieldy, Afieldz, fcor;
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
  double x,y,z,x1,y1, z1, delr, vx, vy, vz, vel, typemax, dismax, sq_error=0, rms_error=0, tot_error=0;
  double delfac=0.9, n0, b0, a0, a0_x, a0_y, n_anal, b_anal;
  double nh, nhmax, mass, mmax, dis, xmax, ymax, zmax, sl, masstot, temp, tmax, h2, h2max, gam, gammin;
  double sinkposx, sinkposy, sinkposz, disAU, vrad, vrotx, vroty, vrotz, vrot, mh_mass;
  double disx, disy, disz;
  FILE *outfile, *outfile2, *infile;



  sprintf(path, "/work/00863/minerva");
  sprintf(pathout, "/work/00863/minerva");
  sprintf(path2, "/work/00863/minerva/orion/");

  sprintf(basename2, "homolog");
  sprintf(basename, "homolog");
  sprintf(basenameout, "snaphomolog");

  int arrnum = 0, idnum = 0;
  int snapmin = 6, snapmax = 40, snapbegin = 3, snapend = 98, jinc = 5;
  double width_smallx = 10.e0, width_smally = 10.e0, width_smallz = 10.e0;
  

  if(which_sim == 1)
  {
  sprintf(basename2, "sod");
  sprintf(basename, "sod");
  sprintf(basenameout, "snapsod");

  snapmin = 2, snapmax = 2, snapbegin = 2, snapend = 98, jinc = 1;
  width_smallx = 6.e2;
  width_smally = 1.e1; 
  width_smallz = 1.e1;
  }

  if(which_sim == 2)
  {
  sprintf(basename2, "shear");
  sprintf(basename, "shear");
  sprintf(basenameout, "snapshear");

  snapmin = 2, snapmax = 30, snapbegin = 1, snapend = 98, jinc = 1;
  width_smallx = 2.e1;
  width_smally = 5.e0;
  width_smallz = 2.e1;
  }


for(j=snapbegin;j<=snapbegin;j=j+jinc){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */

  if(j<=999)
  {
  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%03d", pathout, basenameout, snapshot_number);
  }

  if(j>999 || which_sim == 1)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }

  Ngas = load_snapshot(input_fname, files);

  xmax = ymax = zmax = header1.BoxSize/2.0;

  double cen_min = 100;
  if(which_sim == 2) //xmax = ymax = zmax = 0.503049;
  {
    for(i = 0; i < Ngas; i++)
      {
      if(fabs(P[i].Pos[1] - ymax) < cen_min && P[i].Pos[1] > ymax) cen_min = fabs(P[i].Pos[1] - ymax);
      }
  xmax = header1.BoxSize/2.0;
  ymax = header1.BoxSize/2.0 + cen_min;
  zmax = header1.BoxSize/2.0;
  }
  printf("nhmax = %lg, xmax = %lg, ymax = %lg, zmax = %lg\n", nhmax, xmax, ymax, zmax);

  printf("line 117\n");

  for(i = 0; i < Ngas; i++)
     {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3;
     disy = fabs((P[i].Pos[1]-ymax))*1.e3;
     disz = fabs((P[i].Pos[2]-zmax))*1.e3;
     if(disx < width_smallx/2.0 && disy < width_smally/2.0 && disz < width_smallz/2.0)
       idnum++;
     }

  printf("idnum = %d\n", idnum);

  if(!(plist=(struct plist_data *) malloc(idnum*sizeof(struct plist_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }

  n = 0;
  for(i = 0; i < Ngas; i++)
     {
     disx = fabs((P[i].Pos[0]-xmax))*1.e3;
     disy = fabs((P[i].Pos[1]-ymax))*1.e3;
     disz = fabs((P[i].Pos[2]-zmax))*1.e3;
     if(disx < width_smallx/2.0 && disy < width_smally/2.0 && disz < width_smallz/2.0)
        {
        plist[n].Id = P[i].Id;
        n++;
        }
     }


free(P);
} 
 

  for(j=snapmin;j<=snapmax;j=j+1){

  snapshot_number= j;                   /* number of snapshot */
  files=1;                               /* number of files per snapshot */
  arrnum = 0; 

  if(j<=999)
  {
  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  sprintf(output_fname, "%s/%s_%03d", pathout, basenameout, snapshot_number);
  }

  if(j>999 || which_sim == 1)
  {
    sprintf(input_fname, "%s/%s_%04d", path, basename, snapshot_number);
    sprintf(output_fname, "%s/%s_%04d", pathout, basenameout, snapshot_number);
  }


  sprintf(output_fname2, "nh_u_evol");

  Ngas = load_snapshot(input_fname, files);

  /*    reordering();*/ /* call this routine only if your ID's are set properly */
 
  printf("t_Hubble = %lg \n", 5.4e8/pow((1.e0+header1.redshift)/10.e0, 1.5e0));
  printf("hi 1, Ngas = %d, NumPart = %d\n", Ngas, NumPart);

  unit_conversion();  

  printf("hi 2\n");

  outfile=fopen(output_fname, "w");

  ncount = 0;
  ncount2 = 0;

  printf("hi 3\n");

//Read in B-field info!///////////////////////////////////////////////////////////
    sprintf(input_fname2, "%s/%s_bfield_%04d", path2, basename2, snapshot_number);
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
          temp = P[n].Temp;

          if(nh > nhmax)
          //if(mass > mmax)
            {
              //printf("Found the sink!\n");
              nhmax = nh;
              mmax=mass;
              xmax=P[n].Pos[0];
              ymax=P[n].Pos[1];
              zmax=P[n].Pos[2];
              sl = P[n].hsm;
              tmax=P[n].Temp;
              typemax = P[n].Type;
              idmax = P[n].Id;
              vx = P[n].Vel[0];
              vy = P[n].Vel[1];
              vz = P[n].Vel[2];
            }
    }

  xmax = ymax = zmax = header1.BoxSize/2.0;


  double nh_typ=0, b_typ=0, u_typ=0, num_typ=0, bfield, ufield, u_norm, bx_typ=0, by_typ=0, bz_typ=0, ax_typ=0, ay_typ=0, az_typ=0, a_anal=0, a_typ;
  int pfac, pfac2;
  double dwdx, dwdy, dwdz, ax, ay, az;
  double width_currentx, width_currenty, width_currentz;
  int nskip=500, ncheck;
  double rcheck, density_check=0, density_check2, dcor, hcor;

  pfac = pfac2 = 3*(j-snapmin);
  if(which_sim == 2)
   {  
   pfac = 0; 
   delfac = 0.1;
   }
  width_currentx = width_smallx*pow(delfac, j-snapmin);
  width_currenty = width_smally*pow(delfac, j-snapmin);
  width_currentz = width_smallz*pow(delfac, j-snapmin);

/*
  for(n=0;n<Ngas;n++)
      {
      P[n].error = 2. * (P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test);
      disx = fabs((P[n].Pos[0]-xmax))*1.e3;
      disy = fabs((P[n].Pos[1]-ymax))*1.e3;
      disz = fabs((P[n].Pos[2]-zmax))*1.e3;
      if(fabs(P[n].error) > 0 && P[n].nh_test == P[n].nh_test && disx < width_current && disy < width_current && disz < width_current)
        ncheck = n;
      }
*/

   n = 5000;
   for(i=0;i<Ngas;i++)
     {
     if(i == n) continue;
     rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
     if(rcheck > P[n].hsm) continue;
     density_check = density_check + P[i].Mass*calc_kernel_spline(i, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
     }
    density_check2 = 56.*P[n].Mass / (4./3.*PI) * pow(P[n].hsm,-3); 
    dcor = P[n].Density / density_check;
    hcor = pow(dcor, -1./3.);
    printf("density_check = %lg, density_check2 = %lg, P[n].Density = %lg, dcor = %lg, hcor = %lg\n", density_check, density_check2, P[n].Density, dcor, hcor);


  double x_alt;
  for(n=0;n<Ngas;n++)
       {

          P[n].error = 2. * (P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test);
          disx = fabs((P[n].Pos[0]-xmax))*1.e3;
          disy = fabs((P[n].Pos[1]-ymax))*1.e3;
          disz = fabs((P[n].Pos[2]-zmax))*1.e3;

          //if(fabs(P[n].error) > 0 && P[n].nh_test == P[n].nh_test && disx < width_currentx && disy < width_currenty && disz < width_currentz)
          //if(plist[n].disx < width_small/2.0 && plist[n].disy < width_small/2.0 && plist[n].disz < width_small/2.0)
          //^ NOTE!!! Above n indeces do NOT correpsond to same particle IDs from snapshot to snapshot!!!!  
          //Don't use this to look at same set of particles between snapshots!!!
          for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
            {
             arrnum++;
             P[n].Density = 0;

             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {

               x_alt = P[i].Pos[0];
/*
               if(P[n].Pos[0] > header1.BoxSize - P[n].hsm && P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] + header1.BoxSize;
               if(P[n].Pos[0] < P[n].hsm && header1.BoxSize - P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] - header1.BoxSize;
*/
               if(i == n) continue;
               rcheck = pow((P[n].Pos[0] - x_alt)*(P[n].Pos[0] - x_alt) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue;
               P[n].Density = P[n].Density + P[i].Mass * calc_kernel_spline(n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], x_alt,  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
                }

             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {

               x_alt = P[i].Pos[0];
/*
               if(P[n].Pos[0] > header1.BoxSize - P[n].hsm && P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] + header1.BoxSize;
               if(P[n].Pos[0] < P[n].hsm && header1.BoxSize - P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] - header1.BoxSize;
*/
               if(i == n) continue;
               rcheck = pow((P[n].Pos[0] - x_alt)*(P[n].Pos[0] - x_alt) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue;
               P[n].fcor = P[n].fcor + (P[n].hsm /  3. / P[n].Density) * P[i].Mass * del_rho(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], x_alt,  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
                }
             

             if(vpot == 1)
             for(i=0;i<Ngas;i++)
               {

               x_alt = P[i].Pos[0];
/*
               if(P[n].Pos[0] > header1.BoxSize - P[n].hsm && P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] + header1.BoxSize;
               if(P[n].Pos[0] < P[n].hsm && header1.BoxSize - P[i].Pos[0] < P[n].hsm)
                   x_alt = P[i].Pos[0] - header1.BoxSize;
*/
               if(i == n) continue;
               rcheck = pow((P[n].Pos[0] - x_alt)*(P[n].Pos[0] - x_alt) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
               if(rcheck > P[n].hsm) continue; 
               dwdx = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], x_alt, P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[0] - x_alt);
               dwdy = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], x_alt, P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[1] - P[i].Pos[1]);
               dwdz = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], x_alt, P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[2] - P[i].Pos[2]);
               ax = P[n].Afieldx - P[i].Afieldx;
               ay = P[n].Afieldy - P[i].Afieldy;
               az = P[n].Afieldz - P[i].Afieldz;
               //P[n].fcor = 1.0;
               P[n].Bfieldx = P[n].Bfieldx + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(0, ax, ay, az, dwdx, dwdy, dwdz);
               P[n].Bfieldy = P[n].Bfieldy + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(1, ax, ay, az, dwdx, dwdy, dwdz);
               P[n].Bfieldz = P[n].Bfieldz + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(2, ax, ay, az, dwdx, dwdy, dwdz); 
               }

             nh_typ = nh_typ + P[n].nh;
             bfield = pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);
             //bfield =  P[n].Bfieldz;
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
             //if(fabs(P[n].error) > 1.0) printf("P[n].nh = %lg, P[n].nh_test = %lg, P[n].error = %lg\n", P[n].Rho, P[n].nh_test, P[n].error);
             //if(j == snapmax) printf("P[n].nh = %lg, P[n].nh_test = %lg, P[n].error = %lg\n", P[n].Rho, P[n].nh_test, P[n].error);
             printf("n = %d x = %lg, y = %lg, z = %lg, bfieldx = %lg, bfieldy = %lg, bfieldz = %lg  ax = %lg, ay = %lg az = %lg, fcor = %lg\n", n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].Afieldx, P[n].Afieldy, P[n].Afieldz, P[n].fcor);

            }
       }

  printf("arrnum = %d\n", arrnum);

  nh_typ = nh_typ / num_typ;
  b_typ  = b_typ / num_typ;
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

  if(j == snapmin) 
    {
    n0 = nh_typ;
    b0 = b_typ;
    a0 = a_typ;
    a0_x = ax_typ;
    a0_y = ay_typ; 
    }
  n_anal = n0 * pow(delfac, -pfac); 
  b_anal = b0 * pow(n_anal/n0, 2./3.);

  if(which_sim == 1) a_anal = a0 * pow(n_anal/n0, 1./3.);

  //represents change in A Y-COMPONENT for shear flow where del_fac = 0.1
  if(which_sim == 2)  a_anal = a0_y - delfac*a0_x*(j - snapmin); 

  printf("nh_typ = %lg, nh_anal = %lg, b_typ = %lg, b_anal = %lg, u_typ = %lg, u_norm = %lg a_typ = %lg, a_anal = %lg\n", nh_typ, n_anal, b_typ, b_anal, u_typ, u_norm, a_typ, a_anal);

  outfile2=fopen(output_fname2, "a");
  fprintf(outfile2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
  nh_typ, n_anal, b_typ, b_anal, a_typ, a_anal, u_typ, u_norm, bx_typ, by_typ, bz_typ, ax_typ, ay_typ, az_typ); 
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

         unit_conversion();
 
         //for(n = 1; n <= NumPart; n++)
         for(n = 0; n < Ngas; n++)
             {

             P[n].Vel[0] = P[n].Vel[0] - vx;
             P[n].Vel[1] = P[n].Vel[1] - vy;
             P[n].Vel[2] = P[n].Vel[2] - vz;

             x = P[n].Pos[0];
             y = P[n].Pos[1];
             z = P[n].Pos[2];

             P[n].error = 2. * (P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test); 

             dis = pow(((P[n].Pos[0]-sinkposx)*(P[n].Pos[0]-sinkposx) + (P[n].Pos[1]-sinkposy)*(P[n].Pos[1]-sinkposy) + (P[n].Pos[2]-sinkposz)*(P[n].Pos[2]-sinkposz)), 0.5);
             //dis = pow(((P[n].Pos[0]-xmax)*(P[n].Pos[0]-xmax) + (P[n].Pos[1]-ymax)*(P[n].Pos[1]-ymax) + (P[n].Pos[2]-zmax)*(P[n].Pos[2]-zmax)), 0.5);
             if(dis > dismax)
                dismax = dis;
             dis=dis*1.e3*Time/(hubble_param);
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


    random = rand();
    int to_print = 0;
    if(random*(1.e0/RAND_MAX)< 0.1e0) to_print++;
    if(fabs(P[n].Bfieldx) > 0 || fabs(P[n].Bfieldy) > 0 || fabs(P[n].Bfieldz) > 0) to_print++;

    //if((random*(1.e0/RAND_MAX)< 0.05e0 && P[n].nh > 3.e-2) || P[n].nh > 2.e1 /*&& P[n].sink < 0*/ /* ||  P[n].sink > 0.5*/ ) {
    if(j == snapmax /* && to_print > 0*/) 
      {
      fprintf(outfile,"%15.13g %12d %15.13g %15.13g %15.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g\n",
      P[n].error, P[n].Id,x,y,z,P[n].Temp,P[n].nh, P[n].Rho, 1.0, 1.0, 5./3., disAU, vrad, vrot, P[n].hsm*Time/0.7, P[n].Mass, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz, P[n].nh_test);
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
  HubbleParam= hubble_param;


  for(i=0; i<NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
 
       MeanWeight=1.2195;

          MeanWeight=MeanWeight*PROTONMASS;

	  //MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  gamma= 5.0/3.0;

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

printf("U = %lg\n", P[1000].U);

	  SKIP;
	  for(n=0, pc_sph=pc; n<header1.npart[0];n++)
	    {
	      fread(&P[pc_sph].Density, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

printf("Density = %lg\n", P[1000].Density);

         SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].hsm, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

printf("hsm = %lg\n", P[1000].hsm);
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

printf("HeIII = %lg\n", P[1000].HeIII);

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].gam, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

printf("gam = %lg\n", P[1000].gam);

          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
              fread(&P[pc_sph].sink, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;
printf("sink = %lg\n", P[1000].sink);
*/
	}

      fclose(fd);
    }


  Time= header1.time;
  zred= header1.redshift;

  Time = 1.0;

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
    dkernel = (8./PI/pow(hsm,3)) * ( -12*pow(ratio,1)*(1./hsm)*drdu + 18*pow(ratio,2)*(1./hsm)*drdu );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = (8./PI/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(-1./hsm)*drdu;
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
    dkernel = -3.*(8./PI/pow(hsm,4)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3))
             + (8./PI/pow(hsm,3)) * ( - 12.*pow(ratio,1)*(-rad/hsm/hsm) +  18*pow(ratio,2)*(-rad/hsm/hsm) );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = -3.*(8./PI/pow(hsm,4)) * 2.*pow(1. - ratio, 3) + (8./PI/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(rad/hsm/hsm);
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


  











