#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXREF 20
#define PI 3.1415927
#define m_H 1.6726e-24
#define k_B 1.3806e-16
#define Hubble 3.2407789e-18
#define G 6.672e-8
#define X 0.76
#define yr 3.1536e7
#define gamma 5.0/3.0
#define unit_length 1.0
//#define unit_mass 1.67e-24
#define unit_mass 1.0
#define unit_energy 1.989e53
#define hubble_param 1.0

#define readB 0
#define convert_vpot 0
#define which_sim 6

#define do_vel 1

int load_snapshot(char *fname, int files);
int allocate_memory(void);
int reordering(void);

//functions for calculating bfield from vpot
double calc_kernel_spline_vpot(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double cosmo_fac);
double del_kernel(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double del_rho(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double curl(int dir, double x1, double y1, double z1, double x2, double y2, double z2);

struct io_header
  {
    int npart[MAXREF];
    double mass[MAXREF];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    int npartTotal[MAXREF];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
  } header;

int idnum=1;
struct particle_data
  {
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int Type;
  int Id;

  //double  U, Temp, 
  double nh, Density, Rho, Hsml;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, gam;
  double sink;
  //double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;

  double nh_test, Bfieldx, Bfieldy, Bfieldz, Bmag, error;
  double Afieldx, Afieldy, Afieldz, fcor;

  double dummy;
  int convert; 
  } *P;
//int *Id;

struct plist_data
{
  int Id;
  double disx, disy, disz;
}*plist;

int N_part;
int N_gas;

double  Time, zred;


int main(int argc, char **argv)
  {
    int i = 0;
    int j = 0;
    int n = 0;
    int m = 0;
    int i_min = 0;
    int i_max = 0;
    int j_min = 0;
    int j_max= 0;
    int flag_i = 0;
    int flag_j = 0;
    int N_dm = 0;
    int N_grid = 256;

    int N_begin, N_end;

    if(which_sim == 1)
    {
    N_begin = 2267;
    N_end = 2267;
    }

    if(which_sim == 2)
    {
    N_begin = 4713;
    N_end = 4713;
    }

    if(which_sim == 3)
    {
    //N_begin = 3300;
    //N_end = 3300;    
    N_begin = 7098;
    N_end = 7098;
    }

    if(which_sim == 4)
    {
    N_begin = 50;
    N_end = 50;
    }

    if(which_sim == 6)
    {
    N_begin = 85;
    N_end = 85;
    }

//    int n_grid1[N_grid][N_grid];
    double n_grid1[N_grid][N_grid];
//    int n_grid2[N_grid][N_grid];
    double n_grid2[N_grid][N_grid];

    double center_x = 500.0;  
    double center_y = 500.0;
    double center_z = 500.0;

    double width = 1000.;
    double slice = 1000.;
    double entries = 0.0;
    double temp_pos = 0.0;
    double min = 1.e-1;     //hires density range  
    double max = 5.e-1;
    //double min = 0.2;  
    //double max = 0.3;
    if(do_vel == 1) {min = 1.e-1; max = 1.5e0;}
    double abs_min = min;
    double abs_max = max;
    double time_in_Myr = 0.0;
    double grid1[N_grid][N_grid];
    double grid2[N_grid][N_grid];
    double grid3[1][N_grid];
    char dir[500], dir_out[500];
    char buf[500], endtype[400];
    FILE *infile;
    FILE *outfile;
    int files;
    int N_gas;
    double hsml_factor=1.e0;
    double center_i, center_j, W_x, weight, h, x2, x;
    double nh, nhmax, xmax, ymax, zmax, vx, vy, vz, xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;

    files=1;

    double IDmass=0.0;
    int IDsink, rotate=0;

    sprintf(dir, "/scratch/00863/minerva/ot");
    sprintf(dir_out, "/scratch/00863/minerva/ot");

    if(which_sim == 0) 
      {
      sprintf(dir, "/work/00863/minerva/ot_lowres/ot");
      sprintf(dir_out, "/work/00863/minerva/ot");
      }
    if(which_sim == 1) 
      {
      sprintf(dir, "/scratch/00863/minerva/ot_twodim/ot"); 
      sprintf(dir_out, "/work/00863/minerva/ot");
      }
    if(which_sim == 2) 
      {
      sprintf(dir, "/scratch/00863/minerva/ot_twodim_midres/ot");
      sprintf(dir_out, "/work/00863/minerva/ot");
      }
    if(which_sim == 3) 
      {
      sprintf(dir, "/scratch/00863/minerva/ot_twodim_highres/ot");
      sprintf(dir_out, "/work/00863/minerva/ot");
      }

    sprintf(endtype, ".nh");
    if(rotate == 1) sprintf(endtype, ".nhxz");

   if(do_vel == 1)
   {
    sprintf(endtype, ".vel");
    if(rotate == 1) sprintf(endtype, ".velxz");
   }

    for(m = N_begin; m < N_end+1; m=m+10)
      {
        if(m < 10)
          {
            sprintf(buf, "%s_000%d", dir, m);

            printf("reading 000%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_000%d%s", dir_out, m, endtype);

            printf("processing 000%d...\n", m);
          }
        else if(m > 9 && m < 100)
          {
            sprintf(buf, "%s_00%d", dir, m);

            printf("reading 00%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_00%d%s", dir_out, m, endtype);

            printf("processing 0%d...\n", m);
          }
         else if(m > 99 && m < 1000)
          {
            sprintf(buf, "%s_0%d", dir, m);

            printf("reading 0%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_0%d%s", dir_out, m, endtype);

            printf("processing 0%d...\n", m);
          }
        else
          {
            sprintf(buf, "%s_%d", dir, m);

            printf("reading %d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_%d%s", dir_out, m, endtype);

            printf("processing %d...\n", m);
          }

//Read in B-field info!///////////////////////////////////////////////////////////

  char path[200], pathout[200], path2[200], basename2[200], input_fname[200], input_fname2[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  FILE *outfile, *outfile2, *infile;
  int snapshot_number, Ngas;
  snapshot_number = m;
  Ngas = N_gas;

  sprintf(path2, "/work/00863/minerva/orion/");
  if(which_sim == 0) sprintf(path2, "/work/00863/minerva/orion/ot_lowres");
  if(which_sim == 1) sprintf(path2, "/work/00863/minerva/orion/ot_twodim");
  if(which_sim == 2) sprintf(path2, "/work/00863/minerva/orion/ot_twodim_midres");
  if(which_sim == 3) sprintf(path2, "/work/00863/minerva/orion/ot_twodim_highres");
  sprintf(basename2, "ot");

  if(which_sim == 4)
    {
    sprintf(path2, "/work/00863/minerva/orion/vort");
    sprintf(basename2, "vort");
    }

  if(which_sim == 0)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
   sprintf(basenameout, "snapbin_zoom10_new_cut");
   }

  if(convert_vpot == 1) sprintf(path2, "/work/00863/minerva/orion/bfield_comp_vpot_fwards");

if(readB == 1)
{
    sprintf(input_fname2, "%s/%s_bfield_%04d", path2, basename2, snapshot_number);
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }

    printf("reading bfield file %s!\n", input_fname2);

  if(convert_vpot != 1)
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
  if(convert_vpot == 1)
    {
    for(n=0;n<Ngas;n++)
        P[n].Bfieldx = P[n].Bfieldy = P[n].Bfieldz = P[n].Density = P[n].convert = 0;
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
}
      
//////////////////////////////////////////////////////////////////////////////////
        double time_fac = 1.0;
        double width_cm = 1.e8;
 
        width =  width_cm*header.HubbleParam/(time_fac)/1.e8;

        slice= width;
        slice= (1./32)*width;

        IDmass = 0.0;
        IDsink = 0;

//test our vortex algorithm
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   double x_center, y_center, z_center;
   double x, y, rad, velx, vely, phi, phi_new;
   double xnew, ynew, rnew, half_box;

   half_box = header.BoxSize / 2.;
   x_center = header.BoxSize / 2.;
   y_center = header.BoxSize / 2.;
   z_center = header.BoxSize / 2.;

   double r1 = .45*half_box, r2 = 0.9*half_box, del_phi = 0.01*2.*PI;

   for(i = 0; i < N_gas; i++)
      {

      x = P[i].Pos[0] - x_center;
      y = P[i].Pos[1] - y_center;
      rad = pow(x*x + y*y,  0.5);

      phi = atan2(y,x);

      phi_new = phi;

      if(rad < r1)
         phi_new = phi + del_phi;

      xnew = rad * cos(phi_new);
      ynew = rad * sin(phi_new);
      rnew =  pow(xnew*xnew + ynew*ynew,  0.5);

/*
      if(fabs(rad-rnew) > 0.001*rnew)
        printf("RAD CHANGE id = %d, rho = %lg, phi = %lg, phi_new = %lg, x = %lg, xnew = %lg, y = %lg, ynew = %lg, r = %lg, rnew = %lg\n", P[i].Id, P[i].Density, phi, phi_new, x, xnew, y, ynew, rad, rnew);

      if(fabs(x-xnew) > 0.1*half_box)
        printf("XNEW CHANGE id = %d, rho = %lg, phi = %lg, phi_new = %lg, x = %lg, xnew = %lg, y = %lg, ynew = %lg, r = %lg, rnew = %lg\n", P[i].Id, P[i].Density, phi, phi_new, x, xnew, y, ynew, rad, rnew);

      if(fabs(y-ynew) > 0.1*half_box)
        printf("YNEW CHANGE id = %d, rho = %lg, phi = %lg, phi_new = %lg, x = %lg, xnew = %lg, y = %lg, ynew = %lg, r = %lg, rnew = %lg\n", P[i].Id, P[i].Density, phi, phi_new, x, xnew, y, ynew, rad, rnew);

      if(fabs(xnew) < 0.01*half_box && fabs(ynew) < 0.01*half_box)
        printf("XNEW CHANGE id = %d, rho = %lg, phi = %lg, phi_new = %lg, x = %lg, xnew = %lg, y = %lg, ynew = %lg, r = %lg, rnew = %lg\n", P[i].Id, P[i].Density, phi, phi_new, x, xnew, y, ynew, rad, rnew);
*/
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        double xtemp, ytemp, ztemp;

        for(n = 0; n < N_gas; n++)
          {
            //use g/cm^3 here instead of cm^-3
            P[n].nh = P[n].Rho*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(time_fac,3.0);
            if(P[n].nh > 0.35) { /*printf("high nh = %lg\n", P[n].nh);*/  P[n].sink = 1;}
            if (n%10000 ==0) printf("%lg\n", P[n].nh);

            xtemp = P[n].Pos[0];
            ytemp = P[n].Pos[1];
            ztemp = P[n].Pos[2];

            if(rotate == 1)
            {
            P[n].Pos[0] = ytemp;
            P[n].Pos[1] = ztemp;
            P[n].Pos[2] = xtemp;
            } 
         }
 
         xCOM = yCOM = zCOM = vxCOM = vyCOM = vzCOM = ncount_doub = nhmax = 0;
         for(n=0;n < N_gas; n++)
            {
             nh = P[n].nh;
             if(nh > nhmax)
               {
               nhmax = nh;
               xmax=P[n].Pos[0];
               ymax=P[n].Pos[1];
               zmax=P[n].Pos[2];
               vx = P[n].Vel[0];
               vy = P[n].Vel[1];
               vz = P[n].Vel[2];
               }
            }

         for(n=0;n < N_gas; n++)
            {
            nh = P[n].nh;
            if(nh > nhmax/1.5)
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

         vx = vxCOM/ncount_doub;
         vy = vyCOM/ncount_doub;
         vz = vzCOM/ncount_doub;

         center_x=xCOM/ncount_doub;
         center_y=yCOM/ncount_doub;
         center_z=zCOM/ncount_doub;

         center_x = center_y = center_z = header.BoxSize/2.0;
         center_z = 0;
         printf("center_x = %lg, center_y = %lg, center_z = %lg, width = %lg, nhmax = %lg, hubble = %lg\n", center_x, center_y, center_z, width, nhmax, header.HubbleParam);

////////////////////////////////////////////////////////////////////////////////////////////////////////
// convert vector potential to b-field
#if (convert_vpot)

  double disx, disy, disz;

  for(n=0;n<Ngas;n++)
     {
     disx = ((P[n].Pos[0]-center_x))*1.e3*Time/(hubble_param);
     disy = ((P[n].Pos[1]-center_y))*1.e3*Time/(hubble_param);
     disz = ((P[n].Pos[2]-center_z))*1.e3*Time/(hubble_param);

     if(fabs(disx) < width_pc && fabs(disy) < width_pc && fabs(disz) < width_pc)
       {
       P[n].convert = 1;
       idnum++;
       }
     }

    printf("final idnum = %d, Ngas = %d\n", idnum, Ngas);

    if(!(plist=(struct plist_data *) malloc(idnum*sizeof(struct plist_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }
    printf("plist memory allocated\n");

    n=0;
    for(i = 0; i < Ngas; i++)
      {
      if(P[i].convert == 1)
        {
        plist[n].Id = P[i].Id;
        n++;
        }
      if(i % 1000 == 0 && n > 1)
        printf("i = %d, rho = %lg, Id = %d, id = %d\n", i, P[i].Rho, P[i].Id, plist[n-1].Id);
      }

  double dwdx, dwdy, dwdz, ax, ay, az, rcheck;
  double bfac = pow(Time / hubble_param, -1);
  double cosmo_fac = Time/hubble_param;

  for(n=0;n<Ngas;n++)
    {
    P[n].Pos[0] = P[n].Pos[0]*cosmo_fac;
    P[n].Pos[1] = P[n].Pos[1]*cosmo_fac;
    P[n].Pos[2] = P[n].Pos[2]*cosmo_fac;
    P[n].Hsml = P[n].Hsml*cosmo_fac;
    }

  int k;
  double nh_min = 1.e-21;
  struct particle_data psave, psource;
  printf("start reordering\n");
  for(n=0;n<Ngas;n++)
      if(P[n].nh > nh_min)
      for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
             {
             psource = P[n];
             psave   = P[k];
             P[n] = psave;
             P[k] = psource;
             }
  printf("finished reordering\n");

  for(n=0;n<Ngas;n++)
     {
     if(P[n].convert == 1)
       {
       for(i=0;i<idnum;i++)
         {
         if(i == n) continue;
         rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
         if(rcheck > P[n].Hsml) continue;
         P[n].Density = P[n].Density + P[i].Mass * calc_kernel_spline_vpot(n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].Hsml, 0);
          }
        if(n%1000 == 0)  printf("n = %d, Density = %lg\n", n, P[n].Density);
        }
     }

  printf("finished density calculation\n");

  for(n=0;n<Ngas;n++)
     {
     if(P[n].convert == 1)
        {
        for(i=0;i<idnum;i++)
           {
           if(i == n) continue;
           rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
           if(rcheck > P[n].Hsml) continue;
           P[n].fcor = P[n].fcor + (P[n].Hsml /  3. / P[n].Density) * P[i].Mass * del_rho(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].Hsml, 0);
           }
         if(n%1000 == 0)  printf("n = %d, fcor = %lg\n", n, P[n].fcor);
         }   //end if(convert_vpot == 1) statement
      }  //end for(n=0;n < Ngas; n++) loop

  printf("finished fcor calculation\n");

  for(n=0;n<Ngas;n++)
     {
     if(P[n].convert == 1)
       {
        for(i=0;i<idnum;i++)
           {
           if(i == n) continue;
           rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
           if(rcheck > P[n].Hsml) continue;
           dwdx = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].Hsml, P[n].Pos[0] - P[i].Pos[0]);
           dwdy = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].Hsml, P[n].Pos[1] - P[i].Pos[1]);
           dwdz = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].Hsml, P[n].Pos[2] - P[i].Pos[2]);
           ax = P[n].Afieldx - P[i].Afieldx;
           ay = P[n].Afieldy - P[i].Afieldy;
           az = P[n].Afieldz - P[i].Afieldz;
           P[n].Bfieldx = P[n].Bfieldx + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(0, ax, ay, az, dwdx, dwdy, dwdz);
           P[n].Bfieldy = P[n].Bfieldy + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(1, ax, ay, az, dwdx, dwdy, dwdz);
           P[n].Bfieldz = P[n].Bfieldz + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(2, ax, ay, az, dwdx, dwdy, dwdz);
           }
         if(n%2000 == 0)
        printf("n = %d, dwdx = %lg, dwdy = %lg, dwdz = %lg, ax = %lg, ay = %lg, az = %lg\n", n, dwdx, dwdy, dwdz, ax, ay, az);
         }   //end if(convert_vpot == 1) statement
      }  //end for(n=0;n < Ngas; n++) loop

  for(n=0;n<Ngas;n++)
    {
    P[n].Pos[0] = P[n].Pos[0]/cosmo_fac;
    P[n].Pos[1] = P[n].Pos[1]/cosmo_fac;
    P[n].Pos[2] = P[n].Pos[2]/cosmo_fac;
    P[n].Hsml = P[n].Hsml/cosmo_fac;
    }

#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////// 

        time_in_Myr = 2.0/3.0/Hubble/header.HubbleParam/sqrt(header.Omega0)*pow(header.time,3.0/2.0)/1.0e6/yr;

        flag_i = 0;
        flag_j = 0;

        for(i = 0; i < N_grid; i++)
          {
            for(j = 0; j < N_grid; j++)
              {
                n_grid1[i][j] = 0;

                grid1[i][j] = 0.0;
              }
          }

        if(readB == 1)
          for(n=0;n<Ngas;n++)
            P[n].Bmag = pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);

        double KE = 0;
        for(n=0;n<Ngas;n++)
          KE = KE + 0.5*P[n].Mass*(P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1] + P[n].Vel[2]*P[n].Vel[2]); 
        printf("KE total = %lg, mass = %lg\n", KE, P[10].Mass);

        int n_include;

        for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > -1 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < slice/2.0)
              {

n_include++;
//if(n%1000 == 0) printf("n_include = %d, nh_include = %lg, hsm_include = %lg\n", n_include, P[n].nh, P[n].Hsml);

//h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0);
//weight = P[n].nh * P[n].nh;
h = hsml_factor * P[n].Hsml;
weight = pow(h, -3);

if(P[n].sink > 0)
 //h = hsml_factor * 3.57101e-07;
 h = hsml_factor * 2.5696e-08;

                i = 0;
                j=0;

                i_min = int((P[n].Pos[0]-h-center_x+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[0]+h-center_x+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[1]-h-center_y+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[1]+h-center_y+width/2.0)/width*N_grid);

                if(j_max > N_grid-1)
                  {
                    j_max = N_grid-1;
                  }

                do
                  {
                    if(i >= i_min && i <= i_max)
                      {
                        flag_i = 1;

center_i = center_x - width / 2.0 + (i + 0.5) * width / (double) N_grid;
                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

center_j = center_y - width / 2.0 + (j + 0.5) * width / (double) N_grid;

x2 = ((P[n].Pos[0] - center_i) * (P[n].Pos[0] - center_i) + (P[n].Pos[1] - center_j) * (P[n].Pos[1] - center_j)) / h / h;

if(x2 <= 1.0)
  {
    x = sqrt(x2);

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

   if(do_vel == 1) grid1[i][j] += weight * pow(P[n].Vel[0]*P[n].Vel[0] + P[n].Vel[1]*P[n].Vel[1], 0.5) * W_x;
   else   grid1[i][j] += weight * P[n].nh * W_x;

    n_grid1[i][j] += weight * W_x;
   }


                              }
                            else if(j > j_max)
                              {
                                flag_j = 2;
                              }
                            else
                              {
                                flag_j = 0;
                              }

                             j++;
                           }
                         while((flag_j == 0 || flag_j == 1) && j < N_grid);

                         j = 0;
                       }
                     else if(i > i_max)
                       {
                         flag_i = 2;
                       }
                     else
                       {
                         flag_i = 0;
                       }

                     i++;
                   }
                 while((flag_i == 0 || flag_i == 1) && i < N_grid);

                 i = 0;
               }
           }

        for(i = 0; i < N_grid; i++)
          {
            for(j = 0; j < N_grid; j++)
              {
                if(n_grid1[i][j] > 0)
                  {
                    grid1[i][j] /= double(n_grid1[i][j]);
                  }

                if(grid1[i][j] < min)
                  {
                    grid1[i][j] = min;
                  }

                if(grid1[i][j] > max)
                  {
                    grid1[i][j] = max;
                  }

                grid1[i][j] = log10(grid1[i][j]);
              }
           }

        grid1[0][0] = log10(abs_min);
        grid1[N_grid-1][N_grid-1] = log10(abs_max);

        for(i = 0; i < N_grid; i++)
          for(j = 0; j < N_grid; j++)
             //if(grid1[i][j] > log10(abs_min))  printf("grid[%d][%d] = %lg\n", i, j, grid1[i][j]); 
             if(i%100 == 0 and j%100 == 0)  printf("grid[%d][%d] = %lg, n_grid = %lg\n", i, j, grid1[i][j], n_grid1[i][j]);

        flag_i = 0;
        flag_j = 0;

        for(i = 0; i < N_grid; i++)
          {
            for(j = 0; j < N_grid; j++)
              {
                n_grid2[i][j] = 0;

                grid2[i][j] = 0.0;
              }
          }

        for(n = 0; n < N_gas; n++)
          {
            if(double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < slice/2.0)
              {

//h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0); 
hsml_factor * P[n].Hsml;
weight = pow(h, -3);

                i = 0;
                j=0;

                i_min = int((P[n].Pos[0]-h-center_x+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[0]+h-center_x+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[1]-h-center_y+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[1]+h-center_y+width/2.0)/width*N_grid);

                if(j_max > N_grid-1)
                  {
                    j_max = N_grid-1;
                  }

                do
                  {
                    if(i >= i_min && i <= i_max)
                      {
                        flag_i = 1;

center_i = center_x - width / 2.0 + (i + 0.5) * width / (double) N_grid;

                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

center_j = center_y - width / 2.0 + (j + 0.5) * width / (double) N_grid;

x2 = ((P[n].Pos[0] - center_i) * (P[n].Pos[0] - center_i) + (P[n].Pos[1] - center_j) * (P[n].Pos[1] - center_j)) / h / h;

if(x2 <= 1.0)
  {
    x = sqrt(x2);

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

    grid2[i][j] += weight * pow(P[n].nh - pow(10., grid1[i][j]), 2) * W_x;

    n_grid2[i][j] += weight * W_x;
  }


                              }
                            else if(j > j_max)
                              {
                                flag_j = 2;
                              }
                            else
                              {
                                flag_j = 0;
                              }

                             j++;
                           }
                         while((flag_j == 0 || flag_j == 1) && j < N_grid);

                         j = 0;
                       }
                     else if(i > i_max)
                       {
                         flag_i = 2;
                       }
                     else
                       {
                         flag_i = 0;
                       }

                     i++;
                   }
                 while((flag_i == 0 || flag_i == 1) && i < N_grid);

                 i = 0;
               }
           }

        for(i = 0; i < N_grid; i++)
          {
            for(j = 0; j < N_grid; j++)
              {
                if(n_grid2[i][j] > 0)
                  {
                    grid2[i][j] /= double(n_grid2[i][j]);
                  }

                grid2[i][j] = pow(grid2[i][j],0.5);

                if(grid2[i][j] < min)
                  {
                    grid2[i][j] = min;
                  }

                if(grid2[i][j] > max)
                  {
                    grid2[i][j] = max;
                  }

                grid2[i][j] = log10(grid2[i][j]);
              }
           }

        grid2[0][0] = log10(abs_min);
        grid2[N_grid-1][N_grid-1] = log10(abs_max);

        for(i = 0; i < N_grid; i++)
          for(j = 0; j < N_grid; j++)
            if(i%100 == 0 and j%100 == 0)  printf("grid2[%d][%d] = %lg, n_grid2 = %lg\n", i, j, grid2[i][j], n_grid2[i][j]);

        for(j = 0; j < N_grid; j++)
          {
            grid3[0][j] = log10(min)+(2.0*j+1.0)*(log10(max)-log10(min))/double(N_grid)/2.0;
          }

        grid3[0][0] = log10(abs_min);

        outfile = fopen(buf, "w");

        fwrite(grid1,sizeof(grid1),1,outfile);
        fwrite(grid2,sizeof(grid2),1,outfile);
        fwrite(grid3,sizeof(grid3),1,outfile);
        fwrite(&header.redshift,sizeof(double),1,outfile);
        fwrite(&time_in_Myr,sizeof(double),1,outfile);

        entries = 0.0;

        for(n = 0; n < N_gas; n++)
          {
           if(P[n].sink > 0.5 && P[n].Mass > IDmass)
             {
             IDmass = P[n].Mass;
             IDsink = P[n].Id;
             }
          }


        for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > 0.5 /*Id[n]*/ /*P[n].Id == IDsink*/  && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                entries++;
              }
          }

        fwrite(&entries,sizeof(double),1,outfile);

        for(n = 0; n < N_gas; n++)
          {
            if(/*P[n].flag_split > 0.5*/  /*Id[n]*/ P[n].Id == IDsink  && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[0]-center_x+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }


        for(n = 0; n < N_gas; n++)
          {
            if(/*P[n].flag_split > 0.5*/  /*Id[n]*/ P[n].Id == 4333881  && P[n].sink > 0.5 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[0]-center_x+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

        for(n = 0; n < N_gas; n++)
          {
             if(P[n].sink > 0.5 && /*Id[n]*/  P[n].Id != IDsink  &&  P[n].Id != 4333881 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[0]-center_x+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }     


        for(n = 0; n < N_gas; n++)
          {
            if(/*P[n].flag_split > 0.5 &&*/ /*Id[n]*/ P[n].Id == IDsink  && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[1]-center_y+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

       for(n = 0; n < N_gas; n++)
          {
            if(/*P[n].flag_split > 0.5 &&*/ /*Id[n]*/ P[n].Id == 4333881  && P[n].sink > 0.5 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[1]-center_y+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

       for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > 0.5 && /*Id[n]*/ P[n].Id != IDsink  && P[n].Id != 4333881 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[1]-center_y+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }


        for(n = 0; n < N_gas; n++)
          {
            if(/*P[n].flag_split > 0.5 &&*/  /*Id[n]*/ P[n].Id == IDsink  &&  double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[2]-center_z+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

        for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > 0.5 &&  /*Id[n]*/ P[n].Id != IDsink  && P[n].Id != 4333881 &&  double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[2]-center_z+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

        fclose(outfile);

        free(P);
      }

    printf("done!\n");
  }

int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=(struct particle_data *) malloc(N_part*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
//  P--;   /* start with offset 1 */

 /* 
  if(!(Id=(int *) malloc(N_part*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   // start with offset 1 
*/
  printf("allocating memory...done\n");
}

int load_snapshot(char *fname, int files)
  {
  FILE *fd;
  char   buf[200];
  int    l,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc=0,pc_new=0,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
//#define SKIP1 fread(&dummy1, sizeof(dummy1), 1, fd);

  for(l=0, pc=0; l<files; l++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,l);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, N_part=0, ntot_withmasses=0; k<5; k++)
	    N_part+= header.npart[k];
	  N_gas= header.npart[0];
	}
      else
	{
	  for(k=0, N_part=0, ntot_withmasses=0; k<5; k++)
	    N_part+= header.npartTotal[k];
	  N_gas= header.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header.mass[k]==0)
	    ntot_withmasses+= header.npart[k];
	}

      if(l==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
      printf("N_gas= %6d \n",N_gas); 
      printf("pos= %lg \n", &P[pc_new].Pos[0]);
/*      printf("N_DM= %6d \n",N_part-N_gas); */

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header.npart[k];n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header.npart[k];n++)
	    {
	      //fread(&Id[pc_new], sizeof(int), 1, fd);
              fread(&P[pc_new].Id, sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;



      if(ntot_withmasses>0)
	SKIP;
      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      
	      if(header.mass[k]==0)
	      	fread(&P[pc_new].Mass, sizeof(double), 1, fd);
	      else
		P[pc_new].Mass= header.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
	SKIP;

      printf("line 320\n");
      printf("Npart= %6d \n",N_part);
      printf("M_B= %15.6e \n",header.mass[0]);
      printf("M_DM= %15.6e \n",header.mass[1]);
      
      printf("line 324\n");

      if(header.npart[0]>0)
	{
	  SKIP;
	  for(n=0, pc_sph=pc; n<header.npart[0];n++)
	    {
	      fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header.npart[0];n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

         SKIP;
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
            {
              fread(&P[pc_sph].Hsml, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;


/*          SKIP;
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
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
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
            {
              fread(&P[pc_sph].dummy, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

          SKIP;
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
            {
              fread(&P[pc_sph].sink, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;
*/
	  
	}

      fclose(fd);
  }

  Time= header.time;
  zred= header.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.4e \n",Time);
  printf("L= %6.2f \n",header.BoxSize);
  return(N_gas);
 }

double calc_kernel_spline_vpot(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double grid_size_half)
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
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].Hsml);
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

  for(i=0; i<N_part; i++)
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

  printf("space for particle ID freed\n");

  return(0);
}

