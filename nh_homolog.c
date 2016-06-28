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
#define unit_length 3.085678e21
#define unit_mass 1.989e43
#define unit_energy 1.989e53
#define hubble_param 0.7

#define readB 0
#define convert_vpot 0
#define BFF 1
#define which_sim 3

int load_snapshot(char *fname, int files);
int allocate_memory(void);

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
  //double  Vel[3];
  double  Mass;
  int Type;
  int Id;

  double nh, Rho, Hsml;
  //double H2I, HII, DII, HDI, HeII, HeIII, gam,
  double sink;
  double Bfieldx, Bfieldy, Bfieldz, Bmag;
#if(readB)
  double Density;
  double Afieldx, Afieldy, Afieldz, fcor;
#endif
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
    int N_grid = 200;
    int N_begin = 0;  //which_sim==3  
    int N_end = 2;

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

    double min = 1.e-9;
    double max = 1.e-2;

    double abs_min = min;
    double abs_max = max;
    double time_in_Myr = 0.0;
    double grid1[N_grid][N_grid];
    double grid2[N_grid][N_grid];
    double grid3[1][N_grid];
    char dir[500];
    char buf[500], endtype[400];
    FILE *infile;
    FILE *outfile;
    int files;
    int N_gas;
    //double hsml_factor=3.0;
    double hsml_factor=1.7;
    //double hsml_factor=1.0;
    double center_i, center_j, W_x, weight, h, x2, x;
    double nh, nhmax, xmax, ymax, zmax, vx, vy, vz, xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;

    files=1;

    double IDmass=0.0;
    int IDsink, rotate=0;

    //sprintf(dir, "/work/00863/minerva/homolog");
    sprintf(dir, "/work/00863/minerva/homologR2");
    sprintf(endtype, ".nh");
    if(rotate == 1) sprintf(endtype, ".nhxzxz");

    for(m = N_begin; m < N_end+1; m=m+1)
      {
        if(m < 10)
          {
            sprintf(buf, "%s_000%d", dir, m);

            printf("reading 000%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_000%d%s", dir, m, endtype);

            printf("processing 000%d...\n", m);
          }
        else if(m > 9 && m < 100)
          {
            sprintf(buf, "%s_00%d", dir, m);

            printf("reading 00%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_00%d%s", dir, m, endtype);

            printf("processing 00%d...\n", m);
          }
        else if(m > 99 && m < 1000)
          {
            sprintf(buf, "%s_0%d", dir, m);

            printf("reading 0%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_0%d%s", dir, m, endtype);

            printf("processing 0%d...\n", m);
          }
        else
          {
            sprintf(buf, "%s_%d", dir, m);

            printf("reading %d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_%d%s", dir, m, endtype);

            printf("processing %d...\n", m);
          }

//Read in B-field info!///////////////////////////////////////////////////////////

  char path[200], pathout[200], path2[200], basename2[200], input_fname[200], input_fname2[200], output_fname[200], output_fname2[200], basename[200], basenameout[200];
  FILE *outfile, *outfile2, *infile;
  int snapshot_number, Ngas;
  snapshot_number = m;
  Ngas = N_gas;

  sprintf(path2, "/work/00863/minerva/orion/");
  sprintf(basename2, "bin_zoom10");

  if(which_sim == 0)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
   sprintf(basenameout, "snapbin_zoom10_new_cut");
   }
  if(which_sim == 3)
   {
   sprintf(path, "/work/00863/minerva");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "/bin_zoom_cut/bin_zoom10_new_cut");
   sprintf(basenameout, "snapbin_zoom10_new_cut");
   }
  if(which_sim == 4)
   {
   sprintf(path, "/nobackup/astacy");
   sprintf(path2, "/work/00863/minerva/orion/bin_ideal");
   sprintf(pathout, "/work/00863/minerva");
   sprintf(basename, "bin_HR10_ideal");
   sprintf(basename2,"bin_HR10_ideal");
   sprintf(basenameout, "snapbin_HR10_ideal");
   }

  if(convert_vpot == 1) sprintf(path2, "/work/00863/minerva/orion/bfield_comp_vpot_fwards");
  //printf(path2, "/work/00863/minerva/orion/bin_zoom/mid_tstep");

#if(readB)
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
#endif 
//////////////////////////////////////////////////////////////////////////////////
        width =  header.BoxSize;

        //slice= width;
        slice= (1./20)*width;

        IDmass = 0.0;
        IDsink = 0;

        double xtemp, ytemp, ztemp;

        for(n = 0; n < N_gas; n++)
          {
            //P[n].nh = P[n].Rho*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0)*X/m_H;
            //use g/cm^3 here instead of cm^-3
            P[n].nh = P[n].Rho*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0);
            if (n%10000 ==0) printf("%lg\n", P[n].Rho);

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

/*
         center_x=xCOM/ncount_doub;
         center_y=yCOM/ncount_doub;
         center_z=zCOM/ncount_doub;
*/

         center_x = header.BoxSize/2.0;
         center_y = header.BoxSize/2.0;
         center_z = header.BoxSize/2.0;

         printf("nhmax = %lg, center_x = %lg, center_y = %lg, center_z = %lg\n", nhmax, center_x, center_y, center_z);


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

        double Bmag_anal, Btot=0, ntot, Bmean=0, Brms=0, Bmax = 0;
        ntot = 0;
        for(n=0;n<Ngas;n++)
          {
          P[n].nh = P[n].nh / 1.67e-24;
          P[n].Bmag = pow(P[n].Bfieldx*P[n].Bfieldx + P[n].Bfieldy*P[n].Bfieldy + P[n].Bfieldz*P[n].Bfieldz, 0.5);

          if(P[n].convert == 1 && P[n].Bmag < 1.e-5) {Btot = Btot + P[n].Bmag; ntot++;}
          //if(P[n].convert == 1) printf("n = %d, Bmag = %lg\n", n, P[n].Bmag);

          if(P[n].Bmag > Bmax && P[n].Bmag < 1.e-5) Bmax = P[n].Bmag;
          }            
          
         Bmean = Btot / ntot;
         printf("Bmean = %lg Bmax = %lg\n", Bmean, Bmax);

//find a fitting function to B versus n

        double Afit, Bfit, lnx=0, lny=0, lnx2=0, lnxlny=0, nfit=0, nskip=0, hifac;
        for(n=0;n<Ngas;n++)
          {
          //if(P[n].Bmag > 1.e-19)
          if(P[n].convert == 1)
            {
            nfit++;
            lnx = lnx + log(P[n].nh);
            lnx2 = lnx2 + log(P[n].nh)*log(P[n].nh);
            lny = lny + log(P[n].Bmag);
            lnxlny = lnxlny + log(P[n].nh)*log(P[n].Bmag);
            }
          }
        Bfit = (nfit*lnxlny - lnx*lny)/(nfit*lnx2 - lnx*lnx);
        Afit = (lny - Bfit*lnx)/nfit;
        Afit = pow(2.71818, Afit);

        printf("nfit = %lg, lnx = %lg, lnx2 = %lg, lny = %lg \n", nfit, lnx, lnx2, lny);
        printf("Afit = %lg, Bfit = %lg \n", Afit, Bfit);            

        for(n = 0; n < N_gas; n++)
          {
            if(P[n].Bmag > 2.e-20)
              if(P[n].Bmag > 1.e2 * Afit * pow(P[n].nh, Bfit) || P[n].Bmag < 1.e-3 * Afit * pow(P[n].nh, Bfit)) 
                 {
                 nskip ++;
                 hifac = P[n].Bmag / (Afit * pow(P[n].nh, Bfit));
                 //P[n].Bfieldx = P[n].Bfieldx / hifac; 
                 //P[n].Bfieldy = P[n].Bfieldy / hifac; 
                 //P[n].Bfieldz = P[n].Bfieldz / hifac; 
                 //P[n].Bmag = P[n].Bmag / hifac; 
                 //printf("nskip = %lg\n", nskip); 
                 //continue;
                 }

            if(P[n].sink > -1 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < slice/2.0)
              {

//P[n].Bmag = Afit * pow(P[n].nh, Bfit);

//h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0);
h = hsml_factor * P[n].Hsml; 

//weight = P[n].nh * P[n].nh;
//weight = pow(h, -3);
weight = 1.;

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

//x2 = ((P[n].Pos[0] - center_i) * (P[n].Pos[0] - center_i) + (P[n].Pos[1] - center_j) * (P[n].Pos[1] - center_j) + (P[n].Pos[2] - center_z) * (P[n].Pos[2] - center_z) );
//x2 = sqrt(x2) / h;

if(x2 <= 1.0)
  {
    x = sqrt(x2);
    //x = x2;

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

    grid1[i][j] += weight * P[n].Rho * W_x;

    n_grid1[i][j] += weight * W_x;
   }


               //                 n_grid1[i][j] += 1;

                 //               grid1[i][j] += P[n].nh;
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
            if(double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[0]-center_x)) < slice/2.0)
              {

h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0); 
//weight = P[n].nh * P[n].nh;
weight = pow(h, -3);

                i = 0;
                j=0;

                i_min = int((P[n].Pos[1]-h-center_y+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[1]+h-center_y+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[2]-h-center_z+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[2]+h-center_z+width/2.0)/width*N_grid);

                if(j_max > N_grid-1)
                  {
                    j_max = N_grid-1;
                  }

                do
                  {
                    if(i >= i_min && i <= i_max)
                      {
                        flag_i = 1;

center_i = center_y - width / 2.0 + (i + 0.5) * width / (double) N_grid;

                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

center_j = center_z - width / 2.0 + (j + 0.5) * width / (double) N_grid;

x2 = ((P[n].Pos[1] - center_i) * (P[n].Pos[1] - center_i) + (P[n].Pos[2] -
center_j) * (P[n].Pos[2] - center_j)) / h / h;

if(x2 <= 1.0)
  {
    x = sqrt(x2);

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

    grid2[i][j] += weight * P[n].Rho * W_x;

    n_grid2[i][j] += weight * W_x;
  }


                                //n_grid2[i][j] += 1;

                                //grid2[i][j] += P[n].nh;
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
  int    l,j,k,p,dummy,ntot_withmasses;
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
	      //fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
              for(p=0;p<3;p++) fread(&P[pc_new].dummy, sizeof(double), 1, fd);
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


          SKIP;
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

#if(BFF)
          SKIP;
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
            {
              fread(&P[pc_sph].Bfieldx, sizeof(double), 1, fd);
              fread(&P[pc_sph].Bfieldy, sizeof(double), 1, fd);
              fread(&P[pc_sph].Bfieldz, sizeof(double), 1, fd);
              pc_sph++;
            }
          SKIP;

          printf("Bfieldx = %lg\n", P[1000].Bfieldx);

#endif
	  
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


