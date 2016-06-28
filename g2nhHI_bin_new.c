#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define float double

#define MAXREF 20
//define N_grid 700
#define pi 3.1415927
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

int load_snapshot(char *fname, int files);
int allocate_memory(void);

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

struct particle_data
  {
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int Type;
  int Id;

  double  U, Temp, nh, Density, Hsml;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
  double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy; 
  } *P;
//int *Id;

int N_part;
int N_gas;

double  Time, zred;


int main(int argc, char **argv)
  {
    int halo = 12;
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
    int N_grid = 500;
    //int N_begin = 5980;
    //int N_end =5980;
    int N_begin = 15781;
    int N_end =15781;
	  
//    int n_grid1[N_grid][N_grid];
    double n_grid1[N_grid][N_grid];
//    int n_grid2[N_grid][N_grid];
    double n_grid2[N_grid][N_grid];

    double center_x = 50.339339267;
    double center_y = 50.121980449;
    double center_z = 49.486589836;

    if(halo == 1)
      {
      printf("hihihihihi\n");
      center_x = 501.10838567;   //bin_HR1 center
      center_y = 502.63651214;
      center_z = 500.57892754;
      }
    else if(halo == 2)
      {
      center_x = 504.74726968;  //bin_HR2 center
      center_y = 502.54626581;
      center_z = 498.3606391;
      }
    else if(halo == 3)
      {
      center_x = 497.9768553;  //bin_HR3 center
      center_y = 478.97606043;
      center_z = 533.38845866;
      }
    else if(halo == 4)
      {
      center_x = 492.79147686;  //bin_HR4 center
      center_y = 481.6921907;
      center_z = 488.42287753;
      }
    else if(halo == 5)
      {
      center_x = 500.79219936;  //bin_HR5 center
      center_y = 499.63973988;
      center_z = 498.79510558;
      }
    else if(halo == 6)
      {
      center_x = 491.44465306;  //bin_HR6 center
      center_y = 524.86928186;
      center_z = 504.57007015;
      }
    else if(halo == 7)
	  {		           
	  center_x = 500.88455606;  //bin_HR7 center
	  center_y = 499.5879108;
	  center_z = 504.99181437;
	  } 	  
    else if(halo == 8)
  	  {		          
		center_x = 506.01033587;  //bin_HR8 center
		center_y = 501.56205197;
		center_z = 503.11039489;
	  } 	  	  
	else if(halo == 9)
	  {
		center_x = 496.96029045;  //bin_HR9 center
		center_y = 497.98100903;
		center_z = 514.16695849;
	  }
    else if(halo == 10)
	  {
		center_x = 512.65800039;  //bin_HR10 center
		center_y = 505.21837386;
		center_z = 497.39347151;
	  }
	  
    else
      {
      double center_x = 0;  //nothing
      double center_y = 0;
      double center_z = 0;
      }

    printf("halo = %d, x = %lg, y = %lg, z = %lg\n", halo, center_x, center_y, center_z);

    double width = 140.;
    double slice = 140.;
    double entries = 0.0;
    double temp_pos = 0.0;
    //double min = 1.e-3;       //bin_HR density range
    //double max = 1.e3;
    double min = 1.e6;
    double max = 1.e8;
    double abs_min = min;
    double abs_max = max;
    double time_in_Myr = 0.0;
    double grid1[N_grid][N_grid];
    double grid2[N_grid][N_grid];
    double grid3[1][N_grid];
    char dir[500];
    char buf[500];
    FILE *infile;
    FILE *outfile;
    int files;
    int N_gas;
    double hsml_factor=1.7;
    double center_i, center_j, W_x, weight, h, x2, x;
    double nh, nhmax, xmax, ymax, zmax, vx, vy, vz, xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;
	  
    files=1;

    double IDmass=0.0;
    int IDsink;

    if(halo == 1)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR1");
    if(halo == 2)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR2");
    if(halo == 3)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR3");
    if(halo == 4)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR4");
    if(halo == 5)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR5");
    if(halo == 6)
      sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR6");
    if(halo == 7)
	  sprintf(dir, "/work/00863/minerva/bin_zoom7");
	if(halo == 8)
	  sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR8");
	if(halo == 9)
	  sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHR/bin_HR9");
	if(halo == 10)
	  sprintf(dir, "/work/00863/minerva/bin_zoom10");
        if(halo == 11)
          sprintf(dir, "/work/00863/minerva/bin_zoom_cut/bin_zoom10_new_cut");
        if(halo == 12)
          sprintf(dir, "/scratch/00863/minerva/bin_zoom10_new_cut_ref2");	 
        if(halo == 13)
          sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHRL/bin_HRL7");
        if(halo == 14)
          sprintf(dir, "/Volumes/Apps_and_Docs/astacy/SimData/binHRL/bin_HRL10_alt");	 
 
 

    for(m = N_begin; m < N_end+1; m=m+1)
      {
        if(m < 10)
          {
            sprintf(buf, "%s_00%d", dir, m);

            printf("reading 00%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_00%d.nh", dir, m);

            printf("processing 00%d...\n", m);
          }
        else if(m > 9 && m < 100)
          {
            sprintf(buf, "%s_0%d", dir, m);

            printf("reading 0%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_0%d.nh", dir, m);

            printf("processing 0%d...\n", m);
          }
        else
          {
            sprintf(buf, "%s_%d", dir, m);

            printf("reading %d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_%d.nh", dir, m);

            printf("processing %d...\n", m);
          }

        //width= 50000.0*header.HubbleParam/(1.e3*header.time);  //50,000 pc
        //slice= 50000.0*header.HubbleParam/(1.e3*header.time);

        width= 0.2*header.HubbleParam/(1.e3*header.time);  //0.2 pc
        slice= 0.2*header.HubbleParam/(1.e3*header.time);       

        //width= 8.*25.*4.85e-4*header.HubbleParam/(1.e3*header.time);  //20,000 AU
        //slice= 8.*25.*4.85e-4*header.HubbleParam/(1.e3*header.time);

        //width= 100.0*4.85e-4*header.HubbleParam/(1.e3*header.time);  //10,000 AU
        //slice= 100.0*4.85e-4*header.HubbleParam/(1.e3*header.time);


        IDmass = 0.0;
        IDsink = 0;


        for(n = 0; n < N_gas; n++)
          {
            P[n].nh = P[n].Density*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0)*X/m_H;
            if (n%100000 ==0) printf("%lg\n", P[n].nh);
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
			  if(nh > nhmax/10.0)
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

        for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > -1 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < slice/2.0)
              {

h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0);
weight = P[n].nh * P[n].nh;

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

    grid1[i][j] += weight * P[n].nh * W_x;

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
weight = P[n].nh * P[n].nh;

                i = 0;
                j=0;

                i_min = int((P[n].Pos[2]-h-center_y+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[2]+h-center_y+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[1]-h-center_z+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[1]+h-center_z+width/2.0)/width*N_grid);

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

    grid2[i][j] += weight * P[n].nh * W_x;

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

//      for(n = 0; n < N_gas; n++)
//          {
 //           if(/*P[n].flag_split > 0.5 &&*/  /*Id[n]*/ P[n].Id == 4333881  && P[n].sink > 0.5 &&  double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
//              {
//                temp_pos = double(int((P[n].Pos[2]-center_z+width/2.0)/width*double(N_grid)))/double(N_grid);

//                fwrite(&temp_pos,sizeof(double),1,outfile);
//              }
//          }


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
      //printf("pos= %lg \n", P[pc_new].Pos[0]);
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
	      fread(&P[pc_sph].U, sizeof(double), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n=0, pc_sph=pc; n<header.npart[0];n++)
	    {
	      fread(&P[pc_sph].Density, sizeof(double), 1, fd);
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

/*
          SKIP;
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
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
          for(n=0, pc_sph=pc; n<header.npart[0];n++)
            {
              fread(&P[pc_sph].gam, sizeof(double), 1, fd);
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

