#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXREF 20
#define pi 3.1415927
#define m_H 1.6726e-24
#define k_B 1.3806e-16
#define Hubble 3.2407789e-18
#define G 6.672e-8
#define X 0.76
#define yr 3.1536e7
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
    double Pos[3];
    double Vel[3];
    double Mass;
    int Type;
    int Id;

    double  EgySpec, Temp, nh, Density, Hsml;
    double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
    double dummy;

  } *P;  
//int *Id;

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
    int N_grid = 700;
    int N_begin = 152;
    int N_end = 152;
 //   int n_grid1[N_grid][N_grid];
    double n_grid1[N_grid][N_grid];
    int n_grid2[N_grid][N_grid];

    double center_x = 70.285976684;  //hires_nf center
    double center_y = 70.266798846;
    double center_z = 69.644070258;

    //double center_x = 50.339339267;
    //double center_y = 50.121980449;
    //double center_z = 49.486589836;

    //  double center_x = 50.339354975;
    //  double center_y = 50.121967594;
    //  double center_z = 49.48659644;

    //double center_x = 50.339110394;
    //double center_y = 50.122129369;
    //double center_z = 49.486652645;
  
    //double center_x = 70.0;
    //double center_y = 70.0;
    //double center_z = 70.0;

    double width = 140.0;
    double slice = 140.0;
    double entries = 0.0;
    double temp_pos = 0.0;
    double min = 1.0e2;
    double max = 1.0e4;
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
   
    double IDmass=0.0;
    int IDsink; 
    double mu, u, h2frac, muh2, muh2in, gamma;   
    int files;
 
    double hsml_factor=1.7;
    double center_i, center_j, W_x, weight, h, x2, x;

    files=1;

    //sprintf(dir, "/work/00863/minerva/ds_nr6e");
    //sprintf(dir, "/work/00863/minerva/tgridsink_ifront/tgridsink_ifront_nr6e");
    sprintf(dir, "/nobackupp1/astacy/hires_test6");
    //sprintf(dir, "/nobackupp1/astacy/hires_nf");
 
    for(m = N_begin; m < N_end+1; m=m+10)
      {
        if(m < 10)
          {
            sprintf(buf, "%s_00%d", dir, m);

            printf("reading 00%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_00%d.Tempxz", dir, m);

            printf("processing 00%d...\n", m);
          }
        else if(m > 9 && m < 100)
          {
            sprintf(buf, "%s_0%d", dir, m);

            printf("reading 0%d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_0%d.Tempxz", dir, m);

            printf("processing 0%d...\n", m);
          }
        else
          {
            sprintf(buf, "%s_%d", dir, m);

            printf("reading %d...\n", m);

            N_gas=load_snapshot(buf, files);

            sprintf(buf, "%s_%d.Tempxz", dir, m);

            printf("processing %d...\n", m);
          }
        
	printf("Point 1\n");

        width= 5*4.85e-4*header.HubbleParam/(1.e3*header.time);  //500 AU
        slice= 5*4.85e-4*header.HubbleParam/(1.e3*header.time);

        //width= 2.4e-2*header.HubbleParam/(1.e3*header.time);
        //slice= 2.4e-2*header.HubbleParam/(1.e3*header.time);

        //width= 2e0*header.HubbleParam/(1.e3*header.time);   
        //slice= 2e0*header.HubbleParam/(1.e3*header.time);  


         IDmass = 0.0;
         IDsink = 0;

        for(n = 0; n < N_gas; n++)
          {
              //P[n].nh = P[n].Density*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0)*X/m_H;
              //mu = 4.0/(3.0*X+1.0+4.0*X*P[n].elec);
              mu=1.2195;
              h2frac=2.0*P[n].H2I;

              muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
              muh2=pow(muh2in, -1.0);

              if(muh2 >= 1.22){
                 mu=muh2;
                }

              gamma = P[n].gam;
   
              P[n].Temp = P[n].EgySpec*unit_energy/unit_mass*mu*(gamma-1.0)*m_H/k_B;
              P[n].nh = P[n].Density*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0)*X/m_H;             
    
              if(n%100000 == 0)
                printf("Temp= %lg \n", P[n].Temp);


          }


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
	printf("Point 2\n");
        for(n = 0; n < N_gas; n++)
          {
            if(fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_y) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0 && fabs(P[n].Pos[2]-center_z) < slice/2.0)
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

center_i = center_x - width / 2.0 + (i + 0.5) * width / (double) N_grid;

                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

center_j = center_z - width / 2.0 + (j + 0.5) * width / (double) N_grid;

x2 = ((P[n].Pos[0] - center_i) * (P[n].Pos[0] - center_i) + (P[n].Pos[2] -
center_j) * (P[n].Pos[2] - center_j)) / h / h;

if(x2 <= 1.0)
  {
    x = sqrt(x2);

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

    grid1[i][j] += weight * P[n].Temp * W_x;

    n_grid1[i][j] += weight * W_x;
  }


                 //               n_grid1[i][j] += 1;

                   //             grid1[i][j] += P[n].Temp;
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
	printf("Point 3\n");
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
	printf("Point 4\n");
        for(i = 0; i < N_grid; i++)
          {
            for(j = 0; j < N_grid; j++)
              {
                n_grid2[i][j] = 0;

                grid2[i][j] = 0.0;
              }
          }
	printf("Point 5\n");
        for(n = 0; n < N_gas; n++)
          {
            if(fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_x) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0 && fabs(P[n].Pos[0]-center_x) < slice/2.0)
              {
                i = 0;
                j=0;

                i_min = int((P[n].Pos[2]-h-center_z+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[2]+h-center_z+width/2.0)/width*N_grid);

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

                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

                                n_grid2[i][j] += 1;

                                grid2[i][j] += P[n].Temp;
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
	printf("Point 6\n");
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
	printf("Point 7\n");

	/*
       for(n = 1; n <= N_gas; n++)
          {
            if(P[n].flag_split > 0.5 (Id[n] == 4025757 || Id[n] == 3410514 || Id[n] == 4014428) && fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_y) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0)
              {
                entries++;
              }
          }

        fwrite(&entries,sizeof(double),1,outfile);

        for(n = 1; n <= N_gas; n++)
          {
            if(P[n].flag_split > 0.5 (Id[n] == 4025757 || Id[n] == 3410514 || Id[n] == 4014428) && fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_y) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[0]-center_x+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

        for(n = 1; n <= N_gas; n++)
          {
            if(P[n].flag_split > 0.5 (Id[n] == 4025757 || Id[n] == 3410514 || Id[n] == 4014428) && fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_y) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[1]-center_y+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }

        for(n = 1; n <= N_gas; n++)
          {
            if(P[n].flag_split > 0.5 (Id[n] == 4025757 || Id[n] == 3410514 || Id[n] == 4014428)  && fabs(P[n].Pos[0]-center_x) < width/2.0 && fabs(P[n].Pos[1]-center_y) < width/2.0 && fabs(P[n].Pos[2]-center_z) < width/2.0)
              {
                temp_pos = double(int((P[n].Pos[2]-center_z+width/2.0)/width*double(N_grid)))/double(N_grid);

                fwrite(&temp_pos,sizeof(double),1,outfile);
              }
          }
   */

        for(n = 0; n < N_gas; n++)
          {
           if(P[n].sink > 0.5 && P[n].Mass > IDmass)
             {
             IDmass = P[n].Mass;
             IDsink = P[n].Id;
             }
          } 


     for(n = 0; n < N_gas; n++)
          {            if(P[n].sink > 0.5 /*Id[n]*/ /*P[n].Id == IDsink*/  && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
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
            if(/*P[n].flag_split > 0.5 &&*/  /*Id[n]*/ P[n].Id == 4333881  && P[n].sink > 0.5 &&  double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
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

	printf("Point 10\n");
        fclose(outfile);
	printf("Point 11\n");
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
  
  printf("allocating memory...done\n");
}

int load_snapshot(char *fname, int files)
  {
  FILE *fd;
  char   buf[200];
  int    l,j,k,dummy,ntot_withmasses;
  int    t,n,off,pc=0,pc_new=0,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

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
      printf("pos= %g \n", &P[pc_new].Pos[0]);
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
              fread(&P[pc_sph].EgySpec, sizeof(double), 1, fd);
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

