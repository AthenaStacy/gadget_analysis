#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
int allocate_memory(int NumPart);

struct particle_data
  {
    double Pos[3];
    double Vel[3];
    double Mass;
    int Type;
    double EgySpec;
    double Density;
    int Id;
    double pres_x;
    double pres_y;
    double pres_z;
    double visc_x;
    double visc_y;
    double visc_z;
    double Hsml;
    double sink;
    double nh;
    double vischeat;
    double time;
    double dt;
    //double Pressure;
    double mu;
    double gam;
    double DtDens;
    double adiab;
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
    int N_begin = 160;
    int N_end = 160;
    int n_grid1[N_grid][N_grid];
    int n_grid2[N_grid][N_grid];
    //double center_x = 50.339110552;
    //double center_y = 50.122129498;
    //double center_z = 49.486652504;
    double center_x = 50.339102858;
    double center_y = 50.122070744;
    double center_z = 49.486640197;
    //double center_x = 50.0;
    //double center_y = 50.0;
    //double center_z = 50.0;
    double width = 100.;
    double slice = 100.;
    double entries = 0.0;
    double temp_pos = 0.0;
    double min = 1.e-16;
    double max = 1.e-10;
    //double min = 1.e0;
    //double max = 1.e3;
    double abs_min = min;
    double abs_max = max;
    double time_in_Myr = 0.0;
    double grid1[N_grid][N_grid];
    double grid2[N_grid][N_grid];
    double grid3[1][N_grid];
    char dir[500];
    char dir2[500];
    char buf[500];
    FILE *infile;
    FILE *infile2;
    FILE *outfile;
    FILE *outfile2;
    int files;
    int N_gas;

    files=1;

    double IDmass=0.0;
    int IDsink;
    double time, HubbleParam=0.7, redshift;
    char input_fname[200], input_fname2[200];
    double Omega0=0.3;

    double parttime, dt, Mass, Pos0, Pos1, Pos2, nh, pres_x, pres_y, pres_z, visc_x, visc_y, visc_z, dummy, vischeat;
    int Id;
    double Hsml, sink, mu; 
    double hubble_a, AllHub, OmegaLambda, dtconv;
    double UnitEnergy_in_cgs, UnitTime_in_s, UnitMass_in_g, UnitVelocity_in_cm_per_s, UnitLength_in_cm, UnitDensity_in_cgs;
    double Dens, a3inv;

    UnitLength_in_cm= 3.085678e21;
    UnitVelocity_in_cm_per_s= 1.0e5;
    UnitMass_in_g= 1.989e43;
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
    OmegaLambda = 0.7;
    UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);

    //N_gas=29431;
    //N_gas=26518;
    //N_gas=26897; 
    //N_gas=25579;
    //N_gas=27426;
    //N_gas=25895;
    N_gas=25359;

    sprintf(dir, "/work/utexas/ao/minerva/torque");
    sprintf(dir2,"/work/utexas/ao/minerva/dens");

    for(m = N_begin; m < N_end+1; m++)
      {
   
            sprintf(buf, "%s_00%d", dir, m);

            sprintf(input_fname, "%s_%03d", dir, m);
            sprintf(input_fname2, "%s_%03d", dir2, m);
            
            printf("reading `%s' ...\n",input_fname);

            infile = fopen(input_fname, "r");
            for(i = 0; i < N_gas; i++)
	      {

               if(i==0)
                allocate_memory(N_gas);

             fscanf(infile, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d %lg %lg %lg %lg\n",   
                        &P[i].time, &P[i].dt, &P[i].Mass,
                        &P[i].Pos[0], &P[i].Pos[1], &P[i].Pos[2],
                        &P[i].nh,
                        &P[i].pres_x,
                        &P[i].pres_y,
                        &P[i].pres_z,
                        &P[i].visc_x,
                        &P[i].visc_y,
                        &P[i].visc_z,
                        &P[i].dummy,
                        &P[i].dummy,
                        &P[i].dummy,
                        &P[i].vischeat,
                        &P[i].Id,
                        &P[i].Hsml,
                        &P[i].sink,
                        &P[i].gam,
                        &P[i].mu); 
	        }
            fclose(infile);


           infile2 = fopen(input_fname2, "r"); 
           for(i = 0; i < N_gas; i++)
	      {
              fscanf(infile2, "%lg %lg \n", &P[i].DtDens, &P[i].EgySpec); 
	      }
            fclose(infile2);

          time = P[0].time; 
          redshift = (1.0/time) - 1.0;
          AllHub = Hubble*UnitTime_in_s;

          hubble_a = AllHub * sqrt(Omega0 / (time * time * time)
                                   + (1 - Omega0 - OmegaLambda) / (time * time) +
                                   OmegaLambda);
          dtconv = (1.0/(hubble_a*time))*UnitTime_in_s/HubbleParam;
          a3inv = 1.0 / (time * time * time);

          if(m < 10)
            {
            sprintf(buf, "%s_00%d.adiab", dir, m);
            }
          else if(m > 9 && m < 100)
            {
            sprintf(buf, "%s_0%d.adiab", dir, m);
            }
          else
            {
            sprintf(buf, "%s_%d.adiab", dir, m); 
            }

            printf("processing %d...\n", m);
 

	width= 2.4e-2*HubbleParam/(1.e3*time);   
        slice= 2.4e-2*HubbleParam/(1.e3*time);  

        IDmass = 0.0;
        IDsink = 0;

        outfile2 = fopen("/work/utexas/ao/minerva/visc_160", "a");
        for(n = 0; n < N_gas; n++)
          {
            //P[n].nh = P[n].Density*unit_mass/pow(unit_length,3.0)*pow(HubbleParam,2.0)/pow(time,3.0)*X/m_H;
            P[n].vischeat = P[n].vischeat*1.e10;  //convert from km^2/s^2 per time  to cm^2/s^2 per time
            P[n].vischeat = P[n].vischeat*P[n].mu*m_H*(P[n].gam-1.0); //convert from cm^2/s^2 per time to erg per time
            P[n].vischeat = P[n].vischeat*P[n].nh/dtconv;      //convert from erg per time to erg/(cm^3 s)

            Dens = (P[n].mu * m_H * P[n].nh)/(a3inv * UnitDensity_in_cgs * HubbleParam * HubbleParam); 
            //P[n].vischeat = (P[n].vischeat)*(P[n].gam-1.0) / pow(P[n].mu*m_H*P[n].nh, (P[n].gam - 1.0));   //use dEntropy/dt instead

            P[n].adiab = P[n].EgySpec*1.e10;  //convert from km^2/s^2 to cm^2/s^2
            P[n].adiab = P[n].adiab*fabs(P[n].DtDens)*(P[n].gam-1.0)/P[n].nh;  //convert from cm^2/s^2 to cm^2/s^2 per time
            P[n].adiab = P[n].adiab*P[n].mu*m_H*(P[n].gam-1.0); //convert from cm^2/s^2 per time to erg per time
            P[n].adiab = P[n].adiab*P[n].nh/dtconv;      //convert from erg per time to erg/(cm^3 s)

            P[n].EgySpec = P[n].EgySpec* UnitEnergy_in_cgs/ UnitMass_in_g;
            P[n].EgySpec = P[n].mu * m_H /k_B * (P[n].gam - 1.0)*P[n].EgySpec; 
            fprintf(outfile2, "%lg %lg %lg \n", P[n].nh, P[n].EgySpec, P[n].vischeat);

           if (n%100 ==0) printf("%lg %lg %lg %lg\n", P[n].nh, P[n].gam, P[n].mu, P[n].adiab);
          }
          fclose(outfile2);

        time_in_Myr = 2.0/3.0/Hubble/HubbleParam/sqrt(Omega0)*pow(time,3.0/2.0)/1.0e6/yr;

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
            if(double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < slice/2.0)
              {
                i = 0;
                j=0;

                i_min = int((P[n].Pos[0]-P[n].Hsml-center_x+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[0]+P[n].Hsml-center_x+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[1]-P[n].Hsml-center_y+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[1]+P[n].Hsml-center_y+width/2.0)/width*N_grid);

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

                                n_grid1[i][j] += 1;

                                grid1[i][j] += P[n].adiab;
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
                i = 0;
                j=0;

                i_min = int((P[n].Pos[1]-P[n].Hsml-center_y+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[1]+P[n].Hsml-center_y+width/2.0)/width*N_grid);

                if(i_max > N_grid-1)
                  {
                    i_max = N_grid-1;
                  }

                j_min = int((P[n].Pos[2]-P[n].Hsml-center_z+width/2.0)/width*N_grid);

                if(j_min < 0)
                  {
                    j_min = 0;
                  }

                j_max = int((P[n].Pos[2]+P[n].Hsml-center_z+width/2.0)/width*N_grid);

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

                                grid2[i][j] += P[n].adiab;
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
        fwrite(&redshift,sizeof(double),1,outfile);
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

       printf("IDsink = %d\n", IDsink);

        for(n = 0; n < N_gas; n++)
          {
            if(P[n].sink > 0.5 /*Id[n]*/ /*P[n].Id == IDsink*/  && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
              {
                entries++;
              }
          }

        //entries = 5.0;
        fwrite(&entries,sizeof(double),1,outfile);

        printf("entries %lg\n", entries);

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

        fclose(outfile);

        free(P);
      }

    printf("done!\n");
  }

/* this routine allocates the memory for the
 * particle data.
 */
int allocate_memory(int NumPart)
{
  printf("allocating memory...\n");

  if(!(P=(struct particle_data *) malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
}

