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
    double t0, dt_grav;
    int Id;
    double pres_x;
    double pres_y;
    double pres_z;
    double visc_x;
    double visc_y;
    double visc_z;
    double grav_x, grav_y, grav_z;
    double Hsml;
    double sink;
    double nh;
    double vischeat;
    double Temp;
    //double Pressure;
    double gam;
    double dummy;
    double H2I, HII; 
    double tgrav, tpres, tvisc, tgravx, tgravy, tgravz, tpresx, tpresy, tpresz, tviscx, tviscy, tviscz;
    double tx, ty, tz;
  } *P;
//int *Id;

int N_part;

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
    int N_begin = 4338;
    int N_end = 4338;
    double n_grid1[N_grid][N_grid];
    double n_grid2[N_grid][N_grid];
    double center_x = 50.339102858;
    double center_y = 50.122070744;
    double center_z = 49.486640197;
    double width = 100.;
    double slice = 100.;
    double entries = 0.0;
    double temp_pos = 0.0;
    double min = 1.e10;
    double max = 1.e13;
    //double min = 1.e0;
    //double max = 1.e3;
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
    FILE *outfile2;
    int files;
    int N_gas;

    files=1;

    double IDmass=0.0;
    int IDsink;
    double time, HubbleParam=0.7, redshift;
    char input_fname[200];
    double Omega0=0.3;

    double parttime, dt, Mass, Pos0, Pos1, Pos2, nh, pres_x, pres_y, pres_z, visc_x, visc_y, visc_z, dummy, vischeat;
    int Id, ID_sink1;
    double Hsml, sink, mu; 
    double hubble_a, AllHub, OmegaLambda, dtconv, fac1, fac2;
    double UnitEnergy_in_cgs, UnitTime_in_s, UnitMass_in_g, UnitVelocity_in_cm_per_s, UnitLength_in_cm, UnitDensity_in_cgs;
    double UnitAcc_in_cgs, a3inv;
    double sinkposx1, sinkposy1, sinkposz1, x1, y1, z1, vx1, vy1, vz1;
    double xtimesm, ytimesm, ztimesm, vxtimesm, vytimesm, vztimesm;
    double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, axCOM, ayCOM, azCOM, jzcent, diskmass;    
    double ncount_doub, nhmax, xmax, ymax, zmax, vx, vy, vz;
    double center_i, center_j, W_x, weight, h, x2, x;
    double hsml_factor=1.7;

    UnitLength_in_cm= 3.085678e21;
    UnitVelocity_in_cm_per_s= 1.0e5;
    UnitMass_in_g= 1.989e43;
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
    OmegaLambda = 0.7;
    UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
    UnitAcc_in_cgs = UnitLength_in_cm/pow(UnitTime_in_s,2);

    N_gas=79982;

    sprintf(dir, "/nobackupp1/astacy/torque");

    for(m = N_begin; m < N_end+1; m++)
      {
   
            sprintf(buf, "%s_00%d", dir, m);

            sprintf(input_fname, "%s_%03d", dir, m);
           
            if(m >= 1000)
              sprintf(input_fname, "%s_%04d", dir, m);
 
            printf("reading `%s' ...\n",input_fname);

            infile = fopen(input_fname, "r");

            for(i = 0; i < N_gas; i++)
	      {

               if(i==0)
                allocate_memory(N_gas);
            fscanf(infile, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
                        &P[i].t0, &P[i].dt_grav, &P[i].Mass,
                        &P[i].Pos[0], &P[i].Pos[1], &P[i].Pos[2],
                        &P[i].nh,
                        &P[i].pres_x, &P[i].pres_y, &P[i].pres_z,
                        &P[i].visc_x, &P[i].visc_y, &P[i].visc_z,
                        &P[i].grav_x, &P[i].grav_y, &P[i].grav_z,
                        &P[i].dummy,
                        &P[i].Id,
                        &P[i].Hsml,
                        &P[i].sink,
                        &P[i].gam,
                        &P[i].Temp,
                        &P[i].Vel[0], &P[i].Vel[1], &P[i].Vel[2],
                        &P[i].H2I,
                        &P[i].HII);
	        }

            fclose(infile);

          time = P[0].t0; 
          redshift = (1.0/time) - 1.0;
          AllHub = Hubble*UnitTime_in_s;

          hubble_a = AllHub * sqrt(Omega0 / (time * time * time)
                                   + (1 - Omega0 - OmegaLambda) / (time * time) +
                                   OmegaLambda);
          dtconv = (1.0/(hubble_a*time))*UnitTime_in_s/HubbleParam;
          a3inv = 1.0 / (time * time * time);
          printf("time = %lg, redshift = %lg\n", time, redshift);

          fac1 = 1 / (time * time);

          if(m < 10)
            {
            sprintf(buf, "%s_00%d.torque_nega_y", dir, m);
            }
          else if(m > 9 && m < 100)
            {
            sprintf(buf, "%s_0%d.torque_nega_y", dir, m);
            }
          else
            {
            sprintf(buf, "%s_%d.torque_nega_y", dir, m); 
            }

            printf("processing %d...\n", m);
 

	//width= 2.4e-2*HubbleParam/(1.e3*time);   
        //slice= 2.4e-2*HubbleParam/(1.e3*time);  

        width= 5*4.85e-4*HubbleParam/(1.e3*time);  //500 AU   
        slice= 5*4.85e-4*HubbleParam/(1.e3*time);

        IDmass = 0.0;
        IDsink = 0;

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
             //if(P[n].sink > 0.5)
             //  printf("Found a sink!, mass = %lg\n", P[n].Mass);
            }

         printf("nhmax = %lg\n", nhmax);

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

         sinkposx1 = center_x;
         sinkposy1 = center_y;
         sinkposz1 = center_z;

        printf("ncount_doub = %lg, center_x = %lg, center_y = %lg, center_z = %lg\n", ncount_doub, center_x, center_y, center_z);

        for(i = 0; i < N_gas; i++)
          {

            fac2 = 1 / pow(time, 3 * P[i].gam - 2);             

            //P[i].nh = P[i].Density*unit_mass/pow(unit_length,3.0)*pow(HubbleParam,2.0)/pow(time,3.0)*X/m_H;
            P[i].vischeat = P[i].vischeat*1.e10;  //convert from km^2/s^2 per time  to cm^2/s^2 per time
            P[i].vischeat = P[i].vischeat*1.22*m_H; //convert from cm^2/s^2 per time to erg per time
            P[i].vischeat = P[i].vischeat*P[i].nh/dtconv;      //convert from erg per time to erg/(cm^3 s)

            x1 = P[i].Pos[0];
            y1 = P[i].Pos[1];
            z1 = P[i].Pos[2];

/*
            P[i].pres_x = (P[i].pres_x)*pow(time, 0.5)*1.e5/dtconv;  //convert from code units to cm/s per time
            P[i].pres_y = (P[i].pres_y)*pow(time, 0.5)*1.e5/dtconv;  //also convert from cm/s per time to cm/s^2
            P[i].pres_z = (P[i].pres_z)*pow(time, 0.5)*1.e5/dtconv;  /for GADGET1 "Accel[k]" units!
            P[i].visc_x = (P[i].visc_x)*pow(time, 0.5)*1.e5/dtconv;
            P[i].visc_y = (P[i].visc_y)*pow(time, 0.5)*1.e5/dtconv;
            P[i].visc_z = (P[i].visc_z)*pow(time, 0.5)*1.e5/dtconv;
            P[i].grav_x = (P[i].grav_x)*pow(time, 0.5)*1.e5/dtconv;
            P[i].grav_y = (P[i].grav_y)*pow(time, 0.5)*1.e5/dtconv;
            P[i].grav_z = (P[i].grav_z)*pow(time, 0.5)*1.e5/dtconv;
*/

//NOTE: In Gadget2, the "Accel[k]" values have a different conversion to get physical acceleration units than used in Gadget1 

            P[i].pres_x = (P[i].pres_x)*fac2*UnitAcc_in_cgs;  //convert from code units to cm/s per time
            P[i].pres_y = (P[i].pres_y)*fac2*UnitAcc_in_cgs;  //also convert from cm/s per time to cm/s^2
            P[i].pres_z = (P[i].pres_z)*fac2*UnitAcc_in_cgs;
            P[i].visc_x = (P[i].visc_x)*fac2*UnitAcc_in_cgs;
            P[i].visc_y = (P[i].visc_y)*fac2*UnitAcc_in_cgs;
            P[i].visc_z = (P[i].visc_z)*fac2*UnitAcc_in_cgs;
            P[i].grav_x = (P[i].grav_x)*fac1*UnitAcc_in_cgs;
            P[i].grav_y = (P[i].grav_y)*fac1*UnitAcc_in_cgs;
            P[i].grav_z = (P[i].grav_z)*fac1*UnitAcc_in_cgs;

            P[i].tpresx = (y1-sinkposy1)*P[i].pres_z  - (z1-sinkposz1)*P[i].pres_y;
            P[i].tpresx = P[i].tpresx*3.086e18*1.e3*time/(0.7);   //convert to cm*(cm/s^2) 
            P[i].tpresy = (z1-sinkposz1)*P[i].pres_x - (x1-sinkposx1)*P[i].pres_z;
            P[i].tpresy = P[i].tpresy*3.086e18*1.e3*time/(0.7);
            P[i].tpresz = (x1-sinkposx1)*P[i].pres_y  -  (y1-sinkposy1)*P[i].pres_x;
            P[i].tpresz = P[i].tpresz*3.086e18*1.e3*time/(0.7);
            P[i].tpres = pow((P[i].tpresx*P[i].tpresx + P[i].tpresy*P[i].tpresy + P[i].tpresz*P[i].tpresz),0.5);

            P[i].tviscx = (y1-sinkposy1)*P[i].visc_z  - (z1-sinkposz1)*P[i].visc_y;
            P[i].tviscx = P[i].tviscx*3.086e18*1.e3*time/(0.7);
            P[i].tviscy = (z1-sinkposz1)*P[i].visc_x - (x1-sinkposx1)*P[i].visc_z;
            P[i].tviscy = P[i].tviscy*3.086e18*1.e3*time/(0.7);
            P[i].tviscz = (x1-sinkposx1)*P[i].visc_y  -  (y1-sinkposy1)*P[i].visc_x;
            P[i].tviscz = P[i].tviscz*3.086e18*1.e3*time/(0.7);
            P[i].tvisc = pow((P[i].tviscx*P[i].tviscx + P[i].tviscy*P[i].tviscy + P[i].tviscz*P[i].tviscz),0.5);

            P[i].tgravx = (y1-sinkposy1)*P[i].grav_z  - (z1-sinkposz1)*P[i].grav_y;
            P[i].tgravx = P[i].tgravx*3.086e18*1.e3*time/(0.7);
            P[i].tgravy = (z1-sinkposz1)*P[i].grav_x - (x1-sinkposx1)*P[i].grav_z;
            P[i].tgravy = P[i].tgravy*3.086e18*1.e3*time/(0.7);
            P[i].tgravz = (x1-sinkposx1)*P[i].grav_y  -  (y1-sinkposy1)*P[i].grav_x;
            P[i].tgravz = P[i].tgravz*3.086e18*1.e3*time/(0.7);
            P[i].tgrav = pow((P[i].tgravx*P[i].tgravx + P[i].tgravy*P[i].tgravy + P[i].tgravz*P[i].tgravz),0.5);

            P[i].tx = P[i].tgravx + P[i].tpresx + P[i].tviscx;
            P[i].ty = P[i].tgravy + P[i].tpresy + P[i].tviscy;
            P[i].tz = P[i].tgravz + P[i].tpresz + P[i].tviscz;


           if(P[i].tx > 0)
                 P[i].tx = 0;

            if(P[i].ty > 0)
                 P[i].ty = 0;

            if(P[i].tz > 0)
                 P[i].tz = 0;


           if(P[i].tx < 0)
                 P[i].tx = -P[i].tx;

            if(P[i].ty < 0)
                 P[i].ty = -P[i].ty;

            if(P[i].tz < 0)
                 P[i].tz = -P[i].tz;


            if (i%1000 ==0) printf("%lg %lg %lg %lg\n", P[i].nh, P[i].gam, P[i].Temp, P[i].ty);
          }

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

h = fmax(hsml_factor * P[n].Hsml, width / N_grid / 2.0);
weight = P[n].nh * P[n].nh;

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

    grid1[i][j] += weight * P[n].ty * W_x;

    n_grid1[i][j] += weight * W_x;
  }


                                //n_grid1[i][j] += 1;

                                //grid1[i][j] += P[n].ty;
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

                i_min = int((P[n].Pos[2]-P[n].Hsml-center_z+width/2.0)/width*N_grid);

                if(i_min < 0)
                  {
                    i_min = 0;
                  }

                i_max = int((P[n].Pos[2]+P[n].Hsml-center_z+width/2.0)/width*N_grid);

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

center_i = center_z - width / 2.0 + (i + 0.5) * width / (double) N_grid;

                        do
                          {
                            if(j >= j_min && j <= j_max)
                              {
                                flag_j = 1;

if(x2 <= 1.0)
  {
    x = sqrt(x2);

    if(x <= 0.5)
      W_x = 1.0 - 6.0 * x * x + 6.0 * x * x * x;
    else
      W_x = 2.0 * (1.0 - x) * (1.0 - x) * (1.0 - x);

    grid2[i][j] += weight * P[n].ty * W_x;

    n_grid2[i][j] += weight * W_x;
   }

                                //n_grid2[i][j] += 1;

                                //grid2[i][j] += P[n].ty;
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
             if(P[n].sink > 0.5 && /*Id[n]*/  P[n].Id != IDsink  &&  P[n].Id != 4025630 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
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
            if(P[n].sink > 0.5 && /*Id[n]*/ P[n].Id != IDsink  && P[n].Id != 4025630 && double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
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
            if(P[n].sink > 0.5 &&  /*Id[n]*/ P[n].Id != IDsink  && P[n].Id != 4025630 &&  double(fabs(P[n].Pos[0]-center_x)) < width/2.0 && double(fabs(P[n].Pos[1]-center_y)) < width/2.0 && double(fabs(P[n].Pos[2]-center_z)) < width/2.0)
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

