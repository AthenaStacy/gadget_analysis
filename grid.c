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
#define gamma 5.0/3.0
#define unit_length 3.085678e21
#define unit_mass 1.989e43
#define unit_energy 1.989e53

int main(int argc, char **argv)
  {
    int i = 0;
    int j = 0;
    int k = 0;
    int N_grid = 180;
    double boxsize = 140.0;
    //double center_x = 114.4868774;
    //double center_y = 128.1215363;
    //double center_z = 60.31565094;
    //double center_x = 0.9632437069;
    //double center_y = 67.25763737;
    //double center_z = 49.10663532;
    double center_x = 132.972017; 
    double center_y = 97.45263718;
    double center_z = 70.35441062;
    double width_1 = 50.0;
    double width_2 = 45.0;
    double width_3 = 40.0;
    int zoom_1 = int(width_1/boxsize/2.0*N_grid);
    int zoom_2 = int(width_2/boxsize/2.0*N_grid);
    int zoom_3 = int(width_3/boxsize/2.0*N_grid);
    int flag_i_right_1 = 0;
    int flag_i_left_1 = 0;
    int flag_j_right_1 = 0;
    int flag_j_left_1 = 0;
    int flag_k_right_1 = 0;
    int flag_k_left_1 = 0;
    int flag_i_right_2 = 0;
    int flag_i_left_2 = 0;
    int flag_j_right_2 = 0;
    int flag_j_left_2 = 0;
    int flag_k_right_2 = 0;
    int flag_k_left_2 = 0;
    int flag_i_right_3 = 0;
    int flag_i_left_3 = 0;
    int flag_j_right_3 = 0;
    int flag_j_left_3 = 0;
    int flag_k_right_3 = 0;
    int flag_k_left_3 = 0;
    int i_right_1 = int(center_x/boxsize*N_grid)+zoom_1;
    int i_left_1 = int(center_x/boxsize*N_grid)-zoom_1;
    int j_right_1 = int(center_y/boxsize*N_grid)+zoom_1;
    int j_left_1 = int(center_y/boxsize*N_grid)-zoom_1;
    int k_right_1 = int(center_z/boxsize*N_grid)+zoom_1;
    int k_left_1 = int(center_z/boxsize*N_grid)-zoom_1;
    int i_right_2 = int(center_x/boxsize*N_grid)+zoom_2;
    int i_left_2 = int(center_x/boxsize*N_grid)-zoom_2;
    int j_right_2 = int(center_y/boxsize*N_grid)+zoom_2;
    int j_left_2 = int(center_y/boxsize*N_grid)-zoom_2;
    int k_right_2 = int(center_z/boxsize*N_grid)+zoom_2;
    int k_left_2 = int(center_z/boxsize*N_grid)-zoom_2;
    int i_right_3 = int(center_x/boxsize*N_grid)+zoom_3;
    int i_left_3 = int(center_x/boxsize*N_grid)-zoom_3;
    int j_right_3 = int(center_y/boxsize*N_grid)+zoom_3;
    int j_left_3 = int(center_y/boxsize*N_grid)-zoom_3;
    int k_right_3 = int(center_z/boxsize*N_grid)+zoom_3;
    int k_left_3 = int(center_z/boxsize*N_grid)-zoom_3;
    char mesh[N_grid*N_grid*N_grid];

    FILE *outfile;

    printf("writing...\n");

    outfile = fopen("/work/utexas/ao/minerva/testc_180.dat", "w");

    fwrite(&N_grid, sizeof(int), 1, outfile);

    if(i_right_1 >= N_grid)
      {
        i_right_1 -= N_grid;
        flag_i_right_1 = 1;
      }

    if(i_left_1 < 0)
      {
        i_left_1 += N_grid;
        flag_i_left_1 = 1;
      }

    if(j_right_1 >= N_grid)
      {
        j_right_1 -= N_grid;
        flag_j_right_1 = 1;
      }

    if(j_left_1 < 0)
      {
        j_left_1 += N_grid;
        flag_j_left_1 = 1;
      }

    if(k_right_1 >= N_grid)
      {
        k_right_1 -= N_grid;
        flag_k_right_1 = 1;
      }

    if(k_left_1 < 0)
      {
        k_left_1 += N_grid;
        flag_k_left_1 = 1;
      }

    if(i_right_2 >= N_grid)
      {
        i_right_2 -= N_grid;
        flag_i_right_2 = 1;
      }

    if(i_left_2 < 0)
      {
        i_left_2 += N_grid;
        flag_i_left_2 = 1;
      }

    if(j_right_2 >= N_grid)
      {
        j_right_2 -= N_grid;
        flag_j_right_2 = 1;
      }

    if(j_left_2 < 0)
      {
        j_left_2 += N_grid;
        flag_j_left_2 = 1;
      }

    if(k_right_2 >= N_grid)
      {
        k_right_2 -= N_grid;
        flag_k_right_2 = 1;
      }

    if(k_left_2 < 0)
      {
        k_left_2 += N_grid;
        flag_k_left_2 = 1;
      }

    if(i_right_3 >= N_grid)
      {
        i_right_3 -= N_grid;
        flag_i_right_3 = 1;
      }

    if(i_left_3 < 0)
      {
        i_left_3 += N_grid;
        flag_i_left_3 = 1;
      }

    if(j_right_3 >= N_grid)
      {
        j_right_3 -= N_grid;
        flag_j_right_3 = 1;
      }

    if(j_left_3 < 0)
      {
        j_left_3 += N_grid;
        flag_j_left_3 = 1;
      }

    if(k_right_3 >= N_grid)
      {
        k_right_3 -= N_grid;
        flag_k_right_3 = 1;
      }

    if(k_left_3 < 0)
      {
        k_left_3 += N_grid;
        flag_k_left_3 = 1;
      }

    for(i = 0; i < N_grid; i++)
      {
        for(j = 0; j < N_grid; j++)
          {
            for(k = 0; k < N_grid; k++)
              {
                mesh[i*N_grid*N_grid+j*N_grid+k] = 0;
              }
          }
      }

    for(i = 0; i < N_grid; i++)
      {
        for(j = 0; j < N_grid; j++)
          {
            for(k = 0; k < N_grid; k++)
              {
                if((flag_i_right_1 == 0 && flag_i_left_1 == 0 && i <= i_right_1 && i >= i_left_1 || flag_i_right_1 == 1 && flag_i_left_1 == 0 && (i <= i_right_1 || i >= i_left_1) || flag_i_right_1 == 0 && flag_i_left_1 == 1 && (i >= i_left_1 || i <= i_right_1)) && (flag_j_right_1 == 0 && flag_j_left_1 == 0 && j <= j_right_1 && j >= j_left_1 || flag_j_right_1 == 1 && flag_j_left_1 == 0 && (j <= j_right_1 || j >= j_left_1) || flag_j_right_1 == 0 && flag_j_left_1 == 1 && (j >= j_left_1 || j <= j_right_1)) && (flag_k_right_1 == 0 && flag_k_left_1 == 0 && k <= k_right_1 && k >= k_left_1 || flag_k_right_1 == 1 && flag_k_left_1 == 0 && (k <= k_right_1 || k >= k_left_1) || flag_k_right_1 == 0 && flag_k_left_1 == 1 && (k >= k_left_1 || k <= k_right_1)))
                  {
                    mesh[i*N_grid*N_grid+j*N_grid+k] = 1;
                  }
              }
          }
      }

    for(i = 0; i < N_grid; i++)
      {
        for(j = 0; j < N_grid; j++)
          {
            for(k = 0; k < N_grid; k++)
              {
                if((flag_i_right_2 == 0 && flag_i_left_2 == 0 && i <= i_right_2 && i >= i_left_2 || flag_i_right_2 == 1 && flag_i_left_2 == 0 && (i <= i_right_2 || i >= i_left_2) || flag_i_right_2 == 0 && flag_i_left_2 == 1 && (i >= i_left_2 || i <= i_right_2)) && (flag_j_right_2 == 0 && flag_j_left_2 == 0 && j <= j_right_2 && j >= j_left_2 || flag_j_right_2 == 1 && flag_j_left_2 == 0 && (j <= j_right_2 || j >= j_left_2) || flag_j_right_2 == 0 && flag_j_left_2 == 1 && (j >= j_left_2 || j <= j_right_2)) && (flag_k_right_2 == 0 && flag_k_left_2 == 0 && k <= k_right_2 && k >= k_left_2 || flag_k_right_2 == 1 && flag_k_left_2 == 0 && (k <= k_right_2 || k >= k_left_2) || flag_k_right_2 == 0 && flag_k_left_2 == 1 && (k >= k_left_2 || k <= k_right_2)))
                  {
                    mesh[i*N_grid*N_grid+j*N_grid+k] = 2;
                  }
              }
          }
      }

    for(i = 0; i < N_grid; i++)
      {
        for(j = 0; j < N_grid; j++)
          {
            for(k = 0; k < N_grid; k++)
              {
                if((flag_i_right_3 == 0 && flag_i_left_3 == 0 && i <= i_right_3 && i >= i_left_3 || flag_i_right_3 == 1 && flag_i_left_3 == 0 && (i <= i_right_3 || i >= i_left_3) || flag_i_right_3 == 0 && flag_i_left_3 == 1 && (i >= i_left_3 || i <= i_right_3)) && (flag_j_right_3 == 0 && flag_j_left_3 == 0 && j <= j_right_3 && j >= j_left_3 || flag_j_right_3 == 1 && flag_j_left_3 == 0 && (j <= j_right_3 || j >= j_left_3) || flag_j_right_3 == 0 && flag_j_left_3 == 1 && (j >= j_left_3 || j <= j_right_3)) && (flag_k_right_3 == 0 && flag_k_left_3 == 0 && k <= k_right_3 && k >= k_left_3 || flag_k_right_3 == 1 && flag_k_left_3 == 0 && (k <= k_right_3 || k >= k_left_3) || flag_k_right_3 == 0 && flag_k_left_3 == 1 && (k >= k_left_3 || k <= k_right_3)))
                  {
                    mesh[i*N_grid*N_grid+j*N_grid+k] = 3;
                  }
              }
          }
      }

    fwrite(&mesh, sizeof(char), N_grid*N_grid*N_grid, outfile);

    fclose(outfile);

    printf("done\n");
  }
