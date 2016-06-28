
/* Maps Gadget output to grid, together with map.pro  */

/* @ Thomas H. Greif @ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

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

void load_snapshot(char *fname);
void allocate_memory(void);

struct io_header
  {
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    int npartTotal[6];
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
    int ID;
    float Pos[3];
    float Vel[3];
    float Mass;
    float EgySpec;
    float Density;
    float elec;
    float HI;
    float HII;
    float HeI;
    float HeII;
    float HeIII;
    float H2I;
    float H2II;
    float HM;
    float Hsml;
    float DI;
    float DII;
    float HDI;
    float DM;
    float HDII;
    float nh;
    float Temp;
  } *P;

int N_part;
int N_gas;

int main(int argc, char **argv)
  {
    int i, j, k, m, n, i_min, i_max, j_min, j_max, flag_i, flag_j, start, end, local_start, local_end, N_files, N_files_per_task, overhead, Task, N_tasks, N_begin, N_end, slice_flag, box_flag, center_ID, N_grid = 700;

    int tot[3][N_grid][N_grid];

    float grid1[3][N_grid][N_grid], grid2[3][N_grid][N_grid], sink_pos[3];

    double mu, x, x2, W_x, center_i, center_j, width, min1, min2, max1, max2, center[3], proj_width, dummy1, dummy2;

    char dir[500], buf[500];

    FILE *outfile;

/* Free Parameters */

    sprintf(dir, "/work/00025/tgreif/test");			/* file path plus base name */

    N_begin = 0;							/* snapshot start */
    N_end = 999;							/* snapshot end */
    center_ID = 0;							/* center on this particle, only valid for box_flag = 0 */

    slice_flag = 0;							/* [0|1] for projection/slice mode */
    box_flag = 1;							/* [1] for entire box */

    width = 0.0;								/* width of target region in code units, only valid for box_flag = 0, make sure entire region is within box! */

    proj_width = 0.0;							/* width of projected axis in code units, only valid for slice_flag = 0 and box_flag = 0*/

    min1 = 1.0e-4;							/* minimum value for array1 */
    max1 = 10.0;							/* maximum value for array1 */
    min2 = 10.0;							/* minimum value for array2 */
    max2 = 1.0e4;							/* maximum value for array2 */

/* Free Parameters */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &Task);
    MPI_Comm_size(MPI_COMM_WORLD, &N_tasks);

    N_files = N_end-N_begin+1;

    N_files_per_task = N_files/N_tasks;

    overhead = N_files-N_files_per_task*N_tasks;

    start = end = local_start = local_end = 0;

    for(m = 0; m < N_tasks; m++)
      {
        start = end;

        end = start+N_files_per_task;

        if(m < overhead)
          end += 1;

        if(Task == m)
          {
            local_start = N_begin+start;

            local_end = N_begin+end;
          }
      }

    for(m = local_start; m < local_end; m++)
      {
        sprintf(buf, "%s_%03d", dir, m);

        printf("reading %03d...\n", m);

        load_snapshot(buf);

        sprintf(buf, "%s_%03d.map", dir, m);

        printf("processing %03d...\n", m);

        for(n = 0; n < N_part; n++)
          if(P[n].ID == center_ID)
            for(k = 0; k < 3; k++)
              center[k] = P[n].Pos[k];

        if(box_flag == 1)
          {
            width = proj_width = header.BoxSize;

            for(k = 0; k < 3; k++)
              center[k] = header.BoxSize/2.0;
          }

        for(k = 0; k < 3; k++)
          {
            flag_i = flag_j = 0;

            for(i = 0; i < N_grid; i++)
              for(j = 0; j < N_grid; j++)
                grid1[k][i][j] = grid2[k][i][j] = tot[k][i][j] = 0;

            for(n = 0; n < N_gas; n++)
              {
                if(fabs(P[n].Pos[0]-center[0]) <= width/2.0 && fabs(P[n].Pos[1]-center[1]) <= width/2.0 && fabs(P[n].Pos[2]-center[2]) <= width/2.0 && P[n].ID > 0 && (slice_flag == 0 && fabs(P[n].Pos[(k+2)%3]-center[(k+2)%3]) <= proj_width/2.0 || slice_flag == 1 && fabs(P[n].Pos[(k+2)%3]-center[(k+2)%3]) <= 2.0*P[n].Hsml))
                  {
                    if(header.mass[0] != 0.0)
                      P[n].Mass = header.mass[0];

                    mu = 4.0/(3.0*X+1.0+4.0*X*P[n].HII);

                    P[n].nh = P[n].Density*unit_mass/pow(unit_length,3.0)*pow(header.HubbleParam,2.0)/pow(header.time,3.0)*X/m_H;

                    P[n].Temp = P[n].EgySpec*unit_energy/unit_mass*mu*(gamma-1.0)*m_H/k_B;

/* Set both arrays here */

                    dummy1 = P[n].nh;
                    dummy2 = P[n].Temp;

/* Set both arrays here */

                    i = j = 0;

                    i_min = (int)((P[n].Pos[k%3]-P[n].Hsml-center[k%3]+width/2.0)/width*N_grid);

                    if(i_min < 0)
                      i_min = 0;

                    i_max = (int)((P[n].Pos[k%3]+P[n].Hsml-center[k%3]+width/2.0)/width*N_grid)+1;

                    if(i_max > N_grid-1)
                      i_max = N_grid-1;

                    j_min = (int)((P[n].Pos[(k+1)%3]-P[n].Hsml-center[(k+1)%3]+width/2.0)/width*N_grid);

                    if(j_min < 0)
                      j_min = 0;

                    j_max = (int)((P[n].Pos[(k+1)%3]+P[n].Hsml-center[(k+1)%3]+width/2.0)/width*N_grid)+1;

                    if(j_max > N_grid-1)
                      j_max = N_grid-1;

                    do
                      {
                        if(i >= i_min && i <= i_max)
                          {
                            flag_i = 1;

                            center_i = center[k%3]-width/2.0+(i+1.0/2.0)*width/(double)N_grid;

                            do
                              {
                                if(j >= j_min && j <= j_max)
                                  {
                                    flag_j = 1;

                                    center_j = center[(k+1)%3]-width/2.0+(j+1.0/2.0)*width/(double)N_grid;

                                    if(slice_flag == 1)
                                      {
                                        x2 = ((P[n].Pos[k%3]-center_i)*(P[n].Pos[k%3]-center_i)+(P[n].Pos[(k+1)%3]-center_j)*(P[n].Pos[(k+1)%3]-center_j)+(P[n].Pos[(k+2)%3]-center[(k+2)%3])*(P[n].Pos[(k+2)%3]-center[(k+2)%3]))/P[n].Hsml/P[n].Hsml;

                                        if(x2 <= 1.0)
                                          {
                                            x = sqrt(x2);

                                            if(x <= 1.0/2.0)
                                              W_x = 8.0/pi/P[n].Hsml/P[n].Hsml/P[n].Hsml*(1.0-6.0*x2+6.0*x*x2);
                                            else
                                              W_x = 8.0/pi/P[n].Hsml/P[n].Hsml/P[n].Hsml*2.0*(1.0-x)*(1.0-x)*(1.0-x);

                                            grid1[k][i][j] += P[n].Mass*dummy1/P[n].Density*W_x;
                                            grid2[k][i][j] += P[n].Mass*dummy2/P[n].Density*W_x;
                                          }
                                      }
                                    else
                                      {
                                        grid1[k][i][j] += dummy1;
                                        grid2[k][i][j] += dummy2;

                                        tot[k][i][j]++;
                                      }
                                  }
                                else if(j > j_max)
                                  flag_j = 2;
                                else
                                  flag_j = 0;

                                 j++;
                               }
                             while((flag_j == 0 || flag_j == 1) && j < N_grid);

                             j = 0;
                           }
                         else if(i > i_max)
                           flag_i = 2;
                         else
                           flag_i = 0;

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
                    if(tot[k][i][j] > 0)
                      {
                        grid1[k][i][j] /= (double)tot[k][i][j];
                        grid2[k][i][j] /= (double)tot[k][i][j];
                      }

                    if(grid1[k][i][j] < min1)
                      grid1[k][i][j] = min1;

                    if(grid2[k][i][j] < min2)
                      grid2[k][i][j] = min2;

                    if(grid1[k][i][j] > max1)
                      grid1[k][i][j] = max1;

                    if(grid2[k][i][j] > max2)
                      grid2[k][i][j] = max2;

                    grid1[k][i][j] = log10(grid1[k][i][j]);
                    grid2[k][i][j] = log10(grid2[k][i][j]);
                  }
              }

            grid1[k][0][0] = log10(min1);
            grid2[k][0][0] = log10(min2);
            grid1[k][N_grid-1][N_grid-1] = log10(max1);
            grid2[k][N_grid-1][N_grid-1] = log10(max2);
          }

        outfile = fopen(buf, "w");

        fwrite(&N_grid,sizeof(int),1,outfile);
        fwrite(&header.npart[5],sizeof(int),1,outfile);
        fwrite(&header.time,sizeof(double),1,outfile);
        fwrite(&width,sizeof(double),1,outfile);
        fwrite(&min1,sizeof(double),1,outfile);
        fwrite(&max1,sizeof(double),1,outfile);
        fwrite(&min2,sizeof(double),1,outfile);
        fwrite(&max2,sizeof(double),1,outfile);
        fwrite(grid1,sizeof(grid1),1,outfile);
        fwrite(grid2,sizeof(grid2),1,outfile);

        for(n = 0; n < header.npart[5]; n++)
          for(j = 0; j < 3; j++)
            {
              sink_pos[j] = (P[N_part-header.npart[5]+n].Pos[j]-center[j]+width/2.0)/width;

              fwrite(&sink_pos[j],sizeof(float),1,outfile);
            }

        fclose(outfile);

        free(P);
      }

    if(Task == 0)
      printf("done!\n");

    MPI_Finalize();
  }

void allocate_memory(void)
  {
    if(!(P = (struct particle_data *)malloc(N_part*sizeof(struct particle_data))))
      {
        fprintf(stderr, "failed to allocate memory.\n");
        exit(0);
      }
  }

void load_snapshot(char *fname)
  {
    FILE *fd;
    char buf[200];
    int i = 0;
    int n = 0;
    int sizeofheader = 256;
    int pc = 0;
    int pc_new = 0;
    int mass_flag = 0;

    #define SKIP fread(&sizeofheader, sizeof(int), 1, fd);

    sprintf(buf, "%s", fname);

    if(!(fd = fopen(buf, "r")))
      {
        printf("can't open file `%s`\n", buf);
        exit(0);
      }

    SKIP;
    fread(&header, sizeof(header), 1, fd);
    SKIP;

    for(i = 0, N_part = 0; i < 6; i++)
      {
        N_part += header.npart[i];
      }

    N_gas = header.npart[0];

    allocate_memory();

    SKIP;
    for(i = 0, pc_new = pc; i < 6; i++)
      {
        for(n = 0; n < header.npart[i]; n++)
          {
            fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
            pc_new++;
          }
      }
    SKIP;

    SKIP;
    for(i = 0, pc_new = pc; i < 6; i++)
      {
        for(n = 0; n < header.npart[i]; n++)
          {
            fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
            pc_new++;
          }
      }
    SKIP;

    SKIP;
    for(i = 0, pc_new = pc; i < 6; i++)
      {
        for(n = 0; n < header.npart[i]; n++)
          {
            fread(&P[pc_new].ID, sizeof(int), 1, fd);
            pc_new++;
          }
      }
    SKIP;

    for(i = 0; i < 6; i++)
      if(header.mass[i] == 0 && header.npart[i] > 0)
        mass_flag = 1;

    if(mass_flag == 1)
      SKIP;

    for(i = 0, pc_new = pc; i < 6; i++)
      if(header.mass[i] != 0)
        pc_new += header.npart[i];
      else
        for(n = 0; n < header.npart[i]; n++)
          {
            fread(&P[pc_new].Mass, sizeof(float), 1, fd);
            pc_new++;
          }

    if(mass_flag == 1)
      SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].EgySpec, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].Density, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].elec, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HI, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HII, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HeI, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HeII, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HeIII, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].H2I, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].H2II, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HM, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].Hsml, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].DI, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].DII, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HDI, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].DM, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    SKIP;
    for(n = 0, pc_new = pc; n < header.npart[0]; n++)
      {
        fread(&P[pc_new].HDII, sizeof(float), 1, fd);
        pc_new++;
      }
    SKIP;

    fclose(fd);
  }
