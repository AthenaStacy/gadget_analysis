//Adds gamma values of the particles to a new snapshot file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXREF 20
#define width 1.0e0

int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int reordering(void);
int unit_conversion(void);
int do_what_you_want(void);
int allocate_memory(void);

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


struct io_header_old
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
} header_old;


int     NumPart, Ngas, Ngas2;

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  double disx, disy, disz, dis_cm;
  double  Density, Rho, U, Temp, nh,  hsm;
  //double H2I, HII, DII, HDI, HeII, HeIII, gam, 
  double sink;
  double dummy;
} *P;

int *Id;

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
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  j, n, type, snapshot_number, files, Ngas, random, ncount, ncount_doub, ncount2;
  double x,y,z,x1,y1, z1, delr, dis, disx, disy, disz, disAU, nh, nhmax, hubble_param;
  double delx, dely, delz, xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, delvx, delvy, delvz, boxsize;
  FILE *outfile;

  files=1;     
  boxsize = 0.2;
  delx = dely = delz = 0.;
  xCOM = yCOM = zCOM = delvx = delvy = delvz = 0.;

  sprintf(input_fname, "/work/00863/minerva/iso_map_000");  //"nfw" denotes an r^-1 profile!
  sprintf(output_fname, "/work/00863/minerva/iso_map_001");

  Ngas = read_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  //for NON-comoving sims
  Time = hubble_param = 1.0;

  unit_conversion();

  nhmax = 0;
  for(n=0;n<=Ngas;n++)
     {
     nh = P[n].nh;
     if(nh > nhmax)
        {
        nhmax = nh;
        delx = P[n].Pos[0];  //del's based upon location of maximum density
        dely = P[n].Pos[1];
        delz = P[n].Pos[2];
        }
     }

  for(n=0;n <= Ngas; n++)
     {
     dis = pow(((P[n].Pos[0]-delx)*(P[n].Pos[0]-delx) + (P[n].Pos[1]-dely)*(P[n].Pos[1]-dely) + (P[n].Pos[2]-delz)*(P[n].Pos[2]-delz)), 0.5);
     dis=dis*1.e3*Time/(hubble_param);
     disAU=dis*206264.806;

     disx = (P[n].Pos[0]-delx)*1.e3*Time/(hubble_param);
     disy = (P[n].Pos[1]-dely)*1.e3*Time/(hubble_param);
     disz = (P[n].Pos[2]-delz)*1.e3*Time/(hubble_param);
     if(disx < width && disy < width && disz < width)
       {
       vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
       vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
       vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
       ncount_doub = ncount_doub + P[n].Mass;
       }
     }

   delvx = vxCOM/ncount_doub;
   delvy = vyCOM/ncount_doub;
   delvz = vzCOM/ncount_doub;

  printf("hello line 154, nhmax = %lg, delx = %lg, dely = %lg, delz = %lg\n", nhmax, delx, dely, delz);
  printf("delvx = %lg, delvy = %lg, delvz = %lg\n", delvx, delvy, delvz);

  Ngas2 = write_snapshot(output_fname, files, output_fname, delx, dely, delz, boxsize);

  ncount = 0;
  ncount2 = 0;

  printf("ncount = %d.\n", ncount);
  printf("ncount2 = %d.\n", ncount2);

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
  HubbleParam= 0.7e0;

  //For NON-comoving sims
  HubbleParam = 1.0;

  for(i=1; i<=NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */
	  MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  gamma= 5.0/3.0;
          //gamma = P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
          P[i].Rho= P[i].Density * UnitDensity_in_cgs;
	  P[i].nh= P[i].Density * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0)/ MeanWeight;
	  /*  printf("zred = %g", zred);*/
	}
    }
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize)
{
  FILE *fd;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,l,dummy,ntot_withmasses; 
  int    t,n,off,pc=0,pc_new=0,pc_sph;
  int NumPart_new = 0, Ngas_new = 0;
  int Idnew, nnew=1;
  double massnew, hsmnew, x, y, z, dis, nnew_doub=1.0;
  double randomx, randomy, randomz;
  double dmax = 5.0;  //size of 'cut-out' square in pc
  double mass_max = 1.e-10, hubble_param;

  //For NON-comoving setup 
  Time = hubble_param = 1.0;
	
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define WSKIP fwrite(&dummy, sizeof(dummy),1, outfile);

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
            {
	    NumPart+= header1.npart[k];
            printf("NumPart[%d] = %d\n", k, header1.npart[k]);
            }
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

      printf("ntot_withmasses = %d\n", ntot_withmasses);

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(double), 3, fd);
              P[pc_new].disx =fabs(1.e3*Time/(hubble_param)*(P[pc_new].Pos[0] - delx));
              P[pc_new].disy =fabs(1.e3*Time/(hubble_param)*(P[pc_new].Pos[1] - dely));
              P[pc_new].disz =fabs(1.e3*Time/(hubble_param)*(P[pc_new].Pos[2] - delz));
	      pc_new++;
	    }
	}
      SKIP;

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
              fread(&Id[pc_new], sizeof(int), 1, fd);
              pc_new++;
            }
        }
      SKIP;

      if(ntot_withmasses>0)
      {
        SKIP;
     }

      for(k=0, pc_new=pc; k<6; k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              P[pc_new].Type=k;
              if(header1.mass[k]==0)
                {
                fread(&P[pc_new].Mass, sizeof(double), 1, fd);
                }
              else
                P[pc_new].Mass= header1.mass[k];
              pc_new++;
            }
        }
      if(ntot_withmasses>0)
      {
        SKIP;
      }


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
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].Density, sizeof(double), 1, fd);
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
		
            printf("gam = %lg\n",P[100].gam);
	
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].sink, sizeof(double), 1, fd); 
				pc_sph++;
            }
			SKIP;
*/
        }

      fclose(fd);
       
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  return(Ngas);
}

			

//////////write new file!!!!!!!!!!!
	
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize)
{
  FILE *infile;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,l,dummy,ntot_withmasses;
  int    t,n,off,pc=0,pc_new=0,pc_sph=0;
  int NumPart_new = 0, Ngas_new = 0;
  int Idnew, nnew=1;
  int i0, imax= 197090;
  double *pos, *vel, massnew, hsmnew, x, y, z, dis, nnew_doub=1.0;
  double randomx, randomy, randomz;
  double dmax = 1.e10;  //size of 'cut-out' square in pc
  double mass_max = 1.e10, hubble_param;
  double cs_shu, U0, t_shu, r0_shu, rho0, dum1, dum2, dum3, dum4; 
  double r_shu[imax], x_shu[imax], rho_shu[imax], v_shu[imax], u_shu[imax], alpha[imax];  
  double pi=3.14159, grav=6.67e-8, kB = 1.38e-16, mH = 1.67e-24;

  //For NON-cosmological runs
  Time = hubble_param = 1.0;

/* 
 *    interpolate the shu solution for rho(r,t) and infall speed u(r,t)
 *     */
  /* From Shu (1977)
 *      x=r/(at)
 *      v(x) = u(r,t)/a
 *      M(r,t) =m(x)  a**3*t/G
 *      m = x**2 alpha*(x-v)
 *      rho(r,t) = alpha(x)/(4 pi G t**2)
 *      rho(r,0) = (a**2 A)/(4 pi G) r**(-2)*/


  r0_shu = 1.496e18; //length in cm, = 10^5 AU

  for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              P[pc_new].dis_cm = 3.08567758e18*1.e3*Time/(hubble_param)*pow(((P[pc_new].Pos[0]-delx)*(P[pc_new].Pos[0]-delx) + (P[pc_new].Pos[1]-dely)*(P[pc_new].Pos[1]-dely) + (P[pc_new].Pos[2]-delz)*(P[pc_new].Pos[2]-delz)), 0.5);
              P[pc_new].dummy = 1;
              if(P[pc_new].dis_cm > .999*r0_shu && P[pc_new].dis_cm < 1.001*r0_shu)
                 rho0 = P[pc_new].Rho;
              if(P[pc_new].nh > 1.e15 || P[pc_new].dis_cm > width*3.08567758e18)
                 P[pc_new].dummy = -5;
              pc_new++;
            }
        }


  //cs_shu = 180000;
  cs_shu = rho0*4.*pi*grav*pow(r0_shu,2)/2.0;
  cs_shu = pow(cs_shu,0.5);
  U0 = pow(cs_shu,2)*1.e-10;// factor is  UnitMass_in_g / UnitEnergy_in_cgs;
  printf("rho0 = %lg, cs_shu = %lg\n", rho0, cs_shu);

  infile=fopen("/home1/00863/minerva/research_programs/shu.dat","r");
  for(i0=0; i0<imax; i0++)
    {
    fscanf(infile,"%lg %lg %lg %lg\n", &dum1, &dum2, &dum3, &dum4);
    x_shu[i0] = dum1;
    alpha[i0] = dum2;
    v_shu[i0] = dum3;
    }
  fclose(infile);

  //t_shu = r0_shu/x_shu[imax-1]/cs_shu;
  t_shu = 3.e10;
  r0_shu = cs_shu * t_shu;
  printf("shu.dat file has been read. t_shu = %lg, r0_shu = %lg\n", t_shu, r0_shu);

  for(i0=0; i0<imax; i0++)
    {
    r_shu[i0] = x_shu[i0]  * cs_shu * t_shu;
    rho_shu[i0] = alpha[i0] / (4*pi*grav*(pow(t_shu,2)));
    u_shu[i0] = cs_shu * v_shu[i0];
    
    if(i0%10000 == 0)
      printf("x_shu[i0] = %lg, r_shu[i0] = %lg, v_shu[i0] = %lg\n", x_shu[i0], r_shu[i0], v_shu[i0]);
    }

  pos= (double*)malloc(sizeof(double) * 3);
  vel= (double*)malloc(sizeof(double) * 3);

  #define SKIP fread(&dummy, sizeof(dummy), 1, fd);
  #define WSKIP fwrite(&dummy, sizeof(dummy),1, outfile);

  for(i=0, pc=0; i<files; i++, pc=pc_new)
    {

     for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              //if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax)
              if(P[pc_new].dummy > 0)
                NumPart_new++;
              //if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
              if(P[pc_new].dummy > 0)
                Ngas_new++;
              pc_new++;
            }
        }

      for(k=0;k<MAXREF;k++)
        header_old.npart[k] = header1.npart[k];

      header1.BoxSize = boxsize;

      Ngas_new = Ngas_new*nnew;
      header1.npartTotal[0] = Ngas_new;
      header1.npartTotal[1] = 0;
      header1.npart[0] = Ngas_new;
      header1.npart[1] = 0;

      printf("Ngas= %6d \n",Ngas);
      for(k=0;k<6;k++)
        printf("Ngas_new %6d\n",header1.npart[k]);
     for(k=0;k<6;k++)
       printf("Ngas_new_tot %6d\n",header1.npartTotal[k]);
     for(k=0;k<6;k++)
       printf("Ngas_old_tot %6d\n",header1.npartTotal[k]);
      printf("NumPart_new %6d\n",NumPart_new);

    for(k=0, ntot_withmasses=0; k<5; k++)
       header1.mass[k]==0;

      for(k=0, ntot_withmasses=0; k<5; k++)
        {
          if(header1.mass[k]==0)
            ntot_withmasses+= header1.npart[k];
        }

     printf("ntot_withmasses = %d\n", ntot_withmasses);

      outfile=fopen(outname,"w");
		
      fwrite(&dummy, sizeof(dummy), 1, outfile);
      fwrite(&header1, sizeof(header1), 1, outfile);
      fwrite(&dummy, sizeof(dummy), 1, outfile);

      WSKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header_old.npart[k];n++)
            {

            if(P[pc_new].dummy > 0)
              {
//             if(P[n].dis_cm > r0_shu)
//                {
                pos[0]=P[pc_new].Pos[0]-delx + header1.BoxSize/2.;
                pos[1]=P[pc_new].Pos[1]-dely + header1.BoxSize/2.;
                pos[2]=P[pc_new].Pos[2]-delz + header1.BoxSize/2.;

                fwrite(&pos[0], sizeof(double), 3, outfile);
  //              }
/*
              if(P[n].dis_cm <= r0_shu)
                {
		   do{
		      randomx =rand();
		      randomx = randomx/RAND_MAX;
		      randomy =rand();
		      randomy = randomy/RAND_MAX;
		      randomz =rand();
		      randomz = randomz/RAND_MAX;
						
		      pos[0]=(1.0/3.08567758e18/1.e3/Time*hubble_param)*r0_shu*randomx;
		      pos[1]=(1.0/3.08567758e18/1.e3/Time*hubble_param)*r0_shu*randomy;
		      pos[2]=(1.0/3.08567758e18/1.e3/Time*hubble_param)*r0_shu*randomz;

		      dis = 3.08567758e18*1.e3*Time/hubble_param*pow((pos[0])*(pos[0]) + (pos[1])*(pos[1]) + (pos[2])*(pos[2]), 0.5);
		      }while(dis > r0_shu);
						
		    fwrite(&pos[0], sizeof(double), 3, outfile);	
                }
*/
              }
              pc_new++;
            }
        }
      WSKIP;


      WSKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header_old.npart[k];n++)
	    {
            if(P[pc_new].dummy > 0)
              {
              dis = pow((P[pc_new].Pos[0])*(P[pc_new].Pos[0]) + (P[pc_new].Pos[1])*(P[pc_new].Pos[1]) + (P[pc_new].Pos[2])*(P[pc_new].Pos[2]), 0.5);
              vel[0] = vel[1] = vel[2] = 0.0;
              //vel[0] = P[pc_new].Vel[0];
              //vel[1] = P[pc_new].Vel[1];
              //vel[2] = P[pc_new].Vel[2];

/*
             if(P[pc_new].dis_cm < r0_shu)
             {
             printf("dis_cm = %lg, r_shu = %lg, vel[0] = %lg\n", P[pc_new].dis_cm, r_shu[i0], vel[0]);
             vel[0] = P[pc_new].Vel[0];
             vel[1] = P[pc_new].Vel[1];
             vel[2] = P[pc_new].Vel[2];
             }
*/
/*
             for(i0=imax-1; i0>0; i0--)
               {
               //printf("dis_cm = %lg, r_shu = %lg, vel[0] = %lg\n", P[pc_new].dis_cm, r_shu[i0], vel[0]);
               vel[0] = u_shu[i0]*P[pc_new].Pos[0]/dis/1.e5;
               vel[1] = u_shu[i0]*P[pc_new].Pos[1]/dis/1.e5;
               vel[2] = u_shu[i0]*P[pc_new].Pos[2]/dis/1.e5;
	       if(P[pc_new].dis_cm > r_shu[i0])
                  {
                  //printf("i0 = %d, dis_cm = %lg, u_shu = %lg, r_shu = %lg, vel[0] = %lg\n", i0, P[pc_new].dis_cm, u_shu[i0], r_shu[i0], vel[0]);
	          break;
                  }
               } 
*/
             

              fwrite(&vel[0], sizeof(double), 3, outfile);
              }
	      pc_new++;
	    }
	}
      WSKIP;

      WSKIP;
      Idnew = 0;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header_old.npart[k];n++)
	    {
              //if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
              if(P[pc_new].dummy > 0)
                {
		for(l=0;l<nnew;l++)
		   {
		   //Idnew = Id[pc_new]*100 + l;
                   Idnew++;
                   fwrite(&Idnew, sizeof(int), 1, outfile);
	       	   }
		}
	      pc_new++;
	    }
	}
      WSKIP;


      if(ntot_withmasses>0)
      {
        WSKIP;
     }

      for(k=0, pc_new=pc; k<6; k++)
	{
	  for(n=0;n<header_old.npart[k];n++)
	    {
	      P[pc_new].Type=k;
	      if(header1.mass[k]==0)
                {
                //if(P[pc_new].disx < dmax && P[pc_new].disy < dmax && P[pc_new].disz < dmax && k ==0 && P[pc_new].Mass < mass_max)
                if(P[pc_new].dummy > 0)
                  {  
		  for(l=0;l<nnew;l++)
		    {
		    massnew = P[pc_new].Mass/nnew_doub;
                    fwrite(&massnew, sizeof(double), 1, outfile);
		    }
                  }
                }
	      else
		P[pc_new].Mass= header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses>0)
      {
        WSKIP;
      }

      printf("Mass = %lg\n", P[1000].Mass);
      printf("Npart= %6d \n",NumPart);
      printf("M_B= %15.6e \n",header1.mass[0]);
      printf("M_DM= %15.6e \n",header1.mass[1]);
 

   if(header_old.npart[0]>0)
	{
          WSKIP;
	  for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
	    {
              //if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max) 
              if(P[pc_sph].dummy > 0)
                {
                P[pc_sph].U = U0;
		for(l=0;l<nnew;l++)
	         fwrite(&P[pc_sph].U, sizeof(double), 1, outfile);
                }
	      pc_sph++;
	    }
          WSKIP;

          WSKIP;
	  for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
	    {
              //if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
              if(P[pc_sph].dummy > 0)
                {
		for(l=0;l<nnew;l++)
                  fwrite(&P[pc_sph].Density, sizeof(double), 1, outfile);  
                }
	      pc_sph++;
	    }
          WSKIP;

         

//divergence between gadget1 and gadget2 begins here!
          WSKIP;
          for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
            {
              //if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
              if(P[pc_sph].dummy > 0)
                {
		for(l=0;l<nnew;l++)
		  {
		   hsmnew = P[pc_sph].hsm*pow(nnew_doub,-.3333); 
                   fwrite(&hsmnew, sizeof(double), 1, outfile);
	       	  }
                }
              pc_sph++;
            }
          WSKIP;
/*
          WSKIP;
          for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
		for(l=0;l<nnew;l++)
		   {
                    fwrite(&P[pc_sph].H2I, sizeof(double), 1, outfile);
                    fwrite(&P[pc_sph].HII, sizeof(double), 1, outfile);
                    fwrite(&P[pc_sph].DII, sizeof(double), 1, outfile);
                    fwrite(&P[pc_sph].HDI, sizeof(double), 1, outfile);
                    fwrite(&P[pc_sph].HeII, sizeof(double), 1, outfile);
		    fwrite(&P[pc_sph].HeIII, sizeof(double), 1, outfile);
         	   }
                }
              pc_sph++;
            }
          WSKIP;

          WSKIP;
          for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
            {
              if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max)
                {
                 if(n%10000 == 0)
                   printf("gam = %lg\n",P[pc_sph].gam);
                 for(l=0;l<nnew;l++)
                  fwrite(&P[pc_sph].gam, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;
*/
          WSKIP;
          for(n=0, pc_sph=pc; n<header_old.npart[0];n++)
            {
              if(P[pc_sph].dummy > 0)
              //if(P[pc_sph].disx < dmax && P[pc_sph].disy < dmax && P[pc_sph].disz < dmax && P[pc_sph].Mass < mass_max) 
                {
                P[pc_sph].sink = 0;
		for(l=0;l<nnew;l++)
                  fwrite(&P[pc_sph].sink, sizeof(double), 1, outfile);
                }
              pc_sph++;
            }
          WSKIP;

 
	}

      fclose(outfile);
       
    }


  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  return(Ngas2);
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

  
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
//  Id--;   /* start with offset 1 */

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






  











