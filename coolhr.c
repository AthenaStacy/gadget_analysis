//CR cooling is now added
//high density cooling and heating (H2 formation heating and enhanced H2 cooling) are also added

void  hmcool(double nh2, double nh,double nhd, double ncr, double temp,double Tcmb,double *pH2,double *pHD)
/*=======================================================================
** Evaluates H2 and HD cooling function!
=======================================================================*/
{
  /* I've added hdcool declaration below, as compiler said I needed it.*/
  void hdcool(double dens,double temp,double *pw_hd,double dw_hddT);     
  double LH2,LLTE,XNCRIT,LLOW,T3, lamH2, lamHD;
  double LOGLOW,LRLTE,LVLTE;
  double T3C,LOGLOWC,templogc,templog;
  double LLOWC,LRLTEC,LVLTEC,LLTEC,XNCRITC;
  double w_hd,dw_hddT,XHDm,XHDg,w_hdg,xdummy;
  double lamCR, w_CR;
  FILE *outfile;

     
/*=======================================================================
**    include H2 cooling
**    (Galli & Palla 1998)
=======================================================================*/
  T3=temp/1000.e0;
  T3C=Tcmb/1000.e0;
  templog=log10(temp);
  templogc=log10(Tcmb);
  LOGLOW=-103.0 + 97.59*templog - 48.05*pow(templog,2) +
        10.80*pow(templog,3) - 0.9032*pow(templog,4); 
  LOGLOWC=-103.0 + 97.59*templogc - 48.05*pow(templogc,2) +
        10.80*pow(templogc,3) - 0.9032*pow(templogc,4);
  LLOW=pow(10.e0,LOGLOW);
  LLOWC=pow(10.e0,LOGLOWC);
  LRLTE=(9.5e-22*pow(T3,3.76)*exp(-pow(0.13/T3,3))/
        (1.e0+0.12*pow(T3,2.1))+3.e-24*exp(-0.51/T3))/nh;
  LRLTEC=(9.e-22*pow(T3C,3.76)*exp(-pow(0.13/T3C,3))/
        (1.e0+0.12*pow(T3C,2.1))+3.e-24*exp(-0.51/T3C))/nh;
  LVLTE=(6.7e-19*exp(-5.86/T3) + 1.6e-18*exp(-11.7/T3))/nh;
  LVLTEC=(6.7e-19*exp(-5.86e0/T3C) + 
         1.6e-18*exp(-11.7e0/T3C))/nh;
  LLTE=LRLTE + LVLTE;
  LLTEC=LRLTEC + LVLTEC;
  XNCRIT=nh*(LLTE/LLOW);
  XNCRITC=nh*(LLTEC/LLOWC);
  if(nh>0.e0){
  LH2=LLTE/(1.e0+XNCRIT/nh) - LLTEC/(1.e0+XNCRITC/nh);
  lamH2=nh2*nh*LH2;
  }
  else{
    lamH2=1.e-52;
  }
  //lamH2 = lamH2/(nh*nh);
 //if(dens >= 1.e8){
    //if(i%1000==0){
     outfile=fopen("lamH2dat", "a");
     fprintf(outfile, "%g %g\n", nh + nh2, lamH2);
     fclose(outfile);
    //}
     //}

 if((nh+nh2)>1.e10)
 {
  lamH2=lamH2*pow((nh+nh2)/1.e10,-0.2614);
 }

  *pH2=lamH2;

/*======================================================
**    include HD cooling
======================================================*/
  hdcool(nh,temp,&w_hd,dw_hddT);
  hdcool(nh,Tcmb,&w_hdg,xdummy);
  XHDm=nhd*w_hd;
  XHDg=nhd*w_hdg;
  lamHD=XHDm-XHDg; 
  /*lamHD=0.e0;*/
  /* JJ added the below line to have lamHD units match the lamH2 and lamHl units.*/
  //lamHD=lamHD/(nh*nh); 
  //printf("%g %g %g\n", lamHD,XHDm, XHDg);
  *pHD=lamHD;
  return;

}



/* JJ took off the static designation at the beginning of the below line.*/
void hdcool(double densin,double temp,double *pw_hd,double dw_hddT)
/*======================================================
***  Flower et al. (2000) HD cooling function
+======================================================*/
{

  double aa,bb,omeg,phi,c1,c2,d1,d2;
  double dens,x,y,w,dwdx;
  aa=  -26.2982e0;
  bb= -0.215807e0;
  omeg= 2.03657e0;
  phi=  4.63258e0;
  c1= 0.283978e0;
  c2= -1.27333e0;
  d1= -2.08189e0;
  d2= 4.66288e0;

  dens=densin;
  if (densin < 1.e0) dens=1.e0;
  if (densin > 1.e8) dens=1.e8;
  y = log10(dens);
  x = log10(temp);

  w = 0.5e0 * y + aa * pow(x,bb)
       - sqrt( 0.25e0 * y*y
       + (c1*y+c2) * sin(omeg*x+phi) + (d1*y+d2) );

  dwdx =  aa*bb*pow(x,(bb-1.e0))
       -0.5e0/sqrt( 0.25e0 * y*y
       + (c1*y+c2) * sin(omeg*x+phi) + (d1*y+d2) )*
       (c1*y+c2) * cos(omeg*x+phi)*omeg;

  *pw_hd = pow(10.0e0,w);
  dw_hddT=(*pw_hd/temp)*dwdx;
  if (temp < 30.e0 || temp > 3000.e0) {
    *pw_hd = 0.e0;
    dw_hddT=1.e0;
  }
  return;
}

void crcool(double dens, double temp, double *plamCR, double CRheat)
  {
  double w_CR,lamCR;
  //w_CR=-2.e-29/(pow(.323,0.5));
  w_CR=-CRheat;
  //w_CR=0;
  //lamCR=w_CR*pow(dens,0.5);
  lamCR=w_CR*dens;
  *plamCR=lamCR;
  }

void protoheat(double dens, double temp, double *pgampr, double lum)
 {
   double rad, flux, gampr;
   gampr = 0.0;
   rad = 0.09*pow(dens/1.e6,-0.5);   //convert density to radius using r^-2 isothermal profile
                                 //radius is in pc
   rad = rad*3.0857e18;         //convert radius to cm
   lum = 1.0*3.839e33;           //luminosity in erg/s
   flux = lum/(4.0*3.14159*pow(rad, 2.0));
   gampr = lum/pow(rad, 3.0);
   //*pgampr = -gampr;
   *pgampr = 0.0;
 }

//AS has added term to account for H2-formation heating
/*======================================================================
**evaluates H2-formation heating
========================================================================*/
void hrheat(double dens, double temp, double time, double h2frac, double nHI, double nh2p, double nhm,  double *pgamHR)
 {
 double en, rate, gamHR;
 int m=0;
 double k8, k10, ncr, gamhm, gamh2p, yh2, yh;
 FILE *outfile;

 yh2=0.5*h2frac;
 yh=nHI/dens;

 rate=(5.5e-29/temp)*pow(nHI, 3.0) + (1.0/8.0)*(5.5e-29/temp)*(h2frac*0.5*dens)*nHI*nHI;
 en=7.17775e-12;
 gamHR=en*rate;

 k8=4.e-9*pow(temp, -0.17);
 k10=6.0e-10;

ncr=1.e6*pow(temp, -0.5)*pow(1.4*yh2*exp(-12000.0/(temp+1200.0)) + 1.6*yh*exp(-pow(400.0/temp,2.0)),-1.0);

gamhm=k8*nhm*nHI*3.53*1.6e-12/(1.0+(ncr/nHI));
gamh2p=k10*nh2p*nHI*1.83*1.6e-12/(1.0+(ncr/nHI));

gamHR = gamHR + gamhm + gamh2p;

 *pgamHR=-gamHR;
}


//AS has added term to account for excitation of H2 by other H2 molecules at high densities
/*======================================================================
**evaluates extra H2 cooling term
========================================================================*/

void hrcool(double dens, double temp, double Tcmb, double time, double h2frac, double nHI, double *plamHR)
 {
 double LH2,LLTE,XNCRIT,LLOW,T3,lamHR;
 double LOGLOW,LRLTE,LVLTE;
 double T3C,LOGLOWC,templogc,templog;
 double LLOWC,LRLTEC,LVLTEC,LLTEC,XNCRITC;
 double XNCRIT_ROT, XNCRIT_VIB, XNCRIT_ROTC, XNCRIT_VIBC;
 double LLOW_ROT, LLOW_VIB, LLOW_ROTC, LLOW_VIBC;
 double gamma2, gamma3, gamma10, gamma20, gamma2C, gamma3C, gamma10C, gamma20C, e0, e1, e2, e3, e10, e20;
 double k, c, h;
 double we, wexe, weye, weze, Bv, Dv;

 double nh2;
 int m2=0;
 FILE *outfile;

 T3=temp/1000.e0;
 T3C=Tcmb/1000.e0;
 templog=log10(temp);
 templogc=log10(Tcmb);

 nh2=0.5*h2frac*dens;

 k=1.38e-16;
 c=3.e10;
 h=6.625e-27;

 we=4395.2;
 wexe=117.9;
 weye=0.29;
 weze=0.046;
 Bv=60.809-(2.993*0.5);
 Dv=0.04656+(-1.799e-3*0.5);

 gamma2=(3.3e-12 + (6.6e-12*T3))*0.276*pow(2.0,2.0)*exp(-pow((2.0/3.18),1.7));
 gamma3=(3.3e-12 + (6.6e-12*T3))*0.276*pow(3.0,2.0)*exp(-pow((3.0/3.18),1.7));
 gamma10=1.4e-12*pow(temp, 0.5)*exp(-12000.0/(temp+1200.0));
 gamma20=0.0;

 gamma2C=(3.3e-12 + (6.6e-12*T3C))*0.276*pow(2.0,2.0)*exp(-pow((2.0/3.18),1.7));
 gamma3C=(3.3e-12 + (6.6e-12*T3C))*0.276*pow(3.0,2.0)*exp(-pow((3.0/3.18),1.7));
 gamma10C=1.4e-12*pow(Tcmb, 0.5)*exp(-12000.0/(Tcmb+1200.0));
 gamma20C=0.0;

 e0=(we*0.5) - (wexe*pow(0.5,2.0)) + (weye*pow(0.5,3.0)) - (weze*pow(0.5,4.0)); 
 e1=(we*0.5) - (wexe*pow(0.5,2.0)) + (weye*pow(0.5,3.0)) - (weze*pow(0.5,4.0)) + (Bv*1.0*(1.0+1.0)) - (Dv*pow(1.0, 2.0)*pow(1.0+1.0, 2.0));
 e2=(we*0.5) - (wexe*pow(0.5,2.0)) + (weye*pow(0.5,3.0)) - (weze*pow(0.5,4.0)) + (Bv*2.0*(2.0+1.0)) - (Dv*pow(2.0, 2.0)*pow(2.0+1.0, 2.0));
 e3=(we*0.5) - (wexe*pow(0.5,2.0)) + (weye*pow(0.5,3.0)) - (weze*pow(0.5,4.0)) + (Bv*3.0*(3.0+1.0)) - (Dv*pow(3.0, 2.0)*pow(3.0+1.0, 2.0));

 e0=e0*h*c;
 e1=e1*h*c;
 e2=e2*h*c;
 e3=e3*h*c;

 e10=k*5860.0;
 e20=2.0*k*5860.0;

 LLOW_ROT=0.25*(5.0*gamma2*exp(-(e2-e0)/(k*temp))*(e2-e0)) + 0.75*((7.0/3.0)*gamma3*exp(-(e3-e1)/(k*temp))*(e3-e1));
 LLOW_ROTC=0.25*(5.0*gamma2C*exp(-(e2-e0)/(k*Tcmb))*(e2-e0)) + 0.75*((7.0/3.0)*gamma3C*exp(-(e3-e1)/(k*Tcmb))*(e3-e1));

 LLOW_VIB=gamma10*exp(-e10/(k*temp))*e10 + gamma20*exp(-e20/(k*temp))*e20;
 LLOW_VIBC=gamma10C*exp(-e10/(k*Tcmb))*e10 + gamma20C*exp(-e20/(k*Tcmb))*e20;

 LRLTE=(9.5e-22*pow(T3,3.76)*exp(-pow(0.13/T3,3))/
        (1.e0+0.12*pow(T3,2.1))+3.e-24*exp(-0.51/T3))/nh2;
 LRLTEC=(9.e-22*pow(T3C,3.76)*exp(-pow(0.13/T3C,3))/
        (1.e0+0.12*pow(T3C,2.1))+3.e-24*exp(-0.51/T3C))/nh2;
 LVLTE=(6.7e-19*exp(-5.86/T3) + 1.6e-18*exp(-11.7/T3))/nh2;
 LVLTEC=(6.7e-19*exp(-5.86e0/T3C) + 
         1.6e-18*exp(-11.7e0/T3C))/nh2;
 LLTE=LRLTE + LVLTE;
 LLTEC=LRLTEC + LVLTEC;

 XNCRIT=nh2*(LLTE/LLOW);
 XNCRITC=nh2*(LLTEC/LLOWC);
 //LH2=LLTE/(1.e0+XNCRIT/nh2) - LLTEC/(1.e0+XNCRITC/nh2);

 XNCRIT_ROT = nh2*LRLTE/LLOW_ROT;
 XNCRIT_VIB = nh2*LVLTE/LLOW_VIB;
 XNCRIT_ROTC = nh2*LRLTEC/LLOW_ROTC;
 XNCRIT_VIBC = nh2*LVLTEC/LLOW_VIBC;
  
 LH2=LRLTE/(1.e0+XNCRIT_ROT/nh2) + LVLTE/(1.e0+XNCRIT_VIB/nh2) - LRLTEC/(1.e0+XNCRIT_ROTC/nh2) -  LVLTEC/(1.e0+XNCRIT_VIBC/nh2);

 lamHR= nh2*nh2*LH2;
 

 *plamHR=lamHR;
 
 }



void hlcool(double xn,double temp,double zred, double ny[9], double *plamHl)
/**======================================================================
*** Evaluates atomic line cooling due to H and He!
*** (see Cen 1992, ApJS, 78, 341)
***====================================================================*/
{
  double xnh,xnhe,LAM[11], lamHl;
  double gff,XH1,XH2,XH3,XH4,XH5,T3,T5,T6;
  double Tcmb;

  xnh=0.93e0*xn;
  xnhe=0.07e0*xn;
  Tcmb=2.7e0*(1.e0+zred);
      
/* Bremsstrahlung cooling  */
  if (temp > 5.e3) {
    gff=1.1e0+0.34e0*exp(-pow(5.5e0-log10(temp),2.e0)/3.e0);
    LAM[0]=1.42e-27*gff*sqrt(temp)*(ny[1]+ny[7]
                  +4.e0*ny[8])*ny[5];
  } else {
    LAM[0]=1.e-52;
  }

/* Collisional ionization cooling (HI) */
  T5=temp/1.e5;
  XH1=sqrt(temp)/(1.e0+sqrt(T5));
  if (temp > 5.e3) {
    LAM[1]=1.27e-21*XH1*exp(-157809.1e0/temp)*ny[0]*ny[5];
  } else {
    LAM[1]=1.e-52;
  }

/* Collisional ionization cooling (HeI) */
  if (temp > 8.e3) {
    LAM[2]=9.38e-22*XH1*exp(-285335.4e0/temp)*ny[6]*ny[5];
  } else {
    LAM[2]=1.e-52;
  }

/* Collisional ionization cooling (HeII) */
  if (temp > 1.e4) {
    LAM[3]=4.95e-22*XH1*exp(-631515.0e0/temp)*ny[7]*ny[5];
  } else {
    LAM[3]=1.e-52;
  }

/* Recombination cooling (HII) */
  T3=temp/1.e3;
  T6=temp/1.e6;
  XH4=1.e0/(1.e0+pow(T6,0.7e0));
  if (temp > 5.e3) {
    LAM[4]=8.70e-27*XH4*sqrt(temp)*pow(T3,-0.2e0)*ny[1]*ny[5];
  } else {
    LAM[4]=1.e-52;
  }

/* Recombination cooling (HeII) */
  if (temp > 5.e3) {
    LAM[5]=1.55e-26*pow(temp,0.3647e0)*ny[7]*ny[5];
  } else {
    LAM[5]=1.e-52;
  }

/* Recombination cooling (HeIII) */
  if (temp > 5.e3) {
    LAM[6]=3.48e-26*XH4*sqrt(temp)*pow(T3,-0.2e0)*ny[8]*ny[5];
  } else {
    LAM[6]=1.e-52;
  }

/* Dielectronic recombination cooling  */
  XH5=exp(-470000.e0/temp)*(1.e0+0.3e0*exp(-94000.e0/temp));
  if (temp > 5.e3) {
    LAM[7]=1.24e-13*XH5*pow(temp,-1.5e0)*ny[7]*ny[5];
  } else {
    LAM[7]=1.e-52;
  }

/* Collisional excitation cooling (HI) */
  XH2=1.e0/(1.e0+sqrt(T5));
  if (temp > 5.e3) {
    LAM[8]=7.50e-19*XH2*exp(-118348.e0/temp)*ny[0]*ny[5];
  } else {
    LAM[8]=1.e-52;
  }

/* Collisional excitation cooling (HeII) */
  XH3=XH2*pow(temp,-0.397e0);
  if (temp > 8.e3) {
    LAM[9]=5.54e-17*XH3*exp(-473638.e0/temp)*ny[7]*ny[5];
  } else {
    LAM[9]=1.e-52;
  }

/* Compton cooling  */
  LAM[10]=5.4e-36*pow(1.e0+zred,4.e0)*ny[5]*(temp-Tcmb);

  lamHl=LAM[0]+LAM[1]+LAM[2]+LAM[3]+LAM[4]+LAM[5]+LAM[6]+
       LAM[7]+LAM[8]+LAM[9]+LAM[10];
  //lamHl=lamHl/(xnh*xnh);
  *plamHl=lamHl;

}
