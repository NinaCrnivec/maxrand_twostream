#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "twomaxrnd.h"


/*
static int twostream_maxrand (double *dtau_c, double *omega0_c, double *g_c, // cloudy region parameters
			      double *dtau_f, double *omega0_f, double *g_f, // free region parameters
			      double *cf, int nlev, 
			      double S0, double mu0, double Ag, 
			      double Bg, double *B, int delta,
			      double **Edir, double **Edn, double **Eup, double **Lavg);
*/


//int twomaxrnd (float *dtau_org, float *omega0_org, float *g_org,
//	       double *dtau_clr_org, double *omega0_clr_org, double *g_clr_org, // NEW
//	       float *cf, int nlev, 
//	       double S0, double mu0,
//	       double Ag,
//	       int planck,
//	       int delta,
//	       int nzout,
//	       float *zd,
//	       float *temper,
//	       float btemp,
//	       float wvnmlo,
//	       float wvnmhi,
//	       float *zout, 
//	       float *fldn, 
//	       float *flup, 
//	       float *fldir,
//	       float *uavg);
                      

int main (int argc, char **argv)
{

int iStatus;
int ilev;
int ilyr;
int nlev=3; // number of vertical levels in the calculation
int nlyr=nlev-1;

float *dtau_c = calloc(nlyr,sizeof(float));
float *ssa_c = calloc(nlyr,sizeof(float));
float *g_c = calloc(nlyr,sizeof(float));   
double *dtau_f = calloc(nlyr,sizeof(double));
double *ssa_f = calloc(nlyr,sizeof(double));
double *g_f = calloc(nlyr,sizeof(double));
float *cf = calloc(nlyr,sizeof(float)); 
double S0;   // extraterrestrial irradiance [W/m2]
double mu0;  // mu0=cos(sza); 
double Ag;   // ground albedo
int planck;
int delta;  // Flag for delta-scaling; 1=yes; 0=no;
int nzout;
float *zd = calloc(nlev,sizeof(float));
float *temper = calloc(nlev,sizeof(float));
float btemp;
float wvnmlo;
float wvnmhi;
float *zout = calloc(nlev,sizeof(float));
float *Edw = calloc(nlev,sizeof(float));
float *Eup = calloc(nlev,sizeof(float));
float *Edir = calloc(nlev,sizeof(float));
float *uavg = calloc(nlev,sizeof(float));

double tau_c_total;
double tau_f_total;

// Set constants and flags here:
S0=1000.0; 
mu0=0.5; 
Ag=0.1;   
planck=0;
delta=0;
nzout=nlev;
btemp = 1.0; // izmisljena vrednost
wvnmlo = 1.0; // izmisljena vrednost
wvnmhi = 2.0; // izmisljena vrednost

tau_c_total=10.0;
tau_f_total=0.5;

//cf[0] = 0.5;
//cf[1] = 0.0;

for (ilyr=0; ilyr<nlyr; ilyr++){
   dtau_c[ilyr] = tau_c_total/nlyr;
   dtau_f[ilyr] = tau_f_total/nlyr;
   ssa_c[ilyr] = 0.75;
   ssa_f[ilyr] = 0.75;
   g_c[ilyr] = 0.90;
   g_f[ilyr] = 0.10;
   cf[ilyr] = 0.5;
}//e-for

// Example Input for Eddington coeff:
//dtau = 2.0
//omega0 = 0.75
//g = 0.88

// Use this loop if you want to set:
// Optical properties(cloud part) = optical properties(cloud-free part):
/*
for(ilyr=0; ilyr<nlyr; ilyr++){
  dtau_c[ilyr] = dtau_f[ilyr];
  ssa_c[ilyr] = ssa_f[ilyr];
  g_c[ilyr] = g_f[ilyr];
}//e-for
*/

temper[0] = 1.0;
zd[0] = 1.0; 
zout[0] = 1.0; 

for(ilev=1; ilev<nlev; ilev++){
   temper[ilev] = temper[ilev-1] + 1.0;
   zd[ilev] = zd[ilev-1] + 1.0; 
   zout[ilev] = zout[ilev-1] + 1.0; 
}//e-for


fprintf (stderr, "\n");  
fprintf (stderr, "Input parameters:\n");
fprintf (stderr, "\n");  

fprintf (stderr, "nlev = %d \n", nlev);  
fprintf (stderr, "nlyr = %d \n", nlyr);  
fprintf (stderr, "\n");  

fprintf (stderr, "S0 = %.2f \n", S0); 
fprintf (stderr, "mu0 = %.2f \n", mu0); 
fprintf (stderr, "Ag = %.2f \n", Ag); 
fprintf (stderr, "delta = %d \n", delta);  
fprintf (stderr, "planck = %d \n", planck);  
fprintf (stderr, "nzout = %d \n", nzout);  
fprintf (stderr, "btemp = %.2f \n", btemp);
fprintf (stderr, "wvnmlo = %.2f \n", wvnmlo);
fprintf (stderr, "wvnmhi = %.2f \n", wvnmhi);
fprintf (stderr, "\n");  

fprintf (stderr, "Vertical profiles of layer quantities:\n");
fprintf (stderr, "\n");  

for(ilyr=0; ilyr<nlyr; ilyr++){
   fprintf(stderr, "ilyr = %d, cf = %.4f\n", ilyr, cf[ilyr]);
   fprintf(stderr, "ilyr = %d, dtau_c = %.4f, dtau_f = %.4f\n", ilyr, dtau_c[ilyr], dtau_f[ilyr]);
   fprintf(stderr, "ilyr = %d, ssa_c = %.4f, ssa_f = %.4f\n", ilyr, ssa_c[ilyr], ssa_f[ilyr]);
   fprintf(stderr, "ilyr = %d, g_c = %.4f, g_f = %.4f\n", ilyr, g_c[ilyr], g_f[ilyr]);
}//e-for
fprintf (stderr, "\n");

fprintf (stderr, "Vertical profiles of level quantities:\n");
fprintf (stderr, "ilev, temper[ilev], zd[ilev], zout[ilev]\n");

for(ilev=0; ilev<nlev; ilev++){
   fprintf(stderr, "%d %.4f %.4f %.4f\n", ilev, temper[ilev], zd[ilev], zout[ilev]);
}//e-for
fprintf (stderr, "\n");


fprintf (stderr, "\n");  
fprintf (stderr, "Entering function twomaxrnd!\n");
fprintf (stderr, "\n"); 

iStatus = twomaxrnd (dtau_c, ssa_c, g_c, dtau_f, ssa_f, g_f, 
		    cf, nlev, S0, mu0, Ag, planck, delta, nzout, zd, temper, 
		    btemp, wvnmlo, wvnmhi, zout, 
		    Edw, Eup, Edir, uavg);
	        
if (iStatus!=0) {
  fprintf (stderr, "Error %d returned by twomaxrnd()\n", iStatus);
  return iStatus;    
}
 
 
fprintf (stderr, "\n");  
fprintf (stderr, "Entering main again...\n");
fprintf (stderr, "\n");  
 

// Print out the results;
fprintf (stderr, "Final results:\n");
fprintf (stderr, "Vertical profiles of total irradiances [W/m2]:\n");

for (ilev=0; ilev<nlev; ilev++){
   fprintf (stderr, "ilev = %04d, Edir = %10.4lf, Edw = %10.4lf, Eup = %10.4lf\n", ilev, Edir[ilev], Edw[ilev], Eup[ilev]); 
}//e-for

fprintf (stderr, "\n");  
fprintf (stderr, "The END\n");
fprintf (stderr, "\n");  




return 0;
}










