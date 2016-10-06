#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss.h"




//int solve_gauss (double **A, double *b, int n, double **res);


int calcp1p2p3p4 (int nlev, double *cf, double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4);


void delta_scale_hg (double tau, double ssa, double g, 
                    double *tauscale, double *ssascale, double *gscale);


void eddington_coeffc (double dtau, double g, double omega0, double mu0, double *a11, double *a12,
                       double *a13, double *a23, double *a33);


int buildMatrixA (int nlev, double Ag, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                  double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, double **matrixA);


int makeMatrixA2 (int nlev, double **matrixA);


int buildVectorB (int nlev, double Ag, double mu0, double *ar_a13_c, double *ar_a13_f, 
		  double *ar_a23_c, double *ar_a23_f, double *ar_S_c, double *ar_S_f, 
		  double *ar_p1, double *ar_p3, double *vectB);


int makeVectorMinusB (int nlev, double *vectB);

void displayMatrix (int nlev, double **matrixA, char *name);

void displayVector (int nlev, double *vect, char *name);

void freeMemory (int nlev, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                 double *ar_a13_c, double *ar_a13_f, double *ar_a23_c, double *ar_a23_f, 
                 double *ar_a33_c, double *ar_a33_f, 
                 double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, 
                 double *S_c, double *S_f, 
                 double *bb, double *xx, double **AA);



/* Nina's Delta-two-stream method with maximum-random overlap assumption for partial cloudiness;


/*
c = cloudy region;
f = free region (= cloud-free region);
*/


/*
The form of twostream_maxrnd as requested by Bernhard is:
int twostream_maxrnd (double *dtau_org_c, double *omega0_org_c, double *g_org_c,// cloudy region parameters
                      double *dtau_org_f, double *omega0_org_f, double *g_org_f,// free region parameters
                      double *cf, int nlev, double S0, double mu0, double Ag, int delta,
                      double *Edir, double *Edn, double *Eup) 
  
*/


/*
int nina_maxrand(int nlyr,
    double *dtauf, double *w0f,double *gf,
    double *dtauc, double *w0c,double *gc,
    double *cfrac,
    double mu0, double incSolar, double albedo,
    double *Sf, double *Ednf, double*Eupf,
    double *Sc, double *Ednc, double*Eupc) 
*/


int nina_maxrand(int nlyr,
    double *dtau_org_f, double *omega0_org_f, double *g_org_f,// free region parameters
    double *dtau_org_c, double *omega0_org_c, double *g_org_c,// cloudy region parameters
    double *cf,
    double mu0, double S0, double Ag,
    double *Edir_f, double *Edn_f, double *Eup_f,
    double *Edir_c, double *Edn_c, double *Eup_c)

{
  int nlev=nlyr+1;
  int delta=0;
  
  printf("Computing Nina's Maxrand with %d layers \n", nlyr);

  int ilev;
  int ilyr;
  int iStatus;
  
  // Eddington coefficients for cloudy regions:
  double a11_c;
  double a12_c;
  double a13_c;
  double a23_c;
  double a33_c;
    
  // Eddington coefficients for clear-sky regions:
  double a11_f;
  double a12_f;
  double a13_f;
  double a23_f;
  double a33_f;
  

  // CREATE ARRAYS FOR VERTICAL PROFILES OF EDDINGTON COEFFICIENTS:
  
  // Eddington coefficients for cloudy regions:
  double *ar_a11_c=calloc(nlyr,sizeof(double));
  double *ar_a12_c=calloc(nlyr,sizeof(double));
  double *ar_a13_c=calloc(nlyr,sizeof(double));
  double *ar_a23_c=calloc(nlyr,sizeof(double));
  double *ar_a33_c=calloc(nlyr,sizeof(double));
  
  // Eddington coefficients for clear-sky regions:
  double *ar_a11_f=calloc(nlyr,sizeof(double));
  double *ar_a12_f=calloc(nlyr,sizeof(double));
  double *ar_a13_f=calloc(nlyr,sizeof(double));
  double *ar_a23_f=calloc(nlyr,sizeof(double));
  double *ar_a33_f=calloc(nlyr,sizeof(double));  

  // Parameters related to cloud cover of two contiguous layers:
  // in Zdunkowski (pages: 180-183) denoted as: b1,b2,b3,b4;
  // here: p1, p2, p3, p4;   
  double *ar_p1=calloc(nlyr,sizeof(double));
  double *ar_p2=calloc(nlyr,sizeof(double));
  double *ar_p3=calloc(nlyr,sizeof(double));
  double *ar_p4=calloc(nlyr,sizeof(double));
  
  double dtau_c;
  double omega0_c;
  double g_c;
  
  double dtau_f;
  double omega0_f;
  double g_f;

  double **AA=NULL;
  double *bb;
  double *xx; // result vector
  
  double *S_c;  
  double *S_f;
  
  /*
  double *Edir_c;
  double *Edir_f;
  double *Eup_c;
  double *Eup_f;
  double *Edn_c;
  double *Edn_f;
  */
  
  int i, j;
  double value;
  

  /* allocate memory for equation system */
  AA = calloc (4*nlev, sizeof(double*));
   
  for(ilev=0;ilev<4*nlev;ilev++) {
    if ((AA[ilev] = calloc (4*nlev, sizeof(double)))==NULL) {
      fprintf (stderr, "Error allocating memory for AA[%d]\n", ilev);
      return -1;
    }//e-if
  }//e-for
  
  bb = calloc (4*nlev, sizeof(double)); 
  xx = calloc (4*nlev, sizeof(double)); 
  
  S_c = calloc (nlev, sizeof(double));
  S_f = calloc (nlev, sizeof(double));
   
  
  // c = cloudy regions, f = free regions
  /*
  Edir_c = calloc (nlev, sizeof(double));
  Edir_f = calloc (nlev, sizeof(double));
  Eup_c = calloc (nlev, sizeof(double));
  Eup_f = calloc (nlev, sizeof(double));
  Edn_c = calloc (nlev, sizeof(double));
  Edn_f = calloc (nlev, sizeof(double));
  */
 
 
  for (ilyr=0; ilyr<nlyr; ilyr++){
    fprintf (stderr, "CF %d %f\n", ilyr, cf[ilyr]);
  }//e-for
  fprintf (stderr, "\n");
 
  // Calculate vertical profiles of p1, p2, pr3 and p4 from vertical profile of partial cloud cover;
  // INPUT: nlev, cf;
  // OUTPUT: ar_p1, ar_p2, ar_p3, ar_p4
  iStatus = calcp1p2p3p4(nlev, cf, ar_p1, ar_p2, ar_p3, ar_p4);
  if (iStatus != 0)  {
    fprintf (stderr, "Error calculating vertical profiles of p1, p2, p3 and p4 from cloud cover; ERROR=%d \n", iStatus);
    return -1;
  }//e-if
  

  // Calculate vertical profiles of Eddington coefficients
  for(ilyr=0;ilyr<nlyr;ilyr++){
    
    /* delta scaling of optical properties for cloudy regions */
    if (delta) {
      delta_scale_hg (dtau_org_c[ilyr], omega0_org_c[ilyr], g_org_c[ilyr], 
		      &dtau_c, &omega0_c, &g_c);                    
    }else{
      dtau_c   = dtau_org_c[ilyr];
      omega0_c = omega0_org_c[ilyr];
      g_c      = g_org_c[ilyr];
    }//e-if
    
    /* Clear-sky optical properties */
    dtau_f   = dtau_org_f[ilyr];
    omega0_f = omega0_org_f[ilyr];
    g_f      = g_org_f[ilyr];  
    
   
   /* omega0 should not be 1 (avoiding singularity problem) */
   /* restrict omega0 to 0.999999 */
   if (omega0_c > 0.999999) omega0_c = 0.999999;
   if (omega0_f > 0.999999) omega0_f = 0.999999;

   
    /* Calculate Eddington coefficients for a given layer */ 
    // Cloudy:
    eddington_coeffc (dtau_c, g_c, omega0_c, mu0, &a11_c, &a12_c, &a13_c, &a23_c, &a33_c);
    //Clear-sky:
    eddington_coeffc (dtau_f, g_f, omega0_f, mu0, &a11_f, &a12_f, &a13_f, &a23_f, &a33_f);
    
    // Save Eddington coefficients for each layer to arrays:
    ar_a11_c[ilyr] = a11_c;
    ar_a11_f[ilyr] = a11_f;
    
    ar_a12_c[ilyr] = a12_c;
    ar_a12_f[ilyr] = a12_f;
    
    ar_a13_c[ilyr] = a13_c;
    ar_a13_f[ilyr] = a13_f;
    
    ar_a23_c[ilyr] = a23_c;
    ar_a23_f[ilyr] = a23_f;
      
    ar_a33_c[ilyr] = a33_c;
    ar_a33_f[ilyr] = a33_f;
  }//e-for over ilyr;

 
  // Step #A: initialize vectors S_c and S_f
  // S[0]= S0 = S_c[0] + S_f[0] = cf[0]*S0 + (1.0-cf[0])*S0;
  S_c[0]=cf[0]*S0;    
  S_f[0]=(1.0-cf[0])*S0;
   
  for(ilev=1;ilev<nlev;ilev++){
     S_c[ilev]=ar_a33_c[ilev-1]*((1.0-ar_p1[ilev-1])*S_f[ilev-1] + ar_p3[ilev-1]*S_c[ilev-1]); 
     S_f[ilev]=ar_a33_f[ilev-1]*(ar_p1[ilev-1]*S_f[ilev-1] + (1.0-ar_p3[ilev-1])*S_c[ilev-1]); 
  }//e-for
  
  
  // Equation system has the form: xx = AA*xx + bb;
  // Step #B: build matrix AA; 
  // INPUT: nlev, Ag, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, ar_p1, ar_p2, ar_p3, ar_p4;
  // OUTPUT: AA
  iStatus = buildMatrixA (nlev, Ag, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, 
                          ar_p1, ar_p2, ar_p3, ar_p4, AA);
  if (iStatus != 0){
    fprintf (stderr, "buildMatrixA ERROR=%d \n", iStatus);
    //freeMemory (...);
    return -2;
  }//e-if
  
  
  /*
  // Display the nonzero elements of matrix AA to check;
  for (i=0; i<4*nlev; i++){
    for (j=0; j<4*nlev; j++){
       value = AA[i][j];
       if (value > 0.0){
          //fprintf (stderr, "i=%d, j=%d, AA[%2d][%2d] = %f\n", i,j,i,j,(AA[i][j]));
       }//e-if
    }//e-for
  }//e-for
  */
 
  // Display matrix AA to check;
  //displayMatrix (nlev, AA, "A");
  
  
 // Step #C: make matrix A2=AA-IdentityMatrix=AA-II, save in to the same matrix AA;
 // Since the equation system in the form A*x=b is needed for function solve_gauss;
 // xx = AA*xx + bb; AA*xx - xx + bb = 0; (AA-II)*xx + bb = 0;
 // FINAL FORM OF THE SYSTEM SENT IN solve_gauss: (AA-II)*xx = -bb;
  iStatus = makeMatrixA2 (nlev, AA);
  if (iStatus != 0){
    fprintf (stderr, "makeMatrixA2 ERROR=%d \n", iStatus);
    //freeMemory (...);
    return -3;
  }//e-if
 
  // Display matrix A2=AA-II to check;
  //displayMatrix (nlev, AA, "A-I");

  
  // Step #D: build vector bb
  // INPUT: Ag, mu0 (?), ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, S_c, S_f, ar_p1, ar_p3; 
  // OUTPUT: bb
  iStatus = buildVectorB (nlev, Ag, mu0, ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, S_c, S_f, ar_p1, ar_p3, bb);
  if (iStatus != 0){
    fprintf (stderr, "buildVectorB ERROR=%d \n", iStatus);
    //freeMemory (...);
    return -4;
  }//e-if
  
 
  // Display vector B
  //displayVector (nlev, bb, "B");
  
  // STEP #E: Make vector -bb and save it to the same vector bb;
  iStatus = makeVectorMinusB (nlev, bb);
  if (iStatus != 0){
    fprintf (stderr, "makeVectorMinusB ERROR=%d \n", iStatus);
    //freeMemory (...);
    return -5;
  }//e-if

  // Display vector -B
  //displayVector (nlev, bb, "-B");
  
  
  // Step #F: solve AA*xx=bb --> xx = ...
  // AA = quadratic matrix
 //int solve_gauss (double **A, double *b, int n, double **res)
  iStatus = solve_gauss (AA, bb, 4*nlev, &xx);
  if (iStatus != 0){
    fprintf (stderr, "Error %d solving equation system\n", iStatus);
    //freeMemory (...);
    return -6;
  }//e-if
 
  // Display vector X (=xx):
  //displayVector (nlev, xx, "X");
  

 /*
Form of vector X (=xx) for certain ilev:
Eup_f(ilev)
Eup_c(ilev)
Edn_f(ilev)
Edn_c(ilev) 
*/
 
  for(ilev=0;ilev<nlev;ilev++){
     Eup_f[ilev] = xx[4*ilev];
     Eup_c[ilev] = xx[4*ilev+1];
     Edn_f[ilev] = xx[4*ilev+2];
     Edn_c[ilev] = xx[4*ilev+3];
  }//e-for

  
  for(ilev=0;ilev<nlev;ilev++){
    Edir_c[ilev]=S_c[ilev]*mu0;  
    Edir_f[ilev]=S_f[ilev]*mu0;  
  }//e-for

  
  /* Print out irradiances for cloudy and cloud-free regions */
  fprintf(stderr, "\n");
  fprintf(stderr, "Vertical profiles of irradiances [W/m2] for cloudy and cloud-free regions:\n");
  fprintf(stderr, "ilev,    Edir_c,    Edir_f,     Edn_c,    Edn_f,     Eup_c,     Eup_f\n");
  for (ilev=0;ilev<nlev;ilev++){
     fprintf(stderr, "%04d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",ilev, Edir_c[ilev], Edir_f[ilev], Edn_c[ilev], Edn_f[ilev], Eup_c[ilev], Eup_f[ilev]);    
  }//e-for
  
  
  /* SUM UP THE IRRADIANCES FOR CLOUDY AND CLOUD-FREE REGIONS TO OBTAIN THE FINAL RESULT */
  for (ilev=0;ilev<nlev;ilev++){
    //Edir[ilev] = Edir_c[ilev] + Edir_f[ilev];
    //Eup[ilev] = Eup_c[ilev] + Eup_f.1[ilev];  
    //Edn[ilev] = Edn_c[ilev] + Edn_f[ilev];  
  }//e-for  
    
    
  // Free allocated memory before return 0
  freeMemory (nlev, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, 
              ar_a33_c, ar_a33_f, ar_p1, ar_p2, ar_p3, ar_p4, 
              S_c, S_f, bb, xx, AA);

  return 0;
}//e-nina_maxrand





/*  
From Zdunkowski (pages 180-183): 
p1 = (1 - max(cf[i], cf[i-1])) / (1 - cf[i-1])  
p2 = (1 - max(cf[i], cf[i+1]) / (1 - cf[i+1])  
p3 = min(cf[i], cf[i-1]) / cf[i-1]  
p4 = min(cf[i], cf[i+1]) / cf[i+1]
   
REMARK from Zdunkowski:
It is noteworthy that a particular coefficient pj (j=1,2,3,4) is set equal to 1 
if an undetermined expression 0/0 occurs.
This follows from physical reasoning or from applying l'Hopital's rule.
*/


//======================
// FUNCTION calcp1p2p3p4
//======================
// INPUT: nlev, cf;
// OUTPUT: ar_p1, ar_p2, ar_p3, ar_p4;
int calcp1p2p3p4 (int nlev, double *cf, double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4)
{
int ilyr;

int nlyr;

double max_p1; // = max(cf[ilyr],cf[ilyr-1]); appears in the expression for p1;
double min_p3; // = min(cf[ilyr],cf[ilyr-1]); appears in the expression for p3;
double max_p2; // = max(cf[ilyr],cf[ilyr+1]); appears in the expression for p2;
double min_p4; // = min(cf[ilyr],cf[ilyr+1]); appears in the expression for p4;

nlyr=nlev-1;

/*
FORMULAS:
ar_p1[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr-1])) / (1.0 - cf[ilyr-1]);
ar_p2[ilyr] = (1.0 - max(cf[ilyr], cf[ilyr+1]) / (1.0 - cf[ilyr+1]); 
ar_p3[ilyr] = min(cf[ilyr], cf[ilyr-1]) / cf[ilyr-1]; 
ar_p4[ilyr] = min(cf[ilyr], cf[ilyr+1]) / cf[ilyr+1];   
*/


// Calculate vertical profiles of p1, p3
// Special case: ilyr=0;
for (ilyr=0;ilyr<nlyr;ilyr++) {
    
  if (ilyr == 0){
    ar_p1[ilyr] = 1.0; 
    ar_p3[ilyr] = 1.0; 
  }else{  
      
    // Find max_p1:
    if (cf[ilyr] > cf[ilyr-1]){
      max_p1 = cf[ilyr];
    }else{
      max_p1 = cf[ilyr-1];
    }//e-if
    
    // Find_min_p3:
    if (cf[ilyr] < cf[ilyr-1]){
      min_p3 = cf[ilyr];
    }else{
      min_p3 = cf[ilyr-1];
    }//e-if
    
    if (cf[ilyr-1] == 1.0) ar_p1[ilyr] = 1.0;
    else ar_p1[ilyr] = (1.0 - max_p1) / (1.0 - cf[ilyr-1]);
    
    if (cf[ilyr-1] == 0.0) ar_p3[ilyr] = 1.0;
    else ar_p3[ilyr] = min_p3 / cf[ilyr-1];  
    
  }//e-if  
}//e-for over ilyr;
 
 
// Calculate vertical profiles of p2 and p4
// Special case: ilyr = nlyr-1 = nlev-2;
for (ilyr=0;ilyr<nlyr;ilyr++) {
  
  if (ilyr == (nlyr-1)){ 
    ar_p2[ilyr] = 1.0; 
    ar_p4[ilyr] = 1.0; 
  }else{ 
    
    // Find max_p2:
    if (cf[ilyr] > cf[ilyr+1]){
      max_p2 = cf[ilyr];
    }else{
      max_p2 = cf[ilyr+1];
    }//e-if
    
    // Find_min_p4:
    if (cf[ilyr] < cf[ilyr+1]){
      min_p4 = cf[ilyr];
    }else{
      min_p4 = cf[ilyr+1];
    }//e-if
    
    if (cf[ilyr+1] == 1.0) ar_p2[ilyr] = 1.0;
    else ar_p2[ilyr] = (1.0 - max_p2) / (1.0 - cf[ilyr+1]);  
    
    if (cf[ilyr+1] == 0.0) ar_p4[ilyr] = 1.0;
    else ar_p4[ilyr] = min_p4 / cf[ilyr+1];
  }//e-if
}//e-for over ilyr;


// Print to terminal for test;
printf("\n");
printf("nlev = %d\n", nlev);
printf("nlyr = %d\n", nlyr);
printf("\n");
printf("Coefficients p1,p2,p3,p4\n");
printf("\n");
printf("ilyr, cf[ilyr], ar_p1[ilyr],ar_p3[ilyr], ar_p2[ilyr], ar_p4[ilyr]\n");
for (ilyr=0; ilyr<nlyr; ilyr++){
   printf("%d %8.2f %8.2f %8.2f %8.2f %8.2f\n",ilyr,cf[ilyr],ar_p1[ilyr],ar_p3[ilyr],ar_p2[ilyr],ar_p4[ilyr]);
}//e-for
printf("\n");

return 0;
}//e-calcp1p2p3p4






//=========================
// FUNCTION: delta_scale_hg
//=========================
void delta_scale_hg (double tau, double ssa, double g, 
                     double *tauscale, double *ssascale, double *gscale)
{
double f=g*g;
  
*tauscale = (1.0-ssa*f)*tau;
*ssascale = (1.0-f)*ssa/(1.0-ssa*f);
*gscale   = (g-f)/(1.0-f);
}//e-delta_scale_hg



//===========================
// FUNCTION: eddington_coeffc
//===========================
/* calculate Eddington coefficients a11, a12, a13, a23, and a33 from     */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_coeffc (double dtau, double g, double omega0, double mu0,
		       double *a11, double *a12, double *a13, double *a23, double *a33)
{
double alpha1=0, alpha2=0, alpha3=0, alpha4=0, alpha5=0, alpha6=0;
double lambda=0, b=0, A=0;
double denom=0;

alpha1 = (1.0-omega0)+0.75*(1.0-omega0*g);
alpha2 = -(1.0-omega0)+0.75*(1.0-omega0*g);
  
lambda=sqrt(alpha1*alpha1-alpha2*alpha2);

A=1.0/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));
  
*a11=A*2.0*lambda/alpha2;
*a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));
  
b=0.5-0.75*g*mu0;
alpha3=-omega0*b; 
alpha4=omega0*(1-b);
  
denom = (1.0/mu0/mu0-lambda*lambda);
alpha5=((alpha1-1.0/mu0)*alpha3-alpha2*alpha4)/denom;
alpha6=(alpha2*alpha3-(alpha1+1.0/mu0)*alpha4)/denom;
  
*a33=exp(-dtau/mu0);
  
*a13=alpha5*(1.0-(*a11)*(*a33))-alpha6*(*a12);
*a23=-(*a12)*alpha5*(*a33)+alpha6*((*a33)-(*a11));

}//e-eddington_coeffc




//=======================
// FUNCTION: buildMatrixA
//=======================
// INPUT: nlev, Ag, ar_a11_c, ar_a11_f, ar_a12_c, ar_a12_f, ar_p1, ar_p2, ar_p3, ar_p4;
// OUTPUT: matrixA;

int buildMatrixA (int nlev, double Ag, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                  double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, double **matrixA)
{
int i;     // position in levels (e.g. for nlev=21, i in range: 0 - 20)
int iRow;  // row position in matrix A (e.g. for nlev=21, iRow in range: 0 - 83)
   
int nextColEup; // index of column for next Eup data
int nextColEdw; // index of column for next Edw data; nextColEdw = nextColEup - 4   

// Set values for initial four rows of matrix A:
// First row:
matrixA[0][2] = ar_a12_f[0]*ar_p1[0];
matrixA[0][3] = ar_a12_f[0]*(1.0-ar_p3[0]);
matrixA[0][4] = ar_a11_f[0]*ar_p2[0];
matrixA[0][5] = ar_a11_f[0]*(1.0-ar_p4[0]);
// Second row:
matrixA[1][2] = ar_a12_c[0]*(1.0-ar_p1[0]);
matrixA[1][3] = ar_a12_c[0]*ar_p3[0];
matrixA[1][4] = ar_a11_c[0]*(1.0-ar_p2[0]);
matrixA[1][5] = ar_a11_c[0]*ar_p4[0];
// Third row is already zero; (needs to be zero due to upper boundary condition);
// Forth row is already zero; (needs to be zero due to upper boundary condition);  


//Print to the terminal for test:
//fprintf (stderr, "A[0][3]=%f\n", (matrixA[0][3]));
//fprintf (stderr, "ar_a12_f[0]=%f\n", (ar_a12_f[0]));
//fprintf (stderr, "ar_p3[0]=%f\n", (ar_p3[0]));
//fprintf (stderr, "(1.0-ar_p3[0])=%f\n", (1.0-(ar_p3[0])));


nextColEup = 6; // index of column for Eup(level1)
nextColEdw = 2; // index of column for Edw(level1)

for(i=1;i<nlev;i++){	
   iRow = 4*i;
   
   if (i == (nlev-1)){ // lower boundary condition for the forth and third row from bottom up of matrix A;
      matrixA[iRow][nextColEup] = Ag;
      matrixA[iRow+1][nextColEup+1] = Ag;
   }else{
     
      matrixA[iRow][nextColEup]   = ar_a12_f[i]*ar_p1[i];		
      matrixA[iRow][nextColEup+1] = ar_a12_f[i]*(1.0-ar_p3[i]);
      matrixA[iRow][nextColEup+2] = ar_a11_f[i]*ar_p2[i];
      matrixA[iRow][nextColEup+3] = ar_a11_f[i]*(1.0-ar_p4[i]);
      
      matrixA[iRow+1][nextColEup]   = ar_a12_c[i]*(1.0-ar_p1[i]);
      matrixA[iRow+1][nextColEup+1] = ar_a12_c[i]*ar_p3[i];
      matrixA[iRow+1][nextColEup+2] = ar_a11_c[i]*(1.0-ar_p2[i]);
      matrixA[iRow+1][nextColEup+3] = ar_a11_c[i]*ar_p4[i];
   }//e-if
   nextColEup += 4;

   matrixA[iRow+2][nextColEdw]   = ar_a11_f[i-1]*ar_p1[i-1];
   matrixA[iRow+2][nextColEdw+1] = ar_a11_f[i-1]*(1.0-ar_p3[i-1]);
   matrixA[iRow+2][nextColEdw+2] = ar_a12_f[i-1]*ar_p2[i-1];
   matrixA[iRow+2][nextColEdw+3] = ar_a12_f[i-1]*(1.0-ar_p4[i-1]);
   
   matrixA[iRow+3][nextColEdw]   = ar_a11_c[i-1]*(1.0-ar_p1[i-1]);	
   matrixA[iRow+3][nextColEdw+1] = ar_a11_c[i-1]*ar_p3[i-1];
   matrixA[iRow+3][nextColEdw+2] = ar_a12_c[i-1]*(1.0-ar_p2[i-1]);
   matrixA[iRow+3][nextColEdw+3] = ar_a12_c[i-1]*ar_p4[i-1];
   
   nextColEdw += 4;
  
}//e-for
          
return 0;
}//e-buildMatrixA



//======================
// FUNCTION makeMatrixA2
//======================
// Just write -1 to all diagonal elements of matrix A (in the same matrix)
int makeMatrixA2 (int nlev, double **matrixA)
{
int iRow;
for (iRow=0; iRow < 4*nlev; iRow++){
   matrixA[iRow][iRow] = -1.0;
}//e-for
return 0;
}//e-makeMatrixA2




//======================
// FUNCTION buildVectorB
//======================
// INPUT: nlev, Ag, mu0, ar_a13_c, ar_a13_f, ar_a23_c, ar_a23_f, ar_S_c, ar_S_f, ar_p1, ar_p3;
// OUTPUT: vectB;
int buildVectorB (int nlev, double Ag, double mu0, double *ar_a13_c, double *ar_a13_f, 
		  double *ar_a23_c, double *ar_a23_f, double *ar_S_c, double *ar_S_f, 
		  double *ar_p1, double *ar_p3, double *vectB)
{
int i;  // position in levels
int j;  // position in vector B

// Set initial four values: 
vectB[0] = ar_a13_f[0]*ar_S_f[0]; 
vectB[1] = ar_a13_c[0]*ar_S_c[0];
vectB[2] = 0.0; // upper boundary condition
vectB[3] = 0.0; // upper boundary condition


for (i=1; i<(nlev-1);i++){
   j = 4*i;     
   vectB[j]   = ar_a13_f[i]*(ar_p1[i]*ar_S_f[i] + (1.0-ar_p3[i])*ar_S_c[i]);
   vectB[j+1] = ar_a13_c[i]*((1.0-ar_p1[i])*ar_S_f[i] + ar_p3[i]*ar_S_c[i]);
   vectB[j+2] = ar_a23_f[i-1]*(ar_p1[i-1]*ar_S_f[i-1] + (1.0-ar_p3[i-1])*ar_S_c[i-1]);
   vectB[j+3] = ar_a23_c[i-1]*((1.0-ar_p1[i-1])*ar_S_f[i-1] + ar_p3[i-1]*ar_S_c[i-1]); 
}//e-for


// Treat last four values seperately:
vectB[4*nlev-4] = Ag*mu0*ar_S_f[nlev-1]; //Lower boundary condition;
vectB[4*nlev-3] = Ag*mu0*ar_S_c[nlev-1]; //Lower boundary condition;
vectB[4*nlev-2] = ar_a23_f[nlev-2]*(ar_p1[nlev-2]*ar_S_f[nlev-2] + (1.0-ar_p3[nlev-2])*ar_S_c[nlev-2]);
vectB[4*nlev-1] = ar_a23_c[nlev-2]*((1.0-ar_p1[nlev-2])*ar_S_f[nlev-2] + ar_p3[nlev-2]*ar_S_c[nlev-2]); 

return 0;
}//e-buildVectorB



//==========================
// FUNCTION makeVectorMinusB
//==========================
int makeVectorMinusB (int nlev, double *vectB)
{
int iRow;
for (iRow=0; iRow < 4*nlev; iRow++){
   vectB[iRow] = -vectB[iRow];
}//e-for
return 0;
}//e-makeVectorMinusB 




//=====================
// FUNCTION: freeMemory
//=====================
// Function to free all allocated memory
void freeMemory (int nlev, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                 double *ar_a13_c, double *ar_a13_f, double *ar_a23_c, double *ar_a23_f, 
                 double *ar_a33_c, double *ar_a33_f, 
                 double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, 
                 double *S_c, double *S_f, 
                 double *bb, double *xx, double **AA)
{
int i;

// Check if memory is allocated and free it
if (*ar_a11_c != 0) free(ar_a11_c);
if (*ar_a11_f != 0) free(ar_a11_f);
if (*ar_a12_c != 0) free(ar_a12_c);
if (*ar_a12_f != 0) free(ar_a12_f);
if (*ar_a13_c != 0) free(ar_a13_c);
if (*ar_a13_f != 0) free(ar_a13_f);
if (*ar_a23_c != 0) free(ar_a23_c);
if (*ar_a23_f != 0) free(ar_a23_f);
if (*ar_a33_c != 0) free(ar_a33_c);
if (*ar_a33_f != 0) free(ar_a33_f);
if (*ar_p1 != 0) free(ar_p1);
if (*ar_p2 != 0) free(ar_p2);
if (*ar_p3 != 0) free(ar_p3);
if (*ar_p4 != 0) free(ar_p4);
if (*S_c != 0) free(S_c);
if (*S_f != 0) free(S_f);
if (*bb != 0) free(bb);
if (*xx != 0) free(xx);

if (**AA != 0) {
  for (i=0; i < 4*nlev; i++) free(AA[i]);    
}//e-if

}//e-freeMemory




//=======================
// FUNCTION displayVector
//=======================
// Function to display vector B and vector X (to check their form);
void displayVector (int nlev, double *vect, char *name)
{
int iRow;
printf("\n");
printf(" -------------------------------\n");
printf("\n");
for (iRow=0; iRow < 4*nlev; iRow++){
   printf("Vector %s [%02d] = %e \n", name, iRow, vect[iRow]);
}//e-for
printf("\n");
printf(" -------------------------------\n");
printf("\n");
}//e-displayVector




//========================
// FUNCTION: displayMatrix
//========================
// Function to display matrix A or -A
void displayMatrix (int nlev, double **matrixA, char *name)
{
int iRow;
int jCol;
double value;

// Display matrix to check

printf("\n");

// Top header
printf ("Displaying matrix %s - first and then second half:\n", name);
printf("\n");
printf ("Displaying matrix %s - first half:\n", name);
printf("\n");
printf("    ");
for (jCol=0; jCol < 4*nlev; jCol++){
   if (jCol > (2*nlev-1)) continue;
   printf (" [%02d] ",  jCol);	// Header to printf the column number
}//e-for
printf ("\n");

for (iRow = 0; iRow < 4*nlev; iRow++){
   printf ("%02d:", iRow);
   for (jCol = 0; jCol < 4*nlev; jCol++){
      if (jCol > (2*nlev-1)) continue;

      value = matrixA[iRow][jCol];
      if (value > 0.0)  printf (" A%02d%02d", iRow, jCol); // Printf the indices of non-zero elements, later exchange with true values!
      else {
         if (value < 0.0)  printf ("  -1  ");  
         else              printf ("   0  "); 	 
      }//e-if

   }//e-for
   printf ("\n");
}//e-for

// Bottom footer
printf ("\n    ");    // Set spaces to match header and footer line with columns
for (jCol=0;  jCol < 4*nlev; jCol++){
   if (jCol > (2*nlev-1)) continue;
    printf (" [%02d] ", jCol);  // Footer to printf the column number
}//e-for

// Top header
printf ("\n");
printf ("\n");
printf ("Displaying matrix %s - second half:\n", name);
printf ("\n");
printf ("    ");
for (jCol=0; jCol < 4*nlev; jCol++){
   if (jCol < 2*nlev) continue;
   printf (" [%02d] ", jCol);
}//e-for
printf ("\n");

for (iRow = 0; iRow < 4*nlev; iRow++){
   printf ("%02d:", iRow);
   for (jCol = 0; jCol < 4*nlev; jCol++){
      if (jCol < 2*nlev) continue;

      value = matrixA[iRow][jCol];
      
      if (value > 0.0)  printf (" A%02d%02d", iRow, jCol); 
      else {
         if (value < 0.0)  printf ("  -1  ");
         else              printf ("   0  ");
      }//e-if
    
   }//e-for
   printf ("\n");
}//e-for

// Bottom footer
printf ("    ");
for (jCol=0;  jCol < 4*nlev; jCol++){
   if (jCol < 2*nlev) continue;
   printf (" [%02d] ", jCol);
}//e-for
printf ("\n\n");

}//e-displayMatrix





/***********************************************************************************/
/* Function: solve_gauss                                                  @31_30i@ */
/* Description:                                                                    */
/*  Solve a system of n linear equations, A*x = b, using the Gauss algorithm       */
/*  The pivot element is determined using 'relative column maximum strategy'.      */
/*  For a description of the algorithm see H.R.Schwarz: "Numerische Mathematik",   */
/*  pg. 21. Memory for the result vector res is allocated automatically.           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **A:            Matrix[n x n]  (see above).                             */
/*  double *b:             Vector[n] (see above).                                  */
/*  double n:              Number of equations.                                    */
/*  double **res:          Pointer to the result vector[n]; if no unique solution  */ 
/*                         exists, *res will be set to NULL.                       */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if no unique solution.                                          */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i31_30@ */
/***********************************************************************************/

