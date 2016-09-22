#include <stdio.h>

int nina_maxrand(int nlyr,
    double *dtauf, double *w0f,double *gf,
    double *dtauc, double *w0c,double *gc,
    double *cfrac,
    double mu0, double incSolar, double albedo,
    double *Sf, double *Ednf, double*Eupf,
    double *Sc, double *Ednc, double*Eupc) {

  int k;
  printf("Computing Nina's Maxrand with %d layers \n", nlyr);

  for(k=0; k<nlyr+1; k++) {
    Sf[k]   = 100 + k;
    Sc[k]   = 200 + k;
    Ednf[k] = 300 + k;
    Ednc[k] = 400 + k;
    Eupf[k] = 500 + k;
    Eupc[k] = 600 + k;
  }

  return 0;
}
