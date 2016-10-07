/*--------------------------------------------------------------------
 * $Id: twomaxrnd.c 2623 2011-12-23 10:52:38Z nina.crnivec $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

/* full solution of the multi-layer twostream equation solar + thermal */
/* with maximum-random overlap assumption for partial cloudiness */
/* written by Nina Crnivec */

#ifndef __twomaxrnd_h
#define __twomaxrnd_h

#if defined (__cplusplus)
extern "C" {
#endif

int calcp1p2p3p4 (int nlev, double *cf, double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4);

void delta_scale_hg (double tau, double ssa, double g, 
                    double *tauscale, double *ssascale, double *gscale);

void eddington_coeffc (double dtau, double g, double omega0, double mu0, double *a11, double *a12,
                       double *a13, double *a23, double *a33);

void calcThermalComponents(int ilyr, double *B, double dtau, double omega0, double g,
			   double *theComp1, double *theComp2);

int buildMatrixA (int nlev, double Ag, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                  double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, double **matrixA);

int makeMatrixA2 (int nlev, double **matrixA);

int buildVectorBsol (int nlev, double Ag, double mu0, double *ar_a13_c, double *ar_a13_f, 
		     double *ar_a23_c, double *ar_a23_f, double *ar_S_c, double *ar_S_f, 
		     double *ar_p1, double *ar_p3, double *vectB);

int buildVectorBthe (int nlev, double *ar_theComp1_c, double *ar_theComp1_f, 
		     double *ar_theComp2_c, double *ar_theComp2_f, 
		     double *ar_p1, double *ar_p3, double *vectB);

int makeVectorMinusB (int nlev, double *vectB);

void displayMatrix (int nlev, double **matrixA, char *name);

void displayVector (int nlev, double *vect, char *name);

void freeMemory (int nlev, double *ar_a11_c, double *ar_a11_f, double *ar_a12_c, double *ar_a12_f, 
                 double *ar_a13_c, double *ar_a13_f, double *ar_a23_c, double *ar_a23_f, 
                 double *ar_a33_c, double *ar_a33_f, 
                 double *ar_p1, double *ar_p2, double *ar_p3, double *ar_p4, 
                 double *S_c, double *S_f, double *Edir_c, double *Edir_f, 
		 double *Eup_c, double *Eup_f, double *Edn_c, double *Edn_f,
                 double *bb, double *xx, double **AA);

static int twostream_maxrand (double *dtau_c, double *omega0_c, double *g_c, // cloudy region parameters
			      double *dtau_f, double *omega0_f, double *g_f, // free region parameters
			      double *cf, int nlev, 
			      double S0, double mu0, double Ag, 
			      double Bg, double *B, int delta,
			      double **Edir, double **Edn, double **Eup, double **Lavg);
#if defined (__cplusplus)
extern "C" {
#endif

#endif /* _twomaxrnd_h */
