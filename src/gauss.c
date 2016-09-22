/*--------------------------------------------------------------------
 * $Id: equation.c 2623 2011-12-23 10:52:38Z robert.buras $
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

//#include "numeric.h"

#include "gauss.h"

// Copied from equation.h by Nina:
//#define GAUSS_SINGULAR  -20

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

int solve_gauss (double **A, double *b, int n, double **res)
{
  int i=0, j=0, k=0, p=0;
  double div=0.0;
  double *tmp1=NULL;
  double tmp2=0.0;
  double sum=0.0;
  double max = (0.0 - DBL_MAX);
  

  for (k=0; k<n; k++)  {
    /* get pivot element (line p) */
    max = (0.0 - DBL_MAX);
    p = 0;
    
    for (i=k; i<n; i++)  {
      
      sum = 0.0;
      for (j=k; j<n; j++)
	sum += fabs(A[i][j]);
      if (sum == 0.0)
	return GAUSS_SINGULAR;
      
      if ( (tmp2 = fabs (A[i][k]) / sum) > max )  {
	max = tmp2;
	p = i;
      }
    }
    
    /* exchange lines k and p */
    tmp1 = A[k];
    A[k] = A[p];
    A[p] = tmp1;

    tmp2 = b[k];
    b[k] = b[p];
    b[p] = tmp2;
    
    
    if ( (div = A[k][k]) == 0)   /* no definite solution */
      return GAUSS_SINGULAR;
    
    for (i=k; i<n; i++)
      A[k][i] /= div;
    
    b[k] /= div;
    
    for (i=k+1; i<n; i++)  {
      div = A[i][k];
      for (j=k; j<n; j++)
	A[i][j] = A[i][j] - A[k][j] * div;
      b[i] = b[i] - b[k] * div;
    }
  }
  

  /* allocate memory for result vector */
  *res = (double *) calloc (n, sizeof(double));

  for (i=n-1; i>=0; i--)  {
    (*res)[i] += b[i];
    for (k=i+1; k<n; k++)
      (*res)[i] -= A[i][k] * (*res)[k];
  }

  return 0;
}



