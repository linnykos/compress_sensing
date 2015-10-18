/*****************************************************************************

                Implementation of the 
		Primal Dual (i.e. Self Dual) Simplex Method
		R. Vanderbei, 3 October 1994

Solves problems in the form:

	     T
	max c x

	A x  = b
	  x >= 0

A is an m by N matrix (it is convenient to reserve n for 
the difference N-m).  Artificial variables have been 
added (hence N > m).  One may assume that the last 
m columns of A are invertible and hence can be used as 
a starting basis.

******************************************************************************/         
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef QuadPrec
#include "Quad.h"
#define double Quad
#else
#define high(x) (x)
#endif

#include "linalg.h"
#include "lu.h"
#include "myalloc.h"
#include "macros.h"

#define EPS0 1.0e-3
#define EPS1 1.0e-8
#define EPS2 1.0e-12
#define EPS3 1.0e-4
#define MAX_ITER 1000000
#define LEN_GAUSS 100000
#define LEN_IDX 1500
#define GAUSS_OFFSET 10000
#define IDX_OFFSET 150

double sdotprod(double *c, double *x_B, int *basics, int m);

void Nt_times_y( 
    int N, 
    double *at, 
    int *iat, 
    int *kat, 
    int *basicflag,
    double *y, 
    int *iy, 
    int ny, 
    double *yN,
    int *iyN,
    int *pnyN
);

int ratio_test(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y, 
	double *ybar, 
	double mu
);

int solver2(
    int m,		/* number of constraints */
    int n,		/* number of variables */
    int nz,		/* number of nonzeros in sparse constraint matrix */
    int *ia, 		/* array row indices */
    int *ka, 		/* array of indices into ia and a */
    double *a,		/* array of nonzeros in the constraint matrix */
    double *b, 		/* right-hand side */
    double *c,          /* objective coefficients */
    double *c2,         /* objective coefficients */
    double  f, 		/* objective function shift */
    int *basics,
    int *nonbasics,
    int *basicflag,
    int n0,
    double *betahat
    );

main(int argc, char **argv)
{
    int m;		/* number of constraints */
    int n;		/* number of variables */
    int nz;		/* number of nonzeros in sparse constraint matrix */
    int *ia; 		/* array row indices */
    int *ka; 		/* array of indices into ia and a */
    double *a;		/* array of nonzeros in the constraint matrix */
    double *b; 		/* right-hand side */
    double *c;          /* objective coefficients */
    double *c2;         /* objective coefficients */
    int     *basics;	/* list of basic variable indices */
    int     *nonbasics;	/* list of non-basic variable indices */
    int     *basicflag; /* neg if nonbasic, pos if basic */
    double **U, **X, **Z;
    time_t start_time, end_time;
    double elapsed_time;
    int knum = atoi(argv[1]);
	int offset = atoi(argv[2]);
	int counter = offset * GAUSS_OFFSET ;
	int counter2 = offset * IDX_OFFSET ;
 

        int i, j, k, m0, n0, b_idx, n_idx, j1, j2, i1, i2, n1, n2, d1, d2;
    
	double *beta0, *betahat;


        n1 = 33;
        n2 = 34;
        d1 = 141;
        d2 = 142;


 FILE* file = fopen("rng_gauss.csv", "r");
   i = 0;
   j = 0;
   double num = 0;
   double rng_gauss[LEN_GAUSS];

   while(fscanf(file, "%lf,", &num) > 0)
   {
       rng_gauss[i++] = num;
      if (i>=LEN_GAUSS) break;
   }

   fclose(file);

   int num2 = 0;
   int rng_idx[LEN_IDX];
   i = 0;
   file = fopen("rng_idx.csv", "r");
   while(fscanf(file, "%d,", &num2) > 0){
      rng_idx[i] = num2-1;
      i++;
      if (i>=LEN_IDX) break;
   }

   fclose(file);


        m0 = n1*n2;
        n0 = d1*d2;

        m  = m0;
        n  = 2*n0 + 2*m0;
        nz = 2*m0*n0 + 2*m0;
        CALLOC(        a, nz,  double );      
        CALLOC(       ia, nz,   int );      
        CALLOC(       ka, n+1,  int );      
        CALLOC(        b, m,   double );      
        CALLOC(        c, n,   double );      
        CALLOC(       c2, n,   double );      
        CALLOC(   basics, m,   int );      
        CALLOC(nonbasics, n-m, int );      
        CALLOC(basicflag, n,   int );      
        CALLOC(    beta0, n0,  double );      
        CALLOC(  betahat, n0,  double );      
    
	CALLOC(        U, m0,  double *);
        for (i=0; i<m0; i++) {
           CALLOC(        U[i], n0,  double);
        }
	CALLOC(        X, n1,  double *);
        for (i=0; i<n1; i++) {
           CALLOC( X[i], d1,  double);
        }

        CALLOC(        Z, n2,  double *);
        for (i=0; i<n2; i++) {
           CALLOC( Z[i], d2,  double);
           }

        int kk;
        for (kk=0; kk<knum; kk++) {
           j = rng_idx[counter2];
           counter2++;
           if (beta0[j]==0){
              beta0[j] = 1.;
           } else{
              kk--;
           }
		if(counter2>=LEN_IDX) counter2 = 0;
        }

        for (j1=0; j1<d1; j1++) {
           for (i1=0; i1<n1; i1++) {
              X[i1][j1] = rng_gauss[counter++];
			if(counter>=LEN_GAUSS) counter = 0;
           }
        }

        for (j2=0; j2<d2; j2++) {
           for (i2=0; i2<n2; i2++) {
              Z[i2][j2] = rng_gauss[counter++];
			if(counter>=LEN_GAUSS) counter = 0;
           }
        }

        int tmp1;
        int tmp2;
        int tmp3;
        int tmp4;
        for (i=0; i<m0; i++){
           for(j=0; j<n0; j++){
              tmp1 = i%n2; /*for Z (B)*/
              tmp2 = j%d2; /*for Z (B)*/
              tmp3 = i/n2; /*for X (A)*/
              tmp4 = j/d2; /*for X (A)*/
              U[i][j] = Z[tmp1][tmp2]*X[tmp3][tmp4];
              if(counter >= LEN_GAUSS) counter = 0;
           }
        }

    for (i=0; i<n1; i++){
       FREE(X[i]);
    }
    for (i=0; i<n2; i++){
       FREE(Z[i]);
       }

     FREE(X);
        FREE(Z);

        for (i=0; i<m0; i++) {
           b[i] = 0;
           for (j=0; j<n0; j++) {
              b[i] += U[i][j]*beta0[j];
           }
        }

        time(&start_time);

    /*****************************************************************
     * Structure of the problem:
     *
     *                  +          -         +      -
     *              beta       beta       eps    eps
     *
     * zeta    =                         -e^T   -e^T
     * zetabar =  -mu e^T   -mu e^T
     * --------------------------------------------
     *               X        -X          I     -I   = y
     *
     ***************************************************************/

    k=0;
    for (j=0; j<n0; j++) {
	ka[j] = k;
	for (i=0; i<m0; i++) {
	    a[k] = U[i][j];
	    ia[k] = i;
	    k++;
	}
	c [j] =  0;
	c2[j] = -1;
    }
    for (j=n0; j<2*n0; j++) {
	ka[j] = k;
	for (i=0; i<m0; i++) {
	    a[k] = -U[i][j-n0];
	    ia[k] = i;
	    k++;
	}
	c [j] =  0;
	c2[j] = -1;
    }
    i=0;
    for (j=2*n0; j<2*n0+m0; j++) {
	ka[j] = k;
        a[k] = 1;
        ia[k] = i;
        k++;
	i++;
	c [j] = -1;
	c2[j] =  0;
    }
    i=0;
    for (j=2*n0+m0; j<2*n0+2*m0; j++) {
	ka[j] = k;
        a[k] = -1;
        ia[k] = i;
        k++;
	i++;
	c [j] = -1;
	c2[j] =  0;
    }
    ka[n] = k;



    b_idx = 0; n_idx = 0;
    for (j=0; j<2*n0; j++) {
	nonbasics[n_idx] = j; basicflag[j] = -n_idx-1; n_idx++;
    }
    for (i=0; i<m0; i++) {
	j = 2*n0+i;
	if (b[i] >= 0) {
	    basics[b_idx] = j; basicflag[j] = b_idx; b_idx++;
	} else {
	    nonbasics[n_idx] = j; basicflag[j] = -n_idx-1; n_idx++;
	}
    }
    for (i=0; i<m0; i++) {
	j = 2*n0+m0+i;
	if (b[i] < 0) {
	    basics[b_idx] = j; basicflag[j] = b_idx; b_idx++;
	} else {
	    nonbasics[n_idx] = j; basicflag[j] = -n_idx-1; n_idx++;
	}
    }

   /* printf("starting now \n"); fflush(stdout);*/

    solver2(m,n,nz,ia,ka,a,b,c,c2,0,basics,nonbasics,basicflag,n0,betahat);
    time(&end_time);
    elapsed_time = difftime(end_time, start_time);

    /*printf("=================================================\n");*/
	
    double maxerror = 0.0;
    double l1error = 0.0;
    double tmp = 0.;
    int happy = 1;
    double maxbeta0 = 0.0;
    double betahatsum = 0.0;
    for (j=0; j<n0; j++) {
       if (ABS(beta0[j]) > maxbeta0) { maxbeta0 = ABS(beta0[j]); }
    }

    for (j=0; j<n0; j++) {
       betahatsum = betahatsum + betahat[j];
       tmp = ABS(beta0[j]-betahat[j]);
       l1error = l1error+tmp;
       if (tmp>maxerror) maxerror = tmp;
       if (ABS(beta0[j]-betahat[j]) > 1e-5*maxbeta0) {
          happy = 0;
          break;
       }
    }

    /*printf("Betahat sum %f\n", betahatsum);
    printf("Max error: %f\n", maxerror);
    printf("L1 error: %f\n", l1error);

    if (happy) {
       printf("happy: ");
    } else {
       printf("unhappy: ");
    }
    
    printf("Solution time (seconds): %0.2lf \n", elapsed_time);*/


	printf("%d, %f, %f, %f\n",knum, elapsed_time, maxerror, l1error);

    lu_clo();

    for(i=0; i<m0; i++){
       FREE(U[i]);
    }

    FREE(U);
    FREE(beta0);
    FREE(betahat);
    
    FREE(a);
    FREE(ia); 
    FREE(ka);
    FREE(b);
    FREE(c);
    FREE(c2);
    FREE(basicflag);


    return 0;
}

int solver2(
    int m,		/* number of constraints */
    int N,		/* number of variables */
    int nz,		/* number of nonzeros in sparse constraint matrix */
    int *ia, 		/* array row indices */
    int *ka, 		/* array of indices into ia and a */
    double *a,		/* array of nonzeros in the constraint matrix */
    double *b, 		/* right-hand side */
    double *c,          /* objective coefficients */
    double *c2,         /* objective coefficients */
    double  f, 		/* objective function shift */
    int *basics,
    int *nonbasics,
    int *basicflag,
    int n0,
    double *betahat
    )
{
    FILE *fp;

    double  *x_B;	/* primal basics */
    double  *y_N;	/* dual nonbasics */

    double  *xbar_B;	/* primal basic perturbation */
    double  *ybar_N;	/* dual nonbasic perturbation */

    double  *dy_N;	/*  dual  basics step direction - values (sparse) */
    int    *idy_N;	/*  dual  basics step direction - row indices */
    int     ndy_N;	/* number of nonz in dy_N */

    double  *dx_B;	/* primal basics step direction - values (sparse) */
    int    *idx_B;	/* primal basics step direction - row indices */
    int     ndx_B;	/* number of nonz in dx_B */

    double  *at;	/* sparse data structure for A^T */
    int    *iat;
    int    *kat;

    int     col_in;	/* entering column; index in 'nonbasics' */
    int     col_out;	/* leaving column; index in 'basics' */

    int     iter = 0;	/* number of iterations */
    int     i,j,k,n,v=0;

    double  s, t, sbar, tbar, mu=HUGE_VAL, old_mu, primal_obj;

    double  *vec;
    int    *ivec;
    int     nvec;

    int     new_pt = TRUE;

    int	    status;

    int     from_scratch;

    /*******************************************************************
    * For convenience, we put...
    *******************************************************************/

    n = N-m;

    /*******************************************************************
    * Read in the Data and initialize the common memory sites.
    *******************************************************************/

    CALLOC(    x_B, m,   double );      
    CALLOC( xbar_B, m,   double );      
    CALLOC(   dx_B, m,   double );      
    CALLOC(    y_N, n,   double );      
    CALLOC( ybar_N, n,   double );      
    CALLOC(   dy_N, n,   double );      
    CALLOC(    vec, N,   double );
    CALLOC(  idx_B, m,    int );      
    CALLOC(  idy_N, n,    int );      
    CALLOC(   ivec, N,    int );
    CALLOC(     at, nz,  double );
    CALLOC(    iat, nz,   int );
    CALLOC(    kat, m+1,  int );

    /**************************************************************** 
    *  Initialization.              				    *
    ****************************************************************/

    atnum(m,N,ka,ia,a,kat,iat,at);

   
    lufac( m, ka, ia, a, basics, 0 );

    for (j=0; j<n; j++) {
	  y_N[j] = 0;
    }
    nvec = 0;
    for (i=0; i<m; i++) {
	if (c[basics[i]] != 0.0) {
	    vec[nvec] = c[basics[i]];
	    ivec[nvec] = i;
	    nvec++;
	}
    }
    btsolve( m, vec, ivec, &nvec );  		
    Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );
    for (k=0; k<ndy_N; k++) {
	y_N[idy_N[k]] = dy_N[k];
    }
    for (j=0; j<n; j++) {
	y_N[j] -= c[nonbasics[j]];
    }

    for (j=0; j<n; j++) {
	   ybar_N[j] = 0;
    }
    nvec = 0;
    for (i=0; i<m; i++) {
	if (c2[basics[i]] != 0.0) {
	    vec[nvec] = c2[basics[i]];
	    ivec[nvec] = i;
	    nvec++;
	}
    }
    btsolve( m, vec, ivec, &nvec );  		
    Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );
    for (k=0; k<ndy_N; k++) {
	ybar_N[idy_N[k]] = dy_N[k];
    }
    for (j=0; j<n; j++) {
	ybar_N[j] -= c2[nonbasics[j]];
	if (ybar_N[j] < 0) printf("error: ybar_N[%d] = %e \n", j, ybar_N[j]);
    }

    for (i=0; i<m; i++) {
	       x_B[i] = 0;
	    xbar_B[i] = 0; 
    }
    nvec = 0;
    for (i=0; i<m; i++) {
      if ( b[i] != 0.0 ) {
	 vec[nvec] = b[i];
	ivec[nvec] = i;
	nvec++;
      }
    }
    bsolve( m, vec, ivec, &nvec );
    for (i=0; i<nvec; i++) {
	x_B[ivec[i]] = vec[i];
	if (vec[i] < 0) printf("error: x_B[%d] = %e \n", i, vec[i]);
    }
   /*
    printf ("m = %d,n = %d,nz = %d\n",m,N,nz);
    printf(
"---------------------------------------------------------------------------\n"
"          |   Primal      |        |                           arithmetic  \n"
"  Iter    |  Obj Value    |   mu   |   nonz(L)     nonz(U)     operations  \n"
"- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
    );
    fflush(stdout);
*/
    /****************************************************************
    * 	Main loop                                                   *
    ****************************************************************/

    for (iter=0; iter<MAX_ITER; iter++) {

      /*************************************************************
      * STEP 1: Find mu                                            *
      *************************************************************/

      old_mu = mu;
      mu = -HUGE_VAL;
      col_in  = -1;
      for (j=0; j<n; j++) {
		if (ybar_N[j] > EPS2) { 
			if ( mu < -y_N[j]/ybar_N[j] ) {
			     mu = -y_N[j]/ybar_N[j];
			     col_in  = j;
			}
		}
      }
      col_out = -1;
      for (i=0; i<m; i++) {
		if (xbar_B[i] > EPS2) { 
			if ( mu < -x_B[i]/xbar_B[i] ) {
			     mu = -x_B[i]/xbar_B[i];
			     col_out = i;
			     col_in  = -1;
			}
		}
      }

      /*************************************************************
      * STEP 0: Record current portfolio                           *
      *************************************************************/

      primal_obj = sdotprod(c,x_B,basics,m) + f;
      if ( mu <= EPS3 || primal_obj > -EPS0) {	/* OPTIMAL */
	  for (i=0; i<n0; i++) { betahat[i] = 0; }
	  for (i=0; i<m; i++) {
	    if (basics[i] < n0 && x_B[i] > EPS0) {
	      betahat[basics[i]] = x_B[i];
	    }
	    else if (basics[i] < 2*n0 && x_B[i] > EPS0) {
	      betahat[basics[i]-n0] = -x_B[i];
	    }
	  }
	  status = 0;
	  break;
      }

      if ( col_out >= 0 ) {

        /*************************************************************
	*                          -1  T                             *
	* STEP 2: Compute dy  = -(B  N) e                            * 
	*                   N            i			     *
	*         where i = col_out                                  *
        *************************************************************/

	vec[0] = -1.0;
	ivec[0] = col_out;
	nvec = 1;

	btsolve( m, vec, ivec, &nvec );  		

	Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

        /*************************************************************
	* STEP 3: Ratio test to find entering column                 * 
        *************************************************************/

	col_in = ratio_test( dy_N, idy_N, ndy_N, y_N, ybar_N, mu );

	if (col_in == -1) { 	/* INFEASIBLE */
	    status = 2;
	    break;
	}

        /*************************************************************
	*                        -1                                  *
	* STEP 4: Compute dx  = B  N e                               * 
	*                   B         j                              *
	*                                                            *
        *************************************************************/

	j = nonbasics[col_in];
	for (i=0, k=ka[j]; k<ka[j+1]; i++, k++) {
	     dx_B[i] =  a[k];
	    idx_B[i] = ia[k];
	}
	ndx_B = i;

	bsolve( m, dx_B, idx_B, &ndx_B );

      } else {

        /*************************************************************
	*                        -1                                  *
	* STEP 2: Compute dx  = B  N e                               * 
	*                   B         j                              *
        *************************************************************/

	j = nonbasics[col_in];
	for (i=0, k=ka[j]; k<ka[j+1]; i++, k++) {
	     dx_B[i] =  a[k];
	    idx_B[i] = ia[k];
	}
	ndx_B = i;

	bsolve( m, dx_B, idx_B, &ndx_B );

        /*************************************************************
	* STEP 3: Ratio test to find leaving column                  * 
        *************************************************************/

	col_out = ratio_test( dx_B, idx_B, ndx_B, x_B, xbar_B, mu );

	if (col_out == -1) {	/* UNBOUNDED */
	    status = 1;
	    break;
	}

        /*************************************************************
	*                          -1  T                             *
	* STEP 4: Compute dy  = -(B  N) e                            * 
	*                   N            i			     *
	*                                                            *
        *************************************************************/

	 vec[0] = -1.0;
	ivec[0] = col_out;
	nvec = 1;

	btsolve( m, vec, ivec, &nvec );  		

	Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

      }

      /*************************************************************
      *                                                            *
      * STEP 5: Put       t = x /dx                                *
      *                        i   i                               *
      *                   _   _                                    *
      *                   t = x /dx                                *
      *                        i   i                               *
      *                   s = y /dy                                *
      *                        j   j                               *
      *                   _   _                                    *
      *                   s = y /dy                                *
      *                        j   j                               *
      *************************************************************/

      for (k=0; k<ndx_B; k++) if (idx_B[k] == col_out) break;

      t    =    x_B[col_out]/dx_B[k];
      tbar = xbar_B[col_out]/dx_B[k];

      for (k=0; k<ndy_N; k++) if (idy_N[k] == col_in) break;

      s    =    y_N[col_in]/dy_N[k];
      sbar = ybar_N[col_in]/dy_N[k];

      /*************************************************************
      *                                _    _    _                 *
      * STEP 7: Set y  = y  - s dy     y  = y  - s dy              *
      *              N    N       N     N    N       N             *
      *                                _    _                      *
      *             y  = s             y  = s                      *
      *              i                  i                          *
      *             _    _    _                                    *
      *             x  = x  - t dx     x  = x  - t dx              *
      *              B    B       B     B    B       B             *
      *             _    _                                         *
      *             x  = t             x  = t                      *
      *              j                  j                          *
      *************************************************************/

      for (k=0; k<ndy_N; k++) {
		j = idy_N[k];
		y_N[j]    -= s   *dy_N[k];
		ybar_N[j] -= sbar*dy_N[k];

      }
      
      y_N[col_in]    = s;
      ybar_N[col_in] = sbar;

      for (k=0; k<ndx_B; k++) {
		i = idx_B[k];
		x_B[i]    -= t   *dx_B[k];
		xbar_B[i] -= tbar*dx_B[k];

      }


      x_B[col_out]     = t;
      xbar_B[col_out]  = tbar;

      /*************************************************************
      * STEP 8: Update basis                                       * 
      *************************************************************/

      i =    basics[col_out];
      j = nonbasics[col_in];
      basics[col_out]   = j;
      nonbasics[col_in] = i;
      basicflag[i] = -col_in-1;
      basicflag[j] = col_out;

      /*************************************************************
      * STEP 9: Refactor basis and print statistics                *
      *************************************************************/

      from_scratch = refactor( m, ka, ia, a, basics, col_out, v );

      if (from_scratch) {
          primal_obj = sdotprod(c,x_B,basics,m) + f;
/*          printf("%8d   %14.7e %9.2e \n", iter, high(primal_obj), high(mu) );
            fflush(stdout);
*/
      }
    } 

    primal_obj = sdotprod(c,x_B,basics,m) + f;
    /*   printf("\n");
    printf("%8d   %14.7e %9.2e \n", iter, high(primal_obj), high(mu) );
    */

    /****************************************************************
    * 	Transcribe solution to x vector and dual solution to y      *
    ****************************************************************/

    /****************************************************************
    * 	Split out slack variables and shift dual variables.
    ****************************************************************/

  
    /****************************************************************
    * 	Free work space                                             *
    ****************************************************************/

 Nt_times_y(-1, at, iat, kat, basicflag, vec, ivec, nvec, dy_N, idy_N, &ndy_N);

    FREE(at);
    FREE(iat);

    FREE(xbar_B);
    FREE(kat);
    FREE(ybar_N);

    lu_clo();
    btsolve(0, vec, ivec, &nvec);
    bsolve(0, vec, ivec, &nvec);


    FREE(  vec );
    FREE( ivec );
    FREE(  x_B );
    FREE(  y_N );
    FREE( dx_B );
    FREE(idx_B );
    FREE( dy_N );
    FREE(idy_N );
    FREE( nonbasics );
    FREE( basics );

    return status;

}   /* End of solver */

void Display_Solution(int m,int *basics,double *X)
{
	int i;
	
	printf("SOLUTION:\n\n");
	for (i=0;i<m;i++)
		printf("  X[%d] = %lf\n",basics[i], high(X[i]) );
}

void Nt_times_y( 
    int n, 
    double *at, 
    int *iat, 
    int *kat, 
    int *basicflag,
    double *y, 
    int *iy, 
    int ny, 
    double *yN,
    int *iyN,
    int *pnyN
)
{
    int i,j,jj,k,kk;

    static double *a=NULL;
    static int  *tag=NULL;
    static int *link=NULL;
    static int  currtag=1;

	if (n==-1){
       if (a != NULL) FREE(a);
       if (tag != NULL) FREE(tag);
       link--;
       if (link != NULL) FREE(link);
       return;
    }
	

    if (  a  == NULL) CALLOC(  a,  n,   double);
    if ( tag == NULL) CALLOC( tag, n,   int);
    if (link == NULL) {CALLOC(link, n+2, int); link++;}

    jj = -1;
    for (k=0; k<ny; k++) {
	i = iy[k];
	for (kk=kat[i]; kk<kat[i+1]; kk++) {
	    j = iat[kk];
	    if (basicflag[j] < 0) {
		if (tag[j] != currtag) {
		    a[j] = 0.0;
		    tag[j] = currtag;
		    link[jj] = j;
		    jj = j;
		}
		a[j] += y[k]*at[kk];
	    }
	}
    }
    link[jj] = n;
    currtag++;

    k = 0;
    for (jj=link[-1]; jj<n; jj=link[jj]) {
	if ( ABS(a[jj]) > EPS1 ) {
             yN[k] = a[jj];
            iyN[k] = -basicflag[jj]-1;
            k++;
	}
    }
    *pnyN = k;
}

int ratio_test(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y, 
	double *ybar, 
	double mu
)
{
	int j, jj = -1, k, kk;
	double min = HUGE_VAL;

	for (k=0; k<ndy; k++) {
	    if ( dy[k] > EPS1 ) {
	        j = idy[k];
		if ( (y[j] + mu*ybar[j])/dy[k] < min ) {
			min = (y[j] + mu*ybar[j])/dy[k];
			 jj = j;
			 kk = k;
		}
	    }
	}

	return jj;
}

double sdotprod(double *c, double *x_B, int *basics, int m)
{
	int i;
	double prod = 0.0;

	for (i=0; i<m; i++) { prod += c[basics[i]]*x_B[i]; }

	return prod;
}
