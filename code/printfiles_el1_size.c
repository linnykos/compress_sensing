#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "myalloc.h"

#define LEN_GAUSS 100000
#define LEN_IDX 1500
#define GAUSS_OFFSET 100
#define IDX_OFFSET 150

main(int argc, char **argv)
{
    int m;		/* number of constraints */
    int n;		/* number of variables */
    double **U=NULL, *x0=NULL, **A=NULL, **B=NULL;
    int knum = 100;
	int offset = atoi(argv[2]);
	int counter = offset * GAUSS_OFFSET ;
	int counter2 = offset * IDX_OFFSET ;

    int i, j, k, i1, i2, j1, j2, n1, n2, d1, d2, n12, d12, b_idx, n_idx, jstar, m0, n0;

	FILE* fp;

    n1 = 33;
    n2 = 34;
    d1 = 141;
    d2 = atoi(argv[1]);

    m0 = n1*n2;
    n0 = d1*d2;

    MALLOC(        x0, n0,  double );

    CALLOC(        U, m0,  double *);
    for (i=0; i<m0; i++) {
       CALLOC(        U[i], n0,  double);
    }
    CALLOC( A, n1, double *);
    for (i=0; i<n1; i++){
       CALLOC( A[i], d1, double);
    }
    CALLOC( B, n2, double *);
    for (i=0; i<n2; i++){
       CALLOC(B[i], d2, double);
    }

    int num2 = 0;
    int rng_idx[LEN_IDX];
    i = 0;
    FILE* file = fopen("rng_idx.csv", "r");

    while(fscanf(file, "%d,", &num2) > 0){
       rng_idx[i] = num2-1;
       i++;
      if (i>=LEN_IDX) break;
   }

   fclose(file);

   double num = 0;
   double rng_gauss[LEN_GAUSS];
   i = 0;
   file = fopen("rng_gauss.csv", "r");
   while(fscanf(file, "%lf,", &num) > 0)
   {
      rng_gauss[i++] = num;
      if(i>=LEN_GAUSS) break;
   }
   fclose(file);
	
    int kk;  
    for (kk=0; kk<knum; kk++) {
       j = rng_idx[counter2];
       counter2++;
       if (x0[j]==0){
          x0[j] = 1.;
       } else{
          kk--;
       }
	if (counter2>=LEN_IDX) counter2 = 0;
    }

    for (j1=0; j1<d1; j1++){
       for (i1=0; i1<n1; i1++){
          A[i1][j1] = rng_gauss[counter++];
		if (counter>=LEN_GAUSS) counter = 0;
       }
    }
    for (j2=0; j2<d2; j2++){
       for (i2=0; i2<n2; i2++){
          B[i2][j2] = rng_gauss[counter++];
		if (counter>=LEN_GAUSS) counter = 0;
       }
    }

    int tmp1, tmp2, tmp3, tmp4;
    for (i=0; i<m0; i++){
       for (j=0; j<n0; j++){
          tmp1 = i%n2;
          tmp2 = j%d2;
          tmp3 = i/n2;
          tmp4 = j/d2;
          U[i][j] = B[tmp1][tmp2]*A[tmp3][tmp4];
          if(counter >= LEN_GAUSS) counter = 0;
       }
    }
    
    for (i=0; i<n1; i++){
       FREE(A[i]);
    }
    for (i=0; i<n2; i++){
       FREE(B[i]);
    }
    FREE(A);
    FREE(B);

	

        if ( ( fp = fopen("x0.dat","w") ) == NULL ) {
           printf("cannot open %s \n", "x0");
           return 1;
        }
        for (j2=0; j2<n0; j2++) {
           fprintf(fp,"%5.2f ", x0[j2]);
        }

        fclose(fp);

        if ( ( fp = fopen("U.dat","w") ) == NULL ) {
           printf("cannot open %s \n", "A");
           return 1;
        }
        //for (j2=0; j2<n0; j2++) {
        //   fprintf(fp,"%3d ", j2+1);
        //}
        //fprintf(fp,":= \n");
        for (j1=0; j1<m0; j1++) {
           //fprintf(fp,"%3d ", j1+1);
           for (j2=0; j2<n0; j2++) {
              fprintf(fp,"%10.7f ", U[j1][j2]);
           }
           //fprintf(fp,"\n");
        }
        fclose(fp);

	if ( ( fp = fopen("d.dat","w") ) == NULL) {
		printf("cannot open %s \n", "d.dat");
		return 1;
	}
	fprintf(fp,"%3d",n0);
	fclose(fp);

	for (i=0; i<m0; i++){
		FREE(U[i]);
	}
	FREE(U);
	FREE(x0);

    return 0;
}

