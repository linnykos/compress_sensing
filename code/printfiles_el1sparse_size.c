#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "myalloc.h"

#define LEN_GAUSS 1000000
#define LEN_IDX 1500
#define GAUSS_OFFSET 100
#define IDX_OFFSET 150

main(int argc, char **argv)
{

	//printf("hi");

    int m;		/* number of constraints */
    int n;		/* number of variables */
    double **A=NULL, **B=NULL, **X0=NULL;
    int knum = 100;
	int offset = atoi(argv[2]);
	int counter = offset * GAUSS_OFFSET ;
	int counter2 = offset * IDX_OFFSET ;

    int i, j, k, i1, i2, j1, j2, n1, n2, d1, d2, n12, d12, b_idx, n_idx, jstar;

	FILE* fp;

    n1 = 33;
    n2 = 34;
    d1 = 141;
    d2 = atoi(argv[1]);

	//printf("%d",d1);
	//printf("%d",d2);

    MALLOC(        X0, d1,  double *);
    for (i=0; i<d1; i++) {
        CALLOC( X0[i], d2,  double );      
    }

    MALLOC(        A, n1,  double *);
    for (i=0; i<n1; i++) {
        MALLOC( A[i], d1,  double);
    }

    MALLOC(        B, n2,  double *);
    for (i=0; i<n2; i++) {
        MALLOC( B[i], d2,  double);
    }


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
	int rng_i[LEN_IDX];
   int rng_j[LEN_IDX];

   i = 0;
   file = fopen("rng_idx.csv", "r");
   while(fscanf(file, "%d,", &num2) > 0){
      rng_i[i] = (num2-1)/d2;
      rng_j[i] = (num2-1)%d2;

      i++;
      if (i>=LEN_IDX) break;
   }

   fclose(file);


    int kk;  
    for (kk=0; kk<knum; kk++) {


//	printf("%d\n",counter2);
//	printf("%d\n", rng_i[counter2]);
//	printf("%d\n", rng_j[counter2]);

       j1 = rng_i[counter2];
       j2 = rng_j[counter2];
	counter2++;
	if(j1<d1 & j2<d2){
       if (X0[j1][j2]==0){
           X0[j1][j2] = 1.;
       } else{
          kk--;
          }
	if (counter2>=LEN_IDX) counter2 = 0;
    
	}
    }

    for (j1=0; j1<d1; j1++) {
       for (i1=0; i1<n1; i1++) {
          A[i1][j1] = rng_gauss[counter++];
		if (counter>=LEN_GAUSS) counter = 0;
       }
    }
  
    for (j2=0; j2<d2; j2++) {
       for (i2=0; i2<n2; i2++) {
          B[i2][j2] = rng_gauss[counter++];
		if (counter>=LEN_GAUSS) counter = 0;
       }
    }
	

	if ((fp = fopen("d.dat","w")) == NULL){
		printf("cannot open %s \n", "d.dat");
		return 1;
	}
	fprintf(fp,"%d",d2);
	fclose(fp);

//printf("%5.2f",d2);

    if ( ( fp = fopen("X0.dat","w") ) == NULL ) {
	printf("cannot open %s \n", "X0");
	return 1;
    }
    /*for (j2=0; j2<d2; j2++) {
      fprintf(fp,"%3d ", j2+1);
    }
    fprintf(fp,":= \n");*/
    for (j1=0; j1<d1; j1++) {
      //fprintf(fp,"%3d ", j1+1);
      for (j2=0; j2<d2; j2++) {
	  fprintf(fp,"%5.2f ", X0[j1][j2]);
      }
      //fprintf(fp,"\n");
    }
    fclose(fp);

    if ( ( fp = fopen("A.dat","w") ) == NULL ) {
	printf("cannot open %s \n", "A");
	return 1;
    }
   /* for (j2=0; j2<d1; j2++) {
      fprintf(fp,"%3d ", j2+1);
    }
    fprintf(fp,":= \n");*/
   for (j1=0; j1<n1; j1++) {
      //fprintf(fp,"%3d ", j1+1);
      for (j2=0; j2<d1; j2++) {
	  fprintf(fp,"%10.7f ", A[j1][j2]);
      }
      //fprintf(fp,"\n");
    }
    fclose(fp);

    if ( ( fp = fopen("B.dat","w") ) == NULL ) {
	printf("cannot open %s \n", "B");
	return 1;
    }
    /*for (j2=0; j2<d2; j2++) {
      fprintf(fp,"%3d ", j2+1);
    }
    fprintf(fp,":= \n");*/
    for (j1=0; j1<n2; j1++) {
      //fprintf(fp,"%3d ", j1+1);
      for (j2=0; j2<d2; j2++) {
	  fprintf(fp,"%10.7f ", B[j1][j2]);
      }
      //fprintf(fp,"\n");
    }
    fclose(fp);

	for (i=0; i<n1; i++){
		FREE(A[i]);
	}

	for (i=0; i<n2; i++){
		FREE(B[i]);
	}
	FREE(A);
	FREE(B);
	FREE(X0);
   
 return 0;
}

