#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "macros.h"
#include "myalloc.h"
#include "tree.h"
#include "heap.h"
#include "linalg.h"
#include "lu.h"

#define REFRAC_INT 150
#define E_N 200
#define E_NZ 20000
#define LARGE 100000

#define EPS    1.0e-14
#define EPS1   1.0e-14
#define EPSSOL 1.0e-4   /* Zero tolerance for consistent eqns w/dep rows */
#define EPSNUM 1.0e-5
#define NOREORD 0
#define MD  1

static void Gauss_Eta( 
    int m, 
    double *dx_B, 
    int *idx_B, 
    int *pndx_B 
);

static void Gauss_Eta_T( 
    int m, 
    double *vec, 
    int *ivec, 
    int *pnvec 
);

static int     *E_d=NULL;  /* Eta file; new column location */
static double  *E=NULL;    /* Eta file - values (sparse matrix) */
static int    *iE=NULL;    /* Eta file - row indices */
static int    *kE=NULL;    /* Eta file - column start positions */
static int     e_iter = 0; /* number of iterations since last refactorization */
static int     enz=0;


static  int     rank;
static  int     *kL =NULL, *iL =NULL, 
		*kLt=NULL, *iLt=NULL,
		*kU =NULL, *iU =NULL,
		*kUt=NULL, *iUt=NULL;
static  int     *colperm=NULL, *icolperm=NULL, *rowperm=NULL, *irowperm=NULL;
static  double  *L=NULL, *Lt=NULL, *U=NULL, *Ut=NULL, *diagU=NULL;

static  double  cumtime = 0.0;
static  double  ocumtime= 0.0;

struct valind {    /* nonzero entry */
        double d;  /* value */
        int    i;  /* row index */
};
typedef struct valind VALIND;

/*-----------------------------------------------------------------+
| LU factorization.                                                |
| Input:                                                           |
|    m          number of rows (= number of columns)               |
|    kA, iA, A  three array sparse representation of m by n        |
|               matrix A                                           |
|    basis      array of m out of n of the column indices of A     |
|               specify a submatrix B of A                         |
| Output:                                                          |
|    static global variables (only visible to other routines in    |
|    this file:                                                    |
|                                                                  |
|    rank       rank of B                                          |
|    kL, iL, L, three array sparse representation of L             |
|    kUt,iUt,Ut three array sparse representation of U transpose   |
|               without its diagonal                               |
|    diagU      diagonal entries of U                              |
|    colperm, icolperm, rowperm, irowperm                          |
|               column and row permutations and their inverses    */

int refactor(
    int m,
    int *kA,
    int *iA,
    double *A,
    int *basics,
    int col_out,
    int v
)
{
	double starttime, endtime;
	double rffactor = 1.0;
	int    from_scratch;
	int    k;

        /*------------------------------------------------------+
        | Check if it is time to refactor from scratch         */

	if (e_iter > 0) {
	    from_scratch = TRUE;
	    for (k=kE[e_iter]; k<kE[e_iter+1]; k++) {
	        if (iE[k] == col_out) {
		    from_scratch = FALSE;
		    break;
	        }
	    }
	} else {
	    from_scratch = FALSE;
	}
           if (e_iter >= REFRAC_INT || from_scratch == TRUE) {
		ocumtime = 0.0;
		cumtime  = 0.0;
		lufac( m, kA, iA, A, basics, v );
		cumtime  *= rffactor;
		enz = 0;
		e_iter = 0;
		from_scratch = TRUE;
		return from_scratch;
	}

	ocumtime  = cumtime;
	starttime = (double) clock();

	E_d[e_iter] = col_out;
	e_iter++;

	endtime = (double) clock();
	cumtime += endtime - starttime;

	from_scratch = FALSE;
	return from_scratch;
}

void lufac( int m, int *kA, int *iA, double *A, int *basis, int v )
{
        int     kk, kkk, tag, rowdeg, coldeg, row, col, row2, col2;

        int     i, j, k, cnt, lnz, unz, lnzbnd, unzbnd, okey, deg,
                heapnum, cur, method=MD;

        int     *degB=NULL, *degBt=NULL, *hkey=NULL, 
		*heap=NULL, *iheap=NULL, *iwork=NULL, *iwork2=NULL;

        VALIND  tempB, **B=NULL, **Bt=NULL;

        double  narth;

        double starttime, endtime;

	starttime = (double) clock();

        /*---------------------------------------------------------+
        | allocate space for perm and iperm.                      */

        if (colperm == NULL)  { MALLOC( colperm,  m, int ); }
	else                 { REALLOC( colperm,  m, int ); }
	if (icolperm == NULL) { MALLOC( icolperm, m, int ); }
	else		     { REALLOC( icolperm, m, int ); }
	if (rowperm == NULL)  { MALLOC( rowperm,  m, int ); }
	else		     { REALLOC( rowperm,  m, int ); }
	if (irowperm == NULL) { MALLOC( irowperm, m, int ); }
	else		     { REALLOC( irowperm, m, int ); }

        /*---------------------------------------------------------+
        | allocate space for work arrays.                         */

        MALLOC( degB,    m, int   );
        MALLOC( degBt,   m, int   );
        MALLOC( hkey,    m, int   );
        MALLOC( heap,    m, int   );
        MALLOC( iheap,   m, int   );
        MALLOC( iwork,   m, int   );
        MALLOC( iwork2,  m, int   );

        heap--;         /* so that indexing starts from 1 */

        /*---------------------------------------------------------+
        | calculate degrees in B and Bt                           */

        for (i=0; i<m; i++) { degBt[i] = 0; }
        for (i=0; i<m; i++) {
                degB[i] = kA[ basis[i]+1 ] - kA[ basis[i] ];
                for (k=kA[ basis[i] ]; k<kA[ basis[i]+1 ]; k++) {
                        degBt[ iA[k] ]++;
                }
        }

        /*---------------------------------------------------------+
        | calculate initial estimate of number of nonzeros in      |
        | L and Ut                                                */

        lnzbnd = 0;
        for (i=0; i<m; i++) lnzbnd += degB[i];
        lnzbnd = lnzbnd/2;
        unzbnd   = lnzbnd;

        /*---------------------------------------------------------+
        | allocate enough space to store L and Ut                  |
        | (without any fillin)                                    */

	if (kL == NULL)   {  MALLOC(    kL,    m+1,  int    ); }
		     else { REALLOC(    kL,    m+1,  int    ); }
	if (iL == NULL)   {  MALLOC(    iL, lnzbnd,  int    ); }
		     else { REALLOC(    iL, lnzbnd,  int    ); }
	if (L == NULL)    {  MALLOC(     L, lnzbnd,  double ); }
	             else { REALLOC(     L, lnzbnd,  double ); }
	if (kUt == NULL)  {  MALLOC(   kUt,    m+1,  int    ); }
		     else { REALLOC(   kUt,    m+1,  int    ); }
	if (iUt == NULL)  {  MALLOC(   iUt, unzbnd,  int    ); }
		     else { REALLOC(   iUt, unzbnd,  int    ); }
	if (Ut == NULL)   {  MALLOC(    Ut, unzbnd,  double ); }
		     else { REALLOC(    Ut, unzbnd,  double ); }
	if (diagU == NULL){  MALLOC( diagU,      m,  double ); }
		     else { REALLOC( diagU,      m,  double ); }

        MALLOC( B,  m, VALIND * );
        MALLOC( Bt, m, VALIND * );
        for (i=0; i<m; i++) {
		B[i] = NULL;
		Bt[i] = NULL;
                MALLOC( B[i],  degB[i],  VALIND );
                MALLOC( Bt[i], degBt[i], VALIND );
        }

        /*---------------------------------------------------------+
        | initialize B and Bt                                     */

        for (i=0; i<m; i++) { iwork[i] = 0; }
        for (j=0; j<m; j++) {
            kkk = 0;
            for (k=kA[ basis[j] ]; k<kA[ basis[j]+1 ]; k++) {
                row = iA[k];
                kk  = iwork[row];
                B[j][kkk].i = row;
                B[j][kkk].d = A[k];
                Bt[row][kk].i = j;
                Bt[row][kk].d = A[k];
                iwork[row]++;
                kkk++;
            }
        }

        /*---------------------------------------------------------+
        | miscellaneous initializations.                          */

        for (i=0; i<m; i++) { 
            icolperm[i] = -1;
            irowperm[i] = -1;
            iwork[i] = 0; 
            iwork2[i] = -1; 
        }

        rank = m; tag = 0; lnz = 0; unz = 0; kL[0] = 0; kUt[0] = 0;

        /*---------------------------------------------------------+
        | hkey encodes the tie-breaking rule - currently the rule  |
        | is somewhat random.  to make it first occuring minimum,  |
        | change the formula to:                                   |
        |       hkey[node] = degree[node]*m + node;                |
        | warning: with this definition of hkey, there is the      |
        | possibility of integer overflow on moderately large      |
        | problems.                                                |
        |                                                         */

        for (j=0; j<m; j++) {
            if (method == MD) hkey[j] = degB[j];
            else              hkey[j] = j;

            if (hkey[j]==0) hkey[j]=m+1;
        }

        /*---------------------------------------------------------+
        | set up heap structure for quickly accessing minimum.    */

        heapnum = m;
        for (j=m-1; j>=0; j--) {
                cur = j+1;
                iheap[j] = cur;
                heap[cur] = j;
                hfall( heapnum, hkey, iheap, heap, cur );
        }

        /*---------------------------------------------------------+
        | the min degree ordering loop                            */

        for (i=0; i<m; i++) {

                /*------------------------------------------------+
                |  select column with min column degree          */

again:
                col    = heap[1];
                coldeg = degB[col];

                if (coldeg == 0) {
                    printf("singular matrix. rank deficiency = %d\n", m-i);
                    rank = i;
                    goto end;
                }

                /*------------------------------------------------+
                |  select pivot element from this column by       |
                |  choosing nonzero whose row is of minimal       |
                |  degree                                        */

                rowdeg = m+1;
                for (k=0; k<coldeg; k++) {
                    if ( degBt[ B[col][k].i ] < rowdeg 
                         && ABS( B[col][k].d ) > EPSNUM ) {
                        row    = B[col][k].i;
                        rowdeg = degBt[row];
                    }
                }
                if (rowdeg == m+1) {
                    hkey[col]=m+2;
                    hfall( heapnum, hkey, iheap, heap, iheap[col] ); 
                    if (hkey[heap[1]] == m+2) {
                        printf("singular matrix. rank deficiency = %d\n", m-i);
                        rank = i;
                        goto end;
                    } else {
                        goto again;
                    }
                }

                /*------------------------------------------------+
                |  update permutation information                */

                colperm[i] = col;
                icolperm[col] = i;

                rowperm[i] = row;
                irowperm[row] = i;

                /*------------------------------------------------+
                |  reallocate space for L, iL, Ut, and iUt as     |
                |        necessary.                               |
                |                                                 |
                |  lnz stores the number of nonzeros in L so far  |
                |  lnzbnd is an estimate of how many will be in L |
                |  unz stores the number of nonzeros in U so far  |
                |  unzbnd is an estimate of how many will be in U*/

                cnt = lnz + coldeg-1 + coldeg*rowdeg/2;
                if (cnt > lnzbnd) {
                    lnzbnd = cnt;
                    REALLOC(  L, lnzbnd, double );
                    REALLOC( iL, lnzbnd, int );
                }

                cnt = unz + rowdeg-1 + coldeg*rowdeg/2;
                if (cnt > unzbnd) {
                    unzbnd = cnt;
                    REALLOC(  Ut, unzbnd, double );
                    REALLOC( iUt, unzbnd, int );
                }

                /*------------------------------------------------+
                |  copy pivot column into L and pivot row into    |
                |  Ut.                                           */

                kL[i+1] = kL[i] + coldeg-1;

                for (k=0; k<coldeg; k++) {
                    if ( B[col][k].i != row ) {
                        iL[lnz] = B[col][k].i;
                         L[lnz] = B[col][k].d;
                           lnz++;
                    }
                }

                kUt[i+1] = kUt[i] + rowdeg-1;

                for (k=0; k<rowdeg; k++) {
                    if ( Bt[row][k].i != col ) {
                        iUt[unz] = Bt[row][k].i;
                         Ut[unz] = Bt[row][k].d;
                            unz++;
                    } else {
                        diagU[i] = Bt[row][k].d;
                    }
                }

                /*------------------------------------------------+
                |  remove eliminated elements from B and Bt      */

                for (k=0; k<coldeg; k++) {
                    row2 = B[col][k].i;
                    degBt[row2]--;
                    for (kk=0; Bt[row2][kk].i != col; kk++) ;

                    tempB = Bt[row2][ degBt[row2] ];
                    Bt[row2][ degBt[row2] ] = Bt[row2][kk];
                    Bt[row2][kk] = tempB;
                }

                for (k=0; k<rowdeg; k++) {
                    col2 = Bt[row][k].i;
                    degB[col2]--;
                    for (kk=0; B[col2][kk].i != row; kk++) ;

                    tempB = B[col2][ degB[col2] ];
                    B[col2][ degB[col2] ] = B[col2][kk];
                    B[col2][kk] = tempB;
                }
                degB[col] = 0;
                degBt[row] = 0;

                /*------------------------------------------------+
                |  update heap                                   */

                okey = hkey[col];
                heap[1] = heap[heapnum];
                iheap[heap[1]] = 1;
                heapnum--;
                if (okey < hkey[heap[1]]) 
                        hfall(heapnum, hkey, iheap, heap, 1);

                /*------------------------------------------------+
                |  generate fillin and update elements           */

                for (k=kL[i]; k<kL[i+1]; k++) {
                    row2 = iL[k];
                    tag++;
                    for (kk=0; kk<degBt[row2]; kk++) {
                        iwork[ Bt[row2][kk].i] = tag; /* tag these columns */
                        iwork2[Bt[row2][kk].i] = kk;  /* say where they are */
                    }
                    for (kk=kUt[i]; kk<kUt[i+1]; kk++) {
                        col2 = iUt[kk];
                        if ( iwork[col2] == tag ) {
                            Bt[row2][iwork2[col2]].d -= L[k]*Ut[kk]/diagU[i];
                        } else {
                            deg = degBt[row2];
                            REALLOC( Bt[row2], deg+1, VALIND );
                            Bt[row2][deg].i = col2;
                            Bt[row2][deg].d = -L[k]*Ut[kk]/diagU[i];
                            degBt[row2]++;
                        }
                    }
                }

                for (k=kUt[i]; k<kUt[i+1]; k++) {
                    col2 = iUt[k];
                    tag++;
                    for (kk=0; kk<degB[col2]; kk++) {
                        iwork[ B[col2][kk].i] = tag; /* tag these rows */
                        iwork2[B[col2][kk].i] = kk;  /* say where they are */
                    }
                    for (kk=kL[i]; kk<kL[i+1]; kk++) {
                        row2 = iL[kk];
                        if ( iwork[row2] == tag ) {
                            B[col2][iwork2[row2]].d -= L[kk]*Ut[k]/diagU[i];
                        } else {
                            deg = degB[col2];
                            REALLOC( B[col2], deg+1, VALIND );
                            B[col2][deg].i = row2;
                            B[col2][deg].d = -L[kk]*Ut[k]/diagU[i];
                            degB[col2]++;
                        }
                    }
                }

                /*------------------------------------------------+
                |  adjust heap                                   */

                for (k=kUt[i]; k<kUt[i+1]; k++) {
                        col2 = iUt[k];
                        if (method == MD) {
                                hkey[col2] = degB[col2];
                        } else {
                                hkey[col2] = col2;
                        }
                        if (hkey[col2]==0) hkey[col2]=m+1;
                        hrise( hkey, iheap, heap, iheap[col2] );
                        hfall( heapnum, hkey, iheap, heap, iheap[col2] ); 
                }

                /*------------------------------------------------+
                |  free space no longer needed                   */

                /*
                FREE(B[col]);  FREE(Bt[row]);
                */
        }
end:
        /*------------------------------------------------+
        |  process dependent rows/cols                   */

        i = rank;
        for (col=0; col<m; col++) {
            if (icolperm[col] == -1) {
                colperm[i] = col;
                icolperm[col] = i;
                i++;
            }
        }

        i = rank;
        for (row=0; row<m; row++) {
            if (irowperm[row] == -1) {
                rowperm[i] = row;
                irowperm[row] = i;
                i++;
            }
        }

        for (i=rank; i<m; i++) {
                kL[i+1] = kL[i];
                kUt[i+1] = kUt[i];
                diagU[i] = 0.0;
        }

        /*------------------------------------------------+
        |  free up space                                 */

        heap++;
        for (i=0; i<m; i++) { FREE( B[i] ); FREE( Bt[i] ); }
        FREE(degB); FREE(degBt); 
        FREE(hkey); FREE(heap); FREE(iheap);
        FREE(iwork); FREE(iwork2); FREE(B); FREE(Bt);

        /*------------------------------------------------+
        |  update "i" arrays to new indices              */

        for (k=0; k<kL[m]; k++) iL[k] = irowperm[iL[k]];
        for (k=0; k<kUt[m]; k++) iUt[k] = icolperm[iUt[k]];

        /*------------------------------------------------+
        |  divide each column of L by diagonal           */

        for (i=0; i<m; i++) {
            for (k=kL[i]; k<kL[i+1]; k++) {
                L[k] /= diagU[i];
            }
        }

        /*---------------------------------------------------------+
        | calculate and print statistics.                         */

        narth = 0.0e0;
        for (i=0; i<m; i++) {
                k = kL[i+1]-kL[i];   narth += (double) k*k;
                k = kUt[i+1]-kUt[i]; narth += (double) k*k;
        }
        narth = narth + 3*kL[m] + 3*kUt[m] + 2*m;

        lnz    = kL[m];
        unz    = kUt[m];
        if (v) {
                printf("%9d   %9d %15.0f ", lnz, unz, narth);
                fflush(stdout);
        }

	if (  Lt== NULL ) {  MALLOC(  Lt, lnz, double); }
		   else   { REALLOC(  Lt, lnz, double); }
	if ( iLt== NULL ) {  MALLOC( iLt, lnz, int); }
	     	   else   { REALLOC( iLt, lnz, int); }
	if ( kLt== NULL ) {  MALLOC( kLt, m+1, int); }
	     	   else   { REALLOC( kLt, m+1, int); }

	if (  U == NULL ) {  MALLOC(  U,  unz, double); }
		   else   { REALLOC(  U,  unz, double); }
	if ( iU == NULL ) {  MALLOC( iU,  unz, int); }
	     	   else   { REALLOC( iU,  unz, int); }
	if ( kU == NULL ) {  MALLOC( kU,  m+1, int); }
	      	   else   { REALLOC( kU,  m+1, int); }

	atnum(m,m,kL, iL, L, kLt,iLt,Lt);
	atnum(m,m,kUt,iUt,Ut,kU, iU, U );

	if ( E_d == NULL ) {
	    MALLOC( E_d, E_N, int );

	    MALLOC( E, E_NZ, double );
	    MALLOC(iE, E_NZ,    int );
	    MALLOC(kE, E_N+2,    int );
	}
	kE[0] = 0;

	endtime = (double) clock();
	cumtime += endtime - starttime;
}

/*-----------------------------------------------------------------+
| Forward/backward solve using LU factorization                    |
| Input:                                                           |
|    m          dimension of array y                               |
|    y          array containing right-hand side                   |
|                                                                  |
|    static global variables (assumed setup by lufac()):           |
|                                                                  |
|    rank       rank of B                                          |
|    kL, iL, L, three array sparse representation of L             |
|    kUt,iUt,Ut three array sparse representation of U transpose   |
|               without its diagonal                               |
|    diagU      diagonal entries of U                              |
|    colperm, icolperm, rowperm, irowperm                          |
|               column and row permutations and their inverses     |
| Output:                                                          |
|                                          -1                      |
|    y          array containing solution B  y                     |
|                                                                  |
|    integer flag indicating whether system is consistent         */

int     bsolve(
	int m, 
	double *sy,
	int *iy,
	int *pny
)
{
        int i, ny=*pny;
        int k, row, consistent=TRUE;
        double beta;
        double eps;

	static double *y=NULL;
	static int  *tag=NULL;
	static int  currtag=1;

	double starttime, endtime;

	if (m==0) {
	    FREE(y); FREE(tag); currtag=1;
	    Gauss_Eta( 0, sy, iy, &ny);
	    return 0;
	}

        

	starttime = (double) clock();

	if (   y  == NULL) CALLOC(   y, m,   double);
	if ( tag  == NULL) CALLOC( tag, m,   int);

	for (k=0; k<ny; k++) {
	    i = irowperm[iy[k]];
	    y[i] = sy[k];
	    tag[i] = currtag;
	    addtree(i);
	}

        if (rank < m) eps = EPSSOL * maxv(sy,ny);

        /*------------------------------------------------------+
        |               -1                                      |
        |       y  <-  L  y                                    */

        for (i=getfirst(); i < rank && i != -1; i=getnext()) {
                beta = y[i];
                for (k=kL[i]; k<kL[i+1]; k++) {
                        row = iL[k];
			if (tag[row] != currtag) {
			    y[row] = 0.0;
			    tag[row] = currtag;
			    addtree(row);
			}
                        y[row] -= L[k]*beta;
                }
        }

        /*------------------------------------------------------+
        |               -1                                      |
        |       y  <-  U  y                                    */

        for (i=getlast(); i >= rank && i != -1; i=getprev()) {
                if ( ABS( y[i] ) > eps ) consistent = FALSE;
                y[i] = 0.0;
        }
        for ( ; i>=0; i=getprev()) {
                beta = y[i]/diagU[i];
                for (k=kU[i]; k<kU[i+1]; k++) {
			row = iU[k];
			if (tag[row] != currtag) {
			    y[row] = 0.0;
			    tag[row] = currtag;
			    addtree(row);
			}
			y[row] -= U[k]*beta;
                }
                y[i] = beta;
        }

	ny = 0;
	for (i=getfirst(); i != -1; i=getnext()) {
	    if ( ABS(y[i]) > EPS ) {
	        sy[ny] = y[i];
	        iy[ny] = colperm[i];
	        ny++;
	    }
	}

	currtag++;
	killtree();

	Gauss_Eta( m, sy, iy, &ny);

	*pny = ny;

	/*************************************************************
	* Update E and save col_out in E_d[e_iter]                   *
	*************************************************************/

	REALLOC( E, MAX(E_NZ, enz+ny), double );
	REALLOC(iE, MAX(E_NZ, enz+ny),    int );
	for (i=0, k=kE[e_iter]; i<ny; i++, k++) {
	   E[k] = sy[i];
	  iE[k] = iy[i];
	}
	enz = k;
	kE[e_iter+1] = enz;

	endtime = (double) clock();
	cumtime += endtime - starttime;

        return consistent;
}

/*-----------------------------------------------------------------+
| Forward/backward solve using LU factorization                    |
| Input:                                                           |
|    m          dimension of array y                               |
|    y          array containing right-hand side                   |
|                                                                  |
|    static global variables (assumed setup by lufac()):           |
|                                                                  |
|    rank       rank of B                                          |
|    kL, iL, L, three array sparse representation of L             |
|    kUt,iUt,Ut three array sparse representation of U transpose   |
|               without its diagonal                               |
|    diagU      diagonal entries of U                              |
|    colperm, icolperm, rowperm, irowperm                          |
|               column and row permutations and their inverses     |
| Output:                                                          |
|                                          -T                      |
|    y          array containing solution B  y                     |
|                                                                  |
|    integer flag indicating whether system is consistent         */

int     btsolve(
        int m,
        double *sy,
        int *iy,
        int *pny
)
{
        int i, ny=*pny;
        int k, row, consistent=TRUE;
        double beta;
        double eps;

        static double *y=NULL;
	static int  *tag=NULL;
	static int  currtag=1;
	double starttime, endtime;

	if (m==0) { 
	    FREE(y); FREE(tag); currtag=1;
	    Gauss_Eta_T ( 0, sy, iy, &ny );
	    return 0;
	}

        

	starttime = (double) clock();

        if (   y == NULL) CALLOC(   y, m, double);
        if ( tag == NULL) CALLOC( tag, m, int);

	Gauss_Eta_T ( m, sy, iy, &ny );

        for (k=0; k<ny; k++) {
	    i = icolperm[iy[k]];
            y[i] = sy[k];
	    tag[i] = currtag;
	    addtree(i);
        }

        if (rank < m) eps = EPSSOL * maxv(sy,ny);

        /*------------------------------------------------------+
        |               -T                                      |
        |       y  <-  U  y                                    */

        for (i=getfirst(); i < rank && i != -1; i=getnext()) {
                beta = y[i]/diagU[i];
                for (k=kUt[i]; k<kUt[i+1]; k++) {
                        row = iUt[k];
			if (tag[row] != currtag) {
			    y[row] = 0.0;
			    tag[row] = currtag;
			    addtree(row);
			}
                        y[row] -= Ut[k]*beta;
                }
                y[i] = beta;
        }
        for (i=getlast(); i >= rank && i != -1; i=getprev()) {
                if ( ABS( y[i] ) > eps ) consistent = FALSE;
                y[i] = 0.0;
        }

        /*------------------------------------------------------+
        |               -T                                      |
        |       y  <-  L  y                                    */

        for ( ; i>=0; i=getprev()) {
                beta = y[i];
                for (k=kLt[i]; k<kLt[i+1]; k++) {
		    row = iLt[k];
		    if (tag[row] != currtag) {
			y[row] = 0.0;
			tag[row] = currtag;
			addtree(row);
		    }
		    y[row] -= Lt[k]*beta;
                }
        }

	ny = 0;
	for (i=getfirst(); i != -1; i=getnext()) {
	    if ( ABS(y[i]) > EPS ) {
	        sy[ny] = y[i];
	        iy[ny] = rowperm[i];
	        ny++;
	    }
	}
	*pny = ny;

	currtag++;
	killtree();

	endtime = (double) clock();
	cumtime += endtime - starttime;

        return consistent;
}

void lu_clo()
{
        FREE( rowperm ); FREE( irowperm );
        FREE( colperm ); FREE( icolperm );
        FREE( L ); FREE( iL ); FREE( kL ); 
	FREE( U ); FREE( iU ); FREE( kU );
	FREE( Lt ); FREE( iLt ); FREE( kLt );
        FREE( Ut ); FREE( iUt ); FREE( kUt );
        FREE( diagU );
	FREE( E_d );
	FREE(  E );
	FREE( iE );
	FREE( kE );
	e_iter = 0;
	enz=0;
	rank;
	cumtime = 0.0;
	ocumtime= 0.0;
}

/*-----------------------------------------------------------------+
| Forward/backward solve using LU factorization                    |
| Input:                                                           |
|    m          dimension of array y                               |
|    y          array containing right-hand side                   |
|                                                                  |
|    static global variables (assumed setup by lufac()):           |
|                                                                  |
|    rank       rank of B                                          |
|    kL, iL, L, three array sparse representation of L             |
|    kUt,iUt,Ut three array sparse representation of U transpose   |
|               without its diagonal                               |
|    diagU      diagonal entries of U                              |
|    colperm, icolperm, rowperm, irowperm                          |
|               column and row permutations and their inverses     |
| Output:                                                          |
|                                          -1                      |
|    y          array containing solution B  y                     |
|                                                                  |
|    integer flag indicating whether system is consistent         */

int     dbsolve(int m, double *y)
{
        int i;
        int k, row, consistent=TRUE;
        double beta, *dwork;
        double eps;

        double starttime, endtime;

	starttime = (double) clock();

        MALLOC (dwork,m,double);

        if (rank < m) eps = EPSSOL * maxv(y,m);
        for (i=0; i<m; i++) dwork[i] = y[i];
        for (i=0; i<m; i++) y[irowperm[i]] = dwork[i];

        /*------------------------------------------------------+
        |               -1                                      |
        |       y  <-  L  y                                    */

        for (i=0; i<rank; i++) {
                beta = y[i];
                for (k=kL[i]; k<kL[i+1]; k++) {
                        row = iL[k];
                        y[row] -= L[k]*beta;
                }
        }

        /*------------------------------------------------------+
        |               -1                                      |
        |       y  <-  U  y                                    */

        for (i=m-1; i>=rank; i--) {
                if ( ABS( y[i] ) > eps ) consistent = FALSE;
                y[i] = 0.0;
        }
        for (i=rank-1; i>=0; i--) {
                beta = y[i];
                for (k=kUt[i]; k<kUt[i+1]; k++) {
                        beta -= Ut[k]*y[iUt[k]];
                }
                y[i] = beta/diagU[i];
        }

        for (i=0; i<m; i++) dwork[i] = y[i];
        for (i=0; i<m; i++) y[colperm[i]] = dwork[i];

        FREE(dwork);

	endtime = (double) clock();
	cumtime += endtime - starttime;

        return consistent;
}

/*****************************************************************
*  Gaussian elimination for Eta transformations which will solve *
* each system E ... E  d = a for d.				 *
*              1     s                                           *
*****************************************************************/

static void Gauss_Eta( 
    int m, 
    double *dx_B, 
    int *idx_B, 
    int *pndx_B 
)
{
    int i, j, k, col, kcol, ii, ndx_B=*pndx_B;
    double temp;
    static double *a=NULL;
    static int  *tag=NULL;
    static int *link=NULL;
    static int  currtag=1;

    if (m==0) { 
        FREE(a); FREE(tag); link--; FREE(link); 
	currtag=1;
        return;
    }

    if (  a  == NULL) CALLOC(  a,  m,   double);
    if ( tag == NULL) CALLOC( tag, m,   int);
    if (link == NULL) {CALLOC(link, m+2, int); link++;}

    if (e_iter <= 0) return;

    ii = -1;
    for (k=0; k<ndx_B; k++) {
	i = idx_B[k];
	a[i] = dx_B[k];
	tag[i] = currtag;
	link[ii] = i;
	ii = i;
    }
    for (j=0; j<e_iter; j++) {
	col = E_d[j];

	for (k=kE[j]; k<kE[j+1]; k++) {
	    i = iE[k];
	    if (tag[i] != currtag) {
		a[i] = 0.0;
		tag[i] = currtag;
	        link[ii] = i;
	        ii = i;
	    }
	    if (i == col) kcol = k;
	}
    
	temp = a[col]/E[kcol];
	if (temp != 0.0) {
	    for (k=kE[j]; k<kcol; k++) {
	        i = iE[k];
	        a[i] -= E[k] * temp;
	    }
            a[col] = temp;
            for (k=kcol+1; k<kE[j+1]; k++) {
	        i = iE[k];
	        a[i] -= E[k] * temp;	
	    }
	}
    }
    link[ii] = m;
    currtag++;

    k = 0;
    for (i=link[-1]; i<m; i=link[i]) {
	if ( ABS(a[i]) > EPS1 ) {
	     dx_B[k] = a[i];
	    idx_B[k] = i;
	    k++;
	}
    }
    *pndx_B = k;
}


/*****************************************************************
*  Gaussian elimination for Eta transformations which will solve *
* each system B y = c for y                             	 *
*****************************************************************/

static void Gauss_Eta_T( 
    int m, 
    double *vec, 
    int *ivec, 
    int *pnvec 
)
{
    int i, j, k, kk, kkk, col, nvec=*pnvec;
    double temp;
    static double *a=NULL;
    static int  *tag=NULL;
    static int  currtag=1;

    if (m==0) {
        FREE(a); FREE(tag); currtag=1;
        return;
    }

    if (  a  == NULL) CALLOC(  a,  m,   double);
    if ( tag == NULL) CALLOC( tag, m,   int);

    for (j=e_iter-1; j>=0; j--) {
	col = E_d[j];

        for (k=0; k<nvec; k++) {
	    i = ivec[k];
	    if (i == col) kk = k;
	    a[i] = vec[k];
	    tag[i] = currtag;
	}

	if (tag[col] != currtag) {
	    vec[nvec] = 0.0;
	    ivec[nvec] = col;
	    kk = nvec;
	    nvec++;
	    a[col] = 0.0;
	    tag[col] = currtag;
	}
	temp = vec[kk];
	for (k=kE[j]; k<kE[j+1]; k++) {
	    i = iE[k];
	    if (i == col) kkk = k;
	    if (tag[i] == currtag) {
	        if (i != col) {
		    temp -= E[k]*a[i];
		}
	    }
	}
        currtag++;

	vec[kk] = temp/E[kkk];
	*pnvec = nvec;
    }
}
