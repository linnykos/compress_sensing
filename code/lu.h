void lufac( 
    int m, 
    int *kA, 
    int *iA, 
    double *A, 
    int *basics, 
    int v
);

int refactor( 
    int m, 
    int *kA, 
    int *iA, 
    double *A, 
    int *basics, 
    int col_out,
    int v
);

int dbsolve(
    int m, 
    double *y
);

int dbtsolve(
    int m, 
    double *y
);

int btsolve(
    int m,
    double *sy,
    int *iy,
    int *pny
);

int bsolve(
    int m, 
    double *sy,
    int *iy,
    int *pny
);

void lu_clo();
