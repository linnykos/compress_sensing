param ii;
param jj;

param n1 := 33;
param n2 := 34;
param n := n1*n2;
param d1 := 141;
param d2;
param d := d1*d2;

param X1 {i in 1..n1, j in 1..d1};
param X2 {i in 1..n2, j in 1..d2};
param beta0 {j1 in 1..d1, j2 in 1..d2};
param t0 {i1 in 1..n1, j2 in 1..d2} := sum {j1 in 1..d1} X1[i1,j1]*beta0[j1,j2];
param y  {i1 in 1..n1, i2 in 1..n2} := sum {j2 in 1..d2} t0[i1,j2]*X2[i2,j2];

var beta {j1 in 1..d1, j2 in 1..d2};
var absbeta {j1 in 1..d1, j2 in 1..d2} >= 0;
var t {i1 in 1..n1, j2 in 1..d2};

minimize sumabsbeta: sum {j1 in 1..d1, j2 in 1..d2} absbeta[j1,j2];

subject to equality1 {i1 in 1..n1, j2 in 1..d2}:
    sum {j1 in 1..d1} X1[i1,j1]*beta[j1,j2] = t[i1,j2];

subject to equality2 {i1 in 1..n1, i2 in 1..n2}:
    sum {j2 in 1..d2} t[i1,j2]*X2[i2,j2] = y[i1,i2];

subject to absup {j1 in 1..d1, j2 in 1..d2}:  beta[j1,j2] <= absbeta[j1,j2];
subject to absdn {j1 in 1..d1, j2 in 1..d2}: -absbeta[j1,j2] <= beta[j1,j2];

data;

read d2 < d.dat;
read {i in 1..n1, j in 1..d1} X1[i,j] < A.dat;
read {i in 1..n2, j in 1..d2} X2[i,j] < B.dat;
read {j1 in 1..d1, j2 in 1..d2} beta0[j1,j2] < X0.dat;

option loqo_options "verbose=0 timing=1 sigfig=6";
solve;

display beta0, beta;
