param n := 1122;
param d := 20022;

param X {i in 1..n, j in 1..d};
param beta0 {j in 1..d}, default 0;
param y {i in 1..n} := sum {j in 1..d} X[i,j]*beta0[j];

var beta {j in 1..d};
var absbeta {j in 1..d} >= 0;

minimize sumabsbeta: sum {j in 1..d} absbeta[j];

subject to equality {i in 1..n}:
    sum {j in 1..d} X[i,j]*beta[j] = y[i];

subject to absup {j in 1..d}:  beta[j] <= absbeta[j];
subject to absdn {j in 1..d}: -absbeta[j] <= beta[j];

data;

param X: include 'U';
read {j in 1..d} beta0[j] < x0.dat;

option loqo_options "verbose=0 timing=1 sigfig=6";
solve;

display beta0, beta;
