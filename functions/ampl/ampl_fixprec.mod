# We could instead have written this model with a set DOMSTRATA of pairs, such that
# (i, h) is a member of DOMSTRATA.  But then, we need some processing for DOMAINS only.

set DOMAINS;
set STRATA {DOMAINS}; # strata labels for each domain

param N {i in DOMAINS, STRATA[i]} > 0 integer;
param S {i in DOMAINS, STRATA[i]} > 0;
param total {DOMAINS} > 0;
param kappa {DOMAINS} > 0;
	check sum {i in DOMAINS} kappa[i] = 1;
param nmax = sum {i in DOMAINS} (sum {h in STRATA[i]} N[i, h] *  S[i, h])^2 / (sum {h in STRATA[i]} N[i, h] *  S[i, h]^2); 
param n > 0, < nmax integer;

param rho {i in DOMAINS} = total[i] * sqrt(kappa[i]); 
param A {i in DOMAINS, h in STRATA[i]} = (N[i, h] * S[i, h]) / rho[i];
param c {i in DOMAINS} = (1/rho[i]^2) * sum {h in STRATA[i]} N[i, h] * S[i, h]^2;

# Optimization problem

var T >= 0;
var x {i in DOMAINS, h in STRATA[i]} >= 0, <= N[i,h];

minimize base_variance: T;

subject to total_sample_size:
  sum {i in DOMAINS, h in STRATA[i]} x[i, h] = n;
  
subject to Ti {i in DOMAINS}:
  sum {h in STRATA[i]} (A[i, h]^2 / x[i, h]) - c[i] = T;
  