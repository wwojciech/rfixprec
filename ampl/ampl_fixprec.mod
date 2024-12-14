# We could instead have written this model with a set DOMSTRATA of pairs, such that
# (d, h) is a member of DOMSTRATA.  But then, we need some processing for DOMAINS only.

set DOMAINS;
set STRATA {DOMAINS}; # strata labels for each domain
set DOMAINS_J default DOMAINS; # domain labels for which x <= N must be satisfied, subset of DOMAINS

param N {d in DOMAINS, STRATA[d]} > 0 integer;
param S {d in DOMAINS, STRATA[d]} > 0;
param total {DOMAINS} > 0;
param kappa {DOMAINS} > 0;
	check: sum {d in DOMAINS} kappa[d] = 1;

param M {d in DOMAINS, h in STRATA[d]} > 0 <= N[d, h] integer, default N[d, h];
# check {d in DOMAINS, h in STRATA[d]}: M[d, h] <= N[d, h];
	
# param nmax = sum {d in DOMAINS} (sum {h in STRATA[d]} N[d, h] *  S[d, h])^2 / (sum {h in STRATA[d]} N[d, h] *  S[d, h]^2); 
param nmax = sum {d in DOMAINS} (sum {h in STRATA[d]} N[d, h]);
# param n > 0, < nmax integer;
param n > 0, <= nmax integer;

param rho {d in DOMAINS} = total[d] * sqrt(kappa[d]); 
param A {d in DOMAINS, h in STRATA[d]} = (N[d, h] * S[d, h]) / rho[d];
param c {d in DOMAINS} = (1/rho[d]^2) * sum {h in STRATA[d]} N[d, h] * S[d, h]^2;

# Optimization problem

var T >= 0;
#var T;
var x {d in DOMAINS, h in STRATA[d]} >= 0;
#var x {d in DOMAINS, h in STRATA[d]} >= 0, <= M[d,h];

minimize base_variance: T;

subject to total_sample_size:
  sum {d in DOMAINS, h in STRATA[d]} x[d, h] = n;
  
subject to Td {d in DOMAINS}:
  sum {h in STRATA[d]} (A[d, h]^2 / x[d, h]) - c[d] = T;
  
subject to strata_sizes {d in DOMAINS_J, h in STRATA[d]}:
  x[d, h] <= M[d, h];
  