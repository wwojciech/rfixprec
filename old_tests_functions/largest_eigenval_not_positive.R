# Source functions ----

library(rAMPL)
path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow/rfixprec/")
setwd(path)
source(file.path(path, "old_tests_functions/rfixprec_nmax.R"))
source(file.path(path, "ampl/ampl_fixprec.R"))
model <- file.path(path, "ampl/ampl_fixprec.mod")

# Przyklad (nmax < n < sum(N)) ----

J <- c(5, 2)
strata_names <- rep(seq_along(J), times = J)
N <- c(100, 100, 100, 100, 100, 100, 100)
S <- c(154, 178, 213, 134, 124, 102, 12)
names(N) <- strata_names
names(S) <- strata_names
total <- c(13, 2)
kappa <- c(0.8, 0.2)
(rho <- total * sqrt(kappa))
nmax(J, N, S) - 1 # 603
sum(N) # 700
n <- 695 # albo mozna sobie tez sprawdzic 605

x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh
# 100 100 100 100 96.94939 100 98.05060
x_ampl$Topt # 357.8548
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, tol = 10^-2) # 357.8558
(x_fixprec <- fixprec(n, J, N, S, total, kappa, c(1:4, 6)))
# 100 100 100 100 96.94936 100 98.05064
check_kkt(x_fixprec, J, N, S, total, kappa, n, tol = 10^-13) # 357.8601

# (to sa testy recznego sprawdzenie algorytmu zaproponowanego przez Pana Profesora
# dla n > nmax). Ale tutaj jest dla mnie wazne tylko to, ze:
# largest eigenvalue is not strictly positive

# domena I
setNames(fixprec(n, J, N, S, total, kappa), strata_names)
# largest eigenvalue is not strictly positive, ale oczywiscie eigenvector > 0.
# natomiast, largest eigenvalue > 0, jesli n < nmax (nie wiem czy iff, ale chyba tak
# jak sugeruja wyniki ponizszej petli)
for (n in (nmax(J, N, S) - 1):sum(N)) { # 98 elementow
  print(n)
  fixprec(n, J, N, S, total, kappa)
}

# pozostale kroki prowadzace do alokacji optymalnej
setNames(fixprec(n, J, N, S, total, kappa, 1:3), strata_names) # domena I - full take-max
# domena II (j)
fixprec(n - 500, 2, N[-(1:5)], S[-(1:5)], total[-1], kappa[-1])
# domena I
setNames(fixprec(n, J, N, S, total, kappa, 6), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3)), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3, 1)), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3, 1, 4)), strata_names)
# 100 100 100 100 96.94936 100 98.05064
