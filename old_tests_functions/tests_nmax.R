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
n <- 605
x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh
# 89.95402 100.00000 100.00000  78.27168  72.43051 100.00000  64.34379
x_ampl$Topt # 9974.714
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, tol = 10^-2) # 9974.715
x_fixprec <- fixprec(n, J, N, S, total, kappa, c(2, 3, 6))
# 89.95401 100.00000 100.00000  78.27167  72.43050 100.00000  64.34381
check_kkt(x_fixprec, J, N, S, total, kappa, n, tol = 10^-11) # 9974.718

n <- 695
x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh
# 100 100 100 100 96.94939 100 98.05060
x_ampl$Topt # 357.8548
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, tol = 10^-2) # 357.8558
(x_fixprec <- fixprec(n, J, N, S, total, kappa, c(1:4, 6)))
# 100 100 100 100 96.94936 100 98.05064
check_kkt(x_fixprec, J, N, S, total, kappa, n, tol = 10^-13) # 357.8601

# reczne sprawdzenie algorytmu zaproponowanego przez Pana Profesora.
# domena I
setNames(fixprec(n, J, N, S, total, kappa), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, 1:3), strata_names) # domena I - full take-max
# domena II (j)
fixprec(n - 500, 2, N[-(1:5)], S[-(1:5)], total[-1], kappa[-1])
# domena I
setNames(fixprec(n, J, N, S, total, kappa, 6), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3)), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3, 1)), strata_names)
setNames(fixprec(n, J, N, S, total, kappa, c(6, 2, 3, 1, 4)), strata_names)
# 100 100 100 100 96.94936 100 98.05064

# Przyklad (2d, 2+2h, n > nmax) ----

J <- c(2, 2)
N <- c(100, 200, 150, 40)
S <- c(10, 2, 50, 30)
total <- c(2, 3)
kappa <- c(0.4, 0.6)
(rho <- total * sqrt(kappa))
nmax(J, N, S) - 1 # 364
n < sum(N) # 490
n <- 487
x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh # 100.00000 197.04494 150.00000  39.95506
x_ampl$Topt # 7.498376
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, tol = 10^-4) # 7.498926

round(rfixprec(n, J, N, S, total, kappa), 2)

# rfixprec recznie, blokuje cala domene
fixprec(n, J, N, S, total, kappa)
fixprec(n, J, N, S, total, kappa, 1)
fixprec(n, J, N, S, total, kappa, c(1, 2)) # blad algortymu - cala domena zablokowana!
x_opt <- fixprec(n, J, N, S, total, kappa, c(1, 3)) # alokacja optymalna
check_kkt(x_opt, J, N, S, total, kappa, n, tol = 10^-4) # 7.498441

# warunek nie spelniony na poczatku, ale dla alokacji optymalnej spelniony
n < nmax(J, N, S) # 487 < 365, FALSE
# ale
n - sum(N[c(1, 3)]) < nmax(c(1, 1), N[-c(1, 3)], S[-c(1, 3)]) # 237 < 239, TRUE
fixprec(n, J, N, S, total, kappa, c(1, 3))

# Przyklad (2d, 2+2h) ----

J <- c(2, 2)
N <- c(100, 200, 150, 40)
S <- c(10, 2, 50, 30)
total <- c(2, 3)
kappa <- c(0.4, 0.6)
(rho <- total * sqrt(kappa))
nmax(J, N, S) - 1 # 364
n <- 350
x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh # 100.00000  65.35359 150.00000  34.64641
x_ampl$Topt # 1030.138

round(rfixprec(n, J, N, S, total, kappa), 2)
fixprec(n, J, N, S, total, kappa)

# Przyklad (3d, 1>,7>) ----

J <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
N <- c(140, 110, 135, 190, 200, 40, 70)
S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
(n <- nmax(J, N, S) - 1)

x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_dh # 140 103.60437 132.18060 166.39204 195.95002  19.87296  70
x_ampl$Topt # 67.90424

(x <- fixprec(n, J, N, S, total, kappa))
(x_17 <- fixprec(n, J, N, S, total, kappa, active = c(1, 7)))

x_17
x_ampl$n_dh

which(x_ampl$n_dh == N)
which(x_17 == N)

# Przyklad (sctr vs ampl, v0) ----

D <- 4
Hmax <- 5
J <- sample(1:Hmax, D, replace = TRUE)
(N <- ceiling(rchisq(sum(J), 15)) * 10)
(S <- ceiling(sqrt(rf(sum(J), 15, 1))) * 10)
total <- sample(1:100, D)
kappa <- sample(1:100, D)
kappa <- kappa / sum(kappa)
kappa <- round(kappa, 2)
(nm <- nmax(J, N, S) - 1)

kkt_ampl_spctr <- vector("list", nm)
for (n in nm:1) {
  cat("\t n: ", n, "\n")
  x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
  x_spctr <- rfixprec(n, J, N, S, total, kappa)
  kkt_ampl <- check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, tol = 0.1)
  kkt_spctr <- check_kkt(x_spctr, J, N, S, total, kappa, n)
  print(kkt_ampl)
  print(kkt_spctr)
  kkt_ampl_spctr[[n]] <- c(kkt_ampl, kkt_spctr)
}

kkt_ampl_spctr_big <- which(
  abs(sapply(kkt_ampl_spctr, function(i) do.call(`/`, as.list(i))) - 1) > 0.01
)

# Przyklad (spctr vs ampl) ----

H1 <- 6
H2 <- 6
H3 <- 6
J <- c(H1, H2, H3)
N <- c(sample(100, H1), sample(200, H2), sample(500, H3))
S <- round(sqrt(c(rchisq(H1 + H2, 6), sample(50, H3))), 0)
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
(nm <- nmax(J, N, S) - 2)

for (n in nm:1) {
  if (floor(n / 10) == n / 10) {
    cat("\t n: ", n, "\n")
  }
  # alloc
  x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
  x_spctr <- rfixprec(n, J, N, S, total, kappa)

  # active
  active_ampl <- which(x_ampl$n_dh == N)
  active_spctr <- which(x_spctr == N)
  (d <- setdiff(active_spctr, active_ampl))
  x_ampl$n_dh[d]
  N[d]

  # alloc spectr for active ampl
  x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

  # check the KKT
  check_kkt(x_spctr, J, N, S, total, kappa, n, active = active_spctr, tol = 10^-13)
  kkt_spctr_actampl <- check_kkt(x_spctr_actampl, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-13)

  # check_kkt(x_spctr, J, N, S, total, kappa, n, tol = 10^-11)
  if (length(active_spctr) == 0L) {
    break
  }

  # active set diff
  if (!is.na(kkt_spctr_actampl) && length(d) != 0L) { # zbiory active roznia sie a wiezy sepelnione
    cat("\t n: ", n, "\n")
    print("XXXXXXXXX     STOP     XXXXXXXXX")
    break
  }
}

# Przyklad (3d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_3d_2.dat")
data <- ampl_readData(data_ampl, model = model)
J <- data$J
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
n <- data$n

# alloc
x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr <- rfixprec(n, J, N, S, total, kappa)

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr == N)
(ampl_spctr_act_diff <- setdiff(active_spctr, active_ampl))
x_ampl$n_dh[ampl_spctr_act_diff]
N[ampl_spctr_act_diff]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
# x_ampl$n_dh malo dokladne spelnienie wiezow => AMPL tylko do zidentyfikownia zb. active.
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-5) # 157.2966
check_kkt(x_spctr, J, N, S, total, kappa, n, active = active_spctr, tol = 10^-13, details = TRUE) # 157.2966
check_kkt(x_spctr_actampl, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-13, details = TRUE) # 157.2966

# sprawdzenie n
sum(x_spctr) - n # 0
sum(x_spctr_actampl) - n # 0

# spradwdzenie krok po kroku
# domena 3 nr ma takie same zbiory active
active_ampl
# [1]  3  5 10 11 13 14 15 17
active_spctr
# [1]  3  5 10 11 13 14 15 17
act_d3 <- c(13, 14, 15, 17)
which(fixprec(n, J, N, S, total, kappa, active = act_d3) >= N)
# d1
which(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 3)) >= N)
which(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 3, 5)) >= N)
# d2
which(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 10, 11)) >= N)
# d1
which(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 10, 11, 3)) >= N)
which(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 10, 11, 3, 5)) >= N)
all(fixprec(n, J, N, S, total, kappa, active = c(act_d3, 10, 11, 3, 5)) <= N)
setequal(active_ampl, c(act_d3, 10, 11, 3, 5))

# Przyklad (4d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_4d.dat")
data <- ampl_readData(data_ampl, model = model)
J <- data$J
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
n <- data$n

# alloc
x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr <- rfixprec(n, J, N, S, total, kappa)

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr == N)
(ampl_spctr_act_diff <- setdiff(active_spctr, active_ampl))
x_ampl$n_dh[ampl_spctr_act_diff]
N[ampl_spctr_act_diff]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-1) # 25404.1
check_kkt(x_spctr, J, N, S, total, kappa, n, active = active_spctr, tol = 10^-9) # 25404.09
check_kkt(x_spctr_actampl, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-9) # 25404.09

# sprawdzenie n
sum(x_spctr) - n # -9.094947e-13
sum(x_spctr_actampl) - n # -9.094947e-13

# Przyklad (9d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_9d_2.dat")
data <- ampl_readData(data_ampl, model = model)
J <- data$J
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
n <- data$n

# alloc
x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr <- rfixprec(n, J, N, S, total, kappa)

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr == N)
(ampl_spctr_act_diff <- setdiff(active_spctr, active_ampl))
x_ampl$n_dh[ampl_spctr_act_diff]
N[ampl_spctr_act_diff]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
x_ampl$Topt # 36218
check_kkt(x_ampl$n_dh, J, N, S, total, kappa, n, active = active_ampl, tol = 0.1) # 36218
check_kkt(x_spctr, J, N, S, total, kappa, n, active = active_spctr, tol = 10^-11) # 36218.02
check_kkt(x_spctr_actampl, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-11) # 36218.02

# sprawdzenie n
sum(x_ampl$n_dh) - n # -6.726887e-05
sum(x_spctr) - n # 7.275958e-12
sum(x_spctr_actampl) - n # 7.275958e-12
