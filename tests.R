# Source functions ----

library(rAMPL)
path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow/rfixprec/")
setwd(path)
source(file.path(path, "functions/rfixprec.R"))
source(file.path(path, "ampl/ampl_fixprec.R"))
model <- file.path(path, "ampl/ampl_fixprec.mod")

# Przyklad (2d, 2+2h) ----

J <- c(2, 2)
N <- c(100, 200, 150, 40)
S <- c(10, 2, 50, 30)
total <- c(2, 3)
kappa <- c(0.4, 0.6)
(rho <- total * sqrt(kappa))
nmax(J, N, S) - 1 # 449
n <- 350
x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
x_ampl$n_ih
x_ampl$Topt

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
x_ampl$n_ih # 140 103.60437 132.18060 166.39204 195.95002  19.87296  70
x_ampl$Topt # 67.90424

(x <- fixprec(n, J, N, S, total, kappa))
(x_17 <- fixprec(n, J, N, S, total, kappa, active = c(1, 7)))

x_17
x_ampl$n_ih

which(x_ampl$n_ih == N)
which(x_17 == N)

# Przyklad (sctr vs ampl, v0) ----

I <- 4
Hmax <- 5
J <- sample(1:Hmax, I, replace = TRUE)
(N <- ceiling(rchisq(sum(J), 15)) * 10)
(S <- ceiling(sqrt(rf(sum(J), 15, 1))) * 10)
total <- sample(1:100, I)
kappa <- sample(1:100, I)
kappa <- kappa / sum(kappa)
kappa <- round(kappa, 2)
(nm <- nmax(J, N, S) - 1)

kkt_ampl_spctr <- vector("list", nm)
for (n in nm:1) {
  cat("\t n: ", n, "\n")
  x_ampl <- ampl_fixprec(n, J, N, S, total, kappa, model = model)
  x_spctr <- rfixprec(n, J, N, S, total, kappa)
  kkt_ampl <- check_kkt(x_ampl$n_ih, J, N, S, total, kappa, n, tol = 0.1)
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
  active_ampl <- which(x_ampl$n_ih == N)
  active_spctr <- which(x_spctr == N)
  (d <- setdiff(active_spctr, active_ampl))
  x_ampl$n_ih[d]
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
active_ampl <- which(x_ampl$n_ih == N)
active_spctr <- which(x_spctr == N)
(d <- setdiff(active_spctr, active_ampl))
x_ampl$n_ih[d]
N[d]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
# x_ampl$n_ih malo dokladne spelnienie wiezow => AMPL tylko do zidentyfikownia zb. active.
check_kkt(x_ampl$n_ih, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-5) # 157.2966
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
active_ampl <- which(x_ampl$n_ih == N)
active_spctr <- which(x_spctr == N)
(d <- setdiff(active_spctr, active_ampl))
x_ampl$n_ih[d]
N[d]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
check_kkt(x_ampl$n_ih, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-1) # 25404.1
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
active_ampl <- which(x_ampl$n_ih == N)
active_spctr <- which(x_spctr == N)
(d <- setdiff(active_spctr, active_ampl))
x_ampl$n_ih[d]
N[d]

# alloc spectr for active ampl
x_spctr_actampl <- fixprec(n, J, N, S, total, kappa, active = active_ampl)

# check the KKT
x_ampl$Topt # 36218
check_kkt(x_ampl$n_ih, J, N, S, total, kappa, n, active = active_ampl, tol = 0.1) # 36218
check_kkt(x_spctr, J, N, S, total, kappa, n, active = active_spctr, tol = 10^-11) # 36218.02
check_kkt(x_spctr_actampl, J, N, S, total, kappa, n, active = active_ampl, tol = 10^-11) # 36218.02

# sprawdzenie n
sum(x_ampl$n_ih) - n # -6.726887e-05
sum(x_spctr) - n # 7.275958e-12
sum(x_spctr_actampl) - n # 7.275958e-12
