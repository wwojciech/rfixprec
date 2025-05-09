# Source functions ----

library(rAMPL)
path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow/rfixprec/")
setwd(path)
source(file.path(path, "functions/rfixprec.R"))
source(file.path(path, "functions/subsets.R"))
source(file.path(path, "ampl/ampl_fixprec.R"))
model <- file.path(path, "ampl/ampl_fixprec.mod")

# Przyklad (if any(H_dind[U] %in% J) == TRUE) ----

H_counts <- c(2, 2) # two domains with 2 strata each.
N <- c(140, 110, 135, 190)
S <- sqrt(c(180, 20, 5, 4))
total <- c(2, 3)
kappa <- c(0.4, 0.6)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
sum(N) # 575
n <- sum(N) - 1
U <- 1:3
sum(N[U]) < n # TRUE

# optymalna alokacja
rfixprec(n, H_counts, N, S, rho, rho2)
# 140.0000 109.9438 135.0000 189.0562
# U = { 1 3 }

rfixprec(n, H_counts, N, S, rho, rho2, U = U)
# [1] "fixprecact debug: U = { 1 2 3 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# [1] "fixprecact debug: U = { 1 2 3 3 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# [1] "fixprecact debug: U = { 1 2 3 1 2 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# Error in fixprecact(n = n, domain_count = length(H), H_dind = H_dind,  :
#                       U does not meet the precondition

# no error, but x not optimal, even sum(x) != n
rfixprec(n, H_counts, N, S, rho, rho2, U = 1)
# [1] "fixprecact debug: U = { 1 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# [1] "fixprecact debug: U = { 1 3 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# [1] "fixprecact debug: U = { 1 1 }"
# [1] "fixprecact debug: H_ind = ( 1 1 2 2 )"
# [1] 140.00000  97.71259  86.89808 109.38934

# inny trywialny przyklad z 1 domena

H_counts <- 4
H <- H_cnt2glbidx(H_counts)
N <- c(140, 110, 135, 190)
S <- sqrt(c(180, 20, 5, 4))
total <- 5
kappa <- 1
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
sum(N) # 575
n <- sum(N) - 1
U <- 1:3
sum(N[U]) < n # TRUE

# optymalna alokacja
rfixprec(n, H_counts, N, S, rho, rho2)
# 140 110 135 189
# U = { 1 2 3 }

rfixprec(n, H_counts, N, S, rho, rho2, U = U)
# [1] "fixprecact debug: U = { 1 2 3 }"
# [1] "fixprecact debug: H_ind = ( 1 1 1 1 )"
# [1] "fixprecact debug: U = { 1 2 3 1 2 3 }"
# [1] "fixprecact debug: H_ind = ( 1 1 1 1 )"
# Error in fixprecact(n = n, domain_count = length(H), H_dind = H_dind,  :
#                       U does not meet the precondition

# Przyklad (nmax < n < sum(N), J = 1) ----

H_counts <- c(5, 2)
H_names <- rep(seq_along(H_counts), times = H_counts)
N <- c(100, 100, 100, 100, 100, 100, 100)
S <- c(154, 178, 134, 213, 124, 102, 12)
names(N) <- H_names
names(S) <- H_names
total <- c(13, 2)
kappa <- c(0.8, 0.2)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
nmax(H_counts, N, S) - 1 # 603
sum(N) # 700

n <- 623
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model, J = 1)
x_ampl$n_dh
# 100.00000 100.00000 100.00000 100.00000  99.79400 110.23695  12.96905
x_ampl$Topt # 23.47558
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, J = 1, tol_max = -3)
# 23.47645
rfixprec(n, H_counts, N, S, rho, rho2, J = 1)

x_fixprec <- fixprecact(n, H_counts, N, S, rho, rho2, rho2, U = 1:4)
x_fixprec
# 100.00000 100.00000 100.00000 100.00000  99.79400 110.23695  12.96905
check_kkt(x_fixprec, H_counts, N, S, total, kappa, n, J = 1, tol_max = -9)
# 23.476

# Testy dla wszystkich podzbiorow (poza cala domena)

all_subsets <- subsets_domains(H_counts)
kkt <- sapply(all_subsets, function(i) {
  x <- fixprecact(n, H_counts, N, S, rho, rho2, rho2, U = i)
  check_kkt(x, H_counts, N, S, total, kappa, n, J = 1, tol_max = -3)
})
kkt
kkt[which.min(kkt)]
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model, J = 1)
x_ampl$n_dh
x_ampl$Topt

## Algorytm zaproponowany przez P. Profesora.

# a) reczne wyrzucanie warstw - tylko fixprec

# domena I
n <- 640
setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(2, 4)), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(2, 4, 1, 3)), H_names)
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(2, 4, 1, 3, 5)), H_names))
check_kkt(x, H_counts, N, S, total, kappa, n, J = 1, details = TRUE)

(x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model, J = 1))
# Restoration Phase Failed.
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, J = 1, details = TRUE)

# Przyklad (nmax < n < sum(N)) ----

H_counts <- c(5, 2)
H_names <- rep(seq_along(H_counts), times = H_counts)
N <- c(100, 100, 100, 100, 100, 100, 100)
S <- c(154, 178, 134, 213, 124, 102, 12)
names(N) <- H_names
names(S) <- H_names
total <- c(13, 2)
kappa <- c(0.8, 0.2)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
nmax(H_counts, N, S) - 1 # 603
sum(N) # 700

n <- 605
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
x_ampl$n_dh
# 89.95402 100.00000 78.27168 100.00000 72.43051 100.00000  64.34379
x_ampl$Topt # 9974.714
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, tol_max = -1) # 9974.715
(x_fixprec <- fixprecact(n, H_counts, N, S, rho, rho2, U = c(2, 4, 6)))
# 89.95401 100.00000  78.27167 100.00000  72.43050 100.00000  64.34381
check_kkt(x_fixprec, H_counts, N, S, total, kappa, n, tol_max = -11) # 9974.718

n <- 695
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
x_ampl$n_dh
# 100 100 100 100 96.94939 100 98.05060
x_ampl$Topt # 357.8548
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, tol_max = -1)
# 357.8558
(x_fixprec <- fixprecact(n, H_counts, N, S, rho, rho2, U = c(1:4, 6)))
# 100 100 100 100 96.94936 100 98.05064
check_kkt(x_fixprec, H_counts, N, S, total, kappa, n, tol_max = -11)
# 357.8601

## Algorytm zaproponowany przez P. Profesora.

# a) reczne wyrzucanie warstw - tylko fixprec

# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(1:2, 4)), H_names) # domena I - full take-max
# domena II (j)
fixprecact(n - 500, 2, N[-(1:5)], S[-(1:5)], rho[-1], rho2[-1])
# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 6), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4)), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4, 1)), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4, 1, 3)), H_names)
# 100.00000 100.00000 100.00000 100.00000  96.94936 100.00000  98.05064

# b) po zmianach w fixprec, zeby nie trzeba bylo recznie wyrzucac warstw w N, S, ...

# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(1, 2, 4)), H_names) # domena I - full take-max
# domena II
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(1, 2, 4, 3, 5)), H_names)
# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 6), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4)), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4, 1)), H_names)
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2, 4, 1, 3)), H_names)
# 100.00000 100.00000 100.00000 100.00000  96.94936 100.00000  98.05064

microbenchmark::microbenchmark(
  iter = rfixprec_iter(n, H_counts, N, S, rho, rho2),
  rec = rfixprec(n, H_counts, N, S, rho, rho2),
  recmain = rfixprec_main(n, H_counts, N, S, total, kappa),
  times = 100,
  check = "identical",
  unit = "us" # microseconds
)

# Przyklad (3d) ----

H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
H_names <- rep(seq_along(H_counts), times = H_counts)
N <- c(140, 110, 135, 190, 200, 40, 70)
S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
names(N) <- H_names
names(S) <- H_names
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
sum(N) # 885
nmax(H_counts, N, S) - 1
n <- 880
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
x_ampl$n_dh # 140.00000 109.62191 135.00000 188.30800 200.00000  37.07008  70.00000
x_ampl$Topt # 3.793779

# rfixprec rec vs iter
all(sapply(1:(sum(N) - 1), function(n) {
  all(rfixprec(n, H_counts, N, S, rho, rho2) == rfixprec_iter(n, H_counts, N, S, rho, rho2))
}))

# recznie

# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names) >= N
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(1)), H_names) >= N
# domena I full

# domena II
setNames(
  fixprecact(
    n - sum(N[1:2]), c(2, 3), N[-(1:2)], S[-(1:2)], rho[-1], rho2[-1]
  ), H_names[-c(1, 2)]
) >= N[-(1:2)]
setNames(
  fixprecact(n - sum(N[1:2]), c(2, 3), N[-(1:2)], S[-(1:2)], rho[-1], rho2[-1], U = 1),
  H_names[-c(1, 2)]
) >= N[-(1:2)]
# domena II full

# domena III
setNames(
  fixprecact(n - sum(N[1:4]), 3, N[-(1:4)], S[-(1:4)], rho[-(1:2)], rho2[-(1:2)]),
  H_names[-c(1:4)]
) >= N[5:7]

setNames(
  fixprecact(
    n - sum(N[1:4]), 3, N[-(1:4)], S[-(1:4)], rho[-(1:2)], rho2[-(1:2)],
    U = c(1, 3)
  ),
  H_names[-c(1:4)]
) >= N[5:7]
# 1 i 3 z domeny III - active

# domena II
setNames(
  fixprecact(
    n - sum(N[1:2]), c(2, 3), N[-(1:2)], S[-(1:2)], rho[-1], rho2[-1],
    U = c(3, 5)
  ),
  H_names[-c(1, 2)]
) >= N[3:7]
setNames(
  fixprecact(
    n - sum(N[1:2]), c(2, 3), N[-(1:2)], S[-(1:2)], rho[-1], rho2[-1],
    U = c(3, 5, 1)
  ),
  H_names[-c(1, 2)]
) >= N[3:7]
# alokacja w domenach II i III - poprawna

# domena I
setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(3, 5, 1) + 2), H_names)
N[1:2]
(x_alloc <- setNames(
  fixprecact(n, H_counts, N, S, rho, rho2, U = c(c(3, 5, 1) + 2, 1)),
  H_names
)
)
x_alloc <= N
sum(x_alloc)
# wszystkie domeny OK

# Przyklad (spctr vs ampl, v0) ----

D <- 4
Hmax <- 5
H_counts <- sample(1:Hmax, D, replace = TRUE)
(N <- ceiling(rchisq(sum(H_counts), 15)) * 10)
(S <- ceiling(sqrt(rf(sum(H_counts), 15, 1))) * 10)
total <- sample(1:100, D)
kappa <- sample(1:100, D)
kappa <- kappa / sum(kappa)
kappa <- round(kappa, 2)
kappa[D] <- kappa[D] + (1 - sum(kappa)) # so that is sums to 1
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
nmax(H_counts, N, S) - 1
sum(N)

kkt_ampl_spctr <- vector("list", sum(N))
for (n in sum(N):1) {
  cat("\t n: ", n, "\n")
  x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
  x_spctr_iter <- rfixprec_iter(n, H_counts, N, S, rho)
  kkt_ampl <- check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, tol_max = -1)
  kkt_spctr <- check_kkt(x_spctr_iter, H_counts, N, S, total, kappa, n)
  print(kkt_ampl)
  print(kkt_spctr)
  kkt_ampl_spctr[[n]] <- c(kkt_ampl, kkt_spctr)
}

kkt_ampl_spctr_big <- which(
  abs(sapply(kkt_ampl_spctr, function(i) do.call(`/`, as.list(i))) - 1) > 0.01
)
kkt_ampl_spctr[kkt_ampl_spctr_big]

# rfixprec rec vs iter
all(sapply(1:(sum(N) - 1), function(n) {
  all(rfixprec_main(n, H_counts, N, S, total, kappa) == rfixprec_iter(n, H_counts, N, S, rho, rho2))
}))

# Przyklad (spctr vs ampl) ----

H1 <- 6
H2 <- 6
H3 <- 6
H_counts <- c(H1, H2, H3)
N <- c(sample(100, H1), sample(200, H2), sample(500, H3))
S <- round(sqrt(c(rchisq(H1 + H2, 6), sample(50, H3))), 0)
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
nmax(H_counts, N, S) - 1
sum(N)

for (n in (sum(N) - 1):1) {
  if (floor(n / 10) == n / 10) {
    cat("\t n: ", n, "\n")
  }
  # alloc
  x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
  x_spctr <- rfixprec(n, H_counts, N, S, rho, rho2)

  # active
  active_ampl <- which(x_ampl$n_dh == N)
  active_spctr <- which(x_spctr == N)
  (d <- setdiff(active_spctr, active_ampl))
  x_ampl$n_dh[d]
  N[d]

  # alloc spectr for active ampl
  x_spctr_actampl <- fixprecact(n, H_counts, N, S, rho, rho2, U = active_ampl)

  # check the KKT
  check_kkt(x_spctr, H_counts, N, S, total, kappa, n, active = active_spctr, tol_max = -13)
  kkt_spctr_actampl <- check_kkt(
    x_spctr_actampl, H_counts, N, S, total, kappa, n,
    active = active_ampl, tol_max = -13
  )

  # check_kkt(x_spctr, H_counts, N, S, total, kappa, n, tol_max = -11)
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

# rfixprec rec vs iter
all(sapply(((sum(N) - 1):1), function(n) {
  print(n)
  all(rfixprec_main(n, H_counts, N, S, total, kappa) == rfixprec_iter(n, H_counts, N, S, rho, rho2))
}))

# Przyklad (3d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_3d_2.dat")
data <- ampl_readData(data_ampl, model = model)
H_counts <- data$H_counts
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
n <- data$n

# alloc
x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr_iter <- rfixprec_iter(n, H_counts, N, S, rho, rho2)

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr_iter == N)
setequal(active_spctr, active_ampl)

# alloc spectr for active ampl
x_spctr_actampl <- fixprecact(n, H_counts, N, S, rho, rho2, U = active_ampl)

# check the KKT
# x_ampl$n_dh malo dokladne spelnienie wiezow => AMPL tylko do zidentyfikownia zb. active.
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, active = active_ampl, tol_max = -4)
# 157.2966
check_kkt(x_spctr_iter, H_counts, N, S, total, kappa, n, active = active_spctr, details = TRUE)
# 157.2966
check_kkt(x_spctr_actampl, H_counts, N, S, total, kappa, n, active = active_ampl, details = TRUE)
# 157.2966

# sprawdzenie n
sum(x_spctr_iter) - n # 0
sum(x_spctr_actampl) - n # 0

# spradwdzenie krok po kroku
# domena 3 nr ma takie same zbiory active
active_ampl
# [1]  3  5 10 11 13 14 15 17
active_spctr
# [1]  3  5 10 11 13 14 15 17
act_d3 <- c(13, 14, 15, 17)
which(fixprecact(n, H_counts, N, S, rho, rho2, U = act_d3) >= N)
# d1
which(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 3)) >= N)
which(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 3, 5)) >= N)
# d2
which(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 10, 11)) >= N)
# d1
which(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 10, 11, 3)) >= N)
which(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 10, 11, 3, 5)) >= N)
all(fixprecact(n, H_counts, N, S, rho, rho2, U = c(act_d3, 10, 11, 3, 5)) <= N)
setequal(active_ampl, c(act_d3, 10, 11, 3, 5))

# rfixprec rec vs iter
all(sapply(((sum(N) - 1):1), function(n) {
  print(n)
  all(rfixprec_main(n, H_counts, N, S, total, kappa) == rfixprec_iter(n, H_counts, N, S, rho, rho2))
}))

# Przyklad (4d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_4d.dat")
data <- ampl_readData(data_ampl, model = model)
H_counts <- data$H_counts
H <- data$H
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
rho <- data$rho
rho2 <- data$rho2
n <- data$n

# alloc
x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr_iter <- rfixprec_iter(n, H_counts, N, S, rho, rho2)
x_spctr <- rfixprec_main(n, H_counts, N, S, total, kappa)

for (n in (sum(N) - 200):(sum(N) - 400)) {
  print(n)
  print(
    microbenchmark::microbenchmark(
      iter = rfixprec_iter(n, H_counts, N, S, rho, rho2),
      rec = rfixprec_main(n, H_counts, N, S, total, kappa),
      times = 5,
      check = "identical",
      unit = "ms" # microseconds
    )
  )
}

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr == N)
setequal(active_ampl, active_spctr)

# alloc spectr for active ampl
x_spctr_actampl <- fixprecact(n, H_counts, N, S, rho, rho2, U = active_ampl)
x_spctr_actampl1 <- fixprecact(n, H_counts, N, S, rho, rho2, U = active_ampl)

# check the KKT
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, active = active_ampl, tol_max = -1)
# 25404.1
check_kkt(x_spctr_iter, H_counts, N, S, total, kappa, n, active = active_spctr, tol_max = -8)
# 25404.09
check_kkt(x_spctr_actampl, H_counts, N, S, total, kappa, n, active = active_ampl, tol_max = -8)
# 25404.09

# sprawdzenie n
sum(x_spctr_iter) - n # -9.094947e-13
sum(x_spctr_actampl) - n # -9.094947e-13
sum(x_spctr_actampl1) - n # -9.094947e-13

# rfixprec rec vs iter
all(sapply(((sum(N) - 1):1), function(n) {
  print(n)
  all(rfixprec_main(n, H_counts, N, S, total, kappa) == rfixprec_iter(n, H_counts, N, S, rho, rho2))
}))

# Przyklad (9d, spctr vs ampl) ----

data_ampl <- file.path(path, "ampl/data/ampl_fixprec_9d_2.dat")
data <- ampl_readData(data_ampl, model = model)
H_counts <- data$H_counts
H <- data$H
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
rho <- data$rho
rho2 <- data$rho2
data$nmax # 16607.84
n <- data$n # 16494
sum(N) # 44060

# alloc
# n <- 35000
x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)
# x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr_iter <- rfixprec_iter(n, H_counts, N, S, rho, rho2)
x_spctr <- rfixprec_main(n, H_counts, N, S, total, kappa)

# any differences?
sum(x_ampl$n_dh - x_spctr)
sum(x_spctr_iter - x_spctr)

# active
active_ampl <- which(x_ampl$n_dh == N)
active_spctr <- which(x_spctr == N)
setequal(active_spctr, active_ampl)

# alloc spectr for active ampl
x_spctr_actampl <- fixprecact(n, H_counts, N, S, rho, rho2, U = active_ampl)

# check the KKT
x_ampl$Topt # 36218
check_kkt(x_ampl$n_dh, H_counts, N, S, total, kappa, n, active = active_ampl, tol_max = -1)
# 36218
check_kkt(x_spctr, H_counts, N, S, total, kappa, n, active = active_spctr, tol_max = -9)
# 36218.02
check_kkt(x_spctr_actampl, H_counts, N, S, total, kappa, n, active = active_ampl, tol_max = -9)
# 36218.02

# sprawdzenie n
sum(x_ampl$n_dh) - n # -6.726887e-05
sum(x_spctr) - n # 7.275958e-12
sum(x_spctr_iter) - n # 7.275958e-12
sum(x_spctr_actampl) - n # 7.275958e-12

# rfixprec rec vs iter
all(sapply(1:10000, function(n) {
  print(n)
  x_rec <- rfixprec_main(n, H_counts, N, S, total, kappa)
  x_iter <- rfixprec_iter(n, H_counts, N, S, rho, rho2)
  all(x_rec == x_iter)
}))

# sprawdzenie dla wielu n
for (n in (sum(N) - 1):(sum(N) - 1000)) {
  print(n)
  print(
    microbenchmark::microbenchmark(
      iter = round(rfixprec_iter(n, H_counts, N, S, rho, rho2), 1),
      red = round(rfixprec_modreduce(n, H_counts, N, S, rho, rho2), 1),
      rec = round(rfixprec_main(n, H_counts, N, S, total, kappa), 1),
      ampl = round(ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)$n_dh, 1),
      times = 1,
      check = "identical"
      # unit = "ms" # microseconds
    )
  )
}

# Przyklad (3d, basic, microbenchmark) ----

H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
N <- c(140, 110, 135, 190, 200, 40, 70)
S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
n <- 854

for (n in (sum(N) - 1):1) {
  print(n)
  print(
    microbenchmark::microbenchmark(
      iter = round(rfixprec_iter(n, H_counts, N, S, rho, rho2), 1),
      red = round(rfixprec_modreduce(n, H_counts, N, S, rho, rho2), 1),
      rec = round(rfixprec(n, H_counts, N, S, rho, rho2), 1),
      ampl = round(ampl_fixprec(n, H_counts, N, S, total, kappa, model = model)$n_dh, 1),
      times = 1,
      check = "identical",
      unit = "ms" # microseconds
    )
  )
}

# for n = 8, 9 181, ampl -> Restoration Phase Failed.
