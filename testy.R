library(rAMPL)

path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow_proby/fixprec/")
source(file.path(path, "functions/fixprec_WW.R"))
source(file.path(path, "functions/ampl/ampl_fixprec_fun.R"))
source(file.path(path, "functions/subsets.R"))
model <- file.path(path, "functions/ampl/ampl_fixprec.mod")
setwd(path)

# Przyklad 1 (3d, 1>,7>) ----

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

x_17$T_eigenval
x_ampl$Topt
x_17$n_ih
x_ampl$n_ih

which(x_ampl$n_ih == N)
which(x_17$n_ih == N)

# Przyklad 2 (9d, spctr vs ampl) ----

# Inne zbiory active, wiezy spelnione w obu przypadkach ale ampl mniej dokladnie.

data_ampl <- file.path(path, "functions/ampl/ampl_fixprec_9d_2.dat")
data <- ampl_readData(data_ampl, model = model)
J <- data$J
N <- data$N
S <- data$S
total <- data$total
kappa <- data$kappa
n <- data$n

x_ampl <- ampl_fixprec(data = data_ampl, model = model)
x_spctr <- rfixprec(n, J, N, S, total, kappa)

check_kkt(x_ampl$n_ih, J, N, S, total, kappa, n, tol = 0.1)
check_kkt(x_spctr, J, N, S, total, kappa, n, tol = 10^-11)

# Przyklad 3 (wiele domen/warstw) ----

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
  kkt_x <- check_kkt(x_spctr, J, N, S, total, kappa, n)
  print(kkt_ampl)
  print(kkt_x)
  kkt_ampl_spctr[[n]] <- c(kkt_ampl, kkt_x)
}

kkt_ampl_spctr_big <- which(
  abs(sapply(kkt_ampl_spctr, function(i) do.call(`/`, as.list(i))) - 1) > 0.01
)
