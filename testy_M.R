# Funkcje ----

fixprec <- function(n, J, N, S, total, kappa = NULL, M = N, active = NULL, details = TRUE) {
  if (is.null(kappa)) {
    kappa <- rep(1 / length(unique(J)), length(unique(J)))
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa

  A <- N * S / rep(rho, times = table(J))
  M_over_A <- M / A

  J0 <- J
  N0 <- N
  S0 <- S
  n0 <- n
  M0 <- M

  if (!is.null(active)) {
    b_ij <- numeric(length(J))
    b_ij[active] <- N[active]^2 * S[active]^2 / M[active] - N[active] * S[active]^2
    b_i <- tapply(b_ij, J, sum)
    # (b <- 1 / rho2 * b_i - 1 / sum(rho2) * sum(b_i))
    b <- 1 / rho2 * b_i # pomysl P. Profesora

    J <- J[-active]
    N <- N[-active]
    S <- S[-active]
    n <- n - sum(M[active])
  } else {
    b <- 0
  }

  a <- as.matrix(tapply(N * S, J, sum) / rho)
  c_ <- tapply(N * S^2, J, sum) / rho2 - b # (b = 0 if M = N)
  D.matrix <- (a %*% t(a)) / n - diag(c_, nrow = length(c_))

  eigenout <- eigen(D.matrix, symmetric = TRUE)
  lambda <- eigenout$values[1] # largest eigenvalue
  v <- eigenout$vectors[, 1] # corresponding eigenvector
  if (lambda <= 0) {
    stop("Largest eigenvalue is not strictly positive - solution does not exist!")
  }
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s_i <- n * v / as.numeric(t(a) %*% v)
  A <- (N * S) / rep(rho, times = table(J)) # brakets due to finite-prec arithmetic!
  n_opt <- rep(s_i, table(J)) * A
  if (!is.null(active)) {
    M[-active] <- n_opt
    n_opt <- M
  }

  if (details) {
    check_cnstr <- check_kkt(n_opt, M0, J0, N0, S0, total, kappa, n0, active, s_i)
    A <- (N0 * S0) / rep(rho, times = table(J0))
    list(
      T_lambda = lambda, Ti = kappa * lambda, check_cnstr = check_cnstr,
      D.matrix = D.matrix, eigen = eigenout, M_over_A = M0 / A, s_i = s_i, n_opt = n_opt
    )
  } else {
    list(T_lambda = lambda, n_opt = n_opt)
  }
}

check_kkt <- function(n_opt, M, J, N, S, total, kappa, n, active = NULL, s) {
  Topt <- (1 / sum(total^2 * kappa)) * sum((1 / n_opt - 1 / N) * N^2 * S^2)
  Ti <- (1 / total^2) * tapply((1 / n_opt - 1 / N) * N^2 * S^2, J, sum)
  is_Ti_eq_kappa_T <- Ti - kappa * Topt < 10^-12

  # check if mu_ih is non-negative for active constraints
  is_mu_ih_nonneg <- if (is.null(active)) {
    TRUE
  } else {
    rho <- total * sqrt(kappa)
    A <- (N * S) / rep(rho, times = table(J))
    s <- rep(s, times = table(J))
    all((s >= N / A)[active])
  }

  list(
    `T` = Topt,
    Ti = Ti,
    is_Ti_eq_kappa_T = is_Ti_eq_kappa_T,
    is_sum_nopt_eq_n = sum(n_opt) - n < 10^-12,
    n_ih_leq_M = n_opt <= M,
    is_mu_ih_nonneg = is_mu_ih_nonneg
  )
}

# Przyklad 1 (1, 3 >) ----

# J = 1 1 2 2
# w pierwszej iteracji przekroczone: J[1], J[3]
# wystarczy J[1] dac take-max, i w kolejnej iteracji J[3] nie bedzie przekroczone

J <- c(1, 1, 2, 2) # two sub-populations with 2 strata each.
N <- c(100, 150, 200, 300)
S <- sqrt(c(40, 30, 20, 10))
total <- c(2, 3)
kappa <- c(0.4, 0.6)
n <- 500
M <- c(75, 150, 138, 300)

(x <- fixprec(n, J, N, S, total, kappa, M = M))
# 92.59531 120.28484 139.33392 147.78594
# T = 894.7207
(x_1 <- fixprec(n, J, N, S, total, kappa, M = M, active = 1))
# 75.0000 144.2858 136.2254 144.4888
# T = 944.7178
(x_13 <- fixprec(n, J, N, S, total, kappa, M = M, active = c(1, 3)))
# 75.0000 144.2732 138.0000 142.7268
# T = 944.9731

# Przyklad 2 (domena I zablokowana) ----

## 1+2 ----

J <- c(1, 2, 2)
N <- c(30, 20, 40)
S <- c(10, 5, 5)
total <- c(1, 1)
kappa <- c(0.4, 0.6)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
(A <- N * S / rep(rho, times = table(J)))
n <- 20
M <- c(9, N[2:3])
all(M <= N)
n <= sum(M)
(c_ <- 1 / rho2 * tapply(N * S^2, J, sum))
#    1    2
# 7500 2500

# Opcja 0:
(x <- fixprec(n, J, N, S, total, kappa, M))
# n_opt = 10.6727(>9)  3.1091  6.2182
# T_opt = 13581.82

## Opcja 1: h=1->M ----
# Zablokowanie domeny I

(T_opt <- A[1]^2 / M[1] - c_[1]) # wynika z zablokowania domeny I
# 17500

# Reczne "rozrzucenie" pozostalych n - M_1 na n_21, n_22.
n1 <- n - M[1]
wsp_a <- T_opt + c_[2]
wsp_b <- A[2]^2 - A[3]^2 - (T_opt + c_[2]) * n1
wsp_c <- A[3]^2 * n1
(delta <- wsp_b^2 - 4 * wsp_a * wsp_c)

# rozwiazanie a
(n_22_a <- (-wsp_b - sqrt(delta)) / (2 * wsp_a))
(n_21_a <- n1 - n_22_a)
(n_opt_a <- c(M[1], c(n_21_a, n_22_a)))
# 9.000000 7.232588 3.767412
check_kkt(n_opt_a, M, J, N, S, total, kappa, n)

# rozwiazanie b
(n_22_b <- (-wsp_b + sqrt(delta)) / (2 * wsp_a))
(n_21_b <- n1 - n_22_b)
(n_opt_b <- c(M[1], c(n_21_b, n_22_b)))
# 9.000000 1.267412 9.732588
check_kkt(n_opt_b, M, J, N, S, total, kappa, n)

# Pytanie czy n_21 = v2 * A_21 i n_22 = v2 * A_22?
n_opt_a[-1] / A[-1] # rozne,
n_opt_b[-1] / A[-1] # rozne.
# Spelnia KKT z lambda = 0, lambda_2 = 0, mu_2 = 0, mu_1/lambda_1 = A_11^2/M_11^2

## Pozostale opcje dla M ----

M[2] >= n
M[3] >= n
# powyzsze implikuja, ze odpadaja rowniez 1,2; 1,3; 2;3; 1;2;3

## 2+2 ----

J <- c(1, 1, 2, 2) # two sub-populations with 2 strata each.
N <- c(100, 150, 200, 300)
S <- sqrt(c(40, 30, 20, 10))
total <- c(2, 3)
kappa <- c(0.4, 0.6)
n <- 500
M <- c(90, 110, 200, 300)
(x <- fixprec(n, J, N, S, total, kappa, M = M))
# 92.59531 120.28484 139.33392 147.78594
# T = 894.7207

# innny ----

J <- c(1, 2)
N <- c(10, 15)
S <- sqrt(1 / N^2) # dobrane tak, ze N^2 * S^2 = c(1, 1)
total <- c(1, 1)
kappa <- c(0.4, 0.6)
n <- 10
(x <- fixprec(n, J, N, S, total, kappa))
