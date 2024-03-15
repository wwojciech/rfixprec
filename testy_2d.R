library(rAMPL)
path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow_proby/fixprec/")
source(file.path(path, "functions/ampl/ampl_fixprec_fun.R"))
model <- file.path(path, "functions/ampl/ampl_fixprec.mod")

#' List with all possible subsets (of lengths len) of x.
#' @param x source subset
#' @param len vector with sizes of subsets `x`
subsets <- function(x, len = NULL) {
  if (is.null(len)) {
    len <- seq_len(length(x) - 1)
  }
  do.call(c, lapply(len, combn, x = x, simplify = FALSE))
}

#' Subsets of J
#'
#' @examples
#' J <- c(2, 5, 3) # three domains with 2, 5, and 3 strata respectively.
#' subsets_of_J(J)
#'
subsets_of_J <- function(J) {
  indices_i <- lseq_len(J)
  subsets_i <- lapply(indices_i, subsets)
  subsets_i <- subsets_i[lengths(subsets_i) != 0] # remove NULL (occurs if 1s are in J)
  if (is_empty(subsets_i)) {
    return(NULL)
  }

  subsets1 <- do.call(c, subsets_i)
  J <- J[J != 1] # remove domains with 1 stratum
  subsets2 <- if (length(J) == 1L) {
    NULL
  } else {
    crossjoin_J <- subsets(seq_along(J), 2:length(J)) # combination of elements of J.
    subsets2 <- lapply(crossjoin_J, function(x) {
      cj <- do.call(expand.grid, subsets_i[x])
      colnames(cj) <- NULL
      split(cj, seq(nrow(cj)))
    })
    subsets2 <- do.call(c, subsets2)
    setNames(lapply(subsets2, unlist), NULL)
  }

  all_subsets <- c(list(NULL), subsets1, subsets2)

  # verification
  combn_in_i <- 2^J - 2
  if (length(combn_in_i) == 1L) {
    combn_in_i <- c(0, combn_in_i) # since combn() works differently for scalar x.
  }
  no_all_subsets <- sum(
    1,
    sapply(seq_along(J), function(i) sum(combn(combn_in_i, i, prod)))
  )
  if (length(all_subsets) != no_all_subsets) {
    stop("Zle obliczone podzbiory")
  }
  all_subsets
}

#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @param n total sample size
#' @param J vector with domains labels
#' @param N population sizes
#' @param S population standard deviations of surveyed variable
#' @param total totals of surveyed variable in domains
#' @param kappa priority weights for domains
#' @param active user supplied take-max strata set, indices of J.
#' @param short short output
#'
#' @examples
#' J <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(J, N, S) - 1)
#'
#' (x <- fixprec(n, J, N, S, total, kappa))
#' x$n_ih # 162.47371 42.55264 142.98313 179.99052
#' x$check_cnstr
#' x$T_eigenval # 0.8476689
#'
#' (x_opt <- fixprec(n, J, N, S, total, kappa, active = 1))
#' x_opt$n_ih # 140 106.8521 124.4665 156.6814
#' x_opt$check_cnstr
#' x_opt$T_eigenval # 40.50748
#'
fixprec <- function(n, J, N, S, total, kappa = NULL, active = NULL, details = TRUE) {
  if (n >= nmax(J, N, S)) {
    stop("Total sample size n is too large")
  }
  if (is.null(kappa)) {
    kappa <- rep(1 / length(J), length(J))
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa
  J0 <- J
  N0 <- N
  S0 <- S
  n0 <- n

  J <- domain_indicators(J)
  if (!is.null(active)) {
    n <- n - sum(N[active])
    J <- J[-active]
    N <- N[-active]
    S <- S[-active]
  }

  a <- as.matrix(tapply(N * S, J, sum) / rho)
  c_ <- tapply(N * S^2, J, sum) / rho2 # - b (b = 0 if M = N)
  D.matrix <- (a %*% t(a)) / n - diag(c_, nrow = length(c_))

  eigenout <- eigen(D.matrix, symmetric = TRUE)
  T_eigenval <- eigenout$values[1] # largest eigenvalue
  v <- eigenout$vectors[, 1] # corresponding eigenvector
  if (T_eigenval <= 0) {
    stop("Largest eigenvalue is not strictly positive - solution does not exist!")
  }
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s_i <- n * v / as.numeric(t(a) %*% v)
  A <- (N * S) / rep(rho, times = table(J)) # brakets due to finite-prec arithmetic!
  x <- rep(s_i, table(J)) * A
  if (!is.null(active)) {
    N <- N0
    N[-active] <- x
    x <- N
  }

  if (details) {
    check_cnstr <- check_kkt(x, J0, N0, S0, total, kappa, n0, active, s_i, details = TRUE)
    A <- (N0 * S0) / rep(rho, times = J0)
    list(
      T_eigenval = T_eigenval, Ti = kappa * T_eigenval, check_cnstr = check_cnstr,
      D.matrix = D.matrix, eigen = eigenout, N_over_A = N0 / A, s_i = s_i, n_ih = x
    )
  } else {
    list(T_eigenval = T_eigenval, n_ih = x)
  }
}

#' Recursive FIXPREC algorithm (for 2 domains only)
#'
#' @examples
#' J <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(J, N, S) - 1)
#'
#' (x <- rfixprec_2d(n, J, N, S, total, kappa))
#' x$n_ih # 140 106.8521 124.4665 156.6814
#' x$s_i # T = 40.5074790
#'
rfixprec_2d <- function(n, J, N, S, total, kappa = NULL, domain = 1L) {
  W_domains <- lseq_len(J) # list with strata global indices in all domains
  W <- W_domains[[domain]] # strata global indices in recursive domain
  domain2 <- ifelse(domain == 1, 2L, 1L) # zalozenie: sa tylko dwie domeny
  W2 <- W_domains[[domain2]]
  check_for_active2 <- W2

  active <- NULL
  active_temp <- NULL
  iter <- 1
  s_i <- NULL
  repeat{
    # Domena wyrozniona (domyslnie 1)
    check_for_active <- W # potentially active strata in a given domain
    active_temp <- active
    repeat {
      x <- fixprec(n, J, N, S, total, kappa, active_temp)
      s_i <- rbind(s_i, c(iter, x$check_cnstr$T_opt, x$s_i)) # for debug only
      x <- x$n_ih
      active1 <- which(x[check_for_active] >= N[check_for_active])
      if (length(active1) == 0L) {
        break
      } else {
        active_temp <- c(active_temp, check_for_active[active1])
        if (setequal(active1, W)) {
          stop("Wszystkie wiezy na N aktywne w wyroznionej domenie") # wstepne symulacje pokazuje, ze niemozliwe
        }
        check_for_active <- check_for_active[-active1]
      }
    }
    # Pozostala domena
    active2 <- which(x[check_for_active2] >= N[check_for_active2])
    if (length(active2) == 0L) {
      break
    } else {
      active <- c(active, check_for_active2[active2])
      if (setequal(active2, W2)) {
        stop("Wszystkie wiezy na N aktywne w niewyroznionej domenie") # wstepne symulacje pokazuje, ze niemozliwe
      }
      check_for_active2 <- check_for_active2[-active2]
      iter <- iter + 1
    }
  }

  colnames(s_i) <- c("iter", "T", paste0("s_", 1:(ncol(s_i) - 2)))
  list(n_ih = x, s_i = s_i, active = active_temp)
}

#' Check KKT conditions for FIXPREC problem.
#'
#' @return If `details` is `FALSE`, it returns optimal value of T
#'  (if all constraints are satisfied) or `NA` (if any of the constraints is not satisfied).
#'  If `details` is `TRUE`, a detailed constraints check is returned.
#'
#' @examples
#' J <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(J, N, S) - 1)
#'
#' (x <- fixprec(n, J, N, S, total, kappa))
#' check_kkt(x$n_ih, J, N, S, total, kappa, n, active = NULL, s = x$s_i, details = FALSE)
#' check_kkt(x$n_ih, J, N, S, total, kappa, n, active = NULL, s = x$s_i, details = TRUE)
#'
check_kkt <- function(n_opt, J, N, S, total, kappa, n, active, s, tol = 10^-9, details = FALSE) {
  Ti_notscaled <- tapply(N^2 / n_opt * S^2 - N * S^2, domain_indicators(J), sum)
  Ti <- Ti_notscaled / total^2
  T_opt <- sum(Ti_notscaled) / sum(total^2 * kappa)

  is_Ti_eq_kappa_T <- Ti - kappa * T_opt < tol
  is_sum_nopt_eq_n <- sum(n_opt) - n < tol
  n_ih_leq_N <- n_opt <= N

  # check if mu_ih is non-negative for active constraints
  is_mu_ih_nonneg <- if (is.null(active)) {
    TRUE
  } else {
    rho <- total * sqrt(kappa)
    A <- (N * S) / rep(rho, times = J)
    s <- rep(s, times = J)
    all((s >= N / A)[active])
  }

  if (details) {
    list(
      T_opt = T_opt,
      Ti = Ti,
      is_Ti_eq_kappa_T = is_Ti_eq_kappa_T,
      is_sum_nopt_eq_n = is_sum_nopt_eq_n,
      n_ih_leq_N = n_ih_leq_N,
      is_mu_ih_nonneg = is_mu_ih_nonneg
    )
  } else {
    if (all(is_Ti_eq_kappa_T, is_sum_nopt_eq_n, n_ih_leq_N, is_mu_ih_nonneg)) {
      T_opt
    } else {
      NA_real_
    }
  }
}

#' Maximum allowed total sample size
#'
#' @details
#' See (16) from JW 2019. It is usually less than sum(N).
#'
#' @examples
#' J <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' sum(N)
#' nmax(J, N, S)
#'
nmax <- function(J, N, S) {
  J <- domain_indicators(J)
  n_bound <- sum(tapply(N * S, J, sum)^2 / tapply(N * S^2, J, sum))
  n_bound_floor <- floor(n_bound)
  ifelse(n_bound == n_bound_floor, n_bound - 1, n_bound_floor)
}

domain_indicators <- function(J) {
  rep(seq_along(J), times = J)
}

#' Get globally unique strata indices
#'
#' @param J vector with numbers of strata in domains
#'
#' @return A `list` with globally unique strata indices in domains.
#' @examples
#' J <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' lseq_len(J)
#'
lseq_len <- function(J) {
  mapply(
    function(to, len) seq.int(to = to, length.out = len, by = 1L),
    cumsum(J), J,
    SIMPLIFY = FALSE
  )
}

is_non_empty <- function(x) {
  length(x) > 0L
}

is_empty <- function(x) {
  length(x) == 0L
}

# Przyklad 0 ----

J <- c(1, 2)
N <- c(14, 11, 90)
S <- sqrt(c(4950, 50, 2))
total <- c(2, 3)
kappa <- c(0.4, 0.6)
(n <- nmax(J, N, S) - 1)

(x <- fixprec(n, J, N, S, total, kappa, active = 2))

# Przyklad 1 (1>,3>) ----

# J = 2 2
# Przekroczone 1,3,
# -1 => alokacja optymalna
# -3 => alokacja not feasible
# -1, -3 => alokacja feasible, nie optymalna

J <- c(2, 2) # two domains with 2 strata each.
N <- c(140, 110, 135, 190)
S <- sqrt(c(180, 20, 5, 4))
total <- c(2, 3)
kappa <- c(0.4, 0.6)
(n <- nmax(J, N, S) - 1)

ampl_alloc <- ampl_fixprec(n, J = c(2, 2), N, S, total, kappa, model = model)
ampl_alloc$Topt
ampl_alloc$n_ih

(x <- fixprec(n, J, N, S, total, kappa))
# 162.47371  42.55264 142.98313 179.99052 (1>, 3>)
# T = 0.8476689

# Usuniecie 1, naprawia 3.
(x_1 <- fixprec(n, J, N, S, total, kappa, active = 1, details = TRUE))
x_1$T_eigenval
ampl_alloc$Topt
ampl_alloc$n_ih
x_1$n_ih
# 140.0000 106.8521 124.4665 156.6814 (OK)
# T = 40.50748 (optimal)

(x_2 <- fixprec(n, J, N, S, total, kappa, active = 2, details = TRUE))
# T = 43.54576 (negative mu)

(x_3 <- fixprec(n, J, N, S, total, kappa, active = 3, details = TRUE))
# T = 1.511168 (1>)

(x_4 <- fixprec(n, J, N, S, total, kappa, active = 4, details = TRUE))
# T = 1.892899 (1>)

(x_13 <- fixprec(n, J, N, S, total, kappa, active = c(1, 3), details = TRUE))
# T = 42.08127 (negative mu)

(x_14 <- fixprec(n, J, N, S, total, kappa, active = c(1, 4), details = TRUE))
# T = 57.58687 (negative mu)

(x_23 <- fixprec(n, J, N, S, total, kappa, active = c(2, 3), details = TRUE))
# T = 45.72894 (negative mu)

(x_24 <- fixprec(n, J, N, S, total, kappa, active = c(2, 4), details = TRUE))
# T = 65.50605 (negative mu)

# Przyklad 2 ----
# Przekroczone 1,4,
# -1 => alokacja feasible
# -3 => alokacja feasible (optymalna)

# Przyklad 3 ----
# Przekroczone 1, 2
# -1 => alokacja optymalna
# -3 => alokacja feasible

# Przyklad 4 ----

H1 <- 3
H2 <- 4
J <- c(H1, H2)
ind1 <- sample(nrow(stratallo::pop507), H1)
ind2 <- sample(nrow(stratallo::pop507), H2)
N <- c(stratallo::pop507[ind1, "N"], stratallo::pop507[ind2, "N"])
S <- c(stratallo::pop507[ind1, "S"], stratallo::pop507[ind2, "S"])
total <- c(2, 2)
kappa <- c(0.3, 0.3)
(nm <- nmax(J, N, S) - 1)
all_subsets <- subsets_of_J(J)

delt <- logical(nm)
for (n in 1:nm) {
  print(n)
  # recursive FIXPREC
  x_rec <- rfixprec_2d(n, J, N, S, total, kappa, domain = 1L)
  s <- as.vector(tail(x_rec$s_i[, -(1:2), drop = FALSE], 1))
  Topt <- check_kkt(x_rec$n_ih, J, N, S, total, kappa, n, active = x_rec$active, s = s, tol = 0.9)

  # all subsets
  Topt_as <- sapply(all_subsets, function(i) {
    if (n > sum(N[i])) {
      x <- fixprec(n, J, N, S, total, kappa, active = i)
      check_kkt(x$n_ih, J, N, S, total, kappa, n, active = i, s = x$s_i, tol = 0.9)
    } else {
      NA_real_
    }
  })
  Topt_as <- Topt_as[!is.na(Topt_as)]
  if (length(Topt_as) == 0 || length(Topt_as) >= 2) {
    print(n)
  }
  delt[n] <- min(Topt_as) == Topt
}
table(delt, useNA = "ifany")

# 2 duze domeny (monotonicznosc s) ----

H1 <- 10
H2 <- 10

J <- c(H1, H2) # two domains with 2 strata each.
N <- c(sample(1000, H1), sample(500, H2))
S <- sqrt(c(sample(3000, H1), sample(1000, H2)))
total <- c(20, 500)
kappa <- c(0.3, 0.7)
(nm <- nmax(J, N, S) - 1)

al <- rfixprec_2d(nm, J, N, S, total, kappa)

# monotonicznosc s_i
result <- sapply(1:nm, function(n) {
  x <- rfixprec_2d(n, J, N, S, total, kappa)
  if (nrow(x$s_i) >= 2) {
    tapply(1:nrow(x$s_i), x$s_i[, "iter"], function(i) {
      all(diff(x$s_i[i, "s_1", drop = TRUE]) > 0) &&
        all(diff(x$s_i[i, "s_2", drop = TRUE]) < 0)
    })
  } else {
    NA
  }
})
table(unlist(result), useNA = "ifany")

# testy grid (parallel) ----

x <- 2
pop <- expand.grid(20:(20 + x), 30:(30 + x), 20:(20 + x), 40:(40 + x), 2:(2 + x), 12:(12 + x), 21:(21 + x), 30:(30 + x))
J <- c(2, 2)
total <- c(10, 10)
kappa <- c(0.5, 0.5)

library(parallel)
result <- mclapply(
  1:100, #  1:nrow(pop),
  function(i) {
    # print(i / nrow(pop))
    N <- unlist(pop[i, 1:4])
    S <- unlist(pop[i, 5:8])
    nm <- nmax(J, N, S) - 1
    all_subsets <- subsets_of_J(J)

    delt <- logical(nm)
    for (n in 1:nm) {
      # print(n)
      # recursive FIXPREC
      x_rec <- rfixprec_2d(n, J, N, S, total, kappa, domain = 1L)
      s <- as.vector(x_rec$s_i[nrow(x_rec$s_i), -(1:2), drop = FALSE])
      Topt <- check_kkt(x_rec$n_ih, J, N, S, total, kappa, n, active = x_rec$active, s = s, tol = 0.9)

      # all subsets
      Topt_as <- sapply(all_subsets, function(i) {
        if (n > sum(N[i])) {
          x <- fixprec(n, J, N, S, total, kappa, active = i)
          check_kkt(x$n_ih, J, N, S, total, kappa, n, active = i, s = x$s_i, tol = 0.9)
        } else {
          NA_real_
        }
      })
      Topt_as <- Topt_as[!is.na(Topt_as)]
      if (length(Topt_as) == 0 || length(Topt_as) >= 2) {
        stop("Topt not uniqe or missing")
      }
      delt[n] <- min(Topt_as) == Topt
    }
    all(delt)
  },
  mc.cores = 10
)
table(unlist(result), useNA = "ifany")
