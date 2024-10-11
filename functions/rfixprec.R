# ALGORITHMS ----

#' FIXPREC algorithm
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @param n total sample size
#' @param H_ss vector of strata sizes in domains
#' @param N population sizes
#' @param S population standard deviations of surveyed variable
#' @param total totals of surveyed variable in domains
#' @param kappa priority weights for domains
#'
#' @examples
#' H_ss <- c(2, 2, 3) # strata sizes, three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- nmax(H_ss, N, S) - 1)
#'
#' fixprec(n, H_ss, N, S, total, kappa)
#' # 162.46200  42.54957 143.14020 180.18824 203.07887  20.59596  75.98516
#'
fixprec <- function(n, H_ss, N, S, total, kappa = NULL) {
  if (is.null(kappa)) {
    kappa <- rep(1 / length(H_ss), length(H_ss))
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa

  # stratum-domain indicators,
  # e.g. for 2 domains with 4 and 2 strata: 1 1 1 1 2 2
  H_di <- H_domain_indicators(H_ss)

  a.vec <- as.matrix(tapply(N * S, H_di, sum) / rho)
  c.vec <- tapply(N * S^2, H_di, sum) / rho2 # - b (b = 0 if M = N)
  D.matrix <- (a.vec %*% t(a.vec)) / n - diag(c.vec, nrow = length(c.vec))

  eigen_decomp <- eigen(D.matrix, symmetric = TRUE)
  lambda <- eigen_decomp$values[1] # largest eigenvalue
  v <- eigen_decomp$vectors[, 1] # corresponding eigenvector
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s.vec <- n * v / as.numeric(t(a.vec) %*% v)
  A <- (N * S) / rep(rho, table(H_di)) # brackets due to finite-prec arithmetic!
  rep(s.vec, table(H_di)) * A
}

#' FIXPREC algorithm (version that accepts active strata and with debug)
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @inheritParams fixprec
#' @param active global indices of forced take-max strata.
#'   They are effectively removed from the population. It can happen that
#'   whole domain will be removed - then the dim of D matrix is decreased
#'   accordingly.
#' @param details detailed debug output
#'
#' @examples
#' H_ss <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(H_ss, N, S) - 1)
#'
#' (x <- fixprec_act(n, H_ss, N, S, total, kappa, details = TRUE))
#' x$n_dh # 162.47371 42.55264 142.98313 179.99052
#' x$lambda # 0.8476689
#'
#' (x_opt <- fixprec_act(n, H_ss, N, S, total, kappa, active = 1, details = TRUE))
#' x_opt$n_dh # 140 106.8521 124.4665 156.6814
#' x_opt$lambda # 40.50748
#'
fixprec_act <- function(n, H_ss, N, S, total, kappa = NULL, active = NULL, details = FALSE) {
  if (n >= nmax(H_ss, N, S)) {
    # warning("n > nmax, hence optimal T can be negative")
  }
  if (is.null(kappa)) {
    kappa <- rep(1 / length(H_ss), length(H_ss))
  }

  rho <- total * sqrt(kappa)
  rho0 <- rho
  rho2 <- total^2 * kappa
  total0 <- total
  kappa0 <- kappa
  H_di0 <- H_ss
  N0 <- N
  S0 <- S
  n0 <- n

  H_di <- H_domain_indicators(H_ss)
  if (nonempty(active)) {
    n <- n - sum(N[active])
    H_di <- H_di[-active]
    N <- N[-active]
    S <- S[-active]
    # check if there is any domain with all inequality constraints active
    H <- H_s2i(H_di0) # list with strata global indices, list elements - domains
    D_active <- which(sapply(H, function(x) all(x %in% active)))
    if (nonempty(D_active)) {
      paste("Blokada domeny:", paste(1:2, collapse = ","))
      rho <- rho[-D_active]
      rho2 <- rho2[-D_active]
    }
  }

  a.vec <- as.matrix(tapply(N * S, H_di, sum) / rho)
  c.vec <- tapply(N * S^2, H_di, sum) / rho2 # - b (b = 0 if M = N)
  D.matrix <- (a.vec %*% t(a.vec)) / n - diag(c.vec, nrow = length(c.vec))

  eigen_decomp <- eigen(D.matrix, symmetric = TRUE)
  lambda <- eigen_decomp$values[1] # largest eigenvalue
  v <- eigen_decomp$vectors[, 1] # corresponding eigenvector
  if (lambda < 0) {
    # warning("largest eigenvalue is negative")
    # print("lambda < 0")
  }
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s.vec <- n * v / as.numeric(t(a.vec) %*% v)
  A <- (N * S) / rep(rho, table(H_di)) # brackets due to finite-prec arithmetic!
  x <- rep(s.vec, table(H_di)) * A
  if (nonempty(active)) {
    N <- N0
    N[-active] <- x
    x <- N
  }
  x
  if (details) {
    kkt <- check_kkt(
      x, H_di0, N0, S0, total0, kappa0, n0,
      active = active, s = s.vec, details = TRUE
    )
    A <- (N0 * S0) / rep(rho0, times = H_di0)
    list(
      D.matrix = D.matrix, eigen_decomp = eigen_decomp, lambda = lambda,
      N_over_A = N0 / A, s = s.vec, kkt = kkt, n_dh = x
    )
  } else {
    x
  }
}

#' Recursive FIXPREC algorithm
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains. Allocation preserves strata sizes.
#'
#' @inheritParams fixprec
#' @param J vector of domain indices. Specifies domains for which the allocated
#'   samples should preserve strata sizes. For domains other than those specified
#'   in `J`, allocations may exceed strata sizes.
#'
#' @examples
#' H_ss <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- nmax(H_ss, N, S) - 1)
#'
#' rfixprec(n, H_ss, N, S, total, kappa)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
#'
rfixprec <- function(n, H_ss, N, S, total, kappa = NULL, J = NULL) {
  H <- H_s2i(H_ss) # list with strata global indices, list elements - domains
  if (is.null(J)) {
    J <- seq_along(H_ss)
  }

  N0 <- N
  j <- J[1]
  H_j <- H[[j]] # strata global indices for domain j
  j1 <- H_j[1] # global index of first stratum in domain j
  H_j_len <- length(H_j)

  h_tocheck <- H_j
  repeat {
    x <- if (length(J) == 1L) {
      fixprec(n, H_ss, N, S, total, kappa)
    } else {
      rfixprec(n, H_ss, N, S, total, kappa, J[-1])
    }

    h_violated <- which(x[h_tocheck] >= N[h_tocheck])
    if (nonempty(h_violated)) {
      n <- n - sum(N[h_tocheck[h_violated]])
      N <- N[-h_tocheck[h_violated]]
      S <- S[-h_tocheck[h_violated]]

      H_j <- H_j[-h_violated]
      H_ss[j] <- H_ss[j] - length(h_violated)
      if (H_ss[j] == 0L) { # whole domain j is (temporarily) active
        # print("BLOCK DOMAIN")
        H_ss <- H_ss[-j]
        total <- total[-j]
        kappa <- kappa[-j]
        J <- J - 1
        h_tocheck <- integer(0)
      } else {
        h_tocheck <- seq(from = j1, len = H_ss[j])
      }
    } else {
      break
    }
  }
  H[[j]] <- H_j
  N0[unlist(H)] <- x
  N0
}

#' Recursive FIXPREC algorithm (iterative implementation)
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains. Allocation preserves strata sizes.
#'
#' @inheritParams fixprec
#' @param ref_domain reference domain (denoted by j in the paper).
#'
#' @examples
#' H_ss <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- nmax(H_ss, N, S) - 1)
#'
#' rfixprec(n, H_ss, N, S, total, kappa)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
#'
rfixprec_iter <- function(n, H_ss, N, S, total, kappa = NULL, ref_domain = 1L) {
  stopifnot(length(ref_domain) == 1L)
  stopifnot(ref_domain <= length(H_ss))

  H <- H_s2i(H_ss) # list with strata global indices, list elements - domains
  H_ref <- H[[ref_domain]] # reference domain (denoted by j in the paper)
  H_nref <- H[-ref_domain] # remaining domains other than j
  H_nref_active <- vector("list", length(H_nref)) # active strata in remaining domains
  D_nref <- seq_along(H_nref) # indices of remaining domains

  repeat {
    # 1. Allocate in a chosen reference domain j
    H_active <- unlist(H_nref_active)
    repeat {
      x <- fixprec_act(n, H_ss, N, S, total, kappa, H_active)
      violated <- which(x[H_ref] > N[H_ref])
      if (empty(violated)) {
        break
      } else {
        H_active <- c(H_active, H_ref[violated])
        # if (all(H_ref %in% H_active)) { # jesli n < nmax to tak sie nie zdarzy (chyba)
        #   warning("all ineuqality constraints active in the reference domain")
        # }
      }
    }

    # 2. Check for violations in remaining (non recursive) domains.
    for (d in D_nref) {
      H_d <- H_nref[[d]]
      violated <- which(x[H_d] > N[H_d]) # W_d powinno zostac zmniejszone o active temp
      if (nonempty(violated)) {
        if (d >= 2 && nonempty(H_nref_active[1:(d - 1)])) {
          # wyczysc temp active z poprzednich domen
          H_nref_active[1:(d - 1)] <- list(NULL)
        }
        H_nref_active[[d]] <- c(H_nref_active[[d]], H_d[violated])
        break
      }
    }

    if (empty(violated)) {
      break
    }
  }
  x
}

# DIAGNOSTIC FUNCTIONS ----

#' Check KKT conditions for problem of equal-precision optimal allocation in
#' single-stage sampling with domains and strata in domains (allocations
#' do not exceed strata sizes)
#'
#' @inheritParams rfixprec
#' @param active global indices of forced take-max strata.
#' @param s values of functions s_d, for each domain d (optional).
#' @param tol number comparison tolerance.
#' @param details should detailed output be returned?
#'
#' @return If `details` is `FALSE`, it returns optimal value of T
#'  (if all constraints are satisfied) or `NA` (if any of the constraints is not satisfied).
#'  If `details` is `TRUE`, a detailed constraints check is returned.
#'
#' @examples
#' H_ss <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(H_ss, N, S) - 1)
#'
#' (x <- fixprec(n, H_ss, N, S, total, kappa))
#' check_kkt(x, H_ss, N, S, total, kappa, n)
#' check_kkt(x, H_ss, N, S, total, kappa, n, details = TRUE)
#'
check_kkt <- function(x, H_ss, N, S, total, kappa, n, J = seq_along(H_ss),
                      active = NULL, s = NULL, tol = 10^-9, details = FALSE) {
  Td_notscaled <- tapply(N^2 / x * S^2 - N * S^2, H_domain_indicators(H_ss), sum)
  Td <- Td_notscaled / total^2
  T_opt <- sum(Td_notscaled) / sum(total^2 * kappa)
  is_Td_eq_kappa_T <- Td - kappa * T_opt < tol
  is_sum_x_eq_n <- sum(x) - n < tol
  strata_check_N <- unlist(H_s2i(H_ss)[J])
  is_x_leq_N <- x[strata_check_N] <= N[strata_check_N]
  names(is_x_leq_N) <- strata_check_N

  # check if mu_ih is non-negative for active constraints
  if (empty(active)) {
    active <- which(x == N)
  } else {
    stopifnot(setequal(active, which(x == N))) # might be time consuming
  }

  is_mu_dh_nonneg <- if (empty(s)) {
    x[active] >= N[active]
  } else {
    rho <- rep(total * sqrt(kappa), times = H_ss)
    s <- rep(s, times = H_ss)
    s[active] >= rho[active] / S[active]
  }
  names(is_mu_dh_nonneg) <- active

  if (details) {
    list(
      T_opt = T_opt,
      Td = Td,
      is_Td_eq_kappa_T = is_Td_eq_kappa_T,
      is_sum_x_eq_n = is_sum_x_eq_n,
      is_x_leq_N = is_x_leq_N,
      is_mu_dh_nonneg = is_mu_dh_nonneg,
      active = active,
      active_len = length(active)
    )
  } else {
    if (all(is_Td_eq_kappa_T, is_sum_x_eq_n, is_x_leq_N, is_mu_dh_nonneg)) {
      T_opt
    } else {
      NA_real_
    }
  }
}

#' Maximum allowed total sample size so that D is positive
#'
#' @inheritParams fixprec
#' @details
#' See (16) from JW 2019. It is usually less than sum(N).
#'
#' @examples
#' H_ss <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' sum(N)
#' nmax(H_ss, N, S)
#'
nmax <- function(H_ss, N, S) {
  H_di <- H_domain_indicators(H_ss)
  n_bound <- sum(tapply(N * S, H_di, sum)^2 / tapply(N * S^2, H_di, sum))
  n_bound_floor <- floor(n_bound)
  ifelse(n_bound == n_bound_floor, n_bound - 1, n_bound_floor)
}

# HELPERS ----

empty <- function(x) {
  length(x) == 0L
}

nonempty <- function(x) {
  length(x) > 0L
}

H_domain_indicators <- function(H) {
  if (is.list(H)) {
    rep(seq_along(H), times = sapply(H, length))
  } else {
    rep(seq_along(H), times = H)
  }
}

#' Convert strata sizes to globally unique strata indices
#'
#' @param H_ss vector of strata sizes in domains
#'
#' @return A `list` with globally unique strata indices in domains.
#' @examples
#' H_ss <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' H_s2i(H_ss)
#'
H_s2i <- function(H_ss) {
  mapply(
    function(to, len) seq.int(to = to, length.out = len, by = 1L),
    cumsum(H_ss), H_ss,
    SIMPLIFY = FALSE
  )
}
