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
#' check_kkt(x$n_ih, J, N, S, total, kappa, n)
#' check_kkt(x$n_ih, J, N, S, total, kappa, n, details = TRUE)
#'
check_kkt <- function(n_opt, J, N, S, total, kappa, n, active = NULL, s, tol = 10^-9, details = FALSE) {
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

empty <- function(x) {
  length(x) == 0L
}

non_empty <- function(x) {
  length(x) > 0L
}

#' Recursive FIXPREC algorithm
#'
#' @examples
#' J <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- nmax(J, N, S) - 1)
#'
#' x <- rfixprec(n, J, N, S, total, kappa)
#' x # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
#' x$T_eigenval # 67.90425
#'
rfixprec <- function(n, J, N, S, total, kappa = NULL, recursive_domain = 1L) {
  stopifnot(length(recursive_domain) == 1L)
  stopifnot(recursive_domain <= length(J))

  W <- lseq_len(J) # list with strata global indices in all domains
  W_rec <- W[[recursive_domain]] # recursive domain
  W_nrec <- W[-recursive_domain] # non recursive domain
  W_nrec_active <- vector("list", length(W_nrec))
  I_nrec <- seq_along(W_nrec) # indices of non recursive domains

  repeat {
    # 1. Allocate in recursive_domain.
    W_active <- unlist(W_nrec_active)
    repeat {
      x <- fixprec(n, J, N, S, total, kappa, W_active)$n_ih
      violated <- which(x[W_rec] > N[W_rec])
      if (empty(violated)) {
        break
      } else {
        W_active <- c(W_active, W_rec[violated])
        # if (all(W_rec %in% W_active)) {
        #   stop("Wszystkie wiezy aktywne w wyroznionej domenie") # wstepne symulacje pokazuje, ze niemozliwe
        # }
      }
    }

    # 2. Check for violations in remaining (non recursive) domains.
    for (i in I_nrec) {
      W_i <- W_nrec[[i]]
      violated <- which(x[W_i] > N[W_i]) # W_i powinno zostac zmniejszone o active temp
      if (non_empty(violated)) {
        if (i >= 2 && empty(W_nrec_active[[i]])) { # wyczysc poprzednie active
          W_nrec_active[1:(i - 1)] <- list(NULL)
        }
        W_nrec_active[[i]] <- c(W_nrec_active[[i]], W_i[violated])
        break
      }
    }

    if (empty(violated)) {
      break
    }
  }
  x
}
