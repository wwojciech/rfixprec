#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @param n total sample size
#' @param H_count vector with domains labels
#' @param N population sizes
#' @param S population standard deviations of surveyed variable
#' @param total totals of surveyed variable in domains
#' @param kappa priority weights for domains
#' @param active user supplied take-max strata set
#' @param short short output
#'
#' @examples
#' H_count <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(H_count, N, S) - 1)
#'
#' (x <- fixprecact0(n, H_count, N, S, total, kappa, details = TRUE))
#' x$n_dh # 162.47371 42.55264 142.98313 179.99052
#' x$T_eigenval # 0.8476689
#'
#' (x_opt <- fixprecact0(n, H_count, N, S, total, kappa, active = 1, details = TRUE))
#' x_opt$n_dh # 140 106.8521 124.4665 156.6814
#' x_opt$T_eigenval # 40.50748
#'
fixprecact0 <- function(n, H_count, N, S, total, kappa = NULL, active = NULL, details = FALSE) {
  if (n >= nmax(H_count, N, S)) {
    # stop("total sample size n is too large")
    # warning("total sample size n is too large")
  }
  if (is.null(kappa)) {
    kappa <- rep(1 / length(H_count), length(H_count))
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa
  N0 <- N
  S0 <- S
  n0 <- n

  H_dind <- domain_indicators(H_count)
  if (nonempty(active)) {
    n <- n - sum(N[active])
    H_dind <- H_dind[-active]
    N <- N[-active]
    S <- S[-active]
  }

  a.vec <- as.matrix(tapply(N * S, H_dind, sum) / rho)
  c.vec <- tapply(N * S^2, H_dind, sum) / rho2 # - b (b = 0 if M = N)
  D.matrix <- (a.vec %*% t(a.vec)) / n - diag(c.vec, nrow = length(c.vec))

  eigenout <- eigen(D.matrix, symmetric = TRUE)
  T_eigenval <- eigenout$values[1] # largest eigenvalue
  v <- eigenout$vectors[, 1] # corresponding eigenvector
  if (T_eigenval <= 0) {
    # stop("largest eigenvalue is not strictly positive - solution does not exist!")
    # warning("largest eigenvalue is not strictly positive")
    # print("largest eigenvalue is not strictly positive")
  }
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s_d <- n * v / as.numeric(t(a.vec) %*% v)
  A <- (N * S) / rep(rho, table(H_dind)) # brakets due to finite-prec arithmetic!
  x <- rep(s_d, table(H_dind)) * A
  if (nonempty(active)) {
    N <- N0
    N[-active] <- x
    x <- N
  }

  if (details) {
    kkt <- check_kkt(x, H_count, N0, S0, total, kappa, n0, active, s_d, details = TRUE)
    A <- (N0 * S0) / rep(rho, times = H_count)
    list(
      D.matrix = D.matrix, eigen = eigenout, T_eigenval = T_eigenval,
      N_over_A = N0 / A, s_d = s_d, kkt = kkt, n_dh = x
    )
  } else {
    x
  }
}

#' Recursive FIXPREC algorithm (for 2 domains only)
#'
#' @examples
#' H_count <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(H_count, N, S) - 1)
#'
#' (x <- rfixprec_2d(n, H_count, N, S, total, kappa))
#' x$n_dh # 140 106.8521 124.4665 156.6814
#' x$s_d # T = 40.5074790
#'
rfixprec_2d <- function(n, H_count, N, S, total, kappa = NULL, domain = 1L) {
  W_domains <- lseq_len(H_count) # list with strata global indices in all domains
  W <- W_domains[[domain]] # strata global indices in recursive domain
  domain2 <- ifelse(domain == 1, 2L, 1L) # zalozenie: sa tylko dwie domeny
  W2 <- W_domains[[domain2]]
  check_for_active2 <- W2

  active <- NULL
  active_temp <- NULL
  iter <- 1
  s_d <- NULL
  repeat{
    # Domena wyrozniona (domyslnie 1)
    check_for_active <- W # potentially active strata in a given domain
    active_temp <- active
    repeat {
      x <- fixprecact0(n, H_count, N, S, total, kappa, active_temp, details = TRUE)
      s_d <- rbind(s_d, c(iter, x$kkt$T_opt, x$s_d)) # for debug only
      x <- x$n_dh
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

  colnames(s_d) <- c("iter", "T", paste0("s_", 1:(ncol(s_d) - 2)))
  list(n_dh = x, s_d = s_d, active = active_temp)
}

#' Recursive FIXPREC algorithm
#'
#' @examples
#' H_count <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- nmax(H_count, N, S) - 1)
#'
#' rfixprec_iter0(n, H_count, N, S, total, kappa)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
#'
rfixprec_iter0 <- function(n, H_count, N, S, total, kappa = NULL, recursive_domain = 1L) {
  stopifnot(length(recursive_domain) == 1L)
  stopifnot(recursive_domain <= length(H_count))

  W <- lseq_len(H_count) # list with strata global indices in all domains
  W_rec <- W[[recursive_domain]] # recursive domain
  W_nrec <- W[-recursive_domain] # non recursive domain
  W_nrec_active <- vector("list", length(W_nrec))
  I_nrec <- seq_along(W_nrec) # indices of non recursive domains

  repeat {
    # 1. Allocate in recursive_domain.
    W_active <- unlist(W_nrec_active)
    repeat {
      x <- fixprecact0(n, H_count, N, S, total, kappa, W_active)
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
    for (d in I_nrec) {
      W_d <- W_nrec[[d]]
      violated <- which(x[W_d] > N[W_d]) # W_d powinno zostac zmniejszone o active temp
      if (nonempty(violated)) {
        if (d >= 2 && nonempty(W_nrec_active[1:(d - 1)])) { # wyczysc temp active z poprzednich domen
          W_nrec_active[1:(d - 1)] <- list(NULL)
        }
        W_nrec_active[[d]] <- c(W_nrec_active[[d]], W_d[violated])
        break
      }
    }

    if (empty(violated)) {
      break
    }
  }
  x
}

#' Check KKT conditions for FIXPREC problem.
#'
#' @return If `details` is `FALSE`, it returns optimal value of T
#'  (if all constraints are satisfied) or `NA` (if any of the constraints is not satisfied).
#'  If `details` is `TRUE`, a detailed constraints check is returned.
#'
#' @examples
#' H_count <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' (n <- nmax(H_count, N, S) - 1)
#'
#' (x <- fixprecact0(n, H_count, N, S, total, kappa))
#' check_kkt(x, H_count, N, S, total, kappa, n)
#' check_kkt(x, H_count, N, S, total, kappa, n, details = TRUE)
#'
check_kkt <- function(n_opt, H_count, N, S, total, kappa, n, active = NULL, s = NULL, tol = 10^-9, details = FALSE) {
  strata_ind <- seq_along(n_opt)

  Td_notscaled <- tapply(N^2 / n_opt * S^2 - N * S^2, domain_indicators(H_count), sum)
  Td <- Td_notscaled / total^2
  T_opt <- sum(Td_notscaled) / sum(total^2 * kappa)
  is_Td_eq_kappa_T <- Td - kappa * T_opt < tol
  is_sum_nopt_eq_n <- sum(n_opt) - n < tol
  is_nopt_leq_N <- n_opt <= N
  names(is_nopt_leq_N) <- strata_ind

  # check if mu_ih is non-negative for active constraints
  if (empty(active)) {
    active <- which(n_opt == N)
  } else {
    stopifnot(setequal(active, which(n_opt == N))) # might be time consuming
  }

  is_mu_dh_nonneg <- if (empty(s)) {
    n_opt[active] >= N[active]
  } else {
    rho <- rep(total * sqrt(kappa), times = H_count)
    s <- rep(s, times = H_count)
    s[active] >= rho[active] / S[active]
  }
  names(is_mu_dh_nonneg) <- active

  if (details) {
    list(
      T_opt = T_opt,
      Td = Td,
      is_Td_eq_kappa_T = is_Td_eq_kappa_T,
      is_sum_nopt_eq_n = is_sum_nopt_eq_n,
      is_nopt_leq_N = is_nopt_leq_N,
      is_mu_dh_nonneg = is_mu_dh_nonneg,
      active = active,
      active_len = length(active)
    )
  } else {
    if (all(is_Td_eq_kappa_T, is_sum_nopt_eq_n, is_nopt_leq_N, is_mu_dh_nonneg)) {
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
#' H_count <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' sum(N)
#' nmax(H_count, N, S)
#'
nmax <- function(H_count, N, S) {
  H_dind <- domain_indicators(H_count)
  n_bound <- sum(tapply(N * S, H_dind, sum)^2 / tapply(N * S^2, H_dind, sum))
  n_bound_floor <- floor(n_bound)
  ifelse(n_bound == n_bound_floor, n_bound - 1, n_bound_floor)
}

domain_indicators <- function(H_count) {
  rep(seq_along(H_count), times = H_count)
}

#' Convert strata sizes to globally unique strata indices
#'
#' @param H_counts vector of strata sizes in domains
#'
#' @return A `list` with globally unique strata indices in domains.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' H_s2i(H_counts)
#'
H_s2i <- function(H_counts) {
  mapply(
    function(to, len) seq.int(to = to, length.out = len, by = 1L),
    cumsum(H_counts), H_counts,
    SIMPLIFY = FALSE
  )
}

#' Get globally unique strata indices
#'
#' @param H_count vector with numbers of strata in domains
#'
#' @return A `list` with globally unique strata indices in domains.
#' @examples
#' H_count <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' lseq_len(H_count)
#'
lseq_len <- function(H_count) {
  mapply(
    function(to, len) seq.int(to = to, length.out = len, by = 1L),
    cumsum(H_count), H_count,
    SIMPLIFY = FALSE
  )
}

empty <- function(x) {
  length(x) == 0L
}

nonempty <- function(x) {
  length(x) > 0L
}
