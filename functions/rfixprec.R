# ALGORITHMS ----

#' FIXPREC algorithm
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @note This is an internal function and should not be used directly by users.
#'   It is optimized for handling a large number of invocations, specifically
#'   the recursive calls from `rfixprec`, and as a result, parameter assertions
#'   are minimal.
#'
#' @param n (`int`)\cr total sample size. It must be that `0 < n <= sum(N)`.
#' @param H_counts (`integerish`)\cr vector of strata counts in domains.
#' @param N (`integerish`)\cr population sizes.
#' @param S (`numeric`)\cr population standard deviations of surveyed variable.
#' @param rho (`numeric`)\cr `rho = total * sqrt(kappa)`, where `total` is
#'   a vector of totals of surveyed variable in domains
#'   and `kappa` is a vector with priority weights for domains.
#' @param rho2 (`numeric`)\cr `rho2 = rho^2`, provided to improve
#'   the loss of precision due to finite-precision arithmetic issues.
#' @param details (`logical`)\cr detailed debug output.
#'
#' @examples
#' H_counts <- c(2, 2) # 2 domains with 2 strata each
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' rho <- total * sqrt(kappa)
#' rho2 <- total^2 * kappa
#' sum(N) # 575
#' n <- 528
#'
#' (x <- fixprec(n, H_counts, N, S, rho, details = TRUE))
#' x$x # 162.47371 42.55264 142.98313 179.99052
#' x$lambda # 0.8476689
#' x$s # 0.1094155 1.1006847
fixprec <- function(n, H_counts, N, S, rho, rho2, details = FALSE) {
  # H_counts as param instead of H_dind to avoid index misorder.
  H_dind <- H_cnt2dind(H_counts)

  # prepare D matrix
  a.vec <- tapply(N * S, H_dind, sum) / rho
  c.vec <- tapply(N * S^2, H_dind, sum) / rho2 # - b (b = 0 if M = N)
  D.matrix <- (a.vec %*% t(a.vec)) / n - diag(c.vec, nrow = length(c.vec))

  # Eigen
  eigen_decomp <- eigen(D.matrix, symmetric = TRUE)
  lambda <- eigen_decomp$values[1] # largest eigenvalue (can be < 0)
  v <- eigen_decomp$vectors[, 1] # eigenvector corresponding to lambda
  if (diffsigns(v)) {
    stop("eigenvector containts entries of a different sign")
  }

  # allocation
  s.vec <- n * v / sum(a.vec * v) # as.numeric(t(a.vec) %*% v)
  A <- (N * S) / rep(rho, H_counts) # brackets due to finite-prec arithmetic!
  x <- rep(s.vec, H_counts) * A

  # prepare return object
  if (details) {
    list(
      D.matrix = D.matrix, eigen_decomp = eigen_decomp, lambda = lambda,
      s = s.vec, x = x
    )
  } else {
    x
  }
}

#' FIXPRECACT algorithm
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains. It allows to force the set of take-max
#' strata whenever feasible (see `U`).
#'
#' @note This is an internal function and should not be used directly by users.
#'   It is optimized for handling a large number of invocations, specifically
#'   the recursive calls from `rfixprec`, and as a result, parameter assertions
#'   are minimal.
#'
#' @inheritParams fixprec
#' @param U global indices of forced take-max strata. They correspond to indices
#'   of `N`, `S` vectors. Strata in `U` are effectively removed from the
#'   population. It can happen that the entire domain will be removed, in which
#'   case the dimension of the D matrix is decreased accordingly.
#'   `U` must be such that `sum(N[U]) < n` for `sum(N) < n`, or
#'   `U = H` (where `H` denotes the set of all strata indices) for `sum(N) = n`.
#'   Otherwise, the function fails.
#'
#' @examples
#' H_counts <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' rho <- total * sqrt(kappa)
#' rho2 <- total^2 * kappa
#' (n <- round(nmax(H_counts, N, S)) - 1)
#'
#' x <- fixprecact(n, H_counts, N, S, rho, details = TRUE)
#' x$x # 162.47371 42.55264 142.98313 179.99052
#' x$x_Uc # 162.47371 42.55264 142.98313 179.99052
#' x$lambda # 0.8476689
#' x$s # 0.1094155 1.1006847
#'
#' x1 <- fixprecact(n, H_counts, N, S, rho, U = 1, details = TRUE)
#' x1$x # 140 106.8521 124.4665 156.6814
#' x1$x_Uc # 106.8521 124.4665 156.6814
#' x1$lambda # 40.50748
#' x1$s # 0.2747486 0.9581436
#'
fixprecact <- function(n, H_counts, N, S, rho, rho2, U = NULL, details = FALSE) {
  NU <- N[U]
  sumN_U <- sum(NU)
  if (!(sumN_U < n || length(NU) == length(N))) {
    stop("U must be such that sum(N[U]) < n or length(N[U]) == length(N)")
  }
  # if (!(sumN_U < n || sumN_U == n && sum(N) == n)) {
  #   stop("U must be such that sum(N[U]) < n or (sum(N[U]) = n && sum(N) == n)")
  # }

  # reduce the model
  if (nonempty(U)) {
    S <- S[-U]
    if (empty(S)) {
      return(N)
    }
    x <- N
    n <- n - sumN_U
    N <- N[-U]
    H_dind <- H_cnt2dind(H_counts)
    H_dind_Uc <- H_dind[-U]
    H_counts <- tabulate(H_dind_Uc, nbins = length(H_counts))
    D_blocked <- which(H_counts == 0)
    if (nonempty(D_blocked)) { # remove blocked domains, if any
      rho <- rho[-D_blocked]
      rho2 <- rho2[-D_blocked]
      H_counts <- H_counts[-D_blocked]
    }
  }

  # fixprec allocation for reduced model
  fixprec_out <- fixprec(
    n = n, H_counts = H_counts, N = N, S = S, rho = rho, rho2 = rho2,
    details = details
  )

  # prepare return object depending on the details flag
  x_Uc <- if (details) {
    fixprec_out$x
  } else {
    fixprec_out
  }

  if (nonempty(U)) {
    x[-U] <- x_Uc
  } else {
    x <- x_Uc
  }

  if (details) {
    names(fixprec_out)[names(fixprec_out) == "x"] <- "x_Uc"
    c(fixprec_out, list(x = x))
  } else {
    x
  }
}

#' Recursive FIXPREC algorithm
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains. Allocation preserves strata sizes.
#'
#' @note This is an internal function and should not be used directly by users.
#'   It is optimized for handling a large number of recursive calls,
#'   and as a result, parameter assertions are minimal.
#'
#' @inheritParams fixprecact
#' @param J vector of domain indices. Specifies domains for which the allocated
#'   samples should preserve strata sizes. For domains other than those specified
#'   in `J`, allocations may exceed strata sizes.
#'   If `J` is `NULL`, then `J` is treated as it would have contain all domains.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' rho <- total * sqrt(kappa)
#' rho2 <- total^2 * kappa
#' (n <- round(nmax(H_counts, N, S)) - 1)
#'
#' rfixprec(n, H_counts, N, S, rho)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
rfixprec <- function(n, H_counts, N, S, rho, rho2, U = NULL, J = NULL) {
  if (is.null(J)) {
    J <- seq_along(H_counts)
  }
  if (any(H_cnt2dind(H_counts)[U] %in% J)) {
    stop("U must not have any strata from domain that is in J.")
  }

  # strata global indices for domain J[1]
  H_j <- H_strt_idx_in_dom(H_counts, J[1])
  repeat {
    x <- if (length(J) == 1L) {
      fixprecact(
        n = n, H_counts = H_counts, N = N, S = S, rho = rho, rho2 = rho2, U = U
      )
    } else {
      rfixprec(
        n = n, H_counts = H_counts, N = N, S = S, rho = rho, rho2 = rho2,
        U = U, J = J[-1]
      )
    }
    U_Hj <- which(x[H_j] >= N[H_j])
    if (nonempty(U_Hj)) {
      U <- c(U, H_j[U_Hj])
      H_j <- H_j[-U_Hj]
    } else {
      break
    }
  }
  return(x)
}

#' Recursive FIXPREC algorithm (user function)
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains. Allocation preserves strata sizes.
#'
#' @inheritParams rfixprec
#' @param total totals of surveyed variable in domains.
#' @param kappa priority weights for domains.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' (n <- round(nmax(H_counts, N, S)) - 1)
#'
#' rfixprec_main(n, H_counts, N, S, total, kappa)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
rfixprec_main <- function(n, H_counts, N, S, total, kappa = NULL) {
  assert_integerish(H_counts, any.missing = FALSE, min.len = 1L)
  H_dcount <- length(H_counts)
  strata_count <- sum(H_counts)
  assert_number(n, finite = TRUE)
  assert_integerish(N, any.missing = FALSE, len = strata_count)
  assert_numeric(S, finite = TRUE, any.missing = FALSE, len = strata_count)
  assert_integerish(total, any.missing = FALSE, len = H_dcount)
  assert_numeric(kappa, finite = TRUE, any.missing = FALSE, len = H_dcount, null.ok = TRUE)
  assert_true(all(N > 0))
  assert_true(all(S > 0))
  assert_true(all(total > 0))
  assert_true(all(kappa > 0))
  assert_true(n > 0 && n <= sum(N))

  if (n == sum(N)) {
    return(N)
  }

  # convert to integer type
  H_counts <- as.integer(H_counts)
  n <- as.integer(n)
  N <- as.integer(N)
  total <- as.integer(total)

  # initiate parameters
  if (is.null(kappa)) {
    kappa <- rep(1 / length(H_counts), length(H_counts))
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa
  J <- seq_along(H_counts)

  rfixprec(
    n = n, H_counts = H_counts, N = N, S = S, rho = rho, rho2 = rho2,
    U = NULL, J = J
  )
}

## Previous versions of rfixprec ----

# Version that reduces the whole model.
rfixprec_modreduce <- function(n, H_counts, N, S, rho, rho2, J = NULL) {
  H <- H_cnt2glbidx(H_counts)
  if (is.null(J)) {
    J <- seq_along(H_counts)
  }
  # if(length(J) == 1 ) {
  #   cat("\n\tIn 1 \n\n")
  # }

  N0 <- N
  j <- J[1]
  H_j <- H[[j]] # strata global indices for domain j
  j1 <- H_j[1] # global index of first stratum in domain j
  H_j_len <- length(H_j)

  h_tocheck <- H_j
  repeat {
    x <- if (length(J) == 1L) {
      fixprec(n, H_counts, N, S, rho, rho2)
    } else {
      rfixprec_modreduce(n, H_counts, N, S, rho, rho2, J[-1])
    }

    # if(length(J) == 1 ) {
    #   ind <- (sum(H_counts[1:8])+1):((sum(H_counts[1:8])+1)+H_counts[4]-1)
    #   nonN <- which(x[ind] != N[ind])
    #   Adf <- (N[ind]*S[ind])/((total * sqrt(kappa))[4])
    #   s <- x[ind][nonN]/Adf[nonN]
    #   print(s[1])
    #   res <- abs(sum(diff(s))) < 1e-8
    # }

    h_violated <- which(x[h_tocheck] >= N[h_tocheck])
    if (nonempty(h_violated)) {
      n <- n - sum(N[h_tocheck[h_violated]])
      N <- N[-h_tocheck[h_violated]]
      S <- S[-h_tocheck[h_violated]]

      H_j <- H_j[-h_violated]
      H_counts[j] <- H_counts[j] - length(h_violated)
      if (H_counts[j] == 0L) { # whole domain j is (temporarily) active
        # print("BLOCK DOMAIN")
        H_counts <- H_counts[-j]
        rho <- rho[-j]
        rho2 <- rho2[-j]
        J <- J - 1
        h_tocheck <- integer(0)
      } else {
        h_tocheck <- seq(from = j1, len = H_counts[j])
      }
    }
    # } else {
    #   break
    # }
    if (empty(h_violated) || empty(H_counts)) {
      break
    }
  }
  H[[j]] <- H_j
  N0[unlist(H)] <- x
  N0
}

#' iterative implementation
#'
#' @param ref_domain reference domain (denoted by j in the paper).
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' N <- c(140, 110, 135, 190, 200, 40, 70)
#' S <- sqrt(c(180, 20, 5, 4, 35, 9, 40))
#' total <- c(2, 3, 5)
#' kappa <- c(0.5, 0.2, 0.3)
#' rho <- total * sqrt(kappa)
#' (n <- nmax(H_counts, N, S) - 1)
#'
#' rfixprec_iter(n, H_counts, N, S, rho)
#' # 140 103.60438 132.18060 166.39204 195.95002  19.87296  70
#'
rfixprec_iter <- function(n, H_counts, N, S, rho, rho2 = NULL, ref_domain = 1L) {
  if (is.null(rho2)) {
    rho2 <- rho^2
  }

  stopifnot(length(ref_domain) == 1L)
  stopifnot(ref_domain <= length(H_counts))

  H <- H_cnt2glbidx(H_counts) # list with strata global indices, list elements - domains
  H_ref <- H[[ref_domain]] # reference domain (denoted by j in the paper)
  H_nref <- H[-ref_domain] # remaining domains other than j
  H_nref_active <- vector("list", length(H_nref)) # active strata in remaining domains
  D_nref <- seq_along(H_nref) # indices of remaining domains

  repeat {
    # 1. Allocate in a chosen reference domain j
    H_active <- unlist(H_nref_active)
    repeat {
      x <- fixprecact(n, H_counts, N, S, rho, rho2, U = H_active)
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

#' Compute the value of the objective function and constraint functions
#' for a given allocation.
#'
#' @inheritParams rfixprec0
#'
#' @examples
#' H_dind <- c(1, 1, 2, 2) # 2 domains with 2 strata each
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' rho <- total * sqrt(kappa)
#' n <- 500
#'
#' (x <- fixprec(n, H_dind, N, S, rho))
#' compute_obj_cnstr(x, H_counts, N, S, total, kappa, n)
#' compute_obj_cnstr(x, H_counts, N, S, total, kappa, n, 1)
#' compute_obj_cnstr(x, H_counts, N, S, total, kappa, n, 2)
#' compute_obj_cnstr(x, H_counts, N, S, total, kappa, n, NULL)
#'
compute_obj_cnstr <- function(x, H_counts, N, S, total, kappa, n, J = numeric(0)) {
  # TODO: describe in the man that J can be empty.

  # Topt (objective)
  Td_notscaled <- tapply(N^2 / x * S^2 - N * S^2, H_cnt2dind(H_counts), sum)
  Topt <- sum(Td_notscaled) / sum(total^2 * kappa)

  # Td related (Td/kappa_d - Topt)
  Td <- Td_notscaled / total^2
  Td_o_kappad_min_T <- Td / kappa - Topt
  # note: Td / kappa = Td_notscaled / rho2, where rho2 <- total^2 * kappa

  # sum_x_n
  sumx_min_n <- sum(x) - n

  rval <- list(Topt = Topt, Td_o_kappad_min_T = Td_o_kappad_min_T, sumx_min_n = sumx_min_n)

  # N
  if (nonempty(J) || is.null(J)) {
    if (is.null(J)) {
      J <- seq_along(H_counts)
    }
    strata_check_N <- unlist(H_cnt2glbidx(H_counts)[J])
    x_min_N <- x[strata_check_N] - N[strata_check_N]
    names_x <- unlist(lapply(seq_along(H_counts), function(d) {
      paste0(H_cnt2glbidx(H_counts)[[d]], " (d=", d, ")")
    }))
    names(x_min_N) <- names_x[strata_check_N]
    rval <- c(rval, list(x_min_N = x_min_N))
  }

  rval
}

#' Compute the value of the objective function and constraint functions
#' for a given allocation.
#'
#' @inheritParams check_obj_cnstr
#' @param tol_max comparison maximum tolerance in adaptive selection of tolerance.
#'  See `is_equal` function.
#'
#' @examples
#' H_dind <- c(1, 1, 2, 2) # 2 domains with 2 strata each
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' rho <- total * sqrt(kappa)
#' n <- 500
#'
#' (x <- fixprec(n, H_dind, N, S, rho))
#' check_obj_cnstr(x, H_counts, N, S, total, kappa, n)
#' check_obj_cnstr(x, H_counts, N, S, total, kappa, n, 1)
#' check_obj_cnstr(x, H_counts, N, S, total, kappa, n, 2)
#' check_obj_cnstr(x, H_counts, N, S, total, kappa, n, NULL)
#'
check_obj_cnstr <- function(x, H_counts, N, S, total, kappa, n, J = numeric(0),
                            tol_max = 0.1) {
  oc <- compute_obj_cnstr(x, H_counts, N, S, total, kappa, n, J)

  is_Td_eq_kappa_T <- is_equal(oc$Td_o_kappad_min_T, 0, tol_max = tol_max)
  is_sum_x_eq_n <- is_equal(oc$sumx_min_n, 0, tol_max = tol_max)

  rval <- list(
    Topt = oc$Topt,
    is_sum_x_eq_n = is_sum_x_eq_n,
    is_Td_eq_kappa_T = is_Td_eq_kappa_T
  )

  if (nonempty(J) || is.null(J)) {
    is_x_leq_N <- oc$x_min_N <= 0
    rval <- c(rval, list(is_x_leq_N = is_x_leq_N))
  }

  rval
}

#' Check KKT conditions for problem of equal-precision optimal allocation in
#' single-stage sampling with domains and strata in domains (allocations
#' do not exceed strata sizes)
#'
#' @inheritParams rfixprec_main
#' @param active global indices of forced take-max strata.
#' @param s values of functions s_d, for each domain d (optional).
#' @param tol_max comparison maximum tolerance in adaptive selection of tolerance.
#'  See `is_equal` function.
#' @param details should detailed output be returned?
#'
#' @return If `details` is `FALSE`, it returns optimal value of T
#'  (if all constraints are satisfied) or `NA` (if any of the constraints is not
#'  satisfied). If `details` is `TRUE`, a detailed constraints check is returned.
#'
#' @examples
#' H_counts <- c(2, 2) # 2 domains with 2 strata each
#' H_dind <- H_cnt2dind(H_counts)
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' total <- c(2, 3)
#' kappa <- c(0.4, 0.6)
#' rho <- total * sqrt(kappa)
#' (n <- nmax(table(H_dind), N, S) - 1)
#'
#' (x <- fixprec(n, H_dind, N, S, rho))
#' check_kkt(x, H_counts, N, S, total, kappa, n)
#' check_kkt(x, H_counts, N, S, total, kappa, n, details = TRUE)
#'
check_kkt <- function(x, H_counts, N, S, total, kappa, n, J = NULL,
                      active = NULL, s = NULL, tol_max = 0.1, details = FALSE) {
  if (is.null(J)) {
    J <- seq_along(H_counts)
  }

  # value of the objective funct. and flags whether the constraints are satisfied
  oc <- check_obj_cnstr(x, H_counts, N, S, total, kappa, n, J, tol_max)

  # check if mu_ih is non-negative for active constraints
  if (empty(active)) {
    active <- which(x == N)
  } else {
    stopifnot(setequal(active, which(x == N))) # might be time consuming
  }

  is_mu_dh_nonneg <- if (empty(s)) {
    x[active] >= N[active]
  } else {
    rho <- rep(total * sqrt(kappa), times = H_counts)
    s <- rep(s, times = H_counts)
    s[active] >= rho[active] / S[active]
  }
  names(is_mu_dh_nonneg) <- active

  if (details) {
    list(
      Topt = oc$Topt,
      is_sum_x_eq_n = oc$is_sum_x_eq_n,
      is_Td_eq_kappa_T = oc$is_Td_eq_kappa_T,
      is_x_leq_N = oc$is_x_leq_N,
      is_mu_dh_nonneg = is_mu_dh_nonneg,
      active = active,
      active_len = length(active)
    )
  } else {
    if (all(oc$is_Td_eq_kappa_T, oc$is_sum_x_eq_n, oc$is_x_leq_N, is_mu_dh_nonneg)) {
      oc$Topt
    } else {
      NA_real_
    }
  }
}

#' Maximum allowed total sample size so that D is positive
#'
#' @inheritParams fixprec
#'
#' @details
#' See (16) from JW 2019. It is less than or equal to sum(N).
#'
#' @examples
#' H_counts <- c(2, 2) # two domains with 2 strata each.
#' N <- c(140, 110, 135, 190)
#' S <- sqrt(c(180, 20, 5, 4))
#' sum(N)
#' nmax(H_counts, N, S)
#'
nmax <- function(H_counts, N, S) {
  H_dind <- H_cnt2dind(H_counts)
  sum(tapply(N * S, H_dind, sum)^2 / tapply(N * S^2, H_dind, sum))
}

# HELPERS ----

#' Is the length of the vector equal to 0?
#'
#' @param x vector.
#'
empty <- function(x) {
  length(x) == 0L
}

#' Is the length of the vector greater than 0?
#'
#' @param x vector.
#'
nonempty <- function(x) {
  length(x) > 0L
}

#' Convert strata counts to stratum-domain indicator vector
#'
#' @param H_counts vector of strata counts in domains.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # 3 domains with 2, 2, and 3 strata respectively
#' H_cnt2dind(H_counts) # 1 1 2 2 3 3 3
#'
H_cnt2dind <- function(H_counts) {
  if (is.list(H_counts)) {
    rep(seq_along(H_counts), times = sapply(H_counts, length))
  } else {
    rep(seq_along(H_counts), times = H_counts)
  }
}

#' Convert strata counts to globally unique strata indices
#'
#' @param H_counts vector of strata counts in domains.
#'
#' @return A `list` with globally unique strata indices in domains. The d-th
#'   list element is a vector of strata global indices in domain d.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' H_cnt2glbidx(H_counts)
#'
H_cnt2glbidx <- function(H_counts) {
  mapply(
    function(to, len) seq.int(to = to, length.out = len, by = 1L),
    cumsum(H_counts), H_counts,
    SIMPLIFY = FALSE
  )
}

#' Get strata global indices in a given domain
#'
#' @param H_counts vector of strata counts in domains.
#' @param d domain index. It must be that `0 < d <= length(H_counts)`.
#'
#' @return A vector with globally unique strata indices in domain `d`.
#'
#' @examples
#' H_counts <- c(2, 2, 3) # three domains with 2, 2, and 3 strata respectively
#' H_strt_idx_in_domain(H_counts, 3)
#'
H_strt_idx_in_dom <- function(H_counts, d) {
  H_dind <- H_cnt2dind(H_counts)
  which(H_dind == d)
}

#' Check if `x1` and `x2` are equal up to some tolerance level.
#' The tolerance level is found adaptive in the range
#' `10^((-19):(tol_max))`
#'
#' @param x1 numeric vector, fist value to compare
#' @param x1 numeric vector, second value to compare
#' @param tol_max integer, specified as a power of `10`
#'
#' @examples
#' is_equal(c(3, 4), c(3, 4))
#' is_equal(c(3, 4), c(3.01, 4.11))
#' is_equal(c(3, 4), c(3.01, 4.11), tol_max = 0)
#'
is_equal <- function(x1, x2, tol_max = -1) {
  tolerances <- 10^((-19):(tol_max))

  abs_diff <- abs(x1 - x2)
  results <- rep(setNames(FALSE, paste0("tol=1e", tol_max)), length(abs_diff))
  to_check <- seq_along(abs_diff)

  for (tol in tolerances) {
    res <- which(abs_diff[to_check] <= tol)
    if (nonempty(res)) {
      results[to_check[res]] <- TRUE
      names(results)[to_check[res]] <- rep(paste0("tol=", tol), length(res))
      to_check <- to_check[-res]
    }
    if (all(results)) {
      break
    }
  }
  results
}

#' Check if all elements of a vector are not of the same sign
#'
#' @param x vector to check.
#'
#' @return `logical` `TRUE` if all not elements are of the same sign, and `FALSE`
#'   otherwise. `0` is neutral.
#'
#' @examples
#' diffsigns(1:5)
#' diffsigns(-(1:5))
#' diffsigns(c(-1, -2, 3))
#' diffsigns(c(0, -1))
#' diffsigns(c(0, 1))
#' diffsigns(c(0, 1, -1))
#'
diffsigns <- function(x) {
  any(x < 0) && any(x > 0)
}

# PROTOTYPES (under testing) ----

#' FIXPREC algorithm for upper bounds M <= N (prototype)
#'
#' Function for equal-precision optimal allocation in single-stage sampling
#' with domains and strata in domains.
#'
#' @inheritParams fixprec
#' @param M sample size upper bounds, defaults to N.
#' @param U take-max strata global indices.
#'
#' @examples
#' H_counts <- c(5, 2)
#' H_names <- rep(seq_along(H_counts), times = H_counts)
#' N <- c(100, 100, 100, 100, 100, 100, 100)
#' S <- c(154, 178, 134, 213, 124, 102, 12)
#' M <- c(80, 90, 70, 40, 10, 90, 100)
#' names(N) <- H_names
#' names(S) <- H_names
#' names(M) <- H_names
#' total <- c(13, 2)
#' kappa <- c(0.8, 0.2)
#' (rho <- total * sqrt(kappa))
#' n <- 150
#'
#' library(rAMPL)
#' source(file.path("functions/rfixprec0.R"))
#' source(file.path("ampl/ampl_fixprec.R"))
#' model <- file.path("ampl/ampl_fixprec.mod")
#' x_ampl <- ampl_fixprec(n, H_counts, N, S, total, kappa, model = model, M = M)
#' x_ampl$n_dh
#' # 12.754880 14.742653 11.098402 17.641490 10.000000 74.945462  8.817113
#'
#' fixprec_M(n, H_counts, N, S, total, kappa, M = M, U = 5)
#' #         1         1         1         1         1         2         2
#' # 12.754880 14.742653 11.098402 17.641490 10.000000 74.945462  8.817113
#'
fixprec_M <- function(n, H_counts, N, S, total, kappa = NULL, M = N, U = NULL) {
  H_dcount <- length(H_counts)
  if (is.null(kappa)) {
    kappa <- rep(1 / H_dcount, H_dcount)
  }
  rho <- total * sqrt(kappa)
  rho2 <- total^2 * kappa

  # stratum-domain indicators,
  # e.g. for 2 domains with 4 and 2 strata: 1 1 1 1 2 2
  H_di <- H_cnt2dind(H_counts)

  # ASSUMPTION: works only if U does not contain whole domain!
  if (is.null(U)) {
    n_adjusted <- n
    a.vec <- as.matrix(tapply(N * S, H_di, sum) / rho)
    c.vec <- tapply(N * S^2, H_di, sum) / rho2
  } else {
    n_adjusted <- n - sum(M[U])
    a.vec <- as.matrix(tapply(N[-U] * S[-U], H_di[-U], sum) / rho)
    c.vec_b <- vector(mode = "numeric", length = H_dcount)
    c.vec_b[H_di[U]] <- tapply(N[U]^2 * S[U]^2 / M[U], H_di[U], sum)
    c.vec <- (tapply(N * S^2, H_di, sum) - c.vec_b) / rho2
  }
  D.matrix <- (a.vec %*% t(a.vec)) / n_adjusted - diag(c.vec, nrow = length(c.vec))

  eigen_decomp <- eigen(D.matrix, symmetric = TRUE)
  lambda <- eigen_decomp$values[1] # largest eigenvalue
  v <- eigen_decomp$vectors[, 1] # corresponding eigenvector
  if (any(diff(sign(v)) != 0)) {
    stop("eigenvector containts entries of a different sign")
  }

  s.vec <- n_adjusted * v / as.numeric(t(a.vec) %*% v)
  A <- (N * S) / rep(rho, table(H_di)) # brackets due to finite-prec arithmetic!
  x <- rep(s.vec, table(H_di)) * A
  x[U] <- M[U]
  x
}
