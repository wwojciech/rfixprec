#
# This file contains R functions implementing the algorithms described
# in the paper: Jacek Wesołowski & Robert Wieczorkowski (2017) An eigenproblem
# approach to optimal equal-precision sample allocation in subpopulations,
# Communications in Statistics - Theory and Methods, 46:5, 2212-2231,
# DOI: 10.1080/03610926.2015.1040501
# To link to this article: http://dx.doi.org/10.1080/03610926.2015.1040501
#
# actualization date: 17.11.2017
#

#
# Function for equal-precision optimal allocation in subpopulations in two-stage
# sampling with Hartley-Rao pi-ps scheme at the first stage and the SRSWOR
# at the second, and with fixed SSU's sample size withins PSU,
# see Theorem 3.6 from the article
#
# Arguments:
#   m - sample size in terms of PSU
#   n - sample size in terms of SSU
#   data - data frame with input information for each stratum, subpopulation
#   and PSU (should contains identification variable for PSU)
#   J - name of the variable denoting subpopulations
#   H - name of the variable denoting strata
#   N_jhi - name of the variable with population sizes in terms of SSU
#   M_jh - name of the variable with population sizes in terms of PSU
#   z_jhi - name of the auxiliary variable \tilde{z}
#   t_jhi - name of the variable with totals of surveyed variable
#   S2_jhi - name of the variable with population variances of surveyed variable
#   omega_jhi - name of the variable omega
#   T_j - name of the variable with totals of surveyed variable
#   gamma_jh - name of the variable gamma
#   mcv - if given: vector of priority weights which modify common precision
#        (coefficient of variation) for subpopulations;
#        mcv is equal the square root of parameter \kappa from article
#
# Values:
#   output data frame contains input data appended with new columns 'm_jh' and 'n_jhi',
#   which give optimal allocation of numbers of PSU and SSU in subpulations,
#   strata and PSU's.
#
fixprec_HR_SRSWOR_n <- function(m, n, data, J = "sub", H = "h", N_jhi = "N_jhi", M_jh = "M_jh",
                                z_jhi = "z_jhi", t_jhi = "t_jhi", S2_jhi = "S2_jhi",
                                omega_jhi = "omega_jhi", T_j = "T_j", gamma_jh = "gamma_jh", mcv = NULL) {
  kappa <- as.data.frame(table(data[[J]]))
  names(kappa) <- c(J, "kappa")
  if (!is.null(mcv)) {
    kappa$kappa <- mcv^2
  } else {
    kappa$kappa <- 1
  }

  stopifnot(all(data[[gamma_jh]] > 0))

  d0 <- unique(data[c(J, H, gamma_jh)])
  A_jh <- data.frame(d0[, c(J, H)], A_jh = d0[[gamma_jh]])

  B_jhi <- data.frame(data[, c(J, H)], B_jhi = (data[[N_jhi]]^2) * data[[S2_jhi]] / ((data[[T_j]]^2) * data[[z_jhi]]))
  B_jh <- aggregate(B_jhi$B_jhi, list(B_jhi[[J]], B_jhi[[H]]), sum)
  names(B_jh) <- c(J, H, "B_jh")
  c_j <- aggregate(data[[omega_jhi]] * data[[z_jhi]] / (data[[T_j]]^2), list(data[[J]]), sum)
  names(c_j) <- c(J, "c_j")

  A_jh <- merge(A_jh, kappa, by = J, sort = FALSE)
  B_jh <- merge(B_jh, kappa, by = J, sort = FALSE)
  c_j <- merge(c_j, kappa, by = J, sort = FALSE)
  A_jh$A_jh <- A_jh$A_jh / A_jh$kappa
  B_jh$B_jh <- B_jh$B_jh / B_jh$kappa
  c_j$c_j <- c_j$c_j / c_j$kappa
  A_jh$kappa <- NULL
  B_jh$kappa <- NULL
  c_j$kappa <- NULL


  M_jh <- unique(data[, c(J, H, M_jh)])
  N_jhi <- data[, c(J, H, PSU, N_jhi)]
  N_jh <- aggregate(N_jhi$N_jhi, list(N_jhi[[J]], N_jhi[[H]]), sum)
  names(N_jh) <- c(J, H, "N_jh")

  nJ <- length(table(c_j[[J]]))
  a_j <- matrix(tapply(sqrt(A_jh$A_jh / m), A_jh[[J]], sum), nJ, 1)
  b_j <- matrix(tapply(sqrt(B_jh$B_jh / n), B_jh[[J]], sum), nJ, 1)

  if (nrow(c_j) > 1) {
    D.matrix <- a_j %*% t(a_j) + b_j %*% t(b_j) - diag(c_j$c_j)
  } else {
    D.matrix <- a_j %*% t(a_j) + b_j %*% t(b_j) - c_j$c_j
  }

  m1 <- eigen(D.matrix, symmetric = TRUE)
  lambda <- m1$values[1] # largest eigenvalue
  if (lambda <= 0) {
    stop("Largest eigenvalue is not strictly positive - solution does not exist!")
  }


  cv_optimal <- m1$values[1] # maximum eigenvalue
  # cat("CV optimal (in %) = ",100*sqrt(cv_optimal),"\n")
  cat("CV optimal (in %) for subpopulations:  \n ")
  print(cbind(kappa[[J]], round(100 * sqrt(kappa$kappa * cv_optimal), 2)))

  v <- (-m1$vectors[, 1]) # corresponding eigenvector
  # sum(v)

  v_j <- data.frame(J = c_j[[J]], v = v)
  names(v_j)[1] <- J
  # v_j[[J]]<-as.factor(v_j[[J]])

  A_jh <- merge(A_jh, v_j, by = J, sort = FALSE)

  tmp <- with(A_jh, v * sqrt(A_jh))
  m_jh <- (m * ((A_jh$v * sqrt(A_jh$A_jh)) / sum(tmp)))
  # m_jh<-round(m_jh)
  sum(m_jh)

  A_jh <- cbind(A_jh, m_jh = m_jh)
  B_jh <- merge(B_jh, A_jh)

  tmp <- with(B_jh, v * sqrt(B_jh))
  n_jh <- ((B_jh$v * sqrt(B_jh$B_jh)) / sum(tmp)) * n / B_jh$m_jh
  n_jh <- round(n_jh)
  n_jh <- pmax(n_jh, 1)
  sum(n_jh)

  A_jh$m_jh <- round(A_jh$m_jh)
  out <- merge(N_jh, M_jh)
  out <- cbind(out, n_jh = n_jh)
  out <- merge(out, A_jh[, c(J, H, "m_jh")])
  out <- out[order(out[[J]], out[[H]]), ]

  return(out)
}


#
# Function for equal-precision optimal allocation in subpopulations in two-stage
# sampling with Hartley-Rao pi-ps scheme at the first stage and the SRSWOR
# at the second, see Theorem 3.5 from the article
#
# Arguments:
#   m - sample size in terms of PSU
#   n - sample size in terms of SSU
#   data - data frame with input information for each stratum, subpopulation
#   and PSU (should contains identification variable for PSU)
#   J - name of the variable denoting subpopulations
#   H - name of the variable denoting strata
#   N_jhi - name of the variable with population sizes in terms of SSU
#   M_jh - name of the variable with population sizes in terms of PSU
#   z_jhi - name of the auxiliary variable \tilde{z}
#   t_jhi - name of the variable with totals of surveyed variable
#   S2_jhi - name of the variable with population variances of surveyed variable
#   omega_jhi - name of the variable omega
#   T_j - name of the variable with totals of surveyed variable
#   gamma_jh - name of the variable gamma
#   mcv - if given: vector of priority weights which modify common precision
#        (coefficient of variation) for subpopulations;
#        mcv is equal the square root of parameter \kappa from article
#
# Values:
#   output data frame contains input data appended with new columns 'm_jh' and 'n_jhi',
#   which give optimal allocation of numbers of PSU and SSU in subpulations,
#   strata and PSU's.
#
fixprec_HR_SRSWOR <- function(m, n, data, J = "sub", H = "h", N_jhi = "N_jhi", M_jh = "M_jh",
                              z_jhi = "z_jhi", t_jhi = "t_jhi", S2_jhi = "S2_jhi",
                              omega_jhi = "omega_jhi", T_j = "T_j", gamma_jh = "gamma_jh", mcv = NULL) {
  kappa <- as.data.frame(table(data[[J]]))
  names(kappa) <- c(J, "kappa")
  if (!is.null(mcv)) {
    kappa$kappa <- mcv^2
  } else {
    kappa$kappa <- 1
  }

  stopifnot(all(data[[gamma_jh]] > 0))

  d0 <- unique(data[c(J, H, gamma_jh)])
  A_jh <- data.frame(d0[, c(J, H)], A_jh = d0[[gamma_jh]])

  alfa_jhi <- data.frame(data[, c(J, H)], alfa_jhi = data[[z_jhi]])
  B_jhi <- data.frame(data[, c(J, H)], B_jhi = (data[[N_jhi]]^2) * data[[S2_jhi]] / ((data[[T_j]]^2) * data[[z_jhi]]))
  c_j <- aggregate(data[[omega_jhi]] * data[[z_jhi]] / (data[[T_j]]^2), list(data[[J]]), sum)
  names(c_j) <- c(J, "c_j")

  A_jh <- merge(A_jh, kappa, by = J, sort = FALSE)
  # A_jh<-A_jh[order(A_jh[[J]]),]

  B_jhi <- merge(B_jhi, kappa, by = J, sort = FALSE)
  # B_jhi<-B_jhi[order(B_jhi[[J]]),]


  c_j <- merge(c_j, kappa, by = J, sort = FALSE)
  # c_j<-c_j[order(c_j[[J]]),]

  A_jh$A_jh <- A_jh$A_jh / A_jh$kappa
  B_jhi$B_jhi <- B_jhi$B_jhi / B_jhi$kappa
  c_j$c_j <- c_j$c_j / c_j$kappa
  A_jh$kappa <- NULL
  B_jhi$kappa <- NULL
  c_j$kappa <- NULL


  M_jh <- unique(data[, c(J, H, M_jh)])
  N_jhi <- data[, c(J, H, PSU, N_jhi)]


  nJ <- length(table(c_j[[J]]))
  a_j <- matrix(tapply(sqrt(A_jh$A_jh / m), A_jh[[J]], sum), nJ, 1)
  b_j <- matrix(tapply(sqrt(alfa_jhi$alfa_jhi * B_jhi$B_jhi / n), alfa_jhi[[J]], sum), nJ, 1)

  if (nrow(c_j) > 1) {
    D.matrix <- a_j %*% t(a_j) + b_j %*% t(b_j) - diag(c_j$c_j)
  } else {
    D.matrix <- a_j %*% t(a_j) + b_j %*% t(b_j) - c_j$c_j
  }

  m1 <- eigen(D.matrix, symmetric = TRUE)
  lambda <- m1$values[1] # largest eigenvalue
  if (lambda <= 0) {
    stop("Largest eigenvalue is not strictly positive - solution does not exist!")
  }


  cv_optimal <- m1$values[1] # maximum eigenvalue
  # cat("CV optimal (in %) = ",100*sqrt(cv_optimal),"\n")
  cat("CV optimal (in %) for subpopulations:  \n ")
  print(cbind(kappa[[J]], round(100 * sqrt(kappa$kappa * cv_optimal), 2)))

  v <- (-m1$vectors[, 1]) # corresponding eigenvector
  # sum(v)

  v_j <- data.frame(J = c_j[[J]], v = v)
  names(v_j)[1] <- J
  # v_j[[J]]<-as.factor(v_j[[J]])

  A_jh <- merge(A_jh, v_j, by = J, sort = FALSE)

  tmp <- with(A_jh, v * sqrt(A_jh))
  m_jh <- (m * ((A_jh$v * sqrt(A_jh$A_jh)) / sum(tmp)))
  # m_jh<-round(m_jh)
  sum(m_jh)

  B_jhi <- cbind(B_jhi, alfa_jhi = alfa_jhi$alfa_jhi)
  A_jh <- cbind(A_jh, m_jh = m_jh)
  B_jhi <- merge(B_jhi, A_jh)

  tmp <- with(B_jhi, v * sqrt(alfa_jhi * B_jhi))
  n_jhi <- ((B_jhi$v * sqrt(B_jhi$B_jhi / B_jhi$alfa_jhi)) / sum(tmp)) * n / B_jhi$m_jh
  n_jhi <- round(n_jhi)
  n_jhi <- pmax(n_jhi, 1)
  sum(n_jhi)

  A_jh$m_jh <- round(A_jh$m_jh)
  out <- merge(N_jhi, M_jh)
  out <- cbind(out, n_jhi = n_jhi)
  out <- merge(out, A_jh[, c(J, H, "m_jh")])
  out <- out[order(out[[J]], out[[H]]), ]

  return(out)
}

#
# Function for equal-precision optimal allocation in single-stage sampling
# with subpopulations and strata, see Theorem 2.3 in the article
#
# Arguments:
#   n - sample size
#   data - data frame with input information for each stratum and subpopulation
#   J - name of the variable denoting subpopulations
#   H - name of the variable denoting strata
#   N_jh - name of the variable with population sizes
#   S2_jh - name of the variable with population variances of surveyed variable
#   t_j - name of the variable with totals of surveyed variable in subpopulations
#   mcv - if given: vector of priority weights which modify common precision
#        (coefficient of variation) for subpopulations;
#        mcv is equal the square root of parameter \kappa from article
#
# Values:
#   output data frame contains input data appended with new column 'n_jh',
#   which gives optimal allocation in subpulations and strata.
#
fixprec_SRSWOR <- function(n, data, J = "sub", H = "h", N_jh = "N_jh", S2_jh = "S2_jh", t_j = "t_j", mcv = NULL) {
  kappa <- as.data.frame(table(data[[J]]))
  names(kappa) <- c(J, "kappa")
  if (!is.null(mcv)) {
    kappa$kappa <- mcv^2
  } else {
    kappa$kappa <- 1
  }


  A_jh <- subset(data, select = c(J, H))
  A_jh <- merge(A_jh, kappa, by = J, sort = FALSE)
  A_jh$A_jh <- ((data[[N_jh]]^2) * data[[S2_jh]] / (data[[t_j]]^2)) / A_jh$kappa
  A_jh$kappa <- NULL

  c_jh <- subset(A_jh, select = c(J, H))
  c_jh$c_jh <- A_jh$A_jh / data[[N_jh]]

  A_j <- as.vector(tapply(sqrt(A_jh$A_jh), data[[J]], sum))
  c <- as.vector(tapply(c_jh$c_jh, data[[J]], sum))


  a <- matrix(A_j, length(A_j), 1)

  if ((nrow(diag(c)) > 1)) {
    D.matrix <- (1 / n) * (a %*% t(a)) - diag(c)
  } else {
    D.matrix <- (1 / n) * (a %*% t(a)) - c
  }

  m1 <- eigen(D.matrix, symmetric = TRUE)
  lambda <- m1$values[1] # largest eigenvalue
  if (lambda <= 0) {
    stop("Largest eigenvalue is not strictly positive - solution does not exist!")
  }

  cv_opt <- m1$values[1] # maximum eigenvalue
  # cat("cV optimal (in %)  = ",100*sqrt(cv_opt),"\n")
  cat("CV optimal (in %) for subpopulations:  \n ")
  print(cbind(kappa[[J]], round(100 * sqrt(kappa$kappa * cv_optimal), 2)))

  v <- (-m1$vectors[, 1]) # corresponding eigenvalue

  v <- v * (n / as.numeric((t(a) %*% v)))
  out <- as.data.frame(table(data[[J]]))
  names(out) <- c(J, "freq")
  out <- cbind(out, n_jh = v)
  out <- merge(data, out, by = J)
  out$n_jh <- round(sqrt(A_jh$A_jh) * out$n_jh) # optimal allocation
  out$freq <- NULL

  return(out)
}
