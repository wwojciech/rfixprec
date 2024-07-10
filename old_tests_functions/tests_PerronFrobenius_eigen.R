# Przyklad 1 ----

# P-F theorem, more than one strictly positive eigenvalue.
D <- structure(
  c(
    22L, 2L, 13L, 11L, 3L,
    2L, -6L, 17L, 12L, 19L,
    13L, 17L, 4L, 16L, 15L,
    11L, 12L, 16L, 9L, 10L,
    3L, 19L, 15L, 10L, 1L
  ),
  dim = c(5L, 5L)
)
eigen(D, symmetric = TRUE)

# D > 0
D <- structure(
  c(
    22L, 2L, 13L, 11L, 3L,
    2L, 6L, 17L, 12L, 19L,
    13L, 17L, 4L, 16L, 15L,
    11L, 12L, 16L, 9L, 10L,
    3L, 19L, 15L, 10L, 1L
  ),
  dim = c(5L, 5L)
)
eigen(D, symmetric = TRUE)

# Przyklad 2 ----

J <- c(5, 2)
strata_names <- rep(seq_along(J), times = J)
N <- c(100, 100, 100, 100, 100, 100, 100)
S <- c(154, 178, 213, 134, 124, 102, 12)
names(N) <- strata_names
names(S) <- strata_names
total <- c(13, 2)
kappa <- c(0.8, 0.2)
(rho <- total * sqrt(kappa))
n <- 695

D <- structure(
  c(-30564.1469073262, 126649.142224682, 126649.142224682, -1084758.99280576),
  dim = c(2L, 2L)
)
isSymmetric(D) # TRUE
eigenout <- eigen(D, symmetric = TRUE)
eigenout$values # -15562.23 -1099760.91
eigenout$vectors # (-0.9930575, -0.1176302), (-0.1176302, 0.9930575)

eigenout <- eigen(D, symmetric = FALSE)
eigenout$values # -1099760.91 -15562.23
eigenout$vectors # (-0.1176302, 0.9930575), (0.9930575, 0.1176302)
