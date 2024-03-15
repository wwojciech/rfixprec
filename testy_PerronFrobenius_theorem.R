# P-F theorem, wiecej niz jedna dodatnia wartosc wlasna. ----
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
