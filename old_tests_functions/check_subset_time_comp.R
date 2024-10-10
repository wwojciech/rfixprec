A <- sample(10^6, 100, replace = TRUE)
B <- sample(10^7, 50, replace = TRUE)

microbenchmark::microbenchmark(
  a = setequal(intersect(A, B), A),
  b = all(B %in% A),
  times = 10,
  unit = "us"
)
