# Przyklad (2 domeny) ----

# cala domena 1 od razu wyrzucona

# Tutaj tez ciekawostka: powiedzmy, ze J = {1}.
# Zaczynam podejrzewac, ze w tym przypadku minimum nie istnieje, gdy n jest takie,
# ze cala pierwsza domena jest przekroczona.
# Troche analizy: |J| >= 1 => T >= 0 (zachodzi z powodu wiezow rownosciowych na T w domenach).
# czyli jesli np |J| = 1 (nasz przypadek) to T >= 0.
# ampl_fixprec(600, H_ss, N, S, total, kappa, model = model, J = 1)

H_ss <- c(2, 2)
H_names <- rep(seq_along(H_ss), times = H_ss)
N <- c(50, 50, 200, 456)
S <- c(254, 178, 4, 1750)
names(N) <- H_names
names(S) <- H_names
total <- c(13, 8)
kappa <- c(0.6, 0.4)
(rho <- total * sqrt(kappa))
nmax(H_ss, N, S) # 553.9116
sum(N) # 850
n <- 750

# domena I
(x <- setNames(fixprec_act(n, H_ss, N, S, total, kappa), H_names))
check_obj_cnstr(x, H_ss, N, S, total, kappa, n)
(x <- setNames(fixprec_act(n, H_ss, N, S, total, kappa, 1), H_names))
check_obj_cnstr(x, H_ss, N, S, total, kappa, n)

(x <- setNames(fixprec_act(n, H_ss, N, S, total, kappa, 1:2), H_names)) # cala domena 1 wyrzucona
check_obj_cnstr(x[-(1:2)], H_ss[-1], N[-(1:2)], S[-(1:2)], total[-1], kappa[-1], n - sum(N[1:2]))

# domena II
# w d = 2, przekroczona h = 2 (ind = 4)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, 4), H_names)

# domena I
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(4, 1)), H_names)

# Obserwacja - jesli wszystkie w domenie 1 sa przekroczone, to
# T < 0 i wtedy w kazdej pozostalej domenie istnieje co najmniej jedna warstawa przekroczona
# Kiedy odejme cala domene, to nadwyzka z pierwszej domeny,
# musi przejsc na pozostale domeny (z powodu wiezu na n). Tylko jak przechodzi?
# Czy na kazda domene czy tylko na jedna?


# Przyklad (3 domeny) ----

H_ss <- c(5, 2, 4)
H_names <- rep(seq_along(H_ss), times = H_ss)
N <- c(50, 50, 50, 50, 50, 100, 100, 100, 100, 100, 100)
S <- c(154, 178, 134, 213, 124, 102, 12, 34, 90, 30, 200)
names(N) <- H_names
names(S) <- H_names
total <- c(13, 6, 8)
kappa <- c(0.5, 0.2, 0.3)
(rho <- total * sqrt(kappa))
nmax(H_ss, N, S) - 1 # 612
sum(N) # 850
n <- 800

# domena I
setNames(fixprec_act(n, H_ss, N, S, total, kappa), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, 1:5), H_names)

# domena II
# w d = 2, przekroczona h = 1 (ind = 6)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, 6), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4, 1)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4, 1, 3, 5)), H_names)
# w d = 2, przekroczona h = 2 (ind = 7)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 7)), H_names) # cala II-ga domena zablokowana
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 7, 1:5)), H_names)

# domena III
# w d = 3, przekroczona h = 4 (ind = 11)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 1, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 1, 2, 4, 3, 5)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 2, 4, 1)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 2, 4, 1, 3)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 2, 4, 1, 3, 5)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 7)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 7, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 7, 2, 4, 1, 3)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 6, 7, 2, 4, 1, 3, 5)), H_names)
# w d = 3, przekroczona h = 2 (ind = 9)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4, 1, 3)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4, 1, 3, 5)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6, 2, 4, 1)), H_names)

ampl_fixprec(n, H_ss, N, S, total, kappa, model = model)$n_dh
rfixprec(n, H_ss, N, S, total, kappa)

## D matrix ----
rho2 <- total^2 * kappa
H_di <- H_domain_indicators(H_ss)
a.vec <- as.matrix(tapply(N * S, H_di, sum) / rho)
c.vec <- tapply(N * S^2, H_di, sum) / rho2 # - b (b = 0 if M = N)
D.matrix <- (a.vec %*% t(a.vec)) / n - diag(c.vec, nrow = length(c.vec))
D.matrix
# 1          2          3
# 1 -55503.22   23195.61   44108.25
# 2  23195.61 -123937.50   42904.34
# 3  44108.25   42904.34 -179643.23
eigen(D.matrix, symmetric = TRUE)
# eigen() decomposition
# $values
# [1]  -27618.68 -122930.11 -208535.15
#
# $vectors
# [,1]       [,2]       [,3]
# [1,] 0.8616444  0.4675550 -0.1973859
# [2,] 0.3647094 -0.8409070 -0.3998282
# [3,] 0.3529249 -0.2725212  0.8950845

D.matrix_d1_0 <- D.matrix
D.matrix_d1_0[, 1] <- 0
D.matrix_d1_0[1, ] <- 0
eigen(D.matrix_d1_0, symmetric = TRUE)
# eigen() decomposition
# $values
# [1]       0.0 -100638.0 -202942.7
#
# $vectors
# [,1]      [,2]       [,3]
# [1,]    1 0.0000000  0.0000000
# [2,]    0 0.8787798 -0.4772275
# [3,]    0 0.4772275  0.8787798

D.matrix_no_d1 <- D.matrix[, -1]
D.matrix_no_d1 <- D.matrix_no_d1[-1, ]
eigen(D.matrix_no_d1, symmetric = TRUE)
# eigen() decomposition
# $values
# [1] -100638.0 -202942.7
#
# $vectors
# [,1]       [,2]
# [1,] -0.8787798 -0.4772275
# [2,] -0.4772275  0.8787798

active <- 1:5
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
D.matrix
D.matrix_no_d1 # rozne bo w D.matrix zmienia sie jeszcze n.

# Przyklad (3 domeny) ----

H_ss <- c(5, 2, 4)
H_names <- rep(seq_along(H_ss), times = H_ss)
N <- c(50, 50, 50, 50, 50, 100, 100, 100, 50, 100, 50)
S <- c(154, 178, 134, 213, 124, 102, 12, 34, 90, 30, 200)
names(N) <- H_names
names(S) <- H_names
total <- c(13, 6, 8)
kappa <- c(0.5, 0.2, 0.3)
(rho <- total * sqrt(kappa))
nmax(H_ss, N, S) - 1 # 529
sum(N) # 750
n <- 700

# domena I
setNames(fixprec_act(n, H_ss, N, S, total, kappa), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, 1:5), H_names)

# domena II
# w d = 2, przekroczona h = 1 (ind = 6)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, 6), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4, 1)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 2, 4, 1, 3, 5)), H_names)
# w d = 2, przekroczona h = 2 (ind = 7)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 7)), H_names) # cala II-ga domena zablokowana
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 7, 1:4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(6, 7, 1:4, 5)), H_names)

# domena III
# w d = 3, przekroczona h = 2,4 (ind = 9,11)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4, 1, 3)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 2, 4, 1, 3, 5)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6, 2, 4)), H_names)
setNames(fixprec_act(n, H_ss, N, S, total, kappa, c(11, 9, 6, 2, 4, 1)), H_names)

ampl_fixprec(n, H_ss, N, S, total, kappa, model = model)$n_dh
rfixprec(n, H_ss, N, S, total, kappa)
