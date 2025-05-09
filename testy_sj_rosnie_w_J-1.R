H_counts <- c(3, 3)
H_names <- rep(seq_along(H_counts), times = H_counts)
N <- c(40, 100, 100, 100, 100, 100)
S <- c(68, 50, 26, 23, 24, 90)
names(N) <- H_names
names(S) <- H_names
total <- c(20, 10)
kappa <- c(0.4, 0.6)
rho <- total * sqrt(kappa)
rho2 <- total^2 * kappa
A <- (N * S) / rep(rho, H_counts)
nmax(H_counts, N, S) # 415.8198
sum(N) # 540
n <- 400

# zaczynam od domena 2
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 6), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt # 888.9765
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 1), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt # 275.5272
# Czyli T sie zmniejszylo

# zaczynam od domena 1
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 1), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt # 275.5272
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 1:2), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt # 283.9586

(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 6), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt # 888.9765
# Czyli T sie zwiekszylo

(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 1), H_names))
check_obj_cnstr(x, H_counts, N, S, total[-1], kappa[-1], n)$Topt

# domena II
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = 6), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt

# domena I
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 2)), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt

# domena II
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 4)), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt

# domena I
(x <- setNames(fixprecact(n, H_counts, N, S, rho, rho2, U = c(6, 4, 2)), H_names))
check_obj_cnstr(x, H_counts, N, S, total, kappa, n)$Topt
