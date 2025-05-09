#' List with all possible subsets (of lengths len) of x.
#' @param x source subset
#' @param len vector with sizes of subsets `x`
subsets_len <- function(x, len = NULL) {
  if (is.null(len)) {
    len <- seq_len(length(x) - 1)
  }
  do.call(c, lapply(len, combn, x = x, simplify = FALSE))
}

#' Cross-domain subsets with strata
#'
#' @examples
#' H_counts <- c(2, 5, 3) # three domains with 2, 5, and 3 strata respectively.
#' subsets_domains(H_counts)
#'
subsets_domains <- function(H_counts) {
  H <- H_cnt2glbidx(H_counts)
  subsets_d <- lapply(H, subsets_len)
  subsets_d <- subsets_d[lengths(subsets_d) != 0] # remove NULL (occurs if 1s are in H_counts)
  if (length(subsets_d) == 0L) {
    return(NULL)
  }

  subsets1 <- do.call(c, subsets_d)
  H_counts <- H_counts[H_counts != 1] # remove domains with 1 stratum
  subsets2 <- if (length(H_counts) == 1L) {
    NULL
  } else {
    crossjoin_J <- subsets_len(seq_along(H_counts), 2:length(H_counts)) # combination of elements of J.
    subsets2 <- lapply(crossjoin_J, function(x) {
      cj <- do.call(expand.grid, subsets_d[x])
      colnames(cj) <- NULL
      split(cj, seq(nrow(cj)))
    })
    subsets2 <- do.call(c, subsets2)
    setNames(lapply(subsets2, unlist), NULL)
  }

  all_subsets <- c(list(NULL), subsets1, subsets2)

  # verification
  combn_in_d <- 2^H_counts - 2
  if (length(combn_in_d) == 1L) {
    combn_in_d <- c(0, combn_in_d) # since combn() works differently for scalar x.
  }
  no_all_subsets <- sum(
    1,
    sapply(seq_along(H_counts), function(d) sum(combn(combn_in_d, d, prod)))
  )
  if (length(all_subsets) != no_all_subsets) {
    stop("Zle obliczone podzbiory")
  }
  all_subsets
}
