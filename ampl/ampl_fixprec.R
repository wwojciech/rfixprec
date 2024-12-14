# MAIN FUNTIONS ----

#' Solve (using AMPL) equal-precision optimal allocation problem in single-stage
#' sampling with domains and strata in domains. Allocations will not exceed
#' specified upper bounds.
#'
#' @param n total sample size
#' @param H_ss vector of strata sizes in domains
#' @param N population sizes
#' @param S population standard deviations of surveyed variable
#' @param total totals of surveyed variable in domains
#' @param kappa priority weights for domains
#' @param M optional upper bounds, defaults to `N`.
#' @param J vector of domain indices. Specifies domains for which the allocated
#'   samples should preserve bounds specified in `M`. For domains other than
#'   those specified in `J`, allocations may exceed `M`.
#' @param data file path to .dat file (instead of `n`, `H_ss`, `N`, `S`, `total`,
#'   `kappa`, `M`.
#' @param model file path to AMPL model specification file
#' @param print_values should parameter values be printed out?
#' @param output_handler AMPL output handler function, defaults to
#' .  `ampl_output_handler_ipopt`
#' @param solver AMPL solver name, defaults to "ipopt".
#' @param solver_options solver options.
#'
#' @return `list` with the solution and AMPL output details
#'
#' @examples
#' library(rAMPL)
#' model <- model <- file.path("ampl/ampl_fixprec.mod")
#'
#' data <- file.path("ampl/data/ampl_fixprec_2d_2_M.dat")
#' x_ampl <- ampl_fixprec(data = data, model = model)
#' x_ampl
#'
#' H_ss <- c(5, 2)
#' N <- c(100, 100, 100, 100, 100, 100, 100)
#' S <- c(154, 178, 134, 213, 124, 102, 12)
#' total <- c(13, 2)
#' kappa <- c(0.8, 0.2)
#' n <- 623
#' x_ampl <- ampl_fixprec(n, H_ss, N, S, total, kappa, model = model)
#' x_ampl
#'
ampl_fixprec <- function(n, H_ss, N, S, total, kappa, M = N, J = NULL,
                         data = NULL,
                         model = "ampl_fixprec.mod",
                         print_values = FALSE,
                         output_handler = ampl_output_handler_ipopt,
                         solver = "ipopt",
                         solver_options = "max_iter=500 print_level=0") {
  # Setup the AMPL env. and read model file.
  env <- new(Environment, "/Applications/AMPL")
  ampl <- new(AMPL, env)
  ampl$read(model)

  if (!is.null(data)) {
    ampl$readData(data)
  } else {
    # Set values of the parameters.
    set_DOMAINS <- ampl$getSet("DOMAINS")
    set_DOMAINS$setValues(seq_along(H_ss))
    # STRATA
    set_STRATA <- ampl$getSet("STRATA")
    set_STRATA <- set_STRATA$getInstances()
    for (d in seq_along(set_STRATA)) {
      set_STRATA[[d]]$setValues(seq_len(H_ss[d]))
    }
    # J
    if (!is.null(J)) {
      set_DOMAINS_J <- ampl$getSet("DOMAINS_J")
      set_DOMAINS_J$setValues(J)
    }
    # N
    param_N <- ampl$getParameter("N")
    param_N$setValues(N)
    # S
    param_S <- ampl$getParameter("S")
    param_S$setValues(S)
    # total
    param_total <- ampl$getParameter("total")
    param_total$setValues(total)
    # kappa
    param_kappa <- ampl$getParameter("kappa")
    param_kappa$setValues(kappa)
    # n
    param_n <- ampl$getParameter("n")
    param_n$setValues(n)
    # M
    param_M <- ampl$getParameter("M")
    param_M$setValues(M)
  }

  if (print_values) {
    print(set_DOMAINS$getValues())
    print(sapply(set_STRATA, function(d) d$getValues()[[1]]))
    print(set_DOMAINS_J$getValues())
    print(cbind(param_N$getValues(), param_M$getValues(), S = param_S$getValues()[, "S"]))
    print(cbind(param_total$getValues(), kappa = param_kappa$getValues()[, "kappa"]))
    cat(sprintf("Total sample size n: %g\n", param_n$getValues()))
  }
  # ampl$exportData("ampl_fixprec_4d.dat")

  # Set options and solve the problem.
  # ampl$setOption("show_stats", 0)
  # ampl$setOption("solver_msg", 1)
  ampl$setOption("solver", solver)
  if (is.character(solver_options)) {
    ampl$setOption(paste0(solver, "_options"), solver_options)
  }
  if (is.function(output_handler)) {
    ampl$setOutputHandler(output_handler)
  }

  # Solve the problem.
  ampl_out <- capture.output(ampl$solve())
  if (length(ampl_out) > 0) {
    cat(ampl_out, sep = "\t\n")
  }

  # Get variables.
  var_base_variance <- ampl$getObjective("base_variance")
  var_x <- ampl$getVariable("x")
  Topt <- var_base_variance$value()
  n_opt <- var_x$getValues()[, "x.val"]

  # Stop the AMPL engine.
  ampl$close()

  list(Topt = Topt, n_dh = n_opt, ampl_out = ampl_out)
}

#' Read data from AMPL .dat file
#'
#' @param data file path to .dat file
#' @param model file path to AMPL model specification file
#'
#' @return `list` with data objects
#'
#' @examples
#' library(rAMPL)
#' model <- model <- file.path("ampl/ampl_fixprec.mod")
#' data <- file.path("ampl/data/ampl_fixprec_2d_2_M.dat")
#' data <- ampl_readData(data, model)
#' data
#'
ampl_readData <- function(data = "ampl_fixprec_9d_2.dat", model = "ampl_fixprec.mod") {
  # Setup the AMPL env. and read model and data files.
  env <- new(Environment, "/Applications/AMPL")
  ampl <- new(AMPL, env)
  ampl$read(model)
  ampl$readData(data)

  # Get parameters.
  set_STRATA <- ampl$getSet("STRATA")
  set_STRATA <- set_STRATA$getInstances()
  param_N <- ampl$getParameter("N")
  param_S <- ampl$getParameter("S")
  param_total <- ampl$getParameter("total")
  param_kappa <- ampl$getParameter("kappa")
  param_n <- ampl$getParameter("n")
  param_M <- ampl$getParameter("M")

  H_ss <- sapply(set_STRATA, function(i) nrow(i$getValues()))
  names(H_ss) <- NULL
  N <- param_N$getValues()[, "N"]
  S <- param_S$getValues()[, "S"]
  total <- param_total$getValues()[, "total"]
  kappa <- param_kappa$getValues()[, "kappa"]
  n <- param_n$getValues()[, "n"]

  M <- tryCatch(
    {
      param_M$getValues()[, "M"] # it throws the error if M not specified
    },
    error = function(e) {
      NULL
    }
  )

  # Stop the AMPL engine.
  ampl$close()

  list(H_ss = H_ss, n = n, M = M, N = N, S = S, total = total, kappa = kappa)
}

# HELPERS ----

#' AMPL Output handler
#'
#' The function handling the AMPL output derived from interpreting user commands.
#' Used by `ampl_fixprec` function. See `AMPL` documentation for more details on
#' setting a new output handler.
#'
ampl_output_handler_ipopt <- function(x) {
  pattern_ipopt <- "^ \\nIpopt ([[:digit:]])+\\.([[:digit:]])+\\.([[:digit:]])+: "
  pattern_solfound <- paste0(pattern_ipopt, "Optimal Solution Found")
  if (any(grepl(pattern_ipopt, x)) && !any(grepl(pattern_solfound, x))) {
    cat(gsub(pattern_ipopt, "", x))
  }
  pattern_check <- "^.*check:"
  check <- grepl(pattern_check, x)
  if (any(check)) {
    cat(gsub(pattern_check, "check failed:", x[check]))
  }
}
