# library(rAMPL)

ampl_fixprec <- function(n, J, N, S, total, kappa,
                         data = NULL,
                         model = "ampl_fixprec.mod",
                         print_values = FALSE, short_ampl_output = TRUE,
                         ipopt_options = "max_iter=500 print_level=0") {
  # Setup ampl env and read model file.
  env <- new(Environment, "/Applications/AMPL")
  ampl <- new(AMPL, env)
  ampl$read(model)

  if (!is.null(data)) {
    ampl$readData(data)
  } else {
    # Set values of the parameters.
    set_DOMAINS <- ampl$getSet("DOMAINS")
    set_DOMAINS$setValues(seq_along(J))
    # STRATA
    set_STRATA <- ampl$getSet("STRATA")
    set_STRATA <- set_STRATA$getInstances()
    for (i in seq_along(set_STRATA)) {
      set_STRATA[[i]]$setValues(seq_len(J[i]))
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
  }

  if (print_values) {
    print(set_DOMAINS$getValues())
    print(sapply(set_STRATA, function(i) i$getValues()[[1]]))
    print(cbind(param_N$getValues(), S = param_S$getValues()[, "S"]))
    print(cbind(param_total$getValues(), kappa = param_kappa$getValues()[, "kappa"]))
    cat(sprintf("Total sample size n: %g\n", param_n$getValues()))
  }
  # ampl$exportData("ampl_fixprec_9d_2.dat")

  # Set options and solve the problem.
  # ampl$setOption("show_stats", 0)
  # ampl$setOption("solver_msg", 1)
  ampl$setOption("solver", "ipopt")
  ampl$setOption("ipopt_options", ipopt_options)

  if (short_ampl_output) {
    output_handler <- function(x) {
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
  nopt <- var_x$getValues()[, "x.val"]

  # Stop the AMPl engine.
  ampl$close()

  list(Topt = Topt, n_ih = nopt, ampl_out = ampl_out)
}

ampl_readData <- function(data = "ampl_fixprec_9d_2.dat", model = "ampl_fixprec.mod") {
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

  J <- sapply(set_STRATA, function(i) nrow(i$getValues()))
  names(J) <- NULL
  N <- param_N$getValues()[, "N"]
  S <- param_S$getValues()[, "S"]
  total <- param_total$getValues()[, "total"]
  kappa <- param_kappa$getValues()[, "kappa"]
  n <- param_n$getValues()[, "n"]

  # Stop the AMPl engine.
  ampl$close()

  list(J = J, n = n, N = N, S = S, total = total, kappa = kappa)
}
