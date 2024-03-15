library(rAMPL)
path <- file.path("/Volumes/Dane/studia/doktoranckie/alokacja_optymalna/R_code/testy_algorytmow_proby/fixprec/")

env <- new(Environment, "/Applications/AMPL")
ampl <- new(AMPL, env)

# Read model and the data.
ampl$read(file.path(path, "functions/ampl/ampl_fixprec.mod"))
ampl$readData(file.path(path, "functions/ampl/ampl_fixprec_3d.dat"))

# Solve it.
ampl$setOption("solver", "ipopt")
ampl$setOption("ipopt_options", "max_iter=500 print_level=0")
ampl$solve()

# Get objective entity by AMPL name.
Topt <- ampl$getObjective("base_variance")
# Print it.
cat(sprintf("Objective is: %g\n", Topt$value()))

# Get the values of the variable x.
nopt <- ampl$getVariable("x")
(nopt <- nopt$getValues()[, "x.val"])

# Stop the AMPl engine.
ampl$close()
