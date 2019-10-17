# preparations
load("mynetwork.RData")
source(file.path(dirx, "libs.R"))
source(file.path(dirx, "opt_cumo_tools.R"))
source(file.path(dirx, "opt_icumo_tools.R")) # uncoment for influx_i use
tmp=sparse2spa(spa)

# doing something useful
# here, we calculate a vector of cost values, one per MC iteration
#cost_mc=apply(free_mc, 2, function(p) cumo_cost(p, labargs))
# do something else ...
