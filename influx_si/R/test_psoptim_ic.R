# test for psoptim_ic
dirR="/home/sokol/dev/sysbio/ftbl2sys/influx_si/R"
source(file.path(dirR, "nlsic.R"))
source(file.path(dirR, "opt_cumo_tools.R"))
source(file.path(dirR, "psoptim_ic.R"))
fn=function(x) {
#browser()
   sum(x**2)+10*sum(sin(x)**2)+10
}

# scalar case
x=seq(-5, 5, length.out=201)
plot(x, sapply(x, fn), t="l")

fit=psoptim_ic(NA, fn)

# 2D case
n=2L
z=outer(x, x, function(x,y) sapply(seq_along(x), function(i) fn(c(x[i],y[i]))))
plot3D::persp3D(z=z)
fit=psoptim_ic(c(NA,NA), fn, mean=10, control=list(trace=1, maxit=30))

# 2D box case x in [-4,-1], y in [2; 6]
n=2L
ui=rbind(diag(n), -diag(n))
ci=c(-4, 2, 1, -6)
x=seq(ci[1L], -ci[1L+n], length.out=101)
y=seq(ci[2L], -ci[2L+n], length.out=101)
z=outer(x, y, function(x,y) sapply(seq_along(x), function(i) fn(c(x[i],y[i]))))
plot3D::persp3D(x, y, z=z)
# the mean+-magnitude is completely outside of the feasibility domain
#set.seed(7)
fit=psoptim_ic(c(NA,NA), fn, mean=100*c(-2.5, 4), magnitude=10, control=list(trace=1, maxit=30, REPORT=5, maxit.stagnate=3, reltol=1.e-2), ui=ui, ci=ci)
