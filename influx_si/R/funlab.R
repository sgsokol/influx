# common functions for variable label input

# step labeling in a linear metabolic pathway with all different exchange rates 'nu' (nu=flux/metab_concentration).
# all metabolites are labeled at 0 at t=0
# return a matrix where one column corresponds to a given time point from tp and a row corresponds to a given metabolite
steplinpath=function(tp, nu) {
    nb_tp=length(tp)
    nb_nu=length(nu)
    # trivial case
    if (nb_tp == 0 || nb_nu == 0)
        return(matrix(NA, nrow=nb_nu, ncol=nb_tp))
    # test for all nu being > 0 and pairwise different
    stopifnot(all(nu >= 0.))
    diff_nu=outer(nu, nu, "!=")
    diff_nu=diff_nu[lower.tri(diff_nu)]
    stopifnot(all(diff_nu))
    # a is a list of coefficient vectors a[[imet]][knu]
    a=lapply(seq_along(nu), function(imet) sapply(seq(imet), function(knu) {v=nu[seq(imet)][-knu]; -prod(v/(v-nu[knu]))}))
    expt=exp(-nu%o%tp)
    res=t(1.+vapply(seq_along(nu), function(imet) {
        v=a[[imet]]
        colSums(v*expt[seq_along(v),, drop=FALSE])
    }, double(nb_tp)))
    dim(res)=c(nb_nu, nb_tp)
    colnames(res)=tp
    structure(res, a=a)
}

# more general than steplinpath2, arbitary initial conditions in 'init' and arbitrary heihgt of input labeling ('height')
steplinpath2=function(tp, nu, init=double(length(nu)), height=1.) {
    nb_tp=length(tp)
    nb_nu=length(nu)
    # trivial case
    if (nb_tp == 0 || nb_nu == 0)
        return(matrix(NA, nrow=nb_nu, ncol=nb_tp))
    # test for all nu being > 0 and pairwise different
    stopifnot(all(nu >= 0.))
    diff_nu=outer(nu, nu, "!=")
    diff_nu=diff_nu[lower.tri(diff_nu)]
    stopifnot(all(diff_nu))
    # a is a list of coefficient vectors a[[imet]][knu]
    a=Reduce(function(li, imet) {av=li[[imet-1]]*(nu[imet]/(nu[imet]-nu[seq_len(imet-1L)])); alast=init[imet]-(height+sum(av)); c(li, list(c(av, alast)))}, seq_len(nb_nu-1L)+1L, list(init[1L]-height))
    expt=exp(-nu%o%tp)
    res=t(height+vapply(seq_along(nu), function(imet) {
        v=a[[imet]]
        colSums(v*expt[seq_along(v),, drop=FALSE])
    }, double(nb_tp)))
    dim(res)=c(nb_nu, nb_tp)
    colnames(res)=tp
    structure(res, a=a)
}

# periodic step pulses with step intervals in 'Tint' and heights in 'Hint'. Initial conditions are in 'init'.
# The full duration of a complete period is sum(Tint)
ppulseslinpath=function(tp, nu, Tint, Hint=rep_len(c(1., 0.), length(Tint)), init=double(length(nu))) {
    stopifnot(length(Tint) == length(Hint))
    nb_tp=length(tp)
    # test for all nu being > 0 and pairwise different
    stopifnot(all(nu >= 0.))
    diff_nu=outer(nu, nu, "!=")
    diff_nu=diff_nu[lower.tri(diff_nu)]
    stopifnot(all(diff_nu))
    # partition tp in periods. Each period is composed of intervals
    intstart=c(0., cumsum(Tint))
    nb_int=length(Tint)
    tperiod=intstart[nb_int+1L]
    nb_period=max(ceiling(max(tp)/tperiod), 1L)
    res=NULL # will accumulate results
    for (iperiod in seq_len(nb_period)) {
        pstart=tperiod*(iperiod-1L)
        for (iint in seq_len(nb_int)) {
            tstart=pstart+intstart[iint]
            tphere=c(tp[tp >= tstart & tp < tstart+Tint[iint]], tstart+Tint[iint])-tstart
            # take init from the previous interval
            inithere=if (iperiod == 1L && iint == 1) init else reshere[,dim(reshere)[2L]]
            reshere=steplinpath2(tphere, nu, inithere, Hint[iint])
            colnames(reshere)=tphere+tstart
            res=cbind(res, reshere[,-length(tphere), drop=FALSE])
            if (iperiod == nb_period && iint == nb_int && ncol(res) < nb_tp)
                res=cbind(res, reshere[,dim(reshere)[2L], drop=FALSE])
        }
    }
    return(res)
}

# linear interpolation at time points tp with linear piece-wise function defined by x='knots' and y='v'
linterp=function(tp, knots, v) {
    rtp=range(tp)
    rk=range(knots)
    nb_kn=length(knots)
    stopifnot(rtp[1L] >= rk[1L] && rtp[2L] <= rk[2L])
    # find index of the closest left knot for all tp
    # knots are supposed being ordered monotonously and increasingly
    il=apply(outer(tp, knots, ">="), 1L, function(v) rev(which(v))[1L])
    il[il==nb_kn]=nb_kn-1L
    ir=il+1L # right knot
    # l(w)=(1-w)*y_left+w*y_right, w in [0,1]
    w=(tp-knots[il])/(knots[ir]-knots[il])
    v[il]+w*(v[ir]-v[il])
}

# periodic rectangular pulses. Duration of each time interval is defined in Tint numeric vector. Pulse amplitudes are in Hint, by defaults it's a repetition of 1 and 0 as many times as there are intervals in Tint.
ppulses=function(tp, Tint, Hint=rep_len(c(1., 0.), length(Tint))) {
    stopifnot(length(Tint) == length(Hint))
    nb_tp=length(tp)
    if (nb_tp == 0)
        return(double(0))
    tol=.Machine$double.eps*2**7
    # partition tp in periods. Each period is composed of intervals
    intstart=c(0., cumsum(Tint))
    nb_int=length(Tint)
    tperiod=intstart[nb_int+1L]
    nb_period=max(ceiling(max(tp)/tperiod), 1L)
    res=NULL # will accumulate results
    for (iperiod in seq_len(nb_period)) {
        #cat("ip=", iperiod, "\n")
        pstart=tperiod*(iperiod-1L)
        for (iint in seq_len(nb_int)) {
            #cat("iint=", iint, "\n")
            tstart=pstart+intstart[iint]
            if (abs(tp[nb_tp]-tstart-Tint[iint]) < tol) {
                tphere=tp[tp >= tstart-tol & tp <= tstart+Tint[iint]+tol]
            } else {
                tphere=tp[tp >= tstart & tp < tstart+Tint[iint]]
            }
            #cat("tph=", tphere, "\n")
            reshere=rep_len(Hint[iint], length(tphere))
            names(reshere)=tphere
            res=c(res, reshere)
            if (length(res) == nb_tp)
                return(res)
        }
    }
    stop("oops. Should not be there.")
}
