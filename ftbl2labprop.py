#!/usr/bin/env python

"""
Transform an ftbl to R code which will solve an optimization of
labeling propagation problem in pulse experiment:
min(S) over \Theta, where S=||Predicted-Observed||^2_\Sigma^2
here \Theta is a vector of free fluxes (net+xch), metabolite pools
and scaling parameters.
Predicted vector is obtained from cumomer vector x(t) (calculated from
free fluxes and divided in chunks according to the cumo weight) by
multiplying it by the measurement matrices and scale factors, boths coming
from ftbl file. Observed values vector xo is extracted from datafile
whose name is taken from ftbl file.
it is composed of cumomer measurements.
Flux measurments are taken from ftbl file.
\Sigma^2, covariance diagonal matrices sigma[flux|mass|label|peak]
are orginated from ftbl
The system Afl*flnx=bfl is created from ftbl file.

usage: ./ftbl2labprop.py organism
where organism is the ftbl informative part of file name
(before .ftbl), e.g. organism.ftbl
After execution a file organism.R will be created.
If it already exists, it will be silently overwritten.

This code is derived from ftbl2optR.py which generated R code
for static problem.

Important python variables:
Collections:
   netan - (dict) ftbl structured content
   tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
   collection
   measures - (dict) exp data
   rAb - (list) reduced linear systems A*x_cumo=b by weight
   scale - unique scale names
   nrow - counts scale names
   o_sc - ordered scale names
   o_meas - ordered measurement types
org - (str) prefix of .ftbl  file like "PPP"
File names (str):
   n_ftbl (descriptor f_ftbl)
   n_R (R code) (f)
   n_fort (fortran code) (ff)
Counts:
   nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total), nb_ffn, nb_ffx (free fluxes)
Index translators:
   fwrv2i - flux names to index in R:fwrv
   cumo2i - cumomer names to index in R:x
   ir2isc - mapping measurement rows indexes on scale index isc[meas]=ir2isc[meas][ir]
Vector names:
   cumos (list) - names of R:x
   o_mcumos - cumomers involved in measurements

Important R variables:
Scalars:
   nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
   nb_ffn, nb_ffx, nb_ff (free fluxes),
   nb_fcn, nb_fcx, nb_fc (constrained fluxes),
   nb_ineq, nb_param, nb_fmn
Name vectors:
   nm_cumo, nm_fwrv, nm_fallnx, nm_fln, nm_flx, nm_fl, nm_par,
   nm_ffn, nm_ffx,
   nm_fcn, nm_fcx,
   nm_mcumo, nm_fmn
Numeric vectors:
   fwrv - all fluxes (fwd+rev)
   x - all cumomers (weight1+weight2+...)
   param - free flux net, free flux xch, scale label, scale mass, scale peak
   fcn, fcx, fc,
   bp - helps to construct the rhs of flux system
   xi -cumomer input vector
   fallnx - complete flux vector (constr+net+xch)
   bc - helps to construct fallnx
   li - inequality vector (mi%*%fallnx>=li)
   ir2isc - measur row to scale vector replicator
   ci - inequalities for param use (ui%*%param-ci>=0)
   measvecti,
   measinvvar,
   imeas,
   fmn
Matrices:
   Afl, qrAfl, invAfl,
   p2bfl - helps to construct the rhs of flux system
   mf, md - help to construct fallnx
   mi - inequality matrix (ftbl content)
   ui - inequality matrix (ready for param use)
   measmat - measmat*(x[imeas];1)=vec of simulated not-yet-pooled and not-yet-scaled measurements
Functions:
   param2fl_usm_eul - translate param to flux and cumomer vector (initial approximation)
   icumo_cost - cost function (khi2)
   cumo_grad - finite difference gradient
   cumo_gradj - implicit derivative gradient
"""

# 2012-03-06 sokol: initial version (derived from ftbl2optR.py)

if __name__ == "__main__":
    import sys
    import os
    import stat
    import time
    import copy
    import getopt

    me=os.path.realpath(sys.argv[0])
    sys.path.append(os.path.dirname(me))
    me=os.path.basename(me)

    from tools_ssg import *
    import C13_ftbl

    try:
        import psyco
        psyco.full()
        #psyco.log()
        #psyco.profile()
    except ImportError:
        #print 'Psyco not installed, the program will just run slower'
        pass

    def usage():
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--DEBUG] network_name[.ftbl]\n")

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "DEBUG"])
    except getopt.GetoptError, err:
        #pass
        sys.stderr.write(str(err)+"\n")
        usage()
        #sys.exit(1)

    fullsys=False
    DEBUG=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--cost":
            cost=True
        elif o=="--fullsys":
            fullsys=True
        elif o=="--DEBUG":
            DEBUG=True
        else:
            #assert False, "unhandled option"
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            pass
    #aff("args", args);##
    if len(args) != 1:
        usage()
        exit(1)
    org=os.path.basename(args[0])
    dirorg=os.path.dirname(args[0]) or '.'

    # cut .ftbl if any
    if org[-5:]==".ftbl":
        org=org[:-5]
    fullorg=os.path.join(dirorg, org)

    #-->
    #DEBUG=True
    import ftbl2code
    ftbl2code.DEBUG=DEBUG
    #org="ex3"
    #org="PPP_exact"
    #DEBUG=True
    if DEBUG:
        import pdb


    n_ftbl=fullorg+".ftbl"
    n_R=fullorg+".R"
    #n_fort=fullorg+".f"
    f_ftbl=open(n_ftbl, "r")
    try:
        os.chmod(n_R, stat.S_IWRITE)
    except:
        pass
    f=open(n_R, "w")

    # parse ftbl
    try:
        ftbl=C13_ftbl.ftbl_parse(f_ftbl)
    except:
        print sys.exc_info()[1]
        sys.exit(1)
    f_ftbl.close()

    # analyse network
    # reload(C13_ftbl)

    netan=C13_ftbl.ftbl_netan(ftbl)

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys)

#    f.write(
    """
# output flux repartition
cat("Dependent fluxes:\\n")
if (nb_fln) {
   print(paste(nm_fln,"net",sep="_"))
}
if (nb_flx) {
   print(paste(nm_flx,"xch",sep="_"))
}
cat("Free fluxes:\\n")
if (nb_ffn) {
   print(paste(nm_ffn,"net",sep="_"))
}
if (nb_ffx) {
   print(paste(nm_ffx,"xch",sep="_"))
}
cat("Constrained fluxes:\\n")
if (nb_fcn) {
   print(paste(nm_fcn,"net",sep="_"))
}
if (nb_fcx) {
   print(paste(nm_fcx,"xch",sep="_"))
}
"""
    f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\\n", sep="")
}
#browser()


# prepare flux index conversion
ifwrv=1:nb_fwrv
names(ifwrv)=nm_fwrv
ifl_in_fw=ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")]
iff_in_fw=ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")]
# index couples for jacobian df_dfl, df_dffd
cfw_fl=crv_fl=cbind(ifl_in_fw, 1:nb_fl)
cfw_ff=crv_ff=cbind(iff_in_fw, 1:nb_ff)
crv_fl[,1]=(nb_fwrv/2)+crv_fl[,1]
crv_ff[,1]=(nb_fwrv/2)+crv_ff[,1]
""")

    f.write("""

## variables for isotopomer cinetics
tstart=0.;
tmax=%(tmax)f;
dt=%(dt)f;

# metabolite pools are : all (poolall) which is divided in free (poolf) and
# constrained (poolc)

# constrained pool
poolc=c(%(poolc)s)
nm_poolc=c(%(nm_poolc)s)
names(poolc)=nm_poolc

# starting values for free pool
poolf=c(%(poolf)s)
nm_poolf=c(%(nm_poolf)s)
names(poolf)=nm_poolf
nb_poolf=length(poolf)
nb_f$nb_poolf=nb_poolf

nm_poolall=c(nm_poolf, nm_poolc)
poolall=c(poolf, poolc)
names(poolall)=nm_poolall

# extend param vector by free pools
if (nb_poolf > 0) {
   param=c(param, log(poolf))
   nm_par=c(nm_par, nm_poolf)
}
nm_list$par=nm_par

nm_list$poolf=nm_poolf
nm_list$poolc=nm_poolc
nm_list$poolall=nm_poolall

#browser()
if (nb_poolf > 0) {
   # extend inequalities ui, ci by log(cupp) >= poolf >= log(clowp)
   nb_row=nrow(ui)
   nb_col=ncol(ui)
   ui=cbind(ui, matrix (0., nrow=nrow(ui), ncol=nb_poolf))
   ui=rbind(ui, matrix (0., nrow=2*nb_poolf, ncol=ncol(ui)))
   ui[nb_row+1:nb_poolf,nb_col+1:nb_poolf]=diag(1., nb_poolf)
   ui[nb_row+nb_poolf+1:nb_poolf,nb_col+1:nb_poolf]=diag(-1., nb_poolf)
   ci=c(ci, rep(log(clowp), nb_poolf))
   ci=c(ci, rep(-log(cupp), nb_poolf))
   nm_i=c(nm_i, paste(nm_poolf, ">=", clowp, sep=""))
   nm_i=c(nm_i, paste(nm_poolf, "<=", cupp, sep=""))
   rownames(ui)=nm_i
   names(ci)=nm_i
   colnames(ui)=nm_par
}

# prepare mapping of metab pools on cumomers
nminvm=matrix(unlist(strsplit(nm_rcumo, ":")), ncol=2, byrow=T)[,1]
ipc2ix=match(paste("pc", nminvm, sep=":"), nm_poolall, nomatch=0)
ipf2ix=match(paste("pf", nminvm, sep=":"), nm_poolall, nomatch=0)
ip2ix=ipc2ix+ipf2ix
nb_f$ip2ix=ip2ix
nb_f$ipf2ix=ipf2ix
#browser()
# prepare mapping of free pools on vector xw
if (nb_f$nb_poolf > 0) {
   nm_px=nm_poolall[ip2ix]
   nbc_cumos=c(0., cumsum(nb_rcumos))
   iparpf2ix=lapply(1:nb_rw, function(iw){
      ipaire=matrix(0, nrow=0, ncol=2)
      lapply(1:nb_poolf, function(ipf) {
         i=which(nm_px[nbc_cumos[iw]+(1:nb_rcumos[iw])]==nm_poolf[ipf])
         if (length(i)) {
            ipaire <<- rbind(ipaire, cbind(i, ipf))
         }
      })
      return(ipaire)
   })
   nb_f$iparpf2ix=iparpf2ix
}
# read measvecti from a file specified in ftbl
flabcin="%(flabcin)s"
if (nchar(flabcin)) {
   measvecti=as.matrix(read.table(flabcin, header=T, row.names=1, sep="\t", check=F, comment=""))
   nm_row=rownames(measvecti)
   # put in the same row order as simulated measurments
   # check if nm_meas are all in rownames
   if (all(nm_meas %%in%% nm_row)) {
      measvecti=measvecti[nm_meas,,drop=F]
   } else {
      # try to strip row number from measure id
      nm_strip=sapply(strsplit(nm_meas, ":"), function(v) {
         v[length(v)]="";
         paste(v, sep="", collapse=":")
      })
      im=pmatch(nm_strip, nm_row)
      ina=is.na(im)
      if (any(ina)) {
         stop(paste("Cannot match the following measurement(s) in the file '", flabcin, "':\\n", paste(nm_meas[ina], sep="", collapse="\\n"), sep="", collapse=""))
      }
      measvecti=measvecti[im,,drop=F]
      #stopifnot(all(!is.na(measvecti)))
      stopifnot(typeof(measvecti)=="double")
   }
   ti=c(tstart, as.numeric(colnames(measvecti)))
   i=which(ti<=tmax)
   ti=ti[i]
   measvecti=measvecti[,i[-1]-1,drop=F]
   stopifnot(all(!is.na(ti)))
   nb_ti=length(ti)
   # divide the first time step in ndiv intervals with 2**x growing of interval length
   ndiv=4
   tmp=cumsum(2**(1:ndiv))
   tifull=c(tstart, (ti[2]/tmp[ndiv])*tmp, ti[-(1:2)])
   
   # divide each interval in 2.
   tifull=c(rbind(head(tifull, -1), head(tifull, -1)+diff(tifull)/2), tail(tifull, 1))
   
   # once again divide each interval in 2.
   #tifull=c(rbind(head(tifull, -1), head(tifull, -1)+diff(tifull)/2), tail(tifull, 1))

   #tifull=head(rep(ti, each=ndiv), -(ndiv-1))
   #for (idiv in 1:(ndiv-1)) {
   #   i=idiv+1+(0:(nb_ti-2))*ndiv
   #   th=idiv/ndiv
   #   tifull[i]=(1.-th)*head(ti, -1)+th*tail(ti, -1)
   #}
} else {
   measvecti=NULL
   ti=seq(tstart, tmax, by=dt)
   tifull=ti
   if (optimize) {
      warning("A fitting is requested but no labeling data are provided by 'file_labcin' option in the ftbl file.
The fitting is ignored as if '--noopt' option were asked.")
      optimize=F
   }
}
nb_ti=length(ti)

# formated output in kvh file
fkvh=file("%(fullorg)s_res.kvh", "w")
"""%{
    "fullorg": escape(fullorg, "\\"),
    "poolf": join(", ", (-p for p in netan["met_pools"].values() if p < 0.)),
    "nm_poolf": join(", ", (n for (n,p) in netan["met_pools"].iteritems() if p < 0.), '"pf:', '"'),
    "poolc": join(", ", (p for p in netan["met_pools"].values() if p > 0.)),
    "nm_poolc": join(", ", (n for (n,p) in netan["met_pools"].iteritems() if p > 0.), '"pc:', '"'),
    "dt": netan["opt"]["dt"],
    "tmax": netan["opt"]["tmax"],
    "flabcin": netan["opt"].get("file_labcin", ""),
})
    # main part: call optimization
    f.write("""
#browser()
names(param)=nm_par
nb_param=length(param)
if (initrand) {
   param[]=runif(nb_param);
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx
}
if (nb_sc && !is.null(measvecti)) {
   # set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
   # for corresponding measurements
   # cjac=F because param is not complete here, it lacks scaling params
   vr=icumo_resid(param, cjac=F, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
   if (!is.null(vr$err) && vr$err) {
      stop(vr$mes)
   }
   ##save(vr, ti, file="vr.Rdata")
   ##stop("aha")
   # unscaled simulated measurements (usm) [imeas, itime]
   #browser()
   inna=which(!is.na(measvecti))
   simvec=vr$usm
   ms=(measvecti*simvec*measinvvar)[inna]
   ss=(simvec*simvec*measinvvar)[inna]
   for (i in nb_ff+1:nb_sc) {
      im=outer(ir2isc==(i+1), rep(T, nb_ti), "&")[inna]
      param[i]=sum(ms[im])/sum(ss[im])
   }
} else if (nb_sc > 0) {
   # we dont have measurements yet, just set all scalings to 1.
   param[nb_ff+1:nb_sc]=1.
}
# check if initial approximation is feasible
ineq=ui%*%param-ci
param_old=param
if (any(ineq <= -1.e-10)) {
   cat("The following ", sum(ineq<= -1.e-10), " ineqalities are not respected at starting point:\\n", sep="")
   print(ineq[ineq<= -1.e-10,1])
   # put them inside
   param=put_inside(param, ui, ci)
   if (any(is.na(param))) {
      if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
         # fatal error occurd
         stop(paste("put_inside: ", attr(param, "mes"), collapse=""))
      }
   } else if (!is.null(attr(param, "err")) && attr(param, "err")==0){
      # non fatal problem
      warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
   }
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx
}

# see if there are any active inequalities at starting point
ineq=ui%*%param-ci
if (any(abs(ineq)<=1.e-10)) {
   cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point:\\n",
      paste(names(ineq[abs(ineq)<=1.e-10,1]), collapse="\\n"), "\\n", sep="")
}
cat("influx_i\\n", file=fkvh)
cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="")
cat("\\tdate of run\\t", date(), "\\n", file=fkvh, sep="")
# save options of command line
cat("\\tR command line\\n", file=fkvh)
obj2kvh(opts, "opts", fkvh, indent=2)

# resume system sizes
obj2kvh(nb_sys, "system sizes", fkvh)

# save initial fluxes
cat("starting point\\n", file=fkvh)
obj2kvh(param, "free parameters", fkvh, indent=1)

fwrv=jx_f$fwrv
obj2kvh(fwrv, "fwd-rev flux vector", fkvh, indent=1)

f=jx_f$fallnx
obj2kvh(f, "net-xch01 flux vector", fkvh, indent=1)

obj2kvh(ui%*%param-ci, "inequality slacks (must be >=0)", fkvh, indent=1)

#browser() # before the first call to param2fl_usm_eul

# starting cost value
if (!is.null(measvecti)) {
   vr=icumo_resid(param, cjac=F, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
   rcost=icumo_cost(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr, poolall, ti, tifull)
   obj2kvh(rcost, "starting cost value", fkvh, indent=1)
}

obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1)
btmp=c(p2bfl%*%head(param, nb_f$nb_ff)+bp)
names(btmp)=dimnames(Afl)[[1]]
obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1)
obj2kvh(measvec, "measurement set", fkvh, indent=1)

names(param)=nm_par
""")
    f.write("""
control_ftbl=list(%(ctrl_ftbl)s)
"""%{
    "ctrl_ftbl": join(", ", (k[8:]+"="+str(v) for (k,v) in netan["opt"].iteritems() if k.startswith("optctrl_"))),
})
    f.write("""
opt_wrapper=function(measvecti, fmn, ctrace=1) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=ctrace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, icumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, ipooled, measvecti, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr, poolall, ti, tifull)
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=ctrace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, icumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, ipooled, measvecti, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr, poolall, ti, tifull)
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=ctrace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, icumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, ipooled, measvecti, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr, poolall, ti, tifull)
   } else if (method == "nlsic") {
      control=list(trace=ctrace, btdesc=0.75, maxit=50, errx=1.e-5,
         ci=list(report=F))
      control[names(control_ftbl)]=control_ftbl
      res=nlsic(param, icumo_resid, 
         ui, ci, control, e=NULL, eco=NULL, flsi=lsi_fun,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar,
         spAbr, poolall, ti, tifull)
      if (res$err || is.null(res$par)) {
         # store res in kvh
         obj2kvh(res, "optimization aborted here", fkvh)
         warning(res$mes)
         res$par=rep(NA, length(param))
         res$cost=NA
      } else if (nchar(res$mes)) {
         warning(res$mes)
      }
   } else if (method == "ipopt") {
      control=list(max_iter=500, print_level=ctrace*5)
      control[names(control_ftbl)]=control_ftbl
      tui=c(t(ui))
      eval_g=function(x, nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, ipooled=ipooled, measvec=measvecti, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr, pool=poolall, ti=ti, tifull=tifull) {
         return(ui%*%x)
      }
      eval_jac_g=function(x, nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, ipooled=ipooled, measvec=measvecti, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr, pool=poolall, ti=ti, tifull=tifull) {
         return(tui)
      }
      ui_row_spars=rep.int(1, ncol(ui))
      res=ipoptr(param, icumo_cost, cumo_gradj,
         lb=NULL, ub=NULL,
         eval_g=eval_g,
         eval_jac_g=eval_jac_g,
         eval_jac_g_structure=lapply(1:nrow(ui), function(i)ui_row_spars),
         constraint_lb=ci,
         constraint_ub=rep(1.e19, length(ci)),
         eval_h=NULL,
         eval_h_structure=NULL,
         opts=control,
         ipoptr_environment=new.env(),
         nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, ipooled=ipooled, measvec=measvecti, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr, pool=poolall, ti=ti, tifull=tifull)
      res$par=res$solution
      names(res$par)=nm_par
      if (res$status != 0) {
         # store res in kvh
         obj2kvh(res$par, "optimization aborted here", fkvh)
         warning(res$message)
      }
   } else {
      stop(paste("Unknown minimization method '", method, "'", sep=""))
   }
   return(res)
}
if (optimize) {
   # check if at starting position all fluxes can be resolved
#browser()
   if (is.null(jx_f$uujac)) {
      # calculate jacobian here
      vr=icumo_resid(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
      if (!is.null(vr$err) && vr$err) {
         stop(vr$mes)
      }
   }
   qrj=qr(jx_f$dr_dff)
   d=diag(qrj$qr)
   qrj$rank=sum(abs(d)>abs(d[1])*1.e-14)
   if (qrj$rank) {
      nm_uns=nm_ff[qrj$pivot[-(1:qrj$rank)]]
   } else {
      nm_uns=nm_ff
   }
   if (qrj$rank < nb_ff && !(least_norm || method!="nlsic")) {
      # Too bad. The jacobian of free fluxes is not of full rank.
      dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx)
      write.matrix(formatC(jx_f$dr_dff, 15), file="dbg_dr_dff_singular.txt", sep="\t")
      stop(paste("Provided measurements (isotopomers and fluxes) are not
sufficient to resolve all free fluxes.
Unsolvable fluxes may be:
", paste(nm_uns, sep=", ", collapse=", "),
         "\nJacobian dr_dff is dumped in dbg_dr_dff_singular.txt", sep=""))
   }
   if (TIMEIT) {
      cat("optim   : ", date(), "\\n", sep="")
   }
#browser()
   # zero crossing strategy
   # inequalities to keep sens of net flux on first call to opt_wrapper()
   # if active after the first optimization, they are inversed on the second
   # call to opt_wrapper()
   mi_zc=NULL
   li_zc=NULL
   if (nb_fn && zerocross) {
      # add lower limits on [df].net >= zc for positive net fluxes
      # and upper limits on [df].net <= -zc for negative net fluxes
      nm_izc=c()
      ipos=names(which(fallnx[grep("^[df]\\\\.n\\\\.", nm_fallnx)]>=0.))
      ineg=names(which(fallnx[grep("^[df]\\\\.n\\\\.", nm_fallnx)]<0.))
      mi_zc=matrix(0, nrow=length(ipos)+length(ineg), ncol=nb_fallnx)
      colnames(mi_zc)=nm_fallnx
      if (length(ipos)) {
         nm_izc=c(nm_izc, paste("zc ", ipos, ">=", zc, sep=""))
         mi_zc[(1:length(ipos)),ipos]=diag(1., length(ipos))
      }
      if (length(ineg)) {
         nm_izc=c(nm_izc, paste("zc ", ineg, "<=", -zc, sep=""))
         mi_zc[length(ipos)+(1:length(ineg)),ineg]=diag(-1., length(ineg))
      }
      rownames(mi_zc)=nm_izc
      li_zc=rep(zc, length(nm_izc))
      ui_zc=cbind(mi_zc%*%(md%*%invAfl%*%p2bfl+mf),
         matrix(0., nrow=nrow(mi_zc), ncol=nb_sc+nb_poolf))
      ci_zc=li_zc-mi_zc%*%mic
      # remove all zero rows in ui_zc (constrained fluxes with fixed values)
      # find zero indexes
      zi=apply(ui_zc,1,function(v){return(max(abs(v))<=1.e-14)});
      ui_zc=ui_zc[!zi,,drop=F];
      mi_zc=mi_zc[!zi,,drop=F];
      ci_zc=ci_zc[!zi];
      nm_izc=nm_izc[!zi];
      # remove redundant/contradictory inequalities
      nb_zc=nrow(ui_zc)
      nb_i=nrow(ui)
      ired=c()
      if (nb_zc > 0) {
         for (i in 1:nb_zc) {
            nmqry=nm_izc[i]
            for (j in 1:nb_i) {
               if ((diff(range(ui[j,]-ui_zc[i,])) < 1.e-10 ||
                  diff(range(ui[j,]+ui_zc[i,])) < 1.e-10) &&
                  abs(ci[j]-ci_zc[i])<1.e-2) {
#browser()
                  # redundancy
                  cat("inequality '", nmqry, "' redundant or contradictory with '", nm_i[j], "' is removed.\\n", sep="")
                  ired=c(ired, i)
                  break
               }
            }
         }
      }
      if (!is.null(ired)) {
         # remove all ired inequalities
         ui_zc=ui_zc[-ired,,drop=F]
         ci_zc=ci_zc[-ired]
         nm_izc=nm_izc[-ired]
         mi_zc=mi_zc[-ired,,drop=F]
      }
      if (nrow(ui_zc)) {
         # add zc inequalities
         ui=rbind(ui, ui_zc)
         ci=c(ci, ci_zc)
         nm_i=c(nm_i, nm_izc)
         mi=rbind(mi, mi_zc)
      }
   #print(ui)
   #print(ci)
   }
   # pass control to the chosen method
   res=opt_wrapper(measvecti, fmn)
   if (any(is.na(res$par))) {
      obj2kvh(res, "optimization process informations", fkvh)
      stop(res$mes)
   }
#browser()
   if (zerocross && !is.null(mi_zc)) {
      # inverse active "zc" inequalities
      nm_inv=names(which((ui%*%res$par-ci)[,1]<=1.e-10))
      i=grep("^zc ", nm_inv, v=T)
      if (length(i) > 0) {
         reopt=TRUE
         i=str2ind(i, nm_i)
         cat("The following inequalities are active after the first stage
of zero crossing strategy and will be inverted:\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="")
         ipos=grep(">=", nm_i[i], fix=T, v=T)
         ineg=grep("<=", nm_i[i], fix=T, v=T)
         ui[i,]=-ui[i,]
         if (length(ipos)) {
            ipzc=str2ind(ipos, nm_izc)
            ipos=str2ind(ipos, nm_i)
            ci[ipos]=zc+mi_zc[ipzc,,drop=F]%*%mic
            nm_i[ipos]=sub(">=", "<=-", nm_i[ipos])
         }
         if (length(ineg)) {
            inzc=str2ind(ineg, nm_izc)
            ineg=str2ind(ineg, nm_i)
            ci[ineg]=zc+mi_zc[inzc,,drop=F]%*%mic
            nm_i[ineg]=sub("<=-", ">=", nm_i[ineg])
         }
         rownames(ui)=nm_i
         names(ci)=nm_i
      } else {
         reopt=FALSE
      }
      # enforce new inequalities
      if (reopt) {
         # reoptimize
         param=put_inside(res$par, ui, ci)
         if (any(is.na(param))) {
            if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
               # fatal error occurd, don't reoptimize
               warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
               param=res$param
               reopt=FALSE
            }
         } else if (!is.null(attr(param, "err")) && attr(param, "err")==0){
            # non fatal problem
            warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
         }
         res=opt_wrapper(measvecti, fmn)
         if (any(is.na(res$par))) {
            obj2kvh(res, "optimization process informations", fkvh)
            stop("Second optimization (after zero crossing) failed.")
         }
      }
   }
   param=res$par
   names(param)=nm_par
   res$jacobian=NULL # sometimes it's too huge to be presented to humans
   obj2kvh(res, "optimization process informations", fkvh)
}
if (TIMEIT) {
   cat("postopt : ", date(), "\\n", sep="")
}
# active constraints
ine=abs(ui%*%param-ci)
ine=ine < 1.e-10 | ine < abs(ci)*0.01
if (any(ine)) {
   obj2kvh(nm_i[ine], "active (up to 1%) inequality constraints", fkvh)
}
#browser()
if (!is.null(measvecti)) {
   if (is.null(jx_f$jacobian)) {
      rres=icumo_resid(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
   } # else use the last calculated jacobian
   rcost=norm2(jx_f$res)
   obj2kvh(rcost, "final cost", fkvh)
   obj2kvh(jx_f$ureslab, "simulated-measured labeling", fkvh)
   obj2kvh(jx_f$reslab, "(simulated-measured)/sd_exp labeling", fkvh)
   if (nb_fmn > 0) {
      obj2kvh(jx_f$uresflu, "simulated-measured fluxes", fkvh)
      obj2kvh(jx_f$resflu, "(simulated-measured)/sd_exp fluxes", fkvh)
   }
   # gradient -> kvh
   gr=2*c(jx_f$res%tmm%jx_f$jacobian)
   names(gr)=nm_par
   obj2kvh(gr, "gradient vector", fkvh)
   #obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh)
} else {
   if (is.null(jx_f$usm)) {
      # simulate measures
      v=param2fl_usm_eul(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, spAbr, poolall, ti, measmat, irmeas, ipooled, tifull)
   }
}

if (fullsys) {
   nm_flist=nm_list
   nm_flist$rcumo=nm_cumo
   v=param2fl_usm_eul(param, cjac=T, nb_f, nm_flist, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAbr_f, poolall, ti, measmat, imeas, ipooled, tifull)
} # else use the last calculated usm

# simulated measurements -> kvh
obj2kvh(jx_f$usm, "simulated unscaled labeling measurements", fkvh)
if (nb_sc > 0) {
   obj2kvh(jx_f$usm*c(1.,param)[ir2isc], "simulated scaled labeling measurements", fkvh)
} else {
   obj2kvh(jx_f$usm, "simulated scaled labeling measurements", fkvh)
}
if (nb_fmn) {
   obj2kvh(cbind(value=jx_f$fallnx[nm_fmn], sd=1./sqrt(invfmnvar)), "simulated flux measurements", fkvh)
}
# simulated cumomers -> kvh
obj2kvh(jx_f$xsim, "simulated cumomers", fkvh)
# simulated derivatives of cumomers -> kvh
obj2kvh(jx_f$xpsim, "simulated derivatives of cumomers", fkvh)

x=jx_f$x
obj2kvh(cumo2mass(x), "MID vector", fkvh)
# replace :i by :x1's
nm_mask=unlist(lapply(strsplit(names(x), ":"), function(v) {
   n=length(v)
   metab=paste(v[-n], collapse=":")
   mask=int2bit(v[n], len=clen[metab])
   c(metab, mask)
}))
nm_mask=t(matrix(nm_mask, nrow=2))
o=order(nm_mask[,1], nm_mask[,2])
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T)
nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#")
names(nm_mask)=names(x)
names(x)=nm_mask
obj2kvh(x[o], "cumomer vector", fkvh)

fwrv=jx_f$fwrv
names(fwrv)=nm_fwrv
#obj2kvh(fwrv, "fwd-rev flux vector", fkvh)

fallnx=jx_f$fallnx
names(fallnx)=nm_fallnx

flnx=jx_f$flnx
names(flnx)=c(nm_fl)

if (sensitive=="mc") {
   if (TIMEIT) {
      cat("monte-ca: ", date(), "\\n", sep="")
   }
   # Monte-Carlo simulation in parallel way
   mc_inst=library(multicore, warn.conflicts=F, verbose=F, logical.return=T)
   invar=c(rep(measinvvar, nb_ti-1), invfmnvar)
   simcumom=jx_f$simvec
   simfmn=fallnx[nm_fmn]
   mc_sim=function(i) {
      # random measurement generation
      if (nb_meas) {
         meas_mc=matrix(rnorm(nb_meas*nb_ti, simcumom, 1./sqrt(measinvvar)), nrow=nrow(simcumom))
      } else {
         meas_mc=c()
      }
      if (nb_fmn) {
         fmn_mc=rnorm(nb_fmn, simfmn, 1./sqrt(invfmnvar))
      } else {
         fmn_mc=c()
      }
      #cat("imc=", i, "\\n", sep="")
      # minimization
      res=opt_wrapper(meas_mc, fmn_mc, trace=0)
      # return the solution
      return(list(cost=norm2(jx_f$res), par=res$par))
   }
   # parallel execution
   if (mc_inst) {
      mc_res=mclapply(1:nmc, mc_sim)
   } else {
      mc_res=lapply(1:nmc, mc_sim)
   }
   # strip failed mc iterations
   free_mc=sapply(mc_res, function(l) {if (class(l)=="character" || is.na(l$cost)) { ret=rep(NA, nb_param+1) } else { ret=c(l$cost, l$par) }; ret })
   if (length(free_mc)==0) {
      warning("Parallel exectution of Monte-Carlo simulations has failed.")
      free_mc=matrix(NA, nb_param+1, 0)
   }
   cost_mc=free_mc[1,]
   free_mc=free_mc[-1,,drop=F]
   # remove failed m-c iterations
   ifa=which(is.na(free_mc[1,]))
   if (length(ifa)) {
      if (ncol(free_mc)>length(ifa)) {
         warning("Some Monte-Carlo iterations failed.")
      }
      free_mc=free_mc[,-ifa,drop=F]
      cost_mc=cost_mc[-ifa]
      nmc=length(cost_mc)
   }
#browser()
   dimnames(free_mc)[[1]]=nm_par
   cat("monte-carlo\\n", file=fkvh)
   indent=1
   obj2kvh(nmc, "sample number", fkvh, indent)
   avaco=multicore:::detectCores()
   obj2kvh(avaco, "detected cores", fkvh, indent)
   avaco=max(1, avaco, na.rm=T)
   obj2kvh(min(avaco, options("cores")$cores, na.rm=T), "used cores", fkvh, indent)
   # cost section in kvh
   cat("\\tcost\\n", file=fkvh)
   indent=2
   obj2kvh(mean(cost_mc), "mean", fkvh, indent)
   obj2kvh(median(cost_mc), "median", fkvh, indent)
   obj2kvh(sd(cost_mc), "sd", fkvh, indent)
   obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent)
   obj2kvh(quantile(cost_mc, c(0.025, 0.975)), "ci", fkvh, indent)
   # free parameters section in kvh
   cat("\\tfree parameters\\n", file=fkvh)
   indent=2
   # param stats
   # mean
   obj2kvh(apply(free_mc, 1, mean), "mean", fkvh, indent)
   # median
   parmed=apply(free_mc, 1, median)
   obj2kvh(parmed, "median", fkvh, indent)
   # covariance matrix
   covmc=cov(t(free_mc))
   obj2kvh(covmc, "covariance", fkvh, indent)
   # sd
   sdmc=sqrt(diag(covmc))
   obj2kvh(sdmc, "sd", fkvh, indent)
   obj2kvh(sdmc*100/abs(param), "rsd (%)", fkvh, indent)
   # confidence intervals
   ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)))
   ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
   dimnames(ci_mc)[[2]][3]="length"
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
   obj2kvh((ci_mc-cbind(param, param, 0))*100/abs(param),
      "relative 95% confidence intervals (%)", fkvh, indent)
   
   # net-xch01 stats
   fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx)
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx
   if (length(fallnx_mc)) {
      dimnames(fallnx_mc)[[1]]=nm_fallnx
      # form a matrix output
      fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0)
      #cat("\\tall net-xch01 fluxes\\n", file=fkvh)
      # mean
#browser()
      fallout=cbind(fallout, mean=apply(fallnx_mc, 1, mean))
      #obj2kvh(apply(fallnx_mc, 1, mean), "mean", fkvh, indent)
      # median
      parmed=apply(fallnx_mc, 1, median)
      fallout=cbind(fallout, median=parmed)
      #obj2kvh(parmed, "median", fkvh, indent)
      # covariance matrix
      covmc=cov(t(fallnx_mc))
      dimnames(covmc)=list(nm_fallnx, nm_fallnx)
      #obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      fallout=cbind(fallout, sd=sdmc)
      #obj2kvh(sdmc, "sd", fkvh, indent)
      fallout=cbind(fallout, "rsd (%)"=sdmc*100/abs(fallnx))
      #obj2kvh(sdmc*100/abs(fallnx), "rsd (%)", fkvh, indent)
      # confidence intervals
      ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
      ci_mc=cbind(ci_mc, (ci_mc-cbind(fallnx, fallnx, 0))*100/abs(fallnx))
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)")
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
      fallout=cbind(fallout, ci_mc)
      #obj2kvh((ci_mc-cbind(fallnx, fallnx, 0))*100/abs(fallnx),
      #   "relative 95% confidence intervals (%)", fkvh, indent)
      o=order(nm_fallnx)
      obj2kvh(fallout[o,,drop=F], "all net-xch01 fluxes", fkvh, indent)
      obj2kvh(covmc[o,o], "covariance of all net-xch01 fluxes", fkvh, indent)
      
      # fwd-rev stats
      fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fwrv)
      dimnames(fwrv_mc)[[1]]=nm_fwrv
      fallout=matrix(0, nrow=nrow(fwrv_mc), ncol=0)
      #cat("\\tforward-reverse fluxes\\n", file=fkvh)
      # mean
      fallout=cbind(fallout, mean=apply(fwrv_mc, 1, mean))
      #obj2kvh(apply(fwrv_mc, 1, mean), "mean", fkvh, indent)
      # median
      parmed=apply(fwrv_mc, 1, median)
      fallout=cbind(fallout, median=parmed)
      #obj2kvh(parmed, "median", fkvh, indent)
      # covariance matrix
      covmc=cov(t(fwrv_mc))
      dimnames(covmc)=list(nm_fwrv, nm_fwrv)
      #obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      fallout=cbind(fallout, sd=sdmc)
      #obj2kvh(sdmc, "sd", fkvh, indent)
      fallout=cbind(fallout, "rsd (%)"=sdmc*100/abs(fwrv))
      #obj2kvh(sdmc*100/abs(fwrv), "rsd (%)", fkvh, indent)
      # confidence intervals
      ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
      ci_mc=cbind(ci_mc, (ci_mc-cbind(fwrv, fwrv, 0))*100/abs(fwrv))
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)")
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
      #obj2kvh((ci_mc-cbind(fwrv, fwrv, 0))*100/abs(fwrv),
      #   "relative 95% confidence intervals (%)", fkvh, indent)
      fallout=cbind(fallout, ci_mc)
      o=order(nm_fwrv)
      obj2kvh(fallout[o,,drop=F], "forward-reverse fluxes", fkvh, indent)
      obj2kvh(covmc[o,o], "covariance of forward-reverse fluxes", fkvh, indent)
   }
} else if (length(sensitive) && nchar(sensitive)) {
   warning(paste("Unknown sensitivity '", sensitive, "' method chosen.", sep=""))
}

if (TIMEIT) {
   cat("linstats: ", date(), "\\n", sep="")
}
# Linear method based on jacobian x_f
# reset fluxes and jacobians according to param
if (is.null(jx_f$jacobian)) {
   rres=icumo_resid(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
} # else use last calculated jacobian
# use last calculated jacobian
if (DEBUG) {
   library(numDeriv); # for numerical jacobian
   # numerical simulation
   rj=function(v, ...) { r=icumo_resid(v, cjac=F, ...); c(r$res) }
   #dr_dpn=jacobian(rj, param, method="simple", method.args=list(),
   #   nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr, poolall, ti, tifull)
   # to compare with jx_f$dr_dp
   gr=grad(icumo_cost, param, method="Richardson", method.args=list(), nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, ipooled, measvecti, ipooled, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr, poolall, ti, tifull)
#browser()
}

# covariance matrix of free fluxes
svj=svd(jx_f$jacobian)
i=svj$d/svj$d[1]<1.e-10
if (all(!i) && svj$d[1]<1.e-10) {
   # we could not find very small d, take just the last
   i[length(i)]=T
}
ibad=apply(svj$v[, i, drop=F], 2, which.contrib)
ibad=unique(unlist(ibad))
if (length(ibad) > 0) {
   warning(paste("Inverse of covariance matrix is numerically singular.\\nStatistically undefined parameter(s) seems to be:\\n",
      paste(nm_par[ibad], collapse="\\n"), "\\nFor more complete list, see sd columns in '/linear stats'\\nin the result file.", sep=""))
}
# "square root" of covariance matrix (to preserve numerical positive definitness)
rtcov=(svj$u)%*%(t(svj$v)/svj$d)
covpar=crossprod(rtcov)
if (nb_ff > 0) {
   covff=covpar[1:nb_ff,1:nb_ff]
   cat("linear stats\\n", file=fkvh)

   # sd free+dependent net-xch01 fluxes
   covfl=crossprod(rtcov[,1:nb_ff,drop=F]%mmt%(rbind(diag(1., nb_ff), dfl_dff)))
   nm_flfd=c(nm_ff, nm_fl)
   dimnames(covfl)=list(nm_flfd, nm_flfd)
   sdfl=sqrt(diag(covfl))
} else {
   sdfl=rep(0., nb_fl)
}
fl=c(head(param, nb_ff), flnx)
nm_flfd=c(nm_ff, nm_fl)
mtmp=cbind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl))
rownames(mtmp)=nm_flfd
o=order(nm_flfd)
obj2kvh(mtmp[o,,drop=F], "net-xch01 fluxes (sorted by name)", fkvh, indent=1)
if (nb_ff > 0) {
   obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1)
}

# sd of all fwd-rev
if (nb_ff > 0) {
   covf=crossprod(rtcov[,1:nb_ff,drop=F]%mmt%(jx_f$df_dff))
   sdf=sqrt(diag(covf))
} else {
   sdf=rep(0., length(fwrv))
}
mtmp=cbind(fwrv, sdf, sdf/abs(fwrv))
dimnames(mtmp)[[2]]=c("value", "sd", "rsd")
o=order(nm_fwrv)
obj2kvh(mtmp[o,], "fwd-rev fluxes (sorted by name)", fkvh, indent=1)
if (nb_ff > 0) {
   obj2kvh(covf, "covariance fwd-rev fluxes", fkvh, indent=1)
}

# pool -> kvh
sdpf=poolall
sdpf[]=0.

if (nb_poolf > 0) {
   # covariance matrix of free pools
   # "square root" of covariance matrix (to preserve numerical positive definitness)
   poolall[nm_poolf]=exp(param[nm_poolf])
   # cov wrt exp(pf)=pool
   covpf=covpar[nb_ff+nb_sc+1:nb_poolf,nb_ff+nb_sc+1:nb_poolf]/(poolall[nm_poolf]**2)
   dimnames(covpf)=list(nm_poolf, nm_poolf)
   sdpf[nm_poolf]=sqrt(diag(covpf))
}
poolall=poolall

mtmp=cbind("value"=poolall, "sd"=sdpf, "rsd"=sdpf/poolall)
rownames(mtmp)=nm_poolall
o=order(nm_poolall)
obj2kvh(mtmp[o,,drop=F], "metabolite pools (sorted by name)", fkvh, indent=1)
if (nb_poolf > 0) {
   o=order(nm_poolf)
   obj2kvh(covpf[o, o], "covariance free pools", fkvh, indent=1)
}
if (!is.null(measvecti)) {
   # goodness of fit (khi2 test)
   khi2test=list("khi2 value"=rcost, "data points"=length(jx_f$res),
      "fitted parameters"=nb_param, "degrees of freedom"=length(jx_f$res)-nb_param)
   khi2test$`khi2 reduced value`=khi2test$`khi2 value`/khi2test$`degrees of freedom`
   khi2test$`p-value, i.e. P(X^2<=value)`=pchisq(khi2test$`khi2 value`, df=khi2test$`degrees of freedom`)
   khi2test$conclusion=if (khi2test$`p-value, i.e. P(X^2<=value)` > 0.95) "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD" else "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD"
   obj2kvh(khi2test, "goodness of fit (khi2 test)", fkvh, indent=1)
}

if (prof) {
   Rprof(NULL)
}
close(fkvh)
""")

    f.write("""
# write edge.netflux property
fedge=file("%(d)s/edge.netflux.%(org)s", "w")
cat("netflux (class=Double)\\n", sep="", file=fedge)
nm_edge=names(edge2fl)
cat(paste(nm_edge, fallnx[edge2fl], sep=" = "), sep="\\n" , file=fedge)
close(fedge)

# write edge.xchflux property
fedge=file("%(d)s/edge.xchflux.%(org)s", "w")
flxch=paste(".x", substring(edge2fl, 4), sep="")
ifl=charmatch(flxch, substring(names(fallnx), 2))
cat("xchflux (class=Double)\\n", sep="", file=fedge)
cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
close(fedge)

# write node.log2pool property
fnode=file("%(d)s/node.log2pool.%(org)s", "w")
cat("log2pool (class=Double)\\n", sep="", file=fnode)
nm_node=substring(names(poolall), 4)
cat(paste(nm_node, log2(poolall), sep=" = "), sep="\\n" , file=fnode)
close(fnode)

if (TIMEIT) {
   cat("rend    : ", date(), "\\n", sep="")
}
"""%{
    "org": escape(org, "\\"),
    "d": escape(dirorg, "\\"),
})

    f.close()
    # make output files just readable to avoid later casual edition
    os.chmod(n_R, stat.S_IREAD)
