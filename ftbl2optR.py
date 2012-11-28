#!/usr/bin/env python

"""
Transform an ftbl to R code which will solve an optimization of flux analysis
problem min(S) over \Theta, where S=||Predicted-Observed||^2_\Sigma^2
and \Theta is a vector of free fluxes (net+xch), scaling parameters and metabolite
pools.
Predicted vector is obtained from cumomer vector x (calculated from
free fluxes and divided in chunks according to the cumo weight) by
multiplying it by the measurement matrices, weighted by metabolite
pools (in case of pooling) and scale factor, boths coming
from ftbl file. Observed values vector xo is extracted from ftbl file.
it is composed of flux and cumomer measurements.
\Sigma^2, covariance diagonal matrices sigma[flux|mass|label|peak]
are orginated from ftbl

usage: ./ftbl2optR.py organism
where organism is the ftbl informative part of file name
(before .ftbl), e.g. organism.ftbl
after execution a file organism.R will be created.
If it already exists, it will be silently overwritten.
The system Afl*flnx=bfl is created from ftbl file.

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
   measvec,
   fmn
Matrices:
   Afl, qrAfl, invAfl,
   p2bfl - helps to construct the rhs of flux system
   mf, md - help to construct fallnx
   mi - inequality matrix (ftbl content)
   ui - inequality matrix (ready for param use)
   measmat - for measmat*x+memaone=vec of simulated not-yet-scaled measurements
Functions:
   param2fl_x - translate param to flux and cumomer vector (initial approximation)
   cumo_cost - cost function (khi2)
   cumo_grad - finite difference gradient
   cumo_gradj - implicit derivative gradient
"""

# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used for sparse matrices.

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
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--emu] [--DEBUG] network_name[.ftbl]\n")

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "emu", "DEBUG"])
    except getopt.GetoptError, err:
        #pass
        sys.stderr.write(str(err)+"\n")
        usage()
        #sys.exit(1)

    fullsys=False
    DEBUG=False
    emu=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--cost":
            cost=True
        elif o=="--fullsys":
            fullsys=True
        elif o=="--emu":
            emu=True
        elif o=="--DEBUG":
            DEBUG=True
        else:
            #assert False, "unhandled option"
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            pass
    #aff("args", args);##
    #aff("opts", opts);##
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
    ftbl=C13_ftbl.ftbl_parse(f_ftbl)
    f_ftbl.close()

    # analyse network
    # reload(C13_ftbl)

    netan=C13_ftbl.ftbl_netan(ftbl)
    # prepare rcumo system
    if emu:
        rAb=C13_ftbl.rcumo_sys(netan, C13_ftbl.ms_frag_gath(netan))
    else:
        rAb=C13_ftbl.rcumo_sys(netan)

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys, emu)

    f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\\n", sep="")
}
#browser()
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
poolall=as.numeric(c(poolf, poolc))
names(poolall)=nm_poolall
pool=poolall

# extend param vector by free pools
if (nb_poolf > 0) {
   param=c(param, log(poolf))
   nm_par=c(nm_par, nm_poolf)
   nb_param=length(param)
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
   ui=cBind(ui, Matrix (0., nrow=nrow(ui), ncol=nb_poolf))
   ui=rBind(ui, Matrix (0., nrow=2*nb_poolf, ncol=ncol(ui)))
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
# prepare metabolite measurements
nb_poolm=%(nb_poolm)d
nb_f$nb_poolm=nb_poolm
nm_poolm=c(%(nm_poolm)s)
nm_list$poolm=nm_poolm

# measured values
vecpoolm=c(%(v_poolm)s)
names(vecpoolm)=nm_poolm

# inverse of variance for flux measurements
poolmdev=c(%(poolmdev)s)

# simulated metabolite measurements are calculated as
# measmatpool*poolall=>poolm
measmatpool=Matrix(0., nrow=nb_poolm, ncol=length(poolall))
dimnames(measmatpool)=list(nm_poolm, nm_poolall)
i=matrix(1+c(%(imeasmatpool)s), ncol=2, byrow=T)
measmatpool[i]=1.

# gather all measurement information
measurements=list(
   vec=list(labeled=measvec, flux=fmn, pool=vecpoolm),
   dev=list(labeled=measdev, flux=fmndev, pool=poolmdev),
   mat=list(labeled=measmat, flux=ifmn, pool=measmatpool),
   one=list(labeled=memaone)
)
"""%{
    "poolf": join(", ", (-netan["met_pools"][m] for m in netan["vpool"]["free"])),
    "nm_poolf": join(", ", netan["vpool"]["free"], '"pf:', '"'),
    "poolc": join(", ", (netan["met_pools"][m] for m in netan["vpool"]["constrained"])),
    "nm_poolc": join(", ", netan["vpool"]["constrained"], '"pc:', '"'),
    "nb_poolm": len(netan["metab_measured"]),
    "nm_poolm": join(", ", netan["metab_measured"].keys(), '"pm:', '"'),
    "v_poolm": join(", ", (item["val"] for item in netan["metab_measured"].values())),
    "poolmdev": join(", ", (item["dev"] for item in netan["metab_measured"].values())),
    "imeasmatpool": join(", ", valval((ir,netan["vpool"]["all2i"][m]) for (ir, ml) in enumerate(netan["metab_measured"].keys()) for m in ml.split("+"))),
})
    f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\\n", sep="")
}

# check if initial approximation is feasible
ineq=as.numeric(ui%*%param-ci)
names(ineq)=rownames(ui)
param_old=param
if (any(ineq <= -1.e-10)) {
   cat("The following ", sum(ineq<= -1.e-10), " ineqalities are not respected at starting point:\\n", sep="")
   print(ineq[ineq<= -1.e-10])
   # put them inside
   param=put_inside(param, ui, ci)
   if (any(is.na(param))) {
      if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
         # fatal error occured
         stop(paste("put_inside: ", attr(param, "mes"), collapse=""))
      }
   } else if (!is.null(attr(param, "err")) && attr(param, "err")==0) {
      # non fatal problem
      warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
   }
}

# zero crossing strategy
# inequalities to keep sens of net flux on first call to opt_wwrapper()
# if active they are removed on the second call to opt_wrapper()
fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx
mi_zc=NULL
li_zc=NULL
if (nb_fn && zerocross) {
   # add lower limits on [df].net >= zc for positive net fluxes
   # and upper limits on [df].net <= -zc for negative net fluxes
   nm_izc=c()
   ipos=names(which(fallnx[grep("[df].n.", nm_fallnx)]>=0.))
   ineg=names(which(fallnx[grep("[df].n.", nm_fallnx)]<0.))
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
   ui_zc=cBind(mi_zc%*%(md%*%invAfl%*%p2bfl+mf),
      Matrix(0., nrow=nrow(mi_zc), ncol=nb_sc+nb_poolf))
   ci_zc=li_zc-mi_zc%*%mic
   # remove redundant/contradictory inequalities
   nb_zc=nrow(ui_zc)
   nb_i=nrow(ui)
   ired=c()
   if (nb_zc > 0) {
      for (i in 1:nb_zc) {
         nmqry=nm_izc[i]
         for (j in 1:nb_i) {
            if ((max(abs(ui[j,]-ui_zc[i,]))<1.e-10 ||
               max(abs(ui[j,]+ui_zc[i,]))<1.e-10) &&
               abs(abs(ci[j])-abs(ci_zc[i]))<=1.e-2) {
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
      ui=rBind(ui, ui_zc)
      ci=c(ci, ci_zc)
      nm_i=c(nm_i, nm_izc)
      mi=rBind(mi, mi_zc)
   }
#print(ui)
#print(ci)
}

# set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
# for corresponding measurements
vr=param2fl_x(param, cjac=FALSE, nb_f, nm_list, nb_x, invAfl, p2bfl, g2bfl, bp, fc, xi, spa, emu, pool, measurements, ipooled)
if (!is.null(vr$err) && vr$err) {
   stop(vr$mes)
}
if (nb_sc > 0) {
   simvec=jx_f$usimcumom
   measinvvar=1./measurements$dev$labeled**2
   ms=measvec*simvec*measinvvar
   ss=simvec*simvec*measinvvar
   for (i in nb_ff+1:nb_sc) {
      im=(ir2isc==(i+1))
      param[i]=sum(ms[im])/sum(ss[im])
   }
}

# prepare flux index conversion
ifwrv=1:nb_fwrv
names(ifwrv)=nm_fwrv
ifl_in_fw=ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")]
iff_in_fw=ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")]
if (nb_fgr > 0) {
   ifg_in_fw=ifwrv[paste("fwd", substring(nm_fgr, 4), sep="")]
} else {
   ifg_in_fw=numeric(0)
}
# index couples for jacobian df_dfl, df_dffd
cfw_fl=crv_fl=cBind(ifl_in_fw, 1:nb_fl)
cfw_ff=crv_ff=cBind(iff_in_fw, 1:nb_ff)
cfw_fg=crv_fg=cBind(ifg_in_fw, tail(1:(nb_ff+nb_fgr), nb_fgr))
crv_fl[,1]=(nb_fwrv/2)+crv_fl[,1]
crv_ff[,1]=(nb_fwrv/2)+crv_ff[,1]
crv_fg[,1]=(nb_fwrv/2)+crv_fg[,1]

# see if there are any active inequalities at starting point
ineq=as.numeric(ui%*%param-ci)
names(ineq)=rownames(ui)
if (any(abs(ineq)<=1.e-10)) {
   cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point:\\n",
      paste(names(ineq[abs(ineq)<=1.e-10]), collapse="\\n"), "\\n", sep="")
}
""")

    f.write("""
if (TIMEIT) {
   cat("kvh init: ", date(), "\\n", sep="")
}

# formated output in kvh file
fkvh=file("%(fullorg)s_res.kvh", "w")
"""%{
    "fullorg": escape(fullorg, "\\"),
})
    # main part: call optimization
    f.write("""
cat("influx\\n", file=fkvh)
cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="")
# save options of command line
cat("\\tR command line\\n", file=fkvh)
obj2kvh(opts, "opts", fkvh, indent=2)

# resume system sizes
obj2kvh(nb_sys, "system sizes", fkvh)

# save initial flux and cumomer distribution
cat("starting point\\n", file=fkvh)
names(param)=nm_par
obj2kvh(param, "starting free parameters", fkvh, indent=1)

x=vr$x
names(x)=nm_x
obj2kvh(cumo2mass(x), "starting MID vector", fkvh, indent=1)
# replace :i by #x1's
nm_mask=t(sapply(strsplit(names(x), "[:+]"), function(v) {
   n=length(v)
   metab=v[1]
   mask=int2bit(v[2], len=clen[metab])
   c(metab, mask, if (n==3) v[3] else NULL)
}))
o=order(nm_mask[,1], nm_mask[,2])
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T)
tmp=paste(nm_mask[,1], nm_mask[,2], sep="#")
if (emu) {
   nm_mask=paste(tmp, nm_mask[,3], sep="+")
} else {
   nm_mask=tmp
}
names(nm_mask)=names(x)
names(x)=nm_mask
obj2kvh(x[o], "starting cumomer vector", fkvh, indent=1)

fwrv=vr$fwrv
n=length(fwrv)
names(fwrv)=nm_fwrv
obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, indent=1)

f=vr$fallnx
n=length(f)
names(f)=nm_fallnx
obj2kvh(f, "starting net-xch01 flux vector", fkvh, indent=1)

rres=cumo_resid(param, cjac=TRUE, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
names(rres$res)=c(nm_meas, nm_fmn, nm_poolm)
o=order(names(rres$res))
obj2kvh(rres$res[o], "starting residuals (simulated-measured)/sd_exp", fkvh, indent=1)

rcost=cumo_cost(param, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurments, ir2isc, spa, emu, pool, ipooled)
obj2kvh(rcost, "starting cost value", fkvh, indent=1)

obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1)
btmp=as.numeric(p2bfl%*%param[1:nb_f$nb_ff]+bp)
names(btmp)=dimnames(Afl)[[1]]
obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1)
names(measvec)=nm_meas
obj2kvh(measvec, "labeled measurement vector", fkvh, indent=1)

#cat("mass vector:\\n")
#print_mass(x)

names(param)=nm_par
""")
    f.write("""
control_ftbl=list(%(ctrl_ftbl)s)
"""%{
    "ctrl_ftbl": join(", ", (k[8:]+"="+str(v) for (k,v) in netan["opt"].iteritems() if k.startswith("optctrl_"))),
})
    f.write("""
opt_wrapper=function(measurements, trace=1) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "nlsic") {
      control=list(trace=trace, btfrac=0.25, btdesc=0.75, maxit=50, errx=1.e-5,
         ci=list(report=F), history=FALSE, adaptbt=TRUE)
      control[names(control_ftbl)]=control_ftbl
      res=nlsic(param, cumo_resid, 
         ui, ci, control, e=NULL, eco=NULL, flsi=lsi_fun,
         nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc,
         spa, emu, pool, ipooled)
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
      control=list(max_iter=500, print_level=trace*5)
      control[names(control_ftbl)]=control_ftbl
      tui=c(t(ui))
      eval_g=function(x, nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa) {
         return(ui%*%x)
      }
      eval_jac_g=function(x, nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa) {
         return(tui)
      }
      ui_row_spars=rep.int(1, ncol(ui))
      res=ipoptr(param, cumo_cost, cumo_gradj,
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
         nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa, emu=emu, pool=pool, ipooled=ipooled)
      res$par=res$solution
      names(res$par)=nm_par
      if(res$status != 0) {
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
   qrj=qr(jx_f$dr_dff, LAPACK=T)
   d=diag(qrj$qr)
   qrj$rank=sum(abs(d)>abs(d[1])*1.e-14)
   if (qrj$rank) {
      nm_uns=nm_ff[qrj$pivot[-(1:qrj$rank)]]
   } else {
      nm_uns=nm_ff
   }
   if (qrj$rank < nb_ff && !(least_norm || method!="nlsic")) {
      library(MASS)
      # Too bad. The jacobian of free fluxes is not of full rank.
      dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx)
      write.matrix(formatC(jx_f$dr_dff, 15), file="dbg_dr_dff_singular.txt", sep="\t")
      stop(paste("Provided measurements (isotopomers and fluxes) are not
   sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:
   ", paste(nm_uns, sep=", ", collapse=", "),
         "\nJacobian dr_dff is dumped in dbg_dr_dff_singular.txt", sep=""))
   }
   if (TIMEIT) {
      cat("optim   : ", date(), "\\n", sep="")
   }
   # pass control to the chosen optimization method
   res=opt_wrapper(measurements)
   if (any(is.na(res$par))) {
      obj2kvh(res, "optimization process informations", fkvh)
      stop("Optimization failed")
   }
#browser()
   if (zerocross && !is.null(mi_zc)) {
      # inverse active "zc" inequalities
      nm_inv=names(which((ui%*%res$par-ci)[,1]<=1.e-10))
      i=grep("^zc ", nm_inv, v=T)
      if (length(i) > 0) {
         i=str2ind(i, nm_i)
         cat("The following inequalities are active after first stage
of zero crossing strategy and will be inverted:\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="")
         ipos=grep(">=", nm_i[i], v=T)
         ineg=grep("<=", nm_i[i], v=T)
         ui[i,]=-ui[i,,drop=F]
         if (length(ipos)) {
            ipzc=str2ind(ipos, nm_izc)
            ipos=str2ind(ipos, nm_i)
            ci[ipos]=as.numeric(zzc+mi_zc[ipzc,,drop=F]%*%mic)
            nm_i[ipos]=sub(">=", "<=-", nm_i[ipos])
         }
         if (length(ineg)) {
            inzc=str2ind(ineg, nm_izc)
            ineg=str2ind(ineg, nm_i)
            ci[ineg]=as.numeric(zc+mi_zc[inzc,,drop=F]%*%mic)
            nm_i[ineg]=sub("<=-", ">=", nm_i[ineg])
         }
         rownames(ui)=nm_i
         names(ci)=nm_i
         # enforce new inequalities
         reopt=TRUE
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
         # reoptimize
         if (reopt) {
            res=opt_wrapper(measurements)
            if (any(is.na(res$par))) {
               obj2kvh(res, "optimization process informations", fkvh)
               stop("Second optimization failed")
            }
         }
      }
   }
   param=res$par
   names(param)=nm_par
   obj2kvh(res, "optimization process informations", fkvh)
}
if (TIMEIT) {
   cat("postopt : ", date(), "\\n", sep="")
}
# active constraints
ine=as.numeric(abs(ui%*%param-ci))<1.e-10
if (any(ine)) {
   obj2kvh(nm_i[ine], "active inequality constraints", fkvh)
}
poolall[nm_poolf]=exp(param[nm_poolf])

#browser()
rcost=cumo_cost(param, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
rres=cumo_resid(param, cjac=T, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
obj2kvh(rcost, "final cost", fkvh)
names(rres$res)=c(nm_meas, nm_fmn, nm_poolm)
o=order(names(rres$res))
obj2kvh(rres$res[o], "(simulated-measured)/sd_exp", fkvh)

# simulated measurements -> kvh
obj2kvh(cBind(value=jx_f$usimcumom,sd=measurements$dev$labeled), "simulated unscaled labeling measurements", fkvh)
if (nb_sc > 0) {
   val=jx_f$usimcumom*c(1.,param)[ir2isc]
   names(val)=nm_meas
   obj2kvh(cBind(value=val,sd=measurements$dev$labeled), "simulated scaled labeling measurements", fkvh)
}
if (nb_fmn) {
   obj2kvh(cBind(value=jx_f$fallnx[nm_fmn], sd=measurements$dev$flux), "simulated flux measurements", fkvh)
}
if (nb_poolm) {
   obj2kvh(cBind(value=measurements$mat$pool%*%poolall, sd=measurements$dev$pool), "simulated metabolite pool measurements", fkvh)
}

# gradient -> kvh
gr=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
names(gr)=nm_par
obj2kvh(gr, "gradient vector", fkvh)
obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh)

if (fullsys) {
   nm_flist=nm_list
   nm_flist$rcumo=nm_cumo
   v=param2fl_x(param, cjac=F, nb_f, nm_flist, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, spAbr_f, emu=F, , pool, measurements, ipooled)
} else {
   v=param2fl_x(param, cjac=F, nb_f, nm_list, nb_x, invAfl, p2bfl, g2bfl, bp, fc, xi, spa, emu, pool, measurements, ipooled)
}
# set final variable depending on param
x=v$x
if (fullsys) {
   names(x)=nm_cumo
} else {
   names(x)=nm_x
}

# write some info in result kvh
obj2kvh(cumo2mass(x), "MID vector", fkvh)
# replace decimal :i by binary :x1's
nm_mask=t(sapply(strsplit(names(x), "[:+]"), function(v) {
   metab=v[1]
   v[2]=int2bit(v[2], len=clen[metab])
   v
}))
o=order(nm_mask[,1], nm_mask[,2])
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T)
if (emu) {
   nm_mask=paste(paste(nm_mask[,1], nm_mask[,2], sep="#"), nm_mask[,3], sep="+M")
} else {
   nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#")
}
names(x)=nm_mask
obj2kvh(x[o], "labeling vector", fkvh)

fwrv=v$fwrv
names(fwrv)=nm_fwrv
#obj2kvh(fwrv, "fwd-rev flux vector", fkvh)

fallnx=v$fallnx
names(fallnx)=nm_fallnx

flnx=v$flnx
names(flnx)=c(nm_fl)

fgr=fallnx[nm_fgr]

if (sensitive=="mc") {
   if (TIMEIT) {
      cat("monte-ca: ", date(), "\\n", sep="")
   }
   # Monte-Carlo simulation in parallel way
   mc_inst=library(multicore, warn.conflicts=F, verbose=F, logical.return=T)
   simcumom=c(1.,param)[ir2isc]*jx_f$usimcumom
   simfmn=f[nm_fmn]
   simpool=as.numeric(measurements$mat$pool%*%poolall)
   mc_sim=function(i) {
      # random measurement generation
      if (nb_meas) {
         meas_mc=rnorm(nb_meas, simcumom, measurements$dev$labeled)
      } else {
         meas_mc=c()
      }
      if (nb_fmn) {
         fmn_mc=rnorm(nb_fmn, simfmn, measurements$dev$flux)
      } else {
         fmn_mc=c()
      }
      if (nb_poolm) {
         poolm_mc=rnorm(nb_fmn, simfmn, measurements$dev$pool)
      } else {
         poolm_mc=c()
      }
      #cat("imc=", i, "\\n", sep="")
      # minimization
      measurements_mc=measurements
      measurements_mc$vec$labeled=meas_mc
      measurements_mc$vec$flux=fmn_mc
      measurements_mc$vec$pool=poolm_mc
      res=opt_wrapper(measurements_mc, trace=0)
      # return the solution
      return(list(cost=sum(jx_f$res*jx_f$res), par=res$par))
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
   if (nb_poolf) {
      # from log to natural concentraion
      free_mc[nb_ff+nb_sc+1:nb_poolf,]=exp(free_mc[nb_ff+nb_sc+1:nb_poolf,])
   }
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
#browser()
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
   ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
   dimnames(ci_mc)[[2]][3]="length"
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
   obj2kvh((ci_mc-cBind(param, param, 0))*100/abs(param),
      "relative 95% confidence intervals (%)", fkvh, indent)
   
   # net-xch01 stats
   fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx)
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx
   if (length(fallnx_mc)) {
      dimnames(fallnx_mc)[[1]]=nm_fallnx
      # form a matrix output
      fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0)
      #cat("\\tall net-xch01 fluxes\\n", file=fkvh)
      # mean
#browser()
      fallout=cBind(fallout, mean=apply(fallnx_mc, 1, mean))
      #obj2kvh(apply(fallnx_mc, 1, mean), "mean", fkvh, indent)
      # median
      parmed=apply(fallnx_mc, 1, median)
      fallout=cBind(fallout, median=parmed)
      #obj2kvh(parmed, "median", fkvh, indent)
      # covariance matrix
      covmc=cov(t(fallnx_mc))
      dimnames(covmc)=list(nm_fallnx, nm_fallnx)
      #obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      fallout=cBind(fallout, sd=sdmc)
      #obj2kvh(sdmc, "sd", fkvh, indent)
      fallout=cBind(fallout, "rsd (%)"=sdmc*100/abs(fallnx))
      #obj2kvh(sdmc*100/abs(fallnx), "rsd (%)", fkvh, indent)
      # confidence intervals
      ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
      ci_mc=cBind(ci_mc, (ci_mc-cBind(fallnx, fallnx, 0))*100/abs(fallnx))
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)")
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
      fallout=cBind(fallout, ci_mc)
      #obj2kvh((ci_mc-cBind(fallnx, fallnx, 0))*100/abs(fallnx),
      #   "relative 95% confidence intervals (%)", fkvh, indent)
      o=order(nm_fallnx)
      obj2kvh(fallout[o,,drop=F], "all net-xch01 fluxes", fkvh, indent)
      obj2kvh(covmc[o,o], "covariance of all net-xch01 fluxes", fkvh, indent)
      
      # fwd-rev stats
      fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fwrv)
      dimnames(fwrv_mc)[[1]]=nm_fwrv
      fallout=matrix(0, nrow=nrow(fwrv_mc), ncol=0)
      #cat("\\tforward-reverse fluxes\\n", file=fkvh)
      # mean
      fallout=cBind(fallout, mean=apply(fwrv_mc, 1, mean))
      #obj2kvh(apply(fwrv_mc, 1, mean), "mean", fkvh, indent)
      # median
      parmed=apply(fwrv_mc, 1, median)
      fallout=cBind(fallout, median=parmed)
      #obj2kvh(parmed, "median", fkvh, indent)
      # covariance matrix
      covmc=cov(t(fwrv_mc))
      dimnames(covmc)=list(nm_fwrv, nm_fwrv)
      #obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      fallout=cBind(fallout, sd=sdmc)
      #obj2kvh(sdmc, "sd", fkvh, indent)
      fallout=cBind(fallout, "rsd (%)"=sdmc*100/abs(fwrv))
      #obj2kvh(sdmc*100/abs(fwrv), "rsd (%)", fkvh, indent)
      # confidence intervals
      ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
      ci_mc=cBind(ci_mc, (ci_mc-cBind(fwrv, fwrv, 0))*100/abs(fwrv))
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)")
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent)
      #obj2kvh((ci_mc-cBind(fwrv, fwrv, 0))*100/abs(fwrv),
      #   "relative 95% confidence intervals (%)", fkvh, indent)
      fallout=cBind(fallout, ci_mc)
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
   rres=cumo_resid(param, cjac=T, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
} # else use last calculated jacobian

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
# standart deviations of free fluxes
cat("linear stats\\n", file=fkvh)

# sd free+dependent+growth net-xch01 fluxes
nm_flfd=c(nm_ff, nm_fgr, nm_fl)
if (nb_ff > 0 || nb_fgr > 0) {
   i=1:nb_param
   i=c(head(i, nb_ff), tail(i, nb_fgr))
   covfl=crossprod(rtcov[, i, drop=F]%mmt%(rBind(diag(nb_ff+nb_fgr), dfl_dffg)%mrv%c(rep.int(1., nb_ff), fgr)))
   dimnames(covfl)=list(nm_flfd, nm_flfd)
   sdfl=sqrt(diag(covfl))
} else {
   sdfl=rep(0., nb_fl)
}
fl=c(head(param, nb_ff), fgr, flnx)
mtmp=cBind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl))
rownames(mtmp)=nm_flfd
o=order(nm_flfd)
obj2kvh(mtmp[o,,drop=F], "net-xch01 fluxes (sorted by name)", fkvh, indent=1)
obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1)

# sd of all fwd-rev
if (nb_ff > 0 || nb_fgr > 0) {
   i=1:nb_param
   i=c(head(i, nb_ff), tail(i, nb_fgr))
   covf=crossprod(tcrossprod(rtcov[,i, drop=F], jx_f$df_dffp%mrv%c(rep.int(1., nb_ff), head(poolall[nm_poolf], nb_fgr))))
   dimnames(covf)=list(nm_fwrv, nm_fwrv)
   sdf=sqrt(diag(covf))
} else {
   sdf=rep(0., length(fwrv))
}
mtmp=cBind(fwrv, sdf, sdf/abs(fwrv))
dimnames(mtmp)[[2]]=c("value", "sd", "rsd")
o=order(nm_fwrv)
obj2kvh(mtmp[o,], "fwd-rev fluxes (sorted by name)", fkvh, indent=1)
if (nb_ff > 0 || nb_fgr > 0) {
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
   covpf=crossprod(rtcov[,nb_ff+nb_sc+1:nb_poolf, drop=F]%mrv%poolall[nm_poolf])
   dimnames(covpf)=list(nm_poolf, nm_poolf)
   sdpf[nm_poolf]=sqrt(diag(covpf))
}
if (length(poolall) > 0) {
   mtmp=cBind("value"=poolall, "sd"=sdpf, "rsd"=sdpf/poolall)
   rownames(mtmp)=nm_poolall
   o=order(nm_poolall)
   obj2kvh(mtmp[o,,drop=F], "metabolite pools (sorted by name)", fkvh, indent=1)
   if (nb_poolf > 0) {
      o=order(nm_poolf)
      obj2kvh(covpf[o, o], "covariance free pools", fkvh, indent=1)
   }
}

# khi2 test for goodness of fit
# goodness of fit (khi2 test)
khi2test=list("khi2 value"=rcost, "data points"=length(jx_f$res),
   "fitted parameters"=nb_param, "degrees of freedom"=length(jx_f$res)-nb_param)
khi2test$`khi2 reduced value`=khi2test$`khi2 value`/khi2test$`degrees of freedom`
khi2test$`p-value, i.e. P(X^2<=value)`=pchisq(khi2test$`khi2 value`, df=khi2test$`degrees of freedom`)
khi2test$conclusion=if (khi2test$`p-value, i.e. P(X^2<=value)` > 0.95) "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD" else "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD"
obj2kvh(khi2test, "goodness of fit (khi2 test)", fkvh, indent=1)

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

# write edge.xchflux property
fedge=file("%(d)s/edge.xchflux.%(org)s", "w")
flxch=paste(".x", substring(edge2fl, 4), sep="")
ifl=charmatch(flxch, substring(names(fallnx), 2))
cat("xchflux (class=Double)\\n", sep="", file=fedge)
cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
close(fedge)

# write node.log2pool property
if (length(poolall)> 0) {
   fnode=file("%(d)s/node.log2pool.%(org)s", "w")
   cat("log2pool (class=Double)\\n", sep="", file=fnode)
   nm_node=substring(names(poolall), 4)
   cat(paste(nm_node, log2(poolall), sep=" = "), sep="\\n" , file=fnode)
   close(fnode)
}
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
