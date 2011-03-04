#!/usr/bin/env python
#
# Transform an ftbl to R code which will solve an optimization of flux analysis
# problem min(S) over \Theta, where S=||Predicted-Observed||^2_\Sigma^2
# and \Theta is a vector of free fluxes (net+xch) and scaling parameters.
# Predicted vector is obtained from cumomer vector x (calculated from
# free fluxes and divided in chunks according to the cumo weight) by
# multiplying it by the measure matrices and scale factor, boths coming
# from ftbl file. Observed values vector xo is extracted from ftbl file.
# it is composed of flux and cumomer measures.
# \Sigma^2, covariance diagonal matrices sigma[flux|mass|label|peak]
# are orginated from ftbl

# usage: ./ftbl2optR.py organism
# where organism is the ftbl informative part of file name
# (before .ftbl), e.g. organism.ftbl
# after execution a file organism.R will be created.
# If it already exists, it will be silently overwritten.
# The system Afl*flnx=bfl is created from ftbl file.
# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used
#   for sparse matrices.

# Important python variables:
# Collections:
#    netan - (dict) ftbl structured content;
#    tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
#    collection;
#    measures - (dict) exp data;
#    rAb - (list) reduced linear systems A*x_cumo=b by weight;
#    scale - unique scale names;
#    nrow - counts scale names;
#    o_sc - ordered scale names;
#    o_meas - ordered measure types;
# org - (str) prefix of .ftbl  file like "PPP"
# File names (str):
#    n_ftbl (descriptor f_ftbl);
#    n_R (R code) (f);
#    n_fort (fortran code) (ff);
# Counts: nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total),
#         nb_ffn, nb_ffx (free fluxes)
# Index translators:
#    fwrv2i - flux names to index in R:fwrv;
#    cumo2i - cumomer names to index in R:x;
#    ir2isc - mapping measure rows indexes on scale index isc[meas]=ir2isc[meas][ir]
# Vector names:
#    cumos (list) - names of R:x;
#    o_mcumos - cumomers involved in measures;

# Important R variables:
# Scalars:
#    nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
#    nb_ffn, nb_ffx, nb_ff (free fluxes),
#    nb_fcn, nb_fcx, nb_fc (constrained fluxes),
#    nb_ineq, nb_param, nb_fmn
# Name vectors:
#    nm_cumo, nm_fwrv, nm_fallnx, nm_fln, nm_flx, nm_fl, nm_par,
#    nm_ffn, nm_ffx,
#    nm_fcn, nm_fcx,
#    nm_mcumo, nm_fmn
# Numeric vectors:
#    fwrv - all fluxes (fwd+rev);
#    x - all cumomers (weight1+weight2+...);
#    param - free flux net, free flux xch, scale label, scale mass, scale peak
#    fcn, fcx, fc,
#    bp - helps to construct the rhs of flux system
#    xi -cumomer input vector
#    fallnx - complete flux vector (constr+net+xch)
#    bc - helps to construct fallnx
#    li - inequality vector (mi%*%fallnx>=li)
#    ir2isc - measur row to scale vector replicator
#    ci - inequalities for param use (ui%*%param-ci>=0)
#    measvec,
#    measinvvar,
#    imeas,
#    fmn
# Matrices:
#    Afl, qrAfl, invAfl,
#    p2bfl - helps to construct the rhs of flux system
#    mf, md - help to construct fallnx
#    mi - inequality matrix (ftbl content)
#    ui - inequality matrix (ready for param use)
#    measmat - measmat*(x[imeas];1)=vec of simulated not-yet-scaled measures
# Functions:
#    param2fl_x - translate param to flux and cumomer vector (initial approximation)
#    cumo_cost - cost function (khi2)
#    cumo_grad - finite difference gradient
#    cumo_gradj - implicit derivative gradient

import sys;
import os;
import time;
import copy;
import getopt;

#sys.path.append("/home/sokol/dev/python");
from tools_ssg import *;
import C13_ftbl;
try:
    import psyco;
    psyco.full();
    #psyco.log();
    #psyco.profile();
except ImportError:
    #print 'Psyco not installed, the program will just run slower'
    pass;

me=os.path.basename(sys.argv[0]);
def usage():
    sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--DEBUG] network_name[.ftbl]\n");

#<--skip in interactive session
try:
    opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "DEBUG"]);
except getopt.GetoptError, err:
    #pass;
    sys.stderr.write(str(err)+"\n");
    usage();
    #sys.exit(1);

fullsys=False;
DEBUG=False;
for o,a in opts:
    if o in ("-h", "--help"):
        usage();
        sys.exit(0);
    elif o=="--cost":
        cost=True;
    elif o=="--fullsys":
        fullsys=True;
    elif o=="--DEBUG":
        DEBUG=True;
    else:
        #assert False, "unhandled option";
        # unknown options can come from shell
        # which passes all options both to python and R scripts
        # so just ignore unknown options
        pass;
#aff("args", args);##
if len(args) != 1:
    usage();
    exit(1);
org=args[0];

# cut .ftbl if any
if org[-5:]==".ftbl":
    org=org[:-5];

#-->
#DEBUG=True;
import ftbl2code;
ftbl2code.DEBUG=DEBUG;
#org="ex3";
#org="PPP_exact";
#DEBUG=True;
if DEBUG:
    import pdb;


n_ftbl=org+".ftbl";
n_R=org+".R";
#n_fort=org+".f";
f_ftbl=open(n_ftbl, "r");
os.system("chmod u+w '%s' 2>/dev/null"%n_R);
#os.system("chmod u+w '%s' 2>/dev/null"%n_fort);
f=open(n_R, "w");
#ff=open(n_fort, "w");

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();

# analyse network
# reload(C13_ftbl);

netan=C13_ftbl.ftbl_netan(ftbl);

# write initialization part of R code
ftbl2code.netan2Rinit(netan, org, f, fullsys);

#f.write(
"""
# output flux repartition
cat("Dependent fluxes:\\n");
if (nb_fln) {
   print(paste(nm_fln,"net",sep="_"));
}
if (nb_flx) {
   print(paste(nm_flx,"xch",sep="_"));
}
cat("Free fluxes:\\n");
if (nb_ffn) {
   print(paste(nm_ffn,"net",sep="_"));
}
if (nb_ffx) {
   print(paste(nm_ffx,"xch",sep="_"));
}
cat("Constrained fluxes:\\n");
if (nb_fcn) {
   print(paste(nm_fcn,"net",sep="_"));
}
if (nb_fcx) {
   print(paste(nm_fcx,"xch",sep="_"));
}
"""
f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\n", sep="");
}
#browser()
# check if initial approximation is feasible
ineq=ui%*%param-ci;
if (any(ineq < 0.)) {
   cat("The following ", sum(ineq<0.), " ineqalities are not respected:\\n", sep="");
   print(ineq[ineq<0.,1]);
   cat("Starting values are put in feasible domain:\\n", sep="");
   #stop("Inequalities violated, cf. log file.");
   # put them inside
   dp=ldp(ui, -ineq);
   names(param)=nm_par;
   tmp=cbind(param[1:nb_ff], param[1:nb_ff]+1.001*dp[1:nb_ff]);
   dimnames(tmp)=list(nm_par[1:nb_ff], c("old", "new"));
   obj2kvh(tmp, "starting free fluxes");
   # move starting point slightly inside of feasible domain
   param=param+1.001*dp;
}

# set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
# for corresponding measures
vr=param2fl_x(param, cjac=FALSE, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, spAb);
if (!is.null(vr$err) && vr$err) {
   stop(vr$mes);
}
simvec=(measmat%*%c(vr$x[irmeas],1.));
if (DEBUG) {
   cat("initial simvec:\\n");
   print(simvec);
}
if (nb_ff < length(param)) {
   ms=measvec*simvec*measinvvar;
   ss=simvec*simvec*measinvvar;
   for (i in (nb_ff+1):length(param)) {
      im=(ir2isc==(i+1));
      param[i]=sum(ms[im])/sum(ss[im]);
   }
}

# prepare flux index conversion
ifwrv=1:nb_fwrv;
names(ifwrv)=nm_fwrv;
ifl_in_fw=ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")];
iff_in_fw=ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")];
# index couples for jacobian df_dfl, df_dffd
cfw_fl=crv_fl=cbind(ifl_in_fw, 1:nb_fl);
cfw_ff=crv_ff=cbind(iff_in_fw, 1:nb_ff);
crv_fl[,1]=(nb_fwrv/2)+crv_fl[,1];
crv_ff[,1]=(nb_fwrv/2)+crv_ff[,1];

# see if there are any active inequalities at starting point
ineq=ui%*%param-ci;
if (any(abs(ineq)<=1.e-10)) {
   cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point:\\n",
      paste(names(ineq[abs(ineq)<=1.e-10,1]), collapse="\\n"), "\\n", sep="");
}
""");

f.write("""
# formated output in kvh file
fkvh=file("%(org)s_res.kvh", "w");
"""%{
    "org": org,
});
# main part: call optimization
f.write("""
# save options of command line
cat("command line\\n", file=fkvh);
obj2kvh(opts, "opts", fkvh, indent=1);

# resume system sizes
obj2kvh(nb_sys, "system size resume", fkvh);

# save initial flux and cumomer distribution
cat("starting point\\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "starting free parameters", fkvh, indent=1);

x=vr$x;
names(x)=nm_rcumo;
obj2kvh(cumo2mass(x), "starting MID vector", fkvh, indent=1);
# replace :i by #x1's
nm_mask=unlist(lapply(strsplit(names(x), ":"), function(v) {
   n=length(v);
   metab=paste(v[-n], collapse=":");
   mask=int2bit(v[n], len=clen[metab])
   c(metab, mask);
}))
nm_mask=t(matrix(nm_mask, nrow=2));
o=order(nm_mask[,1], nm_mask[,2]);
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T);
nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#");
names(nm_mask)=names(x);
names(x)=nm_mask;
obj2kvh(x[o], "starting cumomer vector", fkvh, indent=1);

fwrv=vr$fwrv;
n=length(fwrv);
names(fwrv)=nm_fwrv;
obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, indent=1);

f=vr$fallnx;
n=length(f);
names(f)=nm_fallnx;
obj2kvh(f, "starting net-xch01 flux vector", fkvh, indent=1);

rres=cumo_resid(param, cjac=TRUE, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "starting residuals (simulated-measured)/sd_exp", fkvh, indent=1);

rcost=cumo_cost(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
obj2kvh(rcost, "starting cost value", fkvh, indent=1);

obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1);
btmp=c(p2bfl%*%param[1:nb_f$nb_ff]+bp);
names(btmp)=dimnames(Afl)[[1]];
obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1);
names(measvec)=nm_meas;
obj2kvh(measvec, "measure vector", fkvh, indent=1);

#cat("mass vector:\\n");
#print_mass(x);

names(param)=nm_par;

opt_wrapper=function(measvec, fmn) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAb);
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAb);
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-07,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAb);
   } else if (method == "nlsic") {
      control=list(trace=1, btfrac=0.25, btdesc=0.75, maxit=50, errx=1.e-5);
      res=nlsic(param, cumo_resid, 
         ui, ci, control, e=NULL, eco=NULL,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar,
         spAb);
      if (res$err || is.null(res$par)) {
         # store res in kvh
         obj2kvh(res, "optimization aborted here", fkvh);
         cat(res$mes, file=stderr());
         res$par=rep(NA, length(param))
         res$cost=NA
      }
   } else {
      stop(paste("Unknown minimization method '", method, "'", sep=""));
   }
   return(res);
}
if (optimize) {
   # check if at starting position all fluxes can be resolved
   #browser();
   qrj=qr(jx_f$dr_dff);
   d=diag(qrj$qr)
   qrj$rank=sum(abs(d)>abs(d[1])*1.e-14)
   if (qrj$rank < nb_ff) {
      # Too bad. The jacobian of free fluxes is not of full rank.
      dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx);
      write.matrix(formatC(jx_f$dr_dff, 15), file="dbg_dr_dff_singular.txt", sep="\t");
      stop(paste("Provided measures (isotopomers and fluxes) are not
   sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:
   ", paste(nm_par[qrj$pivot[-(1:qrj$rank)]], sep=", ", collapse=", "),
         "\nJacobian dr_dff is dumped in dbg_dr_dff_singular.txt", sep=""));
   }
   if (TIMEIT) {
      cat("optim   : ", date(), "\n", sep="");
   }
   # pass control to the true method
   res=opt_wrapper(measvec, fmn);
   param=res$par;
   names(param)=nm_par;
   obj2kvh(res, "optimization process informations", fkvh);
   if (any(is.na(res$par))) {
      stop("Optimization failed")
   }
}
if (TIMEIT) {
   cat("postopt : ", date(), "\n", sep="");
}
# active constraints
ine=abs(ui%*%param-ci)<1.e-10;
if (any(ine)) {
   obj2kvh(nm_i[ine], "active inequality constraints", fkvh);
}
#browser();
rcost=cumo_cost(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
rres=cumo_resid(param, cjac=T, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
obj2kvh(rcost, "final cost", fkvh);
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "(simulated-measured)/sd_exp", fkvh);
gr=2*c((t(jx_f$ures*c(measinvvar, invfmnvar))%*%jx_f$udr_dp));
names(gr)=nm_par;
obj2kvh(gr, "gradient vector", fkvh);
obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh);

if (fullsys) {
   v=param2fl_x(param, cjac=F, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAb_f);
} else {
   v=param2fl_x(param, cjac=F, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, spAb);
}
x=v$x;
if (fullsys) {
   names(x)=nm_cumo;
} else {
   names(x)=nm_rcumo;
}
obj2kvh(cumo2mass(x), "MID vector", fkvh);
# replace :i by :x1's
nm_mask=unlist(lapply(strsplit(names(x), ":"), function(v) {
   n=length(v);
   metab=paste(v[-n], collapse=":");
   mask=int2bit(v[n], len=clen[metab])
   c(metab, mask);
}))
nm_mask=t(matrix(nm_mask, nrow=2));
o=order(nm_mask[,1], nm_mask[,2]);
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T);
nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#");
names(nm_mask)=names(x);
names(x)=nm_mask;
obj2kvh(x[o], "cumomer vector", fkvh);

fwrv=v$fwrv;
names(fwrv)=nm_fwrv;
obj2kvh(fwrv, "fwd-rev flux vector", fkvh);

f=v$fallnx;
names(f)=nm_fallnx;
obj2kvh(f, "net-xch01 flux vector", fkvh);

flnx=v$flnx;
names(flnx)=c(nm_fln, nm_flx);

if (sensitive=="mc") {
   if (TIMEIT) {
      cat("monte-ca: ", date(), "\n", sep="");
   }
   # Monte-Carlo simulation in parallel way
   library(multicore, warn.conflicts=F, verbose=F);
   invar=c(measinvvar, invfmnvar);
   simcumom=c(1.,param)[ir2isc]*jx_f$usimcumom;
   simfmn=f[nm_fmn];
   mc_sim=function(i) {
      # random measure generation
      if (nb_meas) {
         meas_mc=rnorm(nb_meas, simcumom, 1./sqrt(measinvvar));
      } else {
         meas_mc=c();
      }
      if (nb_fmn) {
         fmn_mc=rnorm(nb_fmn, simfmn, 1./sqrt(invfmnvar));
      } else {
         fmn_mc=c();
      }
      cat("imc=", i, "\n", sep="");
      # minimization
      res=opt_wrapper(meas_mc, fmn_mc);
      # return the solution
      return(list(cost=sum(jx_f$res*jx_f$res), par=res$par));
   }
   # parallel execution
   mc_res=mclapply(1:nmc, mc_sim);
   free_mc=matrix(unlist(mc_res), ncol=nmc);
   cost_mc=free_mc[1,];
   free_mc=free_mc[-1,,drop=F];
   # remove failed m-c iterations
   ifa=which(is.na(free_mc[1,]))
   if (length(ifa)) {
      warning("Some Monte-Carlo iterations failed.")
      free_mc=free_mc[,-ifa,drop=F]
      cost_mc=cost_mc[-ifa]
      nmc=length(cost_mc)
   }
#browser();
   dimnames(free_mc)[[1]]=nm_par;
   cat("monte-carlo\n", file=fkvh);
   indent=1;
   obj2kvh(nmc, "sample number", fkvh, indent);
   cat("\tcost\n", file=fkvh);
   indent=2;
   obj2kvh(mean(cost_mc), "mean", fkvh, indent);
   obj2kvh(median(cost_mc), "median", fkvh, indent);
   obj2kvh(sd(cost_mc), "sd", fkvh, indent);
   obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent);
   obj2kvh(quantile(cost_mc, c(0.025, 0.975)), "ci", fkvh, indent);
   cat("\tfree parameters\n", file=fkvh);
   indent=2;
   # param stats
   # mean
   obj2kvh(apply(free_mc, 1, mean), "mean", fkvh, indent);
   # meadian
   parmed=apply(free_mc, 1, median);
   obj2kvh(parmed, "median", fkvh, indent);
   # covariance matrix
   covmc=cov(t(free_mc));
   obj2kvh(covmc, "covariance", fkvh, indent);
   # sd
   sdmc=sqrt(diag(covmc));
   obj2kvh(sdmc, "sd", fkvh, indent);
   obj2kvh(sdmc*100/abs(param), "rsd (%)", fkvh, indent);
   # confidence intervals
   ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)));
   ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
   dimnames(ci_mc)[[2]][3]="length";
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
   obj2kvh((ci_mc-cbind(param, param, 0))*100/abs(param),
      "relative 95% confidence intervals (%)", fkvh, indent);
   
   # net-xch01 stats
   fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, invAfl, p2bfl, bp, fc)$fallnx);
   fallnx=param2fl(param, nb_f, invAfl, p2bfl, bp, fc)$fallnx;
   dimnames(fallnx_mc)[[1]]=nm_fallnx;
   cat("\tall net-xch01 fluxes\n", file=fkvh);
   # mean
   obj2kvh(apply(fallnx_mc, 1, mean), "mean", fkvh, indent);
   # meadian
   parmed=apply(fallnx_mc, 1, median);
   obj2kvh(parmed, "median", fkvh, indent);
   # covariance matrix
   covmc=cov(t(fallnx_mc));
   obj2kvh(covmc, "covariance", fkvh, indent);
   # sd
   sdmc=sqrt(diag(covmc));
   obj2kvh(sdmc, "sd", fkvh, indent);
   obj2kvh(sdmc*100/abs(fallnx), "rsd (%)", fkvh, indent);
   # confidence intervals
   ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)));
   ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
   dimnames(ci_mc)[[2]][3]="length";
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
   obj2kvh((ci_mc-cbind(fallnx, fallnx, 0))*100/abs(fallnx),
      "relative 95% confidence intervals (%)", fkvh, indent);

   # fwd-rev stats
   fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, invAfl, p2bfl, bp, fc)$fwrv);
   dimnames(fwrv_mc)[[1]]=nm_fwrv;
   cat("\tforward-reverse fluxes\n", file=fkvh);
   # mean
   obj2kvh(apply(fwrv_mc, 1, mean), "mean", fkvh, indent);
   # meadian
   parmed=apply(fwrv_mc, 1, median);
   obj2kvh(parmed, "median", fkvh, indent);
   # covariance matrix
   covmc=cov(t(fwrv_mc));
   obj2kvh(covmc, "covariance", fkvh, indent);
   # sd
   sdmc=sqrt(diag(covmc));
   obj2kvh(sdmc, "sd", fkvh, indent);
   obj2kvh(sdmc*100/abs(fwrv), "rsd (%)", fkvh, indent);
   # confidence intervals
   ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)));
   ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
   dimnames(ci_mc)[[2]][3]="length";
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
   obj2kvh((ci_mc-cbind(fwrv, fwrv, 0))*100/abs(fwrv),
      "relative 95% confidence intervals (%)", fkvh, indent);
} else if (length(sensitive) && nchar(sensitive)) {
   stop(paste("Unknown sensitivity '", sensitive, "' method chosen.", sep=""));
}

if (TIMEIT) {
   cat("linstats: ", date(), "\n", sep="");
}
# Linear method based on jacobian x_f
# reset fluxes and jacobians according to param
rres=cumo_resid(param, cjac=T, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
#grj=cumo_gradj(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
if (DEBUG) {
   library(numDeriv); # for numerical jacobian
   # numerical simulation
   rj=function(v, ...) {
      r=cumo_resid(v, cjac=F, ...);
      c(r$res);
   }
   #dr_dpn=jacobian(rj, param, method="Richardson", method.args=list(),
   #   nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
   dx_dfn=num_jacob(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
   dr_dp=num_fjacob(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
   # to compare with jx_f$dr_dp
   gr=grad(cumo_cost, param, method="Richardson", method.args=list(), nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
   grn=cumo_grad(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
#browser();
}

# covariance matrix of free fluxes
invcov=t(jx_f$dr_dff)%*%(jx_f$dr_dff*c(measinvvar, invfmnvar));
covff=try(solve(invcov));
if (inherits(covff, "try-error")) {
   # matrix seems to be singular
   svi=svd(invcov);
   i=svi$d/svi$d[1]<1.e-14;
   ibad=apply(svi$u[, i, drop=F], 2, which.max);
   cat("Inverse of covariance matrix is singular.\\nStatistically undefined fluxe(s) seems to be:\\n",
      paste(nm_ff[ibad], collapse="\\n"), "\\n", sep="", file=stderr());
   # regularize and inverse
   covff=solve(invcov+diag(1.e-14*svi$d[1], nrow(invcov)));
}

# standart deviations of free fluxes
sdff=sqrt(diag(covff));
cat("linear stats\\n", file=fkvh);
mtmp=cbind(param[1:nb_ff], sdff, sdff/abs(param[1:nb_ff]));
dimnames(mtmp)[[2]]=c("val", "sd", "rsd");
o=order(mtmp[,"rsd"]);
obj2kvh(mtmp[o,], "val, sd, rsd free fluxes (sorted by rsd)", fkvh, indent=1);
obj2kvh(covff, "covariance free fluxes", fkvh, indent=1);

# sd dependent net-xch01 fluxes
covfl=dfl_dff%*%covff%*%t(dfl_dff);
sdfl=sqrt(diag(covfl));
mtmp=cbind(flnx, sdfl, sdfl/abs(flnx));
dimnames(mtmp)[[2]]=c("val", "sd", "rsd");
o=order(mtmp[,"rsd"]);
obj2kvh(mtmp[o,], "val, sd, rsd dependent net-xch01 fluxes (sorted by rsd)", fkvh, indent=1);
obj2kvh(sdfl, "sd net-xch01 fluxes", fkvh, indent=1);
obj2kvh(covfl, "covariance net-xch01 fluxes", fkvh, indent=1);

# sd of all fwd-rev
covf=jx_f$df_dff%*%covff%*%t(jx_f$df_dff);
sdf=sqrt(diag(covf));
mtmp=cbind(fwrv, sdf, sdf/abs(fwrv));
dimnames(mtmp)[[2]]=c("val", "sd", "rsd");
o=order(mtmp[,"rsd"]);
obj2kvh(mtmp[o,], "val, sd, rsd fwd-rev fluxes", fkvh, indent=1);
obj2kvh(sdf, "sd fwd-rev fluxes", fkvh, indent=1);
obj2kvh(covf, "covariance fwd-rev fluxes", fkvh, indent=1);

# select best defined flux combinations
s=svd(covfl);
# lin comb coeff matrix
l=t(s$u);
iord=apply(l,1,function(v){o=order(abs(v), decreasing=T); n=which(cumsum(v[o]**2)>=0.999); o[1:n[1]]})
dimnames(l)=dimnames(covfl);
nm_fl=dimnames(l)[[1]];
nb_fl=nrow(s$u);
comb=rep("", nb_fl);
for (j in 1:nb_fl) {
   cva=l[j,iord[[j]]]; # abs(factor) ordered j-th lin combination
   cva=sign(cva[1])*cva;
   cva=sprintf("%+.3g*", cva);
   cva[cva=="+1*"]="+";
   cva[cva=="-1*"]="-";
   cva[1]=substring(cva[1], 2);
   comb[j]=paste(cva, nm_fl[iord[[j]]], sep="", collapse="");
}
# from best to worst defined sd
sta=sqrt(abs(diag(l%*%covfl%*%t(l))));
names(sta)=comb;
o=order(sta);
sta=sta[o];
obj2kvh(sta, "sd for uncorrelated net-xch01 linear combinations", fkvh, indent=1);
if (prof) {
   Rprof(NULL);
}
close(fkvh);
if (TIMEIT) {
   cat("rend    : ", date(), "\n", sep="");
}
""");

f.close();
#ff.close();
# make output files just readable to avoid later casual edition
os.system("chmod a-w '%s'"%n_R);
#os.system("chmod a-w '%s'"%n_fort);
