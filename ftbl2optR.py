#!/usr/bin/env python
#
# Transform an ftbl to R code solving an optimization of flux analysis
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
# The generated R code will use organism_sym.R file (A*x=b for cumomers,
# cf. ftbl2symA.py)
# The system Afl*flnx=bfl is created from ftbl file.
# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package

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
    sys.stderr.write("usage: "+me+" network_name[.ftbl]\n");

#<--skip in interactive session
if len(sys.argv) < 2:
    usage();
    sys.exit(1);

# set some python constants
org=sys.argv[1];
# cut .ftbl if any
if org[-5:]==".ftbl":
    org=org[:-5];

DEBUG=True if len(sys.argv) > 2 and sys.argv[2] else False;
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
n_fort=org+".f";
f_ftbl=open(n_ftbl, "r");
os.system("chmod u+w '%s' 2>/dev/null"%n_R);
os.system("chmod u+w '%s' 2>/dev/null"%n_fort);
f=open(n_R, "w");
ff=open(n_fort, "w");

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();

# analyse network
# reload(C13_ftbl);

netan=C13_ftbl.ftbl_netan(ftbl);

# write initialization part of R code
ftbl2code.netan2Rinit(netan, org, f, ff);

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
# set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
# for corresponding measures
vr=param2fl_x(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
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

# check if initial approximation is feasible
ineq=ui%*%param-ci;
if (!all(ineq > 1.e-10)) {
   if (sum(abs(ineq)<=1.e-10)) {
      cat("The following ", sum(ineq==0.), " ineqalities are on the border:\\n", sep="");
      print(ineq[abs(ineq)<=1.e-10,1]);
      #stop("Inequalities on the border, cf. log file.");
   }
   if (sum(ineq<0.)) {
      cat("The following ", sum(ineq<0.), " ineqalities are not respected:\\n", sep="");
      print(ineq[ineq<0.,1]);
      #stop("Inequalities violated, cf. log file.");
   }
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

# save initial flux and cumomer distribution
cat("starting point\\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "starting free parameters", fkvh, indent=1);

x=vr$x;
names(x)=nm_rcumo;
obj2kvh(x, "starting cumomer vector", fkvh, indent=1);

fwrv=vr$fwrv;
n=length(fwrv);
names(fwrv)=nm_fwrv;
obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, indent=1);

f=vr$fallnx;
n=length(f);
names(f)=nm_fallnx;
obj2kvh(f, "starting net-xch flux vector", fkvh, indent=1);

rres=cumo_resid(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, "fwrv2rAbcumo", "frj_rhs");
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "starting residuals", fkvh, indent=1);

rcost=cumo_cost(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo", "frj_rhs");
obj2kvh(rcost, "starting cost value", fkvh, indent=1);

obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1);
btmp=c(p2bfl%*%param[1:nb_f$nb_ff]+bp);
names(btmp)=dimnames(Afl)[[1]];
obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1);
names(measvec)=nm_meas;
obj2kvh(measvec, "cumomer measure vector", fkvh, indent=1);

#cat("mass vector:\\n");
#print_mass(x);

names(param)=nm_par;

# check if at starting position all fluxes can be resolved
qrj=qr(jx_f$dr_dff);
if (qrj$rank < nb_ff) {
   # Too bad. The jacobian of free fluxes is not of full rank.
   stop(paste("Provided measures (isotopomers and fluxes) are not
sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:
", paste(nm_par[qrj$pivot[-(1:qrj$rank)]], sep=", ", collapse=", ")));
}
opt_wrapper=function(measvec, fmn) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-07,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, "fwrv2rAbcumo", "frj_rhs");
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, "fwrv2rAbcumo", "frj_rhs");
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=1);
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-07,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, "fwrv2rAbcumo", "frj_rhs");
   } else if (method == "nlsic") {
      control=list(trace=1);
      res=nlsic(param, cumo_resid, 
         ui, ci, control, e=NULL, eco=NULL,
         nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc,
         irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar,
         "fwrv2rAbcumo", "frj_rhs");
      if (res$err || any(is.na(res$par))) {
         # store res in kvh
         obj2kvh(res, "optimization aborted here", fkvh);
         stop(res$mes);
      }
   } else {
      stop(paste("Unknown minimization method '", method, "'", sep=""));
   }
   return(res);
}
if (optimize) {
   # pass control to the true method
   res=opt_wrapper(measvec, fmn);
   param=res$par;
   names(param)=nm_par;
   obj2kvh(res, "optimization process informations", fkvh);
}
# active constraints
ine=abs(ui%*%param-ci)<1.e-10;
if (any(ine)) {
   obj2kvh(nm_i[ine], "active inequality constraints", fkvh);
}
rcost=cumo_cost(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo", "frj_rhs");
rres=cumo_resid(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, "fwrv2rAbcumo", "frj_rhs");
obj2kvh(rcost, "final cost", fkvh);
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "residual vector", fkvh);
gr=2*c((t(jx_f$res*c(measinvvar, invfmnvar))%*%jx_f$dr_dp));
names(gr)=nm_par;
obj2kvh(gr, "gradient vector", fkvh);
obj2kvh(jx_f$dr_dp, "jacobian dr_dp", fkvh);

v=param2fl_x(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas,
   measmat, measvec, ir2isc, "fwrv2Abcumo", "fj_rhs");
x=v$x;
names(x)=nm_cumo;
o=order(nm_cumo);
obj2kvh(x[o], "cumomer vector", fkvh);

fwrv=v$fwrv;
n=length(fwrv);
names(fwrv)=nm_fwrv;
obj2kvh(fwrv, "fwd-rev flux vector", fkvh);

f=v$fallnx;
n=length(f);
names(f)=nm_fallnx;
obj2kvh(f, "net-xch flux vector", fkvh);

if (sensitive=="mc") {
   # Monte-Carlo simulation
   # random measure generation
   nb_meas=length(measvec);
   invar=c(measinvvar, invfmnvar);
   simcumom=jx_f$simcumom;
   simfmn=f[nm_fmn];
   meas_mc=matrix(rnorm(nmc*nb_meas, simcumom, 1./sqrt(measinvvar)), nb_meas, nmc);
   fmn_mc=matrix(rnorm(nmc*nb_fmn, simfmn, 1./sqrt(invfmnvar)), nb_fmn, nmc);
   free_mc=matrix(0, nb_param, nmc);
   cost_mc=numeric(nmc);
   for (imc in 1:nmc) {
      cat("imc=", imc, "\n", sep="");
      # minimization
      res=opt_wrapper(meas_mc[,imc], fmn_mc[,imc]);
      # store the solution
      free_mc[,imc]=res$par;
      cost_mc[imc]=sum(jx_f$res*jx_f$res*invar);
   }
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
}

# Linear method based on jacobian x_f
# reset fluxes and jacobians according to param
cost=cumo_cost(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo");
#grj=cumo_gradj(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs");
x_f=v$x_f;
dimnames(x_f)=list(nm_cumo, names(fwrv));
if (DEBUG) {
   library(numDeriv); # for numerical jacobian
   # numerical simulation
   rj=function(v, ...) {
      r=cumo_resid(v, ...);
      c(r$res);
   }
   #dr_dpn=jacobian(rj, param, method="Richardson", method.args=list(),
   #   nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, "fwrv2Abcumo", "frj_rhs");
   dx_dfn=num_jacob(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
   dr_dp=num_fjacob(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs");
   # to compare with jx_f$dr_dp
   gr=grad(cumo_cost, param, method="Richardson", method.args=list(), nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo");
   grn=cumo_grad(param, nb_f, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo", fj_rhs=NULL);
   browser();
}

# covariance matrix of free fluxes
invcov=t(jx_f$dr_dff)%*%(jx_f$dr_dff*c(measinvvar, invfmnvar));
covff=try(solve(invcov));
if (inherits(covff, "try-error")) {
   # matrix seems to be singular
   cat("Inverse of covariance matrix is singular => there are unsolved fluxes\\n",
      file=stderr());
   # regularize and inverse
   covff=solve(invcov+diag(1.e-4, nrow(invcov)));
}
# standart deviations of free fluxes
sdff=sqrt(diag(covff));
cat("linear stats\\n", file=fkvh);
obj2kvh(covff, "covariance free fluxes", fkvh, indent=1);
covf=jx_f$df_dff%*%covff%*%t(jx_f$df_dff);
sdf=sqrt(diag(covf));
obj2kvh(covf, "covariance all fluxes", fkvh, indent=1);
obj2kvh(sdff, "SD free fluxes", fkvh, indent=1);
mtmp=cbind(Val=param[1:nb_ff], SD=sdff, CV=sdff/param[1:nb_ff]);
o=order(mtmp[,"CV"]);
obj2kvh(mtmp[o,], "Val, SD, CV free fluxes", fkvh, indent=1);
obj2kvh(sdf, "SD all fluxes", fkvh, indent=1);
# select best defined flux combinations
s=svd(covf);
# lin comb coeff matrix
l=t(s$u);
iord=apply(l,1,function(v){o=order(abs(v), decreasing=T); n=which(cumsum(v[o]**2)>=0.999); o[1:n[1]]})
dimnames(l)=dimnames(covf);
nm_fvrw=dimnames(l)[[1]];
nb_fvrw=nrow(s$u);
comb=rep("", nb_fvrw);
for (j in 1:nb_fvrw) {
   cva=l[j,iord[[j]]]; # abs(factor) ordered j-th lin combination
   cva=sign(cva[1])*cva;
   cva=sprintf("%+.3g*", cva);
   cva[cva=="+1*"]="+";
   cva[cva=="-1*"]="-";
   cva[1]=substring(cva[1], 2);
   comb[j]=paste(cva, nm_fvrw[iord[[j]]], sep="", collapse="");
}
# from best to worst defined sd
sta=sqrt(abs(diag(l%*%covf%*%t(l))));
names(sta)=comb;
o=order(sta);
sta=sta[o];
obj2kvh(sta, "SD for uncorrelated fw-rv linear combinations", fkvh, indent=1);
if (prof) {
   Rprof(NULL);
}
close(fkvh);
""");

f.close();
ff.close();
# make output files just readable to avoid later casual edition
os.system("chmod a-w '%s'"%n_R);
os.system("chmod a-w '%s'"%n_fort);
