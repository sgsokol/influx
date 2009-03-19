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
#    tflcnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
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
#    n_opt (R code) (f);
#    n_fort (fortran code) (ff);
# Counts: no_fln, no_flx, no_fl (dependent fluxes: net, xch, total),
#         no_ffn, no_ffx (free fluxes)
# Index translators:
#    fwrv2i - flux names to index in R:fwrv;
#    cumo2i - cumomer names to index in R:x;
#    ir2isc - mapping measure rows indexes on scale index isc[meas]=ir2isc[meas][ir]
# Vector names:
#    cumos (list) - names of R:x;
#    o_mcumos - cumomers involved in measures;

# Important R variables:
# Scalars:
#    no_w, no_cumos, no_fln, no_flx, no_fl (dependent or unknown fluxes),
#    no_ffn, no_ffx, no_ff (free fluxes),
#    no_fcn, no_fcx, no_fc (constrained fluxes),
#    no_ineq, no_param, no_fmn
# Name vectors:
#    nm_cumo, nm_fwrv, nm_flcnx, nm_fln, nm_flx, nm_fl, nm_par,
#    nm_ffn, nm_ffx,
#    nm_fcn, nm_fcx,
#    nm_mcumo, nm_fmn
# Numeric vectors:
#    fwrv - all fluxes (fwd+rev);
#    x - all cumomers (weight1+weight2+...);
#    param - free flux net, free flux xch, scale label, scale mass, scale peak
#    fcn, fcx, fc,
#    bp - helps to construct the rhs of flux system
#    flcnx - complete flux vector (constr+net+xch)
#    bc - helps to construct flcnx
#    li - inequality vector (mi%*%flcnx>=li)
#    ir2isc - measur row to scale vector replicator
#    ci - inequalities for param use (ui%*%param-ci>=0)
#    measvec,
#    measinvvar,
#    imeas,
#    fmn
# Matrices:
#    Afl, qrAfl, invAfl,
#    p2bfl - helps to construct the rhs of flux system
#    mf, md - help to construct flcnx
#    mi - inequality matrix (ftbl content)
#    ui - inequality matrix (ready for param use)
#    measmat - measmat*(x[imeas];1)=vec of simulated not-yet-scaled measures
# Functions:
#    param2fl_x - translate param to flux and cumomer vector (initial approximation)
#    cumo_cost - cost function (khi2)
#    cumo_grad - finite difference gradient

import sys;
import os;
import time;
import copy;

#sys.path.append("/home/sokol/dev/python");
from tools_ssg import *;
import C13_ftbl;

me=os.path.basename(sys.argv[0]);
def usage():
    sys.stderr.write("usage: "+me+" network_name[.ftbl]");

#<--skip in interactive session
if len(sys.argv) < 2:
    usage();
    sys.exit(1);

# set some python constants
org=sys.argv[1];
# cut .ftbl if any
if org[-5:]==".ftbl":
    org=org[-5:];

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
n_opt=org+".R";
n_fort=org+".f";
f_ftbl=open(n_ftbl, "r");
os.system("chmod u+w '%s' 2>/dev/null"%n_ftbl);
os.system("chmod u+w '%s' 2>/dev/null"%n_fort);
f=open(n_opt, "w");
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
if (no_fln) {
   print(paste(nm_fln,"net",sep="_"));
}
if (no_flx) {
   print(paste(nm_flx,"xch",sep="_"));
}
cat("Free fluxes:\\n");
if (no_ffn) {
   print(paste(nm_ffn,"net",sep="_"));
}
if (no_ffx) {
   print(paste(nm_ffx,"xch",sep="_"));
}
cat("Constrained fluxes:\\n");
if (no_fcn) {
   print(paste(nm_fcn,"net",sep="_"));
}
if (no_fcx) {
   print(paste(nm_fcx,"xch",sep="_"));
}
"""
f.write("""
# set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
# for corresponding measures
vr=param2fl_x(param, no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
simvec=(measmat%*%c(vr$x[irmeas],1.));
if (DEBUG) {
   cat("initial simvec:\\n");
   print(simvec);
}
if (no_ff < length(param)) {
   ms=measvec*simvec*measinvvar;
   ss=simvec*simvec*measinvvar;
   for (i in (no_ff+1):length(param)) {
      im=(ir2isc==(i+1));
      param[i]=sum(ms[im])/sum(ss[im]);
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
# get initial flux and cumomer distribution
cat("initial approximation\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "free parameters", fkvh, ident=1);

x=vr$x;
names(x)=nm_rcumo;
obj2kvh(x, "starting cumomer vector", fkvh, ident=1);

fwrv=vr$fwrv;
n=length(fwrv);
names(fwrv)=paste(nm_fwrv, c(rep("fwd", n/2), rep("rev", n/2)), sep=".");
obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, ident=1);

f=vr$flcnx;
n=length(f);
names(f)=nm_flcnx;
obj2kvh(f, "starting net-xch flux vector", fkvh, ident=1);

rres=cumo_resid(param, no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
obj2kvh(rres$res, "starting cumomer residuals", fkvh, ident=1);

obj2kvh(rres$flcnx[ifmn]-fmn, "flux residual vector", fkvh, ident=1);

rcost=cumo_cost(param, no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo");
obj2kvh(rcost, "starting cost value", fkvh, ident=1);

obj2kvh(Afl, "flux system (Afl)", fkvh, ident=1);
obj2kvh(p2bfl%*%param[1:no_f$no_ff]+bp, "flux system (bfl)", fkvh, ident=1);

#cat("mass vector:\\n");
#print_mass(x);

f=vr$fwrv;
n=length(f);
names(f)=nm_fwrv;
obj2kvh(f, "fwd-rev flux vector", fkvh, ident=1);

# optimize all this
names(param)=nm_par;
if (method == "BFGS") {
   control=list(maxit=500, trace=1);
   res=constrOptim(param, cumo_cost, grad=cumo_grad,
      ui, ci, mu = 1e-04, control,
      method="BFGS", outer.iterations=100, outer.eps=1e-07,
      no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc,
      irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo");
} else if (method == "Nelder-Mead") {
   control=list(maxit=1000, trace=1);
   res=constrOptim(param, cumo_cost, grad=cumo_grad,
      ui, ci, mu = 1e-04, control,
      method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
      no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc,
      irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, "fwrv2rAbcumo");
} else {
   stop(paste("Unknown minimization method '", method, "'", sep=""));
}
param=res$par;
names(param)=nm_par;
""");
f.write("""
obj2kvh(res, "optimization process stats", fkvh);

rres=cumo_resid(param, no_f, no_rw, no_rcumos, invAfl, p2bfl, bp, fc, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
obj2kvh(rres$res, "cumomer residual vector", fkvh);
obj2kvh(rres$flcnx[ifmn]-fmn, "flux residual vector", fkvh);
obj2kvh(measvec, "cumomer measure vector", fkvh);

v=param2fl_x(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, "fwrv2Abcumo");
x=v$x;
names(x)=nm_cumo;
o=order(nm_cumo);
obj2kvh(x[o], "cumomer vector", fkvh);

fwrv=v$fwrv;
n=length(fwrv);
names(fwrv)=paste(nm_fwrv, c(rep("fwd", n/2), rep("rev", n/2)), sep=".");
obj2kvh(fwrv, "fwd-rev flux vector", fkvh);

f=v$flcnx;
n=length(f);
names(f)=nm_flcnx;
obj2kvh(f, "net-xch flux vector", fkvh);
close(fkvh);

if (sensitive=="grad") {
   # sensitivity analysis
   # sensit=df_i/dm_j (jacobian of solution (f) depending on measures (m)
   # perturb mesures 1 by 1 by factor of 1.+pfact and
   # store soltions as columns in sensit
   pfact=0.1;
   sensit=matrix(0, length(nm_fwrv), 0);
   for (i in 1:length(measvec)) {
      # prepare perturbed measures
      measpert=measvec;
      dv=measvec[i]*pfact;
      measpert[i]=measvec[i]+dv;
      # solve perturbed problem
      if (method == "BFGS") {
         control=list(maxit=500, trace=1)
         res=constrOptim(param, cumo_cost, grad=cumo_grad,
            ui, ci, mu = 1e-04, control,
            method="BFGS", outer.iterations=10, outer.eps=1e-05,
            no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
            imeas, measmat, measpert, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      } else if (method == "Nelder-Mead") {
         control=list(maxit=1000, trace=1);
         res=constrOptim(param, cumo_cost, grad=cumo_grad,
            ui, ci, mu = 1e-04, control,
            method="Nelder-Mead", outer.iterations=100, outer.eps=1e-05,
            no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
            imeas, measmat, measpert, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      } else if (method == "SANN") {
         control=list(maxit=10000, trace=1)
         res=constrOptim(param, cumo_cost, grad=cumo_grad,
            ui, ci, mu = 1e-04, control,
            method="SANN", outer.iterations=100, outer.eps=1e-05,
            no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
            imeas, measmat, measpert, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      } else {
         stop(paste("Unknown minimization method '", method, "'", sep=""));
      }
   #print(dv);
   #print(res);
      # store perturbed solution
      v=param2fl_x(res$par, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measpert, ir2isc);
   cat("Perturbed free fluxes:\\n");
   print(res$par);
   cat("Perturbed fluxes:\\n");
   print(v$fwrv);
      sensit=cbind(sensit, (v$fwrv-fwrv)/dv);
   }
   dimnames(sensit)[[1]]=paste(nm_fwrv,c(rep("fwd",length(nm_fwrv)/2),rep("rev",length(nm_fwrv)/2)),sep=".");
   # SD vector for fluxes
   fl_sd=sqrt((sensit**2)%*%(1./measinvvar));
   #names(fl_sd)=dimnames(sensit)[[1]];

   cat("sensitivity matrix:\\n");
   print(sensit);

   cat("fwd-rev flux Standard Deviation (SD):\\n");
   cat(paste(dimnames(sensit)[[1]], fwrv, rep("+-", length(fwrv)), fl_sd), sep="\\n");
} else if (sensitive=="mo") {
   # Monte-Carlo simulation
   nmc=10; # generated measure sample number
   # random measure generation
   no_meas=length(measvec);
   meas_mc=matrix(rnorm(nmc*no_meas, measvec, 1./sqrt(measinvvar)), no_meas, nmc);
   free_mc=matrix(0, no_param, 0);
   for (imc in 1:nmc) {
      print(paste("imc=",imc,sep=""));
      # minimization
      if (method == "BFGS") {
         control=list(maxit=500, trace=0);
         res=constrOptim(param, cumo_cost, grad=cumo_grad,
            ui, ci, mu = 1e-04, control,
            method="BFGS", outer.iterations=10, outer.eps=1e-05,
            no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
            imeas, measmat, meas_mc[,imc], measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      } else if (method == "Nelder-Mead") {
         control=list(maxit=1000, trace=0);
         res=constrOptim(param, cumo_cost, grad=cumo_grad,
            ui, ci, mu = 1e-04, control,
            method="Nelder-Mead", outer.iterations=100, outer.eps=1e-05,
            no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
            imeas, measmat, meas_mc[,imc], measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      } else {
         stop(paste("Unknown minimization method '", method, "'", sep=""));
      }
      # store the solution
      free_mc=cbind(free_mc, res$par);
   }
   dimnames(free_mc)[[1]]=names(param);
""");
f.write("""
   write.table(meas_mc, file="%(org)s_mmc.txt");
   write.table(free_mc, file="%(org)s_fmc.txt");
}
if (prof) {
   Rprof(NULL);
}
""" % {
   "org": org,
});

f.close();
ff.close();
# make output files just readable to avoid later casual edition
os.system("chmod a-w '%s'"%n_ftbl);
os.system("chmod a-w '%s'"%n_fort);
