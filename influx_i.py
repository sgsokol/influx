#!/usr/bin/env python
#
# Transform an ftbl to R code simulating a dynamique of
# isotopomer propagation. Flux values and metabolite pool sizes
# are taken from .ftbl file (dependent fluxes are calculated from
# free and constrained fluxes)

# usage: ./influx_i.py network_name[.ftbl]
# After execution files network_name.R and network_name.f
# will be created. If they already exist, they will be silently overwritten.
# 2009-03-18 sokol

# Important python variables:
# Collections:
#    netan - (dict) ftbl structured content;
#    tflcnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
#    collection;
#    measures - (dict) exp data;
#    scale - unique scale names;
#    nrow - counts scale names;
#    o_sc - ordered scale names;
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
# Vector names:
#    cumos (list) - names of R:x;

# Important R variables:
# Scalars:
#    nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
#    nb_ffn, nb_ffx, nb_ff (free fluxes),
#    nb_fcn, nb_fcx, nb_fc (constrained fluxes),
#    nb_ineq, nb_param, nb_fmn
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
#    fmn
# Matrices:
#    Afl, qrAfl, invAfl,
#    p2bfl - helps to construct the rhs of flux system
#    mf, md - help to construct flcnx
# Functions:
#    param2fl - translate param to fluxes

import sys;
import os;
import time;
import copy;

from tools_ssg import *;
import C13_ftbl;

me=os.path.basename(sys.argv[0]);
def usage():
    sys.stderr.write("usage: "+me+" organism");

#<--skip in interactive session
if len(sys.argv) < 2:
    usage();
    sys.exit(1);

# set some python constants
org=sys.argv[1] if sys.argv[1][-5:] != ".ftbl" else sys.argv[1][:-5];
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
n_R=org+".R";
n_fort=org+".f";
f_ftbl=open(n_ftbl, "r");
os.chmod(n_R,755);
os.chmod(n_fort,755);
f=open(n_R, "w");
ff=open(n_fort, "w");

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();

# analyse network
# reload(C13_ftbl);

netan=C13_ftbl.ftbl_netan(ftbl);

# write initialization part of R code
# flux part (Afl, bfl, ...)
ftbl2code.netan2R_fl(netan, org, f, ff);
# cumomer part (A, b, ...)
ftbl2code.netan2R_cumo(netan, org, f, ff);
ftbl2code.netan2R_ineq(netan, org, f, ff);
# mesure, rcumo system are skipped

f.write("""
# check if inequalities are satisfied
if (! all(ui%*%param>=ci)) {
   cat("The following inequalities are not satisfied. You should change 
initial values of free fluxes in ftbl file.\\n", file=stderr());
   cat(nm_i[ui%*%param<ci], sep="\\n", file=stderr());
}
""");
f.write("""
plotx=function(x, plottype, ti, ...) {
   # prepare and store data for plot
   # x is a full cumomer vector
   # ti is time point (if "row_col", write just name's row)
   # plottype is one of pos_enrich, mass, labeled
   # ... (e.g. f and indent) are passed to obj2kvh()
   plottype=match.arg(plottype, c("pos_enrich", "mass", "labeled"));
   if (plottype=="pos_enrich") {
      # get only cumomers with 1 bit set (<=> power of 2)
      nms=names(x);
      nb=length(x);
      icumo=log2(as.numeric(unlist(strsplit(nms, ":", fixed=TRUE))[2*(1:nb)]));
      i=icumo==round(icumo);
      if (ti == "row_col") {
         obj2kvh(paste(nms[i], collapse="\t"), ti, ...);
      } else {
         obj2kvh(paste(x[i], collapse="\t"), ti, ...);
      }
   } else if (plottype=="mass"){
      # convert cumomers to mass distrib vector
      mass=cumo2mass(x);
      if (ti == "row_col") {
         obj2kvh(paste(names(mass), collapse="\t"), ti, ...);
      } else {
         obj2kvh(paste(mass, collapse="\t"), ti, ...);
      }
   } else if (plottype=="labeled"){
      # convert cumomers to labeled fractions
      lab=cumo2lab(x);
      if (ti == "row_col") {
         obj2kvh(paste(names(lab), collapse="\t"), ti, ...);
      } else {
         obj2kvh(paste(lab, collapse="\t"), ti, ...);
      }
   }
}
## input cumomers xinp
xinp=c(%(xinp)s);
nm_xinp=c(%(nm_xinp)s);
names(xinp)=nm_xinp;
nm2=matrix(unlist(strsplit(nm_xinp, ":", fixed=TRUE)), ncol=2, byrow=TRUE);
o=order(nm2[,1], as.numeric(nm2[,2]));
xinp=xinp[o];

## variables for isotopomer cinetics
tstart=0.;
tmax=%(tmax)f;
dt=%(dt)f;
metab_scale=%(metab_scale)f;
# Metabolite pool vector numbered and replicated to match cumomers x
# and scaled by metab_scale/dt
Metab_dt=(metab_scale/dt)*c(%(met_pools)s);

# get full flux vectors
lf=param2fl(param, nb_f, invAfl, p2bfl, bp, fc);

# cumulated sum of nb_cumo by weight
nbc_cumos=c(0, cumsum(nb_cumos));

# prepare the cumosys cascade
Acumot=list();
for (iw in 1:nb_w) {
   ncumow=nb_cumos[iw];
   A=matrix(0.,ncumow,ncumow);
   #fwrv2Abcumo(fl, nf, x, nx, iw, n, A, b)
   res<-.Fortran("fwrv2Abcumo",
      fl=as.double(lf$fwrv),
      nf=length(lf$fwrv),
      x=as.double(0.),
      nx=as.integer(0),
      iw=as.integer(iw),
      n=as.integer(ncumow),
      A=as.matrix(A),
      b=as.double(0.),
      calcA=as.integer(TRUE),
      calcb=as.integer(FALSE),
      NAOK=TRUE,
      DUP=FALSE);
   # inverse A*x=b;
   Acumot=append(Acumot, list(solve(A+diag(Metab_dt[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], nrow=ncumow))));
}

# formated output in kvh file
fkvh=file("%(fkvh)s", "w");
fplot=file("%(fplot)s", "w");
"""%{
    "org": org,
    "fkvh": escape("%s_ires.kvh"%org, "\\"),
    "fplot": escape("%s_ipl.kvh"%org, "\\"),
    "dt": netan["opt"]["dt"],
    "tmax": netan["opt"]["tmax"],
    "metab_scale": netan["opt"]["metab_scale"],
    "met_pools": join(", ", (netan["met_pools"][cumo.split(":")[0]]
        for w in netan["vcumo"]
        for cumo in w)),
    "xinp": join(", ", netan["cumo_input"].values()),
    "nm_xinp": join(", ", netan["cumo_input"].keys(), '"', '"'),
});
# main part: ode solver
f.write("""
# initial flux and cumomer distribution
cat("flux parameters\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "free parameters", fkvh, ident=1);

fwrv=lf$fwrv;
n=length(fwrv);
names(fwrv)=paste(nm_fwrv, c(rep("fwd", n/2), rep("rev", n/2)), sep=".");
obj2kvh(fwrv, "fwd-rev flux vector", fkvh, ident=1);

f=lf$flcnx;
n=length(f);
names(f)=nm_flcnx;
obj2kvh(f, "net-xch flux vector", fkvh, ident=1);

# alphabetic order of output cumomers
nm2=matrix(unlist(strsplit(nm_cumo, ":", fixed=TRUE)), ncol=2, byrow=TRUE);
o_acumo=order(nm2[,1], as.numeric(nm2[,2]));
# plot options
plottype="%(plottype)s";
plotby=%(plotby)d;
""" % {
   "plottype": netan["opt"].get("plottype", "none"),
   "plotby": netan["opt"].get("plotby", 1),
});
f.write("""
# time 0 cumomers (are all unlabeled)
xold=rep(0., sum(nb_cumos));
names(xold)=nm_cumo;
cat("time course cumomers\n", file=fkvh);
cat("\trow_col", c(nm_xinp, nm_cumo[o_acumo]), file=fkvh, sep="\t");
cat("\n", file=fkvh);
obj2kvh(paste(c(xinp, xold[o_acumo]), collapse="\t"), "0.", fkvh, ident=1);
plotx(c(xinp,xold[o_acumo]), plottype, "row_col", fplot);
plotx(c(xinp,xold[o_acumo]), plottype, 0., fplot);

# time going (implicite Euler)
it=1;
for (ti in seq(dt, tmax, by=dt)) {
   # solve cumomer cascade at current time
   x=numeric(0);
   for (iw in 1:nb_w) {
      nx=length(x);
      # prepare rhs for current weight
      ncumow=nb_cumos[iw];
      b=double(ncumow);
      res<-.Fortran("fwrv2Abcumo",
         fl=as.double(lf$fwrv),
         nf=length(lf$fwrv),
         x=as.double(x),
         nx=as.integer(nx),
         iw=as.integer(iw),
         n=as.integer(ncumow),
         A=as.matrix(0.),
         b=as.double(b),
         calcA=as.integer(FALSE),
         calcb=as.integer(TRUE),
         NAOK=TRUE,
         DUP=FALSE);
      x=c(x, Acumot[[iw]]%*%(b+Metab_dt[(nbc_cumos[iw]+1):nbc_cumos[iw+1]]*xold[(nbc_cumos[iw]+1):nbc_cumos[iw+1]]));
   }
   # store this time step result
   obj2kvh(paste(c(xinp,x[o_acumo]), collapse="\t"), ti, fkvh, ident=1);
   
   # data to plot
   if (!it%%plotby) {
      names(x)=nm_cumo;
      plotx(c(xinp,x[o_acumo]), plottype, ti, fplot);
   }
   
   xold=x;
   it=it+1;
}
""");
f.close();
ff.close();
# make output files just readable to avoid later casual edition
os.chmod(n_R, 666);
os.chmod(n_fort, 666);
