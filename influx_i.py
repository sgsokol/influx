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
#    tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
#    collection;
#    measures - (dict) exp data;
#    scale - unique scale names;
#    nrow - counts scale names;
#    o_sc - ordered scale names;
# org - (str) prefix of .ftbl  file like "PPP"
# File names (str):
#    n_ftbl (descriptor f_ftbl);
#    n_R (R code) (f);
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
#    fmn
# Matrices:
#    Afl, qrAfl, invAfl,
#    p2bfl - helps to construct the rhs of flux system
#    mf, md - help to construct fallnx
# Functions:
#    param2fl - translate param to fluxes

import sys;
import os;
import stat;
#import time;
#import copy;

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
f_ftbl=open(n_ftbl, "r");
try:
    os.chmod(n_R, stat.S_IRUSR | stat.S_IWUSR);
except:
    pass;
f=open(n_R, "w");

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();

# analyse network
# reload(C13_ftbl);

netan=C13_ftbl.ftbl_netan(ftbl);

# write initialization part of R code
ftbl2code.netan2Rinit(netan, org, f, True);
# flux part (Afl, bfl, ...)
ftbl2code.netan2R_fl(netan, org, f);
# cumomer part (A, b, ...)
ftbl2code.netan2R_cumo(netan, org, f);
ftbl2code.netan2R_ineq(netan, org, f);
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
## input cumomers xi
xi=c(%(xi)s);
nm_xi=c(%(nm_xi)s);
names(xi)=nm_xi;
nm2=matrix(unlist(strsplit(nm_xi, ":", fixed=TRUE)), ncol=2, byrow=TRUE);
o=order(nm2[,1], as.numeric(nm2[,2]));
xi=xi[o];

## variables for isotopomer cinetics
tstart=0.;
tmax=%(tmax)f;
dt=%(dt)f;
metab_scale=%(metab_scale)f;
# Metabolite pool vector numbered and replicated to match cumomers x
# and scaled by metab_scale/dt
Metab_dt=(metab_scale/dt)*c(%(met_pools)s);
#browser()
# get full flux vectors
lf=param2fl(param, nb_f, invAfl, p2bfl, bp, fc);

# cumulated sum of nb_cumo by weight
nbc_cumos=c(0, cumsum(nb_cumos));

# prepare the cumosys cascade
Acumot=list();
nb_xi=length(xi);
x=c(1, xi);
for (iw in 1:nb_w) {
   ncumow=nb_cumos[iw];
   lAb=fwrv2Ab(lf$fwrv, spAb_f[[iw]], x);
   A=lAb$A;
   b=lAb$b;
   # inverse A*x=b;
   Acumot=append(Acumot, list(as.matrix(solve(A+diag(Metab_dt[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], nrow=ncumow)))));
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
    "xi": join(", ", netan["cumo_input"].values()),
    "nm_xi": join(", ", netan["cumo_input"].keys(), '"', '"'),
});
# main part: ode solver
f.write("""
# initial flux and cumomer distribution
cat("flux parameters\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "free parameters", fkvh, indent=1);

fwrv=lf$fwrv;
n=length(fwrv);
names(fwrv)=paste(nm_fwrv, c(rep("fwd", n/2), rep("rev", n/2)), sep=".");
obj2kvh(fwrv, "fwd-rev flux vector", fkvh, indent=1);

f=lf$fallnx;
n=length(f);
names(f)=nm_fallnx;
obj2kvh(f, "net-xch flux vector", fkvh, indent=1);

# alphabetic order of output cumomers
nm_incuf=c("one", nm_xi, nm_cumo); # full cumomer list names with input names
nm2=matrix(unlist(strsplit(nm_incuf[-1], ":", fixed=TRUE)), ncol=2, byrow=TRUE);
o_acumo=order(nm2[,1], as.numeric(nm2[,2]));
o_acumo=c(0, o_acumo+1); # 0 is to exclude the leading "1" in incu vector
# plot options
plottype="%(plottype)s";
plotby=%(plotby)d;
""" % {
   "plottype": netan["opt"].get("plottype", "none"),
   "plotby": netan["opt"].get("plotby", 1),
});
f.write("""
# time 0 cumomers (are all unlabeled)
#browser();
xold=rep(0., sum(nb_cumos)+nb_xi+1);
names(xold)=nm_incuf;
cat("time course cumomers\\n", file=fkvh);
cat("\\trow_col", nm_incuf[o_acumo], file=fkvh, sep="\\t");
cat("\\n", file=fkvh);
obj2kvh(paste(xold[o_acumo], collapse="\\t"), "0.", fkvh, indent=1);
plotx(xold[o_acumo], plottype, "row_col", fplot);
plotx(xold[o_acumo], plottype, 0., fplot);

# time going (implicite Euler)
it=1;
for (ti in seq(dt, tmax, by=dt)) {
   # solve cumomer cascade at current time
   x=c(1., xi);
   for (iw in 1:nb_w) {
      nx=length(x);
      # prepare rhs for current weight
      ncumow=nb_cumos[iw];
      fb=rep(0., ncumow);
      b=rep(0., ncumow);
      res<-.Fortran("f2b",
         bfpr=as.matrix(spAb_f[[iw]]$bfpr),
         n=as.integer(ncol(spAb_f[[iw]]$bfpr)),
         fwrv=as.double(fwrv),
         incu=as.double(x),
         fb=as.double(fb),
         b=as.double(b),
         NAOK=TRUE,
         DUP=FALSE
      );
      b=-b;
      x=c(x, Acumot[[iw]]%*%(b+
         Metab_dt[(nbc_cumos[iw]+1):nbc_cumos[iw+1]]*
         xold[1+nb_xi+((nbc_cumos[iw]+1):nbc_cumos[iw+1])]));
   }
   # store this time step result
   obj2kvh(paste(x[o_acumo], collapse="\t"), ti, fkvh, indent=1);
   
   # data to plot
   if (!it%%plotby) {
      names(x)=nm_incuf;
      plotx(x[o_acumo], plottype, ti, fplot);
   }
   
   xold=x;
   it=it+1;
}
""");
f.close();
# make output files just readable to avoid later casual edition
os.chmod(n_R, stat.S_IRUSR);
