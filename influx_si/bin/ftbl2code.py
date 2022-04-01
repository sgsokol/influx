#!/usr/bin/env python3

"""Module for translation of .ftbl file to R code"""

# 2012-02-21 sokol@insa-toulouse.fr : cumomer matrices and rhs from sparse matrices
#                                     (without fortran code)
# 2009-09-14 sokol@insa-toulouse.fr : flux.[net|xch] -> [dfcg].[nx].flux
#                                     flux.[fwd|rev] -> [fwd|rev].flux
# 2008-12-08 sokol@insa-toulouse.fr : added netan2Rinit()
# 2008-11-25 sokol@insa-toulouse.fr : adapted for reduced cumomer list
# 2008-09-19 sokol@insa-toulouse.fr : initial release
# Copyright 2011-2020, INRAE

import time
import copy
import os
import sys
#import pdb
from operator import itemgetter
from itertools import groupby

me=os.path.abspath(os.path.realpath(sys.argv[0]))
dirx=os.path.dirname(me)
sys.path.append(dirx)
if (dirx.endswith("py3")):
    dirx=os.path.split(dirx)[0]

import influx_si
dirr=os.path.join(os.path.dirname(os.path.realpath(influx_si.__file__)), "R")

from tools_ssg import *
import C13_ftbl

def netan2Abcumo_spr(varname, Al, bl, vcumol, minput, f, fwrv2i, incu2i_b1):
    """
    Transform cumomer linear sytems collection (from ftbl file)
    to a R code calculating sparse matrix A and vector b
    in A*x+b=0 for a given weight of fragment iw (index in resulting list)
    Flux vector fl of all fwd. and rev. fluxes are known at R runtime.
    
    Resulting code is a list sprAb indexed by cumomer weight
    (cf. generated R comments for details on sprAb)
    cumomer vector incu=c(1, xi, xl), xi - input cumomers, xl - lighter cumomers.
    
    incu2i_b1 gives i in incu from cumomer name. i=1 corresponds to the constant 1.
    """
    #2012-02-08 sokol
    #2016-09-23 sokol: any number of fused fragments in b (not limited to 2 as before)
    
    nb_cumu=cumsum(len(l) for l in vcumol)
    f.write(
    """
# sparse matrix static parts
# $varname fields:
#  ind_fa - flux index in a_pre$vfwrv[ind_fa]
#  a_pre - sparse matrix whose colsum() gives the a$v vector
#  prodx - dense matrix whose colprod() will give x[ind_x1]*x[ind_x2]*...
#  ind_fb - flux index in b_pre$v=fwrv[ind_fb1]*colprod(prodx)
#  ind_b - dense matrix of indexes for  b_pre$v=f[ind_b[,"indf"]*x[ind_b[,2+1]]*x[ind_b[,2+2]], ...]
#  b_pre - sparse matrix whose colsum gives b@x

#  a - unsigned sparse cumomer A matrix (off-diagonal part)
#  b - unsigned sparse vector of right hand side

if (TIMEIT) {
   cat("spAbr   : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}

nb_fwrv=%(n)d
nb_w=%(nb_w)d
%(var)s=list()
"""%{
    "var": varname,
    "n": len(fwrv2i),
    "nb_w": len(Al),
})
    # base of cumomers in composed vector incu=c(1, input, xcumo)
    # +1 for c(1,...)
    ba_x=len(incu2i_b1) - sum(len(l) for l in vcumol)+1
    ba_xw=ba_x; # base for current weigth cumomer in incu
    ncucumo=0
    for (iwl,A) in enumerate(Al):
        w=iwl+1
        b=bl[iwl]
        cumos=vcumol[iwl]
        ncumo=len(cumos)
        c2i=dict((c,i) for (i,c) in enumerate(cumos))
        #d=[c for c in netan['cumo_sys']['A'][w-1] if not c in cumos]
        if ncumo != len(A):
            raise Exception("wrongCumomerNumber: ncumo=%d, nrow(A)=%d"%(ncumo, len(A)))
        l_ia=[]; # list of non zero off-diagonal elements in A / row
        l_ib=[]; # list of non zero elements in b / row
        nb_maxfa=0; # how many fluxes in an off-diagonal term in a
        nb_maxprod=0 if ncumo == 0 else max(len(li) for cu,rdi in b.items() for fl,d in rdi.items() for i,li in d.items()); # how many cumomer fragments are fused in b
        for irow in range(ncumo):
            cr=cumos[irow]
            row=A[cr]
            # atuple is list of (icumo, list(fluxes))
            atuple=[(c2i[c], [fwrv2i[fl] for fl in row[c]])
                for c in row] #cumos if c in row and c!=cr]
            #if atuple:
            #    nb_maxfa=max(nb_maxfa, max(len(lf) for (ic, lf) in atuple))
            #elif cr not in b:
            #    raise Exception("Empty row in cumomer matrix, weight=%d (base 1), cumo=%s"%(w, cr))
            # btuple is list of [iflux, [icumo1, icumo2, icumo_i,...]]
            if cr in b:
                btuple=[[fwrv2i[fl], [incu2i_b1[v] for v in l]+[1]*(nb_maxprod-len(l))]
                    for (fl, d) in b[cr].items()
                    for (i,l) in d.items()]
                #nb_maxfb=max(nb_maxfb, len(btuple))
            else:
                btuple=[]
            # one list per row
            l_ia.append(atuple)
            #nb_ax+=len(atuple)
            l_ib.append(btuple)
        #print("w=", w, "A=", A, "l_ia=", l_ia, "\n")
        f.write(
"""
if (TIMEIT) {
   cat("weight %(w)d: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
w=%(w)d
nb_c=%(nbc)d
ba_x=%(ba_x)d; # base of cumomer indexes in incu vector
l=new.env()
l$w=w
l$nb_c=nb_c
l$nb_fwrv=nb_fwrv
l$nb_cl=%(ncucumo)d # number of lighter cumomers
maxprod=%(maxprod)d
if (nb_c > 0) {
   # matrix a
   ind_a=matrix(as.integer(c(%(ind_a)s)), ncol=3, byrow=TRUE)
   colnames(ind_a)=c("indf", "ir0", "ic0")
   l$ind_a=ind_a
   
   # vector b
   ind_b=matrix(as.integer(c(%(ind_b)s)), ncol=2+%(maxprod)d, byrow=TRUE)
   colnames(ind_b)=c("indf", "irow", paste("indx", seq_len(%(maxprod)d), sep=""))
   l$ind_b=ind_b
   
   # jacobian b_x
   imaxprod=seq_len(maxprod)
   ind_bx=c()
   for (ix in imaxprod) {
      i=ind_b[,2+ix]>ba_x # exclude from differentiation plain input entries
      tmp=ind_b[i,,drop=FALSE]
      ind_bx=rbind(ind_bx, tmp[,c(1,2,ix+2,2+imaxprod[-ix])]) # move diff var to ic1 place
   }
   if (length(ind_bx)) {
      colnames(ind_bx)=c("indf", "irow", "ic1", sprintf("indx%%d", seq_len(maxprod-1)))
      ind_bx[,"ic1"]=ind_bx[,"ic1"]-ba_x
   }
   l$ind_bx=ind_bx
}
%(var)s[[w]]=l
"""%{
   "var": varname,
   "w": w,
   "nbc": ncumo,
   "ncucumo": ncucumo,
   "ba_x": ba_x,
   "maxprod": nb_maxprod,
   "ind_a": join(", ", valval((ifl, ir, ic)
      for (ir, lt) in enumerate(l_ia)
      for (ic, lf) in lt
      for ifl in lf)),
   "ind_b": join(", ", valval((ifl, ir+1, ", ".join(str(i) for i in ii))
       for (ir, lt) in enumerate(l_ib)
       for (ifl, ii) in lt
   )),
})
        ba_xw+=ncumo
        ncucumo+=ncumo

def netan2Rinit(netan, org, f, fullsys, emu=False, ropts=[]):
    r"""Write R code for initialization of all variables before
cumomer system resolution by chi2 minimization.
:param netan: a collection of parsed ftbl information
:param f: R code output pointer
:param fullsys (logical): write a code for the full or only reduced cumomer system
:param emu (logical): write equations in EMU framework or cumomer (default)
:param ropts: list of items "param=value" to be written as is in R file.

:returns: a dictionnary with some python variables:
    * "measures": measures,
    * "o_mcumos": o_mcumos,
    * "cumo2i": cumo2i,
    * ...
    
"""
    # Important python variables:
    # Collections:
    #    netan - (dict) ftbl structured content
    #    tfallnx - (3-tuple[reac,["d"|"f"|"c"|"g"], ["n"|"x"]] list)- total flux
    #    collection
    #    measures - (dict) exp data
    #    rAb - (list) reduced linear systems A*x_cumo=b by weight
    #    scale - unique scale names
    #    nrow - counts scale names
    #    o_sc - ordered scale names
    #    o_meas - ordered measure types
    # org - (str) prefix of .ftbl  file like "PPP"
    # File names (str):
    #    n_ftbl (descriptor f_ftbl)
    #    n_opt (R code) (f)
    #    n_fort (fortran code) (ff)
    # Counts: nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total),
    #         nb_ffn, nb_ffx (free fluxes)
    # Index translators:
    #    fwrv2i - flux names to index in fwrv 1-based
    #    cumo2i - cumomer names to index in R:x
    #    ir2isc - mapping measure rows indexes on scale index isc[meas]=ir2isc[meas][ir]
    # Vector names:
    #    cumos (list) - names of R:x
    #    o_mcumos - cumomers involved in measures

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
    #    fwrv - all fluxes (fwd+rev)
    #    x - all cumomers (weight1+weight2+...)
    #    param - free flux net, free flux xch, scale label, scale mass, scale peak
    #    fcn, fcx, fc,
    #    bp - helps to construct the rhs of flux system
    #    fallnx - complete flux vector (dep, free, constr, growth:net+xch)
    #    bc - helps to construct fallnx
    #    li - inequality vector (mi%*%fallnx>=li)
    #    ir2isc - measur row to scale vector replicator
    #    ci - inequalities for param use (ui%*%param-ci>=0)
    #    measvec,
    #    measdev,
    #    fmn
    #    nb_sys - system sizes
    # Matrices:
    #    Afl, qrAfl, invAfl,
    #    p2bfl - helps to construct the rhs of flux system from free fluxes
    #    c2bfl - helps to construct the rhs of flux system from constr. fluxes
    #    mf, md, mc, mg - help to construct fallnx
    #    mi - inequality matrix (ftbl content)
    #    ui - inequality matrix (ready for param use)
    #    measmat - measmat*x+memaone=vec of simulated not-yet-pooled and not-yet-scaled measurements
    # Functions:
    #    lab_sim - translate param to flux and cumomer vector (initial approximation)
    #    cumo_cost - cost function (chi2)
    #    cumo_grad - finite difference gradient
    #    fallnx2fwrv - produce fw-rv fluxes from fallnx

    # Main steps:
    #    python var init
    #    R init
    #    R function fallnx2fwrv()
    #    python measures, cumos, cumo2i
    #    fortran code for cumomer systems A*x=b
    #    R var init
    #    R Afl, qr(Afl), invAfl
    #    R param (without scale factors)
    #    R constrained fluxes
    #    R p2bfl, c2bfl, bp
    #    R mf, md, mc
    #    R mi, li
    #    python measure matrix, vector and vars
    #    R ui, ci
    #    R measure matrix, vector, vars
    #    R flux measurements

    nexp=len(netan["iso_input"])
    # header
    f.write("# This is an automatically generated R code. Don't edit.\n")
    f.write("# Generated by \n# "+join(" ", sys.argv)+"\n# at "+time.ctime()+".\n")
    f.write("""
# Copyright 2011-%d, INRAE, France.
"""%time.localtime()[0])
    res=dict()
    ropts="\n".join(ropts)
    if ropts and ropts[0]=='"':
        ropts=ropts[1:-1]
    f.write("""
# working dir
dirw="%(dirw)s"

# installation dir (where influx_si/R/*.R live)
dirr="%(dirr)s"
# short base name of the FTBL (withount '.ftbl')
baseshort="%(org)s"

fcerr=file(file.path(dirw, sprintf("%%s.err", baseshort)), "ab")
fclog=file(file.path(dirw, sprintf("%%s.log", baseshort)), "ab")

if (options()$warn == 0)
    options(warn=1)
options(digits.secs=2)

case_i=%(case_i)s

if (length(find("bitwAnd"))==0L) {
   suppressPackageStartupMessages(library(bitops))
   bitwAnd=bitAnd
}
source(file.path(dirr, "libs.R"))

# define matprod for simple_triplet_matrix
`%%stm%%` = slam::matprod_simple_triplet_matrix

# default options
version=FALSE
noopt=FALSE
noscale=FALSE
meth="nlsic"
fullsys=FALSE
emu=FALSE
irand=FALSE
sens=""
cupx=0.999
cupn=1.e3
cupp=1.e5
clownr=0
cinout=0
clowp=1.e-8
np=0
ln=FALSE
tikhreg=FALSE
sln=FALSE
lim=FALSE
zc=-.Machine$double.xmax
ffguess=FALSE
fdfit=FALSE
addnoise=FALSE
fseries=""
iseries=""
seed=-.Machine$integer.max
excl_outliers=FALSE
TIMEIT=FALSE
prof=FALSE
time_order="1"

# get runtime arguments
%(ropts)s

# synonymous
myver=version
optimize=!noopt
methods=trimws(strsplit(meth, ",")[[1L]])
sensitive=sens
least_norm=ln
initrand=irand

vernum="%(vernum)s"

# sanity check for command line parameters
if (substring(sensitive, 1, 3)=="mc=") {
   # read the mc iteration number
   nmc=as.integer(substring(sensitive, 4))
   sensitive="mc"
} else if (sensitive=="mc") {
   nmc=10
}
# cupx==0 means no upper limit => cupx=1
cupx=ifelse(cupx, cupx, 1)
if (cupx < 0 || cupx > 1) {
   stop_mes("Option '--cupx N' must have N in the interval [0,1]\\n",
      "Instead, the value ", cupx, " is given.", file=fcerr)
}
if (cinout < 0) {
   stop_mes("Option '--cinout N' must have N non negative\\n",
      "Instead, the value ", cinout, " is given.", file=fcerr)
}
# minimization method
validmethods=c("BFGS", "Nelder-Mead", "SANN", "ipopt", "nlsic", "pso")
if (! all(igood <- (methods %%in%% validmethods))) {
   cat(paste("Warning: optimization methods ", paste0(methods[!igood], collapse=", "), " are not implemented. 'nlsic' is used instead."), "\\n", sep="", file=fcerr)
   methods[!igood]="nlsic"
}
if ("ipopt" %%in%% methods) {
   installed=suppressPackageStartupMessages(library(ipoptr, logical.return=TRUE))
   if (!installed) {
      stop_mes("An optimization method ipopt is requested but not available in this R installation", file=fcerr)
   }
}
if (least_norm && sln) {
   stop_mes("Options --ln and --sln cannot be activated simultaniously.", file=fcerr)
}

avaco=try(detectCores(), silent=TRUE)
if (inherits(avaco, "try-error")) {
   avaco=NULL
}
if (np > 0L && np < 1L) {
   np=round(avaco*np)
} else if (np >= 1L) {
   np=round(np)
} else {
   np=avaco
}
if (is.null(np) || np <= 0L) {
   np=1L
}
if (sensitive=="mc") {
   np=min(np, nmc)
}
options(mc.cores=np)

if (least_norm+tikhreg+lim > 1) {
   stop_mes("Options --ln, --lim and --tikhreg cannot be activated simultaneously. Use only one of them at a time.", file=fcerr)
}
lsi_fun=lsi
if (least_norm || sln) {
   lsi_fun=lsi_ln
} else if (tikhreg) {
   lsi_fun=lsi_reg
} else if (lim) {
   suppressPackageStartupMessages(library(limSolve));
   lsi_fun=lsi_lim
}
if (zc==-.Machine$double.xmax) {
   # no zero scrossing to apply
   zerocross=F
} else {
   if (zc < 0.) {
      stop_mes("Zero crossing value ZC must be non negative, instead ", zc, " is given.", file=fcerr)
   }
   zerocross=T
}
if (seed==-.Machine$integer.max) {
   # no seed to apply
   set_seed=F
} else {
   set_seed=T
   set.seed(seed)
}
time_order=gsub("\\\\s", "", time_order) # remove spaces if any
if (!(time_order %%in%% c("1", "2", "1,2"))) {
   stop_mes("time_order must be '1', '2' or '1,2'. Instead got '", time_order, "'", file=fcerr)
}
opts=commandArgs()
# end command line argument proceeding

# get some cumomer tools
source(file.path(dirr, "opt_cumo_tools.R"))
#loadcmp(file.path(dirr, "opt_cumo_tools.Rc"))

lab_resid=cumo_resid
lab_sim=param2fl_x
jx_f=new.env()
"""%{
    "dirw": escape(os.path.abspath(os.path.dirname(f.name)), '\\"'),
    "dirr": escape(dirr, '\\"'),
    "case_i": "TRUE" if case_i else "FALSE",
    "vernum": open(os.path.join(dirr, "..", "influx_version.txt"), "r").read().strip(),
    "org": escape(os.path.basename(f.name[:-2]), '"'),
    "ropts": ropts,
})

    # parse optctrl in netan["opt"]
    # optctrl_maxit=100 goes to list(default=list(maxit=100))
    # optctrl:bfgs:maxit=1000 goes to list(bfgs=list(maxit=1000))
    dctrl={"default": dict()}
    for k,v in netan["opt"].items():
        if not k.startswith("optctrl") or len(k) < 8:
            continue
        k=k[7:] # strip "optctrl" part
        if k[0] == "_":
            dctrl["default"][k[1:]]=str(v)
        elif k[0] == ":":
            li=k[1:].split(":")
            dctrl[li[0]]=dctrl.get(li[0], dict())
            dctrl[li[0]][join(":", li[1:])]=str(v)
    #print(dctrl)
    sep=",\n\t"
    tmp=f"list({sep.join('`'+m+'`=list('+', '.join('`'+kk+'`='+vv for kk,vv in dd.items())+')' for m,dd in dctrl.items())})"
    f.write(f"control_ftbl={tmp}")

    if case_i:
        f.write("""
source(file.path(dirr, "opt_icumo_tools.R"))
#loadcmp(file.path(dirr, "opt_icumo_tools.Rc"))

lab_resid=icumo_resid
lab_sim=param2fl_usm_rich
""")
    f.write("""
if (TIMEIT) {
   cat("rinit   : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}

# R profiling
if (prof) {
   Rprof(sprintf("%s.Rprof", baseshort))
}

nm_list=list()
nb_f=list()
""")
    netan2R_fl(netan, org, f)
    d=netan2R_rcumo(netan, org, f, emu)
    res.update(d)
    rc_keys=list(netan["rcumo_input"][0].keys())
    emu_keys=list(netan["emu_input"][0].keys()) if emu else []
    #import pdb; pdb.set_trace()
    f.write("""
nb_exp=%(nb_exp)d
nm_exp=c(%(nm_exp)s)
nm_list$nm_exp=nm_exp
# input cumomer vectors
xi=c(%(xi)s)
if (length(xi)) {
   dim(xi)=c(length(xi)/nb_exp, nb_exp)
} else {
   stop_mes("No reduced label entry is defined (may be because no measurement defined in FTBL). Cannot continue.", file=fcerr)
}
nm_xi=c(%(nm_xi)s)
rownames(xi)=nm_xi
nm_list$xi=nm_xi
nb_xi=length(nm_xi)
nb_f$xi=nb_xi
nb_cumoi=nb_xi
nm_inp=nm_xi
nm_incu=c("one", nm_xi, nm_rcumo)
nm_inlab=nm_incu
spa=spAbr
nm_x=nm_rcumo
nb_x=nb_rcumos
nb_f$rcumos=nb_rcumos
nb_f$cumoi=nb_cumoi
if (emu) {
   nm_emu=c(%(nm_emu)s)
   nb_emus=nb_rcumos*(seq_len(nb_rw)+1)
   nb_f$emus=nb_emus
   nm_list$emu=nm_emu
   nm_x=nm_emu
   nb_x=nb_emus
   xiemu=matrix(c(%(xiemu)s), ncol=nb_exp)
   nm_xiemu=c(%(nm_xiemu)s)
   nm_list$xiemu=nm_xiemu
   rownames(xiemu)=nm_xiemu
   nb_xiemu=length(nm_xiemu)
   nb_f$xiemu=nb_xiemu
   nb_f$xi=nb_xiemu
   nb_xi=nb_xiemu
   nm_inp=nm_xiemu
   xi=xiemu
   nm_inemu=c("one", nm_xiemu, nm_emu)
   nm_inlab=nm_inemu
   spa=spr2emu(spAbr, nm_incu, nm_inemu, nb_f)
}
# reorder indexes to accelerate sparse matrix construction
spa=sparse2spa(spa)
#browser()
# composite labeling vector incu c(1, xi, xc) names
nm_inlab=c("one", nm_inp, nm_x); # the constant 1 has name "one"
nm_list$x=nm_x
nm_list$inp=nm_inp
nb_f$x=nb_x
"""%{
    "nb_exp": len(netan["iso_input"]),
    "nm_exp": join(", ", netan["exp_names"], '"', '"'),
    "xi": join(", ", [li[k] if li[k]==li[k] else "NA" for li in netan["rcumo_input"] for k in rc_keys]),
    "nm_xi": join(", ", rc_keys, '"', '"'),
    "xiemu": join(", ", [li[k] if li[k]==li[k] else "NA" for li in netan["emu_input"] for k in emu_keys]),
    "nm_xiemu": join(", ", emu_keys, '"', '"'),
    "nm_emu": join(", ", valval(netan.get('vemu', [])), '"', '"'),
})
    if fullsys:
        d=netan2R_cumo(netan, org, f)
        res.update(d)
    else:
        f.write("nm_cumo=NULL\n")
    d=netan2R_meas(netan, org, f, emu)
    res.update(d)
    netan2R_ineq(netan, org, f)
    f.write("""
nb_sys=list(
   reactions=list(
      reversible=%(rrev)s,
      non_reversible=%(rnonrev)s
   ),
   fluxes=list(
      free=%(ff)s,
      dependent=%(fd)s,
      constrained=%(fc)s
   ),
   metabolites=list(
      input=%(minp)s,
      output=%(moutp)s,
      intra=%(mintra)s
   ),
   measurements=list(
      flux=%(meas_f)s,
      mass=%(meas_m)s,
      peak=%(meas_p)s,
      label=%(meas_l)s,
      metab=%(meas_pool)s
   ),
   equations=list(
      equalities=%(eqe)s,
      inequalities=%(eqi)s
   ),
   label_variables=list(
      full=c(%(lncumo)s),
      reduced_cumomers=c(%(lnrcumo)s)
   ),
   parallel_experiments=%(nb_exp)d
)
if (sum(nb_sys$label_variables$full)==0) {
   nb_sys$label_variables$full=NULL
}
if (emu) {
   x=nb_sys$label_variables$reduced_cumomers
   nb_sys$label_variables$reduced_cumomers=NULL
   nb_sys$label_variables$emu=paste(x, "*", seq_len(length(x)), "=", x*seq_len(length(x)))
}
"""%{
    "rrev": len(netan["reac"])-len(netan["notrev"]),
    "rnonrev": len(netan["notrev"]),
    "ff": len(netan["flux_free"]["net"])+len(netan["flux_free"]["xch"]),
    "fd": len(netan["vflux"]["net"])+len(netan["vflux"]["xch"]),
    "fc": len(netan["vflux_constr"]["net"])+len(netan["vflux_constr"]["xch"]),
    "minp": len(netan["input"]),
    "moutp": len(netan["output"]),
    "mintra": len(netan["metabs"])-len(netan["input"])-len(netan["output"]),
    "meas_f": len(netan["vflux_meas"]["net"]),
    "meas_m": sum(len(netan["measures"]["mass"][ili]["vec"]) for ili in range(nexp)),
    "meas_p": sum(len(netan["measures"]["peak"][ili]["vec"]) for ili in range(nexp)),
    "meas_l": sum(len(netan["measures"]["label"][ili]["vec"]) for ili in range(nexp)),
    "meas_pool": len(netan["metab_measured"]),
    "eqe": len(netan["flux_equal"]["net"])+len(netan["flux_equal"]["xch"]),
    "eqi": len(netan["flux_inequal"]["net"])+len(netan["flux_inequal"]["xch"]),
    "lncumo": ",".join(str(len(a)) for a in netan["cumo_sys"]["A"]),
    "lnrcumo": ",".join(str(len(a)) for a in netan["rcumo_sys"]["A"]),
    "nb_exp": len(netan["iso_input"])
    })
    return res

def netan2R_fl(netan, org, f):
    """netan2R_fl(netan, org, f)
    generate R code for flux and pool part
    for more details cf. netan2Rinit()
    """
    # dependent flux counts
    nb_fln=len(netan['vflux']['net'])
    nb_flx=len(netan['vflux']['xch'])
    nb_fl=nb_fln+nb_flx

    # prepare index translator for free fluxes
    # it will be used in bfl expressions where names like flx.net must
    # be mapped on respecting parameter index
    nb_ffn=len(netan['flux_free']['net'])
    nb_ffx=len(netan['flux_free']['xch'])
    nb_fcn=len(netan['flux_constr']['net'])
    nb_fcx=len(netan['flux_constr']['xch'])
    ffn2iprm=dict(("f.n."+f,(i+1))
        for (f,i) in netan['vflux_free']['net2i'].items())
    ffx2iprm=dict(("f.x."+f,(i+1+nb_ffn))
        for (f,i) in netan['vflux_free']['xch2i'].items())

    # prepare fwrv2i
    fwrv2i=dict((f,i+1) for (f,i) in netan["vflux_fwrv"]["fwrv2i"].items())
    nb_fwrv=len(netan["vflux_fwrv"]["fwrv2i"])

    # make tuple for complete flux vector d,f,c
    # (name,"d|f|c|g","n|x")
    tfallnx=list(zip(
            netan["vflux"]["net"]+
            netan["vflux_free"]["net"]+
            netan["vflux_constr"]["net"]+
            netan["vflux_growth"]["net"]+
            netan["vflux"]["xch"]+
            netan["vflux_free"]["xch"]+
            netan["vflux_constr"]["xch"]+
            netan["vflux_growth"]["net"],

            ["d"]*len(netan["vflux"]["net"])+
            ["f"]*len(netan["vflux_free"]["net"])+
            ["c"]*len(netan["vflux_constr"]["net"])+
            ["g"]*len(netan["vflux_growth"]["net"])+
            ["d"]*len(netan["vflux"]["xch"])+
            ["f"]*len(netan["vflux_free"]["xch"])+
            ["c"]*len(netan["vflux_constr"]["xch"])+
            ["g"]*len(netan["vflux_growth"]["net"]),

            ["n"]*len(netan["vflux"]["net"])+
            ["n"]*len(netan["vflux_free"]["net"])+
            ["n"]*len(netan["vflux_constr"]["net"])+
            ["n"]*len(netan["vflux_growth"]["net"])+
            ["x"]*len(netan["vflux"]["xch"])+
            ["x"]*len(netan["vflux_free"]["xch"])+
            ["x"]*len(netan["vflux_constr"]["xch"])+
            ["x"]*len(netan["vflux_growth"]["net"]),
            ))
    netan["f2dfcg_nx_f"]={
       "net": dict((fl, t+".n."+fl) for (fl,t,nx) in tfallnx if nx=="n"),
       "xch": dict((fl, t+".x."+fl) for (fl,t,nx) in tfallnx if nx=="x"),
    }

    f.write("""
if (TIMEIT) {
   cat("r_flux  : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
""")
    # auxiliary dict for edge-flux coupling
    f2edge=dict()
    for (fl,lr) in netan["sto_r_m"].items():
        if len(lr["left"])==1 and len(lr["right"])==1:
           f2edge[fl]=[lr["left"][0][0]+" ("+fl+") "+lr["right"][0][0]]
        else:
           f2edge[fl]=[]
           subs=[m for m,_ in lr["left"]] # substrates
           prods=[m for m,_ in lr["right"]] # products
           same_subs=len(subs)==2 and subs[0]==subs[1]
           same_prods=len(prods)==2 and prods[0]==prods[1]
           for (i, m) in enumerate(subs):
               f2edge[fl].append(m+" ("+fl+(str(i+1) if same_subs else "")+") "+fl)
           for (i, m) in enumerate(prods):
               f2edge[fl].append(fl+" ("+fl+(str(i+1) if same_prods else "")+") "+m)
    f.write("""
# fwd-rev flux names
nm_fwrv=c(%(nm_fwrv)s)

# edge to netflux name translator
edge2fl=c(%(edge2fl)s)
names(edge2fl)=c(%(nedge2fl)s)

# initialize the linear system Afl*flnx=bfl (0-weight cumomers)
# unknown net flux names
nm_fln=c(%(nm_fln)s)
nb_fln=length(nm_fln)
fln=c(%(fln)s)
names(fln)=nm_fln
# unknown xch flux names
nm_flx=c(%(nm_flx)s)
nb_flx=length(nm_flx)
flx=c(%(flx)s)
names(flx)=nm_flx
nm_fl=c(nm_fln, nm_flx)
nb_fl=nb_fln+nb_flx
fl=c(fln, flx)
# gather flux names in a list
nm_list$flnx=nm_fl
nm_list$fwrv=nm_fwrv

# carbon length of metabolites
clen=c(%(clen)s)
names(clen)=c(%(nm_metab)s)

# metabolite pools are : all (poolall) which is divided in free (poolf) and
# constrained (poolc)

# constrained pool
poolc=c(%(poolc)s)
nm_poolc=c(%(nm_poolc)s)
if (length(nm_poolc)) {
   names(nm_poolc)=substring(nm_poolc, 4)
}
names(poolc)=nm_poolc

# starting values for free pool (the same number and the same alphabetic order than free growth fluxes, if present)
poolf=c(%(poolf)s)
nm_poolf=c(%(nm_poolf)s)
if (length(nm_poolf)) {
   names(nm_poolf)=substring(nm_poolf, 4)
}
names(poolf)=nm_poolf
nb_poolf=length(poolf)
nb_f$nb_poolf=nb_poolf

nm_poolall=c(nm_poolf, nm_poolc)
poolall=as.numeric(c(poolf, poolc))
names(poolall)=nm_poolall
pool=poolall
nm_list$poolf=nm_poolf
nm_list$poolc=nm_poolc
nm_list$poolall=nm_poolall

# flux matrix
nb_flr=%(nb_flr)d
if (nb_fl) {
   Afl=matrix(0, nrow=nb_flr, ncol=nb_fl)
"""%{
    "nb_flr": len(netan["Afl"]),
    "nm_fwrv": join(", ", netan["vflux_fwrv"]["fwrv"], '"', '"'),
    "nm_fln": join(", ", netan["vflux"]["net"], '"d.n.', '"'),
    "fln": join(", ", (netan["flux_dep"]["net"][k] for k in netan["vflux"]["net"])),
    "nm_flx": join(", ", netan["vflux"]["xch"], '"d.x.', '"'),
    "flx": join(", ", (netan["flux_dep"]["xch"][k] for k in netan["vflux"]["xch"])),
    "edge2fl": join(", ", ('"'+netan["f2dfcg_nx_f"]["net"][fl]+'"' for (fl,l) in f2edge.items() for e in l)),
    "nedge2fl": join(", ", ('"'+e+'"' for (fl,l) in f2edge.items() for e in l)),
    "clen": join(",", list(netan["Clen"].values())),
    "nm_metab": join(",", list(netan["Clen"].keys()), '"', '"'),
    "poolf": join(", ", (-netan["met_pools"][m] for m in netan["vpool"]["free"])),
    "nm_poolf": join(", ", netan["vpool"]["free"], '"pf:', '"'),
    "poolc": join(", ", (netan["met_pools"][m] for m in netan["vpool"]["constrained"])),
    "nm_poolc": join(", ", netan["vpool"]["constrained"], '"pc:', '"'),
})
    for (i,row) in enumerate(netan["Afl"]):
        f.write(
"""   Afl[%(i)d, c(%(ic)s)]=c(%(v)s)
"""%{
    "i": i+1,
    "ic": join(", ", (i+1 for (i,v) in enumerate(row) if v!=0.)),
    "v": join(", ", (v for v in row if v!=0.)),
})
    f.write(
"""} else {
   Afl=matrix(0., nb_fl, nb_fl)
}
dimnames(Afl)=list(c(%(nm_rows)s), nm_fl)
#browser()
# prepare param (\Theta) vector
# order: free flux net, free flux xch, scale label, scale mass, scale peak
param=numeric(0)
nm_par=c()
# free net fluxes
nb_ffn=%(nb_ffn)d
nm_ffn=c(%(nm_ffn)s)
# starting values for iterations
param=c(param, c(%(ffn)s))
if (nb_ffn) {
   nm_par=c(nm_par, nm_ffn)
}
# free xch fluxes
nb_ffx=%(nb_ffx)d
nm_ffx=c(%(nm_ffx)s)
# starting values for iterations
param=c(param, c(%(ffx)s))
if (nb_ffx) {
   nm_par=c(nm_par, nm_ffx)
}
names(param)=nm_par
ff=param
nm_ff=c(nm_ffn, nm_ffx)
nm_list$ff=nm_ff
nb_param=length(param)
# scaling factors are added to param later

nb_ff=nb_ffn+nb_ffx

# constrained fluxes
# net
nb_fcn=%(nb_fcn)d
nm_fcn=c(%(nm_fcn)s)
fcn=c(%(fcn)s)
# xch
nb_fcx=%(nb_fcx)d
nm_fcx=c(%(nm_fcx)s)
fcx=c(%(fcx)s)
fc=c(fcn, fcx)
nm_fc=c(nm_fcn, nm_fcx)
names(fc)=nm_fc
nb_fc=nb_fcn+nb_fcx

# variable growth fluxes (constant are already accounted in constrained fluxes)
nb_fgr=%(nb_fgr)d
nm_fgr=c(%(nm_fgr)s)
fgr=c(%(fgr)s)
nm_list$fgr=nm_fgr
nb_f$nb_fgr=nb_fgr

# total flux vector fallnx dimension
nb_fallnx=nb_fl+nb_ff+nb_fc+nb_fgr+nb_fgr
nb_fwrv=nb_fallnx

# net dependent and free fluxes
nm_dfn=c(nm_fln, nm_ffn)
names(nm_dfn)=substring(nm_dfn, 5)

# all flux cardinals
nb_f=append(nb_f, list(nb_fln=nb_fln, nb_flx=nb_flx, nb_fl=nb_fl,
   nb_ffn=nb_ffn, nb_ffx=nb_ffx, nb_ff=nb_ff,
   nb_fcn=nb_fcn, nb_fcx=nb_fcx, nb_fc=nb_fc,
   nb_fallnx=nb_fallnx, nb_fwrv=nb_fwrv,
   nb_fgr=nb_fgr,
   include_growth_flux=%(inc_gr_f)s,
   mu=%(mu)s))
"""%{
    "nm_rows": join(", ", netan["vrowAfl"], '"', '"'),
    "nb_ffn": nb_ffn,
    "nb_ffx": nb_ffx,
    "nm_ffn": join(", ", netan["vflux_free"]["net"], '"f.n.', '"'),
    "nm_ffx": join(", ", netan["vflux_free"]["xch"], '"f.x.', '"'),
    "ffn": join(", ", [netan["flux_free"]["net"][fl]
        for fl in netan["vflux_free"]["net"]]),
    "ffx": join(", ", [netan["flux_free"]["xch"][fl]
        for fl in netan["vflux_free"]["xch"]]),
    "nb_fcn": len(netan["flux_constr"]["net"]),
    "nb_fcx": len(netan["flux_constr"]["xch"]),
    "nm_fcn": join(", ", netan["vflux_constr"]["net"], '"c.n.', '"'),
    "nm_fcx": join(", ", netan["vflux_constr"]["xch"], '"c.x.', '"'),
    "fcn": join(", ", [netan["flux_constr"]["net"][fl]
        for fl in netan["vflux_constr"]["net"]]),
    "fcx": join(", ", [netan["flux_constr"]["xch"][fl]
        for fl in netan["vflux_constr"]["xch"]]),
    "inc_gr_f": "TRUE" if netan["opt"].get("include_growth_flux") else "FALSE",
    "mu": str(netan["opt"].get("mu", "NULL")),
    "nb_fgr": len(netan["vflux_growth"]["net"]),
    "nm_fgr": join(", ", netan["vflux_growth"]["net"], '"g.n.', '"'),
    "fgr": join(", ", [netan["flux_vgrowth"]["net"][fl]
        for fl in netan["vflux_growth"]["net"]]),
})
    f.write("""
# prepare p2bfl, c2bfl, g2bfl, cnst2bfl matrices such that p2bfl%*%param[1:nb_ff]+
# c2bfl%*%fc+g2bfl%*%fgr+cnst2bfl=bfl
# replace f.[nx].flx by corresponding param coefficient
p2bfl=simple_triplet_zero_matrix(nrow=nb_flr, ncol=nb_ff)
# replace c.[nx].flx by corresponding fc coefficient
c2bfl=simple_triplet_zero_matrix(nrow=nb_flr, ncol=nb_fc)
# variable growth fluxes
g2bfl=simple_triplet_zero_matrix(nrow=nb_flr, ncol=nb_fgr)
cnst2bfl=numeric(nb_flr); # may be coming from equalities
colnames(p2bfl)=nm_par
colnames(c2bfl)=nm_fc
colnames(g2bfl)=nm_fgr
""")
    row={"f": None, "c": None, "g": None, "cnst": None}
    for (i,item) in enumerate(netan["bfl"]):
        if not item:
            continue
        # split terms in flux types
        row["cnst"]=item.get("")
        row["f"]=dict((k,v) for (k,v) in item.items() if k[0:2]=="f.")
        row["c"]=dict((k,v) for (k,v) in item.items() if k[0:2]=="c.")
        row["g"]=dict((k,v) for (k,v) in item.items() if k[0:2]=="g.")
        f.write("\n")
        if row["f"]:
            f.write("p2bfl[%(i)d, pmatch(c(%(if)s), nm_par)]=c(%(rowf)s);\n"%\
                {"i": i+1,
                "if": join(", ", list(row["f"].keys()), p='"', s='"'),
                "rowf": join(", ", list(row["f"].values())),
                })
        if row["c"]:
            f.write("c2bfl[%(i)d, pmatch(c(%(ic)s), nm_fc)]=c(%(rowc)s);\n"%\
                {"i": i+1,
                "ic": join(", ", list(row["c"].keys()), p='"', s='"'),
                "rowc": join(", ", list(row["c"].values())),
                })
        if row["g"]:
            f.write("g2bfl[%(i)d, pmatch(c(%(ig)s), nm_fgr)]=c(%(rowg)s);\n"%\
                {"i": i+1,
                "ig": join(", ", list(row["g"].keys()), p='"', s='"'),
                "rowg": join(", ", list(row["g"].values())),
                })
        if row["cnst"]:
            f.write("cnst2bfl[%(i)d]=%(rowcnst)s;\n"%{"i": i+1, "rowcnst": row["cnst"],})
    f.write("""
bp=as.numeric(c2bfl%stm%fc+cnst2bfl)
""")

    f.write("""
if (ffguess) {
   # make an automatic guess for free/dependent flux partition
   afd=as.matrix(cBind(Afl, -p2bfl))
   qafd=qr(afd, LAPACK=TRUE)
   d=abs(diag(qafd$qr))
   rank=sum(d > d[1]*1.e-10)
   qrow=qr(t(afd))
   rankr=qrow$rank
   if (rank != rankr)
      stop_mes("Weird error: column and row ranks are not equal.", file=fcerr)
   
   irows=qrow$pivot[seq_len(rankr)]
   if (rank==0) {
      stop_mes("Error: No free/dependent flux partition could be made. Stoichiometric matrix has rank=0.", file=fcerr)
   }
   Afl=afd[irows, qafd$pivot[1L:rank], drop=FALSE]
   ka=kappa(Afl)
   if (ka > 1.e7) {
      mes=sprintf("Error: No working free/dependent flux partition could be proposed. Stoichiometric matrix has condition number %g.\\n", ka)
      stop_mes(mes, file=fcerr)
   }
   p2bfl=-as.simple_triplet_matrix(afd[irows, qafd$pivot[-seq_len(rank)], drop=FALSE])
   c2bfl=c2bfl[irows, , drop=FALSE]
   g2bfl=g2bfl[irows, , drop=FALSE]
   cnst2bfl=cnst2bfl[irows]
   bp=bp[irows]
   
   # replace names
   nm_fl=sub("f.", "d.", colnames(Afl), fixed=TRUE)
   colnames(Afl)=nm_fl # both net and xch
   nm_fln=sort(grep("^d.n.", nm_fl, v=TRUE))
   nm_flx=sort(grep("^d.x.", nm_fl, v=TRUE))
   nm_fl=c(nm_fln, nm_flx)
   Afl=Afl[, nm_fl, drop=FALSE]
   
   nm_ff=sub("d.", "f.", colnames(p2bfl), fixed=TRUE) # both net and xch
   colnames(p2bfl)=nm_ff
   nm_ffn=sort(grep("^f.n.", nm_ff, v=TRUE))
   nm_ffx=sort(grep("^f.x.", nm_ff, v=TRUE))
   nm_ff=c(nm_ffn, nm_ffx)
   p2bfl=p2bfl[, nm_ff, drop=FALSE]
   
   # remake param vector
   if (!fdfit)
      param=c(runif(length(nm_ff)), if (nb_ff == 0) param else param[-seq_len(nb_ff)])
   names(param)[seq(along=nm_ff)]=nm_ff
#browser()
}
nm_list$flnx=nm_fl
nm_fallnx=c(nm_fln, nm_ffn, nm_fcn, nm_fgr, nm_flx, nm_ffx, nm_fcx, sub(".n.", ".x.", nm_fgr, fixed=TRUE))
nm_list$fallnx=nm_fallnx
nm_net=c(nm_fln, nm_ffn, nm_fcn)
names(nm_net)=substring(nm_net, 5)
nm_xch=c(nm_flx, nm_ffx, nm_fcx)
names(nm_xch)=substring(nm_xch, 5)
edge2fl[]=nm_net[substring(edge2fl, 5)]
nm_list$ff=nm_ff

# accounting numbers
nb_flr=nrow(Afl)
nb_param=length(param)
nb_ffn=length(nm_ffn)
nb_ffx=length(nm_ffx)
nb_ff=nb_ffn+nb_ffx
nb_fln=length(nm_fln)
nb_flx=length(nm_flx)
nb_fl=nb_fln+nb_flx
nm_par=names(param)

for (item in c("nb_fln", "nb_flx", "nb_fl", "nb_ffn", "nb_ffx", "nb_ff")) {
   nb_f[item]=get(item)
}
# translation from n-x to fw-rv
sh_fwrv=substring(nm_fwrv[1:(nb_fwrv/2)], 5)
sh_nx=substring(nm_fallnx, 2)
nb_f$inet2ifwrv=pmatch(paste(".n.", sh_fwrv, sep=""), sh_nx)
nb_f$ixch2ifwrv=pmatch(paste(".x.", sh_fwrv, sep=""), sh_nx)
#nb_f$inet2ifwrv=sapply(nm_fwrv[1:(nb_fwrv/2)], function(f) grep(sprintf("^.\\\\.n\\\\.%s$", substring(f, 5)), nm_fallnx))
#nb_f$ixch2ifwrv=sapply(nm_fwrv[1:(nb_fwrv/2)], function(f) grep(sprintf("^.\\\\.x\\\\.%s$", substring(f, 5)), nm_fallnx))

if (TIMEIT) {
   cat("Afl qr(): ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}

qrAfl=qr(Afl, LAPACK=TRUE)
d=abs(diag(qrAfl$qr))
qrAfl$rank=sum(d > d[1]*1.e-10)
rank=qrAfl$rank
aful=as.matrix(cBind(Afl, -p2bfl, -c2bfl))
qrow=qr(t(aful))
rankr=qrow$rank
#browser()
# first check the presence of lindep rows
if (nrow(Afl) > rankr) {
   prop=sprintf("Error: Among %d equations (rows), %d are redundant and must be eliminated by hand.\\n", nrow(Afl), nrow(Afl)-rankr)
   prop=paste(prop, "Candidate(s) for elimination is (are):\\n",
      paste(rownames(Afl)[qrow$pivot[-(1:rankr)]], sep="", collapse="\\n"),
               "\\n", sep="")
   stop_mes(prop, file=fcerr)
}
if (nrow(Afl) != rank || nrow(Afl) != ncol(Afl)) {
   #write.table(Afl)
   mes=NULL
   if (nrow(Afl) <= rank) {
      mes=paste("Candidate(s) for free or constrained flux(es):\\n",
         paste(colnames(Afl)[-qrAfl$pivot[1L:nrow(Afl)]], collapse="\\n"),
         "\\nFor this choice, condition number of stoichiometric matrix will be ",
         kappa(Afl[,qrAfl$pivot[1L:nrow(Afl)],drop=FALSE]), "\\n", sep="")
   } else if (nrow(Afl) > rank) {
      nextra=nrow(Afl)-rank
      comb=combn(c(nm_ffn, colnames(Afl)[-qrAfl$pivot[1L:rank]]), nextra)
      aextra=cBind(Afl[,-qrAfl$pivot[1L:rank],drop=FALSE], -p2bfl)
      colnames(aextra)=c(colnames(Afl)[-qrAfl$pivot[1L:rank]], colnames(p2bfl))
      ara=Afl[,qrAfl$pivot[1L:rank],drop=FALSE]
      i=which.min(apply(comb, 2, function(i) kappa(cBind(ara, aextra[,i]))))[1L]
      nm_tmp=comb[,i]
      ka=kappa(cBind(ara, aextra[,nm_tmp]))
      if (ka < 1.e7) {
         prop=paste("Proposal to declare dependent flux(es) is:\\n",
            paste(nm_tmp, collapse="\\n"), "\\n", sep="")
         if (rank < ncol(Afl)) {
            prop=prop%s+%"While the following dependent flux(es) should be declared free or constrained:\\n"%s+%join("\\n", colnames(Afl)[-qrAfl$pivot[1L:rank]])%s+%"\\n"
         }
         prop=paste(prop, "For this choice, condition number of stoichiometric matrix will be ", ka, "\\n", sep="")
      } else {
         # add constraint fluxes to candidate list
         if (nb_fcn > 0) {
            aextra=as.matrix(cBind(Afl[,-qrAfl$pivot[1L:rank],drop=FALSE], -p2bfl, -c2bfl))
            colnames(aextra)=c(colnames(Afl)[-qrAfl$pivot[1L:rank]], colnames(p2bfl), colnames(c2bfl))
         }
         aextended=aful
         qae=qr(aextended, LAPACK=TRUE)
         d=abs(diag(qae$qr))
         ranke=sum(d > d[1L]*1.e-10)
         if (ranke == nrow(Afl)) {
            prop=paste("Proposal to declare dependent flux(es) is:\\n",
            join("\\n", colnames(aextended)[qae$pivot[1L:ranke]]), "\\n",
            "while free and constrained fluxes should be:\\n",
            join("\\n", colnames(aextended)[-qae$pivot[1L:ranke]]), "\\n",
            sep="")
            ka=kappa(aextended[,qae$pivot[1L:ranke]])
            prop=paste(prop, "For this choice, condition number of stoichiometric matrix will be ", ka, "\\n", sep="")
         } else {
            prop="No proposal for partition dependent/free fluxes could be made.\\n"
         }
      }
      mes=paste("There is (are) probably ", nextra,
         " extra free flux(es) among the following:\\n",
         paste(nm_ffn, collapse="\\n"), "\\n",
         prop,
         sep="")
   }""")
    f.write("""
   stop_mes("Flux matrix is not square or is singular: (", nrow(Afl), "eq x ", ncol(Afl), "unk)\\n",
      "You have to change your choice of free fluxes in the '%(n_ftbl)s' file.\\n",
      mes, file=fcerr)
}

# make sure that free params choice leads to not singular matrix
if (qrAfl$rank != nb_fl) {
   #write.table(Afl)
   # make a suggestion of new free fluxes
   A=cBind(Afl, -p2bfl, -c2bfl)
   colnames(A)=c(colnames(Afl), nm_ff, nm_fc)
   qa=qr(A, LAPACK=TRUE)
   d=diag(qa$qr)
   qa$rank=sum(abs(d)>=abs(d[1]*1.e-10))
   
   mes=paste("Error: Dependent flux matrix is singular.\\n",
      "Change your partition on free/dependent/constrained fluxes in the '%(n_ftbl)s' file.\\n",
      "Can not resolve dependent fluxe(s):\\n",
      paste(colnames(Afl)[-qrAfl$pivot[(1:qrAfl$rank)]], collapse="\\n"),
      sep="")
   if (qa$rank==nb_fl) {
      mes=paste(mes,
      "\\n\\nSuggested dependent fluxes:\\n",
      paste(colnames(A)[qa$pivot[(1:qa$rank)]], collapse="\\n"),
      "\\n\\nWhich would give the following free and constrained fluxes:\\n",
      paste(colnames(A)[-qa$pivot[(1:qa$rank)]], collapse="\\n"), "\\n",
      sep="")
   } else {
      mes=paste(mes, "\\nNo suggested free fluxes could be found", sep="")
   }
   stop_mes(mes, file=fcerr)
}

# inverse flux matrix
invAfl=solve(qrAfl)
""" % {
    "n_ftbl": escape(org+".ftbl", "\\"),
    })
    f.write("""
if (fdfit) {
#browser()
   # choose free fluxe values such that they fit starting values from ftbl
   #dep=invAfl%*%(p2bfl%*%ff+bp)
   #ff=ff
   ff=qr.solve(rbind(invAfl%stm%p2bfl, diag(ncol(p2bfl))), c(fl-invAfl%*%bp, param))
}
# intermediate jacobian
if (TIMEIT) {
   cat("dfl_dffg: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}

dfl_dffg=invAfl%stm%p2bfl
if (nb_fgr > 0L) {
   dfl_dffg=cBind(dfl_dffg, invAfl%stm%g2bfl)
}
dimnames(dfl_dffg)=list(nm_fl, c(nm_ff, nm_fgr))
dfl_dffg[abs(dfl_dffg) < 1.e-14]=0.
nb_f$dfl_dffg=as.simple_triplet_matrix(dfl_dffg)

# prepare mf, md, mc and mg matrices
# such that mf%*%ff+md%*%fl+mc%*%fc+mg%*%fgr gives fallnx
# here ff free fluxes (param), fl are dependent fluxes, fc are constrained
# fluxes and fgr are variable growth fluxes
mf=matrix(0., nb_fallnx, nb_ff)
dimnames(mf)=list(nm_fallnx, nm_ff)
md=matrix(0., nb_fallnx, nb_fl)
dimnames(md)=list(nm_fallnx, nm_fl)
mc=matrix(0., nb_fallnx, nb_fc)
dimnames(mc)=list(nm_fallnx, nm_fc)
mg=matrix(0., nb_fallnx, nb_fgr)
dimnames(mg)=list(nm_fallnx, nm_fgr)

if (nb_ff > 0) {
   mf[nm_ff, nm_ff]=diag(1., nb_ff)
}
if (nb_fl > 0) {
   md[nm_fl, nm_fl]=diag(1., nb_fl)
}
if (nb_fc > 0) {
   mc[nm_fc, nm_fc]=diag(1., nb_fc)
}
if (nb_fgr > 0) {
   mg[nm_fgr, nm_fgr]=diag(1., nb_fgr)
}
""")
    netan["fwrv2i"]=fwrv2i
    netan["tfallnx"]=tfallnx

def netan2R_meas(netan, org, f, emu=False):
    """netan2R_meas(netan, org, f)
    generate code for measure treatment
    """
    # prepare python measures
    if "measures" not in netan:
        #print("Calculate measures in netan2R_meas.")
        measures=dict()
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas)
        netan["measures"]=measures
    measures=netan["measures"]
    nexp=len(netan["iso_input"])
    #aff("got measures in netan2R_meas", measures);##
    # get scaling factors and their indexes, measure matrices, and measured cumomer value vector
    scale=[{"label": {}, "mass": {}, "peak": {}} for i in range(nexp)] # for unique scale names
    nrow=[{"label": {}, "mass": {}, "peak": {}} for i in range(nexp)] # for counting scale names
    o_sc=[{"label": {}, "mass": {}, "peak": {}} for i in range(nexp)] # for ordered unique scale names
    o_meas=list(measures.keys()); # ordered measure types
    o_meas.sort()

    ir2isc=[{"label": [], "mass": [], "peak": []} for i in range(nexp)] # for mapping measure rows indexes on scale index
    # we want to use it in python like isc[ili][meas]=ir2isc[ili][meas][ir]
    for meas in o_meas:
        for ili in range(nexp):
            # get unique scaling factors
            # and count rows in each group
            # row["scale"] is "metab;group" (metab name may be fake here)
            for (i,row) in enumerate(measures[meas][ili]["mat"]):
                scale[ili][meas][row["scale"]]=0.
                nrow[ili][meas][row["scale"]]=nrow[ili][meas].get(row["scale"],0.)+1
            # remove groups having only one measure in them
            for (k,n) in list(nrow[ili][meas].items()):
                if n<2:
                    del(scale[ili][meas][k])
            # order scaling factor
            o_sc[ili][meas]=list(scale[ili][meas].keys())
            o_sc[ili][meas].sort()
            # map a measure rows (card:n) on corresponding scaling factor (card:1)
            # if a row has not scale factor it is scaled with factor 1
            # vector having scaling parameters is formed like
            # c(1,param)
            ir2isc[ili][meas]=[-1]*len(measures[meas][ili]["mat"])
            for (i,row) in enumerate(measures[meas][ili]["mat"]):
                if row["scale"] in scale[ili][meas]:
                    ir2isc[ili][meas][i]=o_sc[ili][meas].index(row["scale"])
    
            # measured value vector is in measures[meas]["vec"]
            # measured dev vector is in measures[meas]["dev"]

    # create R equivalent structures with indices for scaling
    f.write("""
if (TIMEIT) {
   cat("measure : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
if (!noscale) {
   # make place for scaling factors
   nb_sc=vector("integer", %d)
"""%nexp)
    for ili in range(nexp):
        f.write("# experiment: %d\n"%(ili+1))
        for meas in o_meas:
            if not o_sc[ili][meas]:
                continue
            f.write("""
   # %(meas)s
   # initial values for scales are set later
   param=c(param,%(sc)s)
   nm_par=c(nm_par,c(%(sc_names)s))
   names(param)=nm_par
""" % {
        "meas": meas,
        "sc": join(", ", ["1"]*len(o_sc[ili][meas])),
        "sc_names": join(", ", o_sc[ili][meas], '"'+str(ili+1)+":"+meas+';', '"'),
        })

        f.write("""
   nb_param=length(param)
   nb_sc[%d]=nb_param-nb_ff # at this moment it is cumulated sum. diff() is taken later
"""%(ili+1))
    f.write("""
   # indices mapping from scaling to measure matrix row
   # c(1,par)[ir2isc[[iexp]]] replicates scale parameters
   # for corresponding rows of measure matrix
   ir2isc=vector("list", %d)
"""%nexp)
    #base_isc=2+len(netan["flux_free"]["net"])+len(netan["flux_free"]["xch"])
    base_isc=2 # the shift by nb_ff is made in R because ffguess can make vary nb_ff in runtime
    #import pdb; pdb.set_trace()
    for ili in range(nexp):
        for meas in o_meas:
            if not ir2isc[ili][meas]:
                continue
            f.write("""
   # %(iexp)d:%(meas)s
   ir2isc[[%(iexp)d]]=c(ir2isc[[%(iexp)d]],c(%(ir2isc)s))
""" % {
        "iexp": ili+1,
        "meas": meas,
        "ir2isc": join(", ", ((str(ir2isc[ili][meas][ir]+base_isc) if ir2isc[ili][meas][ir]>=0 else 1) for ir in range(len(ir2isc[ili][meas]))))
        })
            base_isc=base_isc+len(scale[ili][meas])
        f.write("""
   isc=ir2isc[[%(iexp)d]] != 1
   ir2isc[[%(iexp)d]][isc]=ir2isc[[%(iexp)d]][isc]+nb_ff
"""%{"iexp": ili+1})

    f.write("""
   # cumulated base for nb_sc
   nb_sc_base=c(0, nb_sc[-nb_exp])
   nb_sc=diff(c(0, nb_sc))
   nb_sc_tot=sum(nb_sc)
   nb_f$nb_sc=nb_sc
   nb_f$nb_sc_tot=nb_sc_tot
   nb_f$nb_sc_base=nb_sc_base
#browser()
} else {
   # no scaling
   ir2isc=list()
   nb_sc=integer(nb_exp)
   nb_sc_tot=0
}
nb_f$nb_sc=nb_sc
nb_f$nb_sc_tot=nb_sc_tot
nm_list$par=nm_par
""")
    # get the full dict of non zero cumomers involved in measures
    # cumo=metab:icumo where icumo is in [1;2^Clen]
    # or emu=metab:ifrag+Mi
    f.write("""
# make a list of sparse measurement matrices
# measmat*xr+memaone gives a vector of simulated not-yet-pooled and not-yet-scaled measurements
# all but 0. Coefficients of 0-cumomers (by defenition equal to 1)
# are all regrouped in the memaone.
nm_measmat=nm_meas=nb_meas=nb_measmat=measmat=memaone=measvec=measdev=ipooled=vector("list", %d)
"""%nexp)
    for ili in range(nexp):
        meas_cumos=dict()
        for meas in o_meas:
            for row in measures[meas][ili]["mat"]:
                metab=row["metab"]
                if emu:
                    meas_cumos.update((metab+":"+i, "") for i in list(row["emuco"].keys()))
                else:
                    meas_cumos.update((metab+":"+str(icumo), "") for icumo in list(row["coefs"].keys()) if icumo != 0)

        # order involved cumomers (emu)
        o_mcumos=list(meas_cumos.keys())
        o_mcumos.sort()
        imcumo2i=dict((cumo, i) for (i, cumo) in enumerate(o_mcumos))
        nb_mcumo=len(o_mcumos)
        f.write("""
nm_measmat[[%(ili)d]]=c(%(idmeasmat)s)
nm_meas[[%(ili)d]]=c(%(idmeas)s)
nb_meas[[%(ili)d]]=length(nm_meas[[%(ili)d]])
nb_measmat[[%(ili)d]]=length(nm_measmat[[%(ili)d]])
measmat[[%(ili)d]]=simple_triplet_zero_matrix(nrow=nb_measmat[[%(ili)d]], ncol=%(ncol)d)
dimnames(measmat[[%(ili)d]])=list(nm_measmat[[%(ili)d]], nm_x)
memaone[[%(ili)d]]=numeric(nb_measmat[[%(ili)d]])
measvec[[%(ili)d]]=c(%(vmeas)s)
measdev[[%(ili)d]]=c(%(dev)s)
names(measvec[[%(ili)d]])=nm_meas[[%(ili)d]]
names(measdev[[%(ili)d]])=nm_meas[[%(ili)d]]
ipooled[[%(ili)d]]=list(ishort=pmatch(nm_meas[[%(ili)d]], nm_measmat[[%(ili)d]]))
"""%{
    "ili": ili+1,
    "nrow": len([measures[meas][ili]["vec"] for meas in measures]),
    "ncol": sum(len(l) for l in (netan["vemu"] if emu else netan["vrcumo"])),
    "idmeasmat": join(", ", (row["id"] for row in
        valval(measures[o][ili]["mat"] for o in o_meas)),
        p='"', s='"'),
    "idmeas": join(", ",  valval([v for o in o_meas for v in measures[o][ili]["ids"]]), p='"', s='"'),
    "vmeas": join(", ", valval(measures[o][ili]["vec"] for o in o_meas)).replace("nan", "NA"),
    "dev": join(", ", (sd for sd in valval(measures[o][ili]["dev"]
        for o in o_meas))),
    })
    f.write("""
nm_meas_tot=unlist(nm_meas)
nb_meas=unlist(nb_meas)
nb_meas_cumo=c(0., cumsum(nb_meas[-nb_exp]))
iexp_meas=lapply(seq_len(nb_exp), function(iexp) seq_len(nb_meas[iexp])+nb_meas_cumo[iexp])
nm_list$meas=nm_meas
nm_list$measmat=nm_measmat
nm_list$meas_tot=nm_meas_tot
nb_f$nb_meas=nb_meas
""")

    # get coeffs in the order above with their corresponding indices from total cumomer vector
    for ili in range(nexp):
        base_pooled=0
        for meas in o_meas:
            if not measures[meas][ili]["mat"]:
                continue
            #print("meas="+meas+"; mat="+str(measures[meas]["mat"]));##
            for metpool in measures[meas][ili]["pooled"]:
                f.write("""
# prepare indices of pooled measurements
ipooled[[%(ili)d]][["%(rowid)s"]]=1+%(basep)d+c(%(ind)s)
"""%{
    "ili": ili+1,
    "rowid": metpool[0],
    "ind": join(", ", metpool[1:]),
    "basep": base_pooled,
    }
                )
            base_pooled=base_pooled+len(measures[meas][ili]["vec"])+sum(len(item)-2 for item in measures[meas][ili]["pooled"])
    # preepare measmat indexes and values : ir, ic, val
    f.write("""
names(measvec)=names(measdev)=nm_exp

if (!noscale) {
   for (iexp in seq_len(nb_exp))
      ir2isc[[iexp]]=ir2isc[[iexp]][ipooled[[iexp]]$ishort]

   # prepare indexes of dispatching scale params in jacobian
   if (nb_sc_tot > 0) {
      nb_f$is2m=vector("list", nb_exp)
      for (iexp in seq_len(nb_exp)) {
         ipaire=matrix(0, nrow=0, ncol=2)
         tmp=lapply(seq_len(nb_sc[[iexp]]), function(isc) {
            i=which(ir2isc[[iexp]]==isc+nb_sc_base[iexp]+1+nb_ff)
            ipaire <<- rbind(ipaire, cbind(i, isc+nb_sc_base[iexp]))
            return(NULL)
         })
         nb_f$is2m[[iexp]]=ipaire
         # place holder for scale part of jacobian
         jx_f$dr_dsc[[iexp]]=simple_triplet_zero_matrix(nrow=length(ir2isc[[iexp]]), ncol=nb_sc_tot)
      }
   }
}
# prepare measmat indexes and values : ir, ic, val
""")
    lab2i0=netan["emu2i0" if emu else "rcumo2i0"]
    onelab="0+0" if emu else 0
    fcoef="emuco" if emu else "coefs"
    for ili in range(nexp):
        f.write("""
ind_mema=matrix(c(
""")
        i=0
        for meas in o_meas:
            if not measures[meas][ili]["mat"]:
                continue
            for row in measures[meas][ili]["mat"]:
                i+=1
                metab=row["metab"]
                f.write("""%(iricval)s,
"""%{
    "iricval": join(", ", valval((i, lab2i0[metab+":"+str(k)]+1, v)
        for (k, v) in row[fcoef].items() if k != onelab))
})

        f.write(r"""
NULL), ncol=3, byrow=TRUE); # close ind_mema creation
measmat[[%(iexp)d]][ind_mema[,1:2,drop=FALSE]]=ind_mema[,3]
memaone[[%(iexp)d]]=c(%(memaone)s)
"""%{
    "iexp": ili+1,
    "memaone": join(", ", (row[fcoef].get(onelab, 0.)
        for meas in o_meas
        for row in measures[meas][ili]["mat"])),
})
    f.write(r"""
pwe=ipwe=ip2ipwe=pool_factor=ijpwef=dp_ones=meas2sum=dpw_dpf=ipf_in_ppw=vector("list", nb_exp)
mets_in_res=vector("list", nb_exp)
for (iexp in seq_len(nb_exp)) {
   names(memaone[[iexp]])=nm_measmat[[iexp]]

   # prepare weights of label data for pooled metabs
   # prepare ipwe and ip2ipwe such that pwe[ipwe]=pool[ip2ipwe]
   # gives a good base for weight sum and normalization
   pwe[[iexp]]=double(nb_measmat[[iexp]])+1.
   mets_in_res[[iexp]]=sapply(nm_measmat[[iexp]], function(m) strsplit(m, ":")[[1L]][2L])
   for (po in names(ipooled[[iexp]])) {
      if (po == "ishort") next
      nm_sum=strsplit(po, ":")[[1L]][2L]
      mets=strsplit(nm_sum, "\\+")[[1L]]
      irpo=ipooled[[iexp]][[po]]
      ipwe[[iexp]]=c(ipwe[[iexp]], irpo) # where weighting is
      i=pmatch(mets, names(nm_poolf))
      if (any(!is.na(i))) {
         for (ir in i) {
            if (is.na(ir)) next
            ijpwef[[iexp]]=rbind(ijpwef[[iexp]], cbind(irpo, ir)) # where free pools matter
         }
      }
      ip2ipwe[[iexp]]=c(ip2ipwe[[iexp]], pmatch(mets, names(nm_poolall)))
      mets_in_res[[iexp]][irpo]=mets
   }
   # order ijpwef for sparse matrix ordering
   if (!is.null(ijpwef[[iexp]])) {
      o=order(ijpwef[[iexp]][,2L], ijpwef[[iexp]][,1L])
      ijpwef[[iexp]]=ijpwef[[iexp]][o,,drop=FALSE]
   }
   pool_factor[[iexp]]=as.factor(nm_measmat[[iexp]])
   # free pool in principal pool weight
   ipf_in_ppw[[iexp]]=apply(outer(mets_in_res[[iexp]], names(nm_poolf), "=="), 1, function(v) if(length(w <-which(v))) w else NA)
   ipf_in_ppw[[iexp]][is.na(ipf_in_ppw[[iexp]])]=0L
   dp_ones[[iexp]]=matrix(0., nb_measmat[[iexp]], nb_poolf)
   dp_ones[[iexp]][cbind(ipwe[[iexp]], ipf_in_ppw[[iexp]][ipwe[[iexp]]])]=1.
   
   # matrix for summing weighted measurements
   meas2sum[[iexp]]=simple_triplet_zero_matrix(length(ipooled[[iexp]]$ishort), nb_measmat[[iexp]])
   meas2sum[[iexp]][cbind(pmatch(nm_measmat[[iexp]], nm_measmat[[iexp]][ipooled[[iexp]]$ishort], dup=TRUE),       seq_len(nb_measmat[[iexp]]))]=1.
   dimnames(meas2sum[[iexp]])=list(nm_meas[[iexp]], nm_measmat[[iexp]])
   
   # dpw_dpf - matrix for derivation of pool weights by free pools
   if (nb_poolf > 0L && length(ijpwef[[iexp]]) > 0) {
      # indeed, we'll have to do weight derivation by free pools
      dpw_dpf[[iexp]]=simple_triplet_zero_matrix(nb_measmat[[iexp]], nb_poolf)
      dpw_dpf[[iexp]][ijpwef[[iexp]]]=1.
   }
}

""")
    f.write("""
# prepare flux measurements
nm_fmn=nm_net[c(%(nm_fmn)s)]
nm_list$fmn=nm_fmn
nb_fmn=length(nm_fmn)
nb_f$nb_fmn=nb_fmn

# measured values
fmn=c(%(fmn)s)

# SD for flux measurements
fmndev=c(%(fmndev)s)
if (nb_fmn)
   names(fmndev)=names(fmn)=nm_fmn

# indices for measured fluxes
# fallnx[ifmn]=>fmn, here fallnx is complete net|xch flux vector
# combining unknown (dependent), free, constrainded and groth fluxes
ifmn=match(nm_fmn, nm_fallnx)
"""%{
    "nm_fmn": join(", ", netan["vflux_meas"]["net"], '"', '"'),
    "fmn": join(", ", (netan["flux_measured"][fl]["val"]
        for fl in netan["vflux_meas"]["net"])).replace("nan", "NA"),
    "fmndev": join(", ", (netan["flux_measured"][fl]["dev"]
        for fl in netan["vflux_meas"]["net"])),
    })
    return {
        "o_meas": o_meas,
        "measures": measures,
        "o_mcumos": o_mcumos,
        "imcumo2i": imcumo2i,
    }

def netan2R_rcumo(netan, org, f, emu=False):
    # prepare reduced python systems
    #rAb=C13_ftbl.rcumo_sys(netan, emu)
    rAb=netan["rcumo_sys"] if ("rcumo_sys" in netan) else C13_ftbl.rcumo_sys(netan, emu)
    # full matrix is Ab=netan["cumo_sys"]

    # prune ordered cumomer list in reverse order
    # so that deleted item does not change the index
    # for the rest items to prune
    if "vrcumo" not in netan:
        netan["vrcumo"]=copy.deepcopy(netan["vcumo"])
        for i in range(len(netan["vrcumo"]),len(rAb["A"]),-1):
            # delete extra weight systems
            del(netan["vrcumo"][i-1])
        for (iw,cumol) in enumerate(netan["vrcumo"]):
            for i in range(len(cumol), 0, -1):
                i-=1
                if cumol[i] not in rAb["A"][iw]:
                    #print "prune", i, cumol[i];##
                    del(cumol[i])
    # prepare cumo2i
    # translate cumoname like A:7 to its index in R vector of cumomers
    rcumos=list(valval(netan["vrcumo"]))
    rcumo2i=dict((c,i+1) for (i,c) in enumerate(rcumos))
    emus=list(valval(netan.get("vemu", [])))
    emu2i=dict((c,i+1) for (i,c) in enumerate(emus))
    # composit cumomer vector incu=c(1,xi,xc)
    incu2i_b1=dict((c,i+2) for (i,c) in enumerate(list(netan["rcumo_input"][0].keys())+rcumos))
    # write code for reduced cumomer systems
    #netan2Abcumo_f(rAb["A"], rAb["b"],
    #    netan["vrcumo"], netan["input"], ff, netan["fwrv2i"], incu2i_b1, "fwrv2rAbcumo")
    #netan2Abcumo_sp("spAb_old", rAb["A"], rAb["b"],
    #    netan["vrcumo"], netan["input"], f, netan["fwrv2i"], incu2i_b1)
    #print("rab=", rAb["A"], "\n")
    netan2Abcumo_spr("spAbr", rAb["A"], rAb["b"],
        netan["vrcumo"], netan["input"], f, netan["fwrv2i"], incu2i_b1)
    #netan2j_rhs_f(rAb["A"], rAb["b"],
    #    netan["vrcumo"], netan["input"], ff, netan["fwrv2i"], rcumo2i, incu2i_b1, "frj_rhs")
    # write R constants and names
    f.write("""
# weight count
nb_rw=%(nb_rw)d
# cumomer count by weight
nb_rcumos=c(%(nb_rc)s)
nbc_cumos=c(0, cumsum(nb_rcumos))
# cumo names
nm_rcumo=c(%(nm_rcumo)s)
nm_list$rcumo=nm_rcumo
"""%{
    "nb_rw": len(rAb["A"]),
    "nb_rc": join(", ", (len(a) for a in rAb["A"])),
    "nm_rcumo": join(", ", valval(netan['vrcumo']), '"', '"'),
})
    f.write("""
if (case_i) {
   # check the coherence of metabolites/cumomers
   met_net=unique(matrix(unlist(strsplit(nm_rcumo, ":", fixed=TRUE)), nrow=2)[1,])
   net_pool=sort(setdiff(met_net, names(nm_poolall)))
   if (length(net_pool) > 0) {
      stop_mes("The following metabolites are internal in NETWORK section but not in METABOLITE_POOLS one:\\n", paste(net_pool, collapse="\\n"), file=fcerr)
   }
}
""")
    netan["rcumo2i"]=rcumo2i
    netan["emu2i"]=emu2i
    return {
        "rcumo2i": rcumo2i,
        "emu2i": emu2i,
        "rAb": rAb,
    }

def netan2R_cumo(netan, org, f):
    """netan2R_cumo(netan, org, f)->dict
    generate data structures for full cumomer matrices
    """
    # prepare cumo2i
    # translate cumoname like A:7 to its index in R vector of cumomers
    cumos=list(valval(netan["vcumo"]))
    cumo2i=dict((c,i+1) for (i,c) in enumerate(cumos))
    # composite cumomer vector
    incu2i_b1=dict((c,i+2) for (i,c) in enumerate(list(netan["cumo_input"][0].keys())+cumos))

    netan2Abcumo_spr("spAbr_f", netan["cumo_sys"]["A"], netan["cumo_sys"]["b"],
        netan["vcumo"], netan["input"], f, netan["fwrv2i"], incu2i_b1)
    # write R constants and names
    f.write("""
if (TIMEIT) {
   cat("cumo   : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}

# weight count
nb_w=%(nb_w)d

# cumomer count by weight
nb_cumos=c(%(nb_c)s)

# cumo names
nm_cumo=c(%(nm_cumo)s)
"""%{
    "nb_w": len(netan["cumo_sys"]["A"]),
    "nb_c": join(", ", (len(a) for a in netan["cumo_sys"]["A"])),
    "nm_cumo": join(", ", valval(netan['vcumo']), '"', '"'),
})
    netan["cumo2i"]=cumo2i
    return {
        "cumo2i": cumo2i,
    }

def netan2R_ineq(netan, org, f):
    """netan2R_ineq(netan, org, f)
    generate inequality code
    """
    # ex: netan["flux_inequal"]
    # {'net': [], 'xch': [('0.85', '>=', {'v2': '+1.'})]}
    tfallnx=netan["tfallnx"]
    f2dfcg_nx_f=netan["f2dfcg_nx_f"]
    nb_ineq=len(netan["flux_inequal"]["net"])+len(netan["flux_inequal"]["xch"])
    f.write("""
if (TIMEIT) {
   cat("ineq    : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
# prepare mi matrix and li vector
# such that mi*fallnx>=li corresponds
# to the inequalities given in ftbl file
nb_ineq=%(nb_ineq)s
mi=matrix(0., nrow=nb_ineq, ncol=nb_fallnx)
li=numeric(nb_ineq)
nm_i=c(c(%(nm_in)s), c(%(nm_ix)s))
dimnames(mi)=list(nm_i, nm_fallnx)
""" % {
   "nb_ineq": nb_ineq,
   "nm_in": join(", ", (join("", (ineq[0], ineq[1], join("+",
        ((str(fa)+"*" if fa != 1. else "")+fl
        for (fl,fa) in ineq[2].items()))))
        for ineq in netan["flux_inequal"]["net"]), p='"n:', s='"'),
   "nm_ix": join(", ", (join("", (ineq[0], ineq[1], join("+",
        ((str(fa)+"*" if fa != 1. else "")+fl
        for (fl,fa) in ineq[2].items()))))
        for ineq in netan["flux_inequal"]["xch"]), p='"x:', s='"'),
})
    
    for (i, ineq) in enumerate(netan["flux_inequal"]["net"]):
        f.write(
"""mi[%(i)s, nm_net[c(%(f)s)]]=%(sign)sc(%(coef)s)
li[%(i)s]=%(sign)s%(li)g
"""%{
    # as R inequality is always ">=" we have to inverse the sign for "<=" in ftbl
    "i": i+1,
    "sign": ("" if ineq[1]=="<=" or ineq[1]=="=<" else "-"),
    "f": join(", ", list(ineq[2].keys()), p='"', s='"'),
    "coef": join(", ", list(ineq[2].values())),
    "li": ineq[0],
    })
    for (i, ineq) in enumerate(netan["flux_inequal"]["xch"]):
        f.write(
"""mi[%(i)s, nm_xch[c(%(f)s)]]=%(sign)sc(%(coef)s)
li[%(i)s]=%(sign)s%(li)g
"""%{
    # as R inequality is always ">=" we have to inverse the sign for "<=" in ftbl
    "i": len(netan["flux_inequal"]["net"])+i+1,
    "sign": ("" if ineq[1]=="<=" or ineq[1]=="=<" else "-"),
    "f": join(", ", list(ineq[2].keys()), p='"', s='"'),
    "coef": join(", ", list(ineq[2].values())),
    "li": ineq[0],
    })

    nb_fdx=len(netan["vflux"]["xch"])
    nb_ffx=len(netan["vflux_free"]["xch"])
    f.write("""
# add standard limits on [df].xch [0;cupx]
nb_tmp=nrow(mi)
nb_fx=nb_flx+nb_ffx
if (nb_fx) {
   mi=rbind(mi, matrix(0, nrow=2*nb_fx, ncol=nb_fallnx))
   if (nb_flx)
      nm_i=c(nm_i, paste(nm_flx, ">=0", sep=""))
   if (nb_ffx)
      nm_i=c(nm_i, paste(nm_ffx, ">=0", sep=""))
   if (nb_flx)
      nm_i=c(nm_i, paste(nm_flx, "<=", cupx, sep=""))
   if (nb_ffx)
      nm_i=c(nm_i, paste(nm_ffx, "<=", cupx, sep=""))
   li=c(li, rep(0, nb_fx), rep(-cupx, nb_fx))
   mi[nb_tmp+(1:nb_fx),c(nm_flx, nm_ffx)]=diag(1., nb_fx)
   mi[nb_tmp+nb_fx+(1:nb_fx),c(nm_flx, nm_ffx)]=diag(-1., nb_fx)
}
""")
    
    nb_notrev=len(netan["notrev"])
    f.write("""
nm_inout=grep("^[^c]\\\\.", nm_net[c(%(nm_inout)s)], v=TRUE) # strip out constrained fluxes
nb_inout=length(nm_inout)
if (nb_inout > 0) {
   # add cinout low limits on inout net fluxes
   nb_tmp=nrow(mi)
   # explicit inequalities take precedence over generic ones
   # so eliminate inout fluxes which are already in inequalities
   nm_itmp=paste("n:.+<=", substring(nm_inout, 5), sep="")
   i=sapply(1:length(nm_itmp), function(k) {
      j=grep(nm_itmp[k], nm_i)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n", file=fclog)
      if (length(j)==0) {
         return(0)
      } else {
         return(k)
      }
   })
   i=i[i!=0]
   if (length(i) > 0) {
      nm_tmp=nm_inout[-i]
   } else {
      nm_tmp=nm_inout
   }
   len_tmp=length(nm_tmp)
   if (len_tmp > 0) {
      mi=rbind(mi, matrix(0, nrow=len_tmp, ncol=nb_fallnx))
      nm_i=c(nm_i, paste("inout ", nm_tmp, ">=", cinout, sep=""))
      mi[nb_tmp+(1:len_tmp), nm_tmp]=diag(1., len_tmp)
      li=c(li, rep(cinout, len_tmp))
   }
}
if (clownr!=0.) {
   # add low limits on net >= clownr for not reversible reactions
   nb_tmp=nrow(mi)
   nm_tmp=nm_net[c(%(nm_notrev)s)]
   # explicit inequalities take precedence over generic ones
   # so eliminate notrev fluxes which are already in inequalities
   nm_itmp=paste("n:.+<=", substring(nm_tmp, 5), sep="")
   i=sapply(1:length(nm_itmp), function(k) {
      j=grep(nm_itmp[k], nm_i)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n", file=fclog)
      if (length(j)==0) {
         return(0)
      } else {
         return(k)
      }
   })
   i=i[i!=0]
   if (length(i) > 0) {
      nm_tmp=nm_tmp[-i]
   }
   # search for inout too
   nm_itmp=paste("inout ", nm_tmp, ">=", sep="")
   i=sapply(1:length(nm_itmp), function(k) {
      j=grep(nm_itmp[k], nm_i, fix=TRUE)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n", file=fclog)
      if (length(j)==0) {
         return(0)
      } else {
         return(k)
      }
   })
   i=i[i!=0]
   if (length(i) > 0) {
      nm_tmp=nm_tmp[-i]
   }

   len_tmp=length(nm_tmp)
   if (len_tmp > 0) {
      mi=rbind(mi, matrix(0, nrow=len_tmp, ncol=nb_fallnx))
      nm_i=c(nm_i, paste(nm_tmp, ">=", clownr, sep=""))
      mi[nb_tmp+(1:len_tmp), nm_tmp]=diag(1., len_tmp)
      li=c(li, rep(clownr, len_tmp))
   }
}
nb_fn=nb_fln+nb_ffn
if (cupn != 0 && nb_fn > 0) {
   # add absolute upper limits on -cupn <= [df].net <= cupn for net fluxes
   # explicit inequalities take precedence over generic ones
   # so eliminate net fluxes which are already in inequalities
   ## proceed n:smth>=flux
   nm_tmp=c(nm_ffn, nm_fln) # all not fixed net fluxes
   nm_itmp=paste("n:.+>=", substring(nm_tmp, 5), sep="")
   i=sapply(vgrep(nm_itmp, nm_i), length)
   i=which(i!=0)
   if (length(i) > 0) {
      nm_tmp=nm_tmp[-i]
   }
   len_tmp=length(nm_tmp)
   if (len_tmp > 0) {
      nb_tmp=nrow(mi)
      mi=rbind(mi, matrix(0, nrow=len_tmp, ncol=nb_fallnx))
      nm_i=c(nm_i, paste0(nm_tmp, "<=", cupn))
      li=c(li, rep(-cupn, len_tmp))
      mi[nb_tmp+(1:len_tmp),nm_tmp]=diag(-1., len_tmp)
   }
   ## proceed n:smth<=flux
   nm_tmp=c(nm_ffn, nm_fln) # all not fixed net fluxes
   nm_itmp=paste("n:.+<=", substring(nm_tmp, 5), sep="")
   i=sapply(vgrep(nm_itmp, nm_i), length)
   i=which(i!=0)
   if (length(i) > 0) {
      nm_tmp=nm_tmp[-i]
   }
   len_tmp=length(nm_tmp)
   if (len_tmp > 0) {
      nb_tmp=nrow(mi)
      mi=rbind(mi, matrix(0, nrow=len_tmp, ncol=nb_fallnx))
      nm_i=c(nm_i, paste0(nm_tmp, ">=", -cupn))
      li=c(li, rep(-cupn, len_tmp))
      mi[nb_tmp+(1:len_tmp),nm_tmp]=diag(1., len_tmp)
   }
}
#browser()
"""%{
#   "nb_notrev": len([fli for (fli,t,nxi) in tfallnx
#      if nxi=="n" and t!="c" and fli in netan["notrev"]]),
   "nm_notrev": join(", ", netan["notrev"], p='"', s='"'),
   "nm_inout": join(", ", netan["flux_inout"], p='"', s='"'),
})

    f.write("nb_ineq=NROW(li);\n")
    f.write("""
dimnames(mi)=list(nm_i, nm_fallnx)
names(li)=nm_i
# prepare ui matrix and ci vector for optimisation
# ui%*%param-ci>=0
# it is composed of explicite inequalities from ftbl
# and permanent inequalities 0<=xch<=0.999 and scale>=0

# constraints such that ui%*%param-ci>=0
# first flux part
ui=mi%*%(md%*%invAfl%stm%p2bfl+mf)
mic=(md%*%invAfl%*%(c2bfl%stm%fc+cnst2bfl) + mc%*%fc)
ci=as.numeric(li-mi%*%mic)
""")
    f.write("""
# finaly, metab part
uip=matrix(0., %(nb_ip)d, ncol=nb_poolf)
colnames(uip)=nm_poolf
cip=c()
# ind: irow, metab, coef, rhs, name
uip_ind=c(
"""%{
    "nb_ip": len(netan["metab_inequal"]),
})
    st=""
    for (i, (rhs, comp, d, name)) in enumerate(netan["metab_inequal"]):
        si=-1. if comp == ">=" or comp == "=>" else 1.
        rhs*=si
        for (m, coef) in d.items():
            # ind4 (str): irow, metab, coef, rhs, name
            coef*=si
            if netan["met_pools"][m] > 0:
                rhs-=netan["met_pools"][m]*coef
                m=""
                coef=0.
            st=st+"""   "%d", "%s", "%g", "%g", "%s",
"""%(i+1, m, coef, rhs, name)
            rhs=0.
    f.write("""%s
)
if (length(uip_ind) > 0) {
   uip_ind=matrix(uip_ind, byrow=TRUE, ncol=5L)
} else {
   uip_ind=matrix(0, 0L, 5L)
}
colnames(uip_ind)=c("irow", "metab", "coef", "rhs", "name")
"""%st[:-2]+"\n")

    f.write("""
if (nrow(uip_ind) > 0) {
   # rhs are summed up for the same irow by aggregate()
   irow=as.integer(uip_ind[,"irow"])
   cip=aggregate(as.double(uip_ind[,"rhs"]), by=list(irow), sum)[,"x"]
   for (i in seq_len(nrow(uip_ind))) {
      row=uip_ind[i,]
      if (nchar(row["metab"])==0) {
         next
      }
      uip[irow[i], nm_poolf[row["metab"]]]=as.double(row["coef"])
   }
   if (nrow(uip) > 0) {
      rownames(uip)=paste("m:", uip_ind[pmatch(seq_len(nrow(uip)), irow),"name"], sep="")
   }
}
names(cip)=rownames(uip)

#browser() # before null inequality removing
# remove all zero rows in ui (constrained fluxes with fixed values)
# find zero indexes
#print(dim(ui))
if (ncol(ui)) {
   zi=apply(ui,1,function(v){return(max(abs(v))<=1.e-14)})
} else {
   # remove all flux inequalities as there is no free fluxe
   zi=rep(TRUE, nrow(ui))
}

if (!all(ci[zi]<=1.e-10)) {
   cat("The following constant inequalities are not satisfied:\\n", file=fcerr)
   cat(nm_i[zi][ci[zi]>1.e-10], sep="\\n", file=fcerr)
   cat("They are simply ignored.\\n", file=fcerr)
   #stop_mes("", file=fcerr)
}
ui=ui[!zi,,drop=FALSE]
ci=ci[!zi]
nm_i=nm_i[!zi]

# complete ui by zero columns corresponding to scale params
if (nb_sc_tot) {
   ui=cBind(as.matrix(ui), matrix(0., NROW(ui), nb_sc_tot))
   # complete ui by scales >=0
   ui=rBind(ui, cBind(matrix(0, nb_sc_tot, nb_ff), diag(1, nb_sc_tot)))
   ci=c(ci,rep(0., nb_sc_tot))
   nm_i=c(nm_i, paste(nm_par[nb_ff+seq_len(nb_sc_tot)], ">=0", sep=""))
   rownames(ui)=nm_i
   names(ci)=nm_i
}

# remove redundant inequalities
#browser()
nb_i=nrow(ui)
ired=c()
if (nb_i > 1L) {
   for (i in 1L:(nb_i-1L)) {
      nmref=nm_i[i]
      for (j in setdiff((i+1L):nb_i, ired)) {
         if (all(ui[j,]==ui[i,]) && ci[i]==ci[j]) {
            # redundancy
            cat("inequality '", nm_i[j], "' redundant with '", nmref, "' is removed.\n", sep="", file=fclog)
            ired=c(ired, j)
         }
      }
   }
}
if (!is.null(ired)) {
   # remove all ired inequalities
   ui=ui[-ired,,drop=FALSE]
   ci=ci[-ired]
   nm_i=nm_i[-ired]
}

# metabolite equalities
ep=matrix(0., %(nb_ep)d, ncol=nb_poolf)
cp=c()
colnames(ep)=nm_poolf
# ind: irow, metab, coef, rhs, name
ep_ind=c(
"""%{
    "nb_ep": len(netan["metab_equal"]),
})
    st=""
    for (i, (rhs, d, name)) in enumerate(netan["metab_equal"]):
        for (m, coef) in d.items():
            # ind4 (str): irow, metab, coef, rhs, name
            if netan["met_pools"][m] > 0:
                rhs-=netan["met_pools"][m]*coef
                m=""
                coef=0.
            st=st+"""   "%d", "%s", "%g", "%g", "%s",
"""%(i+1, m, coef, rhs, name)
            rhs=0.
    f.write("""%s
)
if (length(ep_ind) > 0) {
   ep_ind=matrix(ep_ind, byrow=TRUE, ncol=5L)
} else {
   ep_ind=matrix(0, 0L, 5L)
}
colnames(ep_ind)=c("irow", "metab", "coef", "rhs", "name")
"""%st[:-2]+"\n")

    f.write("""
if (nrow(ep_ind) > 0) {
   # rhs are summed up for the same irow by aggregate
   irow=as.integer(ep_ind[,"irow"])
   cp=aggregate(as.double(ep_ind[,"rhs"]), by=list(irow))[,"x"]
   for (i in seq_len(nrow(ep_ind))) {
      row=ep_ind[i,]
      if (nchar(row["metab"])==0) {
         next
      }
      ep[irow[i], nm_poolf[row["metab"]]]=as.double(row["coef"])
   }
   if (nrow(ep) > 0) {
      rownames(ep)=paste("m:", ep_ind[pmatch(seq_len(nrow(ep)), irow),"name"], sep="")
   }
}
ep=as.matrix(ep)
names(cp)=rownames(ep)
""")
