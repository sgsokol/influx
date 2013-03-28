#!/usr/bin/env python

"""Module for translation of .ftbl file to R code"""

# 2012-02-21 sokol@insa-toulouse.fr : cumomer matrices and rhs from sparse matrices
#                                     (without fortran code)
# 2009-09-14 sokol@insa-toulouse.fr : flux.[net|xch] -> [dfcg].[nx].flux
#                                     flux.[fwd|rev] -> [fwd|rev].flux
# 2008-12-08 sokol@insa-toulouse.fr : added netan2Rinit()
# 2008-11-25 sokol@insa-toulouse.fr : adapted for reduced cumomer list
# 2008-09-19 sokol@insa-toulouse.fr : initial release
# Copyright 2011, INRA

import time
import copy
import os
import sys
from operator import itemgetter
from itertools import groupby

global DEBUG
me=os.path.realpath(sys.argv[0])
dirx=os.path.dirname(me)
sys.path.append(dirx)

from tools_ssg import *
import C13_ftbl

def netan2Abcumo_spr(varname, Al, bl, vcumol, minput, f, fwrv2i, incu2i_b1):
    """
    Transform cumomer linear sytems collection (from ftbl file)
    to a R code calculating sparse matrix A and vector b
    in A*x+b=0 for a given weight of fragment iw (index in resulting list)
    Flux vector fl of all fwd. and rev. fluxes is known from an
    R environement.
    
    Resulting code is a list sprAb indexed by cumomer weight
    (cf. generated R comments for details on sprAb)
    cumomer vector incu=c(1, xi, xl), xi - input cumomers, xl - lighter cumomers.
    
    incu2i_b1 gives i in incu from cumomer name. i=1 corresponds to the constant 1.
    Difference wrt netan2Abcumo_sp is that pure R code is used
    (@i, @p and @x slots are those from Matrix::dgCMatrix class).
    No need for Fortran compiler.
    """
    #2012-02-08 sokol
    
    nb_cumu=cumsum(len(l) for l in vcumol)
    f.write(
    """
# sparse matrix static parts
# $varname fields:
#  ind_fa - flux index in a_pre@x=fwrv[ind_fa]
#  a_pre - sparse matrix whose colsum gives the at@x vector
#  ind_fb - flux index in b_pre@x=fwrv[ind_fb1]*x[ind_x1]*x[ind_x2]
#  ind_x1 - cumomer index in b_pre@x=fwrv[ind_fb]*x[ind_x1]*x[ind_x2]
#  ind_x2 - cumomer index in b_pre@x=fwrv[ind_fb]*x[ind_x1]*x[ind_x2]
#  b_pre - sparse matrix whose colsum gives b@x

#  tA - unsigned sparse transpose of cumomer A matrix (off-diagonal part)
#  b - unsigned sparse vector of right hand side
#  
#  f2ax indexes of fluxes to calculate tA@x slot (cf. fortran
#   f2ax() subroutine in cumo.f)
#  bfpr - unsigned sparse rhs of Ax=b.
#   [(irow, iflux, incu_index1, incu_index2),...]
#   to calculate terms like flux*x1*x2. When there is only one
#   term x, x2 is set to 1
#  ta_fx - unsigned sparse matrix t(dA_df*x)
#   @i run through the fluxes contibuting to a row of dA_df*x
#  x2ta_fx - 1-based indexes of x0=c(0,x) for dA_df*x which is
#   comosed of differences (x_diag - x_off-diag) and (x_diag)
#   when flux is in b term. In this case the index of x_off-diag is
#   set to 1 (i.e. 0-value in x0).
#   [irow, ixoff] irow is actually ixdiag
#  tb_f - unsigned sparse matrix t(db_df) with indexes ix1, ix2 for product as in bfpr
#   @i is running along fluxes
#  x2tb_f - 2 row matrix with indexes ix1, ix2 for product as in bfpr
#  tb_x - unsigned sparse transpose of db_dx over x in lighter weights
#   @i runs over lighter cumomers

if (TIMEIT) {
   cat("spAbr   : ", date(), "\\n", sep="")
}

nb_fwrv=%(n)d
%(var)s=list()
"""%{
    "var": varname,
    "n": len(fwrv2i),
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
            raise Exceptiopn("wrongCumomerNumber")
        l_ia=[]; # list of non zero off-diagonal elements in A / row
        l_ib=[]; # list of non zero elements in b / row
        nb_maxfa=0; # how many fluxes in an off-diagonal term in a
        #nb_maxfb=0; # how many fluxes in a term in b
        #nb_ax=0
        for irow in xrange(ncumo):
            cr=cumos[irow]
            row=A[cr]
            # atuple is list of (icumo, list(fluxes))
            atuple=[(c2i[c], [fwrv2i[fl] for fl in row[c]])
                for c in row] #cumos if c in row and c!=cr]
            #if atuple:
            #    nb_maxfa=max(nb_maxfa, max(len(lf) for (ic, lf) in atuple))
            #elif cr not in b:
            #    raise Exception("Empty row in cumomer matrix, weight=%d (base 1), cumo=%s"%(w, cr))
            # btuple is list of (iflux, icumo1, icumo2)
            if cr in b:
                btuple=[(fwrv2i[fl], incu2i_b1[l[0]], (incu2i_b1[l[1]] if len(l)==2 else 1))
                    for (fl, d) in b[cr].iteritems()
                    for (i,l) in d.iteritems()]
                #nb_maxfb=max(nb_maxfb, len(btuple))
            else:
                btuple=[]
            # one list per row
            l_ia.append(atuple)
            #nb_ax+=len(atuple)
            l_ib.append(btuple)
        f.write(
"""
if (TIMEIT) {
   cat("weight %(w)d: ", date(), "\\n", sep="")
}
w=%(w)d
nb_c=%(nbc)d
ba_x=%(ba_x)d; # base of cumomer indexes in incu vector
l=list()
l$w=w
l$nb_c=nb_c
l$nb_fwrv=nb_fwrv
l$nb_cl=%(ncucumo)d # number of lighter cumomers
if (nb_c > 0) {
   # matrix a
   ind_a=matrix(as.integer(c(%(ind_a)s)), ncol=3, byrow=T)
   colnames(ind_a)=c("indf", "ir0", "ic0")
   l$ind_a=ind_a
   
   # vector b
   ind_b=matrix(as.integer(c(%(ind_b)s)), ncol=4, byrow=T)
   colnames(ind_b)=c("indf", "indx1", "indx2", "irow")
   # put the lowest of ix1 and ix2 to ix2
   ima=pmax(ind_b[,"indx1"], ind_b[,"indx2"])
   imi=pmin(ind_b[,"indx1"], ind_b[,"indx2"])
   ind_b[,"indx1"]=ima
   ind_b[,"indx2"]=imi
   l$ind_b=ind_b
   
   # jacobian b_x
   i=ind_b[,"indx2"]!=1 # exclude from derivation plain input entries
   tmp=ind_b[i,,drop=F]
   
   # term of d/d_x1 ( is garanted to be internal, not input cumomer)
   # => indx remain in place in indx2, ind_store remain in column indx1
   
   # term of d/d_x2 (x2 can be an input cumomer => no derivation)
   # => indx is taken from indx1 and goes to indx2, while ind_store goes to indx1
   i=which(tmp[,"indx2"] > ba_x)
   if (length(i)) {
      ind_bx=rbind(tmp, tmp[i,c(1,3,2,4)])
   } else {
      ind_bx=tmp
   }
   colnames(ind_bx)=c("indf", "ic1", "indx", "irow")
   ind_bx[,"ic1"]=ind_bx[,"ic1"]-ba_x
   l$ind_bx=ind_bx
}
%(var)s[[w]]=l
"""%{
   "var": varname,
   "w": w,
   "nbc": ncumo,
   "ncucumo": ncucumo,
   "ba_x": ba_x,
   "ind_a": join(", ", valval((ifl, ir, ic)
      for (ir, lt) in enumerate(l_ia)
      for (ic, lf) in lt
      for ifl in lf)),
   "ind_b": join(", ", valval((ifl, i1, i2, ir+1)
       for (ir, lt) in enumerate(l_ib)
       for (ifl, i1, i2) in lt
   )),
})
        ba_xw+=ncumo
        ncucumo+=ncumo

def netan2Rinit(netan, org, f, fullsys, emu=False, ropts=[]):
    r"""Write R code for initialization of all variables before
    cumomer system resolution by khi2 minimization.
    Args:
     netan: a collection of parsed ftbl information
     f: R code output pointer
     fullsys (logical): write a code for the full or only reduced cumomer system
     emu (logical): write equations in EMU framework or cumomer (default)
     ropts: list of items "param=value" to be written as is in R file.
    Return:
     a dictionnary with some python variables:
        * "measures": measures,
        * "o_mcumos": o_mcumos,
        * "cumo2i": cumo2i,
        * ...
    """
    global DEBUG
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
    #    param2fl_x - translate param to flux and cumomer vector (initial approximation)
    #    cumo_cost - cost function (khi2)
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

    # header
    f.write("# This is an automatically generated R code. Don't edit.\n")
    f.write("# Generated by \n# "+join(" ", sys.argv)+"\n# at "+time.ctime()+".\n")
    f.write("""
# Copyright INRA, France.
""")
    #if DEBUG:
    #    pdb.set_trace()
    res={}
    f.write("""
fcerr=file("%(errfile)s", "ab")
sink(fcerr, split=F, type="message")
sink("%(logfile)s", append=T, split=F, type="output")
options(warn=1)
suppressPackageStartupMessages(library(bitops))
suppressPackageStartupMessages(library(nnls)); # for non negative least square
suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # for sparse matrices
options(Matrix.quiet=TRUE)
suppressPackageStartupMessages(library(parallel))

# get some common tools
source("%(dirx)s/tools_ssg.R")
source("%(dirx)s/nlsic.R")
source("%(dirx)s/kvh.R")

# default options
version=F
noopt=F
noscale=F
meth="nlsic"
fullsys=F
emu=F
irand=F
sens=""
cupx=0.999
cupn=1.e3
cupp=1.e5
clownr=0
cinout=0
clowp=1.e-8
np=0
ln=F
zc=-.Machine$double.xmax
fseries=""
iseries=""
seed=-.Machine$integer.max
excl_outliers=F
DEBUG=F
TIMEIT=F
prof=F

# get runtime arguments
%(ropts)s

# synonymous
myver=version
optimize=!noopt
method=meth
sensitive=sens
least_norm=ln
initrand=irand

vernum="%(vernum)s"
if (myver) {
   cat("%(prog)s %(vernum)s\\n")
   q("no")
}
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
   stop(paste("Option '--cupx N' must have N in the interval [0,1]\n",
      "Instead, the value ", cupx, " si given.", sep=""))
}
if (cinout < 0) {
   stop(paste("Option '--cinout N' must have N non negative\n",
      "Instead, the value ", cinout, " si given.", sep=""))
}
# minimization method
validmethods=list("BFGS", "Nelder-Mead", "SANN", "ipopt", "nlsic")
if (! method %%in%% validmethods) {
   warning(paste("method", method, "is not known. 'nlsic' is used instead."))
   method="nlsic"
}
if ( method == "ipopt") {
   installed=suppressPackageStartupMessages(library(ipoptr, logical.return=T))
}

if (np) {
   options(cores=np)
}
lsi_fun=lsi
if (least_norm) {
   lsi_fun=lsi_ln
}
if (zc==-.Machine$double.xmax) {
   # no zero scrossing to apply
   zerocross=F
} else {
   if (zc < 0.) {
      stop(paste("Zero crossing value ZC must be non negative, instead ", zc, " is given.", sep=""))
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
opts=commandArgs()
# end command line argument proceeding

# get some cumomer tools
source("%(dirx)s/opt_cumo_tools.R")
if (TIMEIT) {
   cat("rinit   : ", date(), "\n", sep="")
}

# R profiling
if (prof) {
   Rprof("%(proffile)s")
}

nm_list=list()
"""%{
    "dirx": escape(dirx, "\\"),
    "vernum": file(os.path.join(dirx, "influx_version.txt"), "r").read().strip(),
    "proffile": escape(f.name[:-1]+"Rprof", "\\"),
    "prog": os.path.basename(f.name),
    "ropts": join("\n", ropts)[1:-1],
    "logfile": escape(f.name[:-1]+"log", "\\"),
    "errfile": escape(f.name[:-1]+"err", "\\"),
})
    netan2R_fl(netan, org, f)
    d=netan2R_rcumo(netan, org, f)
    res.update(d)
    f.write("""
# input cumomer vector
xi=c(%(xi)s)
nm_xi=c(%(nm_xi)s)
nm_list$xi=nm_xi
nb_xi=length(xi)
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
   nb_emus=nb_rcumos*(1:nb_rw+1)
   nb_f$emus=nb_emus
   nm_list$emu=nm_emu
   nm_x=nm_emu
   nb_x=nb_emus
   xiemu=c(%(xiemu)s)
   nm_xiemu=c(%(nm_xiemu)s)
   nm_list$xiemu=nm_xiemu
   nb_xiemu=length(xiemu)
   nb_f$xiemu=nb_xiemu
   nb_f$xi=nb_xiemu
   nb_xi=nb_xiemu
   nm_inp=nm_xiemu
   xi=xiemu
   nm_inemu=c("one", nm_xiemu, nm_emu)
   nm_inlab=nm_inemu
   spa=spr2emu(spAbr, nm_incu, nm_inemu, nb_f)
}
# composite labeling vector incu c(1, xi, xc) names
nm_inlab=c("one", nm_inp, nm_x); # the constant 1 has name "one"
nm_list$x=nm_x
nb_f$x=nb_x
"""%{
    "xi": join(", ", netan["rcumo_input"].values()),
    "nm_xi": join(", ", netan["rcumo_input"].keys(), '"', '"'),
    "xiemu": join(", ", netan["emu_input"].values()),
    "nm_xiemu": join(", ", netan["emu_input"].keys(), '"', '"'),
    "nm_emu": join(", ", valval(netan['vemu']), '"', '"'),
    #"so": ("dll" if sys.platform in ("win32","cygwin") else "dylib" if sys.platform == "darwin" else "so"),
})
    if fullsys:
        d=netan2R_cumo(netan, org, f)
        res.update(d)
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
      label=%(meas_l)s
   ),
   equations=list(
      equalities=%(eqe)s,
      inequalities=%(eqi)s
   ),
   cumomer=list(
      full=c(%(lncumo)s),
      reduced=c(%(lnrcumo)s)
   )
)
# carbon length of metabolites
clen=c(%(clen)s)
names(clen)=c(%(nm_metab)s)
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
    "meas_m": len(netan["measures"]["mass"]["vec"]),
    "meas_p": len(netan["measures"]["peak"]["vec"]),
    "meas_l": len(netan["measures"]["label"]["vec"]),
    "eqe": len(netan["flux_equal"]["net"])+len(netan["flux_equal"]["xch"]),
    "eqi": len(netan["flux_inequal"]["net"])+len(netan["flux_inequal"]["xch"]),
    "lncumo": ",".join(str(len(a)) for a in netan["cumo_sys"]["A"]),
    "lnrcumo": ",".join(str(len(a)) for a in netan["rcumo_sys"]["A"]),
    
    "clen": join(",", netan["Clen"].values()),
    "nm_metab": join(",", netan["Clen"].keys(), '"', '"')
    })
    return res

def netan2R_fl(netan, org, f):
    """netan2R_fl(netan, org, f)
    generate R code for flux part
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
        for (f,i) in netan['vflux_free']['net2i'].iteritems())
    ffx2iprm=dict(("f.x."+f,(i+1+nb_ffn))
        for (f,i) in netan['vflux_free']['xch2i'].iteritems())

    # prepare fwrv2i
    fwrv2i=dict((f,i+1) for (f,i) in netan["vflux_fwrv"]["fwrv2i"].iteritems())
    nb_fwrv=len(netan["vflux_fwrv"]["fwrv2i"])

    # make tuple for complete flux vector d,f,c
    # (name,"d|f|c|g","n|x")
    tfallnx=zip(
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
            )
    netan["f2dfcg_nx_f"]={
       "net": dict((fl, t+".n."+fl) for (fl,t,nx) in tfallnx if nx=="n"),
       "xch": dict((fl, t+".x."+fl) for (fl,t,nx) in tfallnx if nx=="x"),
    }

    f.write("""
if (TIMEIT) {
   cat("r_flux  : ", date(), "\n", sep="")
}

# custom functions
# produce fw-rv fluxes from fallnx
fallnx2fwrv=function(fallnx) {
   n=length(fallnx)
   # extract and reorder in fwrv order
   net=fallnx[c(%(inet2ifwrv)s)]
   xch=fallnx[c(%(ixch2ifwrv)s)]
   # expansion 0;1 -> 0;+inf of xch (second half of fallnx)
   xch=xch/(1-xch)
   # fw=xch-min(-net,0)
   # rv=xch-min(net,0)
   fwrv=c(xch-pmin(-net,0),xch-pmin(net,0))
   if (DEBUG) {
      n=length(fwrv)
      library(MASS)
      names(fwrv)=nm_fwrv
      write.matrix(fwrv, file="dbg_fwrv.txt", sep="\\t")
   }
   return(fwrv)
}
""" % {
#dyn.load("%(sofile)s")
        "inet2ifwrv": join(", ", (1+
        (netan["vflux"]["net2i"][fl[4:]] if fl[4:] in netan["vflux"]["net2i"]
        else nb_fln+netan["vflux_free"]["net2i"][fl[4:]] if fl[4:] in netan["vflux_free"]["net2i"]
        else nb_fln+nb_ffn+netan["vflux_constr"]["net2i"][fl[4:]] if fl[4:] in netan["vflux_constr"]["net2i"]
        else nb_fln+nb_ffn+nb_fcn+netan["vflux_growth"]["net2i"][fl[4:]])
        for fl in netan["vflux_fwrv"]["fwrv"][:nb_fwrv/2])),
        "ixch2ifwrv": join(", ", (1+len(netan["vflux_fwrv"]["fwrv"])/2+
        (netan["vflux"]["xch2i"][fl[4:]] if fl[4:] in netan["vflux"]["xch2i"]
        else nb_flx+netan["vflux_free"]["xch2i"][fl[4:]] if fl[4:] in netan["vflux_free"]["xch2i"]
        else nb_flx+nb_ffx+netan["vflux_constr"]["xch2i"][fl[4:]] if fl[4:] in netan["vflux_constr"]["xch2i"]
        else nb_flx+nb_ffx+nb_fcx+netan["vflux_growth"]["xch2i"][fl[4:]])
        for fl in netan["vflux_fwrv"]["fwrv"][:nb_fwrv/2])),
    })
    # auxiliary dict for edge-flux coupling
    f2edge=dict()
    for (fl,lr) in netan["sto_r_m"].iteritems():
        if len(lr["left"])==1 and len(lr["right"])==1:
           f2edge[fl]=[lr["left"][0]+" ("+fl+") "+lr["right"][0]]
        else:
           f2edge[fl]=[]
           for m in lr["left"]:
               f2edge[fl].append(m+" ("+fl+") "+fl)
           for m in lr["right"]:
               f2edge[fl].append(fl+" ("+fl+") "+m)
    #sys.stderr.write(str(f2edge)+"\n")
    #sys.stderr.write(str(netan["f2dfcg_nx_f"]["net"])+"\n")
    f.write("""
# fwd-rev flux names
nm_fwrv=c(%(nm_fwrv)s)

# net-xch flux names
nm_fallnx=c(%(nm_fallnx)s)

# edge to netflux name translator
edge2fl=c(%(edge2fl)s)

# initialize the linear system Afl*flnx=bfl (0-weight cumomers)
# unknown net flux names
nm_fln=c(%(nm_fln)s)
nb_fln=length(nm_fln)
# unknown xch flux names
nm_flx=c(%(nm_flx)s)
nb_flx=length(nm_flx)
nm_fl=c(nm_fln, nm_flx)
nb_fl=nb_fln+nb_flx
# gather flux names in a list
nm_list$flnx=nm_fl
nm_list$fallnx=nm_fallnx
nm_list$fwrv=nm_fwrv
# flux matrix
nb_flr=%(nb_flr)d
if (nb_fl) {
   Afl=matrix(0, nrow=nb_flr, ncol=nb_fl)
"""%{
    "nb_flr": len(netan["Afl"]),
    "nm_fwrv": join(", ", netan["vflux_fwrv"]["fwrv"], '"', '"'),
    "nm_fallnx": join(", ", (join(".", (t[1],t[2],t[0])) for t in tfallnx), '"', '"'),
    "nm_fln": join(", ", netan["vflux"]["net"], '"d.n.', '"'),
    "nm_flx": join(", ", netan["vflux"]["xch"], '"d.x.', '"'),
    "edge2fl": join(", ", ('"'+e+'"="'+netan["f2dfcg_nx_f"]["net"][fl]+'"' for (fl,l) in f2edge.iteritems() for e in l)),
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
if (DEBUG) {
   library(MASS)
   write.matrix(Afl, file="dbg_Afl.txt", sep="\\t")
}
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
nm_ff=c(nm_ffn, nm_ffx)
nm_list$ff=nm_ff
nb_param=length(param)
if (initrand) {
   param=runif(nb_param)
}
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

# total flux vector fallnx dimension
nb_fallnx=nb_fl+nb_ff+nb_fc+nb_fgr+nb_fgr
nb_fwrv=nb_fallnx

# all flux cardinals
nb_f=list(nb_fln=nb_fln, nb_flx=nb_flx, nb_fl=nb_fl,
   nb_ffn=nb_ffn, nb_ffx=nb_ffx, nb_ff=nb_ff,
   nb_fcn=nb_fcn, nb_fcx=nb_fcx, nb_fc=nb_fc,
   nb_fallnx=nb_fallnx, nb_fwrv=nb_fwrv,
   nb_fgr=nb_fgr,
   include_growth_flux=%(inc_gr_f)s,
   mu=%(mu)s)
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
p2bfl=Matrix(0., nrow=nb_flr, ncol=nb_ff)
# replace c.[nx].flx by corresponding fc coefficient
c2bfl=Matrix(0., nrow=nb_flr, ncol=nb_fc)
# variable growth fluxes
g2bfl=Matrix(0., nrow=nb_flr, ncol=nb_fgr)
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
        row["f"]=dict((k,v) for (k,v) in item.iteritems() if k[0:2]=="f.")
        row["c"]=dict((k,v) for (k,v) in item.iteritems() if k[0:2]=="c.")
        row["g"]=dict((k,v) for (k,v) in item.iteritems() if k[0:2]=="g.")
        f.write("\n")
        if row["f"]:
            f.write("p2bfl[%(i)d, c(%(if)s)]=c(%(rowf)s);\n"%\
                {"i": i+1,
                "if": join(", ", row["f"].keys(), p='"', s='"'),
                "rowf": join(", ", row["f"].values()),
                })
        if row["c"]:
            f.write("c2bfl[%(i)d, c(%(ic)s)]=c(%(rowc)s);\n"%\
                {"i": i+1,
                "ic": join(", ", row["c"].keys(), p='"', s='"'),
                "rowc": join(", ", row["c"].values()),
                })
        if row["g"]:
            f.write("g2bfl[%(i)d, c(%(ig)s)]=c(%(rowg)s);\n"%\
                {"i": i+1,
                "ig": join(", ", row["g"].keys(), p='"', s='"'),
                "rowg": join(", ", row["g"].values()),
                })
        if row["cnst"]:
            f.write("cnst2bfl[%(i)d]=%(rowcnst)s;\n"%{"i": i+1, "rowcnst": row["cnst"],})
    f.write("""
bp=as.numeric(c2bfl%*%fc+cnst2bfl)
""")

    f.write("""
if (TIMEIT) {
   cat("Afl qr(): ", date(), "\\n", sep="")
}

qrAfl=qr(Afl, LAPACK=T)
d=diag(qrAfl$qr)
qrAfl$rank=sum(abs(d)>=abs(d[1]*1.e-10))
#browser()
if (nrow(Afl) != ncol(Afl)) {
   #write.table(Afl)
   if (nrow(Afl) < ncol(Afl)) {
      mes=paste("Candidate(s) for free flux(es):\\n",
         paste(colnames(Afl)[-qrAfl$pivot[1:nrow(Afl)]], collapse="\\n"), sep="")
   } else {
      nextra=nrow(Afl)-ncol(Afl)
      comb=combn(nb_ffn, nextra)
      i=which.min(apply(comb, 2, function(i)kappa(cBind(Afl, p2bfl[,i]))))[1]
      i=comb[,i]
      ka=kappa(cBind(Afl, p2bfl[,i]))
      if (ka!=Inf) {
         prop=paste("Proposal to delcare dependent flux(es) is:\\n",
            paste(nm_ffn[i], collapse="\\n"), "\\n", sep="")
      } else {
         # test constraint candidate
         if (nb_fcn > 0) {
            comb=combn(nb_fcn, nextra)
            i=which.min(apply(comb, 2, function(i)kappa(cBind(Afl, c2bfl[,i]))))[1]
            i=comb[,i]
            ka=kappa(cBind(Afl, c2bfl[,i]))
            if (ka!=Inf) {
               prop=paste("Proposal to delcare dependent flux(es) is:\\n",
               paste(nm_fcn[i], collapse="\\n"), "\\n", sep="")
            } else {
               prop="No proposal for dependent fluxes could be made.\\n"
            }
         } else {
            prop="No proposal for dependent fluxes could be made.\\n"
         }
      }
      prop=paste(prop, "Condition number of stoechiometric matrix would be ", ka, "\\n", sep="")
      mes=paste("There is (are) probably ", nextra,
         " extra free flux(es) among the following:\\n",
         paste(nm_ffn, collapse="\\n"), "\\n",
         prop,
         "\\nAnother option could be an elimination of some equalities.",
         sep="")
   }
   stop(paste("Flux matrix is not square: (", nrow(Afl), "eq x ", ncol(Afl), "unk)\\n",
      "You have to change your choice of free fluxes in the '%(n_ftbl)s' file.\\n",
      mes, sep=""))
}

# make sure that free params choice leads to not singular matrix
if (qrAfl$rank != nb_fl) {
   #write.table(Afl)
   # make a suggestion of new free fluxes
   A=cBind(Afl, -p2bfl, -c2bfl)
   colnames(A)=c(colnames(Afl), nm_ff, nm_fc)
   qa=qr(A, LAPACK=T)
   d=diag(qa$qr)
   qa$rank=sum(abs(d)>=abs(d[1]*1.e-10))
   
   mes=paste("Dependent flux matrix is singular.\\n",
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
   cat(mes)
   stop(mes)
}

# inverse flux matrix
invAfl=solve(qrAfl)

""" % {
    "n_ftbl": escape(org+".ftbl", "\\"),
    })
    f.write("""
# intermediate jacobian
if (TIMEIT) {
   cat("dfl_dffg : ", date(), "\\n", sep="")
}

dfl_dffg=invAfl %*% p2bfl
if (nb_fgr) {
   dfl_dffg=cBind(dfl_dffg, invAfl%*%g2bfl)
}
dimnames(dfl_dffg)=list(nm_fl, c(nm_ff, nm_fgr))

# prepare mf, md, mc and mg matrices
# such that mf%*%ff+md%*%fl+mc%*%fc+mg%*%fgr gives fallnx
# here ff free fluxes (param), fl are dependent fluxes, fc are constrained
# fluxes and fgr are variable growth fluxes
mf=Matrix(0., nb_fallnx, nb_ff)
dimnames(mf)=list(nm_fallnx, nm_ff)
md=Matrix(0., nb_fallnx, nb_fl)
dimnames(md)=list(nm_fallnx, nm_fl)
mc=Matrix(0., nb_fallnx, nb_fc)
dimnames(mc)=list(nm_fallnx, nm_fc)
mg=Matrix(0., nb_fallnx, nb_fgr)
dimnames(mg)=list(nm_fallnx, nm_fgr)

mf[nm_ff, nm_ff]=diag(1., nb_ff)
md[nm_fl, nm_fl]=diag(1., nb_fl)
mc[nm_fc, nm_fc]=diag(1., nb_fc)
if (nb_fgr) {
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
    #aff("got measures in netan2R_meas", measures);##
    # get scaling factors and their indexes, measure matrices, and measured cumomer value vector
    scale={"label": {}, "mass": {}, "peak": {}}; # for unique scale names
    nrow={"label": {}, "mass": {}, "peak": {}}; # for counting scale names
    o_sc={"label": {}, "mass": {}, "peak": {}}; # for ordered unique scale names
    o_meas=measures.keys(); # ordered measure types
    o_meas.sort()

    if DEBUG:
        pdb.set_trace()

    ir2isc={"label": [], "mass": [], "peak": []}; # for mapping measure rows indexes on scale index
    # we want to use it in python like isc[meas]=ir2isc[meas][ir]
    for meas in o_meas:
        # get unique scaling factors
        # and count rows in each group
        # row["scale"] is "metab;group" (metab name may be fake here)
        for (i,row) in enumerate(measures[meas]["mat"]):
            scale[meas][row["scale"]]=0.
            nrow[meas][row["scale"]]=nrow[meas].get(row["scale"],0.)+1
        # remove groups having only one measure in them
        for (k,n) in list(nrow[meas].iteritems()):
            if n<2:
                del(scale[meas][k])
        # order scaling factor
        o_sc[meas]=scale[meas].keys()
        o_sc[meas].sort()
        # map a measure rows (card:n) on corresponding scaling factor (card:1)
        # if a row has not scale factor it is scaled with factor 1
        # vector having scaling parameters is formed like
        # c(1,param)
        ir2isc[meas]=[-1]*len(measures[meas]["mat"])
        for (i,row) in enumerate(measures[meas]["mat"]):
            if row["scale"] in scale[meas]:
                ir2isc[meas][i]=o_sc[meas].index(row["scale"])

        # measured value vector is in measures[meas]["vec"]
        # measured dev vector is in measures[meas]["dev"]

    if DEBUG:
        pdb.set_trace()

    # create R equivalent structures with indices for scaling
    f.write("""
if (TIMEIT) {
   cat("measure : ", date(), "\n", sep="")
}
# make place for scaling factors
""")
    for meas in o_meas:
        if not o_sc[meas]:
            continue
        f.write("""
# %(meas)s
# initial values for scales are set later
param=c(param,%(sc)s)
nm_par=c(nm_par,c(%(sc_names)s))
""" % {
        "meas": meas,
        "sc": join(", ", (scale[meas][sc] for sc in o_sc[meas])),
        "sc_names": join(", ", o_sc[meas], '"'+meas+';', '"'),
        })

    f.write("""
nb_param=length(param)
nb_sc=nb_param-nb_ff
# indices mapping from scaling to measure matrix row
# c(1,par)[ir2isc] replicates scale parameters
# for corresponding rows of measure matrix
ir2isc=numeric(0)
""")
    base_isc=1
    for meas in o_meas:
        f.write("""
# %(meas)s
ir2isc=c(ir2isc,c(%(ir2isc)s))
""" % {
        "nsc_meas": len(scale[meas]),
        "meas": meas,
        "sc_names": join(", ", o_sc[meas], '"', '"'),
        "ir2isc": join(", ", ((str(ir2isc[meas][ir]+base_isc) if ir2isc[meas][ir]>=0 else -1)
            for ir in xrange(len(ir2isc[meas]))))
        })
        base_isc=base_isc+len(scale[meas])

    f.write("""
# shift indices by nb_ff+1
ir2isc[ir2isc>0]=ir2isc[ir2isc>0]+nb_ff+1
ir2isc[ir2isc<=0]=1

if (noscale) {
#browser()
   # remove scaling params from optimization
   param=head(param, nb_ff)
   nm_par=head(nm_par, nb_ff)
   nb_param=length(param)
   ir2isc[]=1
   nb_sc=0
}
nm_list$par=nm_par
nb_f$nb_sc=nb_sc
""")
    # get the full dict of non zero cumomers involved in measures
    # cumo=metab:icumo where icumo is in [1;2^Clen]
    # or emu=metab:ifrag+Mi
    if DEBUG:
        pdb.set_trace()
    meas_cumos={}
    for meas in o_meas:
        for row in measures[meas]["mat"]:
            metab=row["metab"]
            if emu:
                meas_cumos.update((metab+":"+i, "") for i in row["emuco"].keys())
            else:
                meas_cumos.update((metab+":"+str(icumo), "") for icumo in row["coefs"].keys() if icumo != 0)

    # order involved cumomers (emu)
    o_mcumos=meas_cumos.keys()
    o_mcumos.sort()
    imcumo2i=dict((cumo, i) for (i, cumo) in enumerate(o_mcumos))
    nb_mcumo=len(o_mcumos)
    f.write("""
# make a sparse measurement matrix
# measmat*xr+memaone gives a vector of simulated not-yet-pooled and not-yet-scaled measurements
# all but 0. Coefficients of 0-cumomers (by defenition equal to 1)
# are all regrouped in the memaone.
nm_measmat=c(%(idmeasmat)s)
nm_meas=c(%(idmeas)s)
nm_list$meas=nm_meas
nm_list$measmat=nm_measmat
nb_meas=length(nm_meas)
nb_measmat=length(nm_measmat)
measmat=Matrix(0., nb_measmat, %(ncol)d)
memaone=numeric(nb_measmat)
measvec=c(%(vmeas)s)
measdev=c(%(dev)s)
rownames(measmat)=nm_measmat
names(measvec)=nm_meas
names(measdev)=nm_meas
ipooled=list(ishort=pmatch(nm_meas, nm_measmat))
"""%{
    "nrow": len([measures[meas]["vec"] for meas in measures]),
    "ncol": sum(len(l) for l in (netan["vemu"] if emu else netan["vrcumo"])),
    "idmeasmat": join(", ", (row["id"] for row in
        valval(measures[o]["mat"] for o in o_meas)),
        p='"', s='"'),
    "idmeas": join(", ",  valval(measures[o]["ids"] for o in o_meas), p='"', s='"'),
    "vmeas": join(", ", valval(measures[o]["vec"] for o in o_meas)),
    "dev": join(", ", (sd for sd in valval(measures[o]["dev"]
        for o in o_meas))),
    })

    # get coeffs in the order above with their corresponding indices from total cumomer vector
    base_pooled=0
    for meas in o_meas:
        if not measures[meas]["mat"]:
            continue
        #print("meas="+meas+"; mat="+str(measures[meas]["mat"]));##
        for metpool in measures[meas]["pooled"]:
            f.write("""
# prepare indices of pooled measurements
ipooled[["%(rowid)s"]]=1+%(basep)d+c(%(ind)s)
"""%{
    "rowid": metpool[0],
    "ind": join(", ", metpool[1:]),
    "basep": base_pooled,
}
)
        base_pooled=base_pooled+len(measures[meas]["vec"])
    # preepare measmat indexes and values : ir, ic, val
    f.write("""
ir2isc=ir2isc[ipooled$ishort]
# prepare indexes of dispatching scale params in jacobian
if (nb_sc > 0) {
   ipaire=matrix(0, nrow=0, ncol=2)
   tmp=lapply(1:nb_sc, function(isc) {
      i=which(ir2isc==isc+1+nb_ff)
      ipaire <<- rbind(ipaire, cbind(i, isc))
      return(NULL)
   })
   nb_f$is2m=ipaire
}
# prepare measmat indexes and values : ir, ic, val
ind_mema=matrix(c(
""")
    i=0
    lab2i0=netan["emu2i0" if emu else "rcumo2i0"]
    onelab="0+0" if emu else 0
    fcoef="emuco" if emu else "coefs"
    for meas in o_meas:
        if not measures[meas]["mat"]:
            continue
        for row in measures[meas]["mat"]:
            i+=1
            metab=row["metab"]
            f.write("""%(iricval)s,
"""%{
    "iricval": join(", ", valval((i, lab2i0[metab+":"+str(k)]+1, v)
        for (k, v) in row[fcoef].iteritems() if k != onelab))
    #"i": i,
    #"cumos": join(", ", ((metab+":"+str(k))
    #    for k in row["emuco"].keys()), p='"', s='"') if emu else join(", ", ((metab+":"+str(k) if k else "#x*")
    #    for k in row["coefs"].keys()), p='"', s='"'),
    #"coefs": join(", ", row["emuco"].values()) if emu else join(", ", row["coefs"].values()),
})

    f.write("""
NULL), ncol=3, byrow=T); # close ind_mema creation
measmat[ind_mema[,1:2]]=ind_mema[,3]

memaone=c(%(memaone)s)
names(memaone)=nm_measmat
"""%{
    "memaone": join(", ", (row[fcoef].get(onelab, 0.)
        for meas in o_meas
        for row in measures[meas]["mat"]))
})
    f.write("""
# prepare flux measurements
nb_fmn=%(nb_fmn)d
nb_f$nb_fmn=nb_fmn
nm_fmn=c(%(nm_fmn)s)
nm_list$fmn=nm_fmn

# measured values
fmn=c(%(fmn)s)

# inverse of variance for flux measurements
fmndev=c(%(fmndev)s)

# indices for measured fluxes
# fallnx[ifmn]=>fmn, here fallnx is complete net|xch flux vector
# combining unknown (dependent), free, constrainded and groth fluxes
ifmn=c(%(ifmn)s)
"""%{
    "nb_fmn": len(netan["vflux_meas"]["net"]),
    "nm_fmn": join(", ", trd(("n."+f for f in netan["vflux_meas"]["net"]),
        netan["nx2dfcg"]), '"', '"'),
    "fmn": join(", ", (netan["flux_measured"][fl]["val"]
        for fl in netan["vflux_meas"]["net"])),
    "fmndev": join(", ", (netan["flux_measured"][fl]["dev"]
        for fl in netan["vflux_meas"]["net"])),
    "ifmn": join(", ", (1+netan["vflux_compl"]["net2i"][fl]
        for fl in netan["vflux_meas"]["net"])),
    })
    return {
        "o_meas": o_meas,
        "measures": measures,
        "o_mcumos": o_mcumos,
        "imcumo2i": imcumo2i,
    }

def netan2R_rcumo(netan, org, f):
    # prepare reduced python systems
    #rAb=C13_ftbl.rcumo_sys(netan)
    rAb=netan["rcumo_sys"] if ("rcumo_sys" in netan) else C13_ftbl.rcumo_sys(netan)
    # full matrix is Ab=netan["cumo_sys"]

    # prune ordered cumomer list in reverse order
    # so that deleted item does not change the index
    # for the rest items to prune
    if "vrcumo" not in netan:
        netan["vrcumo"]=copy.deepcopy(netan["vcumo"])
        for i in xrange(len(netan["vrcumo"]),len(rAb["A"]),-1):
            # delete extra weight systems
            del(netan["vrcumo"][i-1])
        for (iw,cumol) in enumerate(netan["vrcumo"]):
            for i in xrange(len(cumol), 0, -1):
                i-=1
                if cumol[i] not in rAb["A"][iw]:
                    #print "prune", i, cumol[i];##
                    del(cumol[i])
    # prepare cumo2i
    # translate cumoname like A:7 to its index in R vector of cumomers
    rcumos=list(valval(netan["vrcumo"]))
    rcumo2i=dict((c,i+1) for (i,c) in enumerate(rcumos))
    emus=list(valval(netan["vemu"]))
    emu2i=dict((c,i+1) for (i,c) in enumerate(emus))
    # composit cumomer vector incu=c(1,xi,xc)
    incu2i_b1=dict((c,i+2) for (i,c) in enumerate(netan["rcumo_input"].keys()+rcumos))
    # write code for reduced cumomer systems
    #netan2Abcumo_f(rAb["A"], rAb["b"],
    #    netan["vrcumo"], netan["input"], ff, netan["fwrv2i"], incu2i_b1, "fwrv2rAbcumo")
    #netan2Abcumo_sp("spAb_old", rAb["A"], rAb["b"],
    #    netan["vrcumo"], netan["input"], f, netan["fwrv2i"], incu2i_b1)
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
    netan["rcumo2i"]=rcumo2i
    netan["emu2i"]=emu2i
    return {
        "rcumo2i": rcumo2i,
        "emu2i": emu2i,
        "rAb": rAb,
    }

def netan2R_cumo(netan, org, f):
    """netan2R_cumo(netan, org, f)->dict
    generate data structures for cumoemr matrices
    """
    # prepare cumo2i
    # translate cumoname like A:7 to its index in R vector of cumomers
    cumos=list(valval(netan["vcumo"]))
    cumo2i=dict((c,i+1) for (i,c) in enumerate(cumos))
    # composite cumomer vector
    incu2i_b1=dict((c,i+2) for (i,c) in enumerate(netan["rcumo_input"].keys()+cumos))

    netan2Abcumo_spr("spAbr_f", netan["cumo_sys"]["A"], netan["cumo_sys"]["b"],
        netan["vcumo"], netan["input"], f, netan["fwrv2i"], incu2i_b1)
    # write R constants and names
    f.write("""
if (TIMEIT) {
   cat("cumo   : ", date(), "\n", sep="")
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
    #dict2kvh(dict((i,t) for (i,t) in enumerate(tfallnx)), "tfallnx.kvh");##
    nb_ineq=len(netan["flux_inequal"]["net"])+len(netan["flux_inequal"]["xch"])
    f.write("""
if (TIMEIT) {
   cat("ineq    : ", date(), "\n", sep="")
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
        for (fl,fa) in ineq[2].iteritems()))))
        for ineq in netan["flux_inequal"]["net"]), p='"n:', s='"'),
   "nm_ix": join(", ", (join("", (ineq[0], ineq[1], join("+",
        ((str(fa)+"*" if fa != 1. else "")+fl
        for (fl,fa) in ineq[2].iteritems()))))
        for ineq in netan["flux_inequal"]["xch"]), p='"x:', s='"'),
})
    
    for (i, ineq) in enumerate(netan["flux_inequal"]["net"]):
        f.write(
"""mi[%(i)s, c(%(f)s)]=%(sign)sc(%(coef)s)
li[%(i)s]=%(sign)s%(li)g
"""%{
    # as R inequality is always ">=" we have to inverse the sign for "<=" in ftbl
    "i": i+1,
    "sign": ("" if ineq[1]=="<=" or ineq[1]=="=<" else "-"),
    "f": join(", ", trd(ineq[2].keys(), f2dfcg_nx_f["net"]), p='"', s='"'),
    "coef": join(", ", ineq[2].values()),
    "li": ineq[0],
    })
    for (i, ineq) in enumerate(netan["flux_inequal"]["xch"]):
        f.write(
"""mi[%(i)s, c(%(f)s)]=%(sign)sc(%(coef)s)
li[%(i)s]=%(sign)s%(li)g
"""%{
    # as R inequality is always ">=" we have to inverse the sign for "<=" in ftbl
    "i": len(netan["flux_inequal"]["net"])+i+1,
    "sign": ("" if ineq[1]=="<=" or ineq[1]=="=<" else "-"),
    "f": join(", ", trd(ineq[2].keys(), f2dfcg_nx_f["xch"]), p='"', s='"'),
    "coef": join(", ", ineq[2].values()),
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
nb_inout=%(nb_inout)d
if (nb_inout > 0) {
   # add cinout low limits on inout net fluxes
   nb_tmp=nrow(mi)
   nm_inout=c(%(nm_inout)s)
   # explicit inequalities take precedence over generic ones
   # so eliminate inout fluxes which are already in inequalities
   nm_itmp=paste("n:.+<=", substring(nm_inout, 5), sep="")
   i=sapply(1:length(nm_itmp), function(k) {
      j=grep(nm_itmp[k], nm_i)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n")
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
   nm_tmp=c(%(nm_notrev)s)
   # explicit inequalities take precedence over generic ones
   # so eliminate notrev fluxes which are already in inequalities
   nm_itmp=paste("n:.+<=", substring(nm_tmp, 5), sep="")
   i=sapply(1:length(nm_itmp), function(k) {
      j=grep(nm_itmp[k], nm_i)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n")
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
      j=grep(nm_itmp[k], nm_i, fix=T)
      #cat(nm_itmp[k], "->", nm_i[j], "\\n")
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
      nm_i=c(nm_i, paste(nm_tmp, ">=", clownr))
      mi[nb_tmp+(1:len_tmp), nm_tmp]=diag(1., len_tmp)
      li=c(li, rep(clownr, len_tmp))
   }
}
nb_fn=nb_fln+nb_ffn
if (cupn != 0 && nb_fn) {
   # add upper limits on [df].net <= cupn for net fluxes
   nb_tmp=nrow(mi)
   mi=rbind(mi, matrix(0, nrow=nb_fn, ncol=nb_fallnx))
   if (nb_fln)
      nm_i=c(nm_i, paste(nm_fln, "<=", cupn, sep=""))
   if (nb_ffn)
      nm_i=c(nm_i, paste(nm_ffn, "<=", cupn, sep=""))
   li=c(li, rep(-cupn, nb_fn))
   mi[nb_tmp+(1:nb_fn),c(nm_fln, nm_ffn)]=diag(-1., nb_fn)
}

"""%{
#   "nb_notrev": len([fli for (fli,t,nxi) in tfallnx
#      if nxi=="n" and t!="c" and fli in netan["notrev"]]),
   "nm_notrev": join(", ", (t+"."+nxi+"."+fli
      for (fli,t,nxi) in tfallnx
      if nxi=="n" and t!="c" and fli in netan["notrev"]),
      p='"', s='"'),
   "nb_inout": len([fli for (fli,t,nxi) in tfallnx
      if nxi=="n" and t!="c" and fli in netan["flux_inout"]]),
   "nm_inout": join(", ", (t+"."+nxi+"."+fli
      for (fli,t,nxi) in tfallnx
      if nxi=="n" and t!="c" and fli in netan["flux_inout"]),
      p='"', s='"'),
})

    f.write("nb_ineq=NROW(li);\n")
    f.write("""
dimnames(mi)=list(nm_i, nm_fallnx)
names(li)=nm_i
# prepare ui matrix and ci vector for optimisation
# ui%*%param-ci>=0
# it is composed of explicite inequalities from ftbl
# and permanent inequalities 0<=xch<=0.999 and scale>=0

# constraints such that ui%*%param[1:nb_ff]-ci>=0
ui=mi%*%(md%*%invAfl%*%p2bfl+mf)
mic=(md%*%invAfl%*%(c2bfl%*%fc+cnst2bfl) + mc%*%fc)
ci=as.numeric(li-mi%*%mic)
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

if (all(ci[zi]<=1.e-10)) {
   ui=ui[!zi,,drop=F]
   ci=ci[!zi]
   nm_i=nm_i[!zi]
} else {
   cat("The following constant inequalities are not satisfied:\\n", file=stderr())
   cat(nm_i[zi][ci[zi]>1.e-14], sep="\\n", file=stderr())
   stop("see above.")
}

# complete ui by zero columns corresponding to scale params
if (nb_sc) {
   ui=cBind(as.matrix(ui), matrix(0., NROW(ui), nb_sc))
   # complete ui by scales >=0
   ui=rBind(ui, cBind(matrix(0, nb_sc, nb_ff), diag(1, nb_sc)))
   ci=c(ci,rep(0., nb_sc))
   nm_i=c(nm_i, paste(nm_par[(nb_ff+1):nb_param], ">=0", sep=""))
   dimnames(ui)[[1]]=nm_i
   names(ci)=nm_i
}

# remove redundant inequalities
nb_i=nrow(ui)
ired=c()
if (nb_i > 0) {
   for (i in 1:(nb_i-1)) {
      nmref=nm_i[i]
      for (j in setdiff((i+1):nb_i, ired)) {
         if (all(ui[j,]==ui[i,]) && ci[i]==ci[j]) {
            # redundancy
            cat("inequality ", nm_i[j], " redundant with ", nmref, " is removed.\n", sep="")
            ired=c(ired, j)
         }
      }
   }
}
if (!is.null(ired)) {
   # remove all ired inequalities
   ui=ui[-ired,,drop=F]
   ci=ci[-ired]
   nm_i=nm_i[-ired]
}
""")
