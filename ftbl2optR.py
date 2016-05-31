#!/usr/bin/env python

r"""
Transform an ftbl to R code which will solve an optimization of flux analysis
problem :math:`\arg \min_{\Theta} S`, where :math:`S=||\mathrm{Predicted}-\mathrm{Observed}||^{2}_{\Sigma}`
and :math:`\Theta` is a vector of parameters to fit: free fluxes (net+xch), scaling parameters and metabolite concentrations
pools.
Two variants of R code can be generated: "s" and "i" for stationary and isotopically
nonstationary labeling.
Predicted vector is obtained from cumomer or emu vector x (calculated from
free fluxes and divided in chunks according to the cumo weight) by
multiplying it by the measurement matrices, weighted by metabolite
pools (in case of pooling) and scale factor (for stationary case only),
boths coming from ftbl file. Observed values vector xo is extracted from ftbl
file for "s" case and from special text file for "i" case.
It is composed of flux, label measurements and metabolite pools.
:math:`\Sigma^2`, covariance diagonal matrices sigma[flux|mass|label|peak|metab.pool]
is orginated from the ftbl file.

usage: ./ftbl2optR.py [opts] organism
where organism is the ftbl informative part of file name
(before .ftbl), e.g. organism.ftbl
after execution a file organism.R will be created.
If it already exists, it will be silently overwritten.
The system Afl*flnx=bfl is created from the ftbl file.

Important python variables:
   * case_i - if True, the case is "i" otherwise it is the "s" case

Collections:
   * netan - (dict) ftbl structured content
   * tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
    collection
   * measures - (dict) exp data
   * rAb - (list) reduced linear systems A*x_cumo=b (a system by weight)
   * scale - unique scale names
   * nrow - counts scale names
   * o_sc - ordered scale names
   * o_meas - ordered measurement types
   
File names (str):
   * n_ftbl (descriptor f_ftbl)
   * n_R (R code) (f)
   * n_fort (fortran code) (ff)

Counts:
   * nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total), nb_ffn, nb_ffx (free fluxes)

Index translators:
   * fwrv2i - flux names to index in R:fwrv
   * cumo2i - cumomer names to index in R:x
   * ir2isc - mapping measurement rows indexes on scale index isc[meas]=ir2isc[meas][ir]

Vector names:
   * cumos (list) - names of R:x
   * o_mcumos - cumomers involved in measurements

Important R variables:

Scalars:
   * nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
   * nb_ffn, nb_ffx, nb_ff (free fluxes),
   * nb_fcn, nb_fcx, nb_fc (constrained fluxes),
   * nb_ineq, nb_param, nb_fmn

Name vectors:
   * nm_cumo, nm_fwrv, nm_fallnx, nm_fln, nm_flx, nm_fl, nm_par,
   * nm_ffn, nm_ffx,
   * nm_fcn, nm_fcx,
   * nm_mcumo, nm_fmn

Numeric vectors:
   * fwrv - all fluxes (fwd+rev)
   * x - all cumomers (weight1+weight2+...)
   * param - free flux net, free flux xch, scale label, scale mass, scale peak, metabolite concentrations
   * fcn, fcx, fc - constrained fluxes
   * bp - helps to construct the rhs of flux system
   * xi -cumomer input vector
   * fallnx - complete flux vector (constr+net+xch)
   * bc - helps to construct fallnx
   * li - inequality vector (mi%*%fallnx>=li)
   * ir2isc - measur row to scale vector replicator
   * ci - inequalities for param use (ui%*%param-ci>=0)
   * measvec - measurement vector
   * fmn

Matrices:
   * Afl, qrAfl, invAfl,
   * p2bfl - helps to construct the rhs of flux system
   * mf, md - help to construct fallnx
   * mi - inequality matrix (ftbl content)
   * ui - inequality matrix (ready for param use)
   * measmat - for measmat*x+memaone=vec of simulated not-yet-scaled measurements

Functions:
   * lab_sim - translate param to flux and cumomer vector (initial approximation)
   * cumo_cost - cost function (khi2)
   * cumo_gradj - implicit derivative gradient
"""

# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used for sparse matrices.
# 2014-04-14 sokol:adapted for both "s" and "i" cases

if __name__ == "__main__":
    import sys
    import os
    import stat
    import time
    import copy
    import getopt
    #import pdb

    me=os.path.realpath(sys.argv[0])
    dirx=os.path.dirname(me)
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
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--emu] [--clownr] [--ropts ROPTS] network_name[.ftbl]\n")

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "emu", "clownr", "ropts=", "case_i"])
    except getopt.GetoptError, err:
        #pass
        sys.stderr.write(str(err)+"\n")
        usage()
        sys.exit(1)

    fullsys=False
    emu=False
    clownr=False
    ropts=['""']
    case_i=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--fullsys":
            fullsys=True
        elif o=="--emu":
            emu=True
        elif o=="--clownr":
            clownr=True
        elif o=="--ropts":
            ropts=a.split("; ") if len(a) else ['""']
        elif o=="--case_i":
            case_i=True
        else:
            #assert False, "unhandled option"
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            #pass
            raise
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
    import ftbl2code
    ftbl2code.case_i=case_i
    C13_ftbl.clownr=clownr
    #org="ex3"
    #org="PPP_exact"

    n_ftbl=fullorg+".ftbl"
    n_R=fullorg+".R"
    #n_fort=fullorg+".f"
    try:
        os.chmod(n_R, stat.S_IWRITE)
    except:
        pass
    f=open(n_R, "w")

    # parse ftbl
    ftbl=C13_ftbl.ftbl_parse(n_ftbl)

    # analyse network
    # reload(C13_ftbl)

    netan=C13_ftbl.ftbl_netan(ftbl, emu, fullsys)
    # prepare rcumo system
    rAb=C13_ftbl.rcumo_sys(netan, emu)

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys, emu, ropts)

    ropts_s="\n\t\t".join(ropts)
    if ropts_s and ropts_s[0]=='"':
        ropts_s=ropts_s[1:-1]
    f.write("""
#browser()

# extend param vector by free pools
if (nb_poolf > 0) {
   param=c(param, poolf)
   nm_par=c(nm_par, nm_poolf)
   nb_param=length(param)
}
nm_list$par=nm_par

#browser()
if (nb_poolf > 0) {
   # extend inequalities ui, ci by uip, cip
   nb_row=nrow(ui)
   nb_col=ncol(ui)
   ui=cbind(ui, matrix (0., nrow=nb_row, ncol=nb_poolf)) # add 0-columns
   ui=rbind(ui, cbind(matrix(0., nrow(uip), ncol=nb_col), uip))
   ci=c(ci, cip)
   
   # extend inequalities ui, ci by cupp>= poolf >= clowp
   # but exclude metabolites that are individually set in the uip (FTBL)
   met_low=met_up=c()
   if (nrow(uip) > 0) {
      # number of non zero entries per row in uip
      i_nz=rowSums(abs(uip) != 0.)
      # number of positive coeffs for alone metabs (i.e. low limit is set)
      i_pos=colSums(uip[i_nz > 0,,drop=F] > 0)
      met_low=nm_poolf[i_pos > 0]
      # number of negative coeffs for alone metabs (i.e. upper limit is set)
      i_neg=colSums(uip[i_nz > 0,,drop=F] < 0)
      met_up=nm_poolf[i_neg > 0]
   }

   # add low limit
   nb_add=nb_poolf-length(met_low)
   if (nb_add > 0) {
      nm_add=nm_poolf[!nm_poolf %in% met_low]
      ui_add=matrix(0., nrow=nb_add, ncol=ncol(ui))
      ui_add[,nb_col+pmatch(nm_add, nm_poolf)]=diag(1., nb_add)
      rownames(ui_add)=paste(nm_add, ">=", clowp, sep="")
      if (nrow(ui)) {
         ui=rbind(ui, ui_add)
      } else {
         ui=ui_add
      }
      ci=c(ci, rep(clowp, nb_add))
   }

   # add upper limit
   nb_add=nb_poolf-length(met_up)
   if (nb_add > 0) {
      nm_add=nm_poolf[!nm_poolf %in% met_up]
      ui_add=matrix(0., nrow=nb_add, ncol=ncol(ui))
      ui_add[,nb_col+pmatch(nm_add, nm_poolf)]=diag(-1., nb_add)
      rownames(ui_add)=paste(nm_add, "<=", cupp, sep="")
      ui=rbind(ui, ui_add)
      ci=c(ci, rep(-cupp, nb_add))
   }

   nm_i=names(ci)=rownames(ui)
   colnames(ui)=nm_par
}
# extend the matrix of metabolite equalities
ep=cbind(matrix(0., nrow(ep), nb_param-nb_poolf), ep)
colnames(ep)=nm_par
""")
    f.write("""
# prepare metabolite pools measurements
nb_poolm=%(nb_poolm)d
nb_f$nb_poolm=nb_poolm
nm_poolm=c(%(nm_poolm)s)
nm_list$poolm=nm_poolm

# measured values
vecpoolm=c(%(v_poolm)s)
names(vecpoolm)=nm_poolm

# inverse of variance for pool measurements
poolmdev=c(%(poolmdev)s)

# simulated metabolite measurements are calculated as
# measmatpool*poolall=>poolm
measmatpool=matrix(0., nrow=nb_poolm, ncol=length(poolall))
dimnames(measmatpool)=list(nm_poolm, nm_poolall)
i=matrix(1+c(%(imeasmatpool)s), ncol=2, byrow=T)
measmatpool[i]=1.

"""%{
    "nb_poolm": len(netan["metab_measured"]),
    "nm_poolm": join(", ", netan["metab_measured"].keys(), '"pm:', '"'),
    "v_poolm": join(", ", (item["val"] for item in netan["metab_measured"].values())).replace("nan", "NA"),
    "poolmdev": join(", ", (item["dev"] for item in netan["metab_measured"].values())),
    "imeasmatpool": join(", ", valval((ir,netan["vpool"]["all2i"][m]) for (ir, ml) in enumerate(netan["metab_measured"].keys()) for m in ml.split("+"))),
})
    if case_i:
        f.write("""
## variables for isotopomer kinetics
tstart=0.
tmax=%(tmax)f
if (tmax < 0) {
   stop_mes(sprintf("The parameter tmax must not be negative (tmax=%%g)", tmax), fcerr)
}
dt=%(dt)f
if (dt <= 0) {
   stop_mes(sprintf("The parameter dt must be positive (dt=%%g)", dt), fcerr)
}

# read measvecti from a file specified in ftbl
flabcin="%(flabcin)s"
if (nchar(flabcin)) {
   if (substr(flabcin, 1, 1) == "/")
      flabcin=file.path(flabcin)
   else
      flabcin=file.path(dirw, flabcin)
   measvecti=as.matrix(read.table(flabcin, header=T, row.names=1, sep="\t", check=F, comment="#"))
   nm_row=rownames(measvecti)
   # put in the same row order as simulated measurements
   # check if nm_meas are all in rownames
   if (all(nm_meas %%in%% nm_row)) {
      measvecti=measvecti[nm_meas,,drop=F]
   } else {
      # try to strip row number from measure id
      nm_strip=sapply(strsplit(nm_meas, ":"), function(v) {
         v[length(v)]=""
         paste(v, sep="", collapse=":")
      })
      im=pmatch(nm_strip, nm_row)
      ina=is.na(im)
      if (any(ina)) {
         mes=paste("Cannot match the following measurement(s) in the file '", flabcin, "':\\n", paste(nm_meas[ina], sep="", collapse="\\n"), "\\n", sep="", collapse="")
         stop_mes(mes, file=fcerr)
      }
      measvecti=measvecti[im,,drop=F]
      #stopifnot(all(!is.na(measvecti)))
      stopifnot(typeof(measvecti)=="double" || all(is.na(measvecti)))
   }
   ti=as.double(colnames(measvecti))
   if (any(is.na(ti))) {
      mes=sprintf("All time moments (in column names) could not be converted to real numbers in the file '%%s'\\nConverted times:\\n%%s", flabcin, join("\\n", ti))
      stop_mes(mes, file=fcerr)
   }
   if (length(ti) < 1L) {
      mes=sprintf("No column found in the file '%%s'", flabcin)
      stop_mes(mes, file=fcerr)
   }
   if (!all(diff(ti) > 0.)) {
      mes=sprintf("Time moments (in column names) are not monotonously increasing in the file '%%s'", flabcin)
      stop_mes(mes, file=fcerr)
   }
   if (ti[1L] < 0.) {
      mes=sprintf("The first time moment cannot be negative in the file '%%s'", flabcin)
      stop_mes(mes, file=fcerr)
   }
   if (ti[1L] != 0.) {
      ti=c(tstart, ti)
   }
   i=which(ti<=tmax)
   ti=ti[i]
   measvecti=measvecti[,i[-1]-1,drop=F]
} else {
   measvecti=NULL
   ti=seq(tstart, tmax, by=dt)
   if (optimize) {
      cat("Warning: a fitting is requested but no file with label data is provided by 'file_labcin' option in the ftbl file.
The fitting is ignored as if '--noopt' option were asked.\\n", file=fcerr)
      optimize=F
   }
}
nb_ti=length(ti)
if (nb_ti < 2L) {
   mes=sprintf("After filtering by tmax, only %%d time moments are kept. It is not sufficient.", nb_ti)
   stop_mes(mes, file=fcerr)
}

# divide the first time interval by n1 geometric intervals
n1=1
tmp=cumsum(2**seq_len(n1))
tifull=c(ti[1L], ti[2L]*tmp/tmp[n1], ti[-(1L:2L)])

# divide each time interval by nsubdiv_dt
dt=diff(tifull)
nsubdiv_dt=max(1L, as.integer(%(nsubdiv_dt)s))
dt=rep(dt/nsubdiv_dt, each=nsubdiv_dt)
tifull=c(tifull[1L], cumsum(dt))
nb_tifu=length(tifull)

tifull2=c(tifull[1L], tifull[1L]+cumsum(rep(diff(tifull)/2., each=2L)))

if (length(ijpwef)) {
   # vector index for many time points
   ijpwef=cbind(ijpwef[,1L], rep(seq_len(nb_ti-1L), each=nrow(ijpwef)), ijpwef[,2L])
   dp_ones=matrix(aperm(array(dp_ones, c(dim(dp_ones), nb_ti-1L)), c(1L, 3L, 2L)), ncol=nb_poolf)
}
"""%{
    "dt": netan["opt"]["dt"],
    "tmax": netan["opt"]["tmax"],
    "flabcin": netan["opt"].get("file_labcin", ""),
    "nsubdiv_dt": netan["opt"].get("nsubdiv_dt", "1"),
})

        f.write("""
# prepare mapping of metab pools on cumomers
nminvm=nm_poolall[matrix(unlist(strsplit(nm_rcumo, ":")), ncol=2, byrow=T)[,1L]]
nb_f$ip2ircumo=match(nminvm, nm_poolall)
nb_f$ipf2ircumo=nb_f$ipf2ircumo2=list()
for (iw in seq_len(nb_w)) {
   ix=seq_len(nb_rcumos[iw])
   ipf2ircumo=ipf2ircumo2=match(nminvm[nbc_cumos[iw]+ix], nm_poolf, nomatch=0L)
   dims=c(1L, nb_rcumos[iw], ifelse(emu, iw, 1L), nb_tifu-1L)
   dims2=c(1L, nb_rcumos[iw], ifelse(emu, iw, 1L), (nb_tifu-1L)*2)
   i=as.matrix(ipf2ircumo)
   i2=as.matrix(ipf2ircumo2)
   for (id in 2L:length(dims)) {
      cstr=sprintf("cbind(%srep(seq_len(dims[id]), each=prod(dims[seq_len(id-1L)])))", paste("i[, ", seq_len(id-1L), "], ", sep="", collapse=""))
      i=eval(parse(text=cstr))
   }
   for (id in 2L:length(dims2)) {
      cstr=sprintf("cbind(%srep(seq_len(dims2[id]), each=prod(dims2[seq_len(id-1L)])))", paste("i2[, ", seq_len(id-1L), "], ", sep="", collapse=""))
      i2=eval(parse(text=cstr))
   }
   colnames(i)=c("ipoolf", "ic", "iw", "iti")
   colnames(i2)=c("ipoolf", "ic", "iw", "iti")
   i=i[i[,1L]!=0L,,drop=F]
   i2=i2[i2[,1L]!=0L,,drop=F]
   # put the poolf column last
   nb_f$ipf2ircumo[[iw]]=cbind(i[,-1L,drop=F], i[,1L,drop=F])
   nb_f$ipf2ircumo2[[iw]]=cbind(i2[,-1L,drop=F], i2[,1L,drop=F])
}

# label state at t=0 (by default=0 but later it should be able to be specified by user)
x0=NULL
nb_ti=length(ti)
""")
    f.write("""
# gather all measurement information
measurements=list(
   vec=list(labeled=measvec, flux=fmn, pool=vecpoolm, kin=if (case_i) measvecti else NULL),
   dev=list(labeled=measdev, flux=fmndev, pool=poolmdev, kin=if (case_i) rep(measdev, nb_ti-1) else NULL),
   mat=list(labeled=measmat, flux=ifmn, pool=measmatpool),
   one=list(labeled=memaone)
)
nm_resid=c(if (case_i) outer(nm_meas, ti[-1L], paste, sep=", t=") else nm_meas, nm_fmn, nm_poolm)
nm_list$resid=nm_resid

if (TIMEIT) {
   cat("preopt  : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
#browser()
names(param)=nm_par
# prepare series of starting points
if (nchar(fseries) > 0) {
   pstart=as.matrix(read.table(fseries, header=T, row.n=1, sep="\\t"))
   # skip parameters (rows) who's name is not in nm_par
   i=rownames(pstart) %in% nm_par
   if (!any(i)) {
      stop_mes("Option --fseries is used but no free parameter with known name is found.\\n")
   }
   pstart=pstart[i,,drop=F]
   cat("Using starting values form '", fseries, "' for the following free parameters:\\n", paste(rownames(pstart), collapse="\\n"), "\\n", sep="", file=fclog)
   nseries=ncol(pstart)
   if (initrand) {
      # fill the rest of rows with random values
      i=nm_par %in% rownames(pstart)
      n=sum(!i)
      pstart=rbind(pstart, matrix(runif(n*nseries), n, nseries))
      rownames(pstart)=c(rownames(pstart)[seq_len(nb_param-n)], nm_par[!i])
   }
   if (nchar(iseries) > 0) {
      iseries=unique(as.integer(eval(parse(t="c("%s+%iseries%s+%")"))))
      iseries=iseries[iseries<=nseries]
      # subsample
      pstart=pstart[,iseries, drop=F]
      nseries=ncol(pstart)
   } else {
      iseries=seq_len(nseries)
   }
} else if (nchar(iseries) > 0) {
   # first construct pstart then if needed fill it with random values
   # and only then subsample
   iseries=unique(as.integer(eval(parse(t="c("%s+%iseries%s+%")"))))
   nseries=max(iseries)
   pstart=matrix(rep(param, nseries), nrow=nb_param, ncol=nseries)
   dimnames(pstart)=list(nm_par, paste("V", seq_len(nseries), sep=""))
   if (initrand) {
      # fill pstart with random values
      pstart[]=runif(length(pstart))
   }
   # subsample
   pstart=pstart[,iseries, drop=F]
   nseries=ncol(pstart)
} else {
   iseries=1L
   pstart=as.matrix(param)
   nseries=1L
   if (initrand) {
      # fill pstart with random values
      pstart[]=runif(length(pstart))
   }
}
nm_pseries=rownames(pstart)

if (is.null(nseries) || nseries==0) {
   stop_mes(sprintf("No starting values in the series file '%s' or --iseries is empty.", fseries),  file=fcerr)
}

pres=matrix(NA, nb_param, nseries)
rownames(pres)=nm_par
colnames(pres)=colnames(pstart)
costres=rep.int(NA, nseries)

# prepare flux index conversion
ifwrv=1:nb_fwrv
names(ifwrv)=nm_fwrv
ifl_in_fw=if (nb_fln) ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")] else integer(0)
iff_in_fw=if (nb_ff > 0) ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")] else integer(0)
ifg_in_fw=if (nb_fgr > 0) ifwrv[paste("fwd", substring(nm_fgr, 4), sep="")] else integer(0)

# index couples for jacobian df_dfl, df_dffd
cfw_fl=crv_fl=cbind(ifl_in_fw, seq_len(nb_fl))
cfw_ff=crv_ff=cbind(iff_in_fw, seq_len(nb_ff))
cfw_fg=crv_fg=cbind(ifg_in_fw, nb_ff+seq_len(nb_fgr))
crv_fl[,1L]=(nb_fwrv/2)+crv_fl[,1L]
crv_ff[,1L]=(nb_fwrv/2)+crv_ff[,1L]
crv_fg[,1L]=(nb_fwrv/2)+crv_fg[,1L]

# store it in nb_f
nb_f=append(nb_f, list(cfw_fl=cfw_fl, crv_fl=crv_fl, cfw_ff=cfw_ff,
   crv_ff=crv_ff, cfw_fg=cfw_fg, crv_fg=crv_fg))

nbc_x=c(0, cumsum(nb_x))
nb_f$nbc_x=nbc_x

# fixed part of jacobian (unreduced by SD)
# measured fluxes
dufm_dp=cbind(dufm_dff(nb_f, nm_list), matrix(0, nrow=nb_fmn, ncol=nb_sc+nb_poolf))
dimnames(dufm_dp)=list(nm_fmn, nm_par)

# measured pools
dupm_dp=matrix(0., nb_poolm, nb_ff+nb_sc)
if (nb_poolf > 0L) {
   dupm_dp=cbind(dupm_dp, measurements$mat$pool[,nm_list$poolf, drop=F])
}
dimnames(dupm_dp)=list(rownames(measurements$mat$pool), nm_par)

#browser()
# prepare argument list for passing to label simulating functions
nm_labargs=c("jx_f", "nb_f", "nm_list", "nb_x", "invAfl", "p2bfl", "g2bfl", "bp", "fc", "xi", "spa", "emu", "pool", "measurements", "ipooled", "ir2isc",  "nb_w", "nbc_x", "measmat", "memaone", "dufm_dp", "dupm_dp", "pwe", "ipwe", "ip2ipwe", "pool_factor", "ijpwef", "ipf_in_ppw", "meas2sum", "dp_ones", "clen", "dirx", "use_magma", "case_i")
""")
    if case_i:
        f.write("""nm_labargs=c(nm_labargs, "ti", "tifull", "tifull2", "x0", "time_order")
""")
    f.write("""
if (TIMEIT) {
   cat("labargs : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
labargs=new.env()
tmp=lapply(nm_labargs, function(nm) assign(nm, get(nm), labargs))
#for (nm in nm_labargs) {
#   labargs[[nm]]=get(nm)
#}
labargs[["nm"]]=labargs[["nm_list"]]
labargs[["spAb"]]=labargs[["spa"]]

# formated output in kvh file
fkvh_saved=file.path(dirw, sprintf("%s_res.kvh", baseshort))
""")
    f.write(r"""
retcode=numeric(nseries)
for (irun in seq_len(nseries)) {
   if (TIMEIT) {
      cat(sprintf("run %4d: %s cpu=%g\n", irun, format(Sys.time()), proc.time()[1], "\\n", sep=""), file=fclog)
   }
   param[nm_pseries]=pstart[nm_pseries, irun]
#browser()
   # prepare kvh file name
   if (nseries > 1) {
      runsuf="." %s+% colnames(pstart)[irun]
   } else {
      runsuf=""
   }
   if (length(nseries) > 0) {
      cat("Starting point", runsuf, "\n", sep="", file=fclog)
   }
   fkvh=file(substring(fkvh_saved, 1, nchar(fkvh_saved)-4) %s+% runsuf %s+% ".kvh", "w");

   # remove zc inequalities from previous runs
   izc=grep("^zc ", nm_i)
   if (length(izc)) {
      ui=ui[-izc,,drop=F]
      ci=ci[-izc]
      nm_i=rownames(ui)
   }
   # check if initial approximation is feasible
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   if (any(ineq <= -1.e-10)) {
      cat("The following ", sum(ineq<= -1.e-10), " ineqalities are not respected at starting point", runsuf, ":\n", sep="", file=fclog)
      i=ineq[ineq<= -1.e-10]
      cat(paste(names(i), i, sep="\t", collapse="\n"), "\n", sep="", file=fclog)
      # put them inside
      capture.output(pinside <- put_inside(param, ui, ci), file=fclog)
      if (any(is.na(pinside))) {
         if (!is.null(attr(pinside, "err")) && attr(pinside, "err")!=0) {
            # fatal error occured
            cat("put_inside", runsuf, ": ", attr(pinside, "mes"), "\n",
               file=fcerr, sep="")
            close(fkvh)
            retcode[irun]=attr(pinside, "err")
            next;
         }
      } else if (!is.null(attr(pinside, "err")) && attr(pinside, "err")==0) {
         # non fatal problem
         cat(paste("put_inside: ", attr(pinside, "mes"), collapse=""), "\n", file=fcerr)
      }
      param[]=pinside
   }

   # prepare zero crossing strategy
   # inequalities to keep sens of net flux on first call to opt_wrapper()
   # if active they are removed on the second call to opt_wrapper()
   # and finaly all zc constraints are relaxed on the last call to opt_wrapper()
   fallnx=param2fl(param, labargs)$fallnx
   mi_zc=NULL
   li_zc=NULL
   if (zerocross && length(grep("^[df]\\.n\\.", nm_fallnx))>0) {
      if (TIMEIT) {
         cat("zc ineq : ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
      }
#browser()
      # prepare fluxes that are already in inequalities in alone mode
      ige=names(which(apply(mi, 1L, function(v) diff(range(v))==1 && sum(v)==1) & li>=0))
      ige=nm_dfn[unique(c(
         sub("^n:.+<=(.+)$", "\\1", grep("^n:.+<=.+$", ige, v=T)),
         sub("^[df]\\.n\\.(.+)>=.+$", "\\1", grep("^[df]\\.n\\..+>=.+$", ige, v=T)),
         sub("^inout [df]\\.n\\.(.+)>=.+$", "\\1", grep("^inout [df]\\.n\\..+>=.+$", ige, v=T))
      ))]
      ile=which(apply(mi, 1L, function(v) diff(range(v))==1 && sum(v)==-1)&li>=0)
      ile=nm_dfn[unique(c(
         sub("^n:.+<=(.+)$", "\\1", grep("^n:.+<=.+$", ile, v=T)),
         sub("^[df]\\.n\\.(.+)>=.+$", "\\1", grep("^[df]\\.n\\..+>=.+$", ile, v=T)),
         sub("^inout [df]\\.n\\.(.+)>=.+$", "\\1", grep("^inout [df]\\.n\\..+>=.+$", ile, v=T))
      ))]
      # add lower limits on [df].net >= zc for positive net fluxes
      # and upper limits on [df].net <= -zc for negative net fluxes
      nm_izc=c()
      ipos=setdiff(names(which(fallnx[grep("^[df]\\.n\\.", nm_fallnx)]>=0.)), ige)
      ineg=setdiff(names(which(fallnx[grep("^[df]\\.n\\.", nm_fallnx)]<0.)), ile)
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
      li_zc=rep(zc, length(nm_izc)) # that's ok for both pos and neg constraints
      ui_zc=cbind(mi_zc%*%(md%*%invAfl%*%p2bfl+mf),
         matrix(0., nrow=nrow(mi_zc), ncol=nb_sc))
      if (nb_fgr > 0) {
         ui_zc=cbind(ui_zc, mi_zc%*%((md%*%invAfl%*%g2bfl)+mg*nb_f$mu))
      } else if (nb_poolf > 0) {
         ui_zc=cbind(ui_zc, matrix(0., nrow=nrow(mi_zc), ncol=nb_poolf))
      }
      ci_zc=li_zc-mi_zc%*%mic
      # remove constant inequalities
      if (ncol(ui_zc)) {
         zi=apply(ui_zc,1,function(v){return(max(abs(v))<=1.e-14)})
      } else {
         # remove all flux inequalities as there is no free params
         zi=rep(TRUE, nrow(ui_zc))
      }

      inotsat=ci_zc[zi]>1.e-10
      if (any(inotsat)) {
         cat("Warning: The following constant inequalities are not satisfied:\n", file=fcerr)
         cat(nm_izc[zi][inotsat], sep="\n", file=fcerr)
      }
      ui_zc=ui_zc[!zi,,drop=F]
      ci_zc=ci_zc[!zi]
      nm_izc=nm_izc[!zi]
      mi_zc=mi_zc[!zi,,drop=F]

      # remove redundant/contradictory inequalities
      nb_zc=nrow(ui_zc)
      nb_i=nrow(ui)
      ired=c()
      tui=t(ui)
      uzcd=sapply(seq_len(nb_zc), function(i) apply(abs(tui-ui_zc[i,]), 2L, max))
      uzcs=sapply(seq_len(nb_zc), function(i) apply(abs(tui+ui_zc[i,]), 2L, max))
      czcd=abs(outer(abs(ci), abs(ci_zc), "-"))
      ired=which(apply((uzcd < 1.e-10 | uzcs < 1.e-10) & czcd <= 1.e-2, 2, any))
      
      if (length(ired) > 0L) {
         # remove all ired inequalities
         cat("The following ", length(ired), " zerocross inequalities are redundant and are removed:\n", paste(nm_izc[ired], collapse="\n"), "\n", sep="", file=fclog)
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
      }
      rm(ui_zc, ci_zc, uzcd, uzcs, czcd)
   }
   rres=NULL
""")
    if case_i:
        f.write("""
   if (nb_sc && !is.null(measvecti)) {
      if (TIMEIT) {
         cat("res esti: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
      }
      capture.output(rres <- lab_resid(param, cjac=F, labargs), file=fclog)
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         close(fkvh)
         retcode[irun]=rres$err
         next
      }
      if (sum(is.infinite(rres$res))) {
         cat("Infinite values appeared in residual vector", file=fcerr)
         retcode[irun]=1
         close(fkvh)
         next
      }
      # set initial scale values to sum(measvec*simlab/dev**2)/sum(simlab**2/dev**2)
      # for corresponding measurements
      # unscaled simulated measurements (usm) [imeas, itime]
      #browser()
      inna=which(!is.na(measvecti))
      simlab=jx_f$usm
      measinvvar=1./measurements$dev$labeled**2
      ms=(measvecti*simlab*measinvvar)[inna]
      ss=(simlab*simlab*measinvvar)[inna]
      for (i in nb_ff+1:nb_sc) {
         im=outer(ir2isc==(i+1), rep(T, nb_ti), "&")[inna]
         param[i]=sum(ms[im])/sum(ss[im])
      }
      # prepare is2mti which is extension of is2m ti times
      nb_f$is2mti=array(0., dim=c(dim(nb_f$is2m), nb_ti))
      nb_f$is2mti[]=nb_f$is2m
      nb_f$is2mti[,1L,]=nb_f$is2m[,1L]+rep((0:(nb_ti-1))*nb_meas, each=nrow(nb_f$is2m))
#browser()
   } else if (nb_sc > 0) {
      # we dont have measurements yet, just set all scalings to 1.
      param[nb_ff+1:nb_sc]=1.
   }
""")
    else:
        f.write("""
   # set initial scale values to sum(measvec*simlab/dev**2)/sum(simlab**2/dev**2)
   # for corresponding measurements
   if (nb_sc > 0) {
      if (optimize) {
         if (TIMEIT) {
            cat("res esti: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         capture.output(rres <- lab_resid(param, cjac=F, labargs), file=fclog)
         if (!is.null(rres$err) && rres$err) {
            cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
            close(fkvh)
            retcode[irun]=rres$err
            next
         }
         if (sum(is.infinite(rres$res))) {
            cat("Infinite values appeared in residual vector (at init scale values)", file=fcerr)
            retcode[irun]=1
            close(fkvh)
            next
         }
         simlab=jx_f$usimcumom
         measinvvar=1./measurements$dev$labeled**2
         ms=measvec*simlab*measinvvar
         ss=simlab*simlab*measinvvar
         # get only valid measurements
         iva=!is.na(ms)
         for (i in nb_ff+1:nb_sc) {
            im=(ir2isc==(i+1)) & iva
            if (sum(im) < 2) {
               mes=sprintf("scaling: no sufficient valid data for scaling factor '%s'\\n", nm_par[i])
               stop_mes(mes, fcerr)
            }
            param[i]=sum(ms[im])/sum(ss[im])
         }
      } else {
         # if no optimization, set all scaling params to 1.
         param[nb_ff+1:nb_sc]=1.
      }
   }
""")
    f.write("""
   # see if there are any active inequalities at starting point
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   if (any(abs(ineq)<=1.e-10)) {
      cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point", runsuf, ":\\n",
         paste(names(ineq[abs(ineq)<=1.e-10]), collapse="\\n"), "\\n", sep="", file=fclog)
   }
""")

    f.write("""
   if (TIMEIT) {
      cat("kvh init: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
   }
""")
    # main part: call optimization
    f.write("""
   cat("influx\\n", file=fkvh)
   cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="")
   cat("\\tlabeling\\t", if (case_i) "instationary" else "stationary", "\\n", file=fkvh, sep="")
   # save options of command line
   cat("\\truntime options\\n", file=fkvh)
   cat("\\t\\t%s\\n", file=fkvh)
   """%ropts_s)
    f.write("""
   obj2kvh(R.Version(), "R.Version", fkvh, indent=1)
   cat("\\tR command line\\n", file=fkvh)
   obj2kvh(opts, "opts", fkvh, indent=2)
   cat("\\t\\texecution date\t", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fkvh)

   # resume system sizes
   obj2kvh(nb_sys, "system sizes", fkvh)

   # save initial param
   cat("starting point\\n", file=fkvh)
   names(param)=nm_par
   obj2kvh(param, "starting free parameters", fkvh, indent=1)
#browser()
   if (is.null(rres)) {
      capture.output(rres <- lab_resid(param, cjac=F, labargs), file=fclog)
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         close(fkvh)
         retcode[irun]=rres$err
         next
      }
      if (sum(is.infinite(rres$res))) {
         cat("Infinite values appeared in residual vector (at starting point)", file=fcerr)
         retcode[irun]=1
         close(fkvh)
         next
      }
   }
   rcost=if (is.null(rres$res)) NA else sum(crossprod(rres$res))
   obj2kvh(rcost, "starting cost value", fkvh, indent=1)

   obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1)
   fg=numeric(nb_f$nb_fgr)
   names(fg)=nm_list$fgr
   if (nb_f$nb_fgr > 0) {
      fg[paste("g.n.", substring(nm_list$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm_list$poolf]
   }
   btmp=as.numeric(p2bfl%*%param[seq_len(nb_f$nb_ff)]+bp+g2bfl%*%fg)
   names(btmp)=dimnames(Afl)[[1]]
   obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1)

   #cat("mass vector:\\n", file=fclog)
   #print_mass(x)

   names(param)=nm_par
""")
    f.write("""
   control_ftbl=list(%(ctrl_ftbl)s)
"""%{
    "ctrl_ftbl": join(", ", (k[8:]+"="+str(v) for (k,v) in netan["opt"].iteritems() if k.startswith("optctrl_"))),
})
    f.write("""
#browser()
   if (optimize && nb_ff+nb_poolf > 0L) {
      if (!(least_norm || method!="nlsic")) {
         # check if at starting position all fluxes can be resolved
         if (TIMEIT) {
            cat("check ja: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         rres=lab_resid(param, cjac=T, labargs)
         if (sum(is.infinite(rres$res))) {
            cat("Infinite values appeared in residual vector (at identifiability check)", file=fcerr)
            retcode[irun]=1
            close(fkvh)
            next
         }
         if (sum(is.infinite(rres$jacobian))) {
            cat("Infinite values appeared in Jacobian (at identifiability check)", file=fcerr)
            retcode[irun]=1
            close(fkvh)
            next
         }
         qrj=qr(jx_f$dr_dff, LAPACK=T)
         d=diag(qrj$qr)
         qrj$rank=sum(abs(d)>abs(d[1])*1.e-10)
         if (is.na(qrj$rank)) {
            cat("Rank of starting jacobian could not be estimated.", file=fcerr)
            retcode[irun]=1
            close(fkvh)
            next
         }
         if (qrj$rank) {
            nm_uns=nm_ff[qrj$pivot[-(1:qrj$rank)]]
         } else {
            nm_uns=nm_ff
         }
         if (qrj$rank < nb_ff) {
            # Too bad. The jacobian of free fluxes is not of full rank.
            dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx)
            fname="dbg_dr_dff_singular" %s+% runsuf %s+% ".csv"
            cat(sprintf("Provided measurements (labeling and fluxes) are not sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:\\n%s\\nJacobian dr_dff is written in the result kvh file.\\n",
               paste(nm_uns, sep=", ", collapse=", ")),
               file=fcerr)
            obj2kvh(jx_f$dr_dff, "Jacobian dr_dff", fkvh, indent=0)
            close(fkvh)
            retcode[irun]=1
            next
         }
      }
      if (TIMEIT) {
         cat("optim   : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
      }
      # pass control to the chosen optimization method
      if (time_order=="1,2")
         labargs$time_order="1" # start with order 1, later continue with 2
      capture.output(res <- opt_wrapper(param, measurements, jx_f), file=fclog)
      if ((!is.null(res$err) && res$err) || is.null(res$par)) {
         cat("first optimization pass", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
         res$par=rep(NA, length(param))
         res$cost=NA
      } else if (!is.null(res$mes) && nchar(res$mes)) {
         cat("first optimization pass", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
      }
      if (any(is.na(res$par))) {
         res$retres$jx_f=NULL # to avoid writing of huge data
         obj2kvh(res, "failed first pass optimization process information", fkvh)
         cat("Optimization failed", runsuf, "\\n", file=fcerr, sep="")
         close(fkvh)
         retcode[irun]=max(res$err, 1)
         next
      }
      param=res$par
#browser()
      if (zerocross && !is.null(mi_zc)) {
         if (TIMEIT) {
            cat("secondzc: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         # inverse active "zc" inequalities
         nm_inv=names(which((ui%*%res$par-ci)[,1]<=1.e-10))
         i=grep("^zc ", nm_inv, v=T)
         if (length(i) > 0) {
            i=str2ind(i, nm_i)
            cat("The following inequalities are active after first pass
of zero crossing strategy and will be inverted", runsuf, ":\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="", file=fclog)
            ipos=grep(">=", nm_i[i], v=T)
            ineg=grep("<=", nm_i[i], v=T)
            ui[i,]=-ui[i,,drop=F]
            if (length(ipos)) {
               ipzc=str2ind(ipos, nm_izc)
               ipos=str2ind(ipos, nm_i)
               ci[ipos]=as.numeric(zc+mi_zc[ipzc,,drop=F]%*%mic)
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
            capture.output(pinside <- put_inside(res$par, ui, ci), file=fclog)
            if (any(is.na(pinside))) {
               if (!is.null(attr(pinside, "err")) && attr(pinside, "err")!=0) {
                  # fatal error occured, don't reoptimize
                  cat(paste("put_inside", runsuf, ": ", attr(pinside, "mes"), "\\n", collapse=""), file=fcerr)
                  reopt=FALSE
               }
            } else if (!is.null(attr(pinside, "err")) && attr(pinside, "err")==0){
               # non fatal problem
               cat(paste("put_inside", runsuf, ": ", attr(pinside, "mes"), "\\n", collapse=""), file=fcerr)
            }
            # reoptimize
            if (reopt) {
               cat("Second zero crossing pass", runsuf, "\\n", sep="", file=fclog)
               capture.output(reso <- opt_wrapper(pinside, measurements, new.env()), file=fclog)
               if (reso$err || is.null(reso$par)) {
                  cat("second zero crossing pass: ", reso$mes, "\\n", sep="", file=fcerr)
               } else if (!is.null(reso$mes) && nchar(reso$mes)) {
                  cat("second zero crossing pass", runsuf, ": ", reso$mes, "\\n", sep="", file=fcerr)
               }
               if(!reso$err && !is.null(reso$par) && !any(is.na(reso$par))) {
                  param=reso$par
                  res=reso
                  jx_f=labargs$jx_f
               }
               if (any(is.na(reso$par))) {
                  reso$retres$jx_f=NULL # to avoid writing of huge data
                  obj2kvh(reso, "failed second pass optimization process information", fkvh)
                  cat("Second zero crossing pass failed. Keep free parameters from previous pass", runsuf, "\\n", file=fcerr, sep="")
               }
            }
            # last pass, free all zc constraints
            i=grep("^zc ", nm_i)
            if (length(i) > 0) {
               if (TIMEIT) {
                  cat("last zc : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
               }
               ui=ui[-i,,drop=F]
               ci=ci[-i]
               nm_i=nm_i[-i]
               cat("Last zero crossing pass (free of zc constraints)", runsuf, "\\n", sep="", file=fclog)
               capture.output(reso <- opt_wrapper(param, measurements, new.env()), file=fclog)
               if (reso$err || is.null(reso$par) || (!is.null(res$mes) && nchar(res$mes))) {
                  cat("last zero crossing (free of zc)", runsuf, ": ", reso$mes, "\\n", sep="", file=fcerr)
               }
               if(!res$err && !is.null(res$par) && !any(is.na(res$par))) {
                  param=reso$par
                  res=reso
                  jx_f=labargs$jx_f
               }
               if (any(is.na(res$par))) {
                  res$retres$jx_f=NULL # to avoid writing of huge data
                  obj2kvh(res, "failed last pass optimization process information", fkvh)
                  cat("Last zero crossing pass failed. Keep free parameters from previous passes", runsuf, "\\n", file=fcerr, sep="")
               }
            }
         } else {
            cat("After the first optimization, no zero crossing inequality was activated. So no reoptimization", runsuf, "\\n", sep="", file=fclog)
         }
      } # end if zero crossing
      param=res$par
      names(param)=nm_par
      if (excl_outliers != F) {
         # detect outliers
         if (TIMEIT) {
            cat("outliers: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         iva=!is.na(res$res)
         iout=which(rz.pval.bi(res$res) <= excl_outliers & iva)
         #cat("iout=", iout, "\\n", file=fclog)
         if (length(iout)) {
            measurements$outlier=iout
            cat("Excluded outliers at p-value ", excl_outliers, ":\\n",
               paste(nm_resid[iout], res$res[iout], sep="\\t", collapse="\\n"), "\\n", sep="", file=fclog)
            capture.output(reso <- opt_wrapper(param, measurements, new.env()), file=fclog)
            if (reso$err || is.null(reso$par) || (!is.null(reso$mes) && nchar(reso$mes))) {
               cat("wo outliers: ", reso$mes, "\\n", sep="", file=fcerr)
            }
            if (any(is.na(reso$par))) {
               cat("Optimization with outliers excluded has failed, run= ", runsuf, "\\n", file=fcerr, sep="")
               # continue without outlier exclusion
               measurements$outlier=NULL
            } else {
               res=reso
               param=reso$par
               names(param)=nm_par
               jx_f=labargs$jx_f
               obj2kvh(nm_resid[iout], "excluded outliers", fkvh)
            }
         } else {
            cat("Outlier exclusion at p-value "%s+%excl_outliers%s+%" has been requested but no outlier was detected at this level.", "\\n", sep="", file=fcerr)
         }
      }
      if (time_order=="1,2") {
         if (TIMEIT) {
            cat("order 2 : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         labargs$time_order="2" # continue with the 2-nd order
         capture.output(reso <- opt_wrapper(param, measurements, new.env()), file=fclog)
         if (reso$err || is.null(reso$par) || (!is.null(reso$mes) && nchar(reso$mes))) {
            cat("order2: ", reso$mes, "\\n", sep="", file=fcerr)
         }
         if (any(is.na(reso$par))) {
            cat("Optimization time_order 2 (in '1,2' suite) has failed, run=", runsuf, "\\n", file=fcerr, sep="")
         } else {
            res=reso
            param=reso$par
            names(param)=nm_par
            jx_f=labargs$jx_f
         }
      }
#browser()
      optinfo=list(
         "fitted parameters"=param,
         "last increment before backtracking"=res$lastp,
         "last increment after backtracking"=res$laststep,
         "iteration number"=res$it,
         "convergence history"=res$hist,
         "exit message"=res$mes
      )
      obj2kvh(optinfo, "optimization process information", fkvh)
   }
   if (TIMEIT) {
      cat("postopt : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
   }
   # active constraints
   ine=as.numeric(abs(ui%*%param-ci))<1.e-10
   if (any(ine)) {
      obj2kvh(nm_i[ine], "active inequality constraints", fkvh)
   }
   poolall[nm_poolf]=param[nm_poolf]

#browser()
   if (is.null(jx_f$jacobian)) {
      # final jacobian calculation
      capture.output(rres <- lab_resid(param, cjac=T, labargs), file=fclog)
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         close(fkvh)
         retcode[irun]=rres$err
         next
      }
   }
   rcost=cumo_cost(param, labargs)
   pres[,irun]=param
   costres[irun]=rcost
   obj2kvh(rcost, "final cost", fkvh)
#browser()
   resid=list(
      "labeled data"=jx_f$reslab,
      "measured fluxes"=jx_f$resflu,
      "measured pools"=jx_f$respool
   )
   obj2kvh(resid, "(simulated-measured)/sd_exp", fkvh)

   # simulated measurements -> kvh
   if (case_i) {
      simul=list(
         "labeled data"=jx_f$usm,
         "measured fluxes"=jx_f$simfmn,
         "measured pools"=jx_f$simpool
      )
   } else {
      simul=list(
         "labeled data (unscaled)"=jx_f$usimcumom,
         "measured fluxes"=jx_f$simfmn,
         "measured pools"=jx_f$simpool
      )
      if (nb_sc > 0) {
         simul=append(simul, list("labeled data (scaled)"=jx_f$usimcumom*c(1.,param)[ir2isc]), after=1)
      }
   }
   obj2kvh(simul, "simulated measurements", fkvh)

   # gradient -> kvh
   if (!is.null(jx_f$res)) {
      gr=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
      names(gr)=nm_par
      obj2kvh(gr, "gradient vector", fkvh)
   }
   colnames(jx_f$udr_dp)=nm_par
   obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh)
   
   # generalized inverse of non reduced jacobian
   svj=svd(jx_f$udr_dp)
   invj=svj$v%*%(t(svj$u)/svj$d)
   dimnames(invj)=rev(dimnames(jx_f$udr_dp))
   obj2kvh(invj, "generalized inverse of jacobian dr_dp (without 1/sd_exp)", fkvh)

   if (fullsys) {
      nm_flist=nm_list
      nm_flist$rcumo=nm_cumo
      nm_flist$rcumo_in_cumo=match(nm_rcumo, nm_cumo)
      nb_f$cumos=nb_cumos""")
    f.write("""
      nm_xi_f=c(%s)
      xi_f=c(%s)"""%(join(", ", netan["cumo_input"].keys(), '"', '"'),
      join(", ", netan["cumo_input"].values())))
    f.write("""
      names(xi_f)=nm_xi_f
      nm_flist$xi=nm_xi_f
      labargs$emu=F
      v=lab_sim(param, cjac=F, labargs)
      labargs$emu=emu
   } else {
      v=lab_sim(param, cjac=F, labargs)
   }
   # set final variable depending on param
   x=as.matrix(v$x)
   if (fullsys) {
      rownames(x)=nm_cumo
   } else {
      rownames(x)=nm_x
   }

   # write some info in result kvh
   obj2kvh(cumo2mass(x), "MID vector", fkvh)

   fwrv=v$lf$fwrv
   fallnx=v$lf$fallnx
   flnx=v$lf$flnx
   fgr=fallnx[nm_fgr]

   # keep last jx_f in jx_f_last
   while (sensitive=="mc") {
      if (TIMEIT) {
         cat("monte-ca: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
      }
      if(set_seed) {
         set.seed(seed)
      }
      # reference simulation corresponding to the final param
      refsim=new.env()
      for (nm_it in c("simlab", "simfmn", "simpool", "usm")) {
         assign(nm_it, jx_f[[nm_it]], envir=refsim)
      }
      
      # Monte-Carlo simulation in parallel way (if asked and possible)
      if (np > 1L) {
         # prepare cluster
         cl_type=ifelse(.Platform$OS.type=="unix", "FORK", "SOCK")
         if (cl_type=="FORK") {
            nodes=np
         } else {
            snow_here=suppressPackageStartupMessages(require(snow))
            if (snow_here) {
               nodes=rep("localhost", np)
            } else {
               cat("Monte-Cralo warning: 'snow' package is not installed =>\\n running Monte-Carlo simulations in sequential mode.\\nIf you don't want/can install snow package and want to avoid this message,\\nrun influx_s with an option '--np=1'.\\n", file=fcerr)
               np=1
               cl_type="None (sequential mode)"
            }
         }
      } else {
         cl_type="None (sequential mode)"
      }
      if (np > 1L) {
         # parallel execution
         cl=makeCluster(nodes, cl_type)
         if (cl_type=="SOCK") {
            if (TIMEIT) {
               cat("cl init : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
            }
            clusterEvalQ(cl, c(require(bitops), require(nnls)))
            if (TIMEIT) {
               cat("cl expor: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
            }
            clusterExport(cl, c("fcerr", "fclog", "lsi_fun", "cumo_jacob", "fx2jr", "trisparse_solv", "fwrv2Abr", "Heaviside", "df_dffp", "fallnx2fwrv", "dfcg2fallnx", "param2fl", "lab_sim", "is.diff", "lab_resid", "ui", "ci", "ep", "cp", "nlsic", "control_ftbl", "param", "norm2", "method", "sln", "opt_wrapper", "labargs", "dirx", "refsim"))
            if (TIMEIT) {
               cat("cl sourc: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
            }
            clusterEvalQ(cl, c(source(file.path(dirx, "tools_ssg.R")), source(file.path(dirx, "nlsic.R"))))
            if (TIMEIT) {
               cat("cl optim: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
            }
         }
         clusterSetRNGStream(cl)
         mc_res=parLapply(cl, seq_len(nmc), mc_sim, refsim=refsim, labargs=labargs)
         stopCluster(cl)
      } else {
         mc_res=lapply(1L:nmc, mc_sim, refsim=refsim, labargs=labargs)
      }
      free_mc=sapply(mc_res, function(l) {if (class(l)=="character" || is.na(l$cost) || l$err) { ret=rep(NA, nb_param+3) } else { ret=c(l$cost, l$it, l$normp, l$par) }; ret })
      if (length(free_mc)==0) {
         cat("Parallel exectution of Monte-Carlo simulations has failed.", "\\n", sep="", file=fcerr)
         free_mc=matrix(NA, nb_param+2, 0)
      }
      cost_mc=free_mc[1,]
      nmc_real=nmc-sum(is.na(free_mc[4,]))
      cat("monte-carlo\\n", file=fkvh)
      indent=1
      obj2kvh(cl_type, "cluster type", fkvh, indent)
      obj2kvh(avaco, "detected cores", fkvh, indent)
      avaco=max(1, avaco, na.rm=T)
      obj2kvh(min(avaco, np, na.rm=T), "used cores", fkvh, indent)
      cat("\\tfitting samples\\n", file=fkvh)
      indent=2
      obj2kvh(nmc, "requested number", fkvh, indent)
      obj2kvh(nmc_real, "calculated number", fkvh, indent)
      obj2kvh(nmc-nmc_real, "failed to calculate", fkvh, indent)
      # convergence section in kvh
      indent=1
      mout=rbind(round(free_mc[1:2,,drop=F], 2),
         format(free_mc[3,,drop=F], di=2, sci=T))
      dimnames(mout)=list(c("cost", "it.numb", "normp"), seq_len(ncol(free_mc)))
      obj2kvh(mout, "convergence per sample", fkvh, indent)
      # remove failed m-c iterations
      free_mc=free_mc[-(1:3),,drop=F]
      ifa=which(is.na(free_mc[1,]))
      if (length(ifa)) {
         if (ncol(free_mc) > length(ifa)) {
            cat("Some Monte-Carlo iterations failed.", "\\n", sep="", file=fcerr)
         }
         free_mc=free_mc[,-ifa,drop=F]
         cost_mc=cost_mc[-ifa]
      }
      if (nmc_real <= 1) {
         cat("No sufficient monter-carlo samples were succesfully calculated to do some statistics.", "\\n", sep="", file=fcerr)
         retcode[irun]=1
         break
      }
#browser()
      rownames(free_mc)=nm_par
      
      # cost section in kvh
      cat("\\tcost\\n", file=fkvh)
      indent=2
      obj2kvh(mean(cost_mc), "mean", fkvh, indent)
      obj2kvh(median(cost_mc), "median", fkvh, indent)
      obj2kvh(sd(cost_mc), "sd", fkvh, indent)
      obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent)
      obj2kvh(quantile(cost_mc, c(0.025, 0.975)), "ci", fkvh, indent)
      
      # free parameters section in kvh
      cat("\\tStatistics\\n", file=fkvh)
      mout=c()
      indent=2
      # param stats
      # mean
      parmean=apply(free_mc, 1, mean)
      # median
      parmed=apply(free_mc, 1, median)
#browser()
      # covariance matrix
      covmc=cov(t(free_mc))
      obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      # confidence intervals
      ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
      colnames(ci_mc)=c("CI 2.5%", "CI 97.5%", "CI length")
      mout=cbind(mout, mean=parmean, median=parmed, sd=sdmc,
         "rsd (%)"=sdmc*100/abs(parmean), ci_mc)
      obj2kvh(mout, "free parameters", fkvh, indent)

      # net-xch01 stats
      fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, labargs)$fallnx)
      fallnx=param2fl(param, labargs)$fallnx
      if (length(fallnx_mc)) {
         dimnames(fallnx_mc)[[1]]=nm_fallnx
         # form a matrix output
         fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0)
         # mean
#browser()
         parmean=apply(fallnx_mc, 1, mean)
         # median
         parmed=apply(fallnx_mc, 1, median)
         # covariance matrix
         covmc=cov(t(fallnx_mc))
         dimnames(covmc)=list(nm_fallnx, nm_fallnx)
         # sd
         sdmc=sqrt(diag(covmc))
         # confidence intervals
         ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)))
         ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
         ci_mc=cbind(ci_mc, ci_mc[,3]*100/abs(parmean))
         colnames(ci_mc)=c("CI 2.5%", "CI 97.5%", "CI 95% length", "relative CI (%)")
         fallout=cbind(fallout, mean=parmean, median=parmed, sd=sdmc,
            "rsd (%)"=sdmc*100/abs(fallnx), ci_mc)
         o=order(nm_fallnx)
         obj2kvh(fallout[o,,drop=F], "all net-xch01 fluxes", fkvh, indent)
         obj2kvh(covmc[o,o], "covariance of all net-xch01 fluxes", fkvh, indent)

         # fwd-rev stats
         fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, labargs)$fwrv)
         dimnames(fwrv_mc)[[1]]=nm_fwrv
         fallout=matrix(0, nrow=nrow(fwrv_mc), ncol=0)
         # mean
         parmean=apply(fwrv_mc, 1, mean)
         # median
         parmed=apply(fwrv_mc, 1, median)
         # covariance matrix
         covmc=cov(t(fwrv_mc))
         dimnames(covmc)=list(nm_fwrv, nm_fwrv)
         # sd
         sdmc=sqrt(diag(covmc))
         # confidence intervals
         ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)))
         ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
         ci_mc=cbind(ci_mc, ci_mc[,3]*100/abs(fwrv))
         dimnames(ci_mc)[[2]]=c("CI 2.5%", "CI 97.5%", "CI 95% length", "relative CI (%)")
         fallout=cbind(fallout, mean=parmean, median=parmed, sd=sdmc,
            "rsd (%)"=sdmc*100/abs(parmean), ci_mc)
         o=order(nm_fwrv)
         obj2kvh(fallout[o,,drop=F], "forward-reverse fluxes", fkvh, indent)
         obj2kvh(covmc[o,o], "covariance of forward-reverse fluxes", fkvh, indent)
      }
      break
   }
#browser()
   if (length(sensitive) && nchar(sensitive) && sensitive != "mc") {
      cat(paste("Unknown sensitivity '", sensitive, "' method chosen.", sep=""), "\\n", sep="", file=fcerr)
      retcode[irun]=1
   }

   if (TIMEIT) {
      cat("linstats: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
   }
   # Linear method based on jacobian x_f
   # reset fluxes and jacobians according to param
   if (is.null(jx_f$jacobian)) {
      capture.output(rres <- lab_resid(param, cjac=T, labargs), file=fclog)
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         close(fkvh)
         retcode[irun]=rres$err
         next
      }
   } # else use the last calculated jacobian

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
      cat(paste(if (nchar(runsuf)) runsuf%s+%": " else "", "Inverse of covariance matrix is numerically singular.\\nStatistically undefined parameter(s) seems to be:\\n",
         paste(nm_par[ibad], collapse="\\n"), "\\nFor more complete list, see sd columns in '/linear stats'\\nin the result file.", sep=""), "\\n", sep="", file=fcerr)
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
      covfl=crossprod(rtcov[, i, drop=F]%mmt%(rbind(diag(nb_ff+nb_fgr), dfl_dffg)%mrv%c(rep.int(1., nb_ff), fgr)))
      dimnames(covfl)=list(nm_flfd, nm_flfd)
      sdfl=sqrt(diag(covfl))
   } else {
      sdfl=rep(0., nb_fl)
      covfl=matrix(0., nb_fl, nb_fl)
   }
   fl=c(head(param, nb_ff), fgr, flnx)
   mtmp=cbind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl))
   rownames(mtmp)=nm_flfd
   o=order(nm_flfd)
   obj2kvh(mtmp[o,,drop=F], "net-xch01 fluxes (sorted by name)", fkvh, indent=1)
   obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1)

   # sd of all fwd-rev
   if (nb_ff > 0 || nb_fgr > 0) {
      i=1:nb_param
      i=c(head(i, nb_ff), tail(i, nb_fgr))
      covf=crossprod(tcrossprod_simple_triplet_matrix(rtcov[,i, drop=F], jx_f$df_dffp%mrv%c(rep.int(1., nb_ff), head(poolall[nm_poolf], nb_fgr))))
      dimnames(covf)=list(nm_fwrv, nm_fwrv)
      sdf=sqrt(diag(covf))
   } else {
      sdf=rep(0., length(fwrv))
   }
   mtmp=cbind(fwrv, sdf, sdf/abs(fwrv))
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
      poolall[nm_poolf]=param[nm_poolf]
      # cov poolf
      covpf=crossprod(rtcov[,nb_ff+nb_sc+1:nb_poolf, drop=F])
      dimnames(covpf)=list(nm_poolf, nm_poolf)
      sdpf[nm_poolf]=sqrt(diag(covpf))
   }
   if (length(poolall) > 0) {
      mtmp=cbind("value"=poolall, "sd"=sdpf, "rsd"=sdpf/poolall)
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
   if (!is.null(jx_f$res)) {
      nvres=sum(!is.na(jx_f$res))
      if (nvres >= nb_param) {
         khi2test=list("khi2 value"=rcost, "data points"=nvres,
            "fitted parameters"=nb_param, "degrees of freedom"=nvres-nb_param)
         khi2test$`khi2 reduced value`=khi2test$`khi2 value`/khi2test$`degrees of freedom`
         khi2test$`p-value, i.e. P(X^2<=value)`=pchisq(khi2test$`khi2 value`, df=khi2test$`degrees of freedom`)
         khi2test$conclusion=if (khi2test$`p-value, i.e. P(X^2<=value)` > 0.95) "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD" else "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD"
         obj2kvh(khi2test, "goodness of fit (khi2 test)", fkvh, indent=1)
      } else {
         cat(sprintf("khi2: Measurement number %d is lower than parameter number %d. Khi2 test cannot be done.\\n", nvres, nb_param), sep="", file=fcerr)
      }
   }
   if (prof) {
      Rprof(NULL)
   }
   close(fkvh)
   # write edge.netflux property
   fedge=file(file.path(dirw, sprintf("edge.netflux.%s%s.attrs", baseshort,  runsuf)), "w")
   cat("netflux (class=Double)\\n", sep="", file=fedge)
   nm_edge=names(edge2fl)
   cat(paste(nm_edge, fallnx[edge2fl], sep=" = "), sep="\\n" , file=fedge)
   close(fedge)

   # write edge.xchflux property
   fedge=file(file.path(dirw, sprintf("edge.xchflux.%s%s.attrs", baseshort,  runsuf)), "w")
   flxch=paste(".x", substring(edge2fl, 4), sep="")
   ifl=charmatch(flxch, substring(names(fallnx), 2))
   cat("xchflux (class=Double)\\n", sep="", file=fedge)
   cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
   close(fedge)

   # write node.log2pool property
   if (length(poolall)> 0) {
      fnode=file(file.path(dirw, sprintf("edge.xchflux.%s%s.attrs", baseshort,  runsuf)), "w")
      cat("log2pool (class=Double)\\n", sep="", file=fnode)
      nm_node=substring(names(poolall), 4)
      cat(paste(nm_node, log2(poolall), sep=" = "), sep="\\n" , file=fnode)
      close(fnode)
   }
}

pres=rbind(cost=costres, pres)
fco=file(file.path(dirw, sprintf("%s.pres.csv", baseshort)), open="w")
cat("row_col\t", file=fco)
write.table(file=fco, pres, row.n=T, quot=F, sep="\\t")
close(fco)
if (TIMEIT) {
   cat("rend    : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
""")
    f.write("""
# source files from FTBL/posttreat_R
postlist=strsplit("%(postlist)s", " *; *")[[1]]
for (post in postlist) {
   fpostR=file.path(dirw, post)
   if (!isTRUE(file.info(fpostR)$isdir)) {
      if (file.exists(fpostR)) {
         source(fpostR)
      } else {
         cat(sprintf("Posttreatment R file '%%s' does not exist. Ignored.\\n", fpostR), file=fcerr)
      }
   }
}
xgc=gc(verbose=FALSE) # to avoid the message "Error in (function (x)  : tentative d'appliquer un objet qui n'est pas une fonction"
close(fclog)
close(fcerr)
retcode=max(retcode)
if (format(parent.frame()) == format(.GlobalEnv)) {
   q("no", status=retcode)
}
"""%{
    "postlist": escape(netan["opt"].get("posttreat_R", ""), '\\"'),
})

    f.close()
    # try to make output files just readable to avoid later casual edition
    try:
        os.chmod(n_R, stat.S_IREAD)
    except:
        pass
