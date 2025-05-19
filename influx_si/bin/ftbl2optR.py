#!/usr/bin/env python3

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
   * nb_w, nb_rw, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
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
   * ir2isc - measure row to scale vector replicator
   * ci - inequalities for param use (ui%*%param-ci>=0)
   * measvec - measurement vector
   * fmn - measured net fluxes

Matrices:
   * Afl, qrAfl, invAfl,
   * p2bfl - helps to construct the rhs of flux system
   * mf, md - help to construct fallnx
   * mi - inequality matrix (ftbl content)
   * ui - inequality matrix (ready for param use)
   * measmat - for measmat*x+memaone=vec of simulated not-yet-scaled measurements

Functions:
   * lab_sim - translate param to flux and cumomer vector (initial approximation)
   * cumo_cost - cost function (chi2)
   * cumo_gradj - implicit derivative gradient
"""

# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used for sparse matrices.
# 2014-04-14 sokol:adapted for both "s" and "i" cases

#import pdb

import sys
import os
import stat
import time
import copy
import getopt
import math
from pathlib import Path
import influx_si
    
from tools_ssg import *
import C13_ftbl
import ftbl2code

def main(argv=sys.argv[1:], wout=sys.stdout.write, werr=sys.stderr.write):
    me=os.path.basename(__file__)
    dirbin=os.path.join(os.path.dirname(influx_si.__file__), "bin")
    sys.path.append(dirbin)

    def usage():
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--emu] [--clownr] [--tblimit[=0]] [--ropts ROPTS] network_name[.ftbl]\n")

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(argv, "h", ["help", "fullsys", "emu", "clownr", "tblimit=", "ropts=", "case_i", "ffguess", "dirres="])
    except getopt.GetoptError as err:
        #pass
        sys.stderr.write(str(err)+"\n")
        usage()
        return 1

    fullsys=False
    emu=False
    clownr=False
    ffguess=False
    ropts=['""']
    case_i=False
    tblimit=0
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            return 0
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
        elif o=="--tblimit":
            tblimit=int(a)
        elif o=="--ffguess":
            ffguess=True
        elif o=="--dirres":
           dirres=a
        else:
            #assert False, "unhandled option"
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            #pass
            #pdb.set_trace()
            raise Exception("unhandled option '%s'"%o)
    #aff("args", args);##
    #aff("opts", opts);##
    if len(args) != 1:
        usage()
        return 1
    org=os.path.basename(args[0])
    dirorg=os.path.dirname(args[0]) or '.'
    sys.tracebacklimit=tblimit

    # cut .ftbl if any
    if org[-5:]==".ftbl":
        org=org[:-5]
    fullorg=os.path.join(dirorg, org)

    #-->
    ftbl2code.case_i=case_i
    C13_ftbl.clownr=clownr
    C13_ftbl.ffguess=ffguess
    C13_ftbl.werr=werr
    C13_ftbl.wout=wout
    
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
    ftbl=C13_ftbl.ftbl_parse(n_ftbl, wout=wout, werr=werr)

    # analyse network
    # reload(C13_ftbl)

    netan=dict();
    C13_ftbl.ftbl_netan(ftbl, netan, emu, fullsys, case_i, wout=wout, werr=werr)
    wout.__self__.flush()
    werr.__self__.flush()
    # prepare rcumo system
    rAb=C13_ftbl.rcumo_sys(netan, emu)

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys, emu, ropts, dirres)

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
      i_pos=colSums(uip[i_nz > 0,,drop=FALSE] > 0)
      met_low=nm_poolf[i_pos > 0]
      # number of negative coeffs for alone metabs (i.e. upper limit is set)
      i_neg=colSums(uip[i_nz > 0,,drop=FALSE] < 0)
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
names(poolmdev)=nm_poolm

# simulated metabolite measurements are calculated as
# measmatpool*poolall=>poolm
measmatpool=matrix(0., nrow=nb_poolm, ncol=length(poolall))
dimnames(measmatpool)=list(nm_poolm, nm_poolall)
i=matrix(1+c(%(imeasmatpool)s), ncol=2, byrow=TRUE)
measmatpool[i]=1.

labargs=new.env()
labargs$labargs2=new.env()
    """%{
    "nb_poolm": len(netan["metab_measured"]),
    "nm_poolm": join(", ", list(netan["metab_measured"].keys()), '"pm:', '"'),
    "v_poolm": join(", ", (item["val"] for item in list(netan["metab_measured"].values()))).replace("nan", "NA"),
    "poolmdev": join(", ", (item["dev"] for item in list(netan["metab_measured"].values()))),
    "imeasmatpool": join(", ", valval((ir,netan["vpool"]["all2i"][m]) for (ir, ml) in enumerate(netan["metab_measured"].keys()) for m in ml.split("+"))),
    })
    if case_i:
        f.write("""
## variables for isotopomer kinetics
tstart=0.
tmax=c(%(tmax)s)
dt=c(%(dt)s)

# funlab list
funlabli=structure(list(%(funlabli)s), names=nm_exp)
funlabli=structure(lapply(nm_exp, function(nm) {
   structure(lapply(names(funlabli[[nm]]), function(met) {
      structure(lapply(names(funlabli[[nm]][[met]]), function(ni) {
         rcode=funlabli[[nm]][[met]][[ni]]
         v=try(parse(text=rcode), silent=TRUE)
         if (inherits(v, "try-error"))
            stop_mes("Error in parsing R code '", rcode, "' for label '", met, "#", ni, "' in '", nm, "':\n", v, file=fcerr)
         v
      }), names=names(funlabli[[nm]][[met]]))
   }), names=names(funlabli[[nm]]))
}), names=nm_exp)
if (inherits(funlabli, "try-error"))
   stop_mes(funlabli, file=fcerr)
funlabR=paste(dirw, c(%(funlabR)s), sep="/")
funlabR[funlabR == paste0(dirw, "/")]="" # set to "" where file names are empty

# read measvecti from file(s) specified in ftbl(s)
flabcin=c(%(flabcin)s)
measvecti=ti=tifull=tifull2=vector("list", nb_exp)
nb_ti=nb_tifu=nb_tifu2=integer(nb_exp)
nsubdiv_dt=pmax(1L, as.integer(c(%(nsubdiv_dt)s)))
nb_f$ipf2ircumo=nb_f$ipf2ircumo2=list()
nminvm=nm_poolall[matrix(unlist(strsplit(nm_rcumo, ":", fixed=TRUE)), ncol=2L, byrow=TRUE)[,1L]]
if (fullsys) {
   nminvmf=nm_poolall[matrix(unlist(strsplit(nm_cumo, ":", fixed=TRUE)), ncol=2L, byrow=TRUE)[,1L]]
   # if simulations are in emu, prepare transition matrix from full cumos to simulated emus
   if (emu) {
      tmp=emu_in_cumo(nm_list$emu, nm_list$cumo)
      labargs$memu_in_cumo=tmp$m
      labargs$cemu_in_cumo=tmp$cnst
      if (case_i) {
         labargs$labargs2$memu_in_cumo=tmp$m
         labargs$labargs2$cemu_in_cumo=tmp$cnst
      }
   }
}
#browser()
for (iexp in seq_len(nb_exp)) {
   if (tmax[iexp] < 0) {
      stop_mes(sprintf("The parameter tmax must not be negative (tmax=%%g in '%%s.ftbl')", tmax[iexp], nm_exp[iexp]), file=fcerr)
   }
   if (dt[iexp] <= 0) {
      stop_mes(sprintf("The parameter dt must be positive (dt=%%g in '%%s.ftbl')", dt[iexp], nm_exp[iexp]), file=fcerr)
   }
   if (nchar(flabcin[iexp])) {
      if (substr(flabcin[iexp], 1, 1) == "/")
         flabcin[iexp]=file.path(flabcin[iexp])
      else
         flabcin[iexp]=file.path(dirw, flabcin[iexp])
      if (!file.exists(flabcin[iexp])) {
         # last resort, see in _res/tmp
         ftmp=file.path(dirw, paste0(nm_exp[[iexp]],"_res"), "tmp", basename(flabcin[iexp]))
         if (file.exists(ftmp))
            flabcin[iexp]=ftmp
      }
      measvecti[[iexp]]=try(as.matrix(read.table(flabcin[iexp], header=TRUE, row.names=1, sep="\t", check=FALSE, comment="#", strip.white=TRUE)), silent=TRUE)
      if (inherits(measvecti[[iexp]], "try-error")) {
         # try with comment '//'
         tmp=try(kvh::kvh_read(flabcin[iexp], comment_str = "//", strip_white = FALSE, skip_blank = TRUE, split_str = "\t", follow_url = FALSE), silent=TRUE)
         if (inherits(tmp, "try-error"))
            stop_mes("Error while reading '", flabcin[iexp], "' from '", nm_exp[iexp], "':\n", tmp, file=fcerr)
         nb_col=sapply(tmp, length)
         if (any(ibad <- nb_col != nb_col[1]))
            stop_mes("Column number varies in '", flabcin[iexp], "'. First row has ", nb_col[1], " columns while the following rows differ:\n\t", paste(c("row", which(ibad)), c("col_nb", nb_col[ibad]), sep="\t", collapse="\n\t"))
         tmp=lapply(tmp, trimws)
         tmp=do.call(rbind, tmp)
         tmp=structure(tmp[-1L,, drop=FALSE], dimnames=list(rownames(tmp)[-1L], tmp[1L,]))
         suppressWarnings(storage.mode(tmp) <- "double")
         measvecti[[iexp]]=tmp
      }
      nm_row=rownames(measvecti[[iexp]])
      # put in the same row order as simulated measurements
      # check if nm_meas are all in rownames
#browser()
      if (all(nm_meas[[iexp]] %%in%% nm_row)) {
         measvecti[[iexp]]=measvecti[[iexp]][nm_meas[[iexp]],,drop=FALSE]
      } else {
         # try to strip row number from measure id
         nm_strip=sapply(strsplit(nm_meas[[iexp]], ":", fixed=TRUE), function(v) {
            paste(c(v[-length(v)], ""), sep="", collapse=":")
         })
         im=pmatch(nm_strip, nm_row)
         ina=is.na(im)
         if (any(ina)) {
            mes=paste("Cannot match the following measurement(s) in the file '", flabcin[iexp], "':\\n", paste(nm_meas[[iexp]][ina], sep="", collapse="\\n"), "\\n", sep="", collapse="")
            stop_mes(mes, file=fcerr)
         }
         measvecti[[iexp]]=measvecti[[iexp]][im,,drop=FALSE]
         #stopifnot(all(!is.na(measvecti)))
         if (typeof(measvecti[[iexp]])!="double") {
            # check for weird  entries
            tmp=measvecti[[iexp]]
            suppressWarnings(storage.mode(tmp) <- "double")
            if (any(ibad <- is.na(tmp) & !is.na(measvecti[[iexp]]))) {
               ibad=which(ibad)[1L]
               stop_mes("This entry '", measvecti[[iexp]][ibad], "' could not be converted to real number (", flabcin[iexp], ")", file=fcerr)
            } else if (!noopt) {
               stop_mes("Entries in file '", flabcin[iexp], "' could not be converted to real numbers", file=fcerr)
            }
         }
         if (!noopt && all(is.na(measvecti[[iexp]]))) {
            stop_mes("All entries in file '", flabcin[iexp], "' are NA (non available).", file=fcerr)
         }
      }
      ti[[iexp]]=as.double(colnames(measvecti[[iexp]]))
      if (any(is.na(ti[[iexp]]))) {
         mes=sprintf("Some time moments (in column names) could not be converted to real numbers in the file '%%s'\\nConverted times:\\n%%s", flabcin[[iexp]], join("\\n", ti[[iexp]]))
         stop_mes(mes, file=fcerr)
      }
      if (length(ti[[iexp]]) < 1L) {
         mes=sprintf("No column found in the file '%%s'", flabcin[[iexp]])
         stop_mes(mes, file=fcerr)
      }
      if (!all(diff(ti[[iexp]]) > 0.)) {
         mes=sprintf("Time moments (in column names) are not monotonously increasing in the file '%%s'", flabcin[[iexp]])
         stop_mes(mes, file=fcerr)
      }
      if (ti[[iexp]][1L] <= 0.) {
         mes=sprintf("The first time moment cannot be negative or 0 in the file '%%s'", flabcin[[iexp]])
         stop_mes(mes, file=fcerr)
      }
      if (ti[[iexp]][1L] != 0.) {
         ti[[iexp]]=c(tstart, ti[[iexp]])
      }
      i=which(ti[[iexp]]<=tmax[[iexp]])
      ti[[iexp]]=ti[[iexp]][i]
      if (tmax[[iexp]] == Inf) {
         tmax[[iexp]]=max(ti[[iexp]])
      }
      measvecti[[iexp]]=measvecti[[iexp]][,i[-1]-1,drop=FALSE]
   } else {
      if (tmax[[iexp]] == Inf) {
         stop_mes(sprintf("Maximal value for time is Inf (probably 'tmax' field is not set in OPTIONS section: '%%s.ftbl')", nm_exp[[iexp]]), file=fcerr)
      }
      ti[[iexp]]=seq(tstart, tmax[[iexp]], by=dt[iexp])
      if (optimize) {
         cat(sprintf("***Warning: a fitting is requested but no file with label data is provided by 'file_labcin' option in '%%s.ftbl' file.
   The fitting is ignored as if '--noopt' option were asked.\\n", nm_exp[[iexp]]), file=fclog)
         optimize=FALSE
      }
      # create measvecti with NA
      measvecti[[iexp]]=structure(matrix(NA, nrow=length(measvec[[iexp]]), ncol=length(ti[[iexp]])-1L), dimnames=list(nm_list$meas[[iexp]], ti[[iexp]][-1L]))
   }
   # recalculate nb_exp from measvecti
   nb_meas=sapply(measvecti, NROW)
   nb_meas_cumo=c(0., cumsum(nb_meas[-nb_exp]))
   iexp_meas=lapply(seq_len(nb_exp), function(iexp) seq_len(nb_meas[iexp])+nb_meas_cumo[iexp])
   nb_f$nb_meas=nb_meas

   nb_ti[iexp]=length(ti[[iexp]])
   if (nb_ti[iexp] < 2L) {
      mes=sprintf("After filtering by tmax, only %%d time moments are kept for experiment '%%s'. It is not sufficient.", nb_ti[iexp], nm_exp[iexp])
      stop_mes(mes, file=fcerr)
   }
   
   # divide the first time interval by n1 geometric intervals
   tifull[[iexp]]=ti[[iexp]]
   
   # divide each time interval by nsubdiv_dt
   dtmp=diff(tifull[[iexp]])
   dtmp=rep(dtmp/nsubdiv_dt[iexp], each=nsubdiv_dt[iexp])
   tifull[[iexp]]=c(tifull[[iexp]][1L], cumsum(dtmp))
   nb_tifu[iexp]=length(tifull[[iexp]])
   
   tifull2[[iexp]]=c(tifull[[iexp]][1L], tifull[[iexp]][1L]+cumsum(rep(diff(tifull[[iexp]])/2., each=2L)))
   nb_tifu2[iexp]=length(tifull2[[iexp]])
   
   if (length(ijpwef[[iexp]])) {
      # vector index for many time points
      ijpwef[[iexp]]=cbind(ijpwef[[iexp]][,1L], rep(seq_len(nb_ti[[iexp]]-1L), each=nrow(ijpwef[[iexp]])), ijpwef[[iexp]][,2L])
      dp_ones[[iexp]]=matrix(aperm(array(dp_ones[[iexp]], c(dim(dp_ones[[iexp]]), nb_ti[[iexp]]-1L)), c(1L, 3L, 2L)), ncol=nb_poolf)
   }

      """%{
          "dt": join(", ", netan["opt"]["dt"]),
          "tmax": join(", ", ["Inf" if math.isinf(v) else v for v in netan["opt"]["tmax"]]),
          "flabcin": join(", ", netan["opt"]["file_labcin"], '"', '"'),
          "nsubdiv_dt": join(", ", netan["opt"]["nsubdiv_dt"]),
          "funlabli": join(", ", (C13_ftbl.mkfunlabli(v) for v in netan["funlab"])),
          "funlabR": join(", ", netan["opt"]["funlabR"], '"', '"')
         })

        f.write("""
   # prepare mapping of metab pools on cumomers
   nb_f$ipf2ircumo[[iexp]]=nb_f$ipf2ircumo2[[iexp]]=list()
   for (iw in seq_len(nb_rw)) {
      ix=seq_len(nb_rcumos[iw])
      ipf2ircumo=ipf2ircumo2=match(nminvm[nbc_rcumos[iw]+ix], nm_poolf, nomatch=0L)
      dims=c(1L, nb_rcumos[iw], ifelse(emu, iw, 1L), nb_tifu[iexp]-1L)
      dims2=c(1L, nb_rcumos[iw], ifelse(emu, iw, 1L), nb_tifu2[iexp]-1L)
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
      i=i[i[,1L]!=0L,,drop=FALSE]
      i2=i2[i2[,1L]!=0L,,drop=FALSE]
      # put the poolf column last
      nb_f$ipf2ircumo[[iexp]][[iw]]=i[, c("ic", "iw", "ipoolf", "iti"), drop=FALSE]
      nb_f$ipf2ircumo2[[iexp]][[iw]]=i2[, c("ic", "iw", "ipoolf", "iti"), drop=FALSE]
   }
   if (fullsys) {
      # prepare mapping of metab pools on cumomers for full system (emu is FALSE here)
      nb_f$ipf2icumo[[iexp]]=nb_f$ipf2icumo2[[iexp]]=list()
      for (iw in seq_len(nb_w)) {
         ix=seq_len(nb_cumos[iw])
         ipf2icumo=ipf2icumo2=match(nminvmf[nbc_cumos[iw]+ix], nm_poolf, nomatch=0L)
         dims=c(1L, nb_cumos[iw], 1L, nb_tifu[iexp]-1L)
         dims2=c(1L, nb_cumos[iw], 1L, nb_tifu2[iexp]-1L)
         i=as.matrix(ipf2icumo)
         i2=as.matrix(ipf2icumo2)
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
         i=i[i[,1L]!=0L,,drop=FALSE]
         i2=i2[i2[,1L]!=0L,,drop=FALSE]
         # put the poolf column last
         nb_f$ipf2icumo[[iexp]][[iw]]=i[, c("ic", "iw", "ipoolf", "iti"), drop=FALSE]
         nb_f$ipf2icumo2[[iexp]][[iw]]=i2[, c("ic", "iw", "ipoolf", "iti"), drop=FALSE]
      }
   }
}
xil=vector("list", nb_exp)
xi2=vector("list", nb_exp)
for (iexp in seq(nb_exp)) {
   fli=funlabli[[iexp]]
   if (length(fli) == 0) {
      # replicate first column in xi as many times as there are time points
      if (time_order == "2" || time_order == "1,2")
         xi2[[iexp]]=matrix(xi[[iexp]], nrow=length(xi[[iexp]]), ncol=nb_tifu2[[iexp]])
      xil[[iexp]]=matrix(xi[[iexp]], nrow=length(xi[[iexp]]), ncol=nb_tifu[[iexp]])
   }  else {
      # use funlab
      envfunlab=new.env() # funlab code won't see influx's variables
      if (nchar(funlabR[[iexp]])) {
         if (file.exists(funlabR[[iexp]])) {
            es=try(source(funlabR[[iexp]], local=envfunlab, echo=FALSE), outFile=fcerr)
         } else {
            stop_mes("funlab script R '", funlabR[[iexp]], "' from '", nm_exp[[iexp]],"' does not exist.", file=fcerr)
         }
      }
      xil[[iexp]]=funlab(tifull[[iexp]], nm_inp, fli, envfunlab, emu, nm_exp[[iexp]], fcerr)
      if (time_order == "2" || time_order == "1,2")
         xi2[[iexp]]=funlab(tifull2[[iexp]], nm_inp, fli, envfunlab, emu, nm_exp[[iexp]], fcerr)
   }
   rownames(xil[[iexp]])=if (emu) nm_xiemu else nm_xi
}
xi=xil

nb_f$ip2ircumo=match(nminvm, nm_poolall)
if (fullsys)
   nb_f$ip2icumo=match(nminvmf, nm_poolall)
nb_f$tifu=nb_tifu
nb_f$tifu2=nb_tifu2

# label state at t=0 (by default=0 but later it should be able to be specified by user)
x0=NULL
nb_f$ti=nb_ti
    """)
    f.write("""
# gather all measurement information
measurements=list(
   vec=list(labeled=measvec, flux=fmn, pool=vecpoolm, kin=if (case_i) measvecti else NULL),
   dev=list(labeled=measdev, flux=fmndev, pool=poolmdev, kin=if (case_i) lapply(structure(seq(nb_exp), names=names(measdev)), function(i) {v=measdev[[i]]; nbc=if (is.null(measvecti[[i]])) 0 else ncol(measvecti[[i]]); suppressWarnings(matrix(v, nrow=length(v), ncol=nbc))}) else NULL),
   mat=list(labeled=measmat, flux=ifmn, pool=measmatpool),
   one=list(labeled=memaone)
)
nm_resid=c(if (case_i) unlist(lapply(seq_len(nb_exp), function(iexp) {m=outer(rownames(measvecti[[iexp]]), ti[[iexp]][-1L], paste, sep=", t="); if (length(m) > 0L) paste(iexp, m, sep=":", recycle0=TRUE) else character(0L)})) else unlist(lapply(seq_len(nb_exp), function(iexp) paste(iexp, nm_meas[[iexp]], sep=":"))), nm_fmn, nm_poolm)
nm_list$resid=nm_resid

if (TIMEIT) {
   cat("preopt  : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
#browser()
names(param)=nm_par
# prepare series of starting points
if (nchar(fseries) > 0) {
   pstart=as.matrix(read.table(file.path(dirw, fseries), header=TRUE, row.names=1L, sep="\\t"))
   # skip parameters (rows) who's name is not in nm_par
   i=rownames(pstart) %in% nm_par
   if (!any(i)) {
      stop_mes("Option --fseries is used but no free parameter with known name is found.\\n", file=fcerr)
   }
   pstart=pstart[i,,drop=FALSE]
   cat("Using starting values form '", fseries, "' for the following free parameters:\\n", paste(rownames(pstart), collapse="\\n"), "\\n", sep="", file=fclog)
   nseries=ncol(pstart)
   if (initrand) {
      # fill the rest of rows with random values
      i=nm_par %in% rownames(pstart)
      n=sum(!i)
      pstart=rbind(pstart, structure(matrix(runif(n*nseries), n, nseries), dimnames=list(NULL, sprintf(paste0("V%0", ceiling(log10(nseries+1)), "d"), seq(nseries)))))
      rownames(pstart)=c(rownames(pstart)[seq_len(nb_param-n)], nm_par[!i])
   }
   if (nchar(iseries) > 0) {
      iseries=unique(as.integer(eval(parse(t="c("%s+%iseries%s+%")"))))
      iseries=iseries[iseries<=nseries]
      # subsample
      pstart=pstart[,iseries, drop=FALSE]
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
   dimnames(pstart)=list(nm_par, sprintf(paste0("V%0", ceiling(log10(nseries+1)), "d"), seq(nseries)))
   if (initrand) {
      # fill pstart with random values
      pstart[]=runif(length(pstart))
   }
   # subsample
   pstart=pstart[,iseries, drop=FALSE]
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
nb_f=as.environment(nb_f)

nbc_x=c(0, cumsum(nb_x))
nb_f$nbc_x=nbc_x

# fixed part of jacobian (unreduced by SD)
# measured fluxes
dufm_dp=cbind(dufm_dff(nb_f, nm_list), matrix(0, nrow=nb_fmn, ncol=nb_sc_tot+nb_poolf))
dimnames(dufm_dp)=list(nm_fmn, nm_par)

# measured pools
dupm_dp=matrix(0., nb_poolm, nb_ff+nb_sc_tot)
if (nb_poolf > 0L) {
   dupm_dp=cbind(dupm_dp, measurements$mat$pool[,nm_list$poolf, drop=FALSE])
}
dimnames(dupm_dp)=list(rownames(measurements$mat$pool), nm_par)

#browser()
# prepare argument list for passing to label simulating functions
nm_labargs=c("jx_f", "nb_f", "nm_list", "nb_x", "invAfl", "p2bfl", "g2bfl", "bp", "fc", "xi", "spa", "spaf", "emu", "pool", "measurements", "ipooled", "ir2isc",  "nb_w", "nb_rw", "nbc_x", "measmat", "memaone", "dufm_dp", "dupm_dp", "pwe", "ipwe", "ip2ipwe", "pool_factor", "ijpwef", "ipf_in_ppw", "meas2sum", "dp_ones", "clen", "dirr", "dirw", "dirres", "baseshort", "case_i", "nb_exp", "noscale", "dpw_dpf", "fclog", "fcerr")
""")
    if case_i:
        f.write('nm_labargs=c(nm_labargs, "ti", "tifull", "tifull2", "x0", "time_order")')
    f.write("""
if (TIMEIT) {
   cat("labargs : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
tmp=lapply(nm_labargs, function(nm) assign(nm, get(nm), labargs))
#for (nm in nm_labargs) {
#   labargs[[nm]]=get(nm)
#}
#labargs[["nm"]]=labargs[["nm_list"]]

# prepare labargs2 if time_order includes 2
if (case_i && (time_order == "2" || time_order == "1,2")) {
   list2env(as.list(labargs), envir=labargs$labargs2)
   labargs$labargs2$tifull=tifull2
   labargs$labargs2$xi=xi2
   labargs$labargs2$jx_f=new.env()
   labargs$labargs2$nb_f=new.env()
   list2env(as.list(labargs$nb_f), envir=labargs$labargs2$nb_f)
   labargs$labargs2$nb_f$ipf2ircumo=nb_f$ipf2ircumo2
   labargs$labargs2$nb_f$tifu=nb_f$tifu2
}

# formated output in kvh file
if (write_res && wkvh) {
   fkvh_saved=file.path(dirres, "tmp", sprintf("%s_res.kvh", baseshort))
} else {
   fkvh_saved=NULL
}
""")
    f.write(r"""
retcode=numeric(nseries)
cl_type="PSOCK"
cl=NULL
if ((case_i && (time_order %in% c("1,2", "2"))) || sensitive == "mc") {
   if (np > 1L) {
      # prepare cluster
      nodes=if (sensitive == "mc" && !parR) np else 2

      cl=makeCluster(nodes, cl_type) #)
      #cl=makeCluster(1L, cl_type, manual=TRUE, outfile="")
#cat("make cluster=")
#print(cl[[1]])
      labargs[["cl"]]=cl
      nodes=length(cl)
      if (TIMEIT) {
         cat("cl expor: ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
      }
      clusterExport(cl, c("lsi_fun", "df_dffp", "lab_sim", "is.diff", "lab_resid", "ui", "ci", "ep", "cp", "control_ftbl", "methods", "sln", "labargs", "dirr", "emu", "%stm%", "case_i", "time_order", "fullsys", "nm_inp"))
      if (TIMEIT) {
         cat("cl sourc: ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
      }
      clusterEvalQ(cl, {
         #idth=myinfo$id
         suppressPackageStartupMessages(library(nlsic))
         suppressPackageStartupMessages(library(slam)) # for quick sparse matrices
         suppressPackageStartupMessages(library(Rcpp))
         suppressPackageStartupMessages(library(RcppArmadillo))
         suppressPackageStartupMessages(library(rmumps))
         suppressPackageStartupMessages(library(arrApply)) # for fast apply() on arrays
         suppressPackageStartupMessages(library(multbxxc))
         compiler::enableJIT(0)
         source(file.path(dirr, "tools_ssg.R"))
         #source(file.path(dirr, "nlsic.R"))
         source(file.path(dirr, "opt_cumo_tools.R"))
         source(file.path(dirr, "opt_icumo_tools.R"))
         labargs$spa=sparse2spa(labargs$spa)
         if (case_i && (time_order == "2" || time_order == "1,2")) {
            labargs$labargs2$spa=labargs$spa
         }
         if (fullsys) {
            labargs$spaf=sparse2spa(labargs$spaf)
            if (case_i && (time_order == "2" || time_order == "1,2")) {
               labargs$labargs2$spaf=labargs$spaf
            }
         }
#cat("evalQ idth=", idth, "\n")
#print(labargs)
#print(labargs$labargs2)
#print(labargs$spa)
#print(labargs$labargs2$spa)
#print(list(cre_cl="", labargs=labargs, spa_a1=labargs$spa[[1L]]$a))
         NULL
      })
      clusterSetRNGStream(cl)
      # set worker id
      idw=parLapply(cl, seq_along(cl), function(i) assign("idw", i, envir=.GlobalEnv))

   } else {
      labargs$cl=NULL
   }
}
for (irun in seq_len(nseries)) {
   if (TIMEIT) {
      cat(sprintf("run %4d: %s cpu=%g\n", irun, format(Sys.time()), proc.time()[1]), "\n", sep="", file=fclog)
   }
   param[nm_pseries]=pstart[nm_pseries, irun]
#browser()
   if (nseries > 1) {
      runsuf="." %s+% colnames(pstart)[irun]
   } else {
      runsuf=""
   }
   if (length(nseries) > 0) {
      cat("Starting point", runsuf, "\n", sep="", file=fclog)
   }
   # prepare kvh file name
   if (write_res && wkvh) {
      fkvh=file(substring(fkvh_saved, 1, nchar(fkvh_saved)-4) %s+% runsuf %s+% ".kvh", "w");
   } else {
      fkvh=NULL
   }

   # remove zc inequalities from previous runs
   izc=grep("^zc ", nm_i)
   if (length(izc)) {
      ui=ui[-izc,,drop=FALSE]
      ci=ci[-izc]
      nm_i=rownames(ui)
   }
   # check if initial approximation is feasible
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   # set tolerance for inequality
   tol_ineq=if ("BFGS" %in% methods) 0. else tol
   nbad=sum(ineq <= -tol_ineq)
   if (anyNA(param))
      stop_mes("NA found in initial 'param' values:\n\t",
      paste0(apply(cbind(nm_par, param), 1, paste0, collapse="\t"), collapse="\n\t"))
   if (nbad > 0) {
      if (TIMEIT) {
         cat(sprintf("put_ins : %s cpu=%g\n", format(Sys.time()), proc.time()[1]), "\n", sep="", file=fclog)
      }
      cat("The following ", nbad, " inequalities are not respected at starting point", runsuf, ":\n", sep="", file=fclog)
      i=ineq[ineq<= -tol_ineq]
      cat(paste(names(i), i, sep="\t", collapse="\n"), "\n", sep="", file=fclog)
      # put them inside
#browser()
      if (write_res) {
         capture.output(pinside <- put_inside(param, ui, ci, rcond=1e7), file=fclog)
      } else {
         pinside <- put_inside(param, ui, ci, rcond=1e7)
      }
      if (anyNA(pinside)) {
         if (!is.null(attr(pinside, "err")) && attr(pinside, "err")!=0) {
            # fatal error occured
            cat("Infeasible inequalities at starting point.\n", file=fcerr)
            cat("put_inside", runsuf, ": ", attr(pinside, "mes"), "\n",
               file=fcerr, sep="")
            #close(fkvh)
            retcode[irun]=attr(pinside, "err")
            next;
         }
      } else if (!is.null(attr(pinside, "err")) && attr(pinside, "err")==0) {
         # non fatal problem
         cat(paste("***Warning: put_inside: ", attr(pinside, "mes"), collapse=""), "\n", file=fclog)
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
   if (!noopt && zerocross && length(grep("^[df]\\.n\\.", nm_fallnx))>0) {
      if (TIMEIT) {
         cat("zc ineq : ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
      }
#browser()
      # prepare fluxes that are already in inequalities in alone mode
      ige=names(which(apply(mi, 1L, function(v) diff(range(v))==1 && sum(v)==1) & li>=0))
      ige=nm_dfn[unique(c(
         sub("^n:.+<=(.+)$", "\\1", grep("^n:.+<=.+$", ige, v=TRUE)),
         sub("^[df]\\.n\\.(.+)>=.+$", "\\1", grep("^[df]\\.n\\..+>=.+$", ige, v=TRUE)),
         sub("^inout [df]\\.n\\.(.+)>=.+$", "\\1", grep("^inout [df]\\.n\\..+>=.+$", ige, v=TRUE))
      ))]
      ile=which(apply(mi, 1L, function(v) diff(range(v))==1 && sum(v)==-1)&li>=0)
      ile=nm_dfn[unique(c(
         sub("^n:.+<=(.+)$", "\\1", grep("^n:.+<=.+$", ile, v=TRUE)),
         sub("^[df]\\.n\\.(.+)>=.+$", "\\1", grep("^[df]\\.n\\..+>=.+$", ile, v=TRUE)),
         sub("^inout [df]\\.n\\.(.+)>=.+$", "\\1", grep("^inout [df]\\.n\\..+>=.+$", ile, v=TRUE))
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
      ui_zc=cbind(mi_zc%*%(md%*%invAfl%stm%p2bfl+mf),
         matrix(0., nrow=nrow(mi_zc), ncol=nb_sc_tot))
      if (nb_fgr > 0) {
         ui_zc=cbind(ui_zc, mi_zc%*%((md%*%invAfl%stm%g2bfl)+mg*nb_f$mu))
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

      inotsat=ci_zc[zi]>tol_ineq
      if (any(inotsat)) {
         cat("***Warning: the following constant zc inequalities are not satisfied:\n", file=fclog)
         cat(nm_izc[zi][inotsat], sep="\n", file=fclog)
      }
      ui_zc=ui_zc[!zi,,drop=FALSE]
      ci_zc=ci_zc[!zi]
      nm_izc=nm_izc[!zi]
      mi_zc=mi_zc[!zi,,drop=FALSE]

      # remove redundant/contradictory inequalities
      nb_zc=nrow(ui_zc)
      nb_i=nrow(ui)
      ired=c()
      tui=t(ui)
      uzcd=sapply(seq_len(nb_zc), function(i) apply(abs(tui-ui_zc[i,]), 2L, max))
      uzcs=sapply(seq_len(nb_zc), function(i) apply(abs(tui+ui_zc[i,]), 2L, max))
      czcd=abs(outer(abs(ci), abs(ci_zc), "-"))
      ired=which(apply((uzcd < tol_ineq | uzcs < tol_ineq) & czcd <= 1.e-2, 2, any))
      
      if (length(ired) > 0L) {
         # remove all ired inequalities
         cat("The following ", length(ired), " zerocross inequalities are redundant and are removed:\n", paste(nm_izc[ired], collapse="\n"), "\n", sep="", file=fclog)
         ui_zc=ui_zc[-ired,,drop=FALSE]
         ci_zc=ci_zc[-ired]
         nm_izc=nm_izc[-ired]
         mi_zc=mi_zc[-ired,,drop=FALSE]
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
    if not case_i:
         f.write("""
   # set initial scale values to sum(measvec*simlab/dev**2)/sum(simlab**2/dev**2)
   # for corresponding measurements
   if (nb_sc_tot > 0) {
      if (TRUE) { # always estimate scaling params even for --noopt
         if (TIMEIT) {
            cat("res esti: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         if (write_res) {
            capture.output(rres <- lab_resid(param, cjac=FALSE, labargs), file=fclog)
         } else {
            rres <- lab_resid(param, cjac=FALSE, labargs)
         }
         if (!is.null(rres$err) && rres$err) {
            cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
            #close(fkvh)
            retcode[irun]=rres$err
            next
         }
         if (sum(is.infinite(rres$res))) {
            cat("Infinite values appeared in residual vector (at init scale values)", file=fcerr)
            retcode[irun]=1
            #close(fkvh)
            next
         }
         for (iexp in seq_len(nb_exp)) {
            simlab=jx_f$usimlab[[iexp]]
            measinvvar=1./measurements$dev$labeled[[iexp]]**2
            ms=measvec[[iexp]]*simlab*measinvvar
            ss=simlab*simlab*measinvvar
            # get only valid measurements
            iva=!is.na(ms)
            for (i in nb_ff+nb_sc_base[iexp]+seq_len(nb_sc[[iexp]])) {
               im=(ir2isc[[iexp]]==(i+1)) & iva
               if (sum(im) < 2) {
                  mes=sprintf("scaling: no sufficient valid data for scaling factor '%s'\\n", nm_par[i])
                  stop_mes(mes, file=fcerr)
               }
               param[i]=sum(ms[im])/sum(ss[im])
            }
         }
      } else {
         # if no optimization, set all scaling params to 1.
         param[nb_ff+seq_len(nb_sc_tot)]=1.
      }
   }
""")
    f.write("""
#browser()
   # see if there are any active inequalities at starting point
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   nbad=sum(abs(ineq)<=tol_ineq)
   if (anyNA(param))
      stop_mes("NA introduced in 'param' values during scaling factor estimation:\n\t",
      paste0(apply(cbind(nm_par, param), 1, paste0, collapse="\t"), collapse="\n\t"))
   if (nbad > 0) {
      cat("The following ", nbad, " ineqalitie(s) are active at starting point", runsuf, ":\\n",
         paste(names(ineq[abs(ineq)<=tol_ineq]), collapse="\\n"), "\\n", sep="", file=fclog)
   }
""")

    f.write("""
   if (TIMEIT && wkvh) {
      cat("kvh init: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
   }
""")
    # main part: call optimization
    f.write("""
   if (write_res && wkvh) {
      cat("influx\\n", file=fkvh)
      cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="")
      cat("\\tlabeling\\t", if (case_i) "instationary" else "stationary", "\\n", file=fkvh, sep="")
      # save options of command line
      cat("\\truntime options\\n", file=fkvh)
      cat("\\t\\t%s\\n", file=fkvh)
   }
   """%ropts_s)
    f.write("""
   if (write_res && wkvh) {
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
   }
#browser()
   if (!length(rres)) {
      rres <- lab_resid(param, cjac=FALSE, labargs)
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         #close(fkvh)
         retcode[irun]=rres$err
         next
      }
      if (sum(is.infinite(rres$res))) {
         cat("Infinite values appeared in residual vector (at starting point)", file=fcerr)
         retcode[irun]=1
         #close(fkvh)
         next
      }
   }
   rcost=if (length(rres$res) && !all(ina <- is.na(rres$res))) sum(crossprod(rres$res[!ina])) else NA
   if (write_res && wkvh) {
      obj2kvh(rcost, "starting cost value", fkvh, indent=1)
      obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1)
   }
   fg=numeric(nb_f$nb_fgr)
   names(fg)=nm_list$fgr
   if (nb_f$nb_fgr > 0) {
      fg[paste("g.n.", substring(nm_list$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm_list$poolf]
   }
   if (write_res && wkvh) {
      btmp=as.numeric(p2bfl%stm%param[seq_len(nb_f$nb_ff)]+bp+g2bfl%stm%fg)
      names(btmp)=dimnames(Afl)[[1]]
      obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1)
   }

   #cat("mass vector:\\n", file=fclog)
   #print_mass(x)

   names(param)=nm_par
""")
    f.write("""
#browser()
   if (optimize && nb_ff+nb_poolf > 0L) {
      if (!(least_norm || sln || !"nlsic" %in% methods)) {
         # check if at starting position all fluxes can be resolved
         if (TIMEIT) {
            cat("check ja: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         rres=lab_resid(param, cjac=TRUE, labargs)
         if (sum(is.infinite(rres$res))) {
            cat("Infinite values appeared in residual vector (at identifiability check)", file=fcerr)
            retcode[irun]=1
            #close(fkvh)
            next
         }
         if (any(is.infinite(rres$jacobian))) {
            cat("Infinite values appeared in Jacobian (at identifiability check)", file=fcerr)
            retcode[irun]=1
            #close(fkvh)
            next
         }
         qrj=qr(jx_f$dr_dff, LAPACK=TRUE)
         d=diag(qrj$qr)
         qrj$rank=sum(abs(d)>abs(d[1])*tol)
         if (is.na(qrj$rank)) {
            cat("Rank of starting jacobian could not be estimated.", file=fcerr)
            retcode[irun]=1
            #close(fkvh)
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
            cat(sprintf("Provided measurements (labeling and fluxes) are not sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:\\n%s\\nJacobian dr_dff is written in the result kvh file (if --wkvh is activated).\\n",
               paste(nm_uns, sep=", ", collapse=", ")),
               file=fcerr)
            if (write_res && wkvh) {
               obj2kvh(jx_f$dr_dff, "Jacobian dr_dff", fkvh, indent=0)
            }
            #close(fkvh)
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
#browser()
      for (method in methods) {
         if (write_res) {
            capture.output(res <- opt_wrapper(param, method, measurements, jx_f, labargs), file=fclog)
         } else {
            res <- opt_wrapper(param, method, measurements, jx_f, labargs)
         }
         if ((!is.null(res$err) && res$err) || is.null(res$par)) {
            cat("***Warning: error occured in first optimization pass", runsuf, ": ", res$mes, "\\n", sep="", file=fclog)
            res$par=rep(NA, length(param))
            res$cost=NA
         } else if (!is.null(res$mes) && nchar(res$mes)) {
            cat("***Warning: in first optimization pass in run ", runsuf, ": ", res$mes, "\\n", sep="", file=fclog)
         }
         if (anyNA(res$par)) {
#browser()
            res$retres$jx_f=NULL # to avoid writing of huge data
            if (write_res && wkvh) {
               obj2kvh(res, "failed first pass optimization process information", fkvh)
            }
            cat("Optimization failed", runsuf, ": ", res$mes, "\\n", file=fcerr, sep="")
            #close(fkvh) # some additional information can be written into fkvh
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
            nm_inv=names(which((ui%*%res$par-ci)[,1]<=tol_ineq))
            i=grep("^zc ", nm_inv, v=TRUE)
            if (length(i) > 0) {
               i=str2ind(i, nm_i)
               cat("The following inequalities are active after first pass
   of zero crossing strategy and will be inverted", runsuf, ":\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="", file=fclog)
               ipos=grep(">=", nm_i[i], v=TRUE)
               ineg=grep("<=", nm_i[i], v=TRUE)
               ui[i,]=-ui[i,,drop=FALSE]
               if (length(ipos)) {
                  ipzc=str2ind(ipos, nm_izc)
                  ipos=str2ind(ipos, nm_i)
                  ci[ipos]=as.numeric(zc+mi_zc[ipzc,,drop=FALSE]%*%mic)
                  nm_i[ipos]=sub(">=", "<=-", nm_i[ipos])
               }
               if (length(ineg)) {
                  inzc=str2ind(ineg, nm_izc)
                  ineg=str2ind(ineg, nm_i)
                  ci[ineg]=as.numeric(zc+mi_zc[inzc,,drop=FALSE]%*%mic)
                  nm_i[ineg]=sub("<=-", ">=", nm_i[ineg])
               }
               rownames(ui)=nm_i
               names(ci)=nm_i
               # enforce new inequalities
               reopt=TRUE
               if (write_res) {
                  capture.output(pinside <- put_inside(res$par, ui, ci, rcond=1e7), file=fclog)
               } else {
                  pinside <- put_inside(res$par, ui, ci, rcond=1.e7)
               }
               if (anyNA(pinside)) {
                  if (!is.null(attr(pinside, "err")) && attr(pinside, "err")!=0) {
                     # fatal (for zc) error occured, don't reoptimize
                     cat("***Warning. Infeasible inequalities after zerocross inversion.\n", file=fclog)
                     reopt=FALSE
                  }
               } else if (!is.null(attr(pinside, "err")) && attr(pinside, "err")==0) {
                  # non fatal problem
                  cat(paste("***Warning: put_inside", runsuf, ": ", attr(pinside, "mes"), "\\n", collapse=""), file=fclog)
               }
               # reoptimize
               if (reopt) {
                  cat("Second zero crossing pass", runsuf, "\\n", sep="", file=fclog)
                  if (write_res) {
                     capture.output(reso <- opt_wrapper(pinside, method, measurements, new.env(), labargs), file=fclog)
                  } else {
                     reso <- opt_wrapper(pinside, method, measurements, new.env(), labargs)
                  }
                  if (reso$err || is.null(reso$par)) {
                     cat("***Warning: error in second zero crossing pass: ", reso$mes, "\\n", sep="", file=fclog)
                  } else if (!is.null(reso$mes) && nchar(reso$mes)) {
                     cat("***Warning: second zero crossing pass", runsuf, ": ", reso$mes, "\\n", sep="", file=fclog)
                  }
                  if(!reso$err && !is.null(reso$par) && !any(is.na(reso$par))) {
                     param=reso$par
                     res=reso
                     jx_f=labargs$jx_f
                  }
                  if (any(is.na(reso$par))) {
                     reso$retres$jx_f=NULL # to avoid writing of huge data
                     if (write_res && wkvh)
                        obj2kvh(reso, "failed second pass optimization process information", fkvh)
                     cat("***Warning: second zero crossing pass failed. Keep free parameters from previous pass", runsuf, "\\n", file=fclog, sep="")
                  }
               }
               # last pass, free all zc constraints
               i=grep("^zc ", nm_i)
               if (length(i) > 0) {
                  if (TIMEIT) {
                     cat("last zc : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
                  }
                  ui=ui[-i,,drop=FALSE]
                  ci=ci[-i]
                  nm_i=nm_i[-i]
                  cat("Last zero crossing pass (free of zc constraints)", runsuf, "\\n", sep="", file=fclog)
                  if (write_res) {
                     capture.output(reso <- opt_wrapper(param, method, measurements, new.env(), labargs), file=fclog)
                  } else {
                     reso <- opt_wrapper(param, method, measurements, new.env(), labargs)
                  }
                  if (reso$err || is.null(reso$par) || (!is.null(res$mes) && nchar(res$mes))) {
                     cat("***Warning: last zero crossing (free of zc)", runsuf, ": ", reso$mes, "\\n", sep="", file=fclog)
                  }
                  if(!reso$err && !is.null(reso$par) && !any(is.na(reso$par))) {
                     param=reso$par
                     res=reso
                     jx_f=labargs$jx_f
                  }
                  if (any(is.na(res$par))) {
                     res$retres$jx_f=NULL # to avoid writing of huge data
                     if (write_res && wkvh)
                        obj2kvh(res, "failed last pass optimization process information", fkvh)
                     cat("***Warning: last zero crossing pass failed. Keep free parameters from previous passes", runsuf, "\\n", file=fclog, sep="")
                  }
               }
            } else {
               cat("After the first optimization, no zero crossing inequality was activated. So no reoptimization", runsuf, "\\n", sep="", file=fclog)
            }
         } # end if zero crossing
      } # for method
#browser()
      param=res$par
      names(param)=nm_par
      if (excl_outliers != F) {
         # detect outliers
         if (TIMEIT) {
            cat("outliers: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         iva=!is.na(res$res)
         zpval=rz.pval.bi(res$res)
         iout=which(zpval <= excl_outliers & iva)
         #cat("iout=", iout, "\\n", file=fclog)
         if (length(iout)) {
            measurements$outlier=iout
            outtab=cbind(residual=res$res[iout], `p-value`=zpval[iout])
            row.names(outtab)=nm_resid[iout]
            cat("Excluded outliers at p-value ", excl_outliers, ":\\n", sep="", file=fclog)
            write.table(outtab, file=fclog, append=TRUE, quote=FALSE, sep="\\t", col.names=FALSE)
            
            # optimize with the last method from methods
            if (write_res) {
               capture.output(reso <- opt_wrapper(param, tail(methods, 1L), measurements, new.env(), labargs), file=fclog)
            } else {
               reso <- opt_wrapper(param, tail(methods, 1L), measurements, new.env(), labargs)
            }
            if (reso$err || is.null(reso$par) || (!is.null(reso$mes) && nchar(reso$mes))) {
               cat("***Warning: error without outliers: ", reso$mes, "\\n", sep="", file=fclog)
            }
            if (any(is.na(reso$par))) {
               cat("***Warning: optimization with outliers excluded has failed", runsuf, "\\n", file=fclog, sep="")
               # continue without outlier exclusion
               measurements$outlier=NULL
            } else {
               res=reso
               param=reso$par
               names(param)=nm_par
               jx_f=labargs$jx_f
               labargs$measurements=measurements # store outliers
               if (write_res && wkvh)
                  obj2kvh(outtab, "excluded outliers", fkvh)
            }
         } else {
            cat("***Warning: outlier exclusion at p-value "%s+%excl_outliers%s+%" has been requested but no outlier was detected at this p-value threshold.", "\\n", sep="", file=fclog)
         }
      }
      if (case_i && time_order=="1,2") {
         if (TIMEIT) {
            cat("order 2 : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
         }
         labargs$time_order="2" # continue with the 2-nd order
         if (write_res) {
            capture.output(reso <- opt_wrapper(param, tail(methods, 1L), measurements, new.env(), labargs), file=fclog)
         } else {
            reso <- opt_wrapper(param, tail(methods, 1L), measurements, new.env(), labargs)
         }
         if (reso$err || is.null(reso$par) || (!is.null(reso$mes) && nchar(reso$mes))) {
            cat("***Warning: order2: ", reso$mes, "\\n", sep="", file=fclog)
         }
         if (any(is.na(reso$par))) {
            cat("***Warning: optimization time_order 2 (in '1,2' suite) has failed, run=", runsuf, "\\n", file=fclog, sep="")
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
      if (write_res && wkvh)
         obj2kvh(optinfo, "optimization process information", fkvh)
      rres=res$retres
   } else {
      if (TIMEIT) {
         cat("residjac: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
      }
      rres=lab_resid(param, TRUE, labargs)
   }
   if (TIMEIT) {
      cat("postopt : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
   }
   # active constraints
   if (!all(is.na(param))) {
      ine=as.numeric(abs(ui%*%param-ci))<tol_ineq
      if (any(ine) && write_res && wkvh) {
         obj2kvh(nm_i[ine], "active inequality constraints", fkvh)
      }
   }
   poolall[nm_poolf]=param[nm_poolf]

#browser()
   if (is.null(jx_f$jacobian)) {
      # final jacobian calculation
      if (write_res) {
         capture.output(rres <- lab_resid(param, cjac=TRUE, labargs), file=fclog)
      } else {
         rres <- lab_resid(param, cjac=TRUE, labargs)
      }
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         #if (write_res && wkvh)
         #   close(fkvh)
         retcode[irun]=rres$err
         next
      }
   }
   rcost=cumo_cost(param, labargs, rres)
   pres[,irun]=param
   costres[irun]=rcost
   if (write_res) {
      if (wkvh) obj2kvh(rcost, "final cost", fkvh)
#browser()
      # get z p-values on residual vector
      zpval=rz.pval.bi(rres$res)
      resid=list()
      if (sum(nb_meas)) {
         resid[["labeled data"]]=lapply(seq_len(nb_exp), function(iexp) if (is.matrix(jx_f$reslab[[iexp]])) jx_f$reslab[[iexp]] else cbind(residual=jx_f$reslab[[iexp]], `p-value`=zpval[seq_along(jx_f$reslab[[iexp]])]))
         names(resid[["labeled data"]])=nm_exp
      
         if (case_i) {
            resid[["labeled data p-value"]]=vector("list", nb_exp)
            names(resid[["labeled data p-value"]])=nm_exp
            for (iexp in seq_len(nb_exp)) {
               mtmp=zpval[seq_along(jx_f$reslab[[iexp]])]
               dim(mtmp)=dim(jx_f$reslab[[iexp]])
               dimnames(mtmp)=dimnames(jx_f$reslab[[iexp]])
               resid[["labeled data p-value"]][[iexp]]=mtmp
               rm(mtmp)
            }
         }
      }
      nb_reslab_tot=sum(sapply(jx_f$reslab, length))
      if (length(jx_f$resflu))
         resid[["measured fluxes"]]=cbind(residual=jx_f$resflu, `p-value`=zpval[nb_reslab_tot+seq_along(jx_f$resflu)])
      if (length(jx_f$respool))
         resid[["measured pools"]]=cbind(residual=if (is.matrix(jx_f$respool)) jx_f$respool[,1] else jx_f$respool, `p-value`=zpval[nb_reslab_tot+length(jx_f$resflu)+seq_along(jx_f$respool)])
      if (wkvh)
         obj2kvh(resid, "(simulated-measured)/sd_exp", fkvh)

      # simulated measurements -> out
      simul=list()
#browser()
      if (case_i) {
         if (sum(nb_meas)) {
            if (addnoise) {
               simul[["labeled data"]]=lapply(seq_len(nb_exp), function(iexp) jx_f$usm[[iexp]]+rnorm(length(jx_f$usm[[iexp]]))*measurements$dev$labeled[[iexp]])
               names(simul[["labeled data"]])=nm_exp
            } else {
               # move mass in usm into valid interval [0, 1] and sum=1
               simul[["labeled data"]]=lapply(seq_len(nb_exp), function(iexp) {
                  x=clamp(jx_f$usm[[iexp]], 0, 1)
                  # get unique mass names to sum up to 1
                  nmx=rownames(x)
                  nm_m=nmx[startsWith(nmx, "m:")]
                  if (length(nm_m)) {
                     # get unique fragments
                     fr_u=unique(sapply(strsplit(nm_m, ":", fixed=TRUE), function(v) paste0(c(v[1L:3L], ""), collapse=":")))
                     lapply(fr_u, function(nm) {
                        # get indexes per fragment
                        i=which(startsWith(nmx, nm))
                        mets=strsplit(nm, ":", fixed=TRUE)[[1L]][2L]
                        met1=strsplit(mets, "+", fixed=TRUE)[[1L]][1L]
                        if (length(i) < clen[met1]+1)
                           return(NULL)
                        s=colSums(x[i,,drop=FALSE])
                        x[i,] <<- arrApply::arrApply(x[i,,drop=FALSE], 2, "multv", v=1./s)
                        NULL
                     })
                  }
                  x
               })
               names(simul[["labeled data"]])=nm_exp
            }
            # simul --> .miso.sim
            mlp2LAB=c(m="MS", l="LAB", p="PEAK")
            cnm=c("Id", "Comment", "Specie", "Fragment", "Dataset", "Isospecies", "Value", "SD", "Time", "Residual", "Pvalue")
            for (fnm in nm_exp) {
               rnm=gsub("#", "", rownames(simul[["labeled data"]][[fnm]]), fixed=TRUE)
               mnm=strsplitlim(rnm, ":", fixed=TRUE, lim=NA, strict=TRUE)
               mnm=matrix(unlist(mnm), ncol=length(mnm[[1L]]), byrow=TRUE)
               ct=rep(colnames(simul[["labeled data"]][[fnm]]), each=nrow(simul[["labeled data"]][[fnm]]))
               c3=suppressWarnings(as.integer(mnm[, 3L]))
               df=cbind(
                  "",
                  "",
                  Specie=mnm[, 2L],
                  Fragment=ifelse(mnm[, 1L] == "m", mnm[, 3L], ""),
                  Dataset=paste0(mlp2LAB[mnm[, 1L]], "_", mnm[, 2L], "_", mnm[, 3L]),
                  Isospecies=ifelse(mnm[, 1L] == "m",
                     paste0("M", mnm[, 4L]), # MS: M0, M1, etc
                     ifelse(mnm[, 1L] == "l", mnm[, 3L], # label: 01x+00x etc
                     paste0(mnm[, 3L], "->", # peak: 2->1,3 etc.
                           ifelse(mnm[, 4L] == "S", "",
                           ifelse(mnm[, 4L] == "D-", c3-1,
                           ifelse(mnm[, 4L] == "D+", c3+1,
                           paste0(c3-1, ",", c3+1) # DD
                           )))))),
                  Value=c(simul[["labeled data"]][[fnm]]),
                  SD=measdev[[fnm]],
                  Time=ct,
                  Resid=c(resid[["labeled data"]][[fnm]]),
                  Pvalue=c(resid[["labeled data p-value"]][[fnm]])
               )
               colnames(df)=cnm
               write.table(df, sep="\t", quote=FALSE, row.names=FALSE, fileEncoding="utf8", file=file.path(dirres, paste0(fnm, runsuf, ".miso.sim")))
            }
         }
      } else {
         if (sum(nb_meas)) {
            if (addnoise) {
               simlab=lapply(seq_len(nb_exp), function(iexp) jx_f$simlab[[iexp]]+rnorm(length(jx_f$simlab[[iexp]]))*measurements$dev$labeled[[iexp]])
               names(simlab)=nm_exp
            } else {
               simlab=jx_f$simlab
               names(simlab)=nm_exp
            }
            if (nb_sc_tot > 0) {
               simul[["labeled data (unscaled)"]]=jx_f$usimlab
               simul[["labeled data (scaled)"]]=simlab
            } else {
               simul[["labeled data"]]=simlab
            }
            # simlab --> .miso.sim
            mlp2LAB=c(m="MS", l="LAB", p="PEAK")
            cnm=c("Id", "Comment", "Specie", "Fragment", "Dataset", "Isospecies", "Value", "SD", "Time", "Residual", "Pvalue")
            for (fnm in nm_exp) {
               rnm=gsub("#", "", names(simlab[[fnm]]), fixed=TRUE)
               mnm=strsplitlim(rnm, ":", fixed=TRUE, lim=NA, strict=TRUE)
               mnm=matrix(unlist(mnm), ncol=length(mnm[[1L]]), byrow=TRUE)
               c3=suppressWarnings(as.integer(mnm[, 3L]))
               df=cbind(
                  "",
                  "",
                  Specie=mnm[, 2L],
                  Fragment=ifelse(mnm[, 1L] == "m", mnm[, 3L], ""),
                  Dataset=paste0(mlp2LAB[mnm[, 1L]], "_", mnm[, 2L], "_", mnm[, 3L]),
                  Isospecies=ifelse(mnm[, 1L] == "m",
                     paste0("M", mnm[, 4L]), # MS: M0, M1, etc
                     ifelse(mnm[, 1L] == "l", mnm[, 3L], # label: 01x+00x etc
                     paste0(mnm[, 3L], "->", # peak: 2->1,3 etc.
                           ifelse(mnm[, 4L] == "S", "",
                           ifelse(mnm[, 4L] == "D-", c3-1,
                           ifelse(mnm[, 4L] == "D+", c3+1,
                           paste0(c3-1, ",", c3+1) # DD
                           )))))),
                  Value=simlab[[fnm]],
                  SD=measdev[[fnm]],
                  Time="",
                  resid[["labeled data"]][[fnm]]
               )
               colnames(df)=cnm
               write.table(df, sep="\t", quote=FALSE, row.names=FALSE, fileEncoding="utf8", file=file.path(dirres, paste0(fnm, runsuf, ".miso.sim")))
            }
         }
      }
      #browser()
      if (nb_fmn) {
         if (addnoise)
            simul[["measured fluxes"]]=jx_f$simfmn+rnorm(length(jx_f$simfm))*measurements$dev$flux
         else
            simul[["measured fluxes"]]=jx_f$simfmn
         # measured fluxes --> .mflux
         cnm=c("Id", "Comment", "Flux", "Value", "SD", "Residual", "Pvalue")
         df=structure(cbind("", "", substring(names(simul[["measured fluxes"]]), 5), simul[["measured fluxes"]], measurements$dev$flux, resid[["measured fluxes"]]), dimnames=list(NULL, cnm))
         write.table(df, sep="\t", quote=FALSE, row.names=FALSE, file=file.path(dirres, paste0(baseshort, runsuf, ".mflux.sim")))
      }
      if (nb_poolm) {
         if (addnoise)
            simul[["measured pools"]]=jx_f$simpool+rnorm(length(jx_f$simpool))*measurements$dev$pool
         else
            simul[["measured pools"]]=jx_f$simpool
         # measured metabolites --> .mmet
         cnm=c("Id", "Comment", "Specie", "Value", "SD", "Residual", "Pvalue")
         df=structure(cbind("", "", substring(names(simul[["measured pools"]]), 4), simul[["measured pools"]], measurements$dev$pool, resid[["measured pools"]]), dimnames=list(NULL, cnm))
         write.table(df, sep="\t", quote=FALSE, row.names=FALSE, file=file.path(dirres, paste0(baseshort, runsuf, ".mmet.sim")))
      }
      rm(resid, zpval)
      if (wkvh) obj2kvh(simul, "simulated measurements", fkvh)
   
      # SD -> out
      # get index of non null components
      iget=sapply(names(measurements$dev), function(nm) !is.null(measurements$dev[[nm]]) & nm %in% c("labeled", "flux", "pool"))
      if (wkvh) obj2kvh(measurements$dev[iget], "measurement SD", fkvh)

      # gradient -> kvh
      if (length(jx_f$res) && !all(ina <- is.na(jx_f$res)) && wkvh) {
         if (any(ina)) {
            gr=2*as.numeric(crossprod(jx_f$res[!ina], jx_f$jacobian[!ina,,drop=FALSE]))
         } else {
            gr=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
         }
         names(gr)=nm_par
         obj2kvh(gr, "gradient vector", fkvh)
      }
      colnames(jx_f$udr_dp)=nm_par
      if (wkvh) obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh)
      # generalized inverse of non reduced jacobian
      if (length(jx_f$udr_dp) > 0L && wkvh) {
         svj=svd(jx_f$udr_dp)
         svj$d=with(svj, ifelse(d >= d[1L]*tol, d, d[1L]*tol))
         invj=svj$v%*%(t(svj$u)/svj$d)
         dimnames(invj)=rev(dimnames(jx_f$udr_dp))
         obj2kvh(invj, "generalized inverse of jacobian dr_dp (without 1/sd_exp)", fkvh)
      }
   }

   labargs$getx=TRUE
   labargs$labargs2$getx=TRUE
   if (fullsys) {
   """)
    cu_keys=list(netan["cumo_input"][0].keys()) if netan["fullsys"] else []
    f.write("""
#browser()
      nm_xif=c(%s)
      # full label input is the same for all experiments
      xif=rep(list(c(%s)), nb_exp)"""%(join(", ", cu_keys, '"', '"'),
         join(", ", [(li[k] if li[k]==li[k] else "NA") for li in netan["cumo_input"] for k in cu_keys]),
    ))
    f.write("""
      xif=lapply(xif, setNames, nm_xif)
      if (case_i) {
         # prepare xif from funlab
         xilf=vector("list", nb_exp)
         xi2f=vector("list", nb_exp)
         for (iexp in seq(nb_exp)) {
            fli=funlabli[[iexp]]
            if (length(fli) == 0) {
               # replicate first column in xi as many times as there are time points
               if (time_order == "2" || time_order == "1,2")
                  xi2f[[iexp]]=matrix(xif[[iexp]], nrow=length(xif[[iexp]]), ncol=nb_tifu2[[iexp]])
               xilf[[iexp]]=matrix(xif[[iexp]], nrow=length(xif[[iexp]]), ncol=nb_tifu[[iexp]])
            }  else {
               # use funlab
               xilf[[iexp]]=funlab(tifull[[iexp]], nm_xif, fli, envfunlab, emu=FALSE, nm_exp[[iexp]], fcerr)
               if (time_order == "2" || time_order == "1,2")
                  xi2f[[iexp]]=funlab(tifull2[[iexp]], nm_xif, fli, envfunlab, emu=FALSE, nm_exp[[iexp]], fcerr)
            }
         }
#browser()
         xif=xilf
         if (time_order == "2" || time_order == "1,2")
            labargs$labargs2$xif=xi2f
      }
      labargs$nm_list$xif=nm_xif
      labargs$nm_list$inp=nm_xif
      labargs$xif=xif
      labargs$nb_f$xif=length(xif[[1L]])
      labargs$labargs2$nm_list$xif=nm_xif
      labargs$labargs2$nm_list$inp=nm_xif
      labargs$labargs2$nb_f$xif=length(xif[[1L]])
      if (case_i && !is.null(labargs$cl)) {
         clusterExport(labargs$cl, c("labargs"))
         clusterEvalQ(labargs$cl, {
            labargs$spa=sparse2spa(labargs$spa)
            labargs$labargs2$spa=labargs$spa
            labargs$spaf=sparse2spa(labargs$spaf)
            labargs$labargs2$spaf=labargs$spaf
#print(c(lab_new="", labargs=labargs, labargs2=labargs$labargs2, a1_1=labargs$spaf[[1]]$a, a1_2=labargs$labargs2spaf[[1]]$a))
         })
      }
#browser()
      v=lab_sim(param, cjac=FALSE, labargs, fullsys)
      if (identical(v$err, 1L)) {
         #save(v$fA$triplet(), file="singular237_triplet.RData")
         stop_mes("fullsys: weight=", v$iw, "; ", v$mes, file=fcerr)
      }
#browser()
      labargs$nm_list$inp=nm_inp
      if (time_order == "2" || time_order == "1,2") {
         labargs$labargs2$nm_list$inp=nm_inp
      }
      x=if (case_i) v$xf else v$x
   } else {
      v=lab_sim(param, cjac=FALSE, labargs)
      x=if (case_i) v$xf else v$x
   }

   mid=cumo2mass(x)
   if (case_i) {
      mid=lapply(mid, function(m) m[sort(rownames(m)),,drop=FALSE])
   } else if (length(mid)) {
      mid=mid[sort(rownames(mid)),,drop=FALSE]
   }
   # write some info in result kvh
   if (write_res && wkvh) {
      obj2kvh(mid, "MID vector", fkvh)
      
      # constrained fluxes to kvh
      obj2kvh(fallnx[nm_fc], "constrained net-xch01 fluxes", fkvh)
   }
   fwrv=v$lf$fwrv
   fallnx=v$lf$fallnx
   flnx=v$lf$flnx
   fgr=fallnx[nm_fgr]

   # keep last jx_f in jx_f_last
#browser()
   while (sensitive=="mc" && !all(is.na(param))) {
      if (TIMEIT) {
         cat("monte-ca: ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
      }
      if (set_seed) {
         set.seed(seed)
      }
      # reference simulation corresponding to the final param
      refsim=new.env()
      for (nm_it in c("simlab", "simfmn", "simpool", "usm")) {
         assign(nm_it, jx_f[[nm_it]], envir=refsim)
      }
      # Monte-Carlo simulation in parallel way (if asked and possible)
      if (np > 1L && !parR) {
         inmc=seq_len(nmc)
         #spli=splitIndices(nmc, nodes);
         spli=suppressWarnings(split(inmc, seq_len(np)))
         ii=inmc
         ii[unlist(spli)]=inmc # inverse the counting in spli
         clusterExport(cl, c("param", "refsim", "runsuf", "cl_worker", "mc_sim")) #, "spli"))
         ##clusterEvalQ(cl, labargs$spa[[1]]$a <- NULL) # to rebuild sparse matrices on cores ## now, they are build once, at the cluster init
         cl_res=clusterEvalQ(cl, {mc_iter=TRUE; labargs$getx=FALSE;})
         mc_res=parLapply(cl, spli, function(iv) lapply(iv, function(i) cl_worker(funth=mc_sim, argth=list(i))))
         mc_res=Reduce(c, mc_res, NULL)[ii]
         #mc_res=vector(nmc, mode="list")
         #for (i in seq(nodes))
         #   mc_res[spli[[i]]]=cl_res[[i]]
         #mc_res=parLapplyLB(cl, seq_len(nmc), function(imc) cl_worker(funth=mc_sim, argth=list(imc)))
      } else {
         mc_res=lapply(seq_len(nmc), function(imc) cl_worker(funth=mc_sim, argth=list(imc)))
      }
#browser()
      free_mc=sapply(mc_res, function(l) {if (class(l)=="character" || is.null(l) || is.na(l$cost) || l$err) { ret=rep(NA, nb_param+3) } else { ret=c(cost=l$cost, it=l$it, normp=l$normp, l$par) }; ret })
      if (length(free_mc)==0) {
         cat("***Warning: parallel exectution of Monte-Carlo simulations has failed.", "\\n", sep="", file=fclog)
         free_mc=matrix(NA, nb_param+2, 0)
      }
      cost_mc=free_mc[1,]
      nmc_real=nmc-sum(is.na(free_mc[4,]))
      # remove failed m-c iterations
      ifa=which(is.na(free_mc[1,]))
      if (length(ifa)) {
         if (ncol(free_mc) > length(ifa)) {
            cat("***Warning: some Monte-Carlo iterations failed.", "\\n", sep="", file=fclog)
         } else {
            cat("All Monter-Carlo samples failed.", "\\n", sep="", file=fcerr)
            retcode[irun]=1
            break
         }
         free_mc=free_mc[,-ifa,drop=FALSE]
         cost_mc=cost_mc[-ifa]
      }
      if (nmc_real <= 1) {
         cat("No sufficient Monter-Carlo samples were successfully calculated to do any statistics.", "\\n", sep="", file=fcerr)
         retcode[irun]=1
         break
      }
#browser()
      fallnx_mc=apply(free_mc[-(1L:3L),,drop=FALSE], 2, function(p)param2fl(p, labargs)$fallnx)
      if (length(fallnx_mc)) {
         # mean
         fallnx_mc_mean=rowMeans(fallnx_mc)
         # median
         fallnx_mc_median=apply(fallnx_mc, 1, median)
         # confidence intervals
         cinx_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)))
         cinx_mc=cbind(cinx_mc, "CI 95% length"=c(diff(t(cinx_mc))))
         cinx_mc=cbind(cinx_mc, "relative CI (%)"=cinx_mc[,3]*100/abs(fallnx_mc_mean))
      }

      if (write_res && wkvh) {
         cat("monte-carlo\\n", file=fkvh)
         indent=1
         obj2kvh(cl_type, "cluster type", fkvh, indent)
         obj2kvh(avaco, "detected cores", fkvh, indent)
         avaco=max(1, avaco, na.rm=TRUE)
         obj2kvh(min(avaco, np, na.rm=TRUE), "used cores", fkvh, indent)
         cat("\\tfitting samples\\n", file=fkvh)
         indent=2
         obj2kvh(nmc, "requested number", fkvh, indent)
         obj2kvh(nmc_real, "calculated number", fkvh, indent)
         obj2kvh(nmc-nmc_real, "failed to calculate", fkvh, indent)
         # convergence section in kvh
         indent=1
         mout=rbind(round(free_mc[1:2,,drop=FALSE], 2),
            format(free_mc[3,,drop=FALSE], di=2, sci=TRUE))
         dimnames(mout)=list(c("cost", "it.numb", "normp"), seq_len(ncol(free_mc)))
         obj2kvh(mout, "convergence per sample", fkvh, indent)
      }
      free_mc=free_mc[-(1L:3L),,drop=FALSE]
      rownames(free_mc)=nm_par

      # param stats
      # mean
      parmean=rowMeans(free_mc)
      # median
      parmed=apply(free_mc, 1, median)

      # confidence intervals
      ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))))
      colnames(ci_mc)=c("CI 2.5%", "CI 97.5%", "CI length")

      # cost section in kvh
      if (write_res && wkvh) {
         cat("\\tcost\\n", file=fkvh)
         indent=2
         obj2kvh(mean(cost_mc), "mean", fkvh, indent)
         obj2kvh(median(cost_mc), "median", fkvh, indent)
         obj2kvh(sd(cost_mc), "sd", fkvh, indent)
         obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent)
         obj2kvh(quantile(cost_mc, c(0.025, 0.95, 0.975)), "ci", fkvh, indent)
         
         # free parameters section in kvh
         cat("\\tStatistics\\n", file=fkvh)
         mout=c()
         indent=2
#browser()
         # covariance matrix
         covmc=cov(t(free_mc))
         obj2kvh(covmc, "covariance", fkvh, indent)
         # sd
         sdmc=sqrt(diag(covmc))
         mout=cbind(mout, mean=parmean, median=parmed, sd=sdmc,
            "rsd (%)"=sdmc*100/abs(parmean), ci_mc)
         obj2kvh(mout, "free parameters", fkvh, indent)

         # net-xch01 stats
         fallnx=param2fl(param, labargs)$fallnx
         if (length(fallnx_mc)) {
            dimnames(fallnx_mc)[[1]]=nm_fallnx
            # form a matrix output
            fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0)
            # covariance matrix
            covmc=cov(t(fallnx_mc))
            dimnames(covmc)=list(nm_fallnx, nm_fallnx)
            # sd
            sdmc=sqrt(diag(covmc))
            fallout=cbind(fallout, mean=fallnx_mc_mean, median=fallnx_mc_median, sd=sdmc,
               "rsd (%)"=sdmc*100/abs(fallnx_mc_mean), cinx_mc)
            o=order(nm_fallnx)
            obj2kvh(fallout[o,,drop=FALSE], "all net-xch01 fluxes", fkvh, indent)
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
            cif_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)))
            cif_mc=cbind(cif_mc, t(diff(t(cif_mc))))
            cif_mc=cbind(cif_mc, cif_mc[,3]*100/abs(fwrv))
            dimnames(cif_mc)[[2]]=c("CI 2.5%", "CI 97.5%", "CI 95% length", "relative CI (%)")
            fallout=cbind(fallout, mean=parmean, median=parmed, sd=sdmc,
               "rsd (%)"=sdmc*100/abs(parmean), cif_mc)
            o=order(nm_fwrv)
            obj2kvh(fallout[o,,drop=FALSE], "forward-reverse fluxes", fkvh, indent)
            obj2kvh(covmc[o,o], "covariance of forward-reverse fluxes", fkvh, indent)
         }
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
      if (write_res) {
         capture.output(rres <- lab_resid(param, cjac=TRUE, labargs), file=fclog)
      } else {
         rres <- lab_resid(param, cjac=TRUE, labargs)
      }
      if (!is.null(rres$err) && rres$err) {
         cat("lab_resid", runsuf, ": ", rres$mes, "\\n", file=fcerr, sep="")
         #if (write_res && wkvh)
         #   close(fkvh)
         retcode[irun]=rres$err
         next
      }
   } # else use the last calculated jacobian
#browser()
   # covariance matrix of free fluxes
   if (length(jx_f$jacobian) > 0L && !all(is.na(param))) {
      svj=svd(jx_f$jacobian)
      if (svj$d[1] == 0.) {
         i=rep(TRUE, length(svj$d))
      } else {
         i=svj$d/svj$d[1]<tol
         if (all(!i) && svj$d[1]<tol) {
            # we could not find very small d, take just the last
            i[length(i)]=TRUE
         }
      }
      ibad=apply(svj$v[, i, drop=FALSE], 2, which.contrib)
      ibad=unique(unlist(ibad))
      if (length(ibad) > 0) {
         cat(paste(if (nchar(runsuf)) runsuf%s+%": " else "", "***Warning: inverse of covariance matrix is numerically singular.\\nStatistically undefined parameter(s) seems to be:\\n",
            paste(sort(nm_par[ibad]), collapse="\\n"), "\\nFor a more complete list, see SD column in '.tvar.sim' result file.", sep=""), "\\n", sep="", file=fclog)
      }
      # "square root" of covariance matrix (to preserve numerical positive definitness)
      svj$d=with(svj, ifelse(d >= d[1L]*tol, d, d[1L]*tol))
      rtcov=(svj$u)%*%(t(svj$v)/svj$d)
      # standard deviations of free fluxes
      if (write_res) {
         if (wkvh) cat("linear stats\\n", file=fkvh)

         # sd free+dependent+growth net-xch01 fluxes
         nm_flfd=c(nm_ff, nm_fgr, nm_fl)
         if (nb_ff > 0 || nb_fgr > 0) {
            i=1:nb_param
            i=c(head(i, nb_ff), tail(i, nb_fgr))
            covfl=crossprod(rtcov[, i, drop=FALSE]%mmt%(rbind(diag(nb_ff+nb_fgr), dfl_dffg)%mrv%c(rep.int(1., nb_ff), fgr)))
            dimnames(covfl)=list(nm_flfd, nm_flfd)
            sdfl=sqrt(diag(covfl))
         } else {
            sdfl=rep(0., nb_fl)
            covfl=matrix(0., nb_fl, nb_fl)
         }
         fl=c(head(param, nb_ff), fgr, flnx)
         stats_nx=cbind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl))
         rownames(stats_nx)=nm_flfd
         if (wkvh) {
            o=order(nm_flfd)
            obj2kvh(stats_nx[o,,drop=FALSE], "net-xch01 fluxes (sorted by name)", fkvh, indent=1)
            obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1)
         }
         
         # flux, pool --> .tvar
#browser()
         rnm=grep("_gr$", nm_fallnx, invert=TRUE, value=TRUE)
         cnm=c("Id", "Comment", "Name", "Kind", "Type", "Value", "SD", "Struct_identif")
         if (sensitive == "mc") {
            cnm=c(cnm, "Low_mc", "Up_mc")
         }
         nx2suf=c(n="NET", x="XCH")
         fd2cap=c(f="F", d="D", c="C")

         mnm=matrix(unlist(strsplitlim(rnm, ".", fixed=TRUE, lim=3L)), ncol=3L, byrow=TRUE)
         o=natorder(mnm[, 3L])
         rnm=rnm[o]
         mnm=mnm[o,,drop=FALSE]
         vfl=fallnx[rnm]
         vfl=ifelse(mnm[,2L] == "x", clamp(vfl, 0, 1), vfl)
         vsd=sdfl[rnm]
         vsd[is.na(vsd)]=0.
         df=cbind("", "", mnm[, 3L], nx2suf[mnm[, 2L]], fd2cap[mnm[, 1L]], vfl, vsd, ifelse(vsd > 10000, "no", "yes"))
         if (sensitive == "mc") {
            df=cbind(df, cinx_mc[rnm, 1L:2L])
         }

         # sd of all fwd-rev
         if (nb_ff > 0 || nb_fgr > 0) {
            i=1:nb_param
            i=c(head(i, nb_ff), tail(i, nb_fgr))
            covf=crossprod(tcrossprod_simple_triplet_matrix(rtcov[,i, drop=FALSE], jx_f$df_dffp%mrv%c(rep.int(1., nb_ff), head(poolall[nm_poolf], nb_fgr))))
            dimnames(covf)=list(nm_fwrv, nm_fwrv)
            sdf=sqrt(diag(covf))
         } else {
            sdf=rep(0., length(fwrv))
         }
         if (wkvh) {
            mtmp=cbind(fwrv, sdf, sdf/abs(fwrv))
            dimnames(mtmp)[[2]]=c("value", "sd", "rsd")
            o=order(nm_fwrv)
            obj2kvh(mtmp[o,], "fwd-rev fluxes (sorted by name)", fkvh, indent=1)
            if (nb_ff > 0 || nb_fgr > 0) {
               obj2kvh(covf, "covariance fwd-rev fluxes", fkvh, indent=1)
            }
         }
         # pool -> kvh
         sdpf=poolall
         sdpf[]=0.

         if (nb_poolf > 0) {
            # covariance matrix of free pools
            # "square root" of covariance matrix (to preserve numerical positive definitness)
            poolall[nm_poolf]=param[nm_poolf]
            # cov poolf
            covpf=crossprod(rtcov[,nb_ff+nb_sc_tot+seq_len(nb_poolf), drop=FALSE])
            dimnames(covpf)=list(nm_poolf, nm_poolf)
            sdpf[nm_poolf]=sqrt(diag(covpf))
         }
         if (length(poolall) > 0) {
            if (wkvh) {
               mtmp=cbind("value"=poolall, "sd"=sdpf, "rsd"=sdpf/poolall)
               rownames(mtmp)=nm_poolall
               o=order(nm_poolall)
               obj2kvh(mtmp[o,,drop=FALSE], "metabolite pools (sorted by name)", fkvh, indent=1)
               if (nb_poolf > 0) {
                  o=order(nm_poolf)
                  obj2kvh(covpf[o, o], "covariance free pools", fkvh, indent=1)
               }
            }
            rnm=names(poolall)
            mnm=matrix(unlist(strsplitlim(rnm, ":", fixed=TRUE, lim=2L)), ncol=2L, byrow=TRUE)
            o=natorder(mnm[, 2L])
            rnm=rnm[o]
            mnm=mnm[o,,drop=FALSE]
            pfc2cap=c(pf="F", pc="C")
            vsd=sdpf[rnm]
            vsd[is.na(vsd)]=0.
            dfp=cbind("", "", mnm[, 2L], "METAB", pfc2cap[mnm[,1L]], poolall[rnm], vsd, ifelse(vsd >= 10000., "no", "yes"))
            if (sensitive == "mc") {
               mci=ci_mc[nm_poolf, 1L:2L]
               mci=rbind(mci, cbind(poolall[nm_poolc], poolall[nm_poolc]))
               dfp=cbind(dfp, mci[rnm,])
            }
            df=rbind(df, dfp)
         }
         colnames(df)=cnm
         write.table(df, sep="\t", quote=FALSE, row.names=FALSE, file=file.path(dirres, paste0(baseshort, runsuf, ".tvar.sim")))
      }
   }

   # chi2 test for goodness of fit
   # goodness of fit (chi2 test)
   if (length(jx_f$res)) {
      if (is.na(rcost)) {
         cat(sprintf("***Warning: chi2: Reduced cost value is NA. Chi2 test cannot be done.\\n"), sep="", file=fclog)
      } else {
         nvres=sum(!is.na(jx_f$res))
         if (nvres > nb_param) {
            chi2test=list("chi2 value"=rcost, "data points"=nvres,
               "fitted parameters"=nb_param, "degrees of freedom"=nvres-nb_param)
            chi2test$`chi2 reduced value`=chi2test$`chi2 value`/chi2test$`degrees of freedom`
            chi2test$`p-value, i.e. P(X^2<=value)`=pchisq(chi2test$`chi2 value`, df=chi2test$`degrees of freedom`)
            chi2test$conclusion=if (chi2test$`p-value, i.e. P(X^2<=value)` > 0.95) "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD" else "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD"
            if (write_res) {
               if (wkvh) obj2kvh(chi2test, "goodness of fit (chi2 test)", fkvh, indent=1)
               fstat=file(file.path(dirres, sprintf("%s%s.stat", baseshort,  runsuf)), "w")
               df=c(rcost, rcost/(nvres-nb_param), nvres, nb_param, nvres-nb_param, chi2test$`p-value, i.e. P(X^2<=value)`, chi2test$conclusion)
               names(df)=c("chi2_value", "chi2/df", "number_of_measurements", "number_of_parameters", "degrees_of_freedom", "p-value", "conclusion")
               write.table(df, sep="\t", quote=FALSE, row.names=TRUE, file=fstat, col.names=FALSE)
               close(fstat)
            }
         } else {
            cat(sprintf("***Warning: chi2: Measurement number %d is lower or equal to parameter number %d. Chi2 test cannot be done.\\n", nvres, nb_param), sep="", file=fclog)
         }
      }
   }
   if (prof) {
      Rprof(NULL)
   }
   if (write_res) {
      if (wkvh) close(fkvh)
      # write edge.netflux property
      fedge=file(file.path(dirres, "tmp", sprintf("edge.netflux.%s%s.attrs", baseshort,  runsuf)), "w")
      cat("netflux (class=Double)\\n", sep="", file=fedge)
      nm_edge=names(edge2fl)
      cat(paste(nm_edge, fallnx[edge2fl], sep=" = "), sep="\\n" , file=fedge)
      close(fedge)

      # write edge.xchflux property
      fedge=file(file.path(dirres, "tmp",  sprintf("edge.xchflux.%s%s.attrs", baseshort,  runsuf)), "w")
      flxch=paste(".x", substring(edge2fl, 4), sep="")
      ifl=charmatch(flxch, substring(names(fallnx), 2))
      cat("xchflux (class=Double)\\n", sep="", file=fedge)
      cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
      close(fedge)

      # write node.log2pool property
      if (length(poolall)> 0) {
         fnode=file(file.path(dirres, "tmp", sprintf("node.log2pool.%s%s.attrs", baseshort,  runsuf)), "w")
         cat("log2pool (class=Double)\\n", sep="", file=fnode)
         nm_node=substring(names(poolall), 4)
         cat(paste(nm_node, log2(poolall), sep=" = "), sep="\\n" , file=fnode)
         close(fnode)
      }
   }
}
#browser()
if (!is.null(cl)) {
   stopCluster(cl)
   labargs$cl=cl=NULL
}

if (write_res) {
#browser()
   pres=rbind(cost=costres, pres)
   fco=file(file.path(dirres, "tmp", sprintf("%s.pres.csv", baseshort)), open="w")
   cat("row_col\t", file=fco)
   write.table(file=fco, pres, row.names=TRUE, quote=FALSE, sep="\\t")
   close(fco)
   # make simlinks to the lowest cost results
   if (length(costres) > 1L) {
      nmbest=colnames(pres)[which.min(costres)]
      for (f in c(list.files(path=dirres, pattern=sprintf("\\\\.%s\\\\..*\\\\.sim$", nmbest), full.names=TRUE),
            list.files(path=dirres, pattern=sprintf("\\\\.%s\\\\.stat$", nmbest), full.names=TRUE),
            list.files(path=file.path(dirres, "tmp"), pattern=sprintf("\\\\.%s\\\\.attrs$", nmbest), full.names=TRUE))) {
         flink=sub(sprintf(".%s", nmbest), "", f, fixed=TRUE)
         file.copy(to=flink, from=f, overwrite = TRUE, copy.date = TRUE)
      }
   }
}
if (TIMEIT) {
   cat("rend    : ", format(Sys.time()), " cpu=", proc.time()[1], "\\n", sep="", file=fclog)
}
""")
    f.write("""
# source files from FTBL/posttreat_R
postlist=strsplit("%(postlist)s", " *; *")[[1]]
for (post in postlist) {
   fpostR=file.path(dirw, post)
   if (file.exists(fpostR) && !isTRUE(file.info(fpostR)$isdir)) {
#cat("post in dirw=", fpostR, "\\n", file=fclog)
      source(fpostR)
   } else {
      # not found in 'dirw', try 'dirr'
      fpostR=file.path(dirr, post)
      if (file.exists(fpostR) && !isTRUE(file.info(fpostR)$isdir)) {
#cat("post in dirr=", fpostR, "\\n", file=fclog)
         source(fpostR)
      } else {
         cat(sprintf("***Warning: posttreatment R file '%%s' does not exist in working directory neither in influx_si one. Ignored.\\n", post), file=fclog)
      }
   }
}
xgc=gc(verbose=FALSE) # to avoid the message "Error in (function (x)  : tentative d'appliquer un objet qui n'est pas une fonction"
if (write_res) {
   close(fclog)
   close(fcerr)
}
retcode=max(retcode)
if (!interactive() && format(parent.frame()) == format(.GlobalEnv))
   q("no", status=retcode)
"""%{
         "postlist": escape(netan["opt"].get("posttreat_R", ""), '\\"'),
         })

    f.close()
    # try to make output files just readable to avoid later casual edition
    try:
        os.chmod(n_R, stat.S_IREAD)
    except:
        pass
    return 0

if __name__ == "__main__" or __name__ == "influx_si.cli":
    sys.exit(main())
