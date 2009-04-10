#!/usr/bin/env python

# Check if a given set of free fluxes can be resolved
# by a given set of measure data. I.e. try to respond
# wether experimental measures add suffient number of
# new constraints on fluxes to solve them in a unique way.
# Because of a non linear nature of the problem, the response
# can be valid in a neighbourhood of the proposed value of
# free fluxes or even worth, it can be valid only in that
# point if the point is on the boundary of feaseable domain.

# 2008-12-18 sokol
# Copyright 2008, INRA

# imports
import sys;
import os;
import time;
import getopt;

#sys.path.append('/home/sokol/dev/python/pytools');
#sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;
import ftbl2code;

def propag_meas(metabm, measures, netan):
    """propag_meas(metabm, measures, netan)
    Propagate measured values from non convergent metabolite to
    closest upstream convergent metabolite. Measure matrix is modified
    appropreatly to reflect the replacement.
    Return a dictionary dpropag[(metabm,icumom)]=set((metabup,icumoup))
    Limitation: propagation is made only for cases of a unique source.
    So the condensations are not treated and are signaled as errors."""
    # Algorithm:
    # - construct the path from metabm to metabup
    # - for each transfer and each product cumomer
    #   find source cumomer and put it in propagation
    # - update the dirctionary dpropag
    #print "propag: start with", metabm;##
    dpropag={};
    # gather all measured cumomers indexes for metabm
    cicumos=set(i for meas in measures
        for row in measures[meas]["mat"]
        for i in row["coefs"]
        if row["scale"].split(";")[0] == metabm);
    #aff("cicumos", cicumos);##
    curmetab=metabm;
    while True:
        # influx tuples: list of (icumo,(in_cumo, fl, imetab))
        #print "curmetab="+curmetab;##
        t_infl=[(icumo,C13_ftbl.cumo_infl(netan, curmetab+":"+str(icumo)))
            for icumo in cicumos];
        #aff("cicumos", cicumos);##
        #print t_infl;##
        #aff("t", t_infl);##
        # are there condensations? I.e. for one influx there are several metabs+imetab
        fl2metabs=dict();
        for (icumo, in_tl) in t_infl:
            for (in_cumo,fl,imetab,iin_metab) in in_tl:
                fl2metabs[fl]=fl2metabs.get(fl,set());
                fl2metabs[fl].add((in_cumo.split(":")[0],iin_metab));
        for (fl,s) in fl2metabs.iteritems():
            if len(s) > 1:
                raise Exception("The flux "+fl+" is a condensation flux for metabolite "+
                    curmetab+": "+str(s));
        
        # how many in-fluxes for the current metab
        nb_infl=len(set(fl for (cumo,in_tl) in t_infl
            for (in_cumo,fl,imetab,iin_metab) in in_tl));
        if nb_infl > 1:
            # we are at the end of our journey
            # update measures
            #print "endmetab="+curmetab;##
            if curmetab != metabm:
                for meas in measures:
                    for (ir,row) in enumerate(measures[meas]["mat"]):
                        if row["scale"].split(";")[0] != metabm:
                            # this row is not about measured metab, skip it
                            #print "skipping row", row["scale"];##
                            continue;
                        #aff("old row", row);##
                        newrow=dict();
                        newrow["scale"]=(curmetab+";"+
                            row["scale"].split(";", 1)[1]+
                            "p");
                        newrow["coefs"]=dict();
                        for (k,v) in row["coefs"].iteritems():
                            # current set of cumomers for a given source cumomer
                            cs=dpropag[(metabm,k)];
                            nb_cs=float(len(cs));
                            for (cm,ci) in cs:
                                newrow["coefs"][ci]=v/nb_cs;
                        #aff("new row", newrow);##
                        measures[meas]["mat"][ir]=newrow;
            return dpropag;
        # update dpropag and cicumos
        newicumos=set();
        dpropag=dict();
        for (icumo,in_tl) in t_infl:
            for (in_cumo,fl,imetab,iin_metab) in in_tl:
                (in_metab,cicumo)=in_cumo.split(":",1);
                cicumo=int(cicumo);
                dpropag[(metabm,icumo)]=dpropag.get((metabm,icumo),set());
                dpropag[(metabm,icumo)].add((curmetab,cicumo));
                newicumos.add(cicumo);
        cicumos=newicumos;
        # update curmetab
        curmetab=in_metab;
    return dpropag;

# take arguments
#<--skip in interactive session
# get arguments
me=os.path.basename(sys.argv[0]);
def usage():
    sys.stderr.write("usage: "+me+
        """ [-h|--help|--DEBUG] network[.ftbl]
Check if a given set of free fluxes can be resolved
by a given set of measure data.

OPTIONS
-h, --help print this message and exit
--DEBUG enable some debuggin features and output (for advanced users)

PARAMETERS
network - the base of an ftbl file (network.ftbl)

OUTPUT
network_res.kvh:
 - feasibles free flux sets with their values;
 - added equations on fluxes resulting from measures;

NB
Base name of ftbl file ('org' in this example)
is used to create or silently overwrite the following files:
network.R
network.f
network.o
network.so
network.log
network.err
network_res.kvh
""");
try:
    opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "cost", "DEBUG"]);
except getopt.GetoptError, err:
    sys.stderr.write(str(err)+"\n");
    usage();
    sys.exit(1);
cost=False;
DEBUG=False;
for o,a in opts:
    if o in ("-h", "--help"):
        usage();
        sys.exit();
    elif o=="--cost":
        cost=True;
    elif o=="--DEBUG":
        DEBUG=True;
    else:
        assert False, "unhandled option";
#aff("args", args);##
if len(args) != 1:
    usage();
    exit(1);
org=args[0];
#-->
# set some python constants
# org="simple"
# org="PPP"
if DEBUG:
    import pdb;

ftbl2code.DEBUG=DEBUG;

if org[-5:]==".ftbl":
    org=org[:-5];
n_ftbl=org+".ftbl";
n_opt=org+".R";
n_fort=org+".f";
f_ftbl=open(n_ftbl, "r");
try:
    os.system("chmod u+w '%s'"%n_opt);
    os.system("chmod u+w '%s'"%n_fort);
except:
    pass;

f=open(n_opt, "w");
ff=open(n_fort, "w");

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();

# analyse network
# reload(C13_ftbl);

netan=C13_ftbl.ftbl_netan(ftbl);

# prepare measures if needed
if "measures" not in netan:
    measures=dict();
    for meas in ("label", "mass", "peak"):
        measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas);
    netan["measures"]=measures;
measures=netan["measures"];
# propagate measures for non convergent metabs
#aff("meas before propag", measures);
for metabm in set(netan["label_meas"].keys()+netan["mass_meas"].keys()+netan["peak_meas"].keys()):
    try:
        propag_meas(metabm, measures, netan);
    except Exception, e:
        print "Warning: "+str(e);
        print "Mesures ignored";
        pass;
#aff("meas after propag", measures);
# write initialization part of R code
# and get measure dictionnary
dm=ftbl2code.netan2Rinit(netan, org, f, ff);
ff.close();
# compile fortran code
os.system("R CMD SHLIB %s"%n_fort);

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
""");
f.write("""
# get (overdetermined) matrix for measured cumomers balance
# sum v_i*MM*e_i=vout*m
# sum v_i=vout
# which gives
# M*[(MM*e_i-MM*e_n)]*tau=M(m-MM*e_n)
# where the [...] matrix is composed of column vectors MM*e_i-MM*e_n,
# MM are appropriate mapping matrices e_i->m vectors and
# M is restiction of measure matrix to the m components

# KM is global fw-rv flux matrix from measured isotopomers
# KL is just a list collection of individual K's
KL=list();

# names of constrained fluxes
nm_fc=c(paste(nm_fcn, "fwd", sep="."), paste(nm_fcx, "rev", sep="."));
# names of isotopomer measured fluxes
nm_imf=c();

# get initial cumomer concentrations x
x=vr$x;

# init total number of added equations
nb_addeq=0;
""");
# measured metab set
o_mcumos=dm["o_mcumos"];
rAb=dm["rAb"];
mmetab=set(cumo.split(":")[0] for cumo in o_mcumos);
# constrained fwd+rev
c_fr=[f_+".fwd" for f_ in netan["vflux_constr"]["net"]]+\
    [f_+".rev" for f_ in netan["vflux_constr"]["xch"]];
for metab in mmetab:
    # for a given measured metab get measured cumovector (vmcumo),
    # vectors of differences of source cumomers (columns of scumo),
    # measure matrix reduction to these measured cumomers (mr)
    # to get finaly tau matrix mr*scumo_k
    mmeas=[(meas, i, row["coefs"]) for meas in dm["o_meas"]
        for (i,row) in enumerate(dm["measures"][meas]["mat"])
        if metab==row["scale"].split(";")[0]];
    vmcumo=[dm["measures"][meas]["vec"][i]
        for (meas, i, row) in mmeas];
    # cumulated sum of row numbers in measure matrixes
    ics=dict(zip(dm["o_meas"], cumsum(len(dm["measures"][meas]["vec"])
        for meas in dm["o_meas"])));
    # ordered list of involved in measures cumomers for the metab
    o_mmcumos=sorted(set(metab+":"+str(ic)
        for (meas, i, row) in mmeas for ic in row if ic));
    # reduction of balance matrix to the metab
    mrA=dict((cumo,row) for wA in rAb["A"]
        for (cumo, row) in wA.iteritems()
        if cumo in o_mmcumos);
    # reduction of rhs to the metab
    mrb=dict((cumo,row) for wb in rAb["b"]
        for (cumo, row) in wb.iteritems()
        if cumo in o_mmcumos);
    # alphabeticaly ordered a flux list
    o_fl=set(fl for row in mrA.values()
        for (cumo, lfl) in row.iteritems()
        for fl in lfl);
    o_fl.update(fl for row in mrb.values()
        for fl in row.keys());
    o_fl=sorted(o_fl);
    if len(o_fl)==1:
        sys.stderr.write(me+
            """: fix me
 Metab """+metab+""" is not convergent point in the network.
 Propagate its measures upstream till next conv. point.
""");
        continue;
    # main python part: flux matrix construction
    # For a given flux (from involved ones) run through all o_mmcumos and
    # gather all source cumomer contributions
    scumo=dict();
    for fl in o_fl:
        scumo[fl]=dict();
        for cumo in o_mmcumos:
            scumo[fl][cumo]=list();
    # gather contributions for this cumo
    # brought by the flux fl (constrained fluxes are skiped)
    # linear terms are just 1-term lists [cumo]
    # condensation terms are lists of product cumomers [cumo1, cumo2]
    # 1. gather linear terms
    for (cumo, row) in mrA.iteritems():
        for (in_cumo, lfl) in row.iteritems():
            if cumo == in_cumo:
                continue;
            for fl in lfl:
                scumo[fl][cumo].append([in_cumo]);
    # 2. gather product terms
    for (cumo, row) in mrb.iteritems():
        for (fl, dfl) in row.iteritems():
            for l in dfl.values():
                scumo[fl][cumo].append(l);
    #aff("s", scumo);##
    f.write("""
metab="%(metab)s";
# indexes of measure matrix rows corresponding to the metab
imrow=c(%(imrow)s);
# measure values
vmcumo=c(%(vmcumo)s);
# involved cumomers
# %(o_mmcumos)s
# involved fluxes
nm_f=c(%(o_fl)s);
nb_f=length(nm_f);
# cumo matrix
scumo=matrix(c(%(scumo)s), nrow=%(nr_s)d);

# reduced to the metab measure matrix
mr=matrix(c(%(mr)s), ncol=%(nc_mr)d, byrow=TRUE);

# scumo modified by measure matrix and weighted by sigma
# K=sigma*(mr*scumo-v*1^t)
K=sqrt(measinvvar[imrow])*(mr%%*%%scumo-vmcumo);
dimnames(K)=list(c(), nm_f);
cat("In measures on", metab, ",", nb_f,
   "involved fluxes:", paste(nm_f, collapse=", "));
fc=intersect(nm_f,nm_fc);
if (length(fc)) {
   cat(" (constrained:", paste(fc, collapse=", "), ")");
}
cat("\\n");
cat("nb of added equations=", qr(K)$rank-1, "\\n");
nb_addeq=nb_addeq+qr(K)$rank-1;

# collect K matrices
KL=append(KL, list(K));
# collect flux names
nm_imf=union(nm_imf, nm_f);
# 
""" % {
    "metab": metab,
    "imrow": join(", ", (ics[meas]+i+1 for (meas, i, row) in mmeas)),
    "vmcumo": join(", ", vmcumo),
    "o_mmcumos": join(", ", o_mmcumos),
    "o_fl": join(", ", o_fl, '"', '"'),
    "scumo": join(",\n", (
        join(", ", [
        join("+", [
        join("*", 
        trd(p, dm["rcumo2i"], "x[", "]"))
#        p, "x[", "]")
        for p in scumo[fl][cumo]])
        for cumo in o_mmcumos]+["1."])
        for fl in o_fl)),
    "nr_s": len(o_mmcumos)+1,
    "mr": join(",\n", (
        join(", ", [row.get(int(cumo.split(":")[1]), "0.")
        for cumo in o_mmcumos]+[row.get(0, "0.")])
        for (m,i,row) in mmeas)),
    "nc_mr": len(o_mmcumos)+1,
    });

f.write("""
# exclude row/columns with constrained fluxes
for (f in nm_fcn) {
   fwd=paste(f,"fwd",sep=".");
   nm_imf=nm_imf[fwd != nm_imf]
}
for (f in nm_fcx) {
   rv=paste(f,"rev",sep=".");
   nm_imf=nm_imf[rv != nm_imf];
}
# construct KM global measure matrix from individual K's
ncol_km=length(nm_imf);
KM=matrix(0,0,ncol=ncol_km);
dimnames(KM)=list(c(), nm_imf);
for (K in KL) {
   # 0 fill init for new addition
   KM=rbind(KM,matrix(0,NROW(K),ncol_km));
   # start, end rows for this addition
   iend=NROW(KM);
   ista=iend-NROW(K)+1;
   # fill non zero entries
   for (f in dimnames(K)[[2]]) {
      if (f %in% nm_imf) {
         KM[ista:iend,f]=K[,f];
      }
   }
}

# rank of isotope measured fluxes (imf) matrix
cat("imf matrix rows=", nrow(KM), "\\n");
rank=qr(KM)$rank;
cat("imf matrix rank=", rank, "\\n");
cat("imf added eq=", nb_addeq, "\\n");

# rank of directly measured fluxes
cat("directly measured net fluxes=", NROW(fmn), "\\n");

# rank of Afl matrix
qra=qr(Afl)
rank=qra$rank;
cat("Afl matrix rank=", rank, "\\n");


# combine [Afl;p2bfl], meas.net and KM matrix into total matrix KT=[[Afl;p2bfl];KM]^t

# all not constrained fluxes .fwd, .rev
# dep:n,x+free:n,x
""");
f.write("""
nm_notrev=c(%(notrev)s);
"""%{
    "notrev": join(", ", netan["notrev"], '"', '"'),
    });
f.write("""
nm_tf=c();
nb_fln=length(nm_fln);
nb_flx=length(nm_flx);
nb_ffn=length(nm_ffn);
nb_ffx=length(nm_ffx);
nb_fmn=length(nm_fmn);
if (nb_fln) {
   nm_tf=c(nm_tf, paste(nm_fln,"fwd",sep="."));
   nm_tmp=setdiff(nm_fln, nm_notrev);
   if (length(nm_tmp)) {
      nm_tf=c(nm_tf, paste(setdiff(nm_fln, nm_notrev),"rev",sep="."));
   }
}
if (nb_flx) {
   nm_tf=union(nm_tf, paste(nm_flx,"rev",sep="."));
}
if (nb_ffn) {
   nm_tf=c(nm_tf, paste(nm_ffn,"fwd",sep="."));
}
if (nb_ffx) {
   nm_tf=c(nm_tf, paste(nm_ffx,"rev",sep="."));
}
if (nb_fmn) {
   nm_tf=c(nm_tf, paste(nm_fmn,"fwd",sep="."));
}
# make flux names unique
nm_tf=union(nm_tf, NULL);

nb_tf=length(nm_tf);
nb_ra=nrow(Afl);
nb_rm=NROW(fmn);
nb_rk=nrow(KM);
nb_rt=nb_ra+nb_rm+nb_rk;
KT=matrix(0, nrow=nb_rt, ncol=nb_tf);
dimnames(KT)=list(c(), nm_tf);
# adding Afl
# .net part is translated in .fwd-.rev
for (i in 1:nb_fln) {
   if (!nb_fln) {
      break;
   }
   f=nm_fln[i];
   fwd=paste(f, "fwd", sep=".");
   rev=paste(f, "rev", sep=".");
   KT[1:nb_ra,fwd]=Afl[,i];
   if (rev %in% nm_tf) {
      KT[1:nb_ra,rev]=-Afl[,i];
   }
}
# .xch is replaced by .rev
for (i in 1:nb_flx) {
   if (!nb_flx) {
      break;
   }
   f=nm_flx[i];
   rev=paste(f, "rev", sep=".");
   KT[1:nb_ra,rev]=Afl[,nb_fln+i];
}
# adding p2bfl
# .net part is translated in -.fwd+.rev
for (i in 1:nb_ffn) {
   if (!nb_ffn) {
      break;
   }
   f=nm_ffn[i];
   fwd=paste(f, "fwd", sep=".");
   rev=paste(f, "rev", sep=".");
   KT[1:nb_ra,fwd]=-p2bfl[,i];
   if (rev %in% nm_tf) {
      KT[1:nb_ra,rev]=p2bfl[,i];
   }
}
# .xch is replaced by .rev
for (i in 1:nb_ffx) {
   if (!nb_ffx) {
      break;
   }
   f=nm_ffx[i];
   rev=paste(f, "rev", sep=".");
   if (! f %in% nm_fln) {
      KT[1:nb_ra,rev]=p2bfl[,nb_ffn+i];
   }
}

# adding directly measured fluxes
i=0;
for (f in nm_fmn) {
   i=i+1;
   fwd=paste(f, "fwd", sep=".");
   KT[nb_ra+i,fwd]=sqrt(invfmnvar[i]);
}

# adding KM
for (f in nm_imf) {
   KT[nb_ra+nb_rm+(1:nb_rk),f]=KM[,f];
}

# rank of total matrix
qrk=qr(KT);
rank=qrk$rank;
cat("total matrix rows=", nrow(KT), "\\n");
cat("total matrix rank=", rank, "\\n");
nb_indeq=nb_addeq+NROW(fmn)+qra$rank;
if (nb_indeq < nb_tf) {
   # find unsolved fluxes
   i=qrk$pivot[(nb_indeq+1):nb_tf];
   cat("Isotopomer measures are not sufficient to calculate all fluxes.",
      "\\nWe have ",
      nb_tf, " fluxes and only ", nb_indeq,
      " independent equations.\\n", file=stderr());
   cat("Seems to be unsolved fluxes (", nb_tf-nb_indeq, "):",
      paste(nm_tf[i], collapse=", "), "\\n", file=stderr());
} else {
   cat("All fluxes can be solved\\n");
}
""");
f.close();
ff.close();

# preserve from accidental editing
try:
    os.system("chmod a-w '%s'"%n_opt);
    os.system("chmod a-w '%s'"%n_fort);
except:
    pass;
