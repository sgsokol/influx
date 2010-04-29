#!/usr/bin/env python

# Check if a given set of free fluxes can be resolved
# by a given set of measure data. I.e. try to respond
# wether experimental measures add suffient number of
# new constraints on fluxes to solve them in a unique way.
# Because of a non linear nature of the problem, the response
# can be valid in a neighbourhood of the proposed value of
# free fluxes or even worth, it can be valid only in that
# point if the point is on the boundary of feaseable domain.

# usage: ff2strid.py [-h|--help|--DEBUG] network[.ftbl]

# 2009-12-10 sokol
# Copyright 2009, INRA

# imports
import sys;
import os;
import time;
import getopt;
import numpy as np;

#sys.path.append('/home/sokol/dev/python/pytools');
#sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;
import ftbl2code;

def propag_meas(cmetab, measures, netan):
    """propag_meas(cmetab, measures, netan)->set((meas,irow)) of rows to delete
    in measures[meas]["mat","vec","dev"].
    Propagate measured values from metabolite cmetab to
    the closest upstream metabolite in_metab (or to two in_matbs if
    cmetab is a product of condensation reaction).
    Rows in measure matrix for cmetab are replaced by
    row(s) for in_metab(s) to reflect the measure propagation.
    Make only one step up in the network
    so we use recursive calls to run through all the upstair path.
    NB: A call to propag_meas must be followed by deleting
    row from returned list (starting with highest row numbers).
    """
    # Algorithm:
    # - count all input fluxes. If >1 there is nothing todo, we
    #   are already at convergent point
    # - for each carbon transfer calculate new coefficient in measures
    #print("propag: cmetab=", cmetab);##
    to_del=set();
    if len(C13_ftbl.infl(cmetab, netan))>1:
        # noting to do, we are at convergent point
        return(to_del);
    # gather all measured cumomers indexes for cmetab
    # current i-cumomers set from all measures
    cicumos=set(i for meas in measures
        for (ir,row) in enumerate(measures[meas]["mat"])
        for i in row["coefs"]
        if row and row["metab"] == cmetab);

    # for each measured icumo get incoming tuple
    # influx tuples: list of (icumo,(in_cumo, fl, imetab, iin_metab))
    #print "curmetab="+curmetab;##
    d_infl=dict((icumo,C13_ftbl.cumo_infl(netan, cmetab+":"+str(icumo)))
        for icumo in cicumos);
    #aff("cicumos", cicumos);##
    #print d_infl;##
    #aff("d", d_infl);##
    
    # one or two in_metabs?
    tin_metabs=set((in_cumo.split(":")[0], iin_metab)
        for l in d_infl.values()
        for (in_cumo, fl, imetab, iin_metab) in l);
    tin_metabs=list(tin_metabs);
    nb_inm=len(tin_metabs);
    if nb_inm > 2:
        raise Exception("Should not be: metabolite "+cmetab+" has too many (in_metab,iin_metab): "+str(tin_metabs));
    if nb_inm < 1:
        raise Exception("Should not be: metabolite "+cmetab+" has not any (in_metab,iin_metab): "+str(tin_metabs));
    # update measures row by row
    # for each cumomer propagate its coefficient to in_cumo
    # (with collected co-factors in condensations, if any)
    for meas in measures:
        for (ir,row) in enumerate(measures[meas]["mat"]):
            if row["metab"] != cmetab:
                # this row is not about measured metab, skip it
                #print "skipping row", row["scale"];##
                continue;
            print("proceeding row="+str(row));
            newrow=dict();
            for (icumo,coef) in row["coefs"].iteritems():
                t_in=d_infl[icumo];
                if len(t_in) == 1:
                    # simple 1 to 1 propagation
                    if len(t_in) != 1:
                        raise Exception("It must not be: t_in has not exactly one element: "+str(t_in));
                    in_cumo, fl, imetab, iin_metab=t_in[0];
                    in_metab, iin_cumo=in_cumo.split(":");
                    iin_cumo=int(iin_cumo);
                    newrow[in_metab]=newrow.get(in_metab, {});
                    newrow[in_metab]["coefs"]=newrow[in_metab].get("coefs", {});
                    newrow[in_metab]["coefs"][iin_cumo]=newrow[in_metab]["coefs"].get(iin_cumo, 0.);
                    newrow[in_metab]["coefs"][iin_cumo]+=coef;
                    #print "1 t_in="+str(t_in);##
                elif len(t_in) == 2:
                    # if row has already cofactors we cannot yet handle this
                    if "cofact" in row and icumo in row["cofact"]:
                        raise Exception("Not yet implemented", "measure row to propagate has already 'cofactor' field.\n We cannot yet propagate this information.");
                    for (ifact, icofact) in ((0,1), (1,0)):
                        in_cumo, fl, imetab, iin_metab=t_in[ifact];
                        in_metab, iin_cumo=in_cumo.split(":");
                        iin_cumo=int(iin_cumo);
                        newrow[in_metab]=newrow.get(in_metab, {});
                        newrow[in_metab]["coefs"]=newrow[in_metab].get("coefs", {});
                        newrow[in_metab]["coefs"][iin_cumo]=newrow[in_metab]["coefs"].get(iin_cumo, 0.);
                        newrow[in_metab]["cofact"]=newrow[in_metab].get("cofact", {});
                        newrow[in_metab]["cofact"][iin_cumo]=newrow[in_metab]["cofact"].get(iin_cumo,[]);
                        newrow[in_metab]["cofact"][iin_cumo].append((coef,t_in[icofact][0]));
                    #print "2 t_in="+str(t_in);##
                else:
                    raise Exception("Must not be", "len(t_in) > 2 or < 1 (%d)"%len(t_in));
            for in_metab in newrow:
                newrow[in_metab]["scale"]=row["scale"];
                newrow[in_metab]["metab"]=in_metab;
                newrow[in_metab]["bcumos"]=["#"+setcharbit("x"*netan["Clen"][in_metab],"1",i) for i in newrow[in_metab]["coefs"].keys()];
                print("meas="+meas+"; ir="+str(ir)+"; inm="+in_metab+"; nr="+str(newrow[in_metab])+"\n or="+str(row));
                # replace the old row by the new one
                measures[meas]["mat"].append(newrow[in_metab]);
                measures[meas]["vec"].append(measures[meas]["vec"][ir]);
                measures[meas]["dev"].append(measures[meas]["dev"][ir]);
                to_del.add((meas, ir));
                #aff("to_del r", to_del);
    # try to propagate in_metabs
    for in_metab in set(t[0] for t in tin_metabs):
        to_del.update(propag_meas(in_metab, measures, netan));
    return(to_del);

def num_netan(netan):
    """Prepare numerical object for a given set of free fluxes
    return a dict with
        "Afl": Afl,
        "bfl": bfl,
        "dfc": dfc_val,
        "nx": nx_val,
        "fwrv_val": fwrv_val,
        "Abx": -> list of matrices, rhs, cumomer vector and dicts
            calculated on reduced and full cumo set (keys: A,b,x,cumo_val).
    """
    # prepare numpy Afl and bfl for dependent flux solving
    Afl=np.array(netan["Afl"]);
    # prepare dictionary with dependent, free and constraint fluxe values
    # (dependent are added later)
    dfc_val=dict(("f.n."+f, v) for (f,v) in netan["flux_free"]["net"].iteritems());
    dfc_val.update(("f.x."+f, v) for (f,v) in netan["flux_free"]["xch"].iteritems());
    dfc_val.update(("c.n."+f, v) for (f,v) in netan["flux_constr"]["net"].iteritems());
    dfc_val.update(("c.x."+f, v) for (f,v) in netan["flux_constr"]["xch"].iteritems());
    bfl=np.array( [(sum(dfc_val.get(f, 1.)*v
        for (f,v) in row.iteritems()) if row else 0.) for row in netan["bfl"]] );
    # solve Afl*d=bfl;
    d=np.linalg.solve(Afl, bfl);
    dfc_val.update(("d.n."+f, d[i]) for (i,f) in enumerate(netan["vflux"]["net"]));
    dfc_val.update(("d.x."+f, d[i+len(netan["vflux"]["net"])]) for (i,f) in enumerate(netan["vflux"]["xch"]));
    
    # cut dfc prefix -> net-xch fluxes
    nx_val=dict();
    for (f, val) in dfc_val.iteritems():
        (dfc, nx, reac) = f.split(".");
        nx_val[nx+"."+reac] = val;

    # transform net-xch to fwd-rev fluxes
    fwrv_val=dict();
    for (f, val) in nx_val.iteritems():
        (nx, reac) = f.split(".");
        fw="fwd."+reac;
        if fw in fwrv_val:
            continue;
        rv="rev."+reac;
        xch_val=nx_val["x."+reac];
        net_val=nx_val["n."+reac];
        xch_val=xch_val/(1.-xch_val);
        fwrv_val[fw]=xch_val-min(-net_val, 0.);
        fwrv_val[rv]=xch_val-min(net_val, 0.);
    
    # prepare list of matrices and rhs for reduced cumomer set: rA, rb
    nAbx={
        "reduced": {"A": [], "b": [], "x": [], "cumo_val": None},
        "full": {"A": [], "b": [], "x": [], "cumo_val": None}
    };
    for (Ab, rf) in ((C13_ftbl.rcumo_sys(netan), "reduced"), (netan["cumo_sys"], "full")):
        #print("calculating "+rf+" system...");##
        cumo_val=dict();
        nAbx[rf]["cumo_val"]=cumo_val;
        vcumos=netan["vcumo"] if rf=="full" else netan["vrcumo"];
        for (w, A) in enumerate(Ab["A"]):
            cumos=vcumos[w];
            cumo2i=dict((c,i) for (i,c) in enumerate(cumos));
            b=Ab["b"][w];
            #print(b);##
            n=len(cumos);
            nA=np.zeros((n,n));
            nAbx[rf]["A"].append(nA);
            nb=np.zeros((n,));
            nAbx[rf]["b"].append(nb);
            for ir in xrange(n):
                # get terms of A in the current row ir
                cumor=cumos[ir];
                for (cumo, fl) in A[cumor].iteritems():
                    ic=cumo2i[cumo];
                    nA[ir, ic]=(-1. if ir!=ic else 1.)*sum(fwrv_val[f] for f in fl);
                if cumor not in b:
                    continue;
                # get term of b for this row
                row=b[cumor];
                #print("w="+str(w)+"; ir="+str(ir)+"; row="+str(row));##
                nb[ir]=sum(fwrv_val[fl]*sum([np.prod([cumo_val.get(x, x)
                    for x in xl])
                    for (imetab,xl) in grap.iteritems()])
                    for (fl, grap) in row.iteritems()
                    );
            # solve the A*x=b for this weight w
            x=np.linalg.solve(nA, nb);
            nAbx[rf]["x"].append(x);
            cumo_val.update((cumo, x[i]) for (i, cumo) in enumerate(cumos));
    return({
        "Afl": Afl,
        "bfl": bfl,
        "dfc": dfc_val,
        "nx": nx_val,
        "fwrv_val": fwrv_val,
        "Abx": nAbx,
        "cumo_val": cumo_val,
        });

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
Base name of ftbl file ('network' in this example)
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
num=num_netan(netan);
# propagate measures for non convergent metabs
#aff("meas before propag", measures);##
for metabm in set(netan["label_meas"].keys()+netan["mass_meas"].keys()+netan["peak_meas"].keys()):
    try:
        to_del=propag_meas(metabm, measures, netan);
        aff("final to_del", to_del);##
    except Exception, e:
        print "Warning: "+str(e);
        print "Mesures ignored";
        raise;
        pass;
    to_del=list(to_del);
    to_del.sort(reverse=True);
    for (meas, ir) in to_del:
        del(measures[meas]["mat"][ir]);
        del(measures[meas]["vec"][ir]);
        del(measures[meas]["dev"][ir]);
#aff("meas after propag", measures);##
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
nm_fc=c(paste("fwd", substr(nm_fcn,5,1000), sep="."), paste("rev", substr(nm_fcx,5,1000), sep="."));
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
c_fr=["fwd."+f_ for f_ in netan["vflux_constr"]["net"]]+\
    ["rev."+f_ for f_ in netan["vflux_constr"]["xch"]];
for metab in mmetab:
    # for a given measured metab get measured cumovector (vmcumo),
    # vectors of differences of source cumomers (columns of scumo),
    # measure matrix reduction to these measured cumomers (mr)
    # to get finaly tau matrix mr*scumo_k
    mmeas=[(meas, i, row["coefs"], row.get("cofact", {}))
        for meas in dm["o_meas"]
        for (i,row) in enumerate(dm["measures"][meas]["mat"])
        if metab==row["metab"]];
    vmcumo=[dm["measures"][meas]["vec"][i]
        for (meas, i, row, cof) in mmeas];
    # cumulated sum of row numbers in measure matrixes
    ics=dict(zip(dm["o_meas"], cumsum(len(dm["measures"][meas]["vec"])
        for meas in dm["o_meas"])));
    # ordered list of involved in measures cumomers for the metab
    o_mmcumos=sorted(set(metab+":"+str(ic)
        for (meas, i, row, cof) in mmeas for ic in row if ic));
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
        raise Exception("Should not be", me+""": fix me
 Metab """+metab+""" is not convergent point in the network.
 At this point it should be propagate its measures upstream
 till the next convergence point.
""");
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
cat("nb of added equations by this metabolite=", qr(K)$rank, "\\n");
nb_addeq=nb_addeq+qr(K)$rank;

# collect K matrices
KL=append(KL, list(K));
# collect flux names
nm_imf=union(nm_imf, nm_f);
# 
""" % {
    "metab": metab,
    "imrow": join(", ", (ics[meas]+i+1 for (meas, i, row, cof) in mmeas)),
    "vmcumo": join(", ", vmcumo),
    "o_mmcumos": join(", ", o_mmcumos),
    "o_fl": join(", ", o_fl, '"', '"'),
    "scumo": join(",\n", (
        join(", ", [
        join("+", [
        join("*", 
        trd(p, dm["rcumo2i"], "x[", "]", a=None))
#        p, "x[", "]")
        for p in scumo[fl][cumo]])
        for cumo in o_mmcumos]+["1."])
        for fl in o_fl)),
    "nr_s": len(o_mmcumos)+1,
    "mr": join(",\n", (
        join(", ", [str(row.get(int(cumo.split(":")[1]), 0.))+
        ("+("+join("+", (str(t[0])+"*"+str(num["cumo_val"][t[1]])
        for t in cof[int(cumo.split(":")[1])]))+")"
        if cof and int(cumo.split(":")[1]) in cof else "")
        for cumo in o_mmcumos]+[row.get(0, "0.")])
        for (m,i,row,cof) in mmeas)),
    "nc_mr": len(o_mmcumos)+1,
    });

f.write("""
# exclude columns with constrained fluxes
for (f in nm_fcn) {
   fwd=paste("fwd", substr(f,5,1000), sep=".");
   nm_imf=nm_imf[fwd != nm_imf]
}
for (f in nm_fcx) {
   rv=paste("rev", substr(f,5,1000), sep=".");
   nm_imf=nm_imf[rv != nm_imf];
}
# exclude columns with mesured fluxes
for (f in nm_fmn) {
   fwd=paste("fwd", substr(f,5,1000), sep=".");
   nm_imf=nm_imf[fwd != nm_imf]
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
i_rank=qr(KM)$rank;
cat("imf matrix rank=", i_rank, "\\n");

# rank of directly measured fluxes
m_rank=NROW(fmn);
cat("directly measured net fluxes=", m_rank, "\\n");

# rank of Afl matrix
qra=qr(Afl);
a_rank=qra$rank;
cat("Afl matrix rank=", a_rank, "\\n");


# combine [Afl;p2bfl], meas.net and KM matrix into total matrix KT=[[Afl;p2bfl];KM]^t

# all not constrained fluxes fwd., rev.
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
   nm_tf=c(nm_tf, paste("fwd",substr(nm_fln,5,1000),sep="."));
#   nm_tmp=setdiff(substr(nm_fln,5,1000), nm_notrev);
#   if (length(nm_tmp)) {
#      nm_tf=c(nm_tf, paste("rev",setdiff(substr(nm_fln,5,1000), nm_notrev),sep="."));
#   }
}
if (nb_flx) {
   nm_tf=c(nm_tf, paste("rev",substr(nm_flx,5,1000),sep="."));
}
if (nb_ffn) {
   nm_tf=c(nm_tf, paste("fwd",substr(nm_ffn,5,1000),sep="."));
}
if (nb_ffx) {
   nm_tf=c(nm_tf, paste("rev",substr(nm_ffx,5,1000),sep="."));
}
if (nb_fmn) {
   nm_tf=c(nm_tf, paste("fwd",substr(nm_fmn,5,1000),sep="."));
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
# .net part is translated in diff (fwd)-(rev)
for (i in 1:nb_fln) {
   if (!nb_fln) {
      break;
   }
   f=substr(nm_fln[i],5,1000);
   fwd=paste("fwd", f, sep=".");
   rev=paste("rev", f, sep=".");
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
   f=substr(nm_flx[i],5,1000);
   rev=paste("rev", f, sep=".");
   KT[1:nb_ra,rev]=Afl[,nb_fln+i];
}
# adding p2bfl
# .net part is translated in -fwd. + rev.
if (nb_ffn) {
   for (i in 1:nb_ffn) {
      f=substr(nm_ffn[i],5,1000);
      fwd=paste("fwd", f, sep=".");
      rev=paste("rev", f, sep=".");
      KT[1:nb_ra,fwd]=-p2bfl[,i];
      if (rev %in% nm_tf) {
         KT[1:nb_ra,rev]=p2bfl[,i];
      }
   }
}
# .xch is replaced by .rev
if (nb_ffx) {
   for (i in 1:nb_ffx) {
      f=substr(nm_ffx[i],5,1000);
      rev=paste("rev", f, sep=".");
      if (! f %in% nm_fln) {
         KT[1:nb_ra,rev]=p2bfl[,nb_ffn+i];
      }
   }
}

# adding directly measured fluxes
i=0;
for (f in nm_fmn) {
   i=i+1;
   fwd=paste("fwd", substr(f,5,1000), sep=".");
   KT[nb_ra+i,fwd]=sqrt(invfmnvar[i]);
}

# adding KM
for (f in nm_imf) {
   KT[nb_ra+nb_rm+(1:nb_rk),f]=KM[,f];
}

cat("It is necessary but not sufficient that
(Afl rank)+(imf rank)+(nb of measured fluxes) >= nb of fluxes\\n");
iq_sign=">=";
ok_ko="OK";
if (a_rank+i_rank+m_rank < nb_tf) {
   iq_sign="<";
   ok_ko="Bad";
}
cat(a_rank, "+", i_rank, "+", m_rank, "=", a_rank+i_rank+m_rank, iq_sign, nb_tf, " : ", ok_ko, "\\n", sep="");

# rank of total matrix
qrk=qr(KT);
t_rank=qrk$rank;
cat("total matrix rows=", nrow(KT), "\\n");
cat("total matrix columns=", ncol(KT), "\\n");
cat("total matrix rank=", t_rank, "\\n");
nb_indeq=nb_addeq+NROW(fmn)+qra$rank;
if (t_rank < nb_tf) {
   # find unsolved fluxes
   i=qrk$pivot[(t_rank+1):nb_tf];
   cat("Isotopomer measures are not sufficient to calculate all fluxes.",
      "\\nWe have ",
      nb_tf, " fluxes and only ", t_rank,
      " independent equations.\\n", file=stderr(), sep="");
   cat("Seems to be unsolved fluxes (", nb_tf-t_rank, "):\\n",
      paste(nm_tf[i], collapse="\\n"), "\\n", file=stderr(), sep="");
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
