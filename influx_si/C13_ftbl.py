"""- Parse .ftbl
- Analyse ftbl

Restrictions:
 - metabolite name cannot have
   
   ":"
     it's a separator in measure id
   "+"
     in measurements it can be metab1+metab2+...
  
"""
            # # normalize or not?
            # collab=sorted(v for li in labs for v in li) # will be collapsed labels.
            # # the group is normalizable if collapsed labs is composed of only "x"
            # while True:
                # collab=sorted(collab) # if there are two collapsible labs, they will be neighbors
                # found=False
                # for i in range(len(collab)):
                    # if i == len(collab)-1:
                        # break
                    # s=np.array(list(collab[i])) if i == 0 else n # string
                    # n=np.array(list(collab[i+1])) # next
                    # cdif=s != n
                    # if cdif.sum() == 1:
                        # idif=np.where(cdif)[0][0]
                        # found=s[idif] == "0" and n[idif] == "1"
                        # if found:
                            # s[idif]="x"
                            # collab[i]="".join(s)
                            # del(collab[i+1])
                            # break
                # if not found or len(collab) == 1:
                    # break
            # if len(collab) == 1 and collab[0] == "x"*len(collab) and all(val == val):
                # s=sum(val)
                # val=val/s
                # sdv=sdv/s
                # norma=True
            # else:
                # norma=False

# 2008-01-22 sokol: ftbl_parse(f)
# 2008-01-25 sokol: ftbl_netan(ftbl)
# 2008-01-30 sokol: sto_r_m is added to the result of ftbl_netan()
# 2008-03-05 sokol: ftbl_parse(): added labeled atoms transitions in "TRANS" key
# 2008-03-05 sokol: ftbl_netan(): added isotopomers balance equations
# 2008-03-06 sokol: ftbl_netan(): added fwd-rev flux balance matrix
# 2008-03-10 sokol: ftbl_netan(): added carbon length
# 2008-03-?? sokol: ftbl_netan(): added chemical formula (reactions)
# 2008-03-?? sokol: ftbl_netan(): added cumomer balance pattern (cumo_balance_pattern)
# 2008-03-19 sokol: enum_metab(): metabolite list enumerated along to chemical pathways
# 2008-03-20 sokol: ftbl_netan(): added cumomer systems (cumo_sys)
# 2008-03-21 sokol: ftbl_netan(): added carbon transitions (carbotrans)
# 2008-03-25 sokol: iw(): added "index in given weight" iterator
# 2008-03-26 sokol: ftbl_netan(): added Cmax, maximal carbon length
# 2008-03-26 sokol: cumo_path(): enumerate cumomers of a given weight along pathways
# 2008-03-26 sokol: ftbl_netan(): removed metab_paths ordering
# 2008-03-26 sokol: ftbl_netan(): removed cumo_w ordering
# 2008-03-26 sokol: ftbl_netan(): added set of all metabolites (metabs)
# 2008-04-14 sokol: ftbl_netan(): added set of free fluxes (flux_free[net|xch])
# 2008-04-14 sokol: ftbl_netan(): added set of constrained fluxes flux_constr[net|xch])
# 2008-04-16 sokol: ftbl_netan(): added set of measured fluxes (flux_measured)
# 2008-04-18 sokol: ftbl_netan(): changed structure of sto_m_r, sto_r_m, flux_m_r, cumo_balance_pattern and matrix A
# 2008-04-18 sokol: ftbl_netan(): removed cumo_balance_pattern
# 2008-04-18 sokol: ftbl_netan(): cumomers matrix construction rewritten from scratch
# 2008-04-21 sokol: ftbl_netan(): added input cumomers to matrix
# 2008-04-21 sokol: ftbl_netan(): added input cumomers to netan (cumo_input)
# 2008-04-21 sokol: iso2cumo(): calculate cumo fraction from isotopomer ones
# 2008-06-12 sokol: ftbl_netan(): added inequality analysis
# 2008-06-12 sokol: ftbl_netan(): added label measurements analysis
# 2008-06-23 sokol: ftbl_netan(): added peak measurements analysis
# 2008-06-23 sokol: ftbl_netan(): added mass measurements analysis
# 2008-06-26 sokol: ftbl_netan(): added matrices for label (H1) and mass
# 2008-06-27 sokol: ftbl_netan(): added matrix for peaks (C13)
# 2008-07-18 sokol: ftbl_netan(): added cumomers ordered lists (one by weight)
# 2008-07-22 sokol: ftbl_netan(): added net and xch equalities
# 2008-07-23 sokol: ftbl_netan(): added Afl, bfl
# 2008-07-23 sokol: ftbl_netan(): added fluxes ordered lists (net and xch)
# 2008-07-24 sokol: ftbl_netan(): added free fluxes ordered lists
# 2008-07-24 sokol: ftbl_netan(): added constrained fluxes ordered lists
# 2008-07-24 sokol: ftbl_netan(): added measured fluxes ordered lists
# 2008-07-24 sokol: ftbl_netan(): added compleet fluxes ordered lists
# 2008-07-28 sokol: ftbl_netan(): added Afl row names (vrowAfl)
# 2008-07-28 sokol: ftbl_netan(): added in-out xch fluxes to constrained flux list
# 2008-07-30 sokol: ftbl_netan(): added in-out fluxes (flux_in, flux_out)
# 2008-09-01 sokol: labprods(): get labeled product by a given metabolite in a given reaction
# 2008-09-03 sokol: allprods(): get all labeled product by a given metabolite with others labeled isotops in a given reaction
# 2008-09-23 sokol: added mat2graph()
# 2008-10-10 sokol: added rcumo_sys()
# 2009-02-02 sokol: added iin_metab to output of cumo_infl()
# 2009-03-24 sokol: ftbl_netan(): added tmax, dt, metab_scale, met_pools
# 2009-03-27 sokol: ftbl_parse(): added float_conv field set
# 2009-05-28 sokol: added t_iso2m(): transition matrix from isotopomer vector to MID vector
# 2009-05-28 sokol: added t_iso2cumo(): transition matrix from isotopomer vector to cumomer vector
# 2009-05-28 sokol: added t_iso2pos(): transition matrix from isotopomer vector to positional labelling vector
# 2009-07-21 sokol: added conv_mid(): convolution of two mass isotopomer distribution
# 2009-09-14 sokol: Flux names changed flux.[net|xch] -> [dfc].[nx].flux
# 2009-10-19 sokol: ftbl_netan(): added consistancy check on fluxes, metabs and lengthes
# 2009-11-26 sokol: added infl(): set of incoming fluxes
# 2010-02-15 sokol: ftbl_parse(): fixed input/output fluxes with respect to d/f/c characterisation
# 2010-05-05 sokol: ftbl_parse(): if file name is given, check for .ftbl and add it if needed
# 2010-05-31 sokol: ftbl_parse(): non blocking errors for not complete ftbl (e.g. for graph plotting)
# 2010-10-11 sokol: ftbl_parse(): for measurements added row number for later identification
# 2010-10-11 sokol: ftbl_netan(): added "id" field in measures
# 2010-11-23 sokol: ftbl_netan(): added "nx2dfc" dictionary
# 2010-12-02 sokol: ftbl_netan(): added growth flux option
# 2011-04-20 sokol: cumo_infl(),...: not reversible fluxes are replaced by flux_out in cumomer balance
# 2011-04-27 sokol: ftbl_neta(): added flux_inout
# 2011-04-27 sokol: cumo_infl(),...: not reversible fluxes are replaced by flux_inout in cumomer balance
# 2012-05-04 sokol: added pooled metabs (same carbon length) in measurements
# 2012-06-25 sokol: added EMU fragment gathering ms_frag_gath()
# 2012-11-08 sokol: added variable growth fluxes depending on variable
#                   concentrations (flux_vgrowth, vflux_growth.[net])
# 2012-11-23 sokol: added concentration measurements
# 2014-01-22 sokol: added possibility of unknown fluxes in EQUALITY section
import numpy as np
import re
import copy
import os
import sys
from codecs import BOM_UTF8, BOM_UTF16_BE, BOM_UTF16_LE, BOM_UTF32_BE, BOM_UTF32_LE
class oset(dict):
    def __init__(*args, **kwds):
        self, *args = args
        if len(args) == 1:
            for k in args[0]:
                self[k]=1
        elif len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
    def copy(self):
        tmp=oset()
        tmp.update(self)
        return(tmp)
    def add(self, x):
        self[x]=1
    def difference_update(self, x):
        for k in x:
            if k in self:
                del(self[k])
    def difference(self, x):
        return(self - oset(x))
    def update(self, x):
        if x:
            for i in x:
                self[i]=1
    def intersection(self, x):
        return self & oset(x)
    def __sub__(self, x):
        tmp=oset(x)
        return oset(i for i in self if i not in tmp)
    def __and__(self, x):
        tmp=oset(x)
        return oset(i for i in self if i in tmp)
    def __or__(self, x):
        tmp=oset(x)
        tmp.update(self)
        return tmp

BOMS = (
    (BOM_UTF8, "UTF-8"),
    (BOM_UTF32_BE, "UTF-32-BE"),
    (BOM_UTF32_LE, "UTF-32-LE"),
    (BOM_UTF16_BE, "UTF-16-BE"),
    (BOM_UTF16_LE, "UTF-16-LE"),
)

werr=sys.stderr.write
wout=sys.stdout.write
me=os.path.realpath(sys.argv[0])
dirx=os.path.dirname(me)
sys.path.append(dirx)
#sys.tracebacklimit=0

from tools_ssg import *
NaN=float("nan")
NA=NaN
tol=sys.float_info.epsilon*2**7
float_conv=oset((
    "VALUE",
    "DEVIATION",
    "OPT_VALUE",
    "META_SIZE",
    "VALUE_S",
    "VALUE_D-",
    "VALUE_D+",
    "VALUE_DD",
    "VALUE(F/C)",
))
# define legal sections/subsections
defsec={
    "PROJECT": oset(),
    "NETWORK": oset(),
    "NOTRACER_NETWORK": oset(("FLUX_NAME", "EQUATION")),
    "FLUXES": oset(("NET", "XCH")),
    "EQUALITIES": oset(("NET", "XCH", "METAB")),
    "INEQUALITIES": oset(("NET", "XCH", "METAB")),
    "LABEL_INPUT": oset(),
    "FLUX_MEASUREMENTS": oset(),
    "LABEL_MEASUREMENTS": oset(),
    "PEAK_MEASUREMENTS": oset(),
    "MASS_SPECTROMETRY": oset(),
    "METAB_MEASUREMENTS": oset(),
    "OPTIONS": oset(),
    "METABOLITE_POOLS": oset(),
}
# at least one of these fields required to be in prl_exp ftbls
req_prl=(
    "LABEL_INPUT",
    "LABEL_MEASUREMENTS",
    "PEAK_MEASUREMENTS",
    "MASS_SPECTROMETRY"
)
if "ffguess" not in locals():
    ffguess=False
def ftbl_parse(f):
    """ftbl_parse(f) -> dict
    read and parse .ftbl file. The only input parameter f is a stream pointer
    with read permission or a file name.
    This function parses the input and returns a dictionnary
    with items corresponding to sections in .ftbl. One section is added.
    "TRANS" correponds to carbon transitions."""
    import re
    import codecs
    ftbl=dict();    # main dictionary to be returned
    
    #print("f=", f)
    if not isstr(f):
        Exception("parameter 'f' must be a string with FTBL file name")
    if f[-5:].lower() != ".ftbl":
        f=f+".ftbl"
    ftbl["name"]=f
    ftbl["base_name"]=os.path.basename(f)[:-5]
    ftbl["abs_path"]=os.path.abspath(f)
    fc=open(f, "rb")
    raw=fc.read()
    fc.close()
    co=[encoding for bom, encoding in BOMS if raw.startswith(bom)]
    if len(co):
        inp=raw.decode(co[0])
    else:
        for co in ["utf-8", "latin9", "utf-16", "utf-32"]:
            try:
                inp="".join(c for c in raw.decode(co) if c != '\x00')
                break;
            except UnicodeDecodeError:
                pass
    #import pdb; pdb.set_trace()
    inp=inp.encode("utf-8").decode("utf-8-sig")
    lines=inp.splitlines()
    # fc=codecs.open(f, "r", encoding="utf-32")
    # try:
        # lines=fc.readlines()
        # fc.close()
    # except:
        # fc.close()
        # try:
            # fc=codecs.open(f, "r", encoding="utf-16")
            # lines=fc.readlines()
            # fc.close()
        # except:
            # fc.close()
            # try:
                # fc=codecs.open(f, "r", encoding="utf-8-sig")
                # lines=fc.readlines()
                # fc.close()
            # except:
                # fc.close()
                # fc=open(f, "r")
                # lines=fc.readlines()
                # fc.close()
    ftbl["pathway"]=dict()
    
    #print f;##
    reblank=re.compile("^[\t ]*$")
    recomm=re.compile("^[\t ]*//.*$")
    repath=re.compile("^[\t ]*//##[\t ]*(.*)[\t ]*$")
    comm=re.compile("^([^(//)]+|.+)//.*$")
    reading="sec_name"
    col_names=[]
    sec_name=subsec_name=""
    irow=0
    dic=dict()
    pathway=""
    #import pdb; pdb.set_trace()
    for l in lines:
        irow+=1
        #print "raw l="+l;##
        # strip out \r
        l=l.replace("\r", "")
        if len(l) == 0:
            continue
        #print "-ctrl-r l="+l;##
        # strip out double quots
        if l[0] == '"': l=l[1:];    # very begining
        #print "-fq l="+l;##
        if l[-1] == '"': l=l[:-1];    # very end
        #print "-lq l="+l;##
        l=l.replace("\t\"", "\t");    # at the field begining
        #print "-ffq l="+l;##
        l=l.replace("\"\t", "\t");    # at the end
        #print "-lfq l="+l;##
        
        # check for pathway name "//## pathname"
        if repath.match(l):
            pathway=repath.sub(r"\1", l).rstrip()
        # skip comments and emty rows
        if recomm.match(l) or reblank.match(l): continue
        
        # skip the comments at the end of the row
        l=comm.sub(r"\1",l).rstrip()
        
        # split in fields
        flds=l.split("\t")
        
        #print "proceeding:"+l;##
        if len(flds) == 1:
            # new section starts here
            sec_name=flds[0]
            subsec_name=""
            if not sec_name in defsec:
                #import pdb; pdb.set_trace()
                raise Exception("FTBL: Illegal section name '%s' (%s: %d)"%(sec_name, ftbl["name"], irow))
            # prepare storage
            ftbl[sec_name]=[]
            # prepare storage for carbon transitions
            if sec_name=="NETWORK":
               #print "TRANS prepared";##
               ftbl["TRANS"]=[]
            try: del stock
            except NameError: pass
            stock=ftbl[sec_name]
            skiptab=1
            data_count=0
            reading="col_names"
            col_names=[]
            subsec_name=""
            continue
        if len(flds) == 2 and len(flds[0]) == 0:
            # read subsection name or what ?
            if len(sec_name) and sec_name in defsec and len(defsec[sec_name]):
                # we are expecting a subsection
                subsec_name=flds[1]
                if subsec_name not in defsec[sec_name]:
                    raise Exception("A subsection '%s' cannot appear in the section '%s' (%s: %d)."%(subsec_name, sec_name, ftbl["name"], irow))
                # prepare sub-storage
                if not ftbl[sec_name]:
                    # replace an empty list by an empty dictionary
                    ftbl[sec_name]=dict()
                #print (irow, reading)##
                ftbl[sec_name][subsec_name]=[]
                try: del stock
                except NameError: pass
                stock=ftbl[sec_name][subsec_name]
                skiptab=2
                data_count=0
                reading="col_names"
                continue
            else:
                # just a very short line
                # it will fall in plain reading data
                pass
        if reading=="col_names" and len(flds) > 2:
            # read column names
            if len(l) < skiptab or l[:skiptab] != "\t"*skiptab:
                raise Exception("Expected at least %d tabulation(s) at the row beginning. Got '%s' (%s: %d)"%(skiptab, l[:min(skiptab, len(l))], ftbl["name"], irow))
            col_names=l[skiptab:].split("\t")
            if len([ item for item in col_names if re.match("^\s*$", item) ]):
                raise Exception("FTBL: row %d has empty column names:\n%s"%(irow,l))
            reading="data"
            #print "col_names=", col_names;##
            continue
        if reading=="data" or reading=="transitions":
            if len(l) < skiptab or l[:skiptab] != "\t"*skiptab:
                raise Exception("Expected at least %d tabulation(s) at the row beginning. Got '%s' (%s: %d)"%(skiptab, l[:min(skiptab, len(l))], ftbl["name"], irow))
            data=[it.strip() for it in l[skiptab:].split("\t")]
            prevdic=dic
            dic={"irow": str(irow)}
            if reading=="data" and sec_name == "NETWORK" and len(data[0]) != 0:
                reading="transitions"
            elif reading=="transitions":
                # here, we are at carbon transition line (e.g. #ABC -> #AB +#C)
                #print "data_count="+str(data_count), \
                #    "\ndata="+str([l for l in enumerate(data)]), \
                #    "\nstock="+str(stock);##
                if sec_name == "NETWORK" and len(data[0]) != 0:
                    raise Exception("Expected label transitions. Got '%s' (%s: %d)"%(l, ftbl["name"], irow))
                reading="data"
                fl_name=str(stock[data_count-1][col_names[0]]) if data_count else ""
                if col_names[0] not in prevdic or len(prevdic[col_names[0]]) == 0:
                    raise Exception("Carbon transition row '%s' is orphan (%s: %d)."%(l[skiptab:], ftbl["name"], irow))
                for i in range(len(col_names)):
                    item=data[i] if i < len(data) else ""
                    dic[col_names[i]]=item
                    metab=stock[data_count-1][col_names[i]]
                    if i > 0 and ((len(metab) and not len(item)) or (not len(metab) and len(item))):
                        #print "i=%d, co='%s', m='%s', tr='%s';"%(i, col_names[i], metab, item)
                        raise Exception("In the reaction '%s', metabolites are misaligned with carbon transitions (%s: %d)."%(fl_name, ftbl["name"], irow))
                ftbl["TRANS"].append(dic)
                continue
            for i in range(len(col_names)):
                # classic data
                if len(data) > len(col_names):
                    raise Exception("FTBL: data have more columns (%d) than column names (%d) (%s: %d)"%(len(data), len(col_names), ftbl["name"], irow))
                try:
                    # decimal point conversion
                    dic[col_names[i]]=data[i].strip() if i < len(data) else ""
                    if col_names[i] in float_conv:
                        oldval=dic[col_names[i]]
                        val=oldval.replace(",", ".")
                        try:
                            fval=float(val)
                            dic[col_names[i]]=val
                        except:
                            dic[col_names[i]]=oldval
                except IndexError:
                    pass
            if sec_name == "NETWORK" and pathway:
                if pathway not in ftbl["pathway"]:
                    ftbl["pathway"][pathway]=[]
                ftbl["pathway"][pathway]+=[dic["FLUX_NAME"]]
            if sec_name == "FLUXES" and subsec_name in ("NET", "XCH") and dic["FCD"] in ("F", "C"):
                try:
                    val=float(eval(dic["VALUE(F/C)"]))
                except:
                    raise Exception("In the field 'VALUE(F/C)', a float value expected (%s: %d)"%(ftbl["name"], irow))
            stock.append(dic)
            data_count+=1
            #print "sec, subsec=", sec_name, subsec_name, data_count, dic;##
        #print "len(flds)=", len(flds), flds, l, data;##
        #print "keys", ftbl.keys(), (ftbl[sec_name].keys() \
        #        if type(ftbl[sec_name])==type({}) else "");##
    if "NETWORK" not in ftbl:
        return ftbl
    # prepare translator reac -> pathway
    ftbl["reac2path"]=dict((reac,path) for path,li in ftbl["pathway"].items() for reac in li)
    # assemble reactions with the same name in one
    # long_reac={reac: {"left": [(minp1, 1), (minp2, 1), ...], "right": [(mout1, 1), (mout2, 1), ...]}
    # and carbon transitions
    # long_trans={reac: {"left": [minp1, minp2, ...], "right": [mout1, mout2, ...]}
    # get unique reaction names
    nw=ftbl["NETWORK"]
    tr=ftbl["TRANS"]
    if (len(nw) != len(tr)):
        raise Exception("Number of reactions (%d) is not equal to label transition number (%d)"%(len(nw), len(tr)) )
    ureac=oset(row["FLUX_NAME"] for row in ftbl["NETWORK"])
    long_reac=dict()
    long_trans=dict()
    for reac in ureac:
        if reac not in long_reac:
            long_reac[reac]={"left": [], "right": []}
            long_trans[reac]={"left": [], "right": []}
        # get rows
        irows=[i for (i, row) in enumerate(nw) if row["FLUX_NAME"] == reac]
        if len(irows) > 1:
            # check that rows are contiguous
            if not all(r == i+irows[0] for i,r in enumerate(irows)):
                raise Exception("Reaction '%s' is split on non contiguous rows: %s"%(reac, ", ".join(nw[i]["irow"] for i in irows)))
        long_reac[reac]["irow"]=", ".join(str(nw[i]["irow"]) for i in irows)
        long_reac[reac]["left"]+=[mecoparse(m) for i in irows for m in (nw[i]["EDUCT_1"], nw[i]["EDUCT_2"]) if m]
        long_reac[reac]["right"]+=[mecoparse(m) for i in irows for m in (nw[i]["PRODUCT_1"], nw[i]["PRODUCT_2"]) if m]
        long_trans[reac]["irow"]=", ".join(str(tr[i]["irow"]) for i in irows)
        long_trans[reac]["left"]+=[m for i in irows for m in (tr[i]["EDUCT_1"], tr[i]["EDUCT_2"]) if m]
        long_trans[reac]["right"]+=[m for i in irows for m in (tr[i]["PRODUCT_1"], tr[i]["PRODUCT_2"]) if m]
        for lr in ("left", "right"):
            nmet=len(long_reac[reac][lr])
            ncarb=len(long_trans[reac][lr])
            if nmet != ncarb:
                raise Exception("Number of metabolites (%d) and carbon strings (%d) on the %s hand side of reaction 'reac' must be equal (%s: %s)"%(nmet, ncarb, lr, reac, ftbl["name"], long_reac[reac]["irow"]))
    ftbl["long_reac"]=long_reac
    ftbl["long_trans"]=long_trans
    return ftbl

def ftbl_netan(ftbl, netan, emu_framework=False, fullsys=False, case_i=False):
    """
    analyse ftbl dictionary to find
     
     - network inputs (input)
     - network outputs (output)
     - substrates (subs)
     - products (prods)
     - metabolites (metabs)
     - reactions (reacs)
     - not reversible reactions (subset of reacs) (notrev)
       all above items are in named sets
     - stocheometric matrix (sto_r_m)
     - stocheometric matrix (sto_m_r)
     - fwd-rev flux matrix (flux_m_r)
     - cumomer balances (cumo_m_r_m)
     - carbon length (Clen)
     - reaction formula (formula)
     - metabolite network (metab_netw)
     - carbon transitions (carbotrans)
     - free fluxes (flux_free)
     - constrained fluxes (flux_constr)
     - measured fluxes (flux_measured)
     - variable growth fluxes (flux_vgrowth)
     - input isotopomers (iso_input)
     - input isotopomers functions (funlab for case_i=True)
     - input cumomers (cumo_input)
     - input reduced cumomers (rcumo_input)
     - flux inequalities (flux_ineqal)
     - flux equalities (flux_eqal)
     - label measurements, H1 (label_meas)
     - peak measurements, C13 (peak_meas)
     - mass measurements (mass_meas)
     - cumomer ordered lists (vcumo)
     - unknown fluxes ordered lists (vflux)
     - linear problem on fluxes (Afl, bfl)
     - free fluxes ordered lists (vflux_free)
     - fw-rv fluxes ordered lists (vflux_fwrv)
     - row names ordered lists for Afl (vrowAfl)
     - in-out fluxes (flux_in, flux_out)
     - measured concentrations (metab_measured)
    """
    # init named sets
    if type(netan)!=type(dict()):
        raise("netan argument must be a dictionary")
    if not netan:
        netan.update({
            "input":oset(),
            "output":oset(),
            "deadend":oset(),
            "subs":oset(),
            "prods":oset(),
            "left":oset(),
            "right":oset(),
            "metabs":oset(),
            "reac":oset(),
            "notrev":oset(),
            "sto_r_m":dict(),
            "sto_m_r":dict(),
            "flux_m_r":dict(),
            "cumo_balance_pattern":dict(),
            "Clen":dict(),
            "formula":dict(),
            "metab_netw":dict(),
            "cumo_sys":dict(),
            "carbotrans":dict(),
            "Cmax":dict(),
            "flux_dep":dict(),
            "flux_free":dict(),
            "flux_constr":dict(),
            "flux_measured":dict(),
            "iso_input":[],
            "funlab":[],
            "cumo_input":[],
            "rcumo_input":[],
            "emu_input":[],
            "flux_inequal":{"net":[], "xch":[]},
            "flux_equal":{"net":[], "xch":[]},
            "label_meas":[],
            "peak_meas":[],
            "mass_meas":[],
            "vcumo":[],
            "Afl":[],
            "bfl":[],
            "vflux":{"net":[], "xch":[], "net2i":dict(), "xch2i":dict()},
            "vflux_free":{"net":[], "xch":[], "net2i":dict(), "xch2i":dict()},
            "vflux_constr":{"net":[], "xch":[], "net2i":dict(), "xch2i":dict()},
            "vflux_meas":{"net":[], "net2i":dict()},
            "vflux_growth":{"net":[], "net2i":dict()},
            "vflux_compl":{"net":[], "xch":[], "net2i":dict(), "xch2i":dict()},
            "vflux_fwrv":{"fw":[], "rv":[], "fw2i":dict(), "rv2i":dict()},
            "vrowAfl":[],
            "flux_in":oset(),
            "flux_out":oset(),
            "flux_inout":oset(),
            "opt":dict(),
            "met_pools":dict(),
            "nx2dfcg":dict(),
            "metab_measured":dict(),
        })
    netan["emu"]=emu_framework
    res="";     # auxiliary short-cut to current result
    netan["exp_names"]=[ftbl["base_name"]]
    netan["fullsys"]=fullsys
    netan["pathway"]=ftbl["pathway"]
    netan["reac2path"]=ftbl["reac2path"]

    # Isotopomer dynamics and other options
    for row in ftbl.get("OPTIONS",[]):
        try:
            netan["opt"][row["OPT_NAME"]]=eval(row["OPT_VALUE"])
        except:
            netan["opt"][row["OPT_NAME"]]=row["OPT_VALUE"]
    for row in ftbl.get("METABOLITE_POOLS",[]):
        metab=row["META_NAME"]
        if metab in netan["met_pools"]:
            raise Exception("Metabolite '%s' is present more than once in the\nftbl secion METABOLITE_POOLS (second appearance on row %s)"%(metab, row["irow"]))
        netan["met_pools"][metab]=eval(row["META_SIZE"])

    # check the presence of fields "NAME", "FCD" and maybe "VALUE(F/C)"
    for suf in ["NET", "XCH"]:
        for row in ftbl.get("FLUXES", dict()).get(suf, dict()):
            if not "NAME" in row:
                raise Exception("No required field NAME in section FLUX/%s (%s: %s)"%(suf, ftbl["name"], row["irow"]))
            if not "FCD" in row:
                raise Exception("No requied field FCD in section FLUX/%s (%s: %s)"%(suf, ftbl["name"], row["irow"]))
            if (row["FCD"]=="F" or row["FCD"]=="C") and not "VALUE(F/C)" in row:
                raise Exception("For flux '%s' in section FLUX/%s, the field 'VALUE(F/C)' is requiered but is absent (%s: %s)"%(row["NAME"], suf, ftbl["name"], row["irow"]))
    # quick consistency check between net and xch fluxes
    fnet=oset(row["NAME"] for row in ftbl.get("FLUXES", dict()).get("NET", dict()))
    fxch=oset(row["NAME"] for row in ftbl.get("FLUXES", dict()).get("XCH", dict()))
    net_wo_xch=fnet-fxch
    if net_wo_xch:
        raise Exception("Fluxe(s) '%s' are defined in the FLUX/NET section but not in the XCH one."%join(", ", net_wo_xch))
    xch_wo_net=fxch-fnet
    if xch_wo_net:
        raise Exception("Fluxe(s) '%s' are defined in the FLUX/XCH section but not in the NET one."%join(", ", xch_wo_net))

    # quick not reversible reactions for complete subs and prods accounting
    revreac=oset(row["NAME"] for row in ftbl.get("FLUXES", dict()).get("XCH", dict()) if row["FCD"]=="F" or (row["FCD"]=="C" and (eval(row["VALUE(F/C)"]) != 0.)))
    # analyse networks
    netw=ftbl.get("long_reac")
    if not netw:
        raise Exception("No long_reac section in the ftbl parameter")
    row_to_del=[]
    for (reac, row) in netw.items():
        #print "reac="+reac;#
        # corresponding carbon transition row
        crow=ftbl["long_trans"][reac]
        # local substrate (es), product (ps) and metabolites (ms) sets
        es=oset(m for m,_ in row["left"])
        ps=oset(m for m,_ in row["right"])
        ms=es|ps
        #if es&ps:
        #   raise Exception("The same metabolite(s) '%s' are present in both sides of a reaction '%s' (%s: %s)."%(join(", ", es&ps), reac, ftbl["name"], row["irow"]))
        
        # all reactions A+B=C or C=A+B or A+B=C+D
        netan["reac"].add(reac)
        netan["left"].update(es)
        netan["right"].update(ps)
        netan["subs"].update(es)
        netan["prods"].update(ps)
        netan["metabs"].update(ms)
        if reac in revreac:
            netan["subs"].update(ps)
            netan["prods"].update(es)
        #aff("ms for "+reac, ms);##

        # create chemical formula
        netan["formula"][reac]={"left": row["left"], "right": row["right"], "all":ms}

        # Carbon length and transitions
        netan["carbotrans"][reac]={"left": [], "right": []}
        for (m, carb, lr) in [(row[lr][i][0], crow[lr][i], lr) for lr in ("left", "right") for i in range(len(row[lr]))]:
                #print "m="+str(m), "; carb="+str(carb);##
            if carb[0] != "#":
                raise Exception("In carbon string '%s' for metabolite '%s' a starting '#' is missing. (%s: %s)"%(carb, m, ftbl["name"], row["irow"]))
            # carbon transitions
            netan["carbotrans"][reac][lr].append((m,carb[1:])); # strip "#" character

            # carbon length
            if netan["Clen"].get(m, 0) and \
                    netan["Clen"][m] != len(carb)-1:
                raise Exception("CarbonLength", "Metabolite "+m+" has length "+
                        str(netan["Clen"][m])+" but in reaction "+reac+
                        " it has length "+str(len(carb)-1)+" (%s: %s)"%(ftbl["name"], row["irow"]))
            netan["Clen"][m]=len(carb)-1; # don't count '#' character
        # check equal carbon length on left and right sides
        lenl=sum(len(l) for m,l in netan["carbotrans"][reac]["left"])
        lenr=sum(len(l) for m,l in netan["carbotrans"][reac]["right"])
        if lenl!=lenr:
            raise Exception("Carbon patterns have different lengths on both sides of reaction '%s' (%d on left, %d on right, %s: %s)"%(reac, lenl, lenr, ftbl["name"], row["irow"]))
        # check unique presence of each letter on each side
        lets=dict()
        uni=dict()
        for lr in ("left", "right"):
            lets[lr]="".join(l for m,l in netan["carbotrans"][reac][lr])
            uni[lr]=oset(lets[lr])
            if len(uni[lr]) != len(lets[lr]):
                # find repeated letters
                di=dict((l,lets[lr].count(l)) for l in uni[lr])
                for (l,c) in di.items():
                    if c > 1:
                        raise Exception("Character '%s' is present %s on the %s side of carbon transition in reaction '%s' (%s: %s)"%(l, ntimes(c), lr, reac, ftbl["name"], row["irow"]))
        # check for perfect mapping
        lmr=uni["left"]-uni["right"]
        if lmr:
            raise Exception("Letter(s) '%s' are present on the left but not on the right hand side in carbon transitions for reaction '%s' (%s: %s)"%(", ".join(sorted(lmr)), reac, ftbl["name"], row["irow"]))
        rml=uni["right"]-uni["left"]
        if rml:
            raise Exception("Letter(s) '%s' are present on the right but not on the left hand side in carbon transitions for reaction '%s' (%s: %s)"%(", ".join(sorted(lmr)), reac, ftbl["name"], row["irow"]))

        # stocheometric matrix in dictionnary form
        # sto_r_m[reac]['left'|'right']=[(metab, coef)] and
        # sto_m_r[metab]['left'|'right']=list((reac, coef))
        # substrates are in 'left' list
        # and products are in the 'right' one.

        netan["sto_r_m"][reac]={"left":row["left"], "right":row["right"]}
        
        for s,c in netan["sto_r_m"][reac]["left"]:
            #print "sto_m_r s="+str(s);##
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]}
            netan["sto_m_r"][s]["left"].append((reac, c))
        for s,c in netan["sto_r_m"][reac]["right"]:
            #print "sto_m_r s="+str(s)
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]}
            netan["sto_m_r"][s]["right"].append((reac, c))
        # end netan["sto_m_r"]
        #aff("sto_m_r"+str(reac), netan["sto_m_r"]);##

    for row in ftbl.get("NOTRACER_NETWORK", []):
        reac=row["FLUX_NAME"]
        if reac in netan["sto_r_m"]:
            raise Exception("Reaction '%s' from NOTRACER_NETWORK section was already introduced in NETWORK section (%s: %s)."%(reac, ftbl["name"], row["irow"]))
        spl=row["EQUATION"].split("=")
        if len(spl) < 2:
            raise Exception("Not found '=' sign in NOTRACER_NETWORK reaction: '%s' (%s: %s)"%(reac, ftbl["name"], row["irow"]))
        elif len(spl) > 2:
            raise Exception("Too many '=' signs (%d) in NOTRACER_NETWORK reaction: '%s' (%s: %s)"%(reac, len(spl)-1, ftbl["name"], row["irow"]))
        lr=("left", "right")
        di=dict()
        for i,su in enumerate(spl):
            terms=su.split("+")
            try:
                di[lr[i]]=mecoparse(terms)
            except:
                werr("Error occured in '%s': %s\n"%(ftbl["name"], row["irow"]))
                raise
        netan["sto_r_m"][reac]=di
        for s,c in netan["sto_r_m"][reac]["left"]:
            #print "sto_m_r s="+str(s);##
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]}
            netan["sto_m_r"][s]["left"].append((reac, c))
        for s,c in netan["sto_r_m"][reac]["right"]:
            #print "sto_m_r s="+str(s)
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]}
            netan["sto_m_r"][s]["right"].append((reac, c))
        netan["reac"].add(reac)
        es=oset(m for m,co in di["left"])
        ps=oset(m for m,co in di["right"])
        ms=es|ps
        for m in ms:
            if m not in netan["Clen"]:
                netan["Clen"][m]=0
        netan["left"].update(es)
        netan["right"].update(ps)
        netan["subs"].update(es)
        netan["prods"].update(ps)
        netan["metabs"].update(ms)
        if reac in revreac:
            netan["subs"].update(ps)
            netan["prods"].update(es)
        #aff("ms for "+reac, ms);##

        # create chemical formula
        netan["formula"][reac]=di
        netan["formula"][reac]["all"]=ms

    # find input and output metabolites
    #import pdb; pdb.set_trace()
    netan["input"].update(netan["subs"]-netan["prods"])
    netan["output"].update(netan["prods"]-netan["subs"])
    netan["deadend"].update(((netan["left"]-netan["right"]) | (netan["right"]-netan["left"])) - netan["input"]-netan["output"])
    # internal metabs
    netan["metabint"]=netan["metabs"].copy()
    netan["metabint"].difference_update(netan["input"] | netan["output"])

    # all met_pools must be in internal metabolites
    mdif=oset(netan["met_pools"]).difference(netan["metabint"])
    if len(mdif) :
        # unknown metabolite
        raise Exception("Unknown metabolite(s). Metabolite(s) '"+", ".join(mdif)+"' defined in section METABOLITE_POOLS are not internal metabolites in NETWORK section.")
    if case_i:
        # check it other way: all metabint must be in metpools
        mdif=oset(netan["metabint"]).difference(netan["met_pools"])
        if len(mdif) :
            # unknown metabolite
            raise Exception("Unknown metabolite concentration. Metabolite(s) '"+", ".join(mdif)+"' defined in section NETWORK are not defined in METABOLITE_POOLS section.")

    # add growth fluxes if requested
    netan["flux_growth"]={"net":dict()}
    netan["flux_vgrowth"]={"net":dict(), "xch":dict()}; # fluxes depending on variable pools
    if netan["opt"].get("include_growth_flux"):
        if not netan["opt"].get("mu"):
            raise Exception("Parameter include_growth_flux is set to True but the growth parameter mu is absent or zero in OPTIONS section")
        for (m,si) in netan["met_pools"].items():
            mgr=m+"_gr"
            reac=mgr
            if reac in netan["reac"]:
                raise Exception("Cannot add growth reaction "+reac+". It is already in the network")
            netan["reac"].add(reac)
            netan["subs"].add(m)
            netan["prods"].add(mgr)
            netan["metabs"].add(m)
            netan["output"].add(mgr)
            if si > 0 :
               netan["flux_growth"]["net"][reac]=None
            else:
               netan["flux_vgrowth"]["net"][reac]=None
               netan["flux_vgrowth"]["xch"][reac]=0.
            #aff("ms for "+reac, ms);##

            # create chemical formula
            netan["formula"][reac]={"left":[(m, 1.)], "right":[(mgr, 1)], "all":[m,mgr]}
            # stoicheometry
            netan["sto_m_r"][m]["left"].append((reac, 1.))
            netan["sto_m_r"][mgr]={"left": [], "right": [(reac, 1.)]}
            netan["sto_r_m"][reac]={"left":[(m, 1.)], "right":[(mgr, 1.)]}

    # check metab names in pools and network
    # all met_pools must be in internal metabolites
    mdif=oset(netan["met_pools"]).difference(netan["metabint"])
    if len(mdif) :
        # unknown metabolite
        raise Exception("Unknown metabolite(s). Metabolite(s) '"+", ".join(mdif)+"' defined in the section METABOLITE_POOLS are not internal metabolites in the NETWORK section.")

    if netan["met_pools"] and me=="ftbl2labprop.py3":
        # all internal metabs must be also in met_pools section
        mdif=netan["metabint"].difference(netan["met_pools"])
        if len(mdif) :
            # unknown metabolite(s)
            raise Exception("Unknown metabolite(s). Metabolite(s) '"+", ".join(mdif)+"' defined as internal in the section NETWORK, are not defined in the METABOLITE_POOLS section.")
    
    # input and output flux names
    netan["flux_in"]=oset(f for m in netan["input"] for f,_ in netan["sto_m_r"][m]["left"])
    netan["flux_out"]=oset(f for m in netan["output"] for f,_ in netan["sto_m_r"][m]["right"])
    netan["flux_inout"]=netan["flux_in"] | netan["flux_out"]
    #print "fl in + out="+str(netan["flux_in"])+"+"+str( netan["flux_out"]);##
    
    # get possible additional fluxes (i.e. from cofactors) from EQAULITY section (eqflux)
    # flux equalities
    # list of tuples (value,dict) where dict is flux:coef
    # net fluxes
    
    # temporary list of constrained net fluxes
    fcnstr=oset(row["NAME"] for row in ftbl.get("FLUXES", dict()).get("NET", []) if row["FCD"] == "C")
    for row in ftbl.get("EQUALITIES", dict()).get("NET",[]):
        dicf=formula2dict(row["FORMULA"])
        # number of non constrained fluxes in a formula
        nb_nonc=sum(fl not in fcnstr for fl in dicf)
        if nb_nonc==0:
            wout("Warning: in EQUALITIES/NET section, the formula '"+
                row["VALUE"]+"="+row["FORMULA"]+"' involves only constrained flux(es)\n.The equality is ignored as meaningless (%s: %s).\n"%(ftbl["name"], row["irow"]))
            continue
        netan["flux_equal"]["net"].append((
                eval(row["VALUE"]),
                dicf,
                row["FORMULA"]+"="+row["VALUE"]+": "+str(row["irow"])))
    # xch fluxes
    # temporary list of constrained xch fluxes
    fcnstr=oset(row["NAME"] for row in ftbl.get("FLUXES", dict()).get("XCH", []) if row["FCD"] == "C")
    for row in ftbl.get("EQUALITIES", dict()).get("XCH",[]):
        dicf=formula2dict(row["FORMULA"])
        # number of non constrained fluxes in a formula
        nb_nonc=sum(fl not in fcnstr for fl in dicf)
        if nb_nonc==0:
            wout("Warning: in EQUALITIES/XCH section, the formula '"+
                row["VALUE"]+"="+row["FORMULA"]+"' involves only constrained flux(es)\n.The equality is ignored as meaningless (%s: %s).\n"%(ftbl["name"], row["irow"]))
            continue
        netan["flux_equal"]["xch"].append((
                eval(row["VALUE"]),
                dicf,
                row["FORMULA"]+"="+row["VALUE"]+": "+str(row["irow"])))
    netan["eqflux"]=oset(f for row in netan["flux_equal"]["net"] for f in [*row[1].keys()]) | oset(f for row in netan["flux_equal"]["xch"] for f in [*row[1].keys()])
    eqflux=netan["eqflux"]
    # metab EQAULITIES
    netan["metab_equal"]=list()
    for row in ftbl.get("EQUALITIES", dict()).get("METAB",[]):
        #print row;##
        dicf=formula2dict(row["FORMULA"])
        nb_neg=0
        for m in dicf:
            if m not in netan["metabint"]:
                raise Exception("Metabolite `%s` is not internal metabolite (%s: %s)."%(m, ftbl["name"], row["irow"]))
            if m not in netan["met_pools"]:
                raise Exception("Metabolite `%s` is not declared in METABOLITE_POOLS section (%s: %s)."%(m, ftbl["name"], row["irow"]))
            nb_neg+=netan["met_pools"][m]<0
        if nb_neg==0:
            raise Exception("At least one of metabolites '%s' must be declared as variable (i.e. having negative value) in the section METABOLITE_POOLS (%s: %s)."%("', '".join(list(dicf.keys())), ftbl["name"], row["irow"]))
        netan["metab_equal"].append((
                eval(row["VALUE"]),
                dicf,
                row["FORMULA"]+"="+row["VALUE"]+": "+str(row["irow"])))

    # analyse reaction reversibility
    # free fluxes(dictionary flux->init value for simulation or minimization)
    # constrained fluxes (dictionary flux->value)
    # variable growth fluxes (dictionary flux->value)
    # The rest are dependent fluxes
    netan["flux_dep"]=dict({"net":dict(), "xch":dict()})
    netan["flux_free"]=dict({"net":dict(), "xch":dict()})
    netan["flux_constr"]=dict({"net":dict(), "xch":dict()})
    flx=ftbl.get("FLUXES")
    fcd=oset(("F", "C", "D"))
    if flx:
        xch=flx.get("XCH")
        net=flx.get("NET")
##        aff("net", net)
        # check that all fluxes are defined in network section
        allreac=netan["reac"] | netan["flux_inout"]
        unk=[ row["NAME"] for row in net if row["FCD"] in fcd and row["NAME"] not in allreac|eqflux ]
        if len(unk):
            raise Exception("The flux name(s) '%s' from the FLUX/NET section is (are) not defined in the NETWORK neither EQUALITY section."%(", ".join(unk)))
        unk=[ row["NAME"] for row in xch if row["FCD"] in fcd and row["NAME"] not in allreac|eqflux ]
        if len(unk):
            raise Exception("The flux name(s) '%s' from the FLUX/XCH section is (are) not defined in the NETWORK neither EQUALITY section."%(", ".join(unk)))

        # check that all reactions from NETWORK are at least in FLUX/NET section
        fnet=oset( row["NAME"] for row in net if row["NAME"] and row["FCD"] in fcd )
        #print("fnet=", fnet)
        #print("allreac=", allreac)
        unk=allreac-fnet-oset(netan["flux_growth"]["net"])-oset(netan["flux_vgrowth"]["net"])
        if len(unk):
            raise Exception("The flux name(s) '%s' from the NETWORK section is (are) not defined in the FLUX/NET section."%(", ".join(unk)))

        # check that all reactions from EQUALITY are at least in FLUX/NET/XCH section
        fnetxch=oset( row["NAME"] for row in net+xch if row["NAME"] and row["FCD"] in fcd )
        unk=eqflux-fnet-oset(netan["flux_growth"]["net"])-oset(netan["flux_vgrowth"]["net"])
        if len(unk):
            raise Exception("The flux name(s) '%s' from the NETWORK section is (are) not defined in the FLUX/NET section."%(", ".join(unk)))
    else:
        werr(os.path.basename(me)+": netan[FLUXES] is not defined\n")
        #return None
    #print "list reac=", netan["reac"];##
    try:
        #print( netan["reac"] | netan["flux_inout"])
        for reac in netan["reac"] | netan["flux_inout"] | eqflux:
            #print("reac=", reac)
            # get xch condition for this reac
            cond=[row for row in xch if row["NAME"]==reac]
            # get net condition for this reac
            ncond=([row for row in net if row["NAME"]==reac]) if net else []
            # no xch dispatch check for input/output fluxes as they are
            # constrained by definition
            #print "r,c,n=", reac, len(cond), len(ncond);##
            if len(cond) > 1:
                raise Exception("FluxDef", "Reaction `%s` is over defined in exchange fluxes."%reac)
            if len(cond) == 0 and not (reac in (netan["flux_in"] | netan["flux_out"])):
                werr("in+out fluxes="+str(netan["flux_inout"])+"\n");##
                raise Exception("FluxDef", "Reaction `%s` is not defined D/F/C in exchange fluxes."%reac)
            if cond:
                cond=cond[0]
            if len(ncond) > 1:
                raise Exception("FluxDef", "Reaction `%s` is over defined in net fluxes."%reac)
            if len(ncond) == 0 and not (reac in (netan["flux_inout"])):
                raise Exception("FluxDef", "Reaction `%s` is not defined D/F/C in net fluxes."%reac)
            if (len(ncond) == 0 and
                (reac in netan["flux_inout"] and
                not reac in netan["flux_growth"]["net"] and
                not reac in netan["flux_vgrowth"]["net"])):
                # input-output fluxes are by definition not reversible
                netan["notrev"].add(reac)
                netan["flux_constr"]["xch"][reac]=0.
                continue
            if netan["opt"].get("include_growth_flux") and reac in netan["flux_growth"]["net"]:
                m=reac[:-3]
                netan["notrev"].add(reac)
                netan["flux_constr"]["xch"][reac]=0.
                netan["flux_constr"]["net"][reac]=netan["opt"]["mu"]*netan["met_pools"][m]
                continue
            if netan["opt"].get("include_growth_flux") and reac in netan["flux_vgrowth"]["net"]:
                m=reac[:-3]
                netan["notrev"].add(reac)
                netan["flux_vgrowth"]["net"][reac]=netan["opt"]["mu"]*abs(netan["met_pools"][m])
                netan["flux_vgrowth"]["xch"][reac]=0.
                continue
            ncond=ncond[0]
            #aff("cond", cond);##
            #aff("ncond", ncond);##
            # not reversibles are those reaction having xch flux=0 or
            # input/output fluxes
            try:
                xval=float(eval(cond["VALUE(F/C)"]))
            except:
                xval=0.
            try:
                nval=float(eval(ncond["VALUE(F/C)"]))
            except:
                nval=0.
            try:
                if (cond and cond["FCD"] == "C") or (reac in netan["flux_inout"]):
                    if (reac in netan["flux_inout"]) or \
                        (xval == 0.):
                        netan["notrev"].add(reac)
                        netan["flux_constr"]["xch"][reac]=0.
                    else:
                        # not null constrained xch
                        if xval < 0. or xval > 1.:
                            raise Exception("FluxVal", "Constrained exchange flux`%s` must have value in [0; 1] interval.\nInstead, a value %f is given"%(reac, xval))
                        netan["flux_constr"]["xch"][reac]=xval
                elif (cond and cond["FCD"] == "F"):
                    # free xch
                    if xval < 0. or xval > 1.:
                        raise Exception("FluxVal", "Free exchange flux '%s' must have starting value in [0; 1] interval.\nInstead, a value %f is given"%(reac, xval))
                    netan["flux_free"]["xch"][reac]=xval
                elif (cond and cond["FCD"] == "D"):
                    # dependent xch
                    netan["flux_dep"]["xch"][reac]=xval
                if (ncond and ncond["FCD"] == "F"):
                    # free net
                    netan["flux_free"]["net"][reac]=nval
                elif (ncond and ncond["FCD"] == "C"):
                    # constr net
                    netan["flux_constr"]["net"][reac]=nval
                elif (ncond and ncond["FCD"] == "D"):
                    # dependent net
                    netan["flux_dep"]["net"][reac]=nval
            except:
                werr("Error occured on the row %d or %d of ftbl file %s:\n"%(ncond["irow"], cond["irow"], ftbl["name"]))
                raise
    except Exception as inst:
        #werr(join(" ", sys.exc_info())+"\n")
        werr(": ".join(inst)+"\n")
    # complete variable growth exchange fluxes
    netan["flux_vgrowth"]["xch"]=dict((k, 0.) for k in netan["flux_vgrowth"]["net"])
    
    # measured fluxes
    for row in ftbl.get("FLUX_MEASUREMENTS",[]):
        if row["FLUX_NAME"] not in netan["reac"]|eqflux:
            raise Exception("Mesured flux `%s` is not defined in NETWORK section neither in EQUALITIES (%s: %s)."%(row["FLUX_NAME"], ftbl["name"], row["irow"]))
        if row["FLUX_NAME"] not in netan["flux_free"]["net"] and \
            row["FLUX_NAME"] not in netan["flux_dep"]["net"]:
            raise Exception("Mesured flux `%s` must be defined as either free or dependent (%s: %s)."%(row["FLUX_NAME"], ftbl["name"], row["irow"]))
        try:
            val=eval(row["VALUE"])
        except:
            val=NaN
        try:
            sdev=float(eval(row["DEVIATION"]))
        except:
            raise Exception("DEVIATION must evaluate to a real number (%s: %s)."%(ftbl["name"], row["irow"]))
        if sdev <= 0.:
            raise Exception("DEVIATION must be positive (%s: %s)."%(ftbl["name"], row["irow"]))
        
        netan["flux_measured"][row["FLUX_NAME"]]={\
                "val": val, \
                "dev": sdev}
    
    # measured concentartions
    for row in ftbl.get("METAB_MEASUREMENTS",[]):
        metabl=row["META_NAME"].split("+")
        found_neg=False
        for m in metabl:
            if m not in netan["metabint"]:
                raise Exception("Mesured metabolite `%s` is not internal metabolite (%s: %s)."%(m, ftbl["name"], row["irow"]))
            if m not in netan["met_pools"]:
                raise Exception("Mesured metabolite `%s` is not declared in METABOLITE_POOLS section (%s: %s)."%(m, ftbl["name"], row["irow"]))
            found_neg=found_neg or netan["met_pools"][m] < 0.
        if not found_neg:
            werr("Warning: metabolite measurements on `%s` does not have a free metabolite (i.e. being negative in the METABOLITE_POOLS ((%s: %s)).\n"%(row["META_NAME"], ftbl["name"], row["irow"]))
            werr("This measurement is ignored\n")
            continue
        try:
            val=float(eval(row["VALUE"]))
        except:
            val=NaN
        try:
            sdev=float(eval(row["DEVIATION"]))
        except:
            raise Exception("DEVIATION must evaluate to a real positive number (%s: %s)."%(ftbl["name"], row["irow"]))
        if sdev <= 0.:
            raise Exception("DEVIATION must be positive (%s: %s)."%(ftbl["name"], row["irow"]))
        netan["metab_measured"][row["META_NAME"]]={\
                "val": val, \
                "dev": sdev}
    
    # proceed LABEL_INPUT
    proc_label_input(ftbl, netan, case_i)
    # flux inequalities
    # list of tuples (value,comp,dict) where dict is flux:coef
    # and comp is of "<", "<=", ...
    # net fluxes
    for row in ftbl.get("INEQUALITIES", dict()).get("NET",[]):
        #print row;##
        if row["COMP"] not in (">=", "=>", "<=", "=<"):
            raise Exception("COMP field in INEQUALITIES section must be one of '>=', '=>', '<=', '=<' and not '%s' (%s: %s)."%(row["COMP"], ftbl["name"], row["irow"]))
        dicf=formula2dict(row["FORMULA"])
        fl=list(dicf.keys())[0]
        if len(dicf)==1 and fl in netan["flux_constr"]["net"]:
            wout("Warning: Inequalities: in NET section, the formula '"+
                row["VALUE"]+row["COMP"]+row["FORMULA"]+"' involves a constrained flux\n"+
                " having a value "+str(netan["flux_constr"]["net"][fl])+". The inequality is ignored as meaningless ((%s: %s)).\n"%(ftbl["name"], row["irow"]))
            continue
        netan["flux_inequal"]["net"].append((
                eval(row["VALUE"]),
                row["COMP"],
                dicf))
    # xch fluxes
    for row in ftbl.get("INEQUALITIES", dict()).get("XCH",[]):
        #print row;##
        if row["COMP"] not in (">=", "=>", "<=", "=<"):
            raise Exception("COMP field in INEQUALITIES section must be one of '>=', '=>', '<=', '=<' and not '%s' (%s: %s)."%(row["COMP"], ftbl["name"], row["irow"]))
        netan["flux_inequal"]["xch"].append((
                eval(row["VALUE"]),
                row["COMP"],
                formula2dict(row["FORMULA"])))
    for (afftype, ftype) in (("Net", "net"), ("Exchange", "xch")):
        for row in netan["flux_inequal"][ftype]:
            for fl in row[2]:
                if fl not in netan["reac"]|eqflux:
                    raise Exception("%s flux `%s` in the inequality\n%s\nis not defined in NETWORK neither EQUALITY sections."%
                        (afftype, fl, join("", row)))
    # metabolite inequalities (like the flux ones)
    netan["metab_inequal"]=list()
    for row in ftbl.get("INEQUALITIES", dict()).get("METAB",[]):
        dicf=formula2dict(row["FORMULA"])
        if row["COMP"] not in (">=", "=>", "<=", "=<"):
            raise Exception("COMP field in INEQUALITIES section must be one of '>=', '=>', '<=', '=<' and not '%s' (%s: %s)."%(row["COMP"], ftbl["name"], row["irow"]))
        nb_neg=0
        for m in dicf:
            if m not in netan["metabint"]:
                raise Exception("Metabolite `%s` is not internal metabolite (%s: %s)."%(m, ftbl["name"], row["irow"]))
            if m not in netan["met_pools"]:
                raise Exception("Metabolite `%s` is not declared in METABOLITE_POOLS section (%s: %s)."%(m, ftbl["name"], row["irow"]))
            nb_neg+=netan["met_pools"][m]<0
        if nb_neg==0:
            raise Exception("At least one of metabolites '%s' must be declared as variable (i.e. having negative value) in the section METABOLITE_POOLS (%s: %s)."%("', '".join(list(dicf.keys())), ftbl["name"], row["irow"]))
        netan["metab_inequal"].append((
                eval(row["VALUE"]),
                row["COMP"],
                dicf,
                row["VALUE"]+row["COMP"]+row["FORMULA"]+": "+str(row["irow"])))

    # Check that fluxes are all in reactions and eqflux
    # then form nx2dfcg dictionary
    for (affnx, nx, nxsh) in (("net", "net", "n."), ("exchange", "xch", "x.")):
        for (affdfcg, dfcg, dfcgsh) in (("Dependent", "flux_dep", "d."), ("Free", "flux_free", "f."), ("Constrained", "flux_constr", "c."), ("Variable growth", "flux_vgrowth", "g.")):
            #print netan[dfcg][nx];##
            for fl in netan[dfcg][nx]:
                if fl not in netan["reac"]|eqflux:
                    raise Exception("%s %s flux `%s` is not defined in NETWORK neither EQUALITY sections."%
                       (affdfcg, affnx, fl))
                netan["nx2dfcg"][nxsh+fl]=dfcgsh+nxsh+fl
    # labeled measurements
    proc_label_meas(ftbl, netan)
    proc_peak_meas(ftbl, netan)
    proc_mass_meas(ftbl, netan)
    
    # discard empty entries
    for e in netan:
        try:
            netan[e].discard("")
            netan[e].discard(None)
        except AttributeError: pass

    # fwd-rev flux balance
    del(res)
    res=netan["flux_m_r"]
    # res[metab]["in"]=list(fwd.flux|rev.flux)
    # res[metab]["out"]=list(fwd.flux|rev.flux)
    for metab,lr in netan["sto_m_r"].items():
        # lr is dico with 'left' and 'right' entries
        #print "fwd metab="+str(metab);##
        if not metab in res:
            res[metab]={"in":[], "out":[]}
        for reac,_ in lr["left"]:
            # here metabolite is consumed in fwd reaction
            # and produced in the reverse one.
            res[metab]["out"].append("fwd."+reac)
            #if reac not in netan["notrev"]:
            if reac not in netan["flux_inout"]:
                res[metab]["in"].append("rev."+reac)
        for reac,_ in lr["right"]:
            # here metabolite is consumed in rev reaction
            # and produced in the forward one.
            res[metab]["in"].append("fwd."+reac)
            #if reac not in netan["notrev"]:
            if reac not in netan["flux_inout"]:
                res[metab]["out"].append("rev."+reac)

    # metabolite network
    res=netan["metab_netw"]
    # left part or reaction points to right part
    for reac,parts in netan["formula"].items():
        for metab,_ in parts["left"]:
            if not metab in res:
                res[metab]=oset()
            res[metab].update(m for m,_ in parts["right"])
    # order cumomers by weight. For a given weight, cumomers are sorted by
    # metabolites order.
    netan["Cmax"]=max(netan["Clen"].values())
    Cmax=netan["Cmax"]
    #return netan;##
    # cumomers systems A*x=b, one by weight
    # weights are going from 1 to Cmax where Cmax is the maximal
    # carbon string length in all metabolites
    # A is a list of matrices (in weight order)
    # b is a list of right hand parts (still in weight order)
    # the dimensions of various weights in b are not the same
    # too short metabolites are dropped when going to higher weights.
    res=netan["cumo_sys"]
    res["A"]=[{} for i in range(Cmax)]
    res["b"]=[{} for i in range(Cmax)]
    try:
        # run through all reactions and update bilan of involved cumomers
        for (reac,lrdict) in iter(netan["carbotrans"].items()) if fullsys else []:
            # run through metabs
            ## aff("lrdict", lrdict);#
            for (imetab,lr,metab,cstr) in ((imetab,lr,metab,cstr)
                    for (lr,lst) in lrdict.items()
                    for (imetab,(metab,cstr)) in enumerate(lst)):
                # if output metab then influx is set to 1
                # so its cumomer distribution is directly
                # defined by cumodistr of inputs
                Clen=netan["Clen"][metab]
                # input metabolite has fixed value so put it in rhs
                # when it is an influx for some internal cumomer
                if metab in netan["input"] or metab in netan["output"]:
                    continue
                # 'out' part of this metab
                fwd_rev=("fwd." if lr=="left" else "rev.")
                flux=fwd_rev+reac
                #if (fwd_rev=="fwd." or reac not in netan["notrev"]):
                if (fwd_rev=="fwd." or reac not in netan["flux_inout"]):
                    # add this out-flux
                    # run through all cumomers of metab
                    for icumo in range(1,1<<Clen):
                        cumo=metab+":"+str(icumo)
                        w=sumbit(icumo)
                        #print "w,i,clen,metab=", w, icumo,Clen,metab;##
                        if cumo not in res["A"][w-1]:
                            res["A"][w-1][cumo]={cumo:[]}
                        # main diagonal term ('out' part)
                        res["A"][w-1][cumo][cumo].append(flux)
                        ##print 'm,ic,w='+metab, icumo, w;#
                        ##aff("res["A"][w-1][cumo][cumo]", res["A"][w-1][cumo][cumo]);#
                # 'in' part
                fwd_rev=("rev." if lr=="left" else "fwd.")
                flux=fwd_rev+reac
                in_lr=("left" if lr=="right" else "right")
                #if (fwd_rev=="rev." and reac in netan["notrev"]):
                if (fwd_rev=="rev." and reac in netan["flux_inout"]):
                    # this cannot be by definition
                    continue
                # add this in-flux
                for (in_i,(in_metab, in_cstr)) in enumerate(lrdict[in_lr]):
                    # run through all cumomers of metab
                    for icumo in range(1,1<<Clen):
                        cumo=metab+":"+str(icumo)
                        w=sumbit(icumo)
                        # get in_cumo
                        in_icumo=src_ind(in_cstr, cstr, icumo)
                        if in_icumo==None:
                            continue
                        in_cumo=in_metab+":"+str(in_icumo)
                        in_w=sumbit(in_icumo)
                        if cumo not in res["A"][w-1]:
                            res["A"][w-1][cumo]={cumo:[]}
                        if in_w==w:
                            if in_metab in netan["input"]:
                                # put it in rhs
                                if cumo not in res["b"][w-1]:
                                    res["b"][w-1][cumo]=dict()
                                if flux not in res["b"][w-1][cumo]:
                                    res["b"][w-1][cumo][flux]=dict()
                                if imetab not in res["b"][w-1][cumo][flux]:
                                    res["b"][w-1][cumo][flux][imetab]=[]
                                if not netan["cumo_input"] or in_cumo not in netan["cumo_input"][0]:
                                    # put this in_cumo ih the dict
                                    iso2cumo(netan, "cumo_input", in_cumo, in_icumo, in_metab)
                                res["b"][w-1][cumo][flux][imetab].append(in_cumo)
                            else:
                                if in_cumo not in res["A"][w-1][cumo]:
                                    res["A"][w-1][cumo][in_cumo]=[]
                                # matrix: linearized off-diagonal term
                                res["A"][w-1][cumo][in_cumo].append(flux)
                        elif in_w < w:
                            # put lighter cumomer product in rhs list[iterm]
                            if cumo not in res["b"][w-1]:
                                res["b"][w-1][cumo]=dict()
                            if flux not in res["b"][w-1][cumo]:
                                res["b"][w-1][cumo][flux]=dict()
                            if imetab not in res["b"][w-1][cumo][flux]:
                                res["b"][w-1][cumo][flux][imetab]=[]
                            #res["b"][w-1][cumo][flux][imetab].append(
                            #    in_cumo if in_metab not in netan["input"] else
                            #    netan["cumo_input"][in_cumo])
                            res["b"][w-1][cumo][flux][imetab].append(in_cumo)
                            #print "b="+str(res["b"][w-1][cumo][flux]);##
                        # if in_w==0 then in_cumo=1 by definition => ignore here
                        # in_w cannot be > w because of src_ind()
    except Exception as inst:
        werr(": ".join(inst)+"\n")
    # ordered cumomer lists
    for w in range(1,netan["Cmax"]+1):
        # weight 1 equations have all metabolites
        ##aff("A "+str(w), netan["cumo_sys"]["A"][w-1]);#
        ##aff("b "+str(w), netan["cumo_sys"]["b"][w-1]);#
        # order cumos along pathways
        # starts are input cumomers
        starts=[cumo for cumo in netan["cumo_sys"]["A"][w-1] \
            if cumo.split(":")[0] in netan["input"]]
        ##aff("st "+str(w), starts)
        # complete starts by all others cumomers
        starts+=[c for c in netan["cumo_sys"]["A"][w-1] if not c in starts]
        cumo_paths=cumo_path(starts, netan["cumo_sys"]["A"][w-1], oset())
        # order
        netan["vcumo"].append([cumo for cumo in valval(cumo_paths)])

    # ordered unknown flux lists
    # get all reactions which are not constrained, not free and not growth
    netan["vflux"]["net"].extend(reac for reac in netan["reac"]|eqflux
        if reac not in netan["flux_constr"]["net"] and
        reac not in netan["flux_free"]["net"] and
        reac not in netan["flux_vgrowth"]["net"])
    netan["vflux"]["xch"].extend(reac for reac in netan["reac"]|eqflux
        if reac not in netan["flux_constr"]["xch"] and
        reac not in netan["flux_free"]["xch"] and
        reac not in netan["flux_vgrowth"]["xch"])

    # order
    netan["vflux"]["net"].sort()
    netan["vflux"]["xch"].sort()
    
    # easy index finder
    netan["vflux"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux"]["net"]))
    netan["vflux"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux"]["xch"]))


    # ordered free flux lists
    netan["vflux_free"]["net"]=list(netan["flux_free"]["net"].keys())
    netan["vflux_free"]["xch"]=list(netan["flux_free"]["xch"].keys())

    # order
    netan["vflux_free"]["net"].sort()
    netan["vflux_free"]["xch"].sort()
    
    # easy index finder
    netan["vflux_free"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_free"]["net"]))
    netan["vflux_free"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_free"]["xch"]))


    # ordered constrained flux lists
    netan["vflux_constr"]["net"]=list(netan["flux_constr"]["net"].keys())
    netan["vflux_constr"]["xch"]=list(netan["flux_constr"]["xch"].keys())

    # order
    netan["vflux_constr"]["net"].sort()
    netan["vflux_constr"]["xch"].sort()
    
    # easy index finder
    netan["vflux_constr"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_constr"]["net"]))
    netan["vflux_constr"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_constr"]["xch"]))


    # ordered measured flux lists
    netan["vflux_meas"]["net"]=list(netan["flux_measured"].keys())

    # order
    netan["vflux_meas"]["net"].sort()
    
    # easy index finder
    netan["vflux_meas"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_meas"]["net"]))

    # ordered variable growth fluxes
    netan["vflux_growth"]["net"]=list(netan["flux_vgrowth"]["net"].keys())
    netan["vflux_growth"]["xch"]=list(netan["flux_vgrowth"]["xch"].keys())

    # order
    netan["vflux_growth"]["net"].sort()
    netan["vflux_growth"]["xch"].sort()
    
    # easy index finder
    netan["vflux_growth"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_growth"]["net"]))
    netan["vflux_growth"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_growth"]["xch"]))

    # ordered complete flux lists
    netan["vflux_compl"]={
        "net": 
        netan["vflux"]["net"]+
        netan["vflux_free"]["net"]+
        netan["vflux_constr"]["net"]+
        netan["vflux_growth"]["net"],
        "xch":
        netan["vflux"]["xch"]+
        netan["vflux_free"]["xch"]+
        netan["vflux_constr"]["xch"]+
        netan["vflux_growth"]["xch"],
    }
    # easy index finding
    # net
    netan["vflux_compl"]["net2i"]=dict(
        (fl,i) for (i,fl) in enumerate(netan["vflux"]["net"])
    )
    netan["vflux_compl"]["net2i"].update(dict(
        (fl,i+len(netan["vflux"]["net"]))
        for (i,fl) in enumerate(netan["vflux_free"]["net"])
    ))
    netan["vflux_compl"]["net2i"].update(dict(
        (fl,i+len(netan["vflux"]["net"])+len(netan["vflux_free"]["net"]))
        for (i,fl) in enumerate(netan["vflux_constr"]["net"])
    ))
    netan["vflux_compl"]["net2i"].update(dict(
        (fl,i+len(netan["vflux"]["net"])+len(netan["vflux_free"]["net"])+len(netan["vflux_constr"]["net"]))
        for (i,fl) in enumerate(netan["vflux_growth"]["net"])
    ))
    # xch
    netan["vflux_compl"]["xch2i"]=dict(
        (fl,i) for (i,fl) in enumerate(netan["vflux"]["xch"])
    )
    netan["vflux_compl"]["xch2i"].update(dict(
        (fl,i+len(netan["vflux"]["xch"]))
        for (i,fl) in enumerate(netan["vflux_free"]["xch"])
    ))
    netan["vflux_compl"]["xch2i"].update(dict(
        (fl,i+len(netan["vflux"]["xch"])+len(netan["vflux_free"]["xch"]))
        for (i,fl) in enumerate(netan["vflux_constr"]["xch"])
    ))
    netan["vflux_compl"]["xch2i"].update(dict(
        (fl,i+len(netan["vflux"]["xch"])+len(netan["vflux_free"]["xch"])+len(netan["vflux_constr"]["xch"]))
        for (i,fl) in enumerate(netan["vflux_growth"]["xch"])
    ))

    # ordered fwd-rev flux lists
    # fw and rv parts are identical so they match each other
    netan["vflux_fwrv"]["fwrv"]=\
        ["fwd."+fl for fl in netan["vflux"]["net"]]+\
        ["fwd."+fl for fl in netan["vflux_free"]["net"]]+\
        ["fwd."+fl for fl in netan["vflux_constr"]["net"]]+\
        ["fwd."+fl for fl in netan["vflux_growth"]["net"]]+\
        ["rev."+fl for fl in netan["vflux"]["net"]]+\
        ["rev."+fl for fl in netan["vflux_free"]["net"]]+\
        ["rev."+fl for fl in netan["vflux_constr"]["net"]]+\
        ["rev."+fl for fl in netan["vflux_growth"]["xch"]]
    
    # order
    netan["vflux_fwrv"]["fwrv"].sort()
    
    # easy index finder
    netan["vflux_fwrv"]["fwrv2i"]=dict((fl,i) for (i,fl) in
        enumerate(netan["vflux_fwrv"]["fwrv"]))
        
    # ordered metabolite pools
    netan["vpool"]={
        "free": [m for (m,v) in netan["met_pools"].items() if v < 0.],
        "constrained": [m for (m,v) in netan["met_pools"].items() if v >= 0.],
    }
    netan["vpool"]["free"].sort()
    netan["vpool"]["constrained"].sort()
    netan["vpool"]["all"]=netan["vpool"]["free"]+netan["vpool"]["constrained"]
    # easy index finding
    netan["vpool"]["all2i"]=dict((m,i) for (i,m) in
        enumerate(netan["vpool"]["all"]))
    
    
    # linear problem on fluxes Afl*(fl_net;fl_xch)=bfl
    # matrix Afl is composed of the following parts :
    # - stocheometic equations (only .net fluxes are involved)
    # - flux equalities
    # constrained to non zero value fluxes are replaced by their values in rhs
    # Afl is a list of lists (two dim array). Primary list is a matrix row,
    # columns are values in the secondary list and contains matrix coefficients
    # bfl is a list of linear expressions. Each expression is a dict
    # where keys are variable names like "f.n.flx" and values are
    # numeric coefficients
    # If the key is empty "", it means that the value is just a scalar to add
    
    # Full flux names are of the format [dfcg].[nx].<reac>
    # where "d", "f", "c" or "g" correspond to
    # - dependent
    # - free
    # - constrained
    # - growth
    # and <reac> correspond to the reaction name
    
    # stocheometric part
    res=netan["Afl"]
    for (metab,lr) in netan["sto_m_r"].items():
        if metab in netan["input"] or metab in netan["output"]:
            continue
        # calculate coefs (repeated fluxes)
        coefs=dict((rea, []) for rea,co in lr["left"]+lr["right"])
        # 'left' part consumes metab, hence the sign '-' in -co
        _=[coefs[rea].append(-co) for rea,co in lr["left"]]
        # 'right' part produces metab
        _=[coefs[rea].append(co) for rea,co in lr["right"]]
        coefs=dict((rea, sum(li)) for rea,li in coefs.items())
        deps=oset(coefs.keys()).intersection(netan["vflux"]["net"])
        if not deps:
            raise Exception("A balance on metabolite '%s' does not contain any dependent flux.\nAt least one of the following net fluxes %s\nmust be declared dependent in the FLUX/NET section (put letter 'D' in the column 'FCD' for some flux)."%(metab, list(coefs.keys())))
        qry=[coefs.get(fl,0) for fl in netan["vflux"]["net"]]
        qry.extend([0]*len(netan["vflux"]["xch"]))
        mqry=-np.array(qry)
        #if qry==[0]*len(qry): must be included even if all zeros, so an R warning will work
        #    # degenerated equation, skip it
        #    #netan["flux_equal"]["net"].append((0., coefs))
        #    raise Exception("Stocheometric equation is zero for metab "+metab+"\n"+str(lr)+"\n"+str(coefs))
        #    continue
        # check if this line was already entered before
        for (i,row) in enumerate(res):
            if not ffguess and (row == qry or (np.array(row) == mqry).all()):
                wout("Warning: when trying to add a balance equation for metabolite '"+metab+
                    "', got equation redundant with those for '"+netan["vrowAfl"][i]+"'\n")
                wout("metab:\t"+join("\t", netan["vflux"]["net"]+netan["vflux"]["xch"])+"\n")
                wout(netan["vrowAfl"][i]+":\t"+join("\t", row)+"\n")
                wout(metab+":\t"+join("\t", qry)+"\n")
                break
        else:
            # identique row is not found, add it
            res.append(qry)
            netan["vrowAfl"].append(metab)
            # prepare right hand side
            netan["bfl"].append({})
            dtmp=netan["bfl"][-1]
            for fl in coefs:
                #print "bfl: fl="+fl;##
                if fl in netan["flux_free"]["net"]:
                    dtmp["f.n."+fl]=dtmp.get("f.n."+fl,0)-coefs[fl]
                elif fl in netan["flux_constr"]["net"]:
                    dtmp["c.n."+fl]=dtmp.get("c.n."+fl,0)-coefs[fl]
                elif fl in netan["flux_vgrowth"]["net"]:
                    dtmp["g.n."+fl]=dtmp.get("g.n."+fl,0)-coefs[fl]
            #print "dtmp=", dtmp

    # flux equality part
    res=netan["Afl"]
    for (nx, nxl) in (("net", "n"), ("xch", "x")):
        for eq in netan["flux_equal"][nx]:
            qry=[]
            qry.extend(eq[1].get(fl,0) for fl in netan["vflux"][nx])
            if nx == "net":
                # add xch zeros
                qry.extend([0]*len(netan["vflux"]["xch"]))
            else:
                # prepend zeros for net fluxes
                qry[0:0]=[0]*len(netan["vflux"]["net"])
            # check qry
            if qry == [0]*len(qry):
                # degenerated equality
                raise Exception("Equality in "+nx.upper()+" section: "+str(eq)+" must have at least one dependent flux\n")
            # check if this line was already entered before
            if not ffguess:
                mqry=-np.array(qry)
                for row in res:
                    if row == qry or (np.array(row) == mqry).all():
                        raise Exception("An equality in "+nx.upper()+" section is redundant. eq:"+str(eq)+
                            "\nqry="+str(qry)+"\nrow="+str(row))
            res.append(qry)
            netan["vrowAfl"].append("eq "+nx+": "+eq[2])
            netan["bfl"].append({"":eq[0]})
            dtmp=netan["bfl"][-1]
            # pass free fluxes to rhs
            for fl in netan["flux_free"][nx]:
                if fl in eq[1]:
                    dtmp["f."+nxl+"."+fl]=dtmp.get("f."+nxl+"."+fl,0)-float(eq[1][fl])
            # pass constrained fluxes to rhs
            for fl in netan["flux_constr"][nx]:
                if fl in eq[1]:
                    dtmp["c."+nxl+"."+fl]=dtmp.get("c."+nxl+"."+fl,0)-float(eq[1][fl])
            # pass growth fluxes to rhs
            for fl in netan["flux_vgrowth"][nx]:
                if fl in eq[1]:
                    dtmp["g."+nxl+"."+fl]=dtmp.get("c."+nxl+"."+fl,0)-float(eq[1][fl])

    # read parallel experiments if any
    proc_kinopt(ftbl, netan)
    if "opt" in netan and "prl_exp" in netan["opt"] and netan["opt"]["prl_exp"]:
        # parse ftbl files
        fli=re.split("\s*;\s*", netan["opt"]["prl_exp"])
        dirw=os.path.dirname(ftbl["abs_path"])
        for fn in fli:
            fn=fn.strip()
            if not fn:
                continue
            fp=ftbl_parse(os.path.join(dirw, fn))
            # test the presence of at least one required field
            if not any(f in fp for f in req_prl):
                raise Exception("Not found label input and/or measurements in '%s'"%fn)
            netan["exp_names"].append(fp["base_name"])
            proc_label_input(fp, netan, case_i)
            proc_label_meas(fp, netan)
            proc_peak_meas(fp, netan)
            proc_mass_meas(fp, netan)
            proc_kinopt(fp, netan)

def enum_path(starts, netw, outs, visited=oset()):
    """Enumerate metabilites along to reaction pathways.
    Algo: start from an input, follow chemical pathways till an output or
    already visited metabolite. Returns a list of metabolite pathways.
    Each pathways is an ordered list."""
    res=[]
    if (len(visited) == len(list(netw.keys()))):
        # no more metabolites 
        return res
    for start in starts:
        if (start in visited):
           continue
        visited.add(start)
        if start in outs:
            return [[start]]
        paths=enum_path(netw[start]-visited, netw, outs, visited)
        if (paths):
            res.append([start]+paths[0])
            if len(paths) > 1:
                res.extend(paths[1:])
        else:
            res.append([start])
    return res
def cumo_path(starts, A, visited=oset()):
    """Enumerate cumomers along reaction pathways.
    Algo: start from an input, follow chemical pathways till no more
    neighbours or till only visited metabolite rest in network.
    Return a list of cumomer pathways.
    Each pathways is an ordered list."""
    #if [s for s in starts if s.find('.') >= 0]:
    #    aff('s', starts)
    res=[]
    if (len(visited) == len(A)):
        # no more cumomers
        return res
    for start in starts:
        if (start in visited or not start in A):
           continue
        visited.add(start)
        # next_starts are cumos influenced by start
        next_starts=oset(cumo for cumo in A if start in A[cumo])
##        aff('next to '+start, next_starts)
        next_starts-=visited
        if next_starts:
            paths=cumo_path(next_starts, A, visited)
##            aff('paths for '+start, paths)
##            aff('visited', visited)
        else:
            path=[]
            res.append([start])
            continue
        if (paths):
            res.append([start]+paths[0])
            if len(paths) > 1:
                res.extend(paths[1:])
        else:
            res.append([start])
##    aff('r', res)
    return res
    
def src_ind(substrate, product, iprod):
    """
    For a given substrate and product carbon strings (e.g. "abc", "ab")
    calculate substrate index corresponding to product index.
    Return None if no source found.
    Return 0 if iprod==0 and intersection of product and substrate strings
    is not empty"""
    movbit=1
    isubstr=0
    substrate=substrate[::-1]
    product=product[::-1]
    for nb in range(len(product)):
        if (movbit & iprod):
            # algo: if the current bit is set in product search for its origin
            # and set the corresponding bit in substrate
            try:
                isubstr|=1<<substrate.find(product[nb])
            except ValueError:
                pass
        movbit<<=1
    return isubstr if (isubstr or
        (iprod==0 and (oset(product) & oset(substrate)))) else None
def labprods(prods, metab, isostr, strs):
    """labprods(prods, metab, isostr, strs)
    Return a set of tuples (vmetab,visostr) which receive at least
    one labeled carbon from (metab, isostr)
    """
    # run through product part of reaction to get
    # metabs containing labeled letter from isostr
    #print "p=", prods, "m=", metab, "i=", isostr, "s=", strs;##
    res=oset()
    # get labeled letters
    lets="".join(s[i] for (i,carb) in enumerate(isostr) if carb=="1" for s in strs)
    #print "ls=", lets;##
    # find labeled letters in reaction products and set
    # them to "1", while others to "0"
    for (vm,vs) in prods:
        rs=""
        for r in vs:
            rs=rs+("1" if r in lets else "0")
            #print "r=",r, "rs=",rs
        if rs.find("1") > -1:
            res.add((vm,rs))
            #print "append=", vm, rs;##
    #print "res=", res;##
    return res
def allprods(srcs, prods, isos, metab, isostr):
    """allprods(srcs, prods, isos, metab, isostr)
    Return a set of tuples (cmetab, cisostr, vmetab, visostr)
    where cmetab and cisostr describe a contex metabolite
    which combined with metab+isostr produced vmetab+visostr.
    if metab is alone on its reaction part cmetab and cisostr are set to
    an empty string ("").
    The set covers all combination of
    metab+isostr and its co-substrates which
    produce isotopes having at least one labeled carbon from
    metab+isostr.
    Co-substrate isotops are
    in a dictionary isos[cmetab]=list(cisotopes)."""
    # run through all combinations of srcs isotopes
    # to get all products
    # |isostr|*|iso2|[+|iso1|*|isostr| if metab is in both position]
    # If there is only one src than m2=""
    #print "allprods: metab=", metab, "isostr=", isostr
    #print "s=", srcs, "p=", prods, "i=", isos;##
    # find metabolite position(s) in the reaction
    mpos=[i for (i,(m,s)) in enumerate(srcs) if m==metab]
    
    # all source isotop couples (im,ic)
    icouples=[]
    if metab==m1:
        icouples=[(isostr,is2) for is2 in isos.get(m2,oset(("",)))]
    if metab==m2:
        icouples.extend((is1,isostr) for is1 in isos.get(m1,oset(("",))))
    #print "icpls=", list(icouples);##
    #sys.exit(1);##
    # labeled source letters
    #print "s1=", s1, "s2=", s2
    lets=["".join(l for (i,l) in enumerate(s1) if is1[i]=="1")+
        "".join(l for (i,l) in enumerate(s2) if is2[i]=="1")
        for (is1,is2) in icouples]
    #print "lets=", lets;##
    # form all products
    res=oset()
    for inm in mpos:
        (m,s)=srcs[inm]
        mlab="".join(l for (i,l) in enumerate(s) if isostr[i]=="1")
        # co-substrate (if any)
        (cm,cs)=("","") if len(srcs)==1 else srcs[(inm+1)%2]
        # run through isotops of co-substrate
        for coiso in isos.get(cm,oset(("",))):
            lab=mlab+"".join(l for (i,l) in enumerate(cs) if coiso[i]=="1")
            # form product isotops
            for (pm,ps) in prods:
                ip="".join(("1" if l in lab else "0") for l in ps)
                if oset(ps)&oset(mlab):
                    # get only products labeled by metab+isostr
                    res.add((cm,coiso, pm, ip))
    #print "res=", res
    return res
def prod(metab, iso, s, cmetab, ciso, cs, prods):
    """prod(metab, iso, s, cmetab, ciso, cs, prods)->oset()
    get isotops from labeled substrates"""
    #print "prod: m=", metab, "s=", s, "i=", iso, "cm=", cmetab, "cs=", cs, "ci=", ciso
    res=oset()
    if cs and ciso == -1:
        return res
    lab=""
    n=len(s)-1
    for (i,b) in iternumbit(iso):
        lab=lab+(s[n-i] if b else "")
    n=len(cs)-1
    for (i,b) in iternumbit(ciso):
        lab=lab+(cs[n-i] if b else "")
    for (pm,ps) in prods:
        ip=sum((1<<i if l in lab else 0) for (i,l) in enumerate(ps[::-1]))
        res.add((pm, ip))
    #print "res=", res
    return res
def frag_prod(metab, frag, s, cmetab, cfrag, cs, prods):
    """frag_prod(metab, frag, s, cmetab, cfrag, cs, prods)
    Get fragments from labeled substrates"""
    #print "frag_prod: m=", metab, "s=", s, "f=", frag, "cm=", cmetab, "cs=", cs, "cf=", cfrag
    res=oset()
    if cs and not cfrag:
        # no fragment for this metabolite is yet produced
        return res
    lab=""
    lab=("".join(s[i] for (i,c) in enumerate(frag) if c!="0")+
        "".join(cs[i] for (i,c) in enumerate(cfrag) if c!="0"))
    flab=frag.replace("0", "")+cfrag.replace("0", "")
    for (pm,ps) in prods:
        fp="".join((flab[lab.index(l)] if l in lab else "0") for (i,l) in enumerate(ps))
        res.add((pm, fp))
    #print "res=", res
    return res
def cumo_iw(w,nlen):
    """iterator for a given cumomer weight w in the carbon length nlen"""
    #print ("cumo_iw: w=%d,nlen=%d\n" % (w,nlen));##
    if w < 0 or w > nlen:
       raise Exception("CumoWeightOutOfRange")
    if nlen < 0:
       raise Exception("CumoLengthNegative")
    if w == 1:
        movbit=1
        for i in range(nlen):
            yield movbit
            movbit<<=1
    else:
        movbit=1<<(w-1)
        for i in range(nlen-w+1):
            for subi in cumo_iw(w-1,w+i-1):
                yield movbit+subi
            movbit<<=1

def iso2cumo(netan, strin, in_cumo, icumo, in_metab):
    """calculate cumomer fraction from isotopomer ones"""
    #return sum(iso_dic.get(iiso,0.)
    #    for iiso in icumo2iiso(icumo, Clen))
    w=sumbit(icumo)
    n=len(netan["iso_input"])
    if len(netan[strin]) != n:
        # init cumo_input list of dicts
        netan[strin]=[{} for i in range(n)]
    for ili in range(len(netan["iso_input"])):
        d=netan["iso_input"][ili][in_metab]
        netan[strin][ili][in_cumo]=sum(val for (iso, val) in d.items() if sumbit(iso&icumo) == w)

def iso2emu(netan, inmetab, mask, mpi, e):
    """calculate emu fraction from isotopomer dict iso_input.
    The fraction corresponds to a fragment defined by a mask and the mass component mpi.
    Return a real number in [0; 1] interval.
    """
    #return(sum([val for (frag, val) in iso_input[inmetab].iteritems() if sumbit(frag&mask) == mpi]))
    n=len(netan["iso_input"])
    if len(netan["emu_input"]) != n:
        # init emu_input list of dicts
        netan["emu_input"]=[{} for i in range(n)]
    for ili in range(n):
        d=netan["iso_input"][ili][inmetab]
        netan["emu_input"][ili][e]=sum([val for (frag, val) in d.items() if sumbit(frag&mask) == mpi])

def formula2dict(f, pterm=re.compile(r'([+-])'), pflux=re.compile(r'(?P<coef>\d+\.?\d*|^)?\s*\*?\s*(?P<var>[a-zA-Z_[\]()][\w\. -\[\]]*)\W*')):
    """parse a linear combination sum([+|-][a_i][*]f_i) where a_i is a 
    positive number and f_i is a string starting by non-digit and not white
    character (# is allowed). Output is a dict f_i:[+-]a_i"""
    res=dict()
    sign=1
    l=(i.strip() for i in pterm.split(str(f)))
    for term in l:
        try:
            next_sign=next(l)
            next_sign=-1 if next_sign == "-" else 1
        except StopIteration:
            next_sign=0
            pass
        if (len(term) == 0):
            continue
        m=pflux.match(term)
        if m:
            coef=m.group("coef")
            coef=1. if coef==None or not len(coef) else float(coef)
            var=m.group("var")
            sign=-1 if sign==-1 else 1
            res[var]=res.get(var, 0.)+sign*coef
            sign=next_sign
        else:
            raise Exception("Term '"+term+"' in formula '"+f+"' cannot be parsed.")
    return res
def mecoparse(terms, pmeco=re.compile(r'\s*((?P<coef>\d+\.?\d*|^)\s*\*\s*)?(?P<metab>[^*+ ]*)\s*$')):
    """Parse a string term from a list (or a sing string) of chemical equation entries.
    The general form of each term is 'coef*metab'.
    coef (if present) must be separated from metab by '*' and be convertible to float.
    metab can start with a number (e.g. '6PG') so the presence of '*' is mandatory to
    separate coef from metab.If coef is absent, it is considered to be 1. 
    Return a list of (or a single for str) tuples (metab (str), coef (real)).
    """
    single=isstr(terms)
    if single:
        terms=[terms]
    res=[]
    for term in terms:
        m=pmeco.match(term)
        if not m:
            raise Exception("Wrong format for term '%s'"%term)
        coef=m.group("coef")
        coef=1. if coef==None or not len(coef) else float(coef)
        metab=m.group("metab")
        if not metab:
            raise Exception("Metabolite cannot be empty in term '%s'"%term)
        res.append((metab, coef))
    if single:
        res=res.pop()
    return res

def label_meas2matrix_vec_dev(netan):
    """use netan["label_meas"] list to construct a corresponding list of
    measure matrix matx_lab
    such that scale_diag*metab_pool_diag*matx_lab*(cumos_vector,1) corresponds to label_measurements_vector.
    matx_lab is defined as
    list of dict{"scale":scale_name, "coefs":dict{icumo:coef}, "metab": metabolite, "poolid": metabolite pool id if pooled}
    where coef is a contribution of cumo in linear combination for given measure.
    scale_name is of the form "metabs;group". Group number is to group
    measurements of the same measurement set.
    poolid is the index of pool list in pooled where each list regroups
    0-based indexes rows in returned matrix for what has to be pooled together.
    vec is a list of measurements (values in .ftbl)
    dev is a list of deviations.
    Elements in matx_lab, vec and dev are ordered in the same way.
    The returned result is a dict (mat,vec,dev)
    """
    # lab[metab][group]=list of {id:txt, val:x, dev:y, bcumos:list of #bcumo}
    # bcumo is like #1xx01 (0,1 and x are allowed)
    # their sum have to be transformed in cumomer linear combination
    nexp=len(netan["iso_input"])
    mli=list(range(nexp))
    for ili in range(nexp):
        mat=[]
        vec=[]
        dev=[]
        pooled=[]
        ids=[]
        imrow=-1
        for (metabs,groups) in netan["label_meas"][ili].items():
            for (group,rows) in groups.items():
                for row in rows:
                    vec.append(row["val"])
                    dev.append(row["dev"])
                    ids.append(row["id"])
                    res=dict()
                    for cumostr in row["bcumos"]:
                        decomp=bcumo_decomp(cumostr)
                        for icumo in decomp["+"]:
                            res.setdefault(icumo,0)
                            res[icumo]+=1
                        for icumo in decomp["-"]:
                            res.setdefault(icumo,0)
                            res[icumo]-=1
                    emuco=dict((str(ic)+"+"+str(sumbit(ic)), c) for (ic, c) in res.items())
                    if len(row["pooled"]) > 1:
                        # init index list
                        pooled.append([row["id"]])
                    for metab in row["pooled"]:
                        imrow+=1
                        mat.append({"id": row["id"], "scale": metabs+";"+group, "coefs":res,
                            "bcumos": row["bcumos"], "metab": metab, "emuco": emuco})
                        if len(row["pooled"]) > 1:
                            # indeed something is pooled
                            pooled[-1].append(imrow)
        mli[ili]={"ids": ids, "mat": mat, "vec": vec, "dev": dev, "pooled": pooled}
    return mli

def mass_meas2matrix_vec_dev(netan):
    """use netan["mass_meas"] list to construct a corresponding list of
    measure matrix matx_mass
    such that scale_diag*matx_mass*cumos_vector corresponds to mass_measures_vector.
    matx_mass is defined as matx_lab in label_meas2matrix_vec_dev()
    Elements in matx_mass, vec and dev are ordered in the same way.
    scale name is defined as "metab;fragment_mask"
    The returned result is a dict (mat,vec,dev)
    """
    emu=netan["emu"]
    # mass[metab][frag_mask][weight]={val:x, dev:y}
    # weight has to be transformed in cumomer linear combination
    nexp=len(netan["iso_input"])
    mli=list(range(nexp))
    for ili in range(nexp):
        mat=[]
        vec=[]
        dev=[]
        pooled=[]
        ids=[]
        imrow=-1
        for (m_id,frag_masks) in netan["mass_meas"][ili].items():
            #print "mass matx calc for ", metab;##
            (l, metabs, frag, m_irow)=m_id.split(":")
            for (fmask,weights) in frag_masks.items():
                onepos=[b_no for (b_no,b) in iternumbit(fmask) if b]
                nmask=sumbit(fmask)
    #            aff("weights for met,fmask "+", ".join((metab,strbit(fmask))), weights);##
                for (weight,row) in weights.items():
                    vec.append(row["val"])
                    dev.append(row["dev"])
                    ids.append(row["id"])
                    metabl=row["pooled"]
                    metab0=metabl[0]
                    n=netan["Clen"][metab0]
                    str0="0"*n
                    strx="x"*n
                    str1="1"*n
                    fmask0x=setcharbit(strx,"0",fmask)
                    fmask01=setcharbit(str0,"1",fmask)
                    if not emu:
                        # for a given weight construct bcumo sum: #x10x+#x01x+...
                        bcumos=["#"+setcharbit(fmask0x,"1",expandbit(iw,onepos))
                                for iw in range(1<<nmask) if sumbit(iw)==weight]
        #                aff("bcumos for met,fmask,w "+", ".join((metab,strbit(fmask),str(weight))), [b for b in bcumos]);##
                        res=dict()
                        for cumostr in bcumos:
                            #print cumostr;##
                            decomp=bcumo_decomp(cumostr)
                            for icumo in decomp["+"]:
                                res.setdefault(icumo,0)
                                res[icumo]+=1
                            for icumo in decomp["-"]:
                                res.setdefault(icumo,0)
                                res[icumo]-=1
                    else:
                        bcumos=None
                        res=None
                    emuco={str(fmask)+"+"+str(weight): 1}
                    if len(row["pooled"]) > 1:
                        # init index list
                        pooled.append([row["id"]])
                    for metab in row["pooled"]:
                        imrow+=1
                        if len(row["pooled"]) > 1:
                            # indeed something is pooled
                            pooled[-1].append(imrow)
                        mat.append({"id": row["id"], "scale": metabs+";"+fmask01+";"+m_irow, "coefs":res,
                            "bcumos": bcumos, "metab": metab, "emuco": emuco})
        mli[ili]={"ids": ids, "mat": mat, "vec": vec, "dev": dev, "pooled": pooled}
    return mli

def peak_meas2matrix_vec_dev(netan, dmask={"S": 2, "D-": 6, "D+": 3, "T": 7, "DD": 7}):
    """use netan["peak_meas"] list to construct a corresponding list of
    measure matrix matx_peak
    such that scale_diag*matx_peak*cumos_vector corresponds to
    peak_measures_vector.
    dmask is a dictionary with 3 carbon labeling pattern mask for
    various peak types. The middle bit corresponds to the targeted carbon,
    lower bit corresponds to the next neighbour (D+) and higher bit
    corresponds to previous carbon (D-).
    matx_peak is defined as matx_lab in label_meas2matrix_vec_dev()
    Elements in matx_peak, vec and dev are ordered in the same way.
    scale name is defined as "metab;c_no;irow"
    The returned result is a dict (mat,vec,dev)
    """
    # peak[metab][c_no][peak_type]={val:x, dev:y}
    # c_no+peak_type have to be transformed in cumomer linear combination
    # c_no is 1-based and left (just after # sign) started.
    nexp=len(netan["iso_input"])
    mli=list(range(nexp))
    for ili in range(nexp):
        mat=[]
        vec=[]
        dev=[]
        pooled=[]
        ids=[]
        imrow=-1
        for (metabs,irows) in netan["peak_meas"][0].items():
            #print "peak matx calc for ", metab;##
            #import pdb; pdb.set_trace()
            metabl=[*[*irows.values()][0].values()][0]["pooled"]
            #print(metabl)
            #metabl=metabl["pooled"]
            
            metab0=metabl[0]
            n=netan["Clen"][metab0]
            strx="#"+"x"*n
            for (irow,peaks) in irows.items():
                # pmask0x put 0 on three targeted carbons and leave x elsewhere
                c_no=[*peaks.values()][0]["c_no"]
                pmask0x=setcharbit(strx,"0",1<<(n-c_no))
                if c_no > 1:
                    # add left neighbour
                    pmask0x=setcharbit(pmask0x,"0",1<<(n-c_no+1))
                if c_no < n:
                    # add right neighbour
                    pmask0x=setcharbit(pmask0x,"0",1<<(n-c_no-1))
    #            aff("peaks for met,c_no,pmask0x="+", ".join((metab,str(c_no),pmask0x)), peaks);##
                for (peak,row) in peaks.items():
                    vec.append(row["val"])
                    dev.append(row["dev"])
                    ids.append(row["id"])
                    res=dict()
                    # for a given (c_no,peak) construct bcumo sum: #x10x+#x01x+...
                    # shift the 3-bit mask to the right carbon position
                    if c_no < n:
                        pmask=dmask[peak]<<(n-c_no-1)
                    else:
                        pmask=dmask[peak]>>1
                    bcumos=[setcharbit(pmask0x,"1",pmask)]
                    if peak=="D-" and "T" in peaks:
                        # D+===D- => add D+ to the list of measured bcumomers
                        # set the 3-bit mask to the right carbon position
                        if c_no < n:
                            pmask=dmask["D+"]<<(n-c_no-1)
                        else:
                            pmask=dmask["D+"]>>1
                        bcumos.append(setcharbit(pmask0x,"1",pmask))
    #                aff("bcumos for met,fmask,w "+", ".join((metab,strbit(fmask),str(weight))), [b for b in bcumos]);##
                    for cumostr in bcumos:
                        #print cumostr;##
                        decomp=bcumo_decomp(cumostr)
                        for icumo in decomp["+"]:
                            res.setdefault(icumo,0)
                            res[icumo]+=1
                        for icumo in decomp["-"]:
                            res.setdefault(icumo,0)
                            res[icumo]-=1
                    emuco=dict((str(ic)+"+"+str(sumbit(ic)), c) for (ic, c) in res.items())
                    if len(row["pooled"]) > 1:
                        # init index list
                        pooled.append([row["id"]])
                    for metab in row["pooled"]:
                        imrow+=1
                        mat.append({"id": row["id"], "scale": metabs+";"+str(c_no)+";"+row["irow"], "coefs":res,
                            "bcumos": bcumos, "metab": metab, "emuco": emuco})
                        if len(row["pooled"]) > 1:
                            # indeed something is pooled
                            pooled[-1].append(imrow)
                    #print "metab,c_no,peak="+",".join((metab,str(c_no),str(peak)));##
                    #print "res=",str(res);##
        mli[ili]={"ids": ids, "mat": mat, "vec": vec, "dev": dev, "pooled": pooled}
    return mli

def bcumo_decomp(bcumo):
    """bcumo is a string of the form #[01x]+. It has to be decomposed
    in the linear combination of cumomers #[1x]+. The coefficients
    of this linear combination are 1 or -1. So it can be represented as
    sum(cumos_positive)-sum(cumos_negative). The result of this function
    is a dictionary {"+": list of icumos, "-": list of icumos}. icumo is
    an integer whose binary form indicates 1's positions in a cumomer."""
    
    # take zero positions
    zpos=[ i for i in range(len(bcumo)-1) if bcumo[-i-1]=="0" ]
    
    # basic 1's cumomer mask
    icumo_base=0
    for i in range(len(bcumo)-1):
        if bcumo[-i-1]=="1":
            icumo_base|=1<<i
    
    # run through all 2^n zero masks separating positive and negative
    # cumomers
    res={"+": [], "-": []}
    nz=len(zpos)
    for zmask in range(1<<nz):
        # odd or even bit number?
        sign="-" if sumbit(zmask)%2 else "+"
        
        # get full cumomer mask
        icumo=icumo_base|expandbit(zmask,zpos)
        # put in result dict
        res[sign].append(icumo)
    #print bcumo+"=" + "+".join(setcharbit("x"*len(bcumo),"1",i) for i in res["+"])+("-" if res["-"] else "") +"-".join(setcharbit("x"*len(bcumo),"1",i) for i in res["-"]);##
    return res
def mat2graph(A, fp):
    """write digraph file on file pointer fp representing
    links in matrix A given as bi-level dictionnary. A key of
    first level (row index) is influenced by keys of second level
    (column indicies)."""
    fp.write("digraph A {\n")
    for i in A:
        labi=i.replace(":", "_")
        for j in A[i]:
            if i==j:
                continue
            labj=j.replace(":", "_")
            fp.write("   \""+labj+"\" -> \""+labi+"\"\n")
    fp.write("}\n")
def dom_cmp(A, i, j):
    """Compares influances of i-th and j-th lements of A.
    Returns 0 if i and j are mutually influenced, 1 if i in A[j]
    (i influences j) , -1 otherwise"""
    return (0 if i in A[j] and j in A[i] else 1 if i in A[j] else -1 if j in A[i] else 0)
def mat2pbm(A, v, fp):
    """Write an image map of non-zero entries of matrix A to file pointer fp.
    Matrix A is a dictionnary, v is a list ordering keys of A."""
    fp.write("P1\n%d %d\n"% (len(A), len(A)))
    for i in v:
        p=0
        for j in v:
            fp.write("1 " if j in A[i] else "0 ")
            p+=2
            if p >= 69:
                fp.write("\n")
                p=0
        fp.write("\n")
        p=0
def aglom_loop1(A):
    """Agglomerate nodes of A if they mutually influence each other
    i.e. they are in a loop of length 1.
    Return a new dictionary of influence where entries are those of A aglomerated and glued "by" tab symbol"""
    # i->loop_name(=min of all nodes in this loop)
    na=copy.deepcopy(A)
    loop=dict()
    while aglom(na,transpose(na),loop):
        pass
    return({"na": na, "loop":loop})
def aglom(na,ta,loop):
    """new matrix A (na), transpose A (ta) are used to
    aglomerate neigbour mutually influencing nodes in a supernode.
    Aglomerated noeds are put in the loop dictionnary.
    Return False if no nodes were aglomerated.
    """
    found=False
    for i in na:
        for j in na[i]:
            if j==i:
                continue
            if i in na[j]:
                # agglomerate i and j in na
                # which means keep lowest netween i and j
                # lo=min(i,j)
                # hi=max(i,j)
                # and imports influences of hi to lo
                found=True
                lo=min(i,j)
                hi=max(i,j)
                del(na[hi][hi])
                del(ta[hi][hi])
                #print "aglom: "+i+"+"+j+"->"+lo;##
                #print "elim row "+hi+str(na[hi].keys());##
                #print "influenced rows "+str(ta[hi].keys());##
                # update na rows influenced by hi
                for ii in [*ta[hi].keys()]:
                    #print "\nold row "+ii+str(na[ii].keys());##
                    na[ii].update(na[hi])
                    del(na[ii][hi])
                    for jj in na[ii]:
                        ta[jj][ii]=na[ii][jj]
                    #print "new row "+ii+str(na[ii].keys());##
                # remove na[hi], ta[hi] and ta's corresponding to na[hi]
                for ii in na[hi]:
                    del(ta[ii][hi])
                del(na[hi])
                del(ta[hi])
                # update loop dictionary
                loop[lo]=loop.get(lo,oset((lo,)))
                loop[lo].update(loop.get(hi, oset((hi,))))
                if hi in loop:
                    del(loop[hi])
                #print "loop:"+str(loop)
                return(found)
def lowtri(A):
    """Try low triangular ordering of matrix A entries"""
    unsrt=[*A.keys()]
    srt=list()
    # first get inputs (no influences on them)
    srt=[k for k in unsrt if len(A[k])==1]
    for k in srt:
        unsrt.remove(k)
    # now move to lower number keys that influence others
    while unsrt:
        for k in list(unsrt):
            srt.extend(oset(A[k].keys()).difference(oset(srt)))
            for i in [*A[k].keys()]:
                try:
                    unsrt.remove(i)
                except:
                    pass
    return srt
def topo_order(A, tA):
    """Try to sort keys of A in topological order. tA is just a transpose of A"""
    unsrt=oset(A.keys())
    srtin=list()
    srtout=list()
    
    while unsrt:
        # shave inputs and outputs
        inps=oset(k for k in unsrt if len(oset(A[k]).difference(srtin))==1)
        srtin.extend(inps)
        unsrt.difference_update(inps)
        outs=oset(k for k in unsrt if len(oset(tA[k]).difference(unsrt))==1)
        if outs:
            srtout.insert(0, outs)
            unsrt.difference_update(outs)
        if not inps and not outs:
            # this is not a dag
            srtin.extend(unsrt)
            #print "not a dag"
            break
    return(srtin+srtout)
def transpose(A):
    """Transpose a matrix defined as a dict."""
    tA=dict()
    for i in A:
        for j in A[i]:
            tA[j]=tA.get(j, dict())
            tA[j][i]=A[i][j]
    return(tA)
def rcumo_sys(netan, emu=False):
    """Calculate reduced cumomers or EMU systems A*x=b
    we start with observed cumomers (emus) of max weight
    and we include only needed involved cumomers (emus)
    A list of cumomer (emu) lists (by weight) is stored
    in netan["vrcumo"] (netan["vemu"])
    """
    # generate measurements dico if not yet done
    if "measures" not in netan:
        measures=dict()
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("%s_meas2matrix_vec_dev(netan)"%meas)
        netan["measures"]=measures
    measures=netan["measures"]

    # init rcumo_input
    n=len(netan["iso_input"])
    if len(netan["rcumo_input"]) != n:
        # init rcumo_input list of dicts
        netan["rcumo_input"]=[{} for i in range(n)]
    
    # get cumomers involved in measurements
    meas_cumos=oset()
    if emu:
        for meas in measures:
            for item in measures[meas]:
                for row in item["mat"]:
                    metab=row["metab"]
                    meas_cumos.update(metab+":"+i.split("+")[0] for i in [*row["emuco"].keys()] if i[-2:]!="+0")
    else:
        for meas in measures:
            for item in measures[meas]:
                for row in item["mat"]:
                    metab=row["metab"]
                    meas_cumos.update(metab+":"+str(icumo) for icumo in [*row["coefs"].keys()] if icumo != 0)
    # make list of observed weights
    weights=oset(sumbit(int(cumo.split(":")[1])) for cumo in meas_cumos)
    if not weights:
        netan["vrcumo"]=[]
        netan["rcumo2i0"]=dict()
        netan["rcumo_sys"]={"A": [], "b": []}
        return {"A": [], "b": []}
    
    maxw=max(weights)
    weights=list(range(1,maxw+1))
    weights.reverse()
    
    # cumomers to visit are stored in lists by weights
    to_visit=dict((w,[]) for w in weights)
    A=[{} for i in range(maxw)]; # store matrices by weight
    b=[{} for i in range(maxw)]; # store rhs by weight
    used=oset()

    # initialize to_visit, we'll stop when it's empty
    for cumo in meas_cumos:
        (m,w)=cumo.split(":")
        if m in netan["input"]:
            continue
        w=sumbit(int(w))
        to_visit[w].append(cumo)

    # run through the network starting with heaviest cumomers
    for w in weights:
        cumo_i=0
        while cumo_i < len(to_visit[w]):
            for cumo in to_visit[w]:
                cumo_i+=1
                if cumo in used:
                    continue
                # add this cumo to used
                used.add(cumo)
                (metab, icumo)=cumo.split(":")
                if metab in netan["output"]:
                    # no equation for output metabs
                    continue
                icumo=int(icumo)
                A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]})
                # get the influents to cumo of all weights: equal and lower
                # equals go to A and lowers go to b
                infl=cumo_infl(netan, cumo)
                #print("cumo=", cumo, "infl=", infl)
                for (incumo,fl,imetab,iinmetab) in infl:
                    if incumo==cumo:
                        # no equation as no transformation
                        continue
                    (inmetab, inw)=incumo.split(":")
                    inicumo=int(inw)
                    inw=sumbit(inicumo)
                    # input metabolites are to rhs, others are to visit
                    if inmetab not in netan["input"] and inw != 0 and incumo not in to_visit[inw]:
                        to_visit[inw].append(incumo)
                    if inmetab in netan["input"]:
                        if emu:
                            # tuple emu (mask, m+i, string)
                            for mask, mpi, e in ((inicumo, i, incumo+"+"+str(i)) for i in range(inw+1)):
                                iso2emu(netan, inmetab, mask, mpi, e)
                        if not netan["rcumo_input"] or incumo not in netan["rcumo_input"][0]:
                            iso2cumo(netan, "rcumo_input", incumo, inicumo, inmetab)
                        #netan["rcumo_input"][incumo]=netan["cumo_input"][incumo]
                    # main part: write equations
                    if inw==w :
                        # equal weight => A
                        if inmetab in netan["input"]:
                            b[w-1][cumo]=b[w-1].get(cumo, dict())
                            b[w-1][cumo][fl]=b[w-1][cumo].get(fl, dict())
                            b[w-1][cumo][fl][imetab]=b[w-1][cumo][fl].get(imetab,[])
                            #b[w-1][cumo][fl][imetab].append(netan["cumo_input"][incumo])
                            b[w-1][cumo][fl][imetab].append(incumo)
                            continue
                        A[w-1][cumo][incumo]=A[w-1][cumo].get(incumo,[])
                        A[w-1][cumo][incumo].append(fl)
                        #aff("A "+str(w)+cumo, A[w-1][cumo]);##
                    elif inw < w:
                        # lower weight => b
                        b[w-1][cumo]=b[w-1].get(cumo, dict())
                        b[w-1][cumo][fl]=b[w-1][cumo].get(fl, dict())
                        b[w-1][cumo][fl][imetab]=b[w-1][cumo][fl].get(imetab,[])
                        b[w-1][cumo][fl][imetab].append(incumo)
                # gather all influx in diagonal term
                #aff("A "+str(w)+" "+cumo, A[w-1][cumo]);##
                #aff("b "+str(w)+" "+cumo, b[w-1].get(cumo));##
                for incumo in A[w-1][cumo]:
                    if incumo == cumo:
                        continue
                    A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]})
                    #print("adding A:", incumo, "->", cumo, A[w-1][cumo][incumo]);##
                    A[w-1][cumo][cumo]+=A[w-1][cumo][incumo]
                if cumo in b[w-1]:
                    A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]})
                    #print("adding b:", cumo, b[w-1][cumo].keys());##
                    A[w-1][cumo][cumo]+=[*b[w-1][cumo].keys()]
    #aff("to_v", to_visit);##
    # make ordered list for reduced cumomer set
    netan["vrcumo"]=[[*a.keys()] for a in A]
    netan["rcumo2i0"]=dict((cumo,i) for (i, cumo) in enumerate(valval(netan["vrcumo"])))
    #print("A=", A)
    netan["rcumo_sys"]={"A": A, "b": b}
    # make order list for emu_input
    if emu:
        netan["vemu_input"]=[*netan["emu_input"][0].keys()]
        # make order list for internal emus
        netan["vemu"]=copy.deepcopy(netan["vrcumo"])
        for w in range(len(netan["vemu"])):
            l=netan["vemu"][w]
            netan["vemu"][w]=[m+"+"+str(i) for i in range(w+2) for m in l]
        netan["emu2i0"]=dict((emu,i) for (i, emu) in enumerate(valval(netan["vemu"])))
    return netan["rcumo_sys"]

def cumo_infl(netan, cumo):
    """cumo_infl(netan, cumo)->list(tuple(in_cumo, fl, imetab, iin_metab))
    return the list of tuples (in_cumo, fl, imetab, iin_metab):
    input cumomer, flux (fwd.fl or rev.fl), index of metab and index of in_metab
    generating cumo. cumo is in format "metab:icumo".
    Condenstation reaction will give the same flux and icumo but various
    iin_metab.
    Convergent point will give multiple fluxes.
    """
    (metab, icumo)=cumo.split(":")
    icumo=int(icumo)
    res=[]
    # run through input forward fluxes of this metab
    for reac,coef in oset(netan["sto_m_r"][metab]["right"]):
        if reac not in netan["carbotrans"]:
            continue # it can happen because of NOTRACER_NETWORK
        # get all cstr for given metab
        for (imetab,cstr) in ((i,s) for (i,(m,s)) in enumerate(netan["carbotrans"][reac]["right"])
            if m==metab):
            # get all input cumomer in this reaction for this cstr
            for (iin_metab, (in_metab,in_str)) in enumerate(netan["carbotrans"][reac]["left"]):
                in_icumo=src_ind(in_str, cstr, icumo)
                if in_icumo != None:
                    in_cumo=in_metab+":"+str(in_icumo)
                    res.append((in_cumo, "fwd."+reac, imetab, iin_metab))
    # run through input reverse fluxes of this metab
    fluxset=oset(f for f,c in netan["sto_m_r"][metab]["left"])
    if (clownr):
        # non reversible reactions are all positive
        fluxset=fluxset.difference(netan["notrev"])
    else:
        # non reversible reactions can change sens => keep reverse flux just in case
        fluxset=fluxset.difference(netan["flux_inout"])
    for reac in fluxset:
        # get all cstr for given metab
        for (imetab,cstr) in ((i,s) for (i,(m,s)) in enumerate(netan["carbotrans"][reac]["left"])
            if m==metab):
            # get all input cumomer in this reaction for this cstr
            for (iin_metab, (in_metab,in_str)) in enumerate(netan["carbotrans"][reac]["right"]):
                in_icumo=src_ind(in_str, cstr, icumo)
                if in_icumo != None:
                    in_cumo=in_metab+":"+str(in_icumo)
                    res.append((in_cumo, "rev."+reac, imetab, iin_metab))
    return res
def infl(metab, netan):
    """infl(metab, netan)->oset(fluxes)
    List incoming fluxes for this metabolite (fwd.reac, rev.reac, ...)
    """
    # run through input forward fluxes of this metab
    res=oset("fwd."+reac for reac,_ in netan["sto_m_r"][metab]["right"])
    # run through input reverse fluxes of this metab
    res.update("rev."+reac for reac in
        oset(f for f,_ in netan["sto_m_r"][metab]["left"]).difference(netan["flux_inout"]))
#        oset(netan["sto_m_r"][metab]["left"]).difference(netan["notrev"]))
    return(res)

def t_iso2m(n):
    """t_iso2m(n) return transition matrix from isotopomers fractions to MID vector
    n - carbon number
    return numpy array of size (n+1,2**n)
    """
    # isotopomer number
    ni=2**n
    return np.array(
        [[1. if sumbit(ii)==im else 0 for ii in range(ni)] for im in range(n+1)]
    )

def t_iso2cumo(n):
    """t_iso2cumo(n) return transition matrix from isotopomers fractions to cumomer vector
    n - carbon number
    return numpy array of size (2**n,2**n)
    """
    if n <= 0:
        return np.array(1, ndmin=2)
    # recursive call
    m_1=t_iso2cumo(n-1)
    nc1,nc1=m_1.shape
    nc=2*nc1
    m=np.zeros((nc,nc))
    m[0:nc1,0:nc1],m[0:nc1,nc1:],m[nc1:,nc1:]=m_1,m_1,m_1
    return m

def t_iso2pos(n):
    """t_iso2pos(n) return transition matrix from isotopomers fractions to positional
    labelling vector (cumomers of weight 1)
    n - carbon number
    return numpy array of size (n,2**n)
    """
    m=t_iso2cumo(n)
    return m[[1<<i for i in range(n)],:]

def conv_mid(x,y):
    """conv_mid(x,y)->z
    convolute two mid vectors (numpy arrays)
    and return the result as numpy array.
    """
    nx=len(x)
    ny=len(y)
    z=np.zeros(nx+ny-1)
    for i in range(ny):
        z[i:i+nx]+=x*y[i]
    return(z)

def ms_frag_gath(netan):
    """gather metabolite fragments necessary to obtain a given set of data
    observed in MS measurements.
    The fragment mask is encoded in the same way as cumomers, Met:7 <=> Met#(0)111
    """
    frags=oset()
    to_visit=oset()
    for di in netan["mass_meas"]:
        for m_id in di:
            for mask in di[m_id]:
                # just the first weight item is sufficient
                item=[*di[m_id][mask].values()][0]
                to_visit.update(met+":"+str(mask) for met in item["pooled"])
    while to_visit:
        for frag in list(to_visit):
            if frag in frags:
                # already accounted
                to_visit.remove(frag)
                continue
            frags.add(frag)
            to_visit.remove(frag)
            # add its contributors for visiting
            to_visit.update(incumo for (incumo,fl,imetab,iinmetab) in cumo_infl(netan, frag))
            break
    return(frags)
def ntimes(n):
    """Return charcater string 'once' for n=1, 'twice' for n=2 and 'n times' for other n"""
    return("once" if n==1 else "twice" if n==2 else "%d times"%n)
def proc_label_input(ftbl, netan, case_i=False):
    """Proceed LABEL_INPUT section in ftbl and add result to the list netan["iso_input"] and netan["funlab"] (case_i)
    List item is a dict {}metab;{isotop_int_index:fraction}}
    """
    ili=len(netan["iso_input"]) # label_input list index
    netan["iso_input"]+=[{}]
    if case_i:
        netan["funlab"]+=[{}]
        resf=netan["funlab"][ili]
    res=netan["iso_input"][ili]
    # input isotopomers
    for row in ftbl.get("LABEL_INPUT",[]):
        metab=row.get("META_NAME", "") or metab
        if metab not in netan["Clen"]:
            raise Exception("Label input metabolite `%s` is not defined in NETWORK (%s: %s)."%(metab, ftbl["name"], row["irow"]))
        ilen=len(row.get("ISOTOPOMER", ""))-1; # -1 for '#' sign
        if ilen != netan["Clen"][metab]:
            raise Exception("Input isotopomer `%s` is of bad length (%d). A length of %d is expected (%s: %s)."%
                (row.get("META_NAME", "")+row.get("ISOTOPOMER", ""), ilen,  netan["Clen"][metab], ftbl["name"], row["irow"]))
        iiso=strbit2int(row["ISOTOPOMER"])
        res[metab]=res.get(metab, {})
        if case_i:
            resf[metab]=resf.get(metab, {})
            resf[metab][iiso]=row["VALUE"]
            res[metab][iiso]=NA
        else:
            try:
                val=eval(row["VALUE"])
            except:
                val=row["VALUE"]
            if type(val) in (float, int):
                if val < 0. and val >= -tol:
                    val=0.
                if val > 1. and val <= 1.+tol:
                    val=1.
                if val < 0 or val > 1:
                    raise Exception("Input isotopomer `%s` has a value (%g) out of range [0; 1] (%s: %s)"%(row.get("META_NAME", "")+row.get("ISOTOPOMER", ""), val, ftbl["name"], row["irow"]))
            res[metab][iiso]=val
    # check that all isoforms sum up to 1 for all inputs for influx_s
    if not case_i:
        for metab in res:
            le=len(res[metab])
            nfo=2**netan["Clen"][metab]
            su=sum(res[metab].values())
            if su > 1.+1.e-10:
                raise Exception("Input metabolite `%s` has label summing up to %g which is greater than 1 (%s)."%(metab, su, ftbl["name"]))
            if le == nfo:
                # all forms are given => must sum up to 1
                if abs(su-1.) > 1.e-10:
                    raise Exception("Input metabolite `%s` has label summing up to %g instead of 1 (%s)."%(metab, su, ftbl["name"]))
            elif le < nfo-1:
                # many forms are lacking (not just one)
                # if fully unlabeled is lacking, use it to complete to 1
                # otherwise raise an error
                if 0 not in res[metab]:
                    res[metab][0]=1.-su
                elif abs(su-1.) > 1.e-10:
                    raise Exception("Input metabolite `%s` has lacking labels to sum up to 1 (we get sum=%g instead) (%s)."%(metab, su, ftbl["name"]))
            elif le == nfo-1:
                # just one form is lacking
                if su != 1.:
                    # add it to complete to 1
                    la=list(oset(range(nfo))-oset(res[metab]))[0]
                    res[metab][la]=1.-su
    if ili > 0:
        # complete absent inputs by their values from the first ftbl
        for m in oset(netan["iso_input"][0])-oset(res):
            res[m]=netan["iso_input"][0][m];
        if case_i:
            for m in oset(netan["funlab"][0])-oset(resf):
                resf[m]=netan["funlab"][0][m];
    zeroc=oset(m for m,n in netan["Clen"].items() if n == 0)
    lmi=oset(res.keys())-netan["input"] # label input - input
    iml=oset(netan["input"])-zeroc-oset(res.keys()) # input - label input
    if lmi:
        raise Exception("LABEL_INPUT section contains metabolite(s) that are not network input(s): '%s' (%s)\n"%(", ".join(lmi), ftbl["name"]))
    if iml:
        raise Exception("LABEL_INPUT section lacks certain metabolite(s) that are network input(s): '%s' (%s)\n"%(", ".join(iml), ftbl["name"]))

def proc_label_meas(ftbl, netan):
    """Proceed LABEL_MEASUREMENT section of ftbl file, add the result to a list of dicts
    """
    # label measurements
    # [ili]][metab][group]->list of {val:x, dev:y, bcumos:list of #bcumo}
    ili=len(netan["label_meas"])
    netan["label_meas"]+=[{}]
    res=netan["label_meas"][ili]
    
    metabs=""
    for row in ftbl.get("LABEL_MEASUREMENTS", []):
        #print row;##
        # test the cumomer pattern validity
        if (not re.match(r"#[01x]+(\+#[01x]+)*", row["CUM_CONSTRAINTS"])):
            raise Exception("Not valid cumomer's pattern in '"+row["CUM_CONSTRAINTS"]+"' (%s: %s)"%(ftbl["name"], row["irow"]))
        metabs=row["META_NAME"] or metabs
        group=row["CUM_GROUP"] or group
        # metabs can be metab1[+metab2[+...]]
        metabl=metabs.split("+")
        if (len(metabl) > 1):
            # pooling metabolites will need their concentraions
            mdif=oset(metabl).difference(netan["met_pools"])
            if mdif:
                raise Exception("Pooled metabolite(s) '%s' are absent in METABOLITE_POOLS section (%s: %s)."%(join(", ", mdif), ftbl["name"], row["irow"]))

        # check that all metabs are unique
        count=dict((i, metabl.count(i)) for i in oset(metabl))
        notuniq=oset((i for i in count if count[i] > 1))
        if notuniq:
            raise Exception("Metabolite(s) '%s' is present twice or more (%s: %s)"%(join(", ", notuniq), ftbl["name"], row["irow"]))

        metab0=metabl[0] if metabl else ""
        if not metab0 in netan["metabs"]:
            raise Exception("Unknown metabolite name '%s' in LABEL_MEASUREMENTS (%s: %s)"%(metab0, ftbl["name"], row["irow"]))
        clen0=netan["Clen"][metab0] if metab0 else 0
        for metab in metabl:
            if not metab in netan["metabs"]:
                raise Exception("Unknown metabolite name '%s' in LABEL_MEASUREMENTS (%s: %s)"%(metab, ftbl["name"], row["irow"]))
            if metab in netan["output"] or metab in netan["input"]:
                raise Exception("""Measured metabolites have to be internal to network (found in input or output metabolites of '%s').
In case of output, you can add a fictious metabolite in your network immediatly after '"""%ftbl["name"]+metab+"' (seen in LABEL_MEASUREMENTS).")
            mlen=netan["Clen"][metab]
            clen=netan["Clen"][metab]
            if clen!=clen0 :
                raise Exception("Carbon length of '%s' (%d) is different from the length of '%s' (%d) (%s: %s)"%(metab, clen, metab0, clen0, ftbl["name"], row["irow"]))
        if not metabs in res:
            res[metabs]=dict()
        if not group in res[metabs]:
            res[metabs][group]=[]
        # prepare cumos list
        # if bcumos is not empty use it
        # else use group name as carbon number (starting from # character)
        if row.get("CUM_CONSTRAINTS",""):
            bcumos=row["CUM_CONSTRAINTS"].split("+")
        else:
            try:
                # just put "1" in group-th place
                i=int(group)
                bcumos="#"+"x"*netan["Clen"][metabl[0]]
                bcumos[i]="1"
                bcumos=[bcumos]
            except:
                raise Exception("Expected integer CUM_GROUP in LABEL_MEASUREMENTS on row "+row["irow"])
        if not row["VALUE"]:
            raise Exception("The field VALUE is empty (%s: %s)"%(ftbl["name"], row["irow"]))
        if not row["DEVIATION"]:
            raise Exception("The field DEVIATION is empty (%s: %s)"%(ftbl["name"], row["irow"]))
        try:
            val=float(eval(row["VALUE"]))
        except:
            val=NaN
        try:
            sdev=float(eval(row["DEVIATION"]))
        except:
            raise Exception("DEVIATION must evaluate to a real positive number (%s: %s)."%(ftbl["name"], row["irow"]))
        if sdev <= 0.:
            raise Exception("DEVIATION must be positive (%s: %s)."%(ftbl["name"], row["irow"]))
        res[metabs][group].append({
                "val":val,
                "dev":sdev,
                "bcumos":row["CUM_CONSTRAINTS"].split("+"),
                "id":":".join(["l", metabs, row["CUM_CONSTRAINTS"], str(row["irow"])]),
                "pooled":metabl,
        })
        # test the icumomer lengths
        if not all(len(ic)==mlen+1 for ic in 
                res[metabs][group][-1]["bcumos"]):
            raise Exception("Wrong cumomer length for %s in %s (%s: %s)"%(metab, row["CUM_CONSTRAINTS"], ftbl["name"], row["irow"]))
def proc_peak_meas(ftbl, netan):
    """Proceed PEAK_MEASUREMENT section of ftbl file, add the result to a list of dicts
    """
    # peak measurements
    # [ili]][metab][c_no][peak_type in (S,D-,D+,(DD|T))]={val:x, dev:y}
    ili=len(netan["peak_meas"])
    netan["peak_meas"]+=[{}]
    res=netan["peak_meas"][ili]
    # peak measurements
    metabs=""
    for row in ftbl.get("PEAK_MEASUREMENTS",[]):
        #print row;##
        # test the pattern validity
        if (row.get("VALUE_DD","") and row.get("VALUE_T","")):
            raise Exception("Not valid value combination. Only one of DD and T has to be in row "+str(row)+" (%s: %s)"%(ftbl["name"], row["irow"]))
        metabs=row["META_NAME"] or metabs
        metabl=metabs.split("+")
        if (len(metabl) > 1):
            # pooling metabolites will need their concentraions
            mdif=oset(metabl).difference(netan["met_pools"])
            if mdif:
                raise Exception("Pooled metabolite(s) '%s' are absent in METABOLITE_POOLS section (%s: %s)"%(join(", ", mdif), ftbl["name"], row["irow"]))

        # check that all metabs are unique
        count=dict((i, metabl.count(i)) for i in oset(metabl))
        notuniq=oset((i for i in count if count[i] > 1))
        if notuniq:
            raise Exception("Metabolite(s) '%s' is present twice or more (%s: %s)"%(join(", ", notuniq), ftbl["name"], row["irow"]))

        metab0=metabl[0] if metabl else ""
        if not metab0 in netan["metabs"]:
            raise Exception("Unknown metabolite name '%s' in PEAK_MEASUREMENTS (%s: %s)"%(metab0, ftbl["name"], row["irow"]))
        clen0=netan["Clen"][metab0] if metab0 else 0
        for metab in metabl:
            if not metab in netan["metabs"]:
                raise Exception("Unknown metabolite name '%s' in PEAK_MEASUREMENTS (%s: %s)"%(metab, ftbl["name"], row["irow"]))
            if metab in netan["output"] or metab in netan["input"]:
                raise Exception("""Measured metabolites have to be internal to network (seen in input or output metabolites (%s: %d) ).
In case of output, you can add a fictious metabolite in your network immediatly after '"""%(ftbl["name"], row["irow"])+metab+"' (seen in PEAK_MEASUREMENTS).")
            clen=netan["Clen"][metab]
            if clen!=clen0 :
                raise Exception("Carbon length of '%s' (%d) is different from the length of '%s' (%d) (%s: %s)"%(metab, clen, metab0, clen0, ftbl["name"], row["irow"]))
        res.setdefault(metabs, dict())
        for suff in ("S", "D-", "D+", "DD", "T"):
            # get val and dev for this type of peak
            val=row.get("VALUE_"+suff,"")
            if (not val or val=="-"):
                continue
            dev=row.get("DEVIATION_"+suff,"") or row.get("DEVIATION_S")
            if suff in ("DD", "T"):
               dev=row.get("DEVIATION_DD/T","") or row.get("DEVIATION_S")
            if not dev:
                raise Exception("The field DEVIATION_%s is empty (%s: %s)"%(suff, ftbl["name"], row["irow"]))
            try:
                val=float(eval(val))
            except:
                val=NaN
            try:
                dev=float(eval(dev))
            except:
                continue
            # test validity
            if not dev:
                raise Exception("Deviation is not valid for VALUE_"+suff+" (%s: %s)")%(ftbl["name"], row["irow"])
            c_no=int(row["PEAK_NO"])
            if c_no > clen0:
                raise Exception("Carbon number "+str(c_no)+" is greater than carbon length "+str(clen0)+" for metabolite '"+metab0+"' (%s: %s)"%(ftbl["name"], row["irow"]))
            if suff == "D-" and c_no == 1:
                raise Exception("Peak D- cannot be set for metabolite "+metab0+", c_no=1 (%s: %s)"%(ftbl["name"], row["irow"]))
            if suff == "D+" and c_no == clen:
                raise Exception("Peak D+ cannot be set for metabolite "+metab+", c_no="+str(c_no)+" (%s: %s)"%(ftbl["name"], row["irow"]))
            if (suff == "DD" or suff == "T") and (c_no == 1 or c_no == clen or clen < 3):
                raise Exception("Peak DD (or T) cannot be set for metabolite "+metab+", c_no="+str(c_no)+", len="+str(clen)+" (%s: %s)"%(ftbl["name"], row["irow"]))
            res[metabs].setdefault(row["irow"], dict())
            res[metabs][row["irow"]][suff]={
               "val": val,
               "dev": dev,
               "id": ":".join(["p", metabs, row["PEAK_NO"], suff, str(row["irow"])]),
               "irow": str(row["irow"]),
               "c_no": c_no,
               "pooled":metabl,
            }
def proc_mass_meas(ftbl, netan):
    """Proceed PEAK_MEASUREMENT section of ftbl file, add the result to a list of dicts
    """
    ili=len(netan["mass_meas"])
    netan["mass_meas"]+=[{}]
    res=netan["mass_meas"][ili]
    # mass measurements
    # [ili]][metab][frag_mask][weight]={val:x, dev:y}
    metabs=""
    for row in ftbl.get("MASS_SPECTROMETRY",[]):
        #print row;##
        metabs=row["META_NAME"] or metabs
        # metabs can be metab1[+metab2[+...]]
        metabl=metabs.split("+")
        if (len(metabl) > 1):
            # pooling metabolites will need their concentraions
            mdif=oset(metabl).difference(netan["met_pools"])
            if mdif:
                raise Exception("Pooled metabolite(s) '%s' are absent in METABOLITE_POOLS section (%s: %s)"%(join(", ", mdif), ftbl["name"], row["irow"]))

        # check that all metabs are unique
        count=dict((i, metabl.count(i)) for i in oset(metabl))
        #print(count)
        notuniq=oset((i for i in count if count[i] > 1))
        if notuniq:
            raise Exception("Metabolite(s) '%s' is present twice or more (%s: %s)"%(join(", ", notuniq), ftbl["name"], row["irow"]))

        metab0=metabl[0] if metabl else ""
        if not metab0 in netan["metabs"]:
            raise Exception("Unknown metabolite name '%s' in MASS_MEASUREMENTS (%s: %s)"%(metab0, ftbl["name"], row["irow"]))
        clen0=netan["Clen"][metab0] if metab0 else 0
        if row["META_NAME"]:
            irow=int(row["irow"])
            frag=""
            #frag=",".join(str(i+1) for i in range(clen0))
            mask=(1<<clen0)-1
        frag=row["FRAGMENT"] or frag
        for metab in metabl:
            clen=netan["Clen"][metab]
            # test the validity
            if not metab in netan["metabs"]:
                raise Exception("Unknown metabolite name '%s' in MASS_SPECTROMETRY (%s: %s)."%(metab, ftbl["name"], row["irow"]))
            if metab in netan["output"] or metab in netan["input"]:
                raise Exception("""Measured metabolites have to be internal to network.
In case of output, you can add a fictious metabolite following to '"""+metab+"' (seen in MASS_MEASUREMENTS).")
            clen=netan["Clen"][metab]
            if clen!=clen0 :
                raise Exception("Carbon length of '%s' (%d) is different from the length of '%s' (%d) (%s: %s)"%(metab, clen, metab0, clen0, ftbl["name"], row["irow"]))
        if row["FRAGMENT"]:
            # recalculate fragment mask
            mask=0
            for item in frag.split(","):
                try:
                    i=int(item)
                    if i > clen0:
                        raise Exception("An item '"+item+"' of carbon fragment is higher than metabolite '"+metab0+"' length "+str(clen0)+" (%s: %d)"%(ftbl["name"], irow))
                    # add this simple item to the mask
                    mask|=1<<(clen0-i)
                except ValueError:
                    # try the interval
                    try:
                        (start,end)=item.split("~")
                        #print "start,end=%s,%s" % (start,end);##
                        try:
                            for i in range(int(start),int(end)+1):
                                if i > clen0:
                                    raise Exception("End of interval '"+item+"' is higher than metabolite '"+metab0+"' length "+str(clen0)+" (%s: %d)"%(ftbl["name"], irow))
                                mask|=1<<(clen0-i)
                        except:
                            raise Exception("Badly formed fragment interval '"+item+"' (%s: %d)"%(ftbl["name"], irow))
                    except:
                        raise Exception("Badly formed fragment interval '"+item+"' (%s: %d)"%(ftbl["name"], irow))
        m_id=":".join(["m", metabs, frag, str(irow)])
        weight=int(row["WEIGHT"])
        if sumbit(mask) < weight:
            raise Exception("Weight "+str(weight)+" is higher than fragment length "+frag+" (%s: %d)"%(ftbl["name"], irow))
        if clen < sumbit(mask):
            raise Exception("Fragment "+frag+" is longer than metabolite length "+str(clen)+" (%s: %d)"%(ftbl["name"], irow))
        res.setdefault(m_id, dict())
        res[m_id].setdefault(mask, dict())
        if not row["VALUE"]:
            raise Exception("The field VALUE is empty (%s: %s)"%(ftbl["name"], row["irow"]))
        if not row["DEVIATION"]:
            raise Exception("The field DEVIATION is empty (%s: %s)"%(ftbl["name"], row["irow"]))
        try:
            val=float(eval(row["VALUE"]))
        except:
            val=NaN
        try:
            sdev=float(eval(row["DEVIATION"]))
        except:
            raise Exception("DEVIATION must evaluate to a real positive number (%s: %s)."%(ftbl["name"], row["irow"]))
        if sdev <= 0.:
            raise Exception("DEVIATION must be positive (%s: %s)."%(ftbl["name"], row["irow"]))
        res[m_id][mask][weight]={
                "val":val,
                "dev":sdev,
                "id":":".join(["m", metabs, frag, row["WEIGHT"], str(row["irow"])]),
                "irow":str(row["irow"]),
                "pooled":metabl,
        }

def proc_kinopt(ftbl, netan):
    """Proceed label kinetics options from OPTIONS section: file_labcin, dt, tmax, nsubdiv_dt, funlab
    """
    # default values (used when absent in ftbl)
    de={"file_labcin": "", "dt": 1, "tmax": float("inf"), "nsubdiv_dt": 1, "funlabR": ""}
    if "opt" not in netan:
        netan["opt"]=dict((k, []) for k in de)
    # get ftbl OPTIONS -> d
    d=dict()
    for row in ftbl.get("OPTIONS",[]):
        if row["OPT_NAME"][:7] == "funlab:":
            d[row["OPT_NAME"]]=row["OPT_VALUE"]
        else:
            try:
                d[row["OPT_NAME"]]=eval(row["OPT_VALUE"])
            except:
                d[row["OPT_NAME"]]=row["OPT_VALUE"]
    for k in de:
        if k not in netan["opt"]:
            netan["opt"][k]=[]
        elif type(netan["opt"][k]) != list:
            netan["opt"][k]=[] # if already set to a value, reset to a list
        ili=len(netan["opt"][k])
        netan["opt"][k]+=[{}]
        res=netan["opt"][k]
        res[ili]=d.get(k, de[k])
    # append funval entries
    if "funlab" not in netan["opt"]:
        netan["opt"]["funlab"]=[]
    ili=len(netan["opt"]["funlab"])
    netan["opt"]["funlab"]+=[{}]
    netan["opt"]["funlab"][ili]=dict((k[7:],v) for k,v in d.items() if k[:7] == "funlab:")

def mkfunlabli(d):
    "transform 'd' dict to a string representing a body of an R function calculating labeling dependent on time 't'"
    # d is like {'Gluc_1:{32: '{T=0; if (t >= T) 1 else NA}'}}
    return "%(li)s"%{"li": "list("+", ".join("`%(m)s`=list(%(i2f)s)"%{"m": m, "i2f": ", ".join(f"`{iso}`='{f}'" for iso,f in idi.items())} for m,idi in d.items())+")"}
    
