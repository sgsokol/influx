"""- Parse .ftbl
- Analyse ftbl
Restrictions:
 - metabolite name cannot have ":" (it's a separator in measure id)
"""

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
# 2008-06-12 sokol: ftbl_netan(): added label measures analysis
# 2008-06-23 sokol: ftbl_netan(): added peak measures analysis
# 2008-06-23 sokol: ftbl_netan(): added mass measures analysis
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
# 2010-10-11 sokol: ftbl_parse(): for measures added row number for later identification
# 2010-10-11 sokol: ftbl_netan(): added "id" field in measures
# 2010-11-23 sokol: ftbl_netan(): added "nx2dfc" dictionary
# 2010-12-02 sokol: ftbl_netan(): added growth flux option
# 2011-04-20 sokol: cumo_infl(),...: not reversible fluxes are replaced by flux_out in cumomer balance
# 2011-04-27 sokol: ftbl_neta(): added flux_inout
# 2011-04-27 sokol: cumo_infl(),...: not reversible fluxes are replaced by flux_inout in cumomer balance
import numpy as np;
import re;
import copy;
import os;
import sys;

werr=sys.stderr.write;
me=os.path.realpath(sys.argv[0])
dirx=os.path.dirname(me);
sys.path.append(dirx);

from tools_ssg import *;

float_conv=set((
    "VALUE",
    "DEVIATION",
    "OPT_VALUE",
    "META_SIZE",
    "VALUE_S",
    "VALUE_D-",
    "VALUE_D+",
    "VALUE_DD",
    "VALUE(F/C)",
));
def ftbl_parse(f):
    """ftbl_parse(f) -> dict
    read and parse .ftbl file. The only input parameter f is a stream pointer
    with read permission or a file name.
    This function parses the input and returns a dictionnary
    with items corresponding to sections in .ftbl. One section is added.
    "TRANS" correponds to carbon transitions.""";
    import re;
    ftbl={};    # main dictionary to be returned;
    open_here=False;
    #print("f=", f);
    if isstr(f):
        if f[-5:].lower() != ".ftbl":
            f=f+".ftbl";
        f=open(f, "r");
        open_here=True;
    
    #print f;##
    reblank=re.compile("^[\t ]*$");
    recomm=re.compile("^[\t ]*//.*$");
    reading="data";
    irow=0;
    for l in f:
        irow+=1;
        #print "raw l="+l;##
        # strip out \r
        l=l.replace("\r", "");
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
        # skip comments and emty rows;
        if recomm.match(l) or reblank.match(l): continue;
        
        # skip the comments at the end of the row;
        l=re.sub(r"^(.+)//.*$",r"\1",l).rstrip();
        
        # split in fields;
        flds=l.split("\t");
        
        #print "proceeding:"+l;##
        if reading=="data" and len(flds)==1:
            # new section starts here;
            sec_name=flds[0];
            subsec_name="";
            # prepare storage
            ftbl[sec_name]=[];
            # prepare storage for carbon transitions
            if sec_name=="NETWORK":
               #print "TRANS prepared";##
               ftbl["TRANS"]={};
            try: del stock;
            except NameError: pass;
            stock=ftbl[sec_name];
            skiptab=1;
            data_count=0;
            reading="col_names";
            continue;
        if len(flds)==2 and len(flds[0])==0:
            # read subsection name;
            subsec_name=flds[1];
            # prepare sub-storage
            if not ftbl[sec_name]:
                ftbl[sec_name]={};
            ftbl[sec_name][subsec_name]=[];
            try: del stock;
            except NameError: pass;
            stock=ftbl[sec_name][subsec_name];
            skiptab=2;
            data_count=0;
            reading="col_names";
            continue;
        if reading=="col_names" and len(flds)>2:
            # read column names;
            col_names=l[skiptab:].split("\t");
            reading="data";
            #print "col_names=", col_names;##
            continue;
        if reading=="data":
            data=l[skiptab:].split("\t");
            dic={"irow": irow};
            if sec_name == "NETWORK" and len(data[0])==0:
                # we are at carbon transition line here (e.g. #ABC -> #AB +#C)
                #print "data_count="+str(data_count), \
                #    "\ndata="+str([l for l in enumerate(data)]), \
                #    "\nstock="+str(stock);##
                fl_name=str(stock[data_count-1][col_names[0]]);
                for i in range(len(col_names)):
                    try:
                        dic[col_names[i]]=data[i].strip();
                    except IndexError:
                        pass;
                ftbl["TRANS"][fl_name]=dic;
                continue;
            for i in xrange(len(col_names)):
                # classic data
                try:
                    # decimal point conversion
                    dic[col_names[i]]=data[i].strip();
                    if col_names[i] in float_conv:
                       dic[col_names[i]]=dic[col_names[i]].replace(",", ".");
                except IndexError:
                    pass;
            stock.append(dic);
            data_count+=1;
            #print "sec, subsec=", sec_name, subsec_name, data_count, dic;##
        #print "len(flds)=", len(flds), flds, l, data;##
        #print "keys", ftbl.keys(), (ftbl[sec_name].keys() \
        #        if type(ftbl[sec_name])==type({}) else "");##
    if open_here:
        f.close();
    return ftbl;

def ftbl_netan(ftbl):
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
     - input isotopomers (iso_input)
     - input cumomers (cumo_input)
     - flux inequalities (flux_ineqal)
     - flux equalities (flux_eqal)
     - label measures, H1 (label_meas)
     - peak measures, C13 (peak_meas)
     - mass measures (mass_meas)
     - cumomer ordered lists (vcumo)
     - unknown fluxes ordered lists (vflux)
     - linear problem on fluxes (Afl, bfl)
     - free fluxes ordered lists (vflux_free)
     - fw-rv fluxes ordered lists (vflux_fwrv)
     - row names ordered lists for Afl (vrowAfl)
     - in-out fluxes (flux_in, flux_out)
    """
    # init named sets
    netan={
        "input":set(),
        "output":set(),
        "subs":set(),
        "prods":set(),
        "metabs":set(),
        "reac":set(),
        "notrev":set(),
        "sto_r_m":{},
        "sto_m_r":{},
        "flux_m_r":{},
        "cumo_balance_pattern":{},
        "Clen":{},
        "formula":{},
        "metab_netw":{},
        "cumo_sys":{},
        "carbotrans":{},
        "Cmax":{},
        "flux_dep":{},
        "flux_free":{},
        "flux_constr":{},
        "flux_measured":{},
        "iso_input":{},
        "cumo_input":{},
        "flux_inequal":{"net":[], "xch":[]},
        "flux_equal":{"net":[], "xch":[]},
        "label_meas":{},
        "peak_meas":{},
        "mass_meas":{},
        "vcumo":[],
        "Afl":[],
        "bfl":[],
        "vflux":{"net":[], "xch":[], "net2i":{}, "xch2i":{}},
        "vflux_free":{"net":[], "xch":[], "net2i":{}, "xch2i":{}},
        "vflux_constr":{"net":[], "xch":[], "net2i":{}, "xch2i":{}},
        "vflux_meas":{"net":[], "net2i":{}},
        "vflux_compl":{"net":[], "xch":[], "net2i":{}, "xch2i":{}},
        "vflux_fwrv":{"fw":[], "rv":[], "fw2i":{}, "rv2i":{}},
        "vrowAfl":[],
        "flux_in":set(),
        "flux_out":set(),
        "flux_inout":set(),
        "opt": {},
        "met_pools": {},
        "nx2dfc": {},
    };
    res="";     # auxiliary short-cut to current result;

    # Isotopomer dynamics and other options
    for row in ftbl.get("OPTIONS",[]):
        try:
            netan["opt"][row["OPT_NAME"]]=eval(row["OPT_VALUE"]);
        except:
            netan["opt"][row["OPT_NAME"]]=row["OPT_VALUE"];
    for row in ftbl.get("METABOLITE_POOLS",[]):
        netan["met_pools"][row["META_NAME"]]=eval(row["META_SIZE"]);

    # analyse networks
    netw=ftbl.get("NETWORK");
    if not netw:
        raise Exception("No NETWORK section in ftbl");
    row_to_del=[];
    for row in netw:
        reac=row["FLUX_NAME"];
        #print "reac="+reac;#
        e1=row.get("EDUCT_1");
        if not e1: raise Exception("EDUCT_1 must be defined in flux "+reac);
        e2=row.get("EDUCT_2");
        p1=row.get("PRODUCT_1");
        if not p1: raise Exception("PRODUCT_1 must be defined in flux "+reac);
        p2=row.get("PRODUCT_2");
        
        # local substrate (es), product (ps) and metabolites (ms) sets
        es=set((e1,));
        if e2:
            es.add(e2);
        ps=set((p1,));
        if p2:
            ps.add(p2);
        ms=es|ps;
        
        # all reactions A+B=C or C=A+B or A+B=C+D
        netan["reac"].add(reac);
        netan["subs"].update(es);
        netan["prods"].update(ps);
        netan["metabs"].update(ms);
        #aff("ms for "+reac, ms);##

        # create chemical formula
        netan["formula"][reac]={"left":es, "right":ps, "all":ms};

        # Carbon length and transitions
        try:
            trans=ftbl["TRANS"][reac];
            
            netan["carbotrans"][reac]={"left":[], "right":[]};
            for (m,carb,lr) in [(e1,trans.get("EDUCT_1"),"left"),
                    (e2,trans.get("EDUCT_2"),"left"),
                    (p1,trans.get("PRODUCT_1"),"right"),
                    (p2,trans.get("PRODUCT_2"),"right")]:
                    #print "m="+str(m), "; carb="+str(carb);##
                if not m or not carb:
                    continue;
                if carb[0] != "#":
                    raise Exception("In carbon string for metabolite "+m+" a starting '#' is missing."+
                        "\nreaction="+str(row)+"\ncarbons ="+str(trans));
                # carbon transitions
                netan["carbotrans"][reac][lr].append((m,carb[1:])); # strip "#" character

                # carbon length
                if netan["Clen"].get(m, 0) and \
                        netan["Clen"][m] != len(carb)-1:
                    raise Exception("CarbonLength", "Metabolite "+m+" has length "+
                            str(netan["Clen"][m])+" but in reaction "+reac+
                            " it has legth "+str(len(carb)-1));
                netan["Clen"][m]=len(carb)-1; # don't count '#' character
        except KeyError:
            werr("CarbonTrans: No reaction "+reac+" in carbon transitions\n");
        except Exception as inst:
            werr(": ".join(inst)+"\n");

        # stocheometric matrix in dictionnary form
        # sto_r_m[reac]['left'|'right']=list(metab) and
        # sto_m_r[metab]['left'|'right']=list(reac)
        # substrates are in 'left' list
        # and products are in the 'right' one.

        netan["sto_r_m"][reac]={"left":[e1], "right":[p1]};
        if (e2):
            netan["sto_r_m"][reac]["left"].append(e2);
        if (p2):
            netan["sto_r_m"][reac]["right"].append(p2);
##        #print "col keys", col.keys();##
        
        for s in netan["sto_r_m"][reac]["left"]:
           #print "sto_m_r s="+str(s);##
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]};
            netan["sto_m_r"][s]["left"].append(reac);
        for s in netan["sto_r_m"][reac]["right"]:
           #print "sto_m_r s="+str(s);
            if not s in netan["sto_m_r"]:
                netan["sto_m_r"][s]={"left":[], "right":[]};
            netan["sto_m_r"][s]["right"].append(reac);
        # end netan["sto_m_r"]
        #aff("sto_m_r"+str(reac), netan["sto_m_r"]);##
    # add growth fluxes if requested
    if netan["opt"].get("include_growth_flux"):
        if not netan["opt"].get("mu"):
            raise Exception("Growth parameter mu is not set in OPTIONS section");
        netan["flux_growth"]=set();
        for (m,si) in netan["met_pools"].iteritems():
            if not m in netan["metabs"]:
                # unknown metabolite, skip its proceeding
                continue;
            mgr=m+"_gr";
            reac=mgr;
            if reac in netan["reac"]:
                raise Exception("Cannot add growth reaction "+reac+". It is already in network");
            netan["reac"].add(reac);
            netan["subs"].update((m,));
            netan["prods"].update((mgr,));
            netan["metabs"].update((m,mgr));
            netan["flux_growth"].add(reac);
            #aff("ms for "+reac, ms);##

            # create chemical formula
            netan["formula"][reac]={"left":[m], "right":[mgr], "all":[m,mgr]};
            # stoicheometry
            netan["sto_m_r"][m]["left"].append(reac);
            netan["sto_m_r"][mgr]={"left": [], "right": [reac]};
            netan["sto_r_m"][reac]={"left":[m], "right":[mgr]};
    # find input and output metabolites
    netan["input"].update(netan["subs"]-netan["prods"]);
    netan["output"].update(netan["prods"]-netan["subs"]);
    
    # input and output flux names
    netan["flux_in"]=set(valval(netan["sto_m_r"][m]["left"] for m in netan["input"]));
    netan["flux_out"]=set(valval(netan["sto_m_r"][m]["right"] for m in netan["output"]));
    netan["flux_inout"]=netan["flux_in"] | netan["flux_out"]
    #print "fl in + out="+str(netan["flux_in"])+"+"+str( netan["flux_out"]);##

    # analyse reaction reversibility
    # and free fluxes(dictionary flux->init value for simulation or minimization)
    # and constrained fluxes (dictionary flux->value)
    # The rest are dependent fluxes
    netan["flux_free"]={"net":{}, "xch":{}};
    netan["flux_constr"]={"net":{}, "xch":{}};
    netan["flux_dep"]={"net":set(), "xch":set()};
    flx=ftbl.get("FLUXES");
    if flx:
        xch=flx.get("XCH");
        net=flx.get("NET");
##        aff("net", net);
    else:
        werr(os.path.basename(me)+": netan[FLUXES][XCH] is not defined");
        #return None;
    #print "list reac=", netan["reac"];##
    try:
        for reac in netan["reac"] | netan["flux_inout"]:
            # get xch condition for this reac
            cond=[row for row in xch if row["NAME"]==reac];
            # get net condition for this reac
            ncond=([row for row in net if row["NAME"]==reac]) if net else [];
            # no xch dispatch check for input/output fluxes as they are
            # constrained by definition
            #print "r,c,n=", reac, len(cond), len(ncond);##
            if len(cond) > 1:
                raise Exception("FluxDef", "Reaction `%s` is over defined in exchange fluxes."%reac);
            if len(cond) == 0 and not (reac in (netan["flux_in"] | netan["flux_out"])):
                werr("in+out fluxes="+str(netan["flux_inout"]));##
                raise Exception("FluxDef", "Reaction `%s` is not defined D/F/C in exchange fluxes."%reac);
            if cond:
                cond=cond[0];
            if len(ncond) > 1:
                raise Exception("FluxDef", "Reaction `%s` is over defined in net fluxes."%reac);
            if len(ncond) == 0 and not (reac in (netan["flux_inout"])):
                raise Exception("FluxDef", "Reaction `%s` is not defined D/F/C in net fluxes."%reac);
            if (len(ncond) == 0 and
                (reac in netan["flux_inout"] and
                not reac in netan.get("flux_growth", ()))):
                # input-output fluxes are by definition not reversible
                netan["notrev"].add(reac);
                netan["flux_constr"]["xch"][reac]=0.;
                continue;
            if netan["opt"].get("include_growth_flux") and reac in netan["flux_growth"]:
                m=reac[:-3];
                netan["notrev"].add(reac);
                netan["flux_constr"]["xch"][reac]=0.;
                netan["flux_constr"]["net"][reac]=netan["opt"].get("metab_scale", 1.)*netan["opt"]["mu"]*netan["met_pools"][m];
                continue;
            ncond=ncond[0];
            #aff("cond", cond);##
            #aff("ncond", ncond);##
            # not reversibles are those reaction having xch flux=0 or
            # input/output fluxes
            try:
                if (cond and cond["FCD"] == "C") or (reac in netan["flux_inout"]):
                    if (reac in netan["flux_inout"]) or \
                        (eval(cond["VALUE(F/C)"]) == 0.):
                        netan["notrev"].add(reac);
                        netan["flux_constr"]["xch"][reac]=0.;
                    else:
                        # not null constrained xch
                        val=eval(cond["VALUE(F/C)"])
                        if val < 0. or val > 1.:
                            raise Exception("FluxVal", "Constrained exchange flux`%s` must have value in [0; 1] interval.\nInstead, a value %f is given"%(reac, val));
                        netan["flux_constr"]["xch"][reac]=val;
                elif (cond and cond["FCD"] == "F"):
                    # free xch
                    val=eval(cond["VALUE(F/C)"])
                    if val < 0. or val > 1.:
                        raise Exception("FluxVal", "Free exchange flux`%s` must have starting value in [0; 1] interval.\nInstead, a value %f is given"%(reac, val));
                    netan["flux_free"]["xch"][reac]=val;
                elif (cond and cond["FCD"] == "D"):
                    # dependent xch
                    netan["flux_dep"]["xch"].add(reac);
                if (ncond and ncond["FCD"] == "F"):
                    # free net
                    netan["flux_free"]["net"][reac]=eval(ncond["VALUE(F/C)"]);
                elif (ncond and ncond["FCD"] == "C"):
                    # constr net
                    netan["flux_constr"]["net"][reac]=eval(ncond["VALUE(F/C)"]);
                elif (ncond and ncond["FCD"] == "D"):
                    # dependent net
                    netan["flux_dep"]["net"].add(reac);
    ##        aff("free net", netan["flux_free"]["net"]);
    ##        aff("free xch", netan["flux_free"]["xch"]);
    ##        aff("constr net", netan["flux_free"]["net"]);
    ##        aff("free xch", netan["flux_free"]["xch"]);
            except:
                werr("Error occured on the row %d of ftbl file:\n"%cond["irow"]);
                raise;
    except Exception as inst:
        werr(": ".join(inst)+"\n");
    
    # measured fluxes
    for row in ftbl.get("FLUX_MEASUREMENTS",[]):
        if row["FLUX_NAME"] not in netan["reac"]:
            raise Exception("Mesured flux `%s` is not defined."%row["FLUX_NAME"]);
        if row["FLUX_NAME"] not in netan["flux_free"]["net"] and \
            row["FLUX_NAME"] not in netan["flux_dep"]["net"]:
            raise Exception("Mesured flux `%s` must be defined as either free or dependent."%row["FLUX_NAME"]);
        netan["flux_measured"][row["FLUX_NAME"]]={\
                "val": eval(row["VALUE"]), \
                "dev": eval(row["DEVIATION"])};
    # input isotopomers
    for row in ftbl.get("LABEL_INPUT",[]):
        metab=row.get("META_NAME", "") or metab;
        if metab not in netan["Clen"]:
            raise Exception("Input metabolite `%s` is not defined."%metab);
        ilen=len(row.get("ISOTOPOMER", ""))-1; # -1 for '#' sign
        if ilen != netan["Clen"][metab]:
            raise Exception("Input isotopomer `%s` has bad length. A length of %d is expected."%
                (row.get("ISOTOPOMER", ""), netan["Clen"][metab]));
        iiso=strbit2int(row["ISOTOPOMER"]);
        if metab not in netan["iso_input"]:
            netan["iso_input"][metab]={};
        netan["iso_input"][metab][iiso]=eval(row["VALUE"]);
    if set(netan["iso_input"].keys()) != set(netan["input"]):
        werr("LabelInput: Label input section must contain all network input metabolites and only them\n"+
            "label_input: "+str(netan["iso_input"].keys())+"\n"+
            "network input: "+str(netan["input"])+"\n");
    
    # translate input iso to input cumo
    for metab in netan["iso_input"]:
       for icumo in xrange(1,1<<(netan["Clen"][metab])):
           cumo=metab+":"+str(icumo);
           netan["cumo_input"][cumo]=iso2cumo(netan["Clen"][metab], icumo, netan["iso_input"][metab]);
    
    # flux inequalities
    # list of tuples (value,comp,dict) where dict is flux:coef
    # and comp is of "<", "<=", ...
    # net fluxes
    for row in ftbl.get("INEQUALITIES",{}).get("NET",[]):
        #print row;##
        dicf=formula2dict(row["FORMULA"])
        fl=dicf.keys()[0];
        if len(dicf)==1 and fl in netan["flux_constr"]["net"]:
            werr("Inequalities: in NET section, the formula '"+
                row["VALUE"]+row["COMP"]+row["FORMULA"]+"' implies a constrained flux\n"+
                " having a value "+str(netan["flux_constr"]["net"][fl])+".\n The inequality is ignored as meaningless.\n")
        netan["flux_inequal"]["net"].append((
                eval(row["VALUE"]),
                row["COMP"],
                dicf));
    # xch fluxes
    for row in ftbl.get("INEQUALITIES",{}).get("XCH",[]):
        #print row;##
        netan["flux_inequal"]["xch"].append((
                eval(row["VALUE"]),
                row["COMP"],
                formula2dict(row["FORMULA"])));
    for (afftype, ftype) in (("Net", "net"), ("Exchange", "xch")):
        for row in netan["flux_inequal"][ftype]:
            for fl in row[2]:
                if fl not in netan["reac"]:
                   raise Exception("%s flux `%s` in the inequality\n%s\nis not defined."%
                       (afftype, fl, join("", row)));

    # flux equalities
    # list of tuples (value,dict) where dict is flux:coef
    # net fluxes
    for row in ftbl.get("EQUALITIES",{}).get("NET",[]):
        #print row;##
        netan["flux_equal"]["net"].append((
                eval(row["VALUE"]),
                formula2dict(row["FORMULA"])));
    # xch fluxes
    for row in ftbl.get("EQUALITIES",{}).get("XCH",[]):
        #print row;##
        netan["flux_equal"]["xch"].append((
                eval(row["VALUE"]),
                formula2dict(row["FORMULA"])));
    for (afftype, ftype) in (("Net", "net"), ("Exchange", "xch")):
        for row in netan["flux_equal"][ftype]:
            for fl in row[1]:
                if fl not in netan["reac"]:
                   raise Exception("%s flux `%s` in the equality\n%s\nis not defined."%
                       (afftype, fl, join("=", row)));
    # Check that fluxes are all in reactions
    # and form nx2dfc dictionary
    for (affnx, nx, nxsh) in (("net", "net", "n."), ("exchange", "xch", "x.")):
        for (affdfc, dfc, dfcsh) in (("Dependent", "flux_dep", "d."), ("Free", "flux_free", "f."), ("Constrained", "flux_constr", "c.")):
            #print netan[dfc][nx];##
            for fl in netan[dfc][nx]:
                if fl not in netan["reac"]:
                   raise Exception("%s %s flux `%s` is not defined."%
                       (affdfc, affnx, fl));
                netan["nx2dfc"][nxsh+fl]=dfcsh+nxsh+fl;
    # Check that all reactions are dispathed between dependent free and constrained
    # fluxes

    # label measurements
    # [metab][group]->list of {val:x, dev:y, bcumos:list of #bcumo}}
    metab="";
    for row in ftbl.get("LABEL_MEASUREMENTS",[]):
        #print row;##
        # test the cumomer pattern validity
        if (not re.match(r"#[01x]+(\+#[01x]+)*", row["CUM_CONSTRAINTS"])):
            raise Exception("Not valid cumomer's pattern in '"+row["CUM_CONSTRAINTS"]+"'");
        metab=row["META_NAME"] or metab;
        if not metab in netan["metabs"]:
            raise Exception("Unknown metabolite name '"+metab+"' in LABEL_MEASUREMENTS");
        if metab in netan["output"]:
            raise Exception("""Measured metabolites have to be internal to network.
You can add a fictious metabolite following to '"""+metab+"' (seen in LABEL_MEASUREMENTS).");
        mlen=netan["Clen"][metab];
        group=row["CUM_GROUP"] or group;
        if not metab in netan["label_meas"]:
            netan["label_meas"][metab]={};
        if not group in netan["label_meas"][metab]:
            netan["label_meas"][metab][group]=[];
        # prepare cumos list
        # if bcumos is not empty use it
        # else use group name as carbon number (starting from # character)
        if row.get("CUM_CONSTRAINTS",""):
            bcumos=row["CUM_CONSTRAINTS"].split("+");
        else:
            try:
                # just put "1" in group-th place
                i=int(group);
                bcumos="#"+"x"*netan["Clen"][metab];
                bcumos[i-1]="1";
                bcumos=[bcumos];
            except:
                raise Exception("Expected integer CUM_GROUP in LABEL_MEASUREMENTS on row="+str(row));
        netan["label_meas"][metab][group].append({
                "val":eval(row["VALUE"]),
                "dev":eval(row["DEVIATION"]),
                "bcumos":row["CUM_CONSTRAINTS"].split("+"),
                "id":":".join(["l", metab, row["CUM_CONSTRAINTS"], str(row["irow"])])
        });
        # test the icumomer lengths
        if not all(len(ic)==mlen+1 for ic in 
                netan["label_meas"][metab][group][-1]["bcumos"]):
            raise Exception("Wrong cumomer length for "+metab+" in "+row["CUM_CONSTRAINTS"]);
    
    # peak measurements
    # [metab][c_no][peak_type in (S,D-,D+,(DD|T))]={val:x, dev:y}
    metab="";
    for row in ftbl.get("PEAK_MEASUREMENTS",[]):
        #print row;##
        # test the pattern validity
        if (row.get("VALUE_DD","") and row.get("VALUE_T","")):
            raise Exception("Not valid value combination. Only one of DD and T has to be in row "+str(row));
        metab=row["META_NAME"] or metab;
        if not metab in netan["metabs"]:
            raise Exception("Unknown metabolite name ""+metab+"" in PEAK_MEASUREMENTS");
        if metab in netan["output"]:
            raise Exception("""Measured metabolites have to be internal to network.
You can add a fictious metabolite following to '"""+metab+"' (seen in PEAK_MEASUREMENTS).");
        clen=netan["Clen"][metab];
        netan["peak_meas"].setdefault(metab,{});
        for suff in ("S", "D-", "D+", "DD", "T"):
            # get val and dev for this type of peak
            val=row.get("VALUE_"+suff,"");
            if (not val or val=="-"):
                continue;
            dev=row.get("DEVIATION_"+suff,"") or row.get("DEVIATION_S");
            try:
                val=float(eval(val));
                dev=float(eval(dev));
            except:
                continue;
            # test validity
            if not dev:
                raise Exception("Deviation is not determined for VALUE_"+suff+" in row \n"+str(row));
            c_no=int(row["PEAK_NO"]);
            if c_no > clen:
                raise Exception("Carbon number "+str(c_no)+" is gretaer than carbon length "+str(clen)+" for metabolite '"+metab+"' in row \n"+str(row));
            if suff == "D-" and c_no == 1:
                raise Exception("Peak D- cannot be set for metabolite "+metab+", c_no=1 in row \n"+str(row));
            if suff == "D+" and c_no == clen:
                raise Exception("Peak D+ cannot be set for metabolite "+metab+", c_no="+str(c_no)+" in row \n"+str(row));
            if (suff == "DD" or suff == "T") and (c_no == 1 or c_no == clen or clen < 3):
                raise Exception("Peak DD (or T) cannot be set for metabolite "+metab+", c_no="+str(c_no)+", len="+str(clen)+" in row \n"+str(row));
            netan["peak_meas"][metab].setdefault(row["irow"],{});
            netan["peak_meas"][metab][row["irow"]][suff]={
               "val": val,
               "dev": dev,
               "id": ":".join(["p", metab, suff, str(row["irow"])]),
               "irow": str(row["irow"]),
               "c_no": c_no,
            };
    
    # mass measurements
    # dict:[metab][frag_mask][weight]={val:x, dev:y}
    metab="";
    for row in ftbl.get("MASS_SPECTROMETRY",[]):
        #print row;##
        metab=row["META_NAME"] or metab;
        if row["META_NAME"]:
            irow=str(row["irow"])
        clen=netan["Clen"][metab];
        # test the validity
        if not metab in netan["metabs"]:
            raise Exception("Unknown metabolite name '"+metab+"' in MASS_SPECTROMETRY");
        if metab in netan["output"]:
            raise Exception("""Measured metabolites have to be internal to network.
You can add a fictious metabolite following to '"""+metab+"' (seen in MASS_MEASUREMENTS).");
        frag=row["FRAGMENT"] or frag;
        m_id=":".join(["m", metab, frag, irow])
        if row["FRAGMENT"]:
            # recalculate fragment mask
            mask=0;
            for item in frag.split(","):
                try:
                    i=int(item);
                    # add this simple item to the mask
                    mask|=1<<(clen-i);
                except ValueError:
                    # try the interval
                    try:
                        (start,end)=item.split("~");
                        #print "start,end=%s,%s" % (start,end);##
                        try:
                            for i in xrange(int(start),int(end)+1):
                                if i > clen:
                                    raise Exception("End of interval '"+item+"' is higher than metabolite "+metab+" length "+str(netan["Clen"][metab])+". \nMASS_SPECTROMETRY, row="+str(row));
                                mask|=1<<(clen-i);
                        except:
                            raise Exception("Badly formed fragment interval '"+item+"' in MASS_SPECTROMETRY,\nrow="+str(row));
                    except:
                        raise Exception("Badly formed fragment interval '"+item+"' in MASS_SPECTROMETRY,\nrow="+str(row));
        weight=int(row["WEIGHT"]);
        if sumbit(mask) < weight:
            raise Exception("Weight "+str(weight)+" is higher than fragment length "+frag+" in MASS_SPECTROMETRY\nrow="+str(row));
        if clen < sumbit(mask):
            raise Exception("Fragment "+frag+" is longer than metabolite length "+str(clen)+" in MASS_SPECTROMETRY\nrow="+str(row));
        netan["mass_meas"].setdefault(m_id, {});
        netan["mass_meas"][m_id].setdefault(mask, {});
        netan["mass_meas"][m_id][mask][weight]={
                "val":eval(row["VALUE"]),
                "dev":eval(row["DEVIATION"]),
                "id":":".join(["m", metab, frag, row["WEIGHT"], str(row["irow"])]),
                "irow":str(row["irow"]),
        }
    
    # discard empty entries
    for e in netan:
        try:
            netan[e].discard("");
            netan[e].discard(None);
        except AttributeError: pass;

    # fwd-rev flux balance
    del(res);
    res=netan["flux_m_r"];
    # res[metab]["in"]=list(fwd.flux|rev.flux)
    # res[metab]["out"]=list(fwd.flux|rev.flux)
    for metab,lr in netan["sto_m_r"].iteritems():
        # lr is dico with 'left' and 'right' entries
        #print "fwd metab="+str(metab);##
        if not metab in res:
            res[metab]={"in":[], "out":[]};
        for reac in lr["left"]:
            # here metabolite is consumed in fwd reaction
            # and produced in the reverse one.
            res[metab]["out"].append("fwd."+reac);
            #if reac not in netan["notrev"]:
            if reac not in netan["flux_inout"]:
                res[metab]["in"].append("rev."+reac);
        for reac in lr["right"]:
            # here metabolite is consumed in rev reaction
            # and produced in the forward one.
            res[metab]["in"].append("fwd."+reac);
            #if reac not in netan["notrev"]:
            if reac not in netan["flux_inout"]:
                res[metab]["out"].append("rev."+reac);

    # metabolite network
    del(res);
    res=netan["metab_netw"];
    # left part or reaction points to right part
    for reac,parts in netan["formula"].iteritems():
        for metab in parts["left"]:
            if not metab in res:
                res[metab]=set();
            res[metab].update(parts["right"]);
    # order cumomers by weight. For a given weight, cumomers are sorted by
    # metabolites order.
    netan["Cmax"]=max(netan["Clen"].values());
    Cmax=netan["Cmax"];
    #return netan;##
    # cumomers systems A*x=b, one by weight
    # weights are going from 1 to Cmax where Cmax is the maximal
    # carbon string length in all metabolites
    # A is a list of matrices (in weight order)
    # b is a list of right hand parts (still in weight order)
    # the dimensions of various weights in b are not the same
    # too short metabolites are dropped when going to higher weights.
    del(res);
    res=netan["cumo_sys"];
    res["A"]=[{} for i in xrange(Cmax)];
    res["b"]=[{} for i in xrange(Cmax)];
    try:
        # run through all reactions and update bilan of involved cumomers
        for (reac,lrdict) in netan["carbotrans"].iteritems():
            # run through metabs
            ## aff("lrdict", lrdict);#
            for (imetab,lr,metab,cstr) in ((imetab,lr,metab,cstr)
                    for (lr,lst) in lrdict.iteritems()
                    for (imetab,(metab,cstr)) in enumerate(lst)):
                # if output metab then influx is set to 1
                # so its cumomer distribution is directly
                # defined by cumodistr of inputs
                Clen=netan["Clen"][metab];
                # input metabolite has fixed value so put it in rhs
                # when it is an influx for some internal cumomer
                if metab in netan["input"] or metab in netan["output"]:
                    continue;
                # 'out' part of this metab
                fwd_rev=("fwd." if lr=="left" else "rev.");
                flux=fwd_rev+reac;
                #if (fwd_rev=="fwd." or reac not in netan["notrev"]):
                if (fwd_rev=="fwd." or reac not in netan["flux_inout"]):
                    # add this out-flux;
                    # run through all cumomers of metab
                    for icumo in xrange(1,1<<Clen):
                        cumo=metab+":"+str(icumo);
                        w=sumbit(icumo);
                        #print "w,i,clen,metab=", w, icumo,Clen,metab;##
                        if cumo not in res["A"][w-1]:
                            res["A"][w-1][cumo]={cumo:[]};
                        # main diagonal term ('out' part)
                        res["A"][w-1][cumo][cumo].append(flux);
                        ##print 'm,ic,w='+metab, icumo, w;#
                        ##aff("res["A"][w-1][cumo][cumo]", res["A"][w-1][cumo][cumo]);#
                # 'in' part
                fwd_rev=("rev." if lr=="left" else "fwd.");
                flux=fwd_rev+reac;
                in_lr=("left" if lr=="right" else "right");
                #if (fwd_rev=="rev." and reac in netan["notrev"]):
                if (fwd_rev=="rev." and reac in netan["flux_inout"]):
                    # this cannot be by definition
                    continue;
                # add this in-flux;
                for (in_i,(in_metab, in_cstr)) in enumerate(lrdict[in_lr]):
                    # run through all cumomers of metab
                    for icumo in xrange(1,1<<Clen):
                        cumo=metab+":"+str(icumo);
                        w=sumbit(icumo);
                        # get in_cumo
                        in_icumo=src_ind(in_cstr, cstr, icumo);
                        if in_icumo==None:
                            continue;
                        in_cumo=in_metab+":"+str(in_icumo);
                        in_w=sumbit(in_icumo);
                        if cumo not in res["A"][w-1]:
                            res["A"][w-1][cumo]={cumo:[]};
                        if in_w==w:
                            if in_metab in netan["input"]:
                                # put it in rhs
                                if cumo not in res["b"][w-1]:
                                    res["b"][w-1][cumo]={};
                                if flux not in res["b"][w-1][cumo]:
                                    res["b"][w-1][cumo][flux]={};
                                if imetab not in res["b"][w-1][cumo][flux]:
                                    res["b"][w-1][cumo][flux][imetab]=[];
                                if in_cumo in netan["cumo_input"]:
                                    #fract=netan["cumo_input"][in_cumo];
                                    pass;
                                else:
                                    raise Exception("LabelInput", in_cumo+" should be in input section of ftbl file:\n"+
                                        str(netan["input"])+"\n"+str(netan["cumo_input"]));
                                #res["b"][w-1][cumo][flux][imetab].append(fract);
                                res["b"][w-1][cumo][flux][imetab].append(in_cumo);
                            else:
                                if in_cumo not in res["A"][w-1][cumo]:
                                    res["A"][w-1][cumo][in_cumo]=[];
                                # matrix: linearized off-diagonal term
                                res["A"][w-1][cumo][in_cumo].append(flux);
                        elif in_w < w:
                            # put lighter cumomer product in rhs list[iterm]
                            if cumo not in res["b"][w-1]:
                                res["b"][w-1][cumo]={};
                            if flux not in res["b"][w-1][cumo]:
                                res["b"][w-1][cumo][flux]={};
                            if imetab not in res["b"][w-1][cumo][flux]:
                                res["b"][w-1][cumo][flux][imetab]=[];
                            #res["b"][w-1][cumo][flux][imetab].append(
                            #    in_cumo if in_metab not in netan["input"] else
                            #    netan["cumo_input"][in_cumo]);
                            res["b"][w-1][cumo][flux][imetab].append(in_cumo);
                            #print "b="+str(res["b"][w-1][cumo][flux]);##
                        # if in_w==0 then in_cumo=1 by definition => ignore here
                        # in_w cannot be > w because of src_ind();
    except Exception as inst:
        werr(": ".join(inst)+"\n");
    # ordered cumomer lists
    for w in xrange(1,netan["Cmax"]+1):
        # weight 1 equations have all metabolites
        ##aff("A "+str(w), netan["cumo_sys"]["A"][w-1]);#
        ##aff("b "+str(w), netan["cumo_sys"]["b"][w-1]);#
        # order cumos along pathways
        # starts are input cumomers
        starts=[cumo for cumo in netan["cumo_sys"]["A"][w-1] \
            if cumo.split(":")[0] in netan["input"]];
        ##aff("st "+str(w), starts);
        # complete starts by all others cumomers
        starts+=[c for c in netan["cumo_sys"]["A"][w-1] if not c in starts];
        cumo_paths=cumo_path(starts, netan["cumo_sys"]["A"][w-1], set());
        # order
        netan["vcumo"].append([cumo for cumo in valval(cumo_paths)]);

    # ordered unknown flux lists
    # get all reactions which are not constrained neither free
    netan["vflux"]["net"].extend(reac for reac in netan["reac"]
        if reac not in netan["flux_constr"]["net"] and
        reac not in netan["flux_free"]["net"]);
    netan["vflux"]["xch"].extend(reac for reac in netan["reac"]
        if reac not in netan["flux_constr"]["xch"] and
        reac not in netan["flux_free"]["xch"]);

    # order
    netan["vflux"]["net"].sort();
    netan["vflux"]["xch"].sort();
    
    # easy index finder
    netan["vflux"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux"]["net"]))
    netan["vflux"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux"]["xch"]))


    # ordered free flux lists
    netan["vflux_free"]["net"]=netan["flux_free"]["net"].keys();
    netan["vflux_free"]["xch"]=netan["flux_free"]["xch"].keys();

    # order
    netan["vflux_free"]["net"].sort();
    netan["vflux_free"]["xch"].sort();
    
    # easy index finder
    netan["vflux_free"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_free"]["net"]))
    netan["vflux_free"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_free"]["xch"]))


    # ordered constrained flux lists
    netan["vflux_constr"]["net"]=netan["flux_constr"]["net"].keys();
    netan["vflux_constr"]["xch"]=netan["flux_constr"]["xch"].keys();

    # order
    netan["vflux_constr"]["net"].sort();
    netan["vflux_constr"]["xch"].sort();
    
    # easy index finder
    netan["vflux_constr"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_constr"]["net"]))
    netan["vflux_constr"]["xch2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_constr"]["xch"]))


    # ordered measured flux lists
    netan["vflux_meas"]["net"]=netan["flux_measured"].keys();

    # order
    netan["vflux_meas"]["net"].sort();
    
    # easy index finder
    netan["vflux_meas"]["net2i"]=dict((fl,i) for (i,fl) in enumerate(netan["vflux_meas"]["net"]))


    # ordered complete flux lists
    netan["vflux_compl"]={
        "net": 
        netan["vflux"]["net"]+
        netan["vflux_free"]["net"]+
        netan["vflux_constr"]["net"],
        "xch":
        netan["vflux"]["xch"]+
        netan["vflux_free"]["xch"]+
        netan["vflux_constr"]["xch"],
    };
    # easy index finding
    # net
    netan["vflux_compl"]["net2i"]=dict(
        (fl,i) for (i,fl) in enumerate(netan["vflux"]["net"])
    );
    netan["vflux_compl"]["net2i"].update(dict(
        (fl,i+len(netan["vflux"]["net"]))
        for (i,fl) in enumerate(netan["vflux_free"]["net"])
    ));
    netan["vflux_compl"]["net2i"].update(dict(
        (fl,i+len(netan["vflux"]["net"])+len(netan["vflux_free"]["net"]))
        for (i,fl) in enumerate(netan["vflux_constr"]["net"])
    ));
    # xch
    netan["vflux_compl"]["xch2i"]=dict(
        (fl,i) for (i,fl) in enumerate(netan["vflux"]["xch"])
    );
    netan["vflux_compl"]["xch2i"].update(dict(
        (fl,i+len(netan["vflux"]["xch"]))
        for (i,fl) in enumerate(netan["vflux_free"]["xch"])
    ));
    netan["vflux_compl"]["xch2i"].update(dict(
        (fl,i+len(netan["vflux"]["xch"])+len(netan["vflux_free"]["xch"]))
        for (i,fl) in enumerate(netan["vflux_constr"]["xch"])
    ));


    # ordered fwd-rev flux lists
    # fw and rv parts are identical so they match each other
    netan["vflux_fwrv"]["fwrv"]=\
        ["fwd."+fl for fl in netan["vflux"]["net"]]+\
        ["fwd."+fl for fl in netan["vflux_free"]["net"]]+\
        ["fwd."+fl for fl in netan["vflux_constr"]["net"]]+\
        ["rev."+fl for fl in netan["vflux"]["net"]]+\
        ["rev."+fl for fl in netan["vflux_free"]["net"]]+\
        ["rev."+fl for fl in netan["vflux_constr"]["net"]];
    
    # order
    netan["vflux_fwrv"]["fwrv"].sort();
    
    # easy index finder
    netan["vflux_fwrv"]["fwrv2i"]=dict((fl,i) for (i,fl) in
        enumerate(netan["vflux_fwrv"]["fwrv"]))
    
    # linear problem on fluxes Afl*(fl_net;fl_xch)=bfl
    # matrix Afl is composed of the following parts :
    # - stocheometic equations (only .net fluxes are involved);
    # - flux equalities
    # constrained to non zero value fluxes are replaced by their values in rhs
    # Afl is a list of lists (two dim array). Primary list is a matrix row,
    # columns are values in the secondary list and contains matrix coefficients
    # bfl is a list of linear expressions. Each expression is a dict
    # where keys are variable names like "f.n.flx" and values are
    # numeric coefficients
    # If the key is empty "", it means that the value is just a scalar to add
    
    # Full flux names are of the format [dfc].[nx].<R>
    # where "d", "f", or "c" corresponds to
    # - dependent
    # - free
    # - constrained
    # and <R> correspond to the reaction name
    
    # stocheometric part
    res=netan["Afl"];
    for (metab,lr) in netan["sto_m_r"].iteritems():
        if (metab in netan["input"] or metab in netan["output"]):
            continue;
        # calculate coefs (repeated fluxes)
        # 'left' part consumes metab
        coefs=list2count(lr["left"], -1);
        # 'right' part produces metab
        # NB: it is assumed that metab cannot appear
        # both in left and right hand side of reaction
        coefs.update(list2count(lr["right"]));
        qry=[coefs.get(fl,0) for fl in netan["vflux"]["net"]];
        qry.extend([0]*len(netan["vflux"]["xch"]));
        if qry==[0]*len(qry):
            # degenerated equation, skip it
            #netan["flux_equal"]["net"].append((0., coefs));
            #raise Exception("Stocheometric equation is zero for metab "+metab+"\n"+str(lr)+"\n"+str(coefs));
            continue;
        # check if this line was already entered before
        for (i,row) in enumerate(res):
            if row==qry:
                werr("Warning: when trying to add balance for "+metab+
                    " got the same equation as for "+netan["vrowAfl"][i]+"\n");
                werr("metab:\t"+join("\t", netan["vflux"]["net"]+netan["vflux"]["xch"])+"\n");
                werr(netan["vrowAfl"][i]+":\t"+join("\t", row)+"\n");
                werr(metab+":\t"+join("\t", qry)+"\n");
                break;
        else:
            # identique row is not found, add it
            res.append(qry);
            netan["vrowAfl"].append(metab);
            # prepare right hand side
            netan["bfl"].append({});
            dtmp=netan["bfl"][-1];
            for fl in coefs:
                #print "bfl: fl="+fl;##
                if fl in netan["flux_free"]["net"]:
                    dtmp["f.n."+fl]=dtmp.get("f.n."+fl,0)-coefs[fl];
                if fl in netan["flux_constr"]["net"]:
                    dtmp["c.n."+fl]=dtmp.get("c.n."+fl,0)-coefs[fl];
            #print "dtmp=", dtmp;

    # flux equality part
    res=netan["Afl"];
    for (nx,nxl) in (("net","n"), ("xch","x")):
        for eq in netan["flux_equal"][nx]:
            qry=[];
            qry.extend(eq[1].get(fl,0) for fl in netan["vflux"][nx]);
            if nx=="net":
                # add xch zeros
                qry.extend([0]*len(netan["vflux"]["xch"]));
            else:
                # prepend zeros for net fluxes
                qry[0:0]=[0]*len(netan["vflux"]["net"]);
            # check qry
            if qry==[0]*len(qry):
                # degenerated equality
                raise Exception("Equality is zero in "+nx+" section: "+str(eq)+"\n");
                continue;
            # check if this line was already entered before
            for row in res:
                if row==qry:
                    raise Exception("An equality in section "+nx+" is redundant. eq:"+str(eq)+
                        "\nqry="+str(qry)+"\nrow="+str(row));
            res.append(qry);
            netan["vrowAfl"].append("eq "+nx+": "+str(eq[1]));
            netan["bfl"].append({"":eq[0]})
            dtmp=netan["bfl"][-1];
            # pass free fluxes to rhs
            for fl in netan["flux_free"][nx]:
                if fl in eq[1]:
                    dtmp["f."+nxl+"."+fl]=dtmp.get("f."+nxl+"."+fl,0)-float(eq[1][fl]);
            # pass constrained fluxes to rhs
            for fl in netan["flux_constr"][nx]:
                if fl in eq[1]:
                    dtmp["c."+nxl+"."+fl]=dtmp.get("c."+nxl+"."+fl,0)-float(eq[1][fl]);
    return netan;

def enum_path(starts, netw, outs, visited=set()):
    """Enumerate metabilites along to reaction pathways.
    Algo: start from an input, follow chemical pathways till an output or
    already visited metabolite. Returns a list of metabolite pathways.
    Each pathways is an ordered list."""
    res=[];
    if (len(visited) == len(netw.keys())):
        # no more metabolites 
        return res;
    for start in starts:
        if (start in visited):
           continue;
        visited.add(start);
        if start in outs:
            return [[start]];
        paths=enum_path(netw[start]-visited, netw, outs, visited);
        if (paths):
            res.append([start]+paths[0]);
            if len(paths) > 1:
                res.extend(paths[1:]);
        else:
            res.append([start]);
    return res;
def cumo_path(starts, A, visited=set()):
    """Enumerate cumomers along to reaction pathways.
    Algo: start from an input, follow chemical pathways till no more
    neighbours or till only visited metabolite rest in network.
    Return a list of cumomer pathways.
    Each pathways is an ordered list."""
    #if [s for s in starts if s.find('.') >= 0]:
    #    aff('s', starts);
    res=[];
    if (len(visited) == len(A)):
        # no more cumomers
        return res;
    for start in starts:
        if (start in visited or not start in A):
           continue;
        visited.add(start);
        # next_starts are cumos influenced by start
        next_starts=set(cumo for cumo in A if start in A[cumo]);
##        aff('next to '+start, next_starts)
        next_starts-=visited;
        if next_starts:
            paths=cumo_path(next_starts, A, visited);
##            aff('paths for '+start, paths)
##            aff('visited', visited)
        else:
            path=[];
            res.append([start]);
            continue;
        if (paths):
            res.append([start]+paths[0]);
            if len(paths) > 1:
                res.extend(paths[1:]);
        else:
            res.append([start]);
##    aff('r', res)
    return res;
    
def src_ind(substrate, product, iprod):
    """
    For a given substrate and product carbon strings (e.g. "abc", "ab")
    calculate substrate index corresponding to product index.
    Return None if no source found.
    Return 0 if iprod==0 and intersection of product and substrate strings
    is not empty"""
    movbit=1;
    isubstr=0;
    substrate=substrate[::-1];
    product=product[::-1];
    for nb in xrange(len(product)):
        if (movbit & iprod):
            # algo: if the current bit is set in product search for its origin
            # and set the corresponding bit in substrate
            try:
                isubstr|=1<<substrate.find(product[nb]);
            except ValueError:
                pass;
        movbit<<=1;
    return isubstr if (isubstr or
        (iprod==0 and (set(product) & set(substrate)))) else None;
def labprods(prods, metab, isostr, strs):
    """labprods(prods, metab, isostr, strs)
    Return a set of tuples (vmetab,visostr) which receive at least
    one labeled carbon from (metab, isostr)
    """
    # run through product part of reaction to get
    # metabs containing labeled letter from isostr
    #print "p=", prods, "m=", metab, "i=", isostr, "s=", strs;##
    res=set();
    # get labeled letters
    lets="".join(s[i] for (i,carb) in enumerate(isostr) if carb=="1" for s in strs);
    #print "ls=", lets;##
    # find labeled letters in reaction products and set
    # them to "1", while others to "0"
    for (vm,vs) in prods:
        rs="";
        for r in vs:
            rs=rs+("1" if r in lets else "0");
            #print "r=",r, "rs=",rs;
        if rs.find("1") > -1:
            res.add((vm,rs));
            #print "append=", vm, rs;##
    #print "res=", res;##
    return res;
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
    in a dictionary isos[cmetab]=list(cisotopes). It is assumed that no
    more than two metabolites can exist in both part of reaction"""
    # run through all combinations of srcs isotopes
    # to get all products
    # |isostr|*|iso2|[+|iso1|*|isostr| if metab is in both position]
    # If there is only one src than m2=""
    #print "allprods: metab=", metab, "isostr=", isostr;
    #print "s=", srcs, "p=", prods, "i=", isos;##
    # find metabolite position(s) in the reaction
    mpos=[i for (i,(m,s)) in enumerate(srcs) if m==metab];
    
    # all source isotop couples (im,ic)
    icouples=[];
    if metab==m1:
        icouples=[(isostr,is2) for is2 in isos.get(m2,set(("",)))];
    if metab==m2:
        icouples.extend((is1,isostr) for is1 in isos.get(m1,set(("",))));
    #print "icpls=", list(icouples);##
    #sys.exit(1);##
    # labeled source letters
    #print "s1=", s1, "s2=", s2;
    lets=["".join(l for (i,l) in enumerate(s1) if is1[i]=="1")+
        "".join(l for (i,l) in enumerate(s2) if is2[i]=="1")
        for (is1,is2) in icouples];
    #print "lets=", lets;##
    # form all products
    res=set();
    for inm in mpos:
        (m,s)=srcs[inm];
        mlab="".join(l for (i,l) in enumerate(s) if isostr[i]=="1");
        # co-substrate (if any)
        (cm,cs)=("","") if len(srcs)==1 else srcs[(inm+1)%2];
        # run through isotops of co-substrate
        for coiso in isos.get(cm,set(("",))):
            lab=mlab+"".join(l for (i,l) in enumerate(cs) if coiso[i]=="1");
            # form product isotops
            for (pm,ps) in prods:
                ip="".join(("1" if l in lab else "0") for l in ps);
                if set(ps)&set(mlab):
                    # get only products labeled by metab+isostr
                    res.add((cm,coiso, pm, ip));
    #print "res=", res;
    return res;
def prod(metab, iso, s, cmetab, ciso, cs, prods):
    """prod(metab, iso, s, cmetab, ciso, cs, prods)->set()
    get isotops from labeled substrates"""
    #print "prod: m=", metab, "s=", s, "i=", iso, "cm=", cmetab, "cs=", cs, "ci=", ciso;
    res=set();
    if cs and ciso == -1:
        return res;
    lab="";
    n=len(s)-1;
    for (i,b) in iternumbit(iso):
        lab=lab+(s[n-i] if b else "");
    n=len(cs)-1;
    for (i,b) in iternumbit(ciso):
        lab=lab+(cs[n-i] if b else "");
    for (pm,ps) in prods:
        ip=sum((1<<i if l in lab else 0) for (i,l) in enumerate(ps[::-1]));
        res.add((pm, ip));
    #print "res=", res;
    return res;
def frag_prod(metab, frag, s, cmetab, cfrag, cs, prods):
    """frag_prod(metab, frag, s, cmetab, cfrag, cs, prods)
    Get fragments from labeled substrates"""
    #print "frag_prod: m=", metab, "s=", s, "f=", frag, "cm=", cmetab, "cs=", cs, "cf=", cfrag;
    res=set();
    if cs and not cfrag:
        # no fragment for this metabolite is yet produced
        return res;
    lab="";
    lab=("".join(s[i] for (i,c) in enumerate(frag) if c!="0")+
        "".join(cs[i] for (i,c) in enumerate(cfrag) if c!="0"));
    flab=frag.replace("0", "")+cfrag.replace("0", "");
    for (pm,ps) in prods:
        fp="".join((flab[lab.index(l)] if l in lab else "0") for (i,l) in enumerate(ps));
        res.add((pm, fp));
    #print "res=", res;
    return res;
def cumo_iw(w,nlen):
    """iterator for a given cumomer weight w in the carbon length nlen"""
    #print ("cumo_iw: w=%d,nlen=%d\n" % (w,nlen));##
    if w < 0 or w > nlen:
       raise Exception("CumoWeightOutOfRange");
    if nlen < 0:
       raise Exception("CumoLengthNegative");
    if w == 1:
        movbit=1;
        for i in xrange(nlen):
            yield movbit;
            movbit<<=1;
    else:
        movbit=1<<(w-1);
        for i in xrange(nlen-w+1):
            for subi in cumo_iw(w-1,w+i-1):
                yield movbit+subi;
            movbit<<=1;
def iso2cumo(Clen, icumo, iso_dic):
    """calculate cumomer fraction from isotopomer ones"""
    return sum(iso_dic.get(iiso,0.)
        for iiso in icumo2iiso(icumo, Clen));
def formula2dict(f):
    """parse a linear combination sum([+|-][a_i][*]f_i) where a_i is a 
    positive number and f_i is a string starting by non-digit and not white
    character (# is allowed). Output is a dict f_i:[+-]a_i"""
    pterm=re.compile(r'\W*([+-])\W*'); # returns list of match1,sep,match2,...
    pflux=re.compile(r'\W*(?P<coef>\d+\.?\d*|^)?\W*\*?\W*(?P<var>[a-zA-Z_][\w\.-]*)\W*');
    res={};
    sign=1;
    l=(i for i in pterm.split(str(f)));
    for term in l:
        try:
            next_sign=l.next();
            next_sign=-1 if next_sign == "-" else 1;
        except StopIteration:
            next_sign=0;
            pass;
        if (len(term) == 0):
            continue;
        m=pflux.match(term);
        if (m):
            coef=m.group("coef");
            coef=1. if coef==None or not len(coef) else float(coef);
            var=m.group("var");
            sign=-1 if sign==-1 else 1;
            res[var]=sign*coef;
            sign=next_sign;
        else:
            raise Exception("Not parsed term '"+term+"'");
    return res;
def label_meas2matrix_vec_dev(netan):
    """use netan["label_meas"] to construct a corresponding measure matrix matx_lab
    such that scale_diag*matx_lab*cumos_vector corresponds to label_measures_vector.
    matx_lab is defined as
    list of dict{"scale":scale_name, "coefs":dict{icumo:coef}, "metab": metabolite}
    where coef is a contribution of cumo in linear combination for given measure.
    scale_name is of the form "metab;group". Group number is to group
    measures of the same metabolite.
    vec is a list of measures (values in .ftbl)
    dev is a list of deviations.
    Elements in matx_lab, vec and dev are ordered in the same way.
    The returned result is a dict (mat,vec,dev)
    """
    # lab[metab][group]=list of {id:txt, val:x, dev:y, bcumos:list of #bcumo}
    # bcumo is like #1xx01 (0,1 and x are allowed)
    # their sum have to be transformed in cumomer linear combination
    mat=[];
    vec=[];
    dev=[];
    for (metab,groups) in netan.get("label_meas", {}).iteritems():
        for (group,rows) in groups.iteritems():
            for row in rows:
                vec.append(row["val"]);
                dev.append(row["dev"]);
                mat.append({"id": row["id"], "scale": metab+";"+group, "coefs":{},
                    "bcumos": row["bcumos"], "metab": metab});
                res=mat[-1]["coefs"];
                for cumostr in row["bcumos"]:
                    decomp=bcumo_decomp(cumostr);
                    for icumo in decomp["+"]:
                        res.setdefault(icumo,0);
                        res[icumo]+=1;
                    for icumo in decomp["-"]:
                        res.setdefault(icumo,0);
                        res[icumo]-=1;
    return {"mat": mat, "vec": vec, "dev": dev};
def mass_meas2matrix_vec_dev(netan):
    """use netan["mass_meas"] to construct a corresponding measure matrix matx_mass
    such that scale_diag*matx_mass*cumos_vector corresponds to mass_measures_vector.
    matx_mass is defined as matx_lab in label_meas2matrix_vec_dev()
    Elements in matx_mass, vec and dev are ordered in the same way.
    scale name is defined as "metab;fragment_mask"
    The returned result is a dict (mat,vec,dev)
    """
    # mass[metab][frag_mask][weight]={val:x, dev:y}
    # weight has to be transformed in cumomer linear combination
    mat=[];
    vec=[];
    dev=[];
    for (m_id,frag_masks) in netan.get("mass_meas",{}).iteritems():
        #print "mass matx calc for ", metab;##
        (l, metab, frag, m_irow)=m_id.split(":");
        n=netan["Clen"][metab];
        str0="0"*n;
        strx="x"*n;
        str1="1"*n;
        for (fmask,weights) in frag_masks.iteritems():
            onepos=[b_no for (b_no,b) in iternumbit(fmask) if b];
            fmask0x=setcharbit(strx,"0",fmask);
            fmask01=setcharbit(str0,"1",fmask);
            nmask=sumbit(fmask);
#            aff("weights for met,fmask "+", ".join((metab,strbit(fmask))), weights);##
            for (weight,row) in weights.iteritems():
                vec.append(row["val"]);
                dev.append(row["dev"]);
                mat.append({"id": row["id"], "scale": metab+";"+fmask01+";"+m_irow, "coefs":{},
                    "bcumos": None, "metab": metab});
                res=mat[-1]["coefs"];
                # for a given weight construct bcumo sum: #x10x+#x01x+...
                bcumos=["#"+setcharbit(fmask0x,"1",expandbit(iw,onepos))
                        for iw in xrange(1<<nmask) if sumbit(iw)==weight];
#                aff("bcumos for met,fmask,w "+", ".join((metab,strbit(fmask),str(weight))), [b for b in bcumos]);##
                mat[-1]["bcumos"]=bcumos;
                for cumostr in bcumos:
                    #print cumostr;##
                    decomp=bcumo_decomp(cumostr);
                    for icumo in decomp["+"]:
                        res.setdefault(icumo,0);
                        res[icumo]+=1;
                    for icumo in decomp["-"]:
                        res.setdefault(icumo,0);
                        res[icumo]-=1;
    return {"mat": mat, "vec": vec, "dev": dev};

def peak_meas2matrix_vec_dev(netan, dmask={"S": 2, "D-": 6, "D+": 3, "T": 7, "DD": 7}):
    """use netan["peak_meas"] to construct a corresponding measure
    matrix matx_peak
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
    mat=[];
    vec=[];
    dev=[];
    for (metab,irows) in netan.get("peak_meas",{}).iteritems():
        #print "peak matx calc for ", metab;##
        n=netan["Clen"][metab];
        strx="#"+"x"*n;
        for (irow,peaks) in irows.iteritems():
            # pmask0x put 0 on three targeted carbons and leave x elsewhere
            c_no=peaks.values()[0]["c_no"]
            pmask0x=setcharbit(strx,"0",1<<(n-c_no));
            if c_no > 1:
                # add left neighbour
                pmask0x=setcharbit(pmask0x,"0",1<<(n-c_no+1));
            if c_no < n:
                # add right neighbour
                pmask0x=setcharbit(pmask0x,"0",1<<(n-c_no-1));
#            aff("peaks for met,c_no,pmask0x="+", ".join((metab,str(c_no),pmask0x)), peaks);##
            for (peak,row) in peaks.iteritems():
                vec.append(row["val"]);
                dev.append(row["dev"]);
                mat.append({"id": row["id"], "scale": metab+";"+str(c_no)+";"+row["irow"], "coefs":{},
                    "bcumos": None, "metab": metab});
                res=mat[-1]["coefs"];
                # for a given (c_no,peak) construct bcumo sum: #x10x+#x01x+...
                # shift the 3-bit mask to the right carbon position
                if c_no < n:
                    pmask=dmask[peak]<<(n-c_no-1);
                else:
                    pmask=dmask[peak]>>1;
                bcumos=[setcharbit(pmask0x,"1",pmask)];
                if peak=="D-" and "T" in peaks:
                    # D+===D- => add D+ to the list of measured bcumomers
                    # set the 3-bit mask to the right carbon position
                    if c_no < n:
                        pmask=dmask["D+"]<<(n-c_no-1);
                    else:
                        pmask=dmask["D+"]>>1;
                    bcumos.append(setcharbit(pmask0x,"1",pmask));
#                aff("bcumos for met,fmask,w "+", ".join((metab,strbit(fmask),str(weight))), [b for b in bcumos]);##
                mat[-1]["bcumos"]=bcumos;
                for cumostr in bcumos:
                    #print cumostr;##
                    decomp=bcumo_decomp(cumostr);
                    for icumo in decomp["+"]:
                        res.setdefault(icumo,0);
                        res[icumo]+=1;
                    for icumo in decomp["-"]:
                        res.setdefault(icumo,0);
                        res[icumo]-=1;
                #print "metab,c_no,peak="+",".join((metab,str(c_no),str(peak)));##
                #print "res=",str(res);##
    return {"mat": mat, "vec": vec, "dev": dev};

def bcumo_decomp(bcumo):
    """bcumo is a string of the form #[01x]+. It has to be decomposed
    in the linear combination of cumomers #[1x]+. The coefficients
    of this linear combination are 1 or -1. So it can be represented as
    sum(cumos_positive)-sum(cumos_negative). The result of this function
    is a dictionary {"+": list of icumos, "-": list of icumos}. icumo is
    an integer who's binary form indicates 1's positions in a cumomer."""
    
    # take zero positions
    zpos=[ i for i in xrange(len(bcumo)-1) if bcumo[-i-1]=="0" ];
    
    # basic 1's cumomer mask
    icumo_base=0;
    for i in xrange(len(bcumo)-1):
        if bcumo[-i-1]=="1":
            icumo_base|=1<<i;
    
    # run through all 2^n zero masks separating positive and negative
    # cumomers
    res={"+": [], "-": []};
    nz=len(zpos);
    for zmask in xrange(1<<nz):
        # odd or even bit number?
        sign="-" if sumbit(zmask)%2 else "+";
        
        # get full cumomer mask
        icumo=icumo_base|expandbit(zmask,zpos);
        # put in result dict
        res[sign].append(icumo);
    #print bcumo+"=" + "+".join(setcharbit("x"*len(bcumo),"1",i) for i in res["+"])+("-" if res["-"] else "") +"-".join(setcharbit("x"*len(bcumo),"1",i) for i in res["-"]);##
    return res;
def mat2graph(A, fp):
    """write digraph file on file pointer fp representing
    links in matrix A given as bi-level dictionnary. A key of
    first level (row index) is influenced by keys of second level
    (column indicies)."""
    fp.write("digraph A {\n");
    for i in A:
        labi=i.replace(":", "_");
        for j in A[i]:
            if i==j:
                continue;
            labj=j.replace(":", "_");
            fp.write("   \""+labj+"\" -> \""+labi+"\"\n");
    fp.write("}\n");
def dom_cmp(A, i, j):
    """Compares influances of i-th and j-th lements of A.
    Returns 0 if i and j are mutually influenced, 1 if i in A[j]
    (i influences j) , -1 if otherwise"""
    return (0 if i in A[j] and j in A[i] else 1 if i in A[j] else -1 if j in A[i] else 0);
def mat2pbm(A, v, fp):
    """Write an image map of non-zero entries of matrix A to file pointer fp.
    Matrix A is a dictionnary, v is a list ordering keys of A."""
    fp.write("P1\n%d %d\n"% (len(A), len(A)));
    for i in v:
        p=0;
        for j in v:
            fp.write("1 " if j in A[i] else "0 ");
            p+=2;
            if p >= 69:
                fp.write("\n");
                p=0;
        fp.write("\n");
        p=0;
def aglom_loop1(A):
    """Agglomerate nodes of A if they are mutually influence each other
    i.e.they are in a loop of length 1.
    Return a new dictionary of influence where entries are those of A aglomerated
    and glued "by" tab symbol"""
    # i->loop_name(=min of all nodes in this loop)
    na=copy.deepcopy(A);
    loop=dict();
    while aglom(na,transpose(na),loop):
        pass;
    return({"na": na, "loop":loop});
def aglom(na,ta,loop):
    """new matrix A (na), transpose A (ta) are used to
    aglomerate neigbour mutually influencing nodes in a supernode.
    Aglomerated noeds are put in the loop dictionnary.
    Return False if no nodes were aglomerated.
    """
    found=False;
    for i in na:
        for j in na[i]:
            if j==i:
                continue;
            if i in na[j]:
                # agglomerate i and j in na
                # which means keep lowest netween i and j
                # lo=min(i,j)
                # hi=max(i,j)
                # and imports influences of hi to lo
                found=True;
                lo=min(i,j);
                hi=max(i,j);
                del(na[hi][hi]);
                del(ta[hi][hi]);
                print "aglom: "+i+"+"+j+"->"+lo;##
                print "elim row "+hi+str(na[hi].keys());##
                print "influenced rows "+str(ta[hi].keys());##
                # update na rows influenced by hi
                for ii in ta[hi].keys():
                    print "\nold row "+ii+str(na[ii].keys());##
                    na[ii].update(na[hi]);
                    del(na[ii][hi]);
                    for jj in na[ii]:
                        ta[jj][ii]=na[ii][jj];
                    print "new row "+ii+str(na[ii].keys());##
                # remove na[hi], ta[hi] and ta's corresponding to na[hi]
                for ii in na[hi]:
                    del(ta[ii][hi]);
                del(na[hi]);
                del(ta[hi]);
                # update loop dictionary
                loop[lo]=loop.get(lo,set((lo,)));
                loop[lo].update(loop.get(hi, set((hi,))));
                if hi in loop:
                    del(loop[hi]);
                print "loop:"+str(loop);
                return(found);
def lowtri(A):
    """Try low triangular ordering of matrix A entries"""
    unsrt=A.keys();
    srt=list();
    # first get inputs (no influences on them)
    srt=[k for k in unsrt if len(A[k])==1];
    for k in srt:
        unsrt.remove(k);
    # now move to lower number keys that influence others
    while unsrt:
        for k in list(unsrt):
            srt.extend(set(A[k].keys()).difference(set(srt)));
            for i in A[k].keys():
                try:
                    unsrt.remove(i);
                except:
                    pass;
    return srt;
def topo_order(A, tA):
    """Try to sort keys of A in topological order. tA is just a transpose of A"""
    unsrt=set(A.keys());
    srtin=list();
    srtout=list();
    
    while unsrt:
        # shave inputs and outputs
        inps=set(k for k in unsrt if len(set(A[k]).difference(srtin))==1);
        srtin.extend(inps);
        unsrt.difference_update(inps);
        outs=set(k for k in unsrt if len(set(tA[k]).difference(unsrt))==1);
        if outs:
            srtout.insert(0, outs);
            unsrt.difference_update(outs);
        if not inps and not outs:
            # this is not a dag
            srtin.extend(unsrt);
            print "not a dag";
            break;
    return(srtin+srtout);
def transpose(A):
    """Transpose a matrix defined as a dict."""
    tA=dict();
    for i in A:
        for j in A[i]:
            tA[j]=tA.get(j,{});
            tA[j][i]=A[i][j];
    return(tA);
def rcumo_sys(netan):
    """Calculate reduced cumomers systems A*x=b
    we start with observed cumomers of max weight
    and we include only needed involved cumomers
    A list of cumomer lists (by weight) is stored
    in netant["vrcumo"]
    """
    # get cumomers involved in measures
    meas_cumos=set();
    # generate measures dico if not yet done
    if "measures" not in netan:
        measures=dict();
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("%s_meas2matrix_vec_dev(netan)"%meas);
        netan["measures"]=measures;
    measures=netan["measures"];
    
    for meas in measures:
        for row in measures[meas]["mat"]:
            metab=row["metab"];
            meas_cumos.update(metab+":"+str(icumo) for icumo in row["coefs"].keys() if icumo != 0);
    # make list of observed weights
    weights=set(sumbit(int(cumo.split(":")[1])) for cumo in meas_cumos);
    if not weights:
        return {"A": {}, "b": {}}
    
    maxw=max(weights);
    weights=range(1,maxw+1);
    weights.reverse();

    # cumomers to visit are stored in lists by weights
    to_visit=dict((w,[]) for w in weights);
    used_fluxes=set();
    A=[{} for i in xrange(maxw)]; # store matrices by weight
    b=[{} for i in xrange(maxw)]; # store rhs by weight
    used=set();

    # initialize to_visit, we'll stop when it's empty
    for cumo in meas_cumos:
        (m,w)=cumo.split(":");
        w=sumbit(int(w));
        to_visit[w].append(cumo);

    # run through the network starting with heaviest cumomers
    for w in weights:
        cumo_i=0;
        while cumo_i < len(to_visit[w]):
            for cumo in to_visit[w]:
                cumo_i+=1;
                if cumo in used:
                    continue;
                # add this cumo used
                used.add(cumo);
                (metab, icumo)=cumo.split(":");
                if metab in netan["output"]:
                    # no equation for output metabs
                    continue;
                icumo=int(icumo);
                A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]});
                # get the influents to cumo of all weights: equal and lower
                # equals go to A and lowers go to b
                infl=cumo_infl(netan, cumo);
                for (incumo,fl,imetab,iinmetab) in infl:
                    if incumo==cumo:
                        # no equation as no transformation
                        continue;
                    (inmetab, inw)=incumo.split(":");
                    inicumo=int(inw);
                    inw=sumbit(inicumo);
                    # input metabolites are to rhs, others are to visit
                    if inmetab not in netan["input"] and inw != 0 and incumo not in to_visit[inw]:
                        to_visit[inw].append(incumo);
                    # main part: write equations
                    if inw==w :
                        # equal weight => A
                        if inmetab in netan["input"]:
                            b[w-1][cumo]=b[w-1].get(cumo, {});
                            b[w-1][cumo][fl]=b[w-1][cumo].get(fl,{});
                            b[w-1][cumo][fl][imetab]=b[w-1][cumo][fl].get(imetab,[]);
                            #b[w-1][cumo][fl][imetab].append(netan["cumo_input"][incumo]);
                            b[w-1][cumo][fl][imetab].append(incumo);
                            continue;
                        A[w-1][cumo][incumo]=A[w-1][cumo].get(incumo,[]);
                        A[w-1][cumo][incumo].append(fl);
                        #aff("A "+str(w)+cumo, A[w-1][cumo]);##
                    elif inw < w:
                        # lower weight => b
                        b[w-1][cumo]=b[w-1].get(cumo, {});
                        b[w-1][cumo][fl]=b[w-1][cumo].get(fl,{});
                        b[w-1][cumo][fl][imetab]=b[w-1][cumo][fl].get(imetab,[]);
                        b[w-1][cumo][fl][imetab].append(incumo);
                # gather all influx in diagonal term
                #aff("A "+str(w)+" "+cumo, A[w-1][cumo]);##
                #aff("b "+str(w)+" "+cumo, b[w-1].get(cumo));##
                for incumo in A[w-1][cumo]:
                    if incumo == cumo:
                        continue;
                    A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]});
                    #print("adding A:", incumo, "->", cumo, A[w-1][cumo][incumo]);##
                    A[w-1][cumo][cumo]+=A[w-1][cumo][incumo];
                if cumo in b[w-1]:
                    A[w-1][cumo]=A[w-1].get(cumo,{cumo:[]});
                    #print("adding b:", cumo, b[w-1][cumo].keys());##
                    A[w-1][cumo][cumo]+=b[w-1][cumo].keys();
    #aff("to_v", to_visit);##
    # make ordered list for reduced cumomer set
    netan["vrcumo"]=copy.deepcopy(netan["vcumo"]);
    for i in xrange(len(netan["vrcumo"]),len(A),-1):
        del(netan["vrcumo"][i-1]);
    for (iw,cumol) in enumerate(netan["vrcumo"]):
        for i in xrange(len(cumol), 0, -1):
            i-=1;
            if cumol[i] not in A[iw]:
                #print "prune", i, cumol[i];##
                del(cumol[i]);
    return {"A": A, "b": b};

def cumo_infl(netan, cumo):
    """cumo_infl(netan, cumo)->list(tuple(in_cumo, fl, imetab, iin_metab))
    return the list of tuples (in_cumo, fl, imetab, iin_metab):
    input cumomer, flux (fwd.fl or rev.fl), index of metab and index of in_metab
    generating cumo. cumo is in format "metab:icumo".
    Condenstation reaction will give the same flux and icumo but various
    iin_metab.
    Convergent point will give multiple fluxes.
    """
    (metab, icumo)=cumo.split(":");
    icumo=int(icumo);
    res=[];
    # run through input forward fluxes of this metab
    for reac in set(netan["sto_m_r"][metab]["right"]):
        # get all cstr for given metab
        for (imetab,cstr) in ((i,s) for (i,(m,s)) in enumerate(netan["carbotrans"][reac]["right"])
            if m==metab):
            # get all input cumomer in this reaction for this cstr
            for (iin_metab, (in_metab,in_str)) in enumerate(netan["carbotrans"][reac]["left"]):
                in_icumo=src_ind(in_str, cstr, icumo);
                if in_icumo != None:
                    in_cumo=in_metab+":"+str(in_icumo);
                    res.append((in_cumo, "fwd."+reac, imetab, iin_metab));
    # run through input reverse fluxes of this metab
    #for reac in set(netan["sto_m_r"][metab]["left"]).difference(netan["notrev"]):
    for reac in set(netan["sto_m_r"][metab]["left"]).difference(netan["flux_inout"]):
        # get all cstr for given metab
        for (imetab,cstr) in ((i,s) for (i,(m,s)) in enumerate(netan["carbotrans"][reac]["left"])
            if m==metab):
            # get all input cumomer in this reaction for this cstr
            for (iin_metab, (in_metab,in_str)) in enumerate(netan["carbotrans"][reac]["right"]):
                in_icumo=src_ind(in_str, cstr, icumo);
                if in_icumo != None:
                    in_cumo=in_metab+":"+str(in_icumo);
                    res.append((in_cumo, "rev."+reac, imetab, iin_metab));
    return res;
def infl(metab, netan):
    """infl(metab, netan)->set(fluxes)
    List incoming fluxes for this metabolite (fwd.reac, rev.reac, ...)
    """
    # run through input forward fluxes of this metab
    res=set("fwd."+reac for reac in netan["sto_m_r"][metab]["right"]);
    # run through input reverse fluxes of this metab
    res.update("rev."+reac for reac in
        set(netan["sto_m_r"][metab]["left"]).difference(netan["flux_inout"]));
#        set(netan["sto_m_r"][metab]["left"]).difference(netan["notrev"]));
    return(res);

def t_iso2m(n):
    """t_iso2m(n) return transition matrix from isotopomers fractions to MID vector
    n - carbon number
    return numpy array of size (n+1,2**n)
    """
    # isotopomer number
    ni=2**n;
    return np.array(
        [[1. if sumbit(ii)==im else 0 for ii in xrange(ni)] for im in xrange(n+1)]
    );

def t_iso2cumo(n):
    """t_iso2cumo(n) return transition matrix from isotopomers fractions to cumomer vector
    n - carbon number
    return numpy array of size (2**n,2**n)
    """
    if n <= 0:
        return np.array(1, ndmin=2);
    # recursive call
    m_1=t_iso2cumo(n-1);
    nc1,nc1=m_1.shape;
    nc=2*nc1;
    m=np.zeros((nc,nc));
    m[0:nc1,0:nc1],m[0:nc1,nc1:],m[nc1:,nc1:]=m_1,m_1,m_1;
    return m;

def t_iso2pos(n):
    """t_iso2pos(n) return transition matrix from isotopomers fractions to positional
    labelling vector (cumomers of weight 1)
    n - carbon number
    return numpy array of size (n,2**n)
    """
    m=t_iso2cumo(n);
    return m[[1<<i for i in xrange(n)],:];

def conv_mid(x,y):
    """conv_mid(x,y)->z
    convolute two mid vectors (numpy arrays)
    and return the result as numpy array.
    """
    nx=len(x);
    ny=len(y);
    z=np.zeros(nx+ny-1);
    for i in xrange(ny):
        z[i:i+nx]+=x*y[i];
    return(z);
