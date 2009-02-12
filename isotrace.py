#!/usr/bin/python
# 2008-09-01 sokol
# follow the path of labeled carbons starting
# from a given isotopomer
# usage: isotrace.py METAB#0010 < file.ftbl
# metabolite METAB must be present in .ftbl file given on stdinput
# and its isotopomer #0010 must match carbon number defined in .ftbl
# Output of the script is a tab separateb list of isotopomer pathes taken by labeled
# carbon(s) till network output, one row by path

# Copyright 2008, INRA/INSA UMR792, MetaSys

# 2008-09-02 initial version, 1C tracer
# 2008-09-03 all products tracer
# 2008-09-05 isotop by int (not str)
# 2008-09-05 minimal step by isotopomer
# 2008-09-08 fragment tracking
# 2009-02-12 non labeled co-sources are used in last resort

import pdb;
import sys;

def usage(mes=""):
    mes and sys.stderr.write(mes+"\n");
    sys.stderr.write("usage: "+sys.argv[0]+" METAB#0010 < file.ftbl\n");

def trace(netan, metab, isostr, visited=set()):
    """get reactions consuming metab
    add their labeled products to result
    and make recursive call on them
    returns empty list if no possibility to pursue exists"""
    
    visited.add((metab,isostr));
    # get labeled neigbours
    v=set();
    for (fr,src,prd) in (("fwd","left","right"),("rev","right","left")):
        v.update([(m,istr,reac+"."+fr) for reac in netan["sto_m_r"][metab][src]
            if src=="left" or reac not in netan["notrev"]
            for (m,istr) in C13_ftbl.labprods(netan["carbotrans"][reac][prd], metab, isostr,
                [s for (sm,s) in netan["carbotrans"][reac][src] if sm==metab])
            if (m,istr) not in visited or m in netan["output"]]);
    #print "m=", metab, "isostr=", isostr, "v=", v;
    #print "vis=", visited;##
    res=[];
    for (met, istr, r) in v:
        if (met,istr) in visited:
            continue;
        if met not in netan["output"]:
            subres=trace(netan, met, istr, visited);
            # add current metabolit as a head of each path
            [sub.insert(0, (metab,isostr,r)) for sub in subres];
            res.extend(subres);
        else:
            visited.add((met,istr));
            res.append([(metab,isostr,r), (met,istr,"end")]);
        #print "res=", res;
    return res;


def front(netan, paths, visited, isos):
    """Track frontal propagation of labeled atoms
    Get last isotop from every path and add its neighbours if
    there is any eligible"""
    eligeable=True;
    #print "paths=", paths;
    #maxpath=500000;##
    step=0;
    while eligeable:
        eligeable=False;
        res=[];
        #if len(paths) > maxpath:
        #    break;
        step+=1;
        cand=[];
        for (ip,p) in enumerate(paths):
            # get last item of this path
            #print "ip=", ip, "p=", p;##
            (reacfr, metab, iso)=p[-1];
            lset=isos.get(metab,[0]*2**netan["Clen"][metab]);
            if not lset[iso]:
                lset[iso]=step;
            isos[metab]=lset;
            #print "front: m=", metab;##
            if metab in netan["output"]:
                continue;
            cand.append((ip, reacfr, metab, iso));
        for (ip, reacfr, metab, iso) in cand:
            # run through all concerned reactionq and
            # all co-isotops to form candidate couples
            # (m,iso,s, cm,ciso,cs)
            v=set();
            for (fr,src,prd) in (("fwd","left","right"),("rev","right","left")):
                for reac in netan["sto_m_r"][metab][src]:
                    if (src=="right") and (reac in netan["notrev"]):
                        # skip reverse of not reversible reaction
                        continue;
                    srcs=netan["carbotrans"][reac][src];
                    prods=netan["carbotrans"][reac][prd];
                    reacfr=reac+"."+fr;
                    #print "reacfr=", reacfr;##
                    mpos=[i for (i,(m,s)) in enumerate(srcs) if m==metab];
                    for inm in mpos:
                        (m,s)=srcs[inm];
                        (cm,cs)=("","") if len(srcs)==1 else srcs[(inm+1)%2];
                        #for ciso in isos.get(cm,set(("0"*len(cs) or "",))):
                        cisos=[i for (i,stp) in enumerate(isos.get(cm,[])) if stp] or [-1];
                        for ciso in cisos:
                            # skip if already visited
                            if (reacfr, metab, iso, cm, ciso) in visited:
                                continue;
                            eligeable=True;
                            # add coisotop to catalog
                            #if cm and ciso > -1:
                            #    clset=isos.get(cm,);
                            #    clset.add(ciso);
                            #    isos[cm]=clset;
                            # register this visit
                            #visited.append((reacfr, metab, iso, cm, ciso));
                            visited.add((reacfr, metab, iso, cm, ciso));
                            # get products
                            v.update((reacfr,m,ist)
                                for (m,ist) in C13_ftbl.prod(metab,
                                iso, s, cm, ciso, cs, prods));
            #print "m=", metab, "v=\n", join("\n", v);##
            if not v:
                continue;
            # stock this results to demultiply
            # concerned paths
            res.append((ip, v));
        # deletion will start from the end so the following (lower) ips are
        # not concerned
        res.sort()
        res.reverse();
        for (ip,v) in res:
            cp=list(paths[ip]);
            del(paths[ip]);
            #print "ip=", ip, "cp=", cp;
            for item in v:
                tmp=list(cp);
                tmp.append(item);
                paths.insert(ip,tmp);
                #print "item=", item;##
                #print "paths=", paths;##

def front_frag(netan, paths, visited, frags):
    """Track frontal propagation of labeled fragments
    Get last fragment pattern from every path and add its neighbours if
    there is any eligible"""
    eligeable=True;
    #print "paths=", paths;
    #maxpath=500000;##
    step=0;
    #pdb.set_trace();##
    while eligeable:
        eligeable=False;
        res=[];
        #if len(paths) > maxpath:
        #    break;
        step+=1;
        # search for candidates and put fresh meat into frags dict
        cand=[];
        for (ip,p) in enumerate(paths):
            # get last item of this path
            #print "ip=", ip, "p=", p;##
            (reacfr, metab, frag)=p[-1];
            lset=frags.get(metab,{});
            if frag not in lset:
                lset[frag]=step;
            frags[metab]=lset;
            #print "front: m=", metab;##
            if metab in netan["output"]:
                continue;
            cand.append((ip, reacfr, metab, frag));
        for (ip,reacfr, metab, frag) in cand:
            # run through all concerned reactions and
            # all co-fragments to form candidate couples
            # (m,frag,s, cm,cfrag,cs)
            v=set();
            for (fr,src,prd) in (("fwd","left","right"),("rev","right","left")):
                for reac in netan["sto_m_r"][metab][src]:
                    if (src=="right") and (reac in netan["notrev"]):
                        # skip reverse of not reversible reaction
                        continue;
                    srcs=netan["carbotrans"][reac][src];
                    prods=netan["carbotrans"][reac][prd];
                    reacfr=reac+"."+fr;
                    #print "reacfr=", reacfr;##
                    mpos=[i for (i,(m,s)) in enumerate(srcs) if m==metab];
                    for inm in mpos:
                        (m,s)=srcs[inm];
                        (cm,cs)=("","") if len(srcs)==1 else srcs[(inm+1)%2];
                        cfrags=[frg for (frg,stp) in frags.get(cm,{}).iteritems()] or [""];
                        # use unlabeled cfrags in last resort
                        if (not cfrags[0]) and cm and cs:
                            cfrags=["z"*len(cs)];
                        for cfrag in cfrags:
                            # skip if already visited
                            if (reacfr, metab, frag, cm, cfrag) in visited:
                                continue;
                            eligeable=True;
                            visited.add((reacfr, metab, frag, cm, cfrag));
                            # get products
                            v.update((reacfr,m,frg)
                                for (m,frg) in C13_ftbl.frag_prod(metab,
                                frag, s, cm, cfrag, cs, prods));
            #print "m=", metab, "v=\n", join("\n", v);##
            if not v:
                continue;
            # stock this results to demultiply
            # concerned paths
            res.append((ip, v));
        # deletion will start from the end so the following (lower) ips are
        # not concerned
        res.sort()
        res.reverse();
        for (ip,v) in res:
            cp=list(paths[ip]);
            del(paths[ip]);
            #print "ip=", ip, "cp=", cp;
            for item in v:
                tmp=list(cp);
                tmp.append(item);
                paths.insert(ip,tmp);
                #print "item=", item;##
                #print "paths=", paths;##


sys.path.append('/home/sokol/dev/python');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;

#print "n=", __name__;
if sys.argv[0] != "":
    # not interactive session
    ftbl=C13_ftbl.ftbl_parse(sys.stdin);
    # check argument number
    if len(sys.argv) != 2:
        usage();
        exit(1);
    # check if metab is network
    isotop=sys.argv[1];
else:
   fname='essai.ftbl';
   fin=open(fname, "r");
   ftbl=C13_ftbl.ftbl_parse(fin);
   fin.close();
   isotop="PEP#111";

# analyze network
netan=C13_ftbl.ftbl_netan(ftbl);

(metab,isostr)=isotop.split("#");
metab in netan["metabs"] or [usage("Metabolite "+metab+\
    " is not in the network "+str(netan["metabs"])), exit(1)];

# check if its length matchs ftbl
len(isostr) == netan["Clen"][metab] or usage("Isotop "+isotop+" has wrong carbon length"+
    "\nExpecting "+str(netan["Clen"][metab])+", got "+str(len(isostr)));
# convert 01 string to int
#iso=sum(1<<i for (i,c) in enumerate(isostr[::-1]) if c=="1");
# convert "011" in "0bc" fragments
frag="".join((letters[i] if c=="1" else "0") for (i,c) in enumerate(isostr));

#print "isotop=", isotop, "m=", metab, "f=", frag;

# run through network to follow labeled carbons
# add visited isotopomers+context in visited and reac+isotop to path list
#paths=trace(netan, metab, isostr);
#print(paths);

# frontal exploration
#pdb.set_trace();
#paths=[[("", metab, iso)]];
paths=[[("", metab, frag)]];
#visited=list();
visited=set();
isos=dict();
frags=dict();
#front(netan, paths, visited, isos);
front_frag(netan, paths, visited, frags);
# print one path by row
maxl=max(len(p) for p in paths);
#print "maxl=", maxl;##

for p in paths:
    print ("\t"*(maxl-len(p)))+join("\t", p);

#print frags;##
print """\nAll products:\n
Metabolite
\tfragm\tstep
""";
[sys.stdout.write(str(metab)+"\n\t"+
    "\n\t".join([frag+"\t"+str(stp)
    for (frag,stp) in sorted(fset.iteritems())])+"\n")
    for (metab,fset) in sorted(frags.iteritems())];
#print "\nVisited nodes:"
#sys.stdout.write(join("\n", visited));
