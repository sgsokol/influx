#!/usr/bin/python
# 2008-09-01 sokol
# follow the path of labeled carbons starting
# from a given isotopomer
# usage: isotrace.py METAB#0010 < file.ftbl
# metabolite METAB must be present in .ftbl file given on stdinput
# and its isotopomer #0010 must mutch carbon number defined in .ftbl
# Output of the script is a tab separateb list of isotopomer pathes taken by labeled
# carbon(s) till network output, one row by path

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
    maxpath=1000;##
    while eligeable:
        eligeable=False;
        res=[];
        if len(paths) > maxpath:
            break;
        for (ip,p) in enumerate(paths):
            # get last item of this path
            #print "ip=", ip, "p=", p;##
            (reacfr, metab, isostr)=p[-1];
            lset=isos.get(metab,set());
            lset.add(isostr);
            isos[metab]=lset;
            #print "front: m=", metab;##
            if metab in netan["output"]:
                continue;
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
                        for ciso in isos.get(cm,set(("",))):
                            # skip if already visited
                            if (reacfr, metab, isostr, cm, ciso) in visited:
                                continue;
                            eligeable=True;
                            # add coisotop to catalog
                            if cm and ciso:
                                clset=isos.get(cm,set());
                                clset.add(ciso);
                                isos[cm]=clset;
                            # register this visit
                            visited.append((reacfr, metab, isostr, cm, ciso));
                            # get products
                            v.update((reacfr,m,istr)
                                for (m,istr) in C13_ftbl.prod(metab,
                                isostr, s, cm, ciso, cs, prods));
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

# check argument number
if len(sys.argv) != 2:
    usage();
    exit(1);

sys.path.append('/home/sokol/dev/python');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;

#fname='ex5.ftbl';
#fin=open(fname, "r");
#ftbl=C13_ftbl.ftbl_parse(fin);
ftbl=C13_ftbl.ftbl_parse(sys.stdin);
#fin.close();

# analyze network
netan=C13_ftbl.ftbl_netan(ftbl);

# check if metab is network
isotop=sys.argv[1];
#isotop="A#01";
(metab,isostr)=isotop.split("#");
metab in netan["metabs"] or [usage("Metabolite "+metab+\
    " is not in the network "+str(netan["metabs"])), exit(1)];

# check if its length matchs ftbl
len(isostr) == netan["Clen"][metab] or usage("Isotop "+isotop+" has wrong carbon length"+
    "\nExpecting "+str(netan["Clen"][metab])+", got "+str(len(isostr)));

# run through network to follow labeled carbons
# add visited isotopomers in path list
#paths=trace(netan, metab, isostr);
#print(paths);

# frontal exploration
paths=[[("", metab, isostr)]];
visited=list();
isos=dict();
front(netan, paths, visited, isos);
# print one path by row
maxl=max(len(p) for p in paths);
#print "maxl=", maxl;##
for p in paths:
    print ("\t"*(maxl-len(p)))+join("\t", p);
print "\nAll products:"
[sys.stdout.write(str(item)+"\n\t"+
    "\n\t".join(v for v in sorted(val))+"\n")
    for (item,val) in sorted(isos.iteritems())];
print "\nVisited nodes:"
sys.stdout.write(join("\n", visited));
