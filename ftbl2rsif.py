#!/usr/bin/python
# read a .ftbl file from stdin and translate it to .sif file for Cytoscape
# on stdout. Two substrates or two products reactions are represented as
# diamond nodes. Node attributes are written in node.<attribute>.fnetw file
# usage: ./ftbl2rsif.py < fnetw.ftbl > fnetw.sif
# or: ./ftbl2rsif.py fnetw.ftbl (in this case node.<attribute>.fnetw are silently
# rewritten)
# 
# 2008-01-23 sokol

#import pdb;
import sys;
sys.path.append('/home/sokol/dev/python');
import re;
from tools_ssg import *;

from C13_ftbl import *;

# determin colour of metabolite (in, out or plain)
def color(m, netan, nc):
    return nc['i'] if m in netan['input'] \
        else nc['o'] if m in netan['output'] else \
        nc['m'];

#print sys.argv;
#exit();

# Configurable constants
ns={
    'm': 'roundrect',    # metabolite shape
    'r': 'diamond'        # reaction shape
};
nc={
    'm': '255,255,255',    # metabolite colour
    'i': '0,255,0',        # input colour (uptake)
    'o': '255,0,0',        # output colour (escape)
    'r': '204,0,255'    # reaction colour
};
et={
    'nr': 'ARROW',        # not reversible target
    'r': 'ARROW',        # reversible target
};
es={
    'nr': 'NONE',        # not reversible source
    'r': 'DIAMOND',        # reversible source
};
##print 'start'
# decide where to read and write
if len(sys.argv) == 2:
    # input file is an argument
    base=re.sub(".ftbl$", "", sys.argv[1]);
##    print base;
    if base==sys.argv[1]:
        sys.stderr.write(sys.argv[0]+": The only argument must be a .ftbl file name, got \
'"+sys.argv[1]+"'");
        raise "FileTypeError";
##    print 'open files'
    fin=open(sys.argv[1], "r");
    fout=open(base+".sif", "w");
    fns=open("node.shape."+base, "w");
    fnc=open("node.fillColor."+base, "w");
    fes=open("edge.sourceArrowShape."+base, "w");
    fet=open("edge.targetArrowShape."+base, "w");
    fel=open("edge.label."+base, "w");
    # write headers for attributes
    fns.write("node.shape\n");
    fnc.write("node.fillColor\n");
    fes.write("edge.sourceArrowShape\n");
    fet.write("edge.targetArrowShape\n");
    fel.write("edge.label\n");
elif len(sys.argv) == 1:
    # standart input and output are used
    fin=sys.stdin;
    fout=sys.stdout;
    fns=fnc=fes=fet=fel=0;
    
# Parse .ftbl file
##print 'parse'
ftbl=ftbl_parse(fin);
##print ftbl
ent="FLUXES";
#aff("ftbl["+ent+"]", ftbl[ent], f=sys.stderr);

# Analyse the network
netan=ftbl_netan(ftbl);
##print netan;

# create and write the graph
#pdb.set_trace();
##print netan['notrev']
for flux in ftbl['NETWORK']:
    reac=str(flux.get('FLUX_NAME'));
    
    if not len(reac): continue    # skip carbon section
    
    path=re.sub(r'[0-9]+$', '', reac);
    
    # collect substrates and products
    subs=[str(flux['EDUCT_1'])];
    if flux.get('EDUCT_2'):
        subs.append(str(flux['EDUCT_2']));
    prods=[str(flux['PRODUCT_1'])];
    if flux.get('PRODUCT_2'):
        prods.append(str(flux['PRODUCT_2']));
    
    edges=[];
    if len(subs)==1 and len(prods)==1:
        # plain 1-to-1 reaction
        fout.write("%s %s %s\n" % (subs[0], reac, prods[0]));
        edges.append(subs[0]+" ("+reac+") "+prods[0]);
    else:
        # 2-to-1, 1-to-2 or 2-to-2 reactions
        for s in subs:
            fout.write("%s %s %s\n" % (s, path, reac));
            edges.append(s+" ("+path+") "+reac);
        for p in prods:
            fout.write("%s %s %s\n" % (reac, path, p));
            edges.append(reac+" ("+path+") "+p);
    revers="r" if reac not in netan['notrev'] else "nr";
##    print reac, revers
    if fns:
        # add node/edge shape and colour
        for m in subs+prods:
            fns.write("%s = %s\n" % (m, ns['m']));
            fnc.write("%s = %s\n" % (m, color(m, netan, nc)));
        if len(edges)>1:
            fns.write("%s = %s\n" % (reac, ns['r']));
            fnc.write("%s = %s\n" % (reac, nc['r']));
        else:
            fel.write(edges[0]+" = "+reac+"\n");
        for e in edges:
            fet.write(e+" = "+et[revers]+"\n");
            fes.write(e+" = "+es[revers]+"\n");

# close opened streams
fin.close();
fout.close();
if fns:
    fns.close();
    fnc.close();
    fet.close();
    fes.close();
