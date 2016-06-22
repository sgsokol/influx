#!/usr/bin/env python

r"""
read a .ftbl file from stdin and translate it to .sif file for Cytoscape
on stdout. Two substrates or two products reactions are represented as
diamond nodes. Node attributes are written in node.<attribute>.fnetw file

usage: ./ftbl2rsif [-h|--help|--DEBUG] mynetwork[.ftbl]
Take an ftbl file and produce a .cif and attribute files for
cytoscape visualisation of the network defined in the ftbl.

OPTIONS
-h, --help print this message and exit
--DEBUG enable some debuggin features and output (for advanced users)

:param mynetwork: the base of an ftbl file (mynetwork.ftbl)

:returns: mynetwork.sif -- file of the network definition suitable for cytoscape

For attribute files, the names are self explanatory :
 * edge.label.mynetwork
 * edge.sourceArrowColor.mynetwork
 * edge.targetArrowColor.mynetwork
 * node.fillColor.mynetwork
 * node.shape.mynetwork

.. note::
 Base name of ftbl file ('mynetwork' in this example)
 is used to create or silently overwrite all result files list
"""

# 2008-01-23 sokol: initial release
# 2010-05-31 non blocking on ftbl parse errors

if __name__ == "__main__":
    import sys
    import os
    import getopt
    import re

    from tools_ssg import *
    from C13_ftbl import *

    werr=sys.stderr.write

    # determine colour of metabolite (in, out or plain)
    def color(m, netan, nc):
        return nc['i'] if m in netan['input'] \
            else nc['o'] if m in netan['output'] else \
            nc['m']

    #print sys.argv
    #exit()

    # Configurable constants
    ns={
        'm': 'roundrect',     # metabolite shape
        'r': 'diamond'        # reaction shape
    }
    nc={
        'm': '255,255,255',    # metabolite colour
        'i': '0,255,0',        # input colour (uptake)
        'o': '255,0,0',        # output colour (escape)
        'r': '204,0,255'       # reaction colour
    }
    et={
        'nr': 'ARROW',       # not reversible target
        'r': 'ARROW',        # reversible target
    }
    es={
        'nr': 'NONE',        # not reversible source
        'r': 'DIAMOND',      # reversible source
    }
    # arrow colours
    ec={
        'd': '0,0,255',    # dependent flux
        'c': '0,0,0',    # constrained flux
        'f': '0,255,0',    # free flux
    }
    ##print 'start'
    # take arguments
    #<--skip in interactive session
    # get arguments
    me=os.path.basename(sys.argv[0])
    def usage():
        sys.stderr.write(__doc__)
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "DEBUG"])
    except getopt.GetoptError, err:
        sys.stderr.write(str(err)+"\n")
        usage()
        sys.exit(1)
    cost=False
    DEBUG=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o=="--DEBUG":
            DEBUG=True
        else:
            assert False, "unhandled option"
    #aff("args", args);##
    if len(args) != 1:
        usage()
        exit(1)
    base=args[0]
    if base[-5:]==".ftbl":
        base=base[:-5]
    #-->
    # define where to read and write
    # input file is an argument
    fdir=os.path.dirname(base) or "."
    base=os.path.basename(base)
    ##    print base
    ##    print 'open files'
    fin=open(os.path.sep.join((fdir, base+".ftbl")), "r")
    fout=open(os.path.sep.join((fdir, base+".sif")), "w")
    fns=open(os.path.sep.join((fdir, "node.shape."+base)), "w")
    fnc=open(os.path.sep.join((fdir, "node.fillColor."+base)), "w")
    fes=open(os.path.sep.join((fdir, "edge.sourceArrowShape."+base)), "w")
    fet=open(os.path.sep.join((fdir, "edge.targetArrowShape."+base)), "w")
    fel=open(os.path.sep.join((fdir, "edge.label."+base)), "w")
    fesc=open(os.path.sep.join((fdir, "edge.sourceArrowColor."+base)), "w")
    fetc=open(os.path.sep.join((fdir, "edge.targetArrowColor."+base)), "w")
    # write headers for attributes
    fns.write("node.shape\n")
    fnc.write("node.fillColor\n")
    fes.write("edge.sourceArrowShape\n")
    fet.write("edge.targetArrowShape\n")
    fel.write("edge.label\n")
    fesc.write("edge.sourceArrowColor\n")
    fetc.write("edge.targetArrowColor\n")

    # Parse .ftbl file
    ##print 'parse'

    try:
       ftbl=ftbl_parse(fin)
    except Exception as inst:
       werr(str(inst)+"\n")
       raise

    #print ftbl
    ent="FLUXES"
    #aff("ftbl["+ent+"]", ftbl[ent], f=sys.stderr)

    # Analyse the network
    netan=dict()
    try:
        ftbl_netan(ftbl, netan)
    except:
        #werr(str(sys.exc_info()[1])+"\n")
        raise
        #sys.exit(1)
    
    ##print netan

    # create and write the graph
    #pdb.set_trace()
    ##print netan['notrev']
    for flux in ftbl['NETWORK']:
        reac=str(flux.get('FLUX_NAME'))

        if not len(reac): continue    # skip carbon section

        #path=re.sub(r'[0-9]+$', '', reac)

        # collect substrates and products
        subs=[str(flux['EDUCT_1'])]
        if flux.get('EDUCT_2'):
            subs.append(str(flux['EDUCT_2']))
        prods=[str(flux['PRODUCT_1'])]
        if flux.get('PRODUCT_2'):
            prods.append(str(flux['PRODUCT_2']))

        edges=[]
        if len(subs)==1 and len(prods)==1:
            # plain 1-to-1 reaction
            fout.write("%s %s %s\n" % (subs[0], reac, prods[0]))
            edges.append(subs[0]+" ("+reac+") "+prods[0])
        else:
            # 2-to-1, 1-to-2 or 2-to-2 reactions
            for s in subs:
                fout.write("%s %s %s\n" % (s, reac, reac))
                edges.append(s+" ("+reac+") "+reac)
            for p in prods:
                fout.write("%s %s %s\n" % (reac, reac, p))
                edges.append(reac+" ("+reac+") "+p)
        revers="r" if reac not in netan['notrev'] else "nr"
        net_dfc=("c" if reac in netan["flux_constr"]["net"] else
                 "f" if reac in netan["flux_free"]["net"] else
                 "d")
        xch_dfc=("c" if reac in netan["flux_constr"]["xch"] else
                 "f" if reac in netan["flux_free"]["xch"] else
                 "d")
    ##    print reac, revers
        if fns:
            # add node/edge shape and colour
            for m in subs+prods:
                fns.write("%s = %s\n" % (m, ns['m']))
                fnc.write("%s = %s\n" % (m, color(m, netan, nc)))
            if len(edges)>1:
                fns.write("%s = %s\n" % (reac, ns['r']))
                fnc.write("%s = %s\n" % (reac, nc['r']))
            else:
                fel.write(edges[0]+" = "+reac+"\n")
            for e in edges:
                fet.write(e+" = "+et[revers]+"\n")
                fes.write(e+" = "+es[revers]+"\n")
                fesc.write(e+" = "+ec[xch_dfc]+"\n")
                fetc.write(e+" = "+ec[net_dfc]+"\n")

    # close opened streams
    fin.close()
    fout.close()
    fns.close()
    fnc.close()
    fet.close()
    fes.close()
