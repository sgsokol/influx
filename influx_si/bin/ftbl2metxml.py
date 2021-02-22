#!/usr/bin/env python3
"""Parse ftbl file from first parameter
and write SBML in XML format on stdout
usage: ftbl2metxml.py network[.ftbl] [> network.xml]
"""
def check(value, message):
   """If 'value' is None, prints an error message constructed using
   'message' and then exits with status code 1.  If 'value' is an integer,
   it assumes it is a libSBML return status code.  If the code value is
   LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
   prints an error message constructed using 'message' along with text from
   libSBML explaining the meaning of the code, and exits with status code 1.
   """
   if value == None:
     raise Exception('LibSBML returned a null value trying to ' + message + '.')
   elif type(value) is int:
     if value == LIBSBML_OPERATION_SUCCESS:
       return
     else:
       err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + OperationReturnValue_toString(value).strip() + '"'
       raise Exception(err_msg)
   else:
     return
if __name__ == "__main__" or __name__ == "influx_si.cli":
    import sys, os, getopt, stat, re
    from libsbml import *

    sys.tracebacklimit=None
    import influx_si
    import tools_ssg
    import C13_ftbl
    import kvh
    def usage():
        print(__doc__)

    me=os.path.basename(sys.argv[0])
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        else:
            assert False, "unhandled option"
    if not args:
        sys.stderr("Expecting ftbl file name\n")
        usage()
        
    fftbl=args[0]
    if fftbl and fftbl[-5:] != ".ftbl":
        fftbl+=".ftbl"
    if not os.path.exists(fftbl):
        sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n")
        sys.exit(1)

    # what kind of output we have?
    mode=os.fstat(1).st_mode
    f=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(fftbl[:-4]+"xml", "w")
    
    ftbl=C13_ftbl.ftbl_parse(fftbl)
    netan=dict()
    C13_ftbl.ftbl_netan(ftbl, netan, False, False)
    allmetabs=set()
    allmetabs.update(netan["input"])
    allmetabs.update(netan["metabs"])
    allmetabs.update(netan["output"])
    
    # sbml building starts here
    try:
        document = SBMLDocument(3, 1)
    except ValueError:
        raise SystemExit('Could not create SBMLDocumention object')
    model = document.createModel()
    check(model, 'create model')
    
    # get list of compartments, i.e. the last "_smth": (idc,c)
    li_compart=set((("_" if c[0].isdigit() else "")+c,c) for m in allmetabs for li in (m.split("_"),) if len(li) > 1 for c in (li[len(li)-1],))
    li_compart.add(("_def", "_def")) # add default compartment
    for idc,compart in li_compart:
        c1=model.createCompartment()
        check(c1, "create compartement %s"%compart)
        check(c1.setId(str(idc)), "set compartment id %s"%idc)
        check(c1.setName(str(compart)), "set compartment name %s"%compart)
    # get dict of species: {mc: (metab, compart)}
    dmc=dict((mc, (m, c)) for mc in allmetabs for li in (mc.split("_"),) for (m,c) in ((("_".join(li[:len(li)-1]),li[len(li)-1]) if len(li) > 1 else (mc, "_def")),))
    # add id which can be different from mc if starts with a digit
    dmc=dict((mc,(m,c,("_" if mc[0].isdigit() else "")+mc, ("_" if c[0].isdigit() else "")+c) ) for mc,(m,c) in dmc.items())
    for mc,(m,c,idm,idc) in dmc.items():
        s=model.createSpecies()
        check(s, "create specie %s"%mc)
        check(s.setId(str(idm)), "set id specie %s"%idm)
        check(s.setCompartment(str(idc)), "set compartment %s for specie %s"%(idc, mc))
        check(s.setName(str(m)), "set name %s for specie %s"%(m, mc))
    # create reactions
    for reac,di in netan["sto_r_m"].items():
        r=model.createReaction()
        check(r, "create reaction %s"%r)
        check(r.setId(str(reac)), "set reaction id %s"%reac)
        check(r.setName(str(reac)), "set reaction name %s"%reac)
        check(r.setReversible(reac not in netan["notrev"]), "set reversibility in reaction %s"%reac)
        # create reactants (i.e. left part of reaction)
        for mc,co in di["left"]:
            mref=r.createReactant()
            check(mref, "create reactant %s in reaction %s"%(mc,reac))
            check(mref.setSpecies(str(dmc[mc][2])), "assign reactant %s in reaction %s"%(mc, reac))
            check(mref.setStoichiometry(co), "set stoichiometry %f for reactant %s in reaction %s"%(co, mc, reac))
        # create products (i.e. right part of reaction)
        for mc,co in di["right"]:
            mref=r.createProduct()
            check(mref, "create product %s in reaction %s"%(mc,reac))
            check(mref.setSpecies(str(dmc[mc][2])), "assign product %s in reaction %s"%(mc, reac))
            check(mref.setStoichiometry(co), "set stoichiometry %f for product %s in reaction %s"%(co, mc, reac))
    # add pathway tags
    docu=[]
    for l in writeSBMLToString(document).split("\n"):
        docu.append(l)
        g=re.match('( *)<reaction id="([^"]+)"[^>]*>$', l)
        if g:
            pathway=netan["reac2path"].get(g.groups()[1])
            if pathway:
                docu.append('%s  <pathway id="%s" name="%s"/>'%(g.groups()[0], pathway, pathway))
    # write final sbml
    f.write("\n".join(docu)+"\n")
    f.close()
    
    # get fluxes from kvh if present
    fkvh=fftbl[:-5]+"_res.kvh"
    if not os.path.exists(fkvh):
        sys.stderr.write("Warning: file '%s' is not found. No flux value is written in txt files.\n"%fkvh)
        sys.exit(0)
    dkvh=kvh.kvh2dict(fkvh, strip=True)
    if not "linear stats" in dkvh or not "fwd-rev fluxes (sorted by name)" in dkvh["linear stats"]:
        sys.stderr.write("Warning: field 'linear stats/fwd-rev fluxes (sorted by name)' is not found in file '%s'. No flux value is written in txt files.\n"%fkvh)
        sys.exit(0)
    table=dkvh["linear stats"]["fwd-rev fluxes (sorted by name)"]
    fwdfl=dict((fl[4:], v.split("\t")[0]) for fl,v in table.items() if fl[:4] == "fwd.")
    revfl=dict((fl[4:], v.split("\t")[0]) for fl,v in table.items() if fl[:4] == "rev.")
    f=open(fftbl[:-5]+"_net.txt", "w")
    f.write("\n".join(fl+"\t"+str(float(v)-float(revfl[fl])) for fl,v in fwdfl.items())+"\n")
    f.close()
    f=open(fftbl[:-5]+"_fwd.txt", "w")
    f.write("\n".join("\t".join((fl,v)) for fl,v in fwdfl.items())+"\n")
    f.close()
    f=open(fftbl[:-5]+"_rev.txt", "w")
    f.write("\n".join("\t".join((fl,v)) for fl,v in revfl.items())+"\n")
    f.close()
    sys.exit(0)
