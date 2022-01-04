#!/usr/bin/env python3
"""Parse ftbl file from first parameter
and write SBML in XML format on stdout
usage: ftbl2metxml.py network[.ftbl] [> network.xml]
"""
import sys, os, stat, re
from pathlib import Path
import argparse
from libsbml import *

import influx_si
import tools_ssg
import C13_ftbl
import kvh

me=os.path.basename(sys.argv[0])
def warn(mes):
    sys.stderr.write(me+": "+mes+"\n")
def werr(mes):
    raise Exception(me+": "+mes)
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
                 + ' LibSBML returned error code ' + str(value) + ': "' \
                 + OperationReturnValue_toString(value).strip() + '"'
       raise Exception(err_msg)
   else:
     return
def valid_id(s):
    "Replace invalid chars in s by '_' thus producing a valid id name"
    return re.sub(r"[^a-zA-Z0-9_.:-]", "_", str(s))
    
def main(argv=sys.argv[1:]):
    sys.tracebacklimit=None

    parser=argparse.ArgumentParser(description="convert FTBL to SBML. If a counterpart '_res.kvh' file is found, produce TSV files with flux values.")
    parser.add_argument('ftbl', nargs="+", help="file(s) to be converted to SBML. Extension '.ftbl' can be omitted.")
    parser.add_argument('--cstart', nargs=1, help="string delimiting start of compartment name, e.g. '--cstart=_'. If not given, all metabolites will be assigned to a default compartment")
    parser.add_argument('--cend', nargs=1, default="", help="string delimiting the end of compartment name, e.g. '--cstart=\"_[\" --cend=\"]\"', here a compartment name is expected to be in brackets preceded by an underscore. For example, with cited values, a metabolite name 'Glc_[c]' will be interpreted as metabolite 'Glc' in compartment 'c'.")
    parser.add_argument('--cdef', default="_def", help="default compartment name, e.g. '--cdef=cytoplasm'. Default is '_def'")
    parser.add_argument('--scale', type=float, default=1., help="If TSV files are generated, flux values are multiplied with this value.")

    args=parser.parse_args(argv)
    cstart=args.cstart
    if type(cstart) == list:
        cstart=cstart[0]
    cend=args.cend
    if type(cend) == list:
        cend=cend[0]
    cdef=args.cdef
    scale=args.scale
    fli=args.ftbl

    for fftbl in fli:
        if fftbl[-5:].lower() != ".ftbl":
            fftbl+=".ftbl"
        bftbl=fftbl[:-5] # base of ftbl name
        fftbl=Path(fftbl)
        if not fftbl.is_file():
            warn("file '"+str(fftbl)+"' does not exist.\n")
            return 1
        
        ftbl=C13_ftbl.ftbl_parse(str(fftbl))
        netan=dict()
        C13_ftbl.ftbl_netan(ftbl, netan, False, False)
        allmetabs=set()
        allmetabs.update(netan["input"])
        allmetabs.update(netan["metabs"])
        allmetabs.update(netan["output"])
        minout=netan["input"] | netan["output"]
        
        # sbml building starts here
        try:
            sbmlns = SBMLNamespaces(3,1,"groups",1)
            document = SBMLDocument(sbmlns)
        except:
            raise SystemExit('Could not create SBMLDocumention object')
        document.setPkgRequired('groups', False)
        model = document.createModel()
        check(model, 'create model')
        mplugin = model.getPlugin("groups")
        check(mplugin, "create plugin 'groups'")
        # get dict of species: {mc: (metab, compart)}
        dmc={}
        li_compart=set()
        if cstart:
            for mc in allmetabs:
                li=mc.split(cstart,1)
                if len(li) == 2 and li[1].endswith(cend):
                    # we have a metab+compartment
                    m=li[0]
                    c=li[1] if not cend else li[1][:-len(cend)]
                else:
                    m=mc
                    c=cdef
                dmc[mc]=(m, c)
                li_compart.add(c)
        else:
            li_compart.add(cdef)
            dmc=dict((mc,(mc,cdef)) for mc in allmetabs)
        # declare compartments
        for compart in li_compart:
            idc=valid_id(compart)
            c1=model.createCompartment()
            check(c1, "create compartement %s"%compart)
            check(c1.setId(idc), "set compartment 'id' attribute %s"%idc)
            check(c1.setName(str(compart)), "set compartment 'name' attribute %s"%compart)
            check(c1.setConstant(True), "set compartment 'constant' attribute %s"%compart)
        # declare species
        for mc,(m,c) in dmc.items():
            idm=valid_id(mc)
            idc=valid_id(c)
            s=model.createSpecies()
            check(s, "create specie %s"%mc)
            check(s.setId(idm), "set 'id' specie '%s'"%idm)
            check(s.setCompartment(idc), "set 'compartment' %s for specie '%s'"%(idc, mc))
            check(s.setName(str(m)), "set 'name' %s for specie '%s'"%(m, mc))
            check(s.setBoundaryCondition(mc in minout), "set 'bounaryCondition' for specie '%s'"%(mc,))
            check(s.setHasOnlySubstanceUnits(True), "set 'hasOnlySubstanceUnits' for specie '%s'"%(mc,))
            check(s.setConstant(mc in minout), "set 'constant' for specie '%s'"%(mc,))
            check(s.setInitialConcentration(1.), "set 'initialConcentration' for specie '%s'"%(mc,))
        # declare and fulfill pathways as groups
        for ip,(p,rli) in enumerate(netan["pathway"].items()):
            group=mplugin.createGroup()
            check(group, "create group for pathway '%s'"%p)
            check(group.setId("ptw_%d"%(ip+1)), "set id for pathway '%s'"%p)
            check(group.setKind("classification"), "set id for pathway '%s'"%p)
            check(group.setName(p), "set id for pathway '%s'"%p)
            # add reactions
            for r in rli:
                member=group.createMember()
                check(member, "create member for reaction '%s'"%r)
                check(member.setIdRef(valid_id(r)), "set idRef for reaction '%s'"%r)
                
        # declare reactions
        for reac,di in netan["sto_r_m"].items():
            vreac=valid_id(reac)
            r=model.createReaction()
            check(r, "create reaction %s"%reac)
            check(r.setId(vreac), "set reaction 'id' %s"%reac)
            check(r.setName(str(reac)), "set reaction 'name' %s"%reac)
            check(r.setReversible(reac not in netan["notrev"]), "set 'reversible' in reaction %s"%reac)
            check(r.setFast(False), "set reaction 'fast' %s"%reac)
            # create reactants (i.e. left part of reaction)
            for mc,co in di["left"]:
                vmc=valid_id(mc)
                mref=r.createReactant()
                check(mref, "create 'reactant' %s in reaction %s"%(mc,reac))
                check(mref.setSpecies(vmc), "assign 'reactant' %s in reaction %s"%(mc, reac))
                check(mref.setStoichiometry(co), "set 'stoichiometry' %f for reactant %s in reaction %s"%(co, mc, reac))
                check(mref.setConstant(mc in minout), "set 'constant' for reactant %s in reaction %s"%(mc, reac))
            # create products (i.e. right part of reaction)
            for mc,co in di["right"]:
                vmc=valid_id(mc)
                mref=r.createProduct()
                check(mref, "create product %s in reaction %s"%(mc,reac))
                check(mref.setSpecies(vmc), "assign 'product' %s in reaction %s"%(mc, reac))
                check(mref.setStoichiometry(co), "set stoichiometry %f for product %s in reaction %s"%(co, mc, reac))
                check(mref.setConstant(mc in minout), "set 'constant' for product %s in reaction %s"%(mc, reac))

        # write final sbml
        writeSBML(document, str(fftbl.with_suffix(".xml")))
        # get fluxes from kvh if present
        fkvh=Path(bftbl+"_res.kvh")
        if not fkvh.is_file():
            warn("warning: file '%s' is not found. No flux values are written in TSV files.\n"%str(fkvh))
            continue
        dkvh=kvh.kvh2dict(str(fkvh), strip=True)
        if not "linear stats" in dkvh or not "fwd-rev fluxes (sorted by name)" in dkvh["linear stats"]:
            warn("warning: field 'linear stats/fwd-rev fluxes (sorted by name)' is not found in file '%s'. No flux values are written in TSV files.\n"%str(fkvh))
            continue
        table=dkvh["linear stats"]["fwd-rev fluxes (sorted by name)"]
        fwdfl=dict((fl[4:], v.split("\t")[0]) for fl,v in table.items() if fl[:4] == "fwd.")
        revfl=dict((fl[4:], v.split("\t")[0]) for fl,v in table.items() if fl[:4] == "rev.")
        with open(bftbl+"_net.txt", "w") as f:
            f.write("Identifier\t%s\n"%bftbl)
            f.write("\n".join("%s\t%f"%(fl,(float(v)-float(revfl[fl]))*scale) for fl,v in fwdfl.items())+"\n")
        with open(bftbl+"_fwd.txt", "w") as f:
            f.write("Identifier\t%s\n"%bftbl)
            f.write("\n".join("%s\t%f"%(fl,float(v)*scale) for fl,v in fwdfl.items())+"\n")
        with open(bftbl+"_rev.txt", "w") as f:
            f.write("Identifier\t%s\n"%bftbl)
            f.write("\n".join("%s\t%f"%(fl,float(v)*scale) for fl,v in revfl.items())+"\n")
    return 0
if __name__ == "__main__":
    sys.exit(main())
