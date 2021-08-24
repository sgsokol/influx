#!/usr/bin/env python3
"""Parse ftbl file from from first parameter or from stdin (if input file is '-')
and write a series of mtf (multiple TSV files).
The file stem ('network' in 'network.ftbl') is used as file name basis
for produced files, e.g. 'network.miso'. Parameter --out can be used to change it.
If out path includes non existing directories, they are automatically created.
Caution! If an existing output file contains "# Created by 'ftbl2mft" or is empty, it is silently overwritten.
Otherwise, the writing is aborted with a warning. Other files may continue to be created.
To force the overwriting, use '--force'.

Output files will have following extensions/meanings:
 .netw: stoechiometric equations and carbon transitions in the meatabolic network;
 .linp: label input;
 .miso: isotopic measurements (MS, label, peak);
 .mflux: flux measurements;
 .mmet: metabolic concentration measurements;
 .tvar: flux/metabolite types partition (free, dependent, constrained) and starting values;
 .cnstr: constraints (equalities, inequalities for both fluxes and concentrations);
 .opt: options.

Copyright 2021 INRAE, INSA, CNRS
Author: Serguei Sokol (sokol [at] insa-toulouse [dot] fr)

usage: ftbl2mtf.py [--out OUT] [--inst] network[.ftbl]
"""
import sys, os, getopt, stat, io
import datetime
from pathlib import Path
import numpy as np

import influx_si
import tools_ssg as tls
import C13_ftbl
from txt2ftbl import tsv2df

LOCAL_TIMEZONE=datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo;
me=os.path.basename(sys.argv[0] or "ftbl2mtf")

def dtstamp():
    "formatted date-time stamp"
    return datetime.datetime.now(LOCAL_TIMEZONE).strftime('%Y-%m-%d %H:%M:%S %Z %z')
def usage():
    print(__doc__)
def warn(mes):
    sys.stderr.write(f"{me}: "+mes+"\n")
def werr(mes):
    raise Exception(f"{me}: "+mes)
def main(argv=sys.argv[1:]):
    # init values
    #sys.tracebacklimit=None
    fftbl=""
    out=""
    case_i=False
    force=False
    invcomp={">=": "<=", "=>": "<=", "<=": ">=", "=<": ">="}
    # parse options
    opts,args=getopt.getopt(argv, "hio:", ["help", "force", "inst", "out="])
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "--out" or o == "-o":
            out=a
        elif o == "--inst":
            case_i=True;
        elif o == "--force":
            force=True
    if not args:
        usage()
        werr("Expecting ftbl file name")
        
    fftbl=Path(args[0]) if args[0] != "-" else sys.stdin
    if type (fftbl) == type(Path()) and not fftbl.is_file():
        fftbl = fftbl.with_suffix(".ftbl")
        if not fftbl.is_file():
            raise Exception(me+": file '"+str(fftbl)+"' does not exist.\n")
    # prepare out
    if out:
        out=Path(out)
    else:
        if type (fftbl) == type(Path()):
            out=fftbl.parent/fftbl.stem
        else:
            out=Path("mtf")
    ftbl=C13_ftbl.ftbl_parse(str(fftbl))
    #print("ftbl parsed=", ftbl)
    #import pdb; pdb.set_trace()
    netan=dict()
    C13_ftbl.ftbl_netan(ftbl, netan)
    bsl="\\" # backslash
    scre=f"# Created by '{me} {' '.join(v.replace(' ', bsl+' ') for v in argv)}'"
    for suff in (".netw", ".linp", ".miso", ".mflux", ".mmet", ".tvar", ".cnstr", ".opt"):
        cout=out.with_suffix(suff) # current output
        dloc={} # local dictionary
        res=""
        header=""
        if suff == ".netw":
            # metab network
            ltr=ftbl["long_trans"]
            for rnm,reac in ftbl["long_reac"].items():
                for lr in ("left", "right"):
                    dloc[lr]=" + ".join(f"{str(tmet[1])+'*' if tmet[1] != 1 else ''}{tmet[0]} {'('+ltr[rnm][lr][i][1:]+')' if ltr[rnm][lr][i] else ''}" for i,tmet in enumerate(reac[lr]))
                s=f"{rnm}:\t{dloc['left']} {'<' if rnm not in netan['notrev'] else ''}-> {dloc['right']}\n"
                res += s
        elif suff == ".linp":
            # label input
            header="Id\tComment\tMetabolite\tIsotopomer\tValue\n"
            for d in ftbl["LABEL_INPUT"]:
                met=d["META_NAME"] if d["META_NAME"] else met
                res += f"\t\t{met}\t{d['ISOTOPOMER'][1:]}\t{d['VALUE']}\n"
        elif suff == ".miso":
            # label, peak, ms measurements
            header="Id\tComment\tMetabolite\tFragment\tDataset\tSpecies\tValue\tSD\tTime\n"
            if case_i:
                # pick file_labcin from opt
                #import pdb; pdb.set_trace()
                flabcin=[d["OPT_VALUE"] for d in ftbl.get("OPTIONS", []) if d["OPT_NAME"] == "file_labcin"]
                if not flabcin:
                    werr("option ''--inst' was activated but a field 'file_labcin' was not found in OPTIONS in '%s'"%fftbl.name)
                if len(flabcin) > 1:
                    werr("option ''--inst' was activated but a field 'file_labcin' is not unique in OPTIONS in '%s'"%fftbl.name)
                flabcin=fftbl.parent/flabcin[0]
                # read flabcin
                try:
                    # default comment='#'
                    df_cin=tsv2df(flabcin, append_iline=None)
                except:
                    df_cin=tsv2df(flabcin, comment="//", append_iline=None)
                # get netan measurements
                vrc=df_cin.iloc[:, 0].to_numpy().astype(str)
                if "measures" not in netan:
                    measures=dict()
                    for meas in ("label", "mass", "peak"):
                        measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas)
            # label
            if case_i:
                # pick "l:..." in df_cin
                ir=np.where(np.char.startswith(vrc, "l:"))[0]
                last_met=last_cgr=""
                labset=0
                for i in ir:
                    let,met,cgr,li=vrc[i].split(":")
                    if met != last_met or cgr != last_cgr:
                        labset += 1
                        last_met,last_cgr=met,cgr
                    # get sd from netan
                    # im - index of the measurement in netan structures
                    try:
                        im=measures["label"][0]["ids"].index(vrc[i])
                    except:
                        im=[ii for ii,v in enumerate(measures["label$"][0]["ids"]) if v.startswith(f"l:{met}:{cgr}:")]
                        if len(im) != 1:
                            continue # not found corresponding unique FTBL entry, silently skip it
                        im=im[0]
                    sdv=measures["label"][0]["dev"][im]
                    for ti in df_cin.columns[1:]:
                        res += f"\t\t{met}\t\tLAB-{labset}\t{cgr.replace('#', '')}\t{df_cin.loc[i, ti]}\t{sdv}\t{ti}\n"
            else:
                labset=0
                for d in ftbl["LABEL_MEASUREMENTS"]:
                    met=d["META_NAME"] if d["META_NAME"] else met
                    if d["META_NAME"]:
                        labset += 1
                        cgr=d["CUM_GROUP"]
                    elif cgr != d["CUM_GROUP"]:
                        labset += 1
                    res += f"\t\t{met}\t\tLAB-{labset}\t{d['CUM_CONSTRAINTS'].replace('#', '')}\t{d['VALUE']}\t{d['DEVIATION']}\t\n"
            # peak
            if case_i:
                # pick "p:..." in df_cin
                ir=np.where(np.char.startswith(vrc, "p:"))[0]
                last_met=last_atom=""
                pset=0
                for i in ir:
                    let,met,atom,p,li=vrc[i].split(":")
                    ps=int(atom)
                    pdm=str(ps-1)
                    pdp=str(ps+1)
                    ps=str(ps)
                    if ps == "1":
                        frag="1,2"
                    else:
                        frag="%s,%s,%s"%(pdm,ps,pdp) # todo: make sure that pdp is not out of molecule
                    if met != last_met or atom != last_atom:
                        pset += 1
                        last_met,last_atom=met,atom
                    # get sd from netan
                    # im - index of the measurement in netan structures
                    try:
                        im=measures["peak"][0]["ids"].index(vrc[i])
                    except:
                        im=[ii for ii,v in enumerate(measures["peak"][0]["ids"]) if v.startswith(f"p:{met}:{atom}:{p}:")]
                        if len(im) != 1:
                            continue # not found corresponding unique FTBL entry, silently skip it
                        im=im[0]
                    sdv=measures["peak"][0]["dev"][im]
                    spec=ps+"->"+(pdm if p == "D-" else pdp if p == "D+" else "" if p == "S" else pdm+","+pdp)
                    for ti in df_cin.columns[1:]:
                        res += f"\t\t{met}\t{frag}\tPEAK-{pset}\t{spec}\t{df_cin.loc[i, ti]}\t{sdv}\t{ti}\n"
            else:
                pset=0
                for d in ftbl["PEAK_MEASUREMENTS"]:
                    met=d["META_NAME"] if d["META_NAME"] else met
                    pset += 1
                    ps=int(d["PEAK_NO"])
                    pdm=str(ps-1)
                    pdp=str(ps+1)
                    ps=str(ps)
                    if ps == "1":
                        frag="1,2"
                    else:
                        frag="%s,%s,%s"%(pdm,ps,pdp) # todo: make sure that pdp is not out of molecule
                    for p in ("S", "D-", "D+", "DD", "T"):
                        val=d.get("VALUE_"+p)
                        if not val:
                            continue
                        spec=ps+"->"+(pdm if p == "D-" else pdp if p == "D+" else "" if p == "S" else pdm+","+pdp)
                        # sd from dict
                        dsdv=d.get("DEVIATION_"+p) if p not in ("DD", "T") else d.get("DEVIATION_DD/T")
                        # current sd
                        sdv=dsdv or sdv
                        res += f"\t\t{met}\t{frag}\tPEAK-{pset}\t{spec}\t{val}\t{sdv}\t\n"
            # ms
            if case_i:
                # pick "m:..." in df_cin
                ir=np.where(np.char.startswith(vrc, "m:"))[0]
                last_met=last_frag=""
                mset=0
                #import pdb; pdb.set_trace()
                for i in ir:
                    let,met,frag,w,li=vrc[i].split(":")
                    if met != last_met or frag != last_frag:
                        mset += 1
                        last_met,last_frag=met,frag
                    # get sd from netan
                    # im - index of the measurement in netan structures
                    try:
                        im=measures["mass"][0]["ids"].index(vrc[i])
                    except:
                        im=[ii for ii,v in enumerate(measures["mass"][0]["ids"]) if v.startswith(f"m:{met}:{frag}:{w}")]
                        if len(im) != 1:
                            continue # not found corresponding unique FTBL entry, silently skip it
                        im=im[0]
                    sdv=measures["mass"][0]["dev"][im]
                    for ti in df_cin.columns[1:]:
                        res += f"\t\t{met}\t{frag.replace('~', '-')}\tMS-{mset}\tM{w}\t{df_cin.loc[i, ti]}\t{sdv}\t{ti}\n"
            else:
                mset=0
                for d in ftbl["MASS_SPECTROMETRY"]:
                    met=d["META_NAME"] if d["META_NAME"] else met
                    frag=d["FRAGMENT"] if mset == 0 or d["FRAGMENT"] else frag
                    if d["META_NAME"] or d["FRAGMENT"]:
                        mset += 1
                    res += f"\t\t{met}\t{frag.replace('~', '-')}\tMS-{mset}\tM{d['WEIGHT']}\t{d['VALUE']}\t{d['DEVIATION']}\t\n" # todo: case_i -> true time
        elif suff == ".mflux":
            # flux measurements
            header="Id\tComment\tFlux\tValue\tSD\n"
            for d in ftbl["FLUX_MEASUREMENTS"]:
                res += f"\t\t{d['FLUX_NAME']}\t{d['VALUE']}\t{d['DEVIATION']}\n"
        elif suff == ".mmet":
            # metabolite measurements
            header="Id\tComment\tMetabolite\tValue\tSD\n"
            for d in ftbl.get("METAB_MEASUREMENTS", []):
                res += f"\t\t{d['META_NAME']}\t{d['VALUE']}\t{d['DEVIATION']}\n"
        elif suff == ".tvar":
            # type of variable and starting values
            header="Id\tComment\tName\tKind\tType\tValue\n"
            for nx in ("NET", "XCH"):
                for d in ftbl["FLUXES"][nx]:
                    if nx == "XCH" and ((d["FCD"] == "C" and float(d['VALUE(F/C)']) == 0.) or (d["FCD"] == "D" and d["NAME"] in netan["flux_inout"])):
                        continue
                    res += f"\t\t{d['NAME']}\t{nx}\t{d['FCD']}\t{d['VALUE(F/C)']}\n"
            for d in ftbl.get("METABOLITE_POOLS", []):
                val=float(d['META_SIZE'])
                fcd="F" if val < 0 else "C"
                res += f"\t\t{d['META_NAME']}\tMETAB\t{fcd}\t{abs(val)}\n"
        elif suff == ".cnstr":
            # constraints (eqs, ineqs)
            header="Id\tComment\tKind\tFormula\tOperator\tValue\n"
            for it in ("EQUALITIES", "INEQUALITIES"):
                for nx in ("NET", "XCH", "METAB"):
                    for d in ftbl.get(it, {}).get(nx, []):
                        op="==" if it == "EQUALITIES" else invcomp[d['COMP']]
                        res += f"\t\t{nx}\t{d['FORMULA']}\t{op}\t{d['VALUE']}\n"
        elif suff == ".opt":
            # options
            header="Id\tComment\tName\tValue\n"
            for d in ftbl["OPTIONS"]:
                res += f"\t\t{d['OPT_NAME']}\t{d['OPT_VALUE']}\n"
        if res:
            # avoid writing empty files
            if not force and cout.is_file() and cout.stat().st_size > 0:
                # check if we can overwrite
                with cout.open() as fc:
                     if scre[:22] != fc.read(22):
                         warn(f"cannot overwrite '{fc.name}' as not created by this script. Use '--force' to go through.")
                         continue
            if not cout.parent.exists():
                cout.parent.mkdir(parents=True)
            print("cout=", str(cout))
            cout.write_text(f"{scre} at {dtstamp()}\n"+header+res)
    return 0

if __name__ == "__main__" or __name__ == "influx_si.cli":
    main()
