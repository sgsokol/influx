#!/usr/bin/env python3
"""Parse ftbl file from first parameter or from stdin (if input file is '-')
and write a series of mtf (multiple TSV files).
The file stem ('network' in 'network.ftbl') is used as file name basis
for produced files, e.g. 'network.miso'. Parameter --out can be used to change it.
If out path includes non existing directories, they are automatically created.
Caution! If an existing output file starts with a comment
"# Created by 'ftbl2mft ..."
or is empty, it is silently overwritten.
Otherwise, the writing is aborted with a warning. Other files may continue to be created.
To force the overwriting, use '--force'.

Output files will have following extensions/meanings:
 .netw: stoichiometric equations and label transitions in the biochemical network;
 .linp: label input;
 .miso: isotopic measurements (MS, label, peak);
 .mflux: flux measurements;
 .mmet: biochemical specie concentration measurements;
 .tvar: flux/specie types partition (free, dependent, constrained) and starting values;
 .cnstr: constraints (equalities, inequalities for both fluxes and concentrations);
 .opt: options.

Copyright 2022 INRAE, INSA, CNRS
Author: Serguei Sokol (sokol [at] insa-toulouse [dot] fr)
"""
import sys, os, argparse, stat, io
import datetime
from pathlib import Path
import numpy as np

import influx_si
import tools_ssg as tls
import C13_ftbl
from txt2ftbl import tsv2df, try_ext, plain_natural_key

LOCAL_TIMEZONE=datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo;
me=os.path.basename(sys.argv[0] or "ftbl2mtf")

def dtstamp():
    "formatted date-time stamp"
    return datetime.datetime.now(LOCAL_TIMEZONE).strftime('%Y-%m-%d %H:%M:%S %Z %z')
def usage():
    print(__doc__)
def warn(mes):
    sys.stderr.write(f"Warning! {me}: "+mes+"\n")
def werr(mes):
    raise Exception(f"{me}: "+mes)
def ftbl2suff(ftbl, fftbl, case_i, netan, force, out, scre, suffs):
    invcomp={">=": "<=", "=>": "<=", "<=": ">=", "=<": ">="}
    out.parent.mkdir(parents=True, exist_ok=True)
    for suff in suffs:
        cout=out.with_suffix(suff) # current output
        dloc={} # local dictionary
        res=""
        header=""
        if suff == ".netw":
            # metab network
            #import pdb; pdb.set_trace()
            ltr=ftbl["long_trans"]
            for rnm,reac in ftbl["long_reac"].items():
                for lr in ("left", "right"):
                    dloc[lr]=" + ".join(f"{str(tmet[1])+'*' if tmet[1] != 1 else ''}{tmet[0]} {'('+ltr[rnm][lr][i][1:]+')' if ltr[rnm][lr][i] else ''}" for i,tmet in enumerate(reac[lr]))
                s=f"{rnm}:\t{dloc['left']} {'<' if rnm not in netan['notrev'] else ''}-> {dloc['right']}\n"
                res += s
            # reactions without label transitions
            # take all if no growth fluxes otherwise exclude ending with "_gr"
            woltr=sorted(set(netan["sto_r_m"].keys() if netan["opt"].get("include_growth_flux", 0) != 1 else (k for k in netan["sto_r_m"].keys() if not k.endswith("_gr")))-set(ftbl["long_reac"].keys()), key=plain_natural_key)
            for rnm in woltr:
                reac=netan["sto_r_m"][rnm]
                for lr in ("left", "right"):
                    dloc[lr]=" + ".join(f"{str(tmet[1])+'*' if tmet[1] != 1 else ''}{tmet[0]}" for i,tmet in enumerate(reac[lr]))
                s=f"{rnm}:\t{dloc['left']} {'<' if rnm not in netan['notrev'] else ''}-> {dloc['right']}\n"
                res += s
        elif suff == ".linp":
            # label input
            header="Id\tComment\tSpecie\tIsotopomer\tValue\n"
            for d in ftbl["LABEL_INPUT"]:
                met=d["META_NAME"] if d["META_NAME"] else met
                res += f"\t\t{met}\t{d['ISOTOPOMER'][1:]}\t{d['VALUE']}\n"
        elif suff == ".miso":
            # label, peak, ms measurements
            header="Id\tComment\tSpecie\tFragment\tDataset\tIsospecies\tValue\tSD\tTime\n"
            if case_i:
                # pick file_labcin from opt
                #import pdb; pdb.set_trace()
                flabcin=[d["OPT_VALUE"] for d in ftbl.get("OPTIONS", []) if d["OPT_NAME"] == "file_labcin"]
                if not flabcin:
                    warn("option '--inst' was activated but a field 'file_labcin' was not found in OPTIONS in '%s'. Only simulations will be possible (not fitting)."%fftbl.name)
                else:
                    if len(flabcin) > 1:
                        werr("option '--inst' was activated but a field 'file_labcin' is not unique in OPTIONS of '%s'"%fftbl.name)
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
            if case_i and flabcin:
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
                #import pdb; pdb.set_trace()
                labset=0
                for d in ftbl.get("LABEL_MEASUREMENTS", {}):
                    met=d["META_NAME"] if d["META_NAME"] else met
                    if d["META_NAME"]:
                        labset += 1
                        cgr=d["CUM_GROUP"]
                    elif cgr != d["CUM_GROUP"]:
                        labset += 1
                    res += f"\t\t{met}\t\tLAB-{labset}\t{d['CUM_CONSTRAINTS'].replace('#', '')}\t{d['VALUE']}\t{d['DEVIATION']}\t\n"
            # peak
            # ps, pdm, pdp = peak singlet, d-, d+
            if case_i and flabcin:
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
                    elif atom == netan["Clen"][met]:
                        frag="%s,%s"%(pdm, ps)
                    else:
                        frag="%s,%s,%s"%(pdm,ps,pdp)
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
                for d in ftbl.get("PEAK_MEASUREMENTS", {}):
                    met=d["META_NAME"] if d["META_NAME"] else met
                    pset += 1
                    atom=int(d["PEAK_NO"])
                    ps=atom
                    pdm=str(ps-1)
                    pdp=str(ps+1)
                    ps=str(ps)
                    if ps == "1":
                        frag="1,2"
                    elif atom == netan["Clen"][met]:
                        frag="%s,%s"%(pdm, ps)
                    else:
                        frag="%s,%s,%s"%(pdm,ps,pdp)
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
            if case_i and flabcin:
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
                        im=[ii for ii,v in enumerate(measures["mass"][0]["ids"]) if v.startswith(f"m:{met}:{frag}:{w}:")]
                        if len(im) != 1:
                            continue # not found corresponding unique FTBL entry, silently skip it
                        im=im[0]
                    sdv=measures["mass"][0]["dev"][im]
                    for ti in df_cin.columns[1:]:
                        res += f"\t\t{met}\t{frag.replace('~', '-')}\tMS-{mset}\tM{w}\t{df_cin.loc[i, ti]}\t{sdv}\t{ti}\n"
            else:
                mset=0
                for d in ftbl.get("MASS_SPECTROMETRY", {}):
                    met=d["META_NAME"] if d["META_NAME"] else met
                    frag=d["FRAGMENT"] if mset == 0 or d["FRAGMENT"] else frag
                    if d["META_NAME"] or d["FRAGMENT"]:
                        mset += 1
                    res += f"\t\t{met}\t{frag.replace('~', '-')}\tMS-{mset}\tM{d['WEIGHT']}\t{d['VALUE']}\t{d['DEVIATION']}\t\n"
        elif suff == ".mflux":
            # flux measurements
            header="Id\tComment\tFlux\tValue\tSD\n"
            for d in ftbl["FLUX_MEASUREMENTS"]:
                res += f"\t\t{d['FLUX_NAME']}\t{d['VALUE']}\t{d['DEVIATION']}\n"
        elif suff == ".mmet":
            # metabolite measurements
            header="Id\tComment\tSpecie\tValue\tSD\n"
            for d in ftbl.get("METAB_MEASUREMENTS", []):
                res += f"\t\t{d['META_NAME']}\t{d['VALUE']}\t{d['DEVIATION']}\n"
        elif suff == ".tvar":
            # type of variable and starting values
            header="Id\tComment\tName\tKind\tType\tValue\n"
            for nx in ("NET", "XCH"):
                for d in ftbl["FLUXES"][nx]:
                    #if d['NAME'] == "BM":
                    #    import pdb; pdb.set_trace()
                    if nx == "XCH" and ((d["FCD"] == "C" and float(d['VALUE(F/C)']) == 0. and d['NAME'] in netan["sto_r_m"]) or (d["FCD"] == "D" and d["NAME"] in netan["flux_inout"])):
                        continue
                    res += f"\t\t{d['NAME']}\t{nx}\t{d['FCD']}\t{d['VALUE(F/C)']}\n"
            for d in ftbl.get("METABOLITE_POOLS", []):
                try:
                    val=float(d['META_SIZE'])
                except:
                    val=float(eval(d['META_SIZE']))
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
            print(str(cout))
            cout.write_text(f"{scre} at {dtstamp()}\n"+header+res)

def main(argv=sys.argv[1:]):
    # parse options
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inst", help="activate instationary mode", action="store_true")
    parser.add_argument("-f", "--force", help="force overwriting of result files", action="store_true")
    parser.add_argument("-o", "--out", help="path prefix for result files", nargs=1)
    parser.add_argument("ftbl", help="input file to be converted to MTF", nargs=1)
    opts = parser.parse_args(argv)
    #print("opts=", opts)
    force=opts.force
    case_i=opts.inst
    out=opts.out[0] if opts.out else ""
    f=opts.ftbl[0]
    fftbl=Path(f) if f != "-" else sys.stdin
    
    if type (fftbl) == type(Path()) and not fftbl.is_file():
        fftbl = fftbl.with_suffix(".ftbl")
        if not fftbl.is_file():
            werr("file '"+str(fftbl)+"' does not exist.\n")
    # prepare out
    if out:
        out=Path(out)
        if out.is_dir() and type (fftbl) == type(Path()):
            out=out/fftbl.stem
    else:
        if type (fftbl) == type(Path()):
            out=fftbl.parent/fftbl.stem
        else:
            out=Path("mtf")
    #print("out=", out)
    #sys.exit(1)
    ftbl=C13_ftbl.ftbl_parse(str(fftbl))
    if not case_i:
        # search for 'file_labcin' or 'funlabR' in OPTIONS to switch to case_i=True if any
        f=[d["OPT_VALUE"] for d in ftbl.get("OPTIONS", []) if d["OPT_NAME"] in ("file_labcin", "funlabR") and d["OPT_VALUE"]]
        if f:
            warn("Switching to instationary mode as non empty 'file_labcin' is found in OPTIONS")
            case_i=True
    #print("ftbl parsed=", ftbl)
    #import pdb; pdb.set_trace()
    netan=dict()
    C13_ftbl.ftbl_netan(ftbl, netan, case_i=case_i)
    if (not case_i) and ("OPTIONS" in ftbl and netan["opt"].get("file_labcin", [""])[0]):
        warn("we are in stationary case but ftbl file has 'file_labcin' option")
    bsl="\\" # backslash
    scre=f"# Created by '{me} {' '.join(v.replace(' ', bsl+' ') for v in argv)}'"
    ftbl2suff(ftbl, fftbl, case_i, netan, force, out, scre, (".netw", ".linp", ".miso", ".mflux", ".mmet", ".tvar", ".cnstr", ".opt"))
    #import pdb; pdb.set_trace()
    if "OPTIONS" in ftbl and "prl_exp" in netan["opt"]:
        for pftbl in netan["opt"]["prl_exp"].split(";"):
            pftbl=pftbl.strip()
            if not pftbl:
                continue
            fpftbl=try_ext(fftbl.parent/pftbl, ["ftbl"])
            dftbl=C13_ftbl.ftbl_parse(str(fpftbl))
            ftbl2suff(dftbl, fpftbl, case_i, netan, force, out.parent/pftbl, scre, (".linp", ".miso", ".opt"))
    return 0

if __name__ == "__main__" or __name__ == "influx_si.cli":
    main()
