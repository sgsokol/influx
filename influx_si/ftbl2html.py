#!/usr/bin/env python3

r"""
read a .ftbl file from a parameter (or from --prefix/--mtf options) and translate to .html file with mermaid.
Reactions involving more than one substrate and/or more than one product are represented
by an additional node while one-to-one reactions are just edges.
Node and edge attributes are written in respective css classes.
Compatibility: mermaid@11

usage: ftbl2html.py [-h|--help] mynetwork.ftbl|--prefix PREFIX|--mtf MTF [> mynetwork.html]

OPTIONS
-h, --help print this message and exit

:param: mynetwork the base of an ftbl file (mynetwork.ftbl)

:returns: mynetwork.html -- file of the network diagram

Copyright 2025, INRAE/INSA/CNRS, France
Author: Serguei Sokol (sokol at insa-toulouse dot fr)
License: Gnu Public License (GPL) v3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import os
import stat
import getopt
import re

import influx_si
from tools_ssg import *
from C13_ftbl import *
import txt2ftbl

def color(m, netan, nc):
    return nc['i'] if m in netan['input'] \
        else nc['o'] if m in netan['output'] else \
        nc['m']
def usage():
    sys.stderr.write(__doc__)
def htmllab(lab):
    """Format metabolite/reaction label as HTML text 'A_c' -> 'A<sub>c</sub>'"""
    li=str(lab).split("_")
    if len(li) > 1:
        # do formatting
        return "_".join(li[:-1])+"<sub>%s</sub>"%li[-1]
    else:
        return lab # as is
# Configurable constants
ns={
    'm': 'rounded',       # metabolite shape
    'r': 'hex'        # reaction shape
}
nc={
    'm': 'ffffff',     # metabolite colour
    'i': 'BED3C3',      # input colour (uptake)
    'o': 'FB7C72',      # output colour (escape)
    'd': 'C7A070',      # dead-end colour
    'r': '7F7F7F'       # reaction colour
}
html_header=r"""
<!doctype html>
<html lang="en">
  <title>%s</title>
  <body>
    <form>
    <button type="submit" id="submit"">Save as SVG image</button>
    </form>
    <pre class="mermaid" id="graph">
flowchart TB
    subgraph %s
    direction TB
"""
# get arguments
me=os.path.basename(__file__)

def main(argv=sys.argv[1:]):
    werr=sys.stderr.write

    # determine colour of metabolite (in, out or plain)

    try:
        opts,args=getopt.getopt(sys.argv[1:], "hfi", ["help", "force", "prefix=", "mtf=", "inst"])
    except getopt.GetoptError as err:
        sys.stderr.write(str(err)+"\n")
        usage()
        sys.exit(1)
    cost=False
    force=False
    case_i=False
    li_ftbl=[]
    mtf_opts=[]
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--force"):
            force=True
        if o in ("-i", "--inst"):
            case_i=True
            mtf_opts += ["--inst"]
        if o in ("--mtf", "--prefix"):
            mtf_opts += [o, a]
    if "--mtf" in mtf_opts or "--prefix" in mtf_opts:
        txt2ftbl.main(mtf_opts, li_ftbl)
    if not args and not li_ftbl:
        sys.stderr.write("Error: expecting ftbl file name or --prefix/--mtf options\n")
        usage()
        sys.exit(1)
    base=args[0] if args else li_ftbl[0]
    if base[-5:]==".ftbl":
        base=base[:-5]
    path_ftbl=base+".ftbl"
    
    # what kind of output we have?
    mode=os.fstat(1).st_mode
    fout=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(path_ftbl[:-4]+"html", "w")
    breakpoint
    # define where to read and write
    # input file is an argument
    fdir=os.path.dirname(base) or "."
    base=os.path.basename(base)
    short_ftbl=base+".ftbl"

    # Parse .ftbl file
    try:
        ftbl=ftbl_parse(path_ftbl)
    except Exception as inst:
        werr(str(inst)+"\n")
        if not force:
            raise
    # Analyse the network
    netan=dict()
    try:
        ftbl_netan(ftbl, netan, case_i=case_i)
    except:
        if force:
            werr(str(sys.exc_info()[1])+"\n")
        else:
            raise
    # transform metabs,reac in nodes (dict of properties) and reac in edges (dict too)
    mlen=len(netan["metabs"])
    
    # metabs -> nodes
    nodes={"metabs": {}, "reacs": {}}
    nodes["metabs"].update((metab,
        {"label": metab, "id": i+1, "shape": ns["m"],
        "color": nc["i"] if metab in netan["input"] else nc["o"] if metab in netan["output"] else nc["d"] if metab in netan["deadend"] else nc["m"]})
        for (i, metab) in enumerate(netan["metabs"]))
    # reacs -> nodes
    nodes["reacs"].update((reac,
        {"label": reac, "id": i+1+mlen, "shape": ns["r"],
        "color": nc["r"]})
        for (i, (reac, d)) in enumerate((reac, d) for (reac, d) in netan["sto_r_m"].items() if len(d["left"]) > 1 or len(d["right"]) > 1))
    # easy id finders
    metab2id=dict((m, d["id"]) for (m, d) in nodes["metabs"].items())
    reac2id=dict((r, d["id"]) for (r, d) in nodes["reacs"].items())
    id2dic=dict((d["id"], d) for mr in ("metabs", "reacs") for (_, d) in nodes[mr].items())
    
    # write header
    fout.write(html_header%(base, base))
    
    # write reactions
    fout.write("\n")
    for (r, d) in netan["sto_r_m"].items():
        # reversible or not?
        rnr="nr" if r in netan["notrev"] else "r"
        # forward flux is dependent, free or constrained?
        fw_dfcg=netan["nx2dfcg"].get("n."+r, "u")[0]
        # revers flux is dependent, free or constrained?
        rv_dfcg=netan["nx2dfcg"].get("x."+r, "u")[0]
        eds=[] # list of dicts: name, ids of source and target, reac, source and target arrow shape, fw&rv color
        if r in reac2id:
            # complex reaction
            rid=reac2id[r]
            subs=d["left"] # substrates of this reaction
            prods=d["right"] # products of this reaction
            eds.extend({
               "label_s": htmllab(m),
               "label_t": "&#8194;"*6+htmllab(r),
               "class_s": "metab"+("_i" if m in netan["input"] else "_o" if m in netan["output"] else "_d" if m in netan["deadend"] else ""),
               "class_t": "reac",
               "link": "==>" if rnr == "r" else "-->",
               "link_text": "",
               "id_s": metab2id[m],
               "id_t": rid,
               "reac": r,
               } for (isu, (m, c)) in enumerate(subs)
            )
            eds.extend({
               "label_s": "&#8194;"*6+htmllab(r),
               "label_t": htmllab(m),
               "class_s": "reac",
               "class_t": "metab"+("_i" if m in netan["input"] else "_o" if m in netan["output"] else "_d" if m in netan["deadend"] else ""),
               "link": "==>" if rnr == "r" else "---",
               "link_text": "",
               "id_s": rid,
               "id_t": metab2id[m],
               "reac": r,
               } for (ipr, (m, c)) in enumerate(prods)
            )
        else:
            # simple reaction
            s=d["left"][0][0]
            t=d["right"][0][0]
            eds.append({
               "label_s": htmllab(s),
               "label_t": htmllab(t),
               "class_s": "metab"+("_i" if s in netan["input"] else "_o" if s in netan["output"] else "_d" if s in netan["deadend"] else ""),
               "class_t": "metab"+("_i" if t in netan["input"] else "_o" if t in netan["output"] else "_d" if t in netan["deadend"] else ""),
               "link": "<==>" if rnr == "r" else "-->",
               "link_text": "|"+htmllab(r)+"|",
               "id_s": metab2id[s],
               "id_t": metab2id[t],
               "reac": r,
            })
        for st in eds:
            # id1:::metab_i@{label: "A", shape: rounded} -->|reac_name| id2:::metab@{label: "B", shape: rounded}
            shape_s="" if st["class_s"] == "reac" else ", shape: "+ns["m"]
            shape_t="" if st["class_t"] == "reac" else ", shape: "+ns["m"]
            A=f"""id{st["id_s"]}:::{st["class_s"]}@{{label: "{st["label_s"]}"{shape_s}}}"""
            B=f"""id{st["id_t"]}:::{st["class_t"]}@{{label: "{st["label_t"]}"{shape_t}}}"""
            if st["class_t"] == "reac":
                # workaround for reverse string
                # instead of A <-- B we have to write 3 edges:
                # A ~~~ B --> A; B ~~~ A
                s=f"""    {A} ~~~ {B} {st["link"]}{st["link_text"]} {A}
    {B} ~~~ {A}
"""
            else:
                s=f"""    {A} {st["link"]}{st["link_text"]} {B}
"""
            fout.write(s)
            #print("s=%s"%s)
    # print footer
    footer=f"""
    end
    subgraph Legend
        direction LR
        id1000000:::metab@{{label: "Internal metabolite", shape: %s}} ~~~ id1000010:::empty@{{label: " "}}
        id1000001:::metab_i@{{label: "Input metabolite", shape: %s}} ~~~ id1000011:::empty@{{label: " "}}
        id1000002:::metab_o@{{label: "Output metabolite", shape: %s}} ~~~ id1000012:::empty@{{label: " "}}
        id1000003:::metab_d@{{label: "Dead-end metabolite", shape: %s}} ~~~ id1000013:::empty@{{label: " "}}
        id1000004:::empty@{{label: " "}} -->|Non reversible reaction| id1000005:::empty@{{label: " "}}
        id1000006:::empty@{{label: " "}} <==>|Reversible reaction| id1000007:::empty@{{label: " "}}
    end
    {base} ~~~ Legend
classDef metab fill:#{nc["m"]}
classDef metab_i fill:#{nc["i"]}
classDef metab_o fill:#{nc["o"]}
classDef metab_d fill:#{nc["d"]}
classDef reac height:0, width: 0
classDef empty height:0, width: 0
    </pre>
    <script type="module">
      import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.esm.min.mjs';
      const saveFile = async function () {{
        const svg=document.getElementById("graph").innerHTML;
        var myBlob = new Blob([svg+"\\n"], {{type: "image/svg+xml"}});
        var anchor = document.createElement("a");
        anchor.href = window.URL.createObjectURL(myBlob);
        anchor.download = "{base}.svg";
        anchor.click();
      }};
      document.querySelector('#submit').addEventListener('click', saveFile)
    </script>
  </body>
</html>
"""%(ns["m"], ns["m"], ns["m"], ns["m"])
    fout.write(footer)
    fout.close()

if __name__ == "__main__":
    main()
