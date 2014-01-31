#!/usr/bin/env python
"""Parse _res.kvh file from stdin or from first parameter
and write the measurment section of ftbl file on stdout
usage: res2ftbl_meas.py [network_res[.kvh]] > meas.ftbl
"""
if __name__ == "__main__":
    import sys, os, getopt
    import tools_ssg
    import kvh
    import pprint
    def usage():
        print(__doc__)
    me="res2ftbl_meas"
    secdef={
        "m": """
MASS_SPECTROMETRY
	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION""",
        "p": """
PEAK_MEASUREMENTS
	META_NAME	PEAK_NO	VALUE_S	VALUE_D-	VALUE_D+	VALUE_DD	VALUE_T	DEVIATION_S	DEVIATION_D-	DEVIATION_D+	DEVIATION_DD/T""",
    }
    pfield="META_NAME	PEAK_NO	VALUE_S	VALUE_D-	VALUE_D+	VALUE_DD	VALUE_T	DEVIATION_S	DEVIATION_D-	DEVIATION_D+	DEVIATION_DD/T".split("\t")
    sdfield={
        "S": "DEVIATION_S",
        "D-": "DEVIATION_D-",
        "D+": "DEVIATION_D+",
        "DD": "DEVIATION_DD/T",
        "T": "DEVIATION_DD/T",
    }
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "DEBUG"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)
    DEBUG=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--DEBUG":
            DEBUG=True
        else:
            assert False, "unhandled option"
    fkvh=args[0] if len(args) else ""
    if fkvh and fkvh[-4:] != ".kvh":
        fkvh+=".ftbl"
    if fkvh and not os.path.exists(fkvh):
        sys.stderr.write(me+": file '"+fkvh+"' does not exist.\n")
        sys.exit(1)
    fkvh=open(fkvh, "r") if fkvh else sys.stdin

    d=kvh.kvh2tlist(fkvh)
    #pprint.pprint(d)
    
    # flux measurements
    meas=[]
    try:
        meas=[ v for (k,v) in d if k == "simulated flux measurements" ].pop()
    except:
        pass
    if meas:
        print("""FLUX_MEASUREMENTS
	FLUX_NAME	VALUE	DEVIATION""")
        for (k,v) in meas[1:]:
            print("\t%s\t%s"%(k[4:], v))
    # labeling measurments
    meas=dict( t for t in d if t[0] in ("simulated scaled labeling measurements", "simulated unscaled labeling measurements") )
    meas=meas.get("simulated scaled labeling measurements", meas["simulated unscaled labeling measurements"])
    #pprint.pprint(meas)
    
    mtype=""
    dpeak={}
    for m in meas[1:]:
        # split id
        sid=m[0].split(":")
        if mtype != sid[0]:
            # new measurement type starts here
            mtype=sid[0]
            print(secdef[mtype])
            oldmetab=""
            oldfrag=""
            oldweight=-1
            oldrid=""
            oldc_no=0
            # flush previous peak dict
            if dpeak:
                print("\t".join([""] + [dpeak.get(f, "") for f in pfield]))
                dpeak={}
        if mtype=="m":
            ## mass measurement
            metab, frag, weight, rid = sid[1:]
            if metab != oldmetab or frag != oldfrag or int(weight) <= oldweight:
                # new metab starts here
                oldmetab, oldfrag, oldweight = (metab, frag, int(weight))
                print("\t%s\t%s\t%d\t%s"%(metab, frag, int(weight), m[1]))
            else:
                print("\t%s\t%s\t%d\t%s"%("", "", int(weight), m[1]))
        elif mtype=="p":
            ## peak measurement
            metab, c_no, ptype, rid = sid[1:]
            val,sd=m[1].split("\t")
            if oldmetab != metab:
                ## new metab starts here
                oldmetab=metab
                oldrid=""
                # flush previous dict
                if dpeak:
                    print("\t".join([""] + [dpeak.get(f, "") for f in pfield]))
                    dpeak={}
                dpeak["META_NAME"]=metab
            if oldrid != rid:
                # new peak for a given metab starts here
                oldrid=rid
                # flush previous dict
                if "PEAK_NO" in dpeak:
                    print("\t".join([""] + [dpeak.get(f, "") for f in pfield]))
                    dpeak={}
                dpeak["PEAK_NO"]=c_no
            dpeak["VALUE_"+ptype]=val
            dpeak[sdfield[ptype]]=sd
    # flush last peak dict
    if dpeak:
        print("\t".join([""] + [dpeak.get(f, "") for f in pfield]))
        dpeak={}
    sys.exit(0)
