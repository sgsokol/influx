#!/bin/sh
# clean files after runing influx_s[i] in the current (default) or a directory given in first and next parameters

if [ $# != 0 ]; then
   dirs="$@"
else
   dirs="."
fi

for d in $dirs; do
   ( cd "$d" &&
   for f in *.ftbl; do
      #echo fb="${f/.ftbl/}"
      fb="${f/.ftbl/}"
      rm -f "$fb"{.R,.err,.log,.pres.csv,.RData,_res.kvh,_ires.kvh,.netan,.sif,.xgmml,_netan.kvh,_rev.txt,_net.txt,_fwd.txt,.xml,.pdf} *.$fb.attrs {edge,node}.*.$fb*\(.V+\([0-9]\)\) "$fb".mc*.err parallel.R
   done )
done
