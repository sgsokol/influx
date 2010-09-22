#!/bin/sh
# wrapper script for minimizing static fluxes by R script
# - generate .f and .R files by ftbl2optR.py
# - compile .f in .so
# - start generated R script for minimization
#
# usage: ./influx_s.sh network[.ftbl] [...]
# optional extra params [...] are passed as is to R script
if [ $# = 0 ]; then
   echo "usage: ./influx_s.sh network[.ftbl] [opt params to R script]";
   exit 1;
fi

echo "code gen:" $(date);
direx=$(dirname "$0");
f="$1"
shift;
eargs="$@"

$direx/ftbl2optR.py "$f" &&
   echo "compil  :" $(date) &&
   R CMD SHLIB "$f.f" &&
   echo "calcul  :" $(date) &&
   R --no-save --silent $eargs < "$f.R" \
   > "$f.log" 2> "$f.err";
echo "end     :" $(date);
