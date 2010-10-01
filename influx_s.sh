#!/bin/sh
# wrapper script for minimizing static fluxes by R script
# - generate .f and .R files by ftbl2optR.py
# - compile .f in .so
# - start generated R script for minimization
#
# usage:
# ./influx_s.sh [...] network[.ftbl]
# optional extra params [...] are passed as is to R script
n=$#
if [ $n = 0 ]; then
   echo "usage: ./influx_s.sh [opt_params_to_R_script] network[.ftbl]";
   exit 1;
fi

echo "code gen:" $(date);
direx=$(dirname "$0");

# test which argument is ftbl file, first or last
f="${1%%.ftbl}"
if [ -r "$f".ftbl ]; then
   shift;
   eargs="$@";
else
   f="${@:n:1}";
   f="${f%%.ftbl}";
   n=$((n-1));
   eargs="$@:1:n";
fi

$direx/ftbl2optR.py "$f" &&
   echo "compil  :" $(date) &&
   R CMD SHLIB "$f.f" &&
   echo "calcul  :" $(date) &&
   R --no-save --silent $eargs < "$f.R" \
   > "$f.log" 2> "$f.err";
echo "end     :" $(date);
