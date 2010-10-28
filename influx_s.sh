#!/bin/sh
# wrapper script for minimizing static fluxes by R script
# - generate .f and .R files by ftbl2optR.py
# - compile .f in .so
# - start generated R script for minimization
#
# usage:
# ./influx_s.sh [...] network[.ftbl]
# optional extra params [...] are passed as is to R script
par="$@"
n=$#
if [ $n = 0 ]; then
   echo "usage: ./influx_s.sh [opt_params_to_R_script] network[.ftbl]";
   exit 1;
fi

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
   eargs="${@:1:n}";
fi
echo "$0 $par" > "$f.log";
echo "code gen:" $(date) | tee -a "$f.log";

# check if static fortran functions are ready
[ $direx/cumo.f -nt $direx/cumo.so ] && { R CMD SHLIB --clean "$direx/cumo.f" || exit 1 ;}

$direx/ftbl2optR.py "$f" 2> "$f.err" &&
#   echo "compil  :" $(date) | tee -a "$f.log" &&
#   R CMD SHLIB --clean "$f.f" 1> /dev/null 2>> /dev/null &&
   echo "calcul  :" $(date) | tee -a "$f.log" &&
   R --no-save --slave --args $eargs < "$f.R" \
      2>> "$f.err" 1>> "$f.log";
echo "end     :" $(date) | tee -a "$f.log";

[ -s "$f.err" ] && echo "=>Check $f.err" | tee -a "$f.log"
