#!/bin/sh
# wrapper script for minimizing static fluxes by R script
# - generate .f and .R files by ftbl2optR.py
# - compile .f in .so
# - start generated R script for minimization
#
# usage: ./opt.sh network[.ftbl] [...]
# optional extra params [...] are passed as is to R script
if [ $# = 0 ]; then
   echo "usage: ./opt.sh network[.ftbl] [opt params to R script]";
   exit 1;
fi

echo "compil:" $(date);
direx=/home/sokol/insa/sysbio/dev/ftbl2sys;
f="$1"
shift;
eargs="$@"

$direx/ftbl2optR.py "$f" &&
   R CMD SHLIB "$f.f" &&
   echo "calcul:" $(date) &&
   R --no-save --silent $eargs < "$f.R" \
   > "$f.log" 2> "$f.err";
echo "end:   " $(date);
