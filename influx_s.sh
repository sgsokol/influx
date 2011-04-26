#!/bin/sh
# wrapper script for minimizing static fluxes by R script
# - generate .f and .R files by ftbl2optR.py
# - compile .f in .so
# - start generated R script for minimization
#
# usage:
# ./influx_s.sh [...] network[.ftbl]
# optional extra params [...] are passed as is to R script
function usage {
   cat <<EOS
usage: ./influx_s.sh [options] network[.ftbl]
Options:
 -h|--help
 --noopt
    no optimization, just use free fluxes as is.
 --noscale
    no scaling factors to optimize => all scaling factors are 1.
 --meth BFGS|Nelder-Mead|nlsic
    nlsic by default
 --fullsys
    calculate all cumomer set (not just reduced necesary to simulate measurements)
 --irand
    ignore initial approximation from FTBL file and use random values instead
 --sens mc[=N]
    sensitivity method: mc for Monte-Carlo. N is number of Monte-Carlo simulations
    Default is 10
 --cupx N
    upper limit for reverse fluxes. Must be in interval [0; 1].
    Default is 0.99999
 --cupn N
    upper limit for net fluxes.
    Default is 1.e5
 --clownr N
    lower limit for not reversible free and dependent fluxes
    Default is 0, i.e. no lower limit
 --np N
    Number of parallel process used in Monte-Carlo simulations.
    Without this option all available cores in a given node are used.
 --ln
    Least norm solution is proposed when Jacobian is rank deficient.
 --zc
    Apply zero crossing strategy for net fluxes.
For developers:
 --DEBUG
 --prof
EOS
}

par="$@"
n=$#
if [ $n -eq 0 -o \( $n -eq 1 -a \( "$1" = "-h" -o "$1" = "--help" \) \) ]; then
   usage;
   echo hope it\'s helping
   exit 0;
fi

direx=$(dirname "$0");

# test which argument is ftbl file, first or last or none
f="${1%%.ftbl}"
if [ -r "$f".ftbl ]; then
   shift;
   eargs="$@";
else
   f="${@:n:1}";
   f="${f%%.ftbl}"
   if [ -r "$f".ftbl ]; then
      f="${f%%.ftbl}";
      n=$((n-1));
      eargs="${@:1:n}";
   else
      eargs="$@";
   fi
fi
# test if there are options to be passed to python and store them in setpy array
pyopt=("-h" "--help" "--fullsys" "--DEBUG")
setpy=()
ipy=0
for o in ${pyopt[@]}; do
   #echo o=$o
   [ "${eargs[*]/$o/}" = "${eargs[*]}" ] || { setpy[ipy]="$o"; ipy=$((ipy+1)) ;}
done
#echo setpy="${setpy[@]}"

if [ ! -r "$f".ftbl ]; then
   usage 1>&2;
   echo 1>&2;
   echo "error: FTBL file '$f.ftbl' does not exist or is not readable." 1>&2
   exit 1;
fi

echo "$0 $par" > "$f.log";
echo "code gen:" $(date) | tee -a "$f.log";

# check if static fortran functions are ready
[ $direx/cumo.f -nt $direx/cumo.so ] && { R CMD SHLIB --clean "$direx/cumo.f" || exit 1 ;}

$direx/ftbl2optR.py ${setpy[@]} "$f" 2> "$f.err" &&
#   echo "compil  :" $(date) | tee -a "$f.log" &&
#   R CMD SHLIB --clean "$f.f" 1> /dev/null 2>> /dev/null &&
   echo "calcul  :" $(date) | tee -a "$f.log" &&
   R --no-save --slave --args $eargs < "$f.R" \
      2>> "$f.err" 1>> "$f.log";
echo "end     :" $(date) | tee -a "$f.log";

[ -s "$f.err" ] && echo "=>Check $f.err" | tee -a "$f.log"
