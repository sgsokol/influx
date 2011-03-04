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

$direx/influx_i.py "$f" $DEBUG 2>"$f.err" | tee -a "$f.log" &&
   echo "calcul  :" $(date) && R --no-save --silent --args $eargs < "$f.R" \
   1>> "$f.log" 2>> "$f.err";
echo "end     :" $(date) | tee -a "$f.log";

[ -s "$f.err" ] && echo "=>Check $f.err" | tee -a "$f.log"
