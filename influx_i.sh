echo "compil:" $(date);
direx=$(dirname $0);
me=$(basename $0)
DEBUG=""
[ "$2" = "DEBUG" ] && DEBUG="DEBUG"

org="$1";
[ -n "$org" ] || { echo 2>1 "$me: expecting ftbl file name"; exit 1; }
org="${org%.ftbl}"; # strip out .ftbl if any
$direx/influx_i.py "$org" $DEBUG 2>"$org.err" &&
   R CMD SHLIB "$org.f" &&
   echo "calcul:" $(date) && R --no-save --silent --args --meth $DEBUG < "$org.R" \
   > "$org.log" 2>> "$org.err";
echo "end   :" $(date);
