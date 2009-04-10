echo "compil:" $(date);
direx=$(dirname $0);
me=$(basename $0)
DEBUG=""
[ "$2" = "DEBUG" ] && DEBUG="DEBUG"

org="$1";
org="${org%.ftbl}"; # strip out .ftbl if any
$direx/ff2strid.py "$org" $DEBUG &&
   R CMD SHLIB "$org.f" &&
   echo "calcul:" $(date) && R --no-save --slave --args $DEBUG < "$org.R" \
   > "$org.log" 2> "$org.err";
echo "end   :" $(date);
