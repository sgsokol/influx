#!/bin/sh
function usage {
echo 1>&2 "
usage:" $me "[-h|--help|--cost|--DEBUG] param_file.kvh

OPTIONS
-h, --help print this message and exit
--cost request to calculate cost function value at feasible points
--DEBUG enables some debugging features and output (for advanced users)

REMARKS
Base name of mesh parameter file ('param_file' in this example)
is used to create or silently overwrite the following files:
param_file.R
param_file.f
param_file.o
param_file.so
param_file.log
param_file.err
param_file_res.kvh
The last file get the results in kvh format:
 - feasibles free flux sets with their values and
 - if requested, cost values at feasible points;
";
}
me=$(basename $0)
#echo me=$me $#;
SHORTOPTS="h"
LONGOPTS="help,cost,DEBUG"
ORIGOPTS="$@"

ARGS=$(getopt -s bash --options $SHORTOPTS  \
  --longoptions $LONGOPTS --name $me -- "$@" )
[ "$?" = "0" ] || { echo "use '$me --help' for more information"; exit 1 ;}
#echo args="$ARGS"
#echo "rargs=$?"

eval set -- "$ARGS"

while true; do
   case $1 in
      -h|--help)
         usage;
         exit 0;
         ;;
      --cost|--DEBUG)
         ;;
      --)
         shift
         break
         ;;
      *)
         echo 1>&2 $# "Unhandled option '$1'";
         exit 1;
         ;;
   esac
   shift
done
direx=$(dirname $0);
DEBUG="${DEBUG:-}"
[ "$me" = "ffmeshd.sh" ] && DEBUG="--DEBUG"
#echo "nargs=$#"
[ $# -eq 1 ] || { echo 1>&2 "Expecting a kvh file name"; usage; exit 1 ;}
fkvh="$1"
base="${fkvh%.kvh}"
date;
$direx/ffmesh.py $ORIGOPTS $DEBUG &&
   R CMD SHLIB "$base".f &&
   date && R --no-save --silent < "$base".R \
   > "$base".log 2> "$base".err;
date;
