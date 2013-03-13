# parse _res.kvh file for extracting free parameters and setting them
# into corresponding ftbl file. The result is written on stdout.
# usage: ffres2ftbl.sh network_res.kvh [base.ftbl] > new.ftbl

# The modified ftbl file name is optional. If omited it is
# extracted from kvh name, here for example it would be network.ftbl

[ $# -eq 0 -o $# -gt 2 ] && { echo 1>&2 "ffre2ftbl: Wrong parameter number\
\nusage: ffres2ftbl.sh network_res.kvh [base.ftbl] > new.ftbl"; }
fkvh="$1"
# strip "_res.kvh"
if [ $# -eq 2 ]; then
   fftbl="$2"
else
   fftbl="${fkvh%_res.kvh}".ftbl
fi

# find executable dir
direx=$(which "$0")
direx=$(dirname "$direx")
ftmp="tmp_${fftbl%.ftbl}".$$.kvh
kvh_get.py "linear stats" "net-xch01 fluxes (sorted by name)" < "$fkvh" |
   cut -f 1,2 | fgrep -v row_col > "$ftmp"
kvh_get.py "linear stats" "metabolite pools (sorted by name)" < "$fkvh" |
   cut -f 1,2 | fgrep -v row_col >> "$ftmp"
"$direx"/ff2ftbl.py "$ftmp" "$fftbl" && rm -f "$ftmp"
