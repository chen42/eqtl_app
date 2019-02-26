#!/bin/bash

window=2000000
chr=$1
loc=$2


grep -P "$chr\t" ./eqtl_data.tab |awk  -v loc="$2" -v win="$window"  'function abs(value) {return (value<0?-value:value)} { if (abs(loc-$4) < win ) print $0}' | sort -k4  >$chr\_$loc\_win$window.tab






