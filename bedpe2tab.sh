
cat large_sv_calls.bedpe |awk  '{if ($1==$4) print "call", $1, $2, $6, $8, $12}' >tmp
#cat ./large_sv_candidates.bedpe | awk  '{if ($1==$4) print "candidate", $1, $2, $6, $8, $12}' >>tmp

perl process_svtype.pl tmp >svs.tab
 

