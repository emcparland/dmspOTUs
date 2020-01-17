#!/bin/sh
NAME=${1?Error: no name given}
cut -d "," -f15 $NAME | grep -v "edge_name"| awk -F_ '{print $1}'|while read headername

do
line1="grep '$headername' mmetsp_all_lineages.csv  |cut -d "," -f6-13 >> edge_lineage.txt"
eval $line1

done

echo "genus,family,order,class,norank,norank,kingdom,generic" | cat - edge_lineage.txt > edge_lineages.txt
rm edge_lineage.txt
