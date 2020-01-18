#!/bin/sh
touch edge_lineage.txt
NAME=${1?Error: no name given}
cut -d "," -f15 $NAME | grep -v "edge_name"| awk -F_ '{print $1}'|awk -F\| '{print $2}' |while read headername

do
line1="grep '$headername' dmspDB_key.txt  |cut -d " " -f2-3 >> edge_lineage.txt"
eval $line1

done

echo "key generic" | cat - edge_lineage.txt > edge_lineages.txt
rm edge_lineage.txt
