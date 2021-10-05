#!/bin/sh
touch tmp.txt
NAME=${1?Error: no name given}
cut -d "," -f15 $NAME | grep -v "edge_name"| awk -F_ '{print $1}'|awk -F\| '{print $2}' |while read headername

do
line1="grep '$headername' dmspdb_key.csv  |cut -d "," -f2-3 >> tmp.txt"
eval $line1
done

echo "key,generic" | cat - tmp.txt > edge_lineages.txt
rm tmp.txt
