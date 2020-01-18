#!/bin/sh
NAME=${1?Error: no name given}
cut -d "," -f4 $NAME | grep -v "edge_num"| while read headername

do
line1="grep '|$headername ' edge_code.txt |awk '{print $"1"}' >> edge_tr.txt"
eval $line1

done

# add a header
sed -i '1iedge_name' edge_tr.txt
# remove the weird numbers added during the alignment
cut -d "," -f15 edge_tr.txt | awk -F\\ '{print $1}' >edge_tr.c.txt
paste -d "," query.align.csv edge_tr.c.txt > query.align.wkey.csv
