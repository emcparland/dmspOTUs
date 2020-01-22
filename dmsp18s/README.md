# Pplacer of OTU sequences in 18S sequence phylogeny of known DMSP producers to assign putative DMSP synthesis gene carriers.

## 
If you'll be using my previously made phylogeny of known DMSP producers (18S sequences I previously retrieved via ncbi based on the Supplemental Table in McParland and Levine 2019), you can *skip to Step 3*. 

There are n=107 18S sequences of known DMSP producers.
```
grep "^>" dmspdb_all.fa | wc -l
```

## 1. Create alignment of reference sequences with Infernal using the 18S ssu cov model as a reference for alignment.
Name of your fasta file with reference sequences
```
NAME_FA=dmspdb_all.fa
```
Degap if necessary
```
seqmagick mogrifty --ungap $NAME_FA
```
Align
```
cmalign --dna -o align.sto --outformat Pfam eukarya-0p1.conv.cm $NAME_FA
```
Convert to fasta format
```
seqmagick convert align.sto align.fa
```
Deduplicate sequences
```
seqmagick mogrify --deduplicate-sequences align.fa
```
**I then take the alignment offline to manually curate and remove any large gaps in the alignment with Jalview or Geneious.**

## 2. Create a phylogenetic tree of the alignment, reference package from the alignment, tree and stats file with pplacer
- Name of alignment to easily sed any extract characters in alignment that will mess up next steps, and name for file outputs
```
NAMEFA=align.fa 
NAMET=dmsp
```
```
sed -i 's/|/_/g' $NAMEFA
sed -i 's/=/_/g' $NAMEFA
sed -i 's/:/_/g' $NAMEFA
sed -i 's/(/_/g' $NAMEFA
sed -i 's/)/_/g' $NAMEFA
#sed -i 's/-/_/g' $NAMEFA
sed -i 's/\//_/g' $NAMEFA
sed -i 's/\.//g' $NAMEFA
```
Per Bowman tutorial: Taxtastic can't read a RAxML stats file with confidence values, work around this by building tree without scores first than calculate scores separately and add back to your previously generated tree.
First build the tree without confidence scores
```
raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -s $NAMEFA -n $NAMET.tre -f d -p 12345
```
Root the tree with RAxML
```
raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -f I -t RAxML_bestTree.$NAMET.tre -n root_$NAMET.tre
```
Generate the confidence scores for tree separately
```
raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root_$NAMET.tre -n root_conf_$NAMET.tre -s $NAMEFA
```
Create reference package with taxit command
```
taxit create -l 18S_rRNA -P refpkg --aln-fasta $NAMEFA --tree-stats RAxML_info.$NAMET.tre  --tree-file RAxML_fastTreeSH_Support.root_conf_$NAMET.tre
```

## 3. Phylogenetic placement of 18S OTUs onto reference
If you're using my phylogeny, you need the following that I have provided in the repo:
- refpkg (created in steps above for pplacer)
- align.fa (reference alignment to concatenate to you query reads)
- tr_edges.sh (a shell script that translates edge numbers into a reference column so you can translate which OTUs were aligned to which reference sequence)

Align the query reads and reference alignment
First, clean names in OTU file
```
sed -i 's/ /_/g' otus.fa
```
If yours have any funky symbols use this handy command from the Bowman tutorial
```
tr "[ -%,;\(\):=\.\\\*[]\"\']" "_" < otus.fa
```
Concatenate the query reads and reference alignment
```
cat align.fa  otus.fa > query.fa
```
Remove gaps
```
seqmagick mogrify --ungap query.fa
```
Align
```
cmalign --dna -o query.sto --outformat Pfam eukarya-0p1.conv.cm query.fa
```
Convert to fasta
```
seqmagick convert query.sto query.align.fa
```

Place your OTUs into the 18S phylogeny.
*Make sure you have read the pplacer documentation to be sure these flag choices are appropriate for your questions.*
```
pplacer -o query.align.jplace -p --keep-at-most 10 -c refpkg query.align.fa
```
Use guppy to convert the json file to a parsible csv file and a phyloxml tree with edges fattened according to density of OTU placements.
```
guppy to_csv --point-mass --pp -o query.align.csv query.align.jplace
guppy fat --node-numbers --point-mass --pp -o query.align.phyloxml query.align.jplace
```
*Take alignment offline and check your OTUs were placed as expected along the full length 18S sequences (I expected my sequences to start around ~400 bp region),*
*and check that all of your OTUs are present in the results (for mine, 5072 OTUs + 1 header line)*
```
cut -d "," -f4 query.align.csv | wc -l
```
Create edgecode and add a new column to the query.align.csv file that interprets which edge each OTU was placed on
 *I'm sure there is a better way to do this, but the following commands creates edge_code.txt which is a list of the edges associated with each reference sequence.*
```
sed 's/,[A-Z]\|([A-Z]/\n&/g' query.align.jplace |head -n107 |tail -n105 |tr '{' '\n'|awk -F} '{print $1}' |awk -F: '{print $1}' | tr '\n' ' ' | sed 's/,//g' | sed 's/(//g' | sed 's/ / |/g' | sed 's/|[A-Z]/\n&/g' > edge_code.txt
```
Create a new column to add to the to the query.align.csv file that interprets which edge each OTU was placed on, will create a new file called query.align.wkey.csv
```
./tr_edges.sh query.align.csv
```

## 4. Process results
You can use the following provided in the repo to follow along:
- assign_lineage.sh (a shell script to add the lineages of strains from ncbi taxonomy and whether strains is high or low producer)
- dmspDB_key.csv (used by the shell script for reference)

Filter results and only keep OTUs placed with posterior probability of 90% and likelihood < -4,000
**Note that my cutoff for the likelihood differs between the two phylogenies. When I moved the alignment offline to an alignment viewer, I found that the posterior probability filter had removed a majority of poorly aligned sequences. However, a few sequences remained that were not as nicely aligned and they had significantly higher likelihood values than all of the other sequences, which is how I landed on both filters.**
```
awk -F, '$6>=0.9 && $7<=-4000' query.align.wkey.csv > query.align.wkey.filt.csv
echo "origin,name,multiplicity,edge_num,like_weight_ratio,post_prob,likelihood,marginal_like,distal_length,pendant_length,classification,map_ratio,map_overlap,map_identity,edge_name" | cat - query.align.wkey.filt.csv > query.align.wkey.filtered.csv
rm query.align.wkey.filt.csv
```
Add HiDP or LoDP key and simple lineages to the csv file (I used an arbitrary 'generic' taxonomy that are biased by my preferences, however the key file also contains the full taxonomy from ncbi if you prefer a different classification level).
```
./assign_lineage.sh query.align.wkey.filtered.csv
paste -d, query.align.wkey.filtered.csv edge_lineages.txt > query.align.final.csv
```
I translate to a tab delimited file out of preference
```
tr "," "\t" < query.align.final.csv > query.align.final.DPtype.txt
```

## Final product: With the steps above you have created query.align.final.DPtype.txt which contains the name of OTU, the statistics from pplacer associated with each OTU, the edge_name which is the name of the DMSP producer the OTU is most significantly related to, the HiDP or LoDP key of this DMSP producer, and the simplified taxonomy of this DMSP producer. The file has been filtered to only include OTUs that were significantly related to one of the known DMSP producers based on our filtering cut-off choices.