# Putative DMSP production of 18S OTUs with pplacer
## Quick summary
(work in progress)

Question: Which eukaryotic DMSP producers are present in a natural mixed community?

Motivation: 
- We can't quite yet measure the eukaryotic DMSP synthesis genes in mixed natural communities.
- Until ([2018](https://www.nature.com/articles/s41564-018-0119-5)) we didn't know the DMSP synthesis genes and still, there are no universal primers for the in situ eukaryotic community.
- We can assess the diversity and relative abundance of the eukaryotic community in mixed natural communities with the 18S gene (V4 region here).

Method:
- I previously created a database of all known measurements of DMSP production in monocultures ([Supp Table 1 in McParland and Levine 2019 2019](https://aslopubs.onlinelibrary.wiley.com/doi/full/10.1002/lno.11076)). 
- Create a reference alignment of these known DMSP producers' 18S gene, then use pplacer ([Matsen, Kodner and Armbrust 2010](https://matsen.fhcrc.org/papers/matsen2010pplacer.pdf)) to align the environmental 18S OTUs with this reference alignment.
- If an OTU is significantly aligned with one of the reference sequences, assume that the OTU could also contain the synthesis gene

A large portion of this analysis was made possible by an awesome blog post from the [Bowman lab](https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/) and helpful discussions with Harriet Alexander.

- If you want to use the phylogeny I built from the 18S sequences of known DMSP producers, you could skip to step 3 to place your own OTUs within that phylogeny.
- If you want to create your own phylogeny, start with step 1.
- This is my first ever official, public repo! Let me know if anything needs clarity.

## You'll need to install:
- [pplacer](https://matsen.fhcrc.org/pplacer/)
- [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html)
- [seqmagick](http://fhcrc.github.io/seqmagick/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (can skip if starting at step 3)

## You'll need:
- [Rfam](https://rfam.xfam.org/family/RF01960) (RF01960) covariance model for euks small subunit rRNA is (provided here for convenience, can skip if starting at step 3)
```
cmconvert eukarya-0p1.cm  > eukarya-0p1.conv.cm
```
- Your 18S OTUs

Mine happen to be here from a long time ago and I duplicated them into a new file so as not to accidentally ruin the hardwork of Erin 3 years ago:
```
cp ~/hu_tutorial/all_seq/pick_open/rep_set.fna otu_seqs.fa
```
I should have n=5072 OTU sequences
```
grep "^>" otu_seqs.fa |wc -l
```

- The 18S sequences of known DMSP producers (provided here for convenience, I previously searched for respectives sequences on ncbi to create Figure 1 in McParland and Levine 2019).

You should have n=107 18S sequences.
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
cat mmetsp.align.nogap.full.fa  otus.fa > query.fa
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
 *I'm sure there is a better way to do this, but the following commands creates 'edge_code.txt' where there are multiple edges in the reference seqs which will then be assigned to the csv file of OTUs*
```
sed 's/,[A-Z]\|([A-Z]/\n&/g' query.align.jplace |head -n107 |tail -n105 |tr '{' '\n'|awk -F} '{print $1}' |awk -F: '{print $1}' | tr '\n' ' ' | sed 's/,//g' | sed 's/(//g' | sed 's/ / |/g' | sed 's/|[A-Z]/\n&/g' > edge_code.txt
```
Create a new column to add to the to the query.align.csv file that interprets which edge each OTU was placed on, will create a new file called query.align.wkey.csv
```
./tr_edges.sh query.align.csv
```

## 4. Process results
For the following I have provided in the rep:
- assign_lineage.sh (a shell script to add the lineages of strains from ncbi taxonomy and whether strains is high or low producer)
- dmspDB_key.txt (used by the shell script for reference)

Filter results and only keep OTUs placed with posterior probability of 90% and likelihood < -4,000
```
awk -F, '$6>=0.9 && $7<=-4000' query.align.wkey.csv > query.align.wkey.filt.csv
echo "origin,name,multiplicity,edge_num,like_weight_ratio,post_prob,likelihood,marginal_like,distal_length,pendant_length,classification,map_ratio,map_overlap,map_identity,edge_name" | cat - query.align.wkey.filt.csv > query.align.wkey.filtered.csv
rm query.align.wkey.filt.csv
```
Add lineages to the key
*** START HERE ALMOST DONE**
```
./assign_lineage.sh query.align.wkey.filtered.csv
paste -d, query.align.wkey.Filt.csv edge_lineages.txt > query.align.final.csv
```