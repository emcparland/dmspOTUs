# Pplacer of OTU sequences in MMETSP 18S sequence phylogeny to assign putative DMSP synthesis gene carriers.

## 
My reference sequences are full-length 18S sequences from the MMETSP strains, which were obtained from [iMicrobe](https://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/18s/18s.fa). (For previous blast analyses to determine which MMETSP contain DMSP synthesis genes, I used the recently re-assembled MMETSP transcriptomes ([Johnson, Alexander and Brown 2018](https://academic.oup.com/gigascience/article/8/4/giy158/5241890)). If you'll be using my previously made phylogeny of MMETSP 18S sequences, you can *skip to Step 3*. 

You will have n=655 18S sequences from the MMETSP transcriptomes. (Note, this is less than the recently updated n=678 transcriptomes in the Johnson study).
```
grep "^>" 18s.fa | wc -l
```
Every transcriptome has a representative 18S sequence meaning there are duplicate 18S sequences, use the shell script below to keep only one sequence per strain (n=393).
```
grep "^>" 18s.fa | sort -k1.12 | uniq -s11 |wc -l
```
Find and remove duplicates strains
```
grep "^>" 18s.fa > all_id.txt
sort -k1.12 all_id.txt | uniq -s11 > unique_id.txt
./extract_erin.sh unique_id.txt 18s.fa
```
Check that it worked
```grep "^>" 18s_all_mmetsp_unique.fa > tmp.txt
comm -12 tmp.txt unique_id.txt |wc -l
rm tmp.txt
```
Clean up symbols
```
sed -i 's/-/_/g' 18s_all_mmetsp_unique.fa
sed -i 's/|/_/g' 18s_all_mmetsp_unique.fa
sed -i 's/=/_/g' 18s_all_mmetsp_unique.fa 
sed -i 's/:/_/g' 18s_all_mmetsp_unique.fa
sed -i 's/ /_/g' 18s_all_mmetsp_unique.fa
```
Check again
```
grep "^>" 18s_all_mmetsp_unique.fa |wc -l
```

## 1. Create alignment of reference sequences with Infernal using the 18S ssu cov model as a reference for alignment.
Name of your fasta file with reference sequences
```
NAME_FA=18s_all_mmetsp_unique.fa
```
Degap if necessary
```
seqmagick mogrifty --ungap $NAME_FA
```
Use infernal to align
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

Name of alignment to easily sed any extract characters in alignment that will mess up next steps and name header for new tree files
```
NAMEFA=mmetsp.align.nogap.full.fa
NAMET=mmetsp
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

If starting with this step, use the following that I have provided in repo:
- mmetsp.align.nogap.full.fa  (reference alignment)
- mmetsp.refpkg (reference package created above for pplacer)
- tr_edge.sh (shell script to translate which edges are associated with which strain)

Align the query reads and reference alignment
First, clean names in OTU file
```
sed -i 's/ /_/g' otus.fa
```
If yours have any funky symbols you could use this handy command from the Bowman tutorial
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
Use infernal to align
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
pplacer -o query.align.jplace -p --keep-at-most 10 -c mmetsp.refpkg query.align.fa
```
Use guppy to convert the json file to a parsible csv file and a phyloxml tree with edges fattened according to density of OTU placements.
```
guppy to_csv --point-mass --pp -o query.align.csv query.align.jplace
guppy fat --node-numbers --point-mass --pp -o query.align.phyloxml query.align.jplace
```
*Take alignment offline and check your OTUs were placed as expected along the full length 18S sequences (I expected my sequences to start around ~400 bp region),*
*and check that all of your OTUs are present in the results (5072 OTUs + 1 header line)*
```
cut -d "," -f4 query.align.csv | wc -l
```
Create edgecode to interpret which edge each OTU was placed on
 *There are multiple edges in the reference seqs*
```
grep "MMETSP" query.align.jplace | sed 's/MM/ MM/g' | tr ' ' '\n' | tr '{' '\n' | awk -F} '{print $1}' | awk -F: '{print $1}' | tr '\n' ' ' | sed 's/ / |/g' | sed 's/|MME/\nMME/g' | grep -v "((" > edge_code.txt
```
Create a new file with column that interprets which edge each OTU was placed on
```
./tr_edges.sh query.align.csv
```

## 4. Process results
Filter results and only keep OTUs placed with posterior probability of 90% and likelihood < -10,000
**Note that my cutoff for the likelihood differs between the two phylogenies. When I moved the alignment offline to an alignment viewer, I found that the posterior probability filter had removed a majority of poorly aligned sequences. However, a few sequences remained that were not as nicely aligned and they had significantly higher likelihood values than all of the other sequences, which is how I landed on both filters.**

- assign_lineage.sh (a shell script to add the lineages of MMETSP strains)
- mmetsp_all_lineages.csv (reference for taoxnomy of MMETSP strains used in shell script)
- mmetspkey.txt (used to assign any associated DMSP synthesis genes to each MMETSP strain in query file)

```
awk -F, '$6>=0.9 && $7<=-10000' query.align.wkey.csv > query.align.wkey.filt.csv
echo "origin,name,multiplicity,edge_num,like_weight_ratio,post_prob,likelihood,marginal_like,distal_length,pendant_length,classification,map_ratio,map_overlap,map_identity,edge_name" | cat - query.align.wkey.filt.csv > query.align.wkey.filtered.csv
rm query.align.wkey.filt.csv
```
Add lineages to the key (I used an arbitrary 'generic' taxonomy that are biased by my preferences, however the lineages csv file also contains the full taxonomy from ncbi if you prefer a different classification level).
```
./assign_lineage.sh query.align.wkey.filtered.csv
paste -d, query.align.wkey.filtered.csv edge_lineages.txt > query.align.final.csv
```
Add gene assignment
```
length=`(wc -l < query.align.final.csv |awk '{print $1}')`
use="$((length-1))"
tail -n $use query.align.final.csv |cut -f15 -d',' | sed -E 's/^MMETSP[0-9]+_//g' | tr "_" "-" | while read headername
do
	grep $headername mmetspkey.txt > tmp
	length=`(wc -l < tmp |awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		cat tmp >> gene_key
	else
		echo "blank blank blank blank" | tr " " "\t" >> gene_key
	fi
done
echo "DSYB TpMT2 TpMT1 MMETSPstrain" | tr " " "\t" > header
cat header gene_key |awk '{print $1,$2,$3}' | tr " " "\t" > gene_keytmp
tr "," "\t" < query.align.final.csv > tmp
paste gene_keytmp tmp > query.align.final.genes.txt
```
## Final product
With the steps above you have created query.align.final.genes.txt which contains the name of OTU, the statistics from pplacer associated with each OTU, the edge_name which is the name of the MMETSP strain the OTU is most significantly related to, the presence (1) or absence (0) of DSYB, TpMT1 or TpMT2 in the associated transcriptome of the MMETSP strain, and the taxonomy of this MMETSP strain. The file has been filtered to only include OTUs that were significantly related to one of the MMETSP strains based on our filtering cut-off choices.