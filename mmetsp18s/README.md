# Putative DMSP synthesis of 18S OTUs with pplacer
## Quick summary
Question: Which eukaryotic DMSP producers are present in a natural mixed community?
Motivation: 
- We can't quite yet measure the eukaryotic DMSP synthesis genes in mixed communities.
- But we do ([finally](https://www.nature.com/articles/s41564-018-0119-5)) know what those genes are!
- We can assess the diversity and relative abundance of the eukaryotic community in mixed natural communities with the 18S gene (V4 region here).
Method:
- I previously used blast to determine which of the 400+ strains in the MMETSP database ([Keeling et al. 2014](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889))contain one (or more) of the DMSP synthesis genes. 
- Create a reference alignment of the 18S gene of each MMETSP strain, then use pplacer ([Matsen, Kodner and Armbrust 2010](https://matsen.fhcrc.org/papers/matsen2010pplacer.pdf)) to align the environmental 18S OTUs with this reference alignment.
- If an OTU is significantly aligned with one of the reference sequences and the respective reference sequence also contains one of the DMSP synthesis genes, assume that the OTU could also contain the synthesis gene

A large portion of this analysis was made possible by an awesome blog post from the [Bowman lab](https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/) and helpful discussions with Harriet Alexander.

- If you want to use the phylogeny I built from the 18S sequences of MMETSP transcriptomes, you could skip to step 3 to place your own OTUs within that phylogeny.
- If you want to create your own phylogeny, start with step 1.
- This is my first ever official, public repo! Let me know if anything needs clarity.

## You'll need to install:
- [pplacer](https://matsen.fhcrc.org/pplacer/)
- [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html)
- [seqmagick](http://fhcrc.github.io/seqmagick/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (can skip if starting at step 3)

## You'll need:
- Your 18S OTUs
Mine happen to be here from a long time ago and I duplicated them into a new file so as not to accidentally ruin the hardwork of Erin 3 years ago:
```cp ~/hu_tutorial/all_seq/pick_open/rep_set.fna otu_seqs.fa```
I should have n=5072 OTU sequences
```grep "^>" otu_seqs.fa |wc -l```

- The MMETSP 18S sequences
My reference sequences are full-length 18S sequences from the transcriptomes of the MMETSP database. For blast analyses I used the recently re-assembled MMETSP transcriptomes ([Johnson, Alexander and Brown 2018](https://academic.oup.com/gigascience/article/8/4/giy158/5241890)), for the 18S sequences, I obtained the sequences from [iMicrobe](https://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/18s/18s.fa)
You should have n=655 18S sequences from the MMETSP transcriptomes. (Note, this is less than the n=678 you would have from the Johnson study).
```grep "^>" 18s.fa | wc -l ```
Every transcriptome has a representative 18S sequence meaning there are duplicate 18S sequences, use the shell script below to keep only the unique strains (n=393).
```grep "^>" 18s.fa | sort -k1.12 | uniq -s11 |wc -l ```
Find and remove duplicates strains
```grep "^>" 18s.fa > all_id.txt```
```sort -k1.12 all_id.txt | uniq -s11 > unique_id.txt```
```./extract_erin.sh unique_id.txt 18s.fa```
Check that it worked
```grep "^>" 18s_all_mmetsp_unique.fa > tmp.txt```
```comm -12 tmp.txt unique_id.txt |wc -l```
```rm tmp.txt```
Clean up symbols
```sed -i 's/-/_/g' 18s_all_mmetsp_unique.fa```
```sed -i 's/|/_/g' 18s_all_mmetsp_unique.fa```
```sed -i 's/=/_/g' 18s_all_mmetsp_unique.fa```  
```sed -i 's/:/_/g' 18s_all_mmetsp_unique.fa```
```sed -i 's/ /_/g' 18s_all_mmetsp_unique.fa```
Check again
```grep "^>" 18s_all_mmetsp_unique.fa |wc -l
```

## 1. Create alignment of reference sequences with Infernal using the 18S ssu cov model as a reference for alignment.
Name of your fasta file with reference sequences
```
NAME_FA=DMSP_prod_db_18S.fa
```
Degap if necessary
```seqmagick mogrifty --ungap $NAME_FA```
- ??????
```cmconvert eukarya-0p1.cm  > eukarya-0p1.conv.cm```
- Align
```cmalign --dna -o align.sto --outformat Pfam eukarya-0p1.conv.cm $NAME_FA```
- Convert to fasta format
```seqmagick convert align.sto align.fa```
- Deduplicate sequences
```seqmagick mogrify --deduplicate-sequences align.fa```

**I then take the alignment offline to manually curate and remove any large gaps in the alignment with Jalview or Geneious.**

## 2. Create a phylogenetic tree of the alignment, reference package from the alignment, tree and stats file with pplacer
- Name of alignment AND to easily sed any extract characters in alignment that will mess up next steps
```NAMEFA=align.fa```
```NAMET=dmsp```
```sed -i 's/|/_/g' $NAMEFA```
```sed -i 's/=/_/g' $NAMEFA```
```sed -i 's/:/_/g' $NAMEFA```
```sed -i 's/(/_/g' $NAMEFA```
```sed -i 's/)/_/g' $NAMEFA```
```#sed -i 's/-/_/g' $NAMEFA```
```sed -i 's/\//_/g' $NAMEFA```
```sed -i 's/\.//g' $NAMEFA```
Per Bowman tutorial: Taxtastic can't read a RAxML stats file with confidence values, work around this by building tree without scores first than calculate scores separately and add back to your previously generated tree.
- First build the tree without confidence scores
```raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -s $NAMEFA -n $NAMET.tre -f d -p 12345```
- Root the tree with RAxML
```raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -f I -t RAxML_bestTree.$NAMET.tre -n root_$NAMET.tre```
- Generate the confidence scores for tree separately
```raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root_$NAMET.tre -n root_conf_$NAMET.tre -s $NAMEFA```
- Create reference package with taxit command
```taxit create -l 18S_rRNA -P refpkg --aln-fasta $NAMEFA --tree-stats RAxML_info.$NAMET.tre  --tree-file RAxML_fastTreeSH_Support.root_conf_$NAMET.tre```

## 3. Phylogenetic placement of 18S OTUs onto reference
- Align the query reads and reference alignment
- First, clean names in OTU file
```sed -i 's/ /_/g' otus.fa```
- If yours have any funky symbols you could use this handy command from the Bowman tutorial
```tr "[ -%,;\(\):=\.\\\*[]\"\']" "_" < otus.fa ```
- Concatenate the query reads and reference alignment
```cat align.fa otus.fa > query.fa```
- Remove gaps
```seqmagick mogrify --ungap query.fa```
- Align
```cmalign --dna -o query.sto --outformat Pfam eukarya-0p1.conv.cm query.fa```
- Convert to fasta
```seqmagick convert query.sto query.align.fa```

- Place your OTUs into the 18S phylogeny!
*Make sure you have read the pplacer documentation to be sure these flag choices are appropriate for your questions.*
```pplacer -o query.align.jplace -p --keep-at-most 10 -c refpkg query.align.fa```
- Use guppy to convert the json file to a parsible csv file and a phyloxml tree with edges fattened according to density of OTU placements.
```guppy to_csv --point-mass --pp -o query.align.csv query.align.jplace```
```guppy fat --node-numbers --point-mass --pp -o query.align.phyloxml query.align.jplace```
- Take alignment offline and check your OTUs were placed as expected along the full length 18S sequences (I expected my sequences to start around ~400 bp region)
- Check that all of your OTUs are present in the results (5072 OTUs + 1 header line)
``` cut -d "," -f4 query.align.csv | wc -l ```

- Create edgecode and add a new column to the query.align.csv file that interprets which edge each OTU was placed on
 *I'm not sure why there are multiple edges in the reference seqs, but here parsing those out so I can assign the csv file of OTUs to be references*
remove first line which is weird symbols
```grep "MMETSP" query.align.jplace | sed 's/MM/ MM/g' | tr ' ' '\n' | tr '{' '\n' | awk -F} '{print $1}' | awk -F: '{print $1}' | tr '\n' ' ' | sed 's/ / |/g' | sed 's/|MME/\nMME/g' | grep -v "((" > edge_code.txt```
**OR**
```	sed 's/,[A-Z]\|([A-Z]/\n&/g' query.align.jplace |head -n107 |tail -n105 |tr '{' '\n'|awk -F} '{print $1}' |awk -F: '{print $1}' | tr '\n' ' ' | sed 's/,//g' | sed 's/(//g' | sed 's/ / |/g' | sed 's/|[A-Z]/\n&/g' > edge_code.txt```
- Create a new column to add to the to the query.align.csv file that interprets which edge each OTU was placed on, will create a new file called query.align.wkey.csv
```./tr_edges.sh query.align.csv ```

## 4. Process results
- Filter results and only keep OTUs placed with posterior probability of 90% and likelihood < -10,000 (for the alignment using known DMSP producers I used a likelihood of < -4000)
```awk -F, '$6>=0.9 && $7<=-10000' query.align.wkey.csv > query.align.wkey.filt.csv```
```echo "origin,name,multiplicity,edge_num,like_weight_ratio,post_prob,likelihood,marginal_like,distal_length,pendant_length,classification,map_ratio,map_overlap,map_identity,edge_name" | cat - query.align.wkey.filt.csv > query.align.wkey.Filt.csv```
```rm query.align.wkey.filt.csv```
- Add lineages to the key
```./assign_lineage.sh query.align.wkey.filt.csv```
```paste -d, query.align.wkey.Filt.csv edge_lineages.txt > query.align.final.csv```
- Add gene assignment
```length=`(wc -l < query.align.final.csv |awk '{print $1}')````
```use="$((length-1))"```
```tail -n $use query.align.final.csv |cut -f15 -d',' | sed -E 's/^MMETSP[0-9]+_//g' | tr "_" "-" | while read headername
do
	grep $headername mmetspkey.txt > tmp
	length=`(wc -l < tmp |awk '{print $1}')`
	if [ $length -gt 0 ]
	then
		cat tmp >> gene_key
	else
		echo "blank blank blank blank" | tr " " "\t" >> gene_key
	fi
done```
```echo "DSYB TpMT2 TpMT1 MMETSPstrain" | tr " " "\t" > header```
```cat header gene_key |awk '{print $1,$2,$3}' | tr " " "\t" > gene_keytmp```
```tr "," "\t" < query.align.final.csv > tmp```
```paste gene_keytmp tmp > query.align.final.genes.txt```






