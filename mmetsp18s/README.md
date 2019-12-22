# Putative DMSP synthesis of 18S OTUs with pplacer

A large portion of this analysis was made possible by an awesome blog post from the [Bowman lab](https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/) and helpful discussions with Harriet Alexander.

## You'll need to install:

## You'll need:

1. Your 18S OTUs
Mine happen to be here from a long time ago and I'll duplicate them into a new file so as not to accidentally ruin the hardwork of Erin 3 years ago:
cp ~/hu_tutorial/all_seq/pick_open/rep_set.fna otu_seqs.fa

I should have n=5072 OTU sequences
```grep "^>" otu_seqs.fa |wc -l```


2. The MMETSP 18S sequences
My reference sequences are full-length 18S sequences from the transcriptomes of the MMETSP database. For blast analyses I used the most recent MMETSP transcriptomes (), for the 18S sequences, I obtained the sequences from [iMicrobe](https://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/18s/18s.fa)

You should have n=655 18S sequences from the MMETSP transcriptomes. (Again, this will be less than the n=678 you would have from the Johnson re-analysis work).
```grep "^>" 18s.fa | wc -l ```

Every transcriptome has a representative 18S sequence meaning there are duplicate 18S sequences, keep only the unique strains (n=393).
```grep "^>" 18s.fa | sort -k1.12 | uniq -s11 |wc -l ```


## 1. Create alignment of reference sequences with Infernal using the 18S ssu cov model as a reference for alignment.

* Run part1.txt to align with cmalign, convert the output to a fasta output and deduplicate any repeat sequences (though there shouldn't be any).

* I then take the alignment offline to manually curate and remove any large gaps in the alignment with Jalview or Geneious.

## 2. Create a phylogenetic tree of the alignment, reference package from the alignment, tree and stats file with pplacer

Per Bowman lab tutorial: Taxtastic can't read a RAxML stats file with confidence values, work around this by building tree without scores first than calculate scores separately and add back to your previously generated tree.
* First build the tree without confidence scores
* Root the tree with RAxML (or do this yourself manually)
* Generate the confidence scores for tree separately
* Create reference package with taxit command

## 3. Phylogenetic placement of 18S OTUs onto reference
Now place your OTUs into the 18S phylogeny! 
* Use pplacer command to place OTUS. Make sure you have read the pplacer documentation to be sure these flag choices are appropriate for your questions. 
* Then use guppy to convert the json file to a parsible csv file and a phyloxml tree with edges fattened according to density of OTU placements.
* Check that all of your OTUs are present in the results (5072 OTUs + 1 header line)
``` cut -d "," -f4 query.align.csv | wc -l ```
* Create edge_code.txt and add a new column to the query.align.csv file that interprets which edge each OTU was placed on
* Filter results and only keep OTUs placed with posterior probability of 90%








