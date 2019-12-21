# Putative DMSP synthesis of 18S OTUs with pplacer

A large portion of this analysis was made possible by an awesome blog post from the [Bowman lab](https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/) and helpful discussions with Harriet Alexander.

## 1. Create alignment of reference sequences with Infernal.
### My reference sequences are full-length 18S sequences from the transcriptomes of the MMETSP database. For blast analyses I used the most recent MMETSP transcriptomes (), for the 18S sequences, I obtained the sequences from [iMicrobe](https://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/18s/18s.fa)

You should have 655 18S sequences from the MMETSP transcriptomes
```grep "^>" 18s.fa | wc -l ```

## 2. Create a phylogenetic tree of the alignment

## 3. Create a reference package from the alignment, tree and stats file with pplacer

## 4. Phylogenetic placement of 18S OTUs onto reference