# Assigning putative DMSP production/synthesis to environmental 18S OTU sequences

## Quick summary (Work in progress)
(This is my first ever official, public repo! Let me know if anything is confusing.)

Question: Who are the eukaryotic DMSP producers present in situ?

Motivation:
- Prior to the recent publication of the DMSP synthesis genes, our knowledge of the ability to produce DMSP was limited to screening monocultures (see Figure 1 in [McParland and Levine 2019]((https://aslopubs.onlinelibrary.wiley.com/doi/full/10.1002/lno.11076)).
- The recently published DMSP synthesis genes encode for the third step of the synthesis pathway which is a methyltransferase enzyme and have been reported in heterotrophic bacteria ([dsyb)](https://www.nature.com/articles/nmicrobiol20179?platform=oscar&draft=collection)), and eukaryotic protists ([DSYB](https://www.nature.com/articles/s41564-018-0119-5) and [TpMT1/TpMT2](https://www.sciencedirect.com/science/article/abs/pii/S0003986118300080?via%3Dihub)).
- Although these genes have been published, we can't quite yet measure the eukaryotic DMSP synthesis genes in mixed natural communities as there are no universal primers for the in situ eukaryotic community.

Method:
- We measured the diversity and relative abundance of in situ eukaryotic communities using the 18S gene (V4 region).
- I assigned the OTUs a putative DMSP phenotype or genotype using pplacer ([Matsen, Kodner and Armbrust 2010](https://matsen.fhcrc.org/papers/matsen2010pplacer.pdf)) by aligning the OTU sequences with one of two reference alignments:
	1. dmsp18s: I previously created a database of all known measurements of DMSP production in monocultures ([Supp Table 1 in McParland and Levine 2019 2019](https://aslopubs.onlinelibrary.wiley.com/doi/full/10.1002/lno.11076)). I created a reference alignment of this database using the 18S sequences of the measured strains retrieved from ncbi.
	2. mmetsp18s: I previously used blast to determine which of the 400+ strains in the MMETSP database ([Keeling et al. 2014](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889)) contain one (or more) of the eukaryotic DMSP synthesis genes. I created a reference alignment of the 18S sequences of the MMETSP strains.

- If an OTU is significantly aligned with:
	1. one of the known DMSP producer sequences, assume that the OTU can also produce DMSP.
	2. one of the MMETSP sequences and the respective transcriptome contained one of the DMSP synthesis genes, assume that the OTU also contains the synthesis gene.

CREDIT: 
A large portion of this analysis was made possible by an awesome blog post from the [Bowman lab](https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/) and  helpful discussions with Harriet Alexander. Science questions were driven by ideas discussed with Naomi M. Levine. The 18S sequences were compiled previously and will be presented elsewhere, but many thanks to the pipeline of Sarah Hu.


## Pipeline description 

If you want to use one of the phylogenies I have already built as described above, you can *skip to step 3* to place your own OTUs within the phylogeny. (If you want to create your own phylogeny, start with step 1).

## You'll need to install:
- [pplacer](https://matsen.fhcrc.org/pplacer/)
- [guppy](https://matsen.github.io/pplacer/generated_rst/guppy.html)
- [seqmagick](http://fhcrc.github.io/seqmagick/)
- [infernal](http://eddylab.org/infernal/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) (can skip if starting at step 3)

## You'll need:
- [Rfam](https://rfam.xfam.org/family/RF01960) (RF01960) covariance model for euks small subunit rRNA is (provided in repos)
```
cmconvert eukarya-0p1.cm  > eukarya-0p1.conv.cm
```
- Your fasta file of 18S OTU sequences

We used the V4 region, however the phylogeny is built with full length 18S sequences so it could be easily adapted for V9 region primers.

My sequences happen to be here from a long time ago and I duplicated them into a new file so as not to accidentally ruin the hardwork of Erin 3 years ago:
```
cp ~/hu_tutorial/all_seq/pick_open/rep_set.fna otu_seqs.fa
```
I should have n=5072 OTU sequences
```
grep "^>" otu_seqs.fa |wc -l
```
- All other files needed for starting at Step 3 are found and further explained in respective directories (dmsp18s and mmetsp18s)
