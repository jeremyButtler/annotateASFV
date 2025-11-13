# Use

This repository is here as a tool to annotate African
  Swine Fever Virus (ASFV) genomes. It is not the best
  method, but it can work. Annotation is done with a
  Waterman Smith alignment to map known genes from a
  reference genome with to the genome of interest.

Memory usage is light, but it is slow. With GNU time it
  took 8 minutes with 10 Mb of ram to annotate the demo
  ASFV genome with the `-slow` flag. Also, it does not
  handle genes that use exons and introns (ASFV does have
  use these) very well.

The `-fast` (default) form may not be quite as sensitive
  (no idea), but for me GNU time clocked it at 7 seconds
  and under 10 Mb of ram.

This repository is under a dual license, pick the license
  (Unlicense or MIT) that works best for you.

# How to install

```
git clone https://github.com/jeremybuttler/annotateASFV;
cd annotateASFV;
make;

The scripts and needed programs will be saved to the bin
  directory. You can place the contents were you want,
  but you must make sure all scripts and programs are in
  the same directory.

# Running annotateASFV.sh

## How to run

To run this repository do:

```
sh annotateASFV.sh \
    -gene genes.fasta \
    -ref ASFVgenome.fasta \
    -prefix name
```

Were genes.fasta are the genes to annotate your genome
  with. ASFVgenome.fasta is the genome to annotate. And
  name is the prefix to add to the output files.

A example is:

```
bash path/to/annotateASFV.sh \
  -gene demo/OP605386-features-rmDup.fa \
  -ref demo/OP605386.fa \
  -prefix demo/demo-OP605386;
```

## Making a fasta file with genes

annotateASFV.sh has a couple rules about the fasta file
  with genes. The first rule is that there can be no
  spaces in the header (line with a ">").

Here is how I prepare a Genbank file for using with
  annotateASFV.sh. The first step is to extract the
  features from the Genbank (.gbk) file. I have used
   [https://www.bioinformatics.org/sms2/genbank_feat.html](
    https://www.bioinformatics.org/sms2/genbank_feat.html)
  to extract the features. I ignore the first feature from
  Genbank extract, because it is the entire genome. Also,
  the first few lines are often meta data. Next I remove
  the duplicate sequences.
 
After I extract the features I use rmDupFeatures.awk
  `awk -f rmDupFeatures.awk features.fasta > cleaned.fasta`
  to remove the duplicate features, set up the names, and
  convert the sequences to one line.

# Output from annotateASFV.sh

I used the prefix demo/demo-delete. There are several files
  made by this program. Each file is there to provide some
  kind of information to you.

# Output:

## prefix-scores.tsv

Is a table of the main alignments and possible alternative
  alignments for each gene. I would recommend viewing this
  file in a spread sheet editor.

As you can see the repeat regions are a bit messy and map
  to both ends of the genome.

## prefix-highScores.tsv

This is like prefix-scores, except these passed the
  filters.

## prefix-featureTable.tsv

Is the feature table that you would submit to Genbank. It
  will need some editing however. In some cases you may
  have a `WARNING` entries. This make something that needs
  to be resolved. the two current warnings are duplicates
  and incomplete reading frames (missing start/stop).

## prefix-featureErrors.tsv

Is a tsv file with the problems that were found in
  prefix-featureTable.tsv. It is a hit list of things to
  correct or double check.

## prefix-aaSequences.fa

Has the ammino acid sequences for each gene or feature
  that was not marked `non-coding`.

## prefix-ntSequences.fa

Has the nucleotides sequences for each feature.
