# Use

This repository is here to allow the annotation of 
  African Swine Fever Virus (ASFV) genomes without the use
  of Geneious or Glimmer. It annotates a genome by using
  a Waterman Smith alignment to map known genes from a
  reference genome with to the genome of interest. This
  takes about 20 minutes per ASFV genome on DDR4 ram.
  You should also double check the results from this
  programs output. It will make mistakes, especially in
  the inverted repeat region.

This repository is under a dual license, pick the license
  (CC0 or MIT) that works best for you.

# How to install

This program is dependent on revseq and transeq from
  EMBOSS, which can be installed with your package manager.
  For Debian or Ubuntu it would be
  `sudo apt-get install emboss`. For Mac, you will likely
  need to install homebrew. After installing EMBOSS run
  `make` to build alnSeq. You can move the scripts and
  compiled alnSer programs to a directory in your PATH if
  you want to or you could just call the annotateASFV.sh
  script from your downloaded directory.

```
apt-get install emboss;
git clone https://github.com/jeremybuttler/annotateASFV;
cd annotateASFV;
make;
```

Install files in /usr/local/bin (optional). I would
  recommend not doing this step.

```
sudo cp rmDupFeatures.awk /usr/local/bin/;
sudo cp detectIncompleteORF.awk /usr/local/bin;
sudo cp annotateASFV.sh /usr/local/bin;
sudo cp alnSeqMem alnSeqScan alnSeq /usr/local/bin;

sudo chmod a+x rmDupFeatures.awk detectIncompleteORF.awk;
sudo chmod a+x annotateASFV.sh alnSeqMem alnSeqScan alnSeq;
```

# Running annotateASFV.sh

## How to run

To run this repository just do
  `bash annotateASFV.sh genes.fasta ASFVgenome.fasta name`.
  Were genes.fasta are the genes to annotate your genome
  with. ASFVgenome.fasta is the genome to annotate. And
  name is the prefix to add to the output files. An example
  command would be

```
bash path/to/annotateASFV.sh\
  demo/OP605386-features-rmDup.fa\
  demo/OP605386.fa\
  demo/demo-OP605386;

# How I ran it for the demo data
cd annotateASFV;

bash annotateASFV.sh\
  demo/OP605386-features-rmDup.fa\
  demo/OP605386.fa\
  demo/demo-delete;
```

## Making a fasta file with genes

annotateASFV.sh has a couple rules about the fasta file
  with genes. The first rule is that there can be no
  spaces in the header (line with a ">"). The second rule
  is that each sequence takes up only one line.

Here is how I prepare a Genbank file for using with
  annotateASFV.sh. The first step is to extract the
  features from the Genbank (.gbk) file. I have used
   [https://www.bioinformatics.org/sms2/genbank_feat.html](
    https://www.bioinformatics.org/sms2/genbank_feat.html)
  to extract the features. I ignore the first feature from
  Genbank extract, because it is the entire genome. Also,
  the first few lines are often meta data.
 
Warning alnSeq will run out of memory if you leave the full
  ASFV genome in the file. You can get around this with
  `mv alnSeqMem alnSeq`. However, annotateASFV.sh will take   
  twice as long, but it will only take 9 Gb instead of
  40 Gb to align two complete ASFV genomes. The better
  option is just to make sure the full ASFV genome from the
  Genbank file is removed. The 7000 base CP2475L gene takes
  around 1.4 Gb of ram to map to the full ASFV genome. So,
  this program should use less than 2 Gb of ram to map each
  genes to the ASFV genome.

  
After I extract the features I use rmDupFeatures.awk
  `awk -f rmDupFeatures.awk features.fasta > cleaned.fasta`
  to remove the duplicate features, set up the names, and
  convert the sequences to one line entires. After that you 
  are ready to use the new features in annotate.sh.


# Output from annotateASFV.sh

I used the prefix demo/demo-delete. There are several files
  made by this program. Each file is there to provide some
  kind of information to you.

## prefix-scores.tsv

The prefix-scores.tsv file has a table of the main
  alignments and possible alternative alignments for each
  gene. I would recommend viewing this file in a spread
  sheet editor. I like sc, which can be installed by
  `sudo apt-get install sc`. To use sc do
  `cat demo/demo-delete-scores.tsv | psc | sc`. It has the
  same key bindings as less; q to quite, page up ctrl-b to
  move down a page, page down or ctrl-f to move down a
  page, and arrow keys or k/j/h/l to move
  up/down/left/right.

As you can see the repeat regions are a bit messy and map
  to both ends of the genome.

## prefix-feat.tbl

The next file is the prefix-feat.tbl. This is the feature
  table that you would submit to Genbank. It will need some
  editing however. In some cases you may have an ALT-FIX
  entry, were you need to make a decision on how to fix
  an entry. In another case you may have to correct an
  entry. From demo/demo-delete-feat.tbl you can see that
  the repeat and genes in the repeats need to be corrected.
  In this case you need to make a choice of alternate
  entries.

## prefix-errors.tsv

The next file is prefix-errors.tsv, which is a tsv file
  with the problems that were printed out at the end of the
  program. I would recommend again viewing this with a
  spread sheet editor like sc
  `cat demo/demo-delete-errors.tsv`. In this case there
  are a few genes without start codons. Since this is an
  Illumina reference, I am going to assume these are real.

## prefix-aa.fa

The next file is prefix-aa.fa (demo/demo-delete-aa.f). This
  file holds the ammino acid sequences for each gene or
  feature. It would be a good idea to ignore the inverted
  repeats, since these are non-coding. The ALT-FIX ending
  means that it was an alternative alignment.

## prefix-map.aln

The next file is prefix-map.aln (demo/demo-delete-map.aln).
  This file holds the alignments for each gene or each
  alternate genes. ALT-FIX does not appear in this file,
  however the ALT-FIX alignments are in this file.

## prefix-seq.fa

The final output file is prefix-seq.fa
  (demo/demo-delete-seq.fa). This file holds the sequences
  for each alignment in the feature table. In this case
  all the sequences are printed in open reading frame one.
  So, the start codon is the first base. The head for each
  sequence ends with number-number-1or0. Were the first
  number is the first reference base the gene mapped to.
  The second number is the last reference base the gene
  mapped to. The 1or0 marks if it is a forward reading
  frame (1) or a reverse complement reading frame (0). You
  will also notice that ALT-FIX also appears in this file
  as well.
