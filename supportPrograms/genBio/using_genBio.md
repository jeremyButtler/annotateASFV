# Goal:

This is to provide a quick table of contents for
  the `using_<name>.md` guides in genBio. It also includes
  a breif description of the goal of each library.

- adjCoords (see using\_adjCoords.md):
  - adjusing sam file reference position from local
    (mapped to a gene) to global (full genome coordinates)
- ampDepth (see using\_amdDepth.md):
  - gets read depths for a reference or if have gene
    coordinates mean and median read depths for genes
- cigToEqx (see using\_cigToEqx.md):
  - converts regular cigar entries ("MDISH") to eqx cigar
    entries ("=XDISH")
- codonFun (see using\_codonFun.md):
  - converts nucleotide sequences to amino acid sequences
- codonTbl (see using\_codonTbl.md):
  - has table to convert a set of three nucleotides to an
    amino acid
- edDist (see using\_edDist.md):
  - finds edit distance (excluding small indels) between
    reference and read in a sam file
  - can also compare two reads in a sam file, but that is
    much less accurate
- geneCoord (see using\_geneCoord.md):
  - using the geneCoord to store and get coordinates for
    genes in a genome, you must input a gene coordaintes
    structure
- gzSeqST (see using\_gzSeqST.md):
  - used to read sequences from fastq, fastq.gz, fasta,
    and fasta.gz files
- kmerBit (see using\_kmerBit.md):
  - has table and macros to convert nucleotide sequences
    to kmers (max of 15 nt)
- kmerCnt (see using\_kmerCnt.md):
  - functions to get number of kmers in a sequence
- maskPrim (see using\_maskPrim.md):
  - functions to mask primers in sam file reads by input
    coordinates
- ntTo2Bit (see using\_ntTo2Bit.md):
  - convert nucleotides to two bit values, which allows
    storing 4 nucleotides per byte
- revNtTo2Bit (see using\_ntTo2Bit.md):
  - convert nucleotides to two bit values, which allows
    storing 4 nucleotides per byte
  - this library also reverese complements
- ntTo5Bit (see using\_ntTo2Bit.md):
  - convert nucleotides to five bit values, not super
    efficent, but allows quick checking of anonymous bases
- rmHomo (see using\_rmHomo.md):
  - left align indels homopolymers or remove small indels
    in homopolymers
- samEntry (see using\_samEntry.md):
  - has stucture to read sam files, print sam file entries
    to store sam file entries, or convert sam file entries
  - think sam file manipulation
- samRef (see using\_samRef.md):
  - get reference sequence names and lengths from sam file
    headers (has its own structure)
- seqST (see using\_seqST.md):
  - supports reading sequences from fastq or fasta files
    into its own internal structure
  - also allows for reverse complementing sequences
- tbCon (see using\_tbCon.md):
  - convert sam file mapped reads into a reference based
    consensus (does not handle insertions well)
- trimSam (see using\_trimSam.md):
  - trim reads in a sam file
  - includes: start/end soft masking, coordinates, or
    set number at ends.
