# Use

AlnSeq uses a Smith Waterman and Needleman Wunsch alignment
  that, depending on flags used when compiling can run
  with memory usage of O(n \* m / 4) Bytes to O(n \* m)
  Bytes. However, the O(n \* m / 4) is twice as slow as a
  traditional Smith Waterman and Needleman Wunsch
  alignment. AlnSeq also includes an Hirschberg alignment.

AlnSeq is a standalone program that can also be compiled
  as a python library.

There are faster and less memory hungry Waterman Smith
  implementations than alnSeq. One example is the stripped
  Waterman Smith alignment, which I think reduces both
  scoring and direction matrix to just a few rows. See 
  [https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
  for an very fast striped Waterman Smith aligner. This
  one does have an issue when the alignment exceeds 33,000
  bases. However, I think there are alternatives without
  this issues.

This program is dual licensed for MIT and CC0. Pick the
  license which works best for you.

# Current work

Right now I am back to work on this project. Things will
  go a bit slowly, but I hope to get this finshed.

Currently I need to fix what I broke in my 20231022
  update. The broken programs are the query/ref scans, the
  two bit Hirschberg (only used if you compliled with
  -DHIRSCHTWOBIT), and the two bit memory efficent Smith
  Waterman (only used if you compiled with -DTWOBITMSW).
  To get around this I included the last version of alnSeq
  (20230908) in this repository, which does work.

The 20231022 update did make the Hirschberg faster, and
  provided a small speed boost to the other programs.

# Building and running alnSeq

## C standalone

### How to complie alnSeq as a standalone

```
# How to install this program

# Install at /usr/local/bin (need root privilage)
make
sudo make install

# Use a different install path
make
make install PREFIX=/path/to/install

# Manual install
make
mv alnSeq /path/to/install
chmod a+x /path/to/install/alnSeq

# Alternative make commands
make fast
  # This command disables the opening gap penalty,
  # directional selection, and two bit arrys. This makes
  # alnSeq behave somewhat similar to bio-alignment.
make mid
  # This command disables directional selection and two
  # bit arrays, but keeps the gap opening penalty. I would
  # recomend using this command.
```

### Standalone, extra build options (standalone only)

The flags alnSeq can be compiled with are:

- -DNOGAPOPEN
  - disable gap opening penalty
  - This is faster, but may produce lower quality
    alignments.
- -DBYTEMATRIX (default option)
  - Waterman and Needleman use a byte matrix instead of
    a two bit array matrix
  - Alignment takes half the time, but also takes 4x more
    memory.
- -DHIRSCHTWOBIT
  - Have the Hirschberg use two bit arrays instead of byte
    arrays for directions (only relevant if not using
    -DNOGAPOPEN).
  - This doubles the time to make an alignment for a very
    minor decrease in memory usage. 
- -DNOSEQCNVT
  - This prevents the conversion of each base in the query
    and reference sequences to an index for alignment.
  - This option slows down the alignment slightly. The only
    reason to use this option would be if the input case
    of a sequence matters.
- -DTWOBITMSW
  - Compiles the more memory efficient smith waterman with
    two bit arrays. This will have little effect on memory,
    but will slow it down.
- These options force alnSeq to prefer only one direction
  and disables all other options. This does speed up alnSeq
  slightly.
  - It only applies when the insertion (ins), snp, or 
    deletion (del) scores for a base pair are equal.
  - -DSNPINSDEL
  - -DSNPDELINS
  - -DINSSNPDEL
  - -DDELSNPINS
  - -DINSDELSNP (on by default)
  - -DDELINSSNP

You can compile with these flags using
  `make CFLAGS="flag"`. You can also compile multiple 
  flags with `make CFLAGS="falg1 flag2"`. One example is
  the `make fast` command, which uses make
  `CFLAGS="-DNOGAPOPEN -DINSSNPDEL -DBYTEMATRIX"`.

### How to run standalone alnSeq

```
# help message
alnSeq -h | less

# Get flags set when compiling alnSeq
alnSeq -flags | less

# Alignment formats 

## For a global alignment (Needleman Wunsch)
alnSeq -query query.fasta -ref ref.fasta > alignment.aln

## For a Hirschberg global alignment (slow, but low memory
# usage)
alnSeq -use-hirschberg -query query.fasta -ref ref.fasta > alignment.aln

## For a single local alignment (Waterman Smith)
alnSeq -use-water -query query.fasta -ref ref.fasta > out.aln

## For a very slow, but more memory efficent Waterman
alnSeq -use-mem-water -query query.fasta -ref ref.fasta > out.aln

# File formatting

## Output an EMBOSS like file
alnSeq -format-emboss -query query.fasta -ref ref.fasta -out out.aln

## Output a clustal file
alnSeq -format-clustal -query query.fasta -ref ref.fasta -out out.aln

## Output a fasta file
alnSeq -format-fasta -query query.fasta -ref ref.fasta -out out.aln

## Trim sequences
alnSeq -print-aligned -query query.fasta -ref ref.fasta -out out.aln

## Print positions for non EMBOSS files
alnSeq -print-positions -query query.fasta -ref ref.fasta -out out.aln
```

## Python library

## How to complie alnSeq as a python library

For python alnSeq is compiled with -DBYTEMATRIX and
  -DINSDELSNP. This install expects a gcc compile
  (minigw or cygin for windows).

You can change the compiler used with CC, but nothing else.

```
# Global install (need root permision)
make python

# Local install (requires pip)
make python local
```

### Using the python functions

The python wrapper functions take in an reference and
  query sequence and returns an aligned reference and query
  sequence. The returned item is a list, which has the
  reference sequence (ret[0]), the query sequence (ret[1]),
  and a score (ret[2]; is always 0 for the Hirschberg).

Required arguments are a reference (first argument or
  ref = sequence) and query sequence (second argument or
  query = sequence).

Optional arguments include:

- gap opening score (gapOpen = -10)
- gap extension score (gapExtend = -1)
- reference start position (refStart = 0)
- reference ending position (refEnd = length(seq) - 1)
- query start position (queryStart = 0)
- query ending position (queryEnd = length(seq) - 1)
- The path to a scoring matrix
  (scoreMatrix = /path/to/matrix.txt)
   - see scoring-matrix.txt for an example.
- The Waterman Smith alignment prints out only the aligned
  regions. This can be disabled with fullAln=True.

`alignedRef,alignedQuery = alnSeqfunction(refSeq,querySeq)`

Function names:

- alnSeqHirsch (Hirschberg)
- alnSeqNeedle (Needleman Wunsch)
- alnSeqWater (Waterman Smith)

# Explaining alnSeq

## What is alnSeq?

AlnSeq is an sequence alignment program that uses either 
  a Needleman Wunsch, Hirschberg, or Smith Waterman
  alignment algorithm. AlnSeq is unique from a traditional
  Needleman Wunsch or Smith Waterman in how it handles the
  scoring matrix and the direction matrix.

One thing I do want to point out is that there are very
  memory efficient algorithms, such as the striped Waterman
  Smith alignment, for optimal local alignments.

## The direction matrix

AlnSeq stores each direction in two bits for the Needleman
  and Waterman alignments. These bits are packed
  into an 8 bit integer. This reduces the directional
  matrix size by 4, but comes at the cost of slower speed.
  This also means that only one direction is stored,
  instead of all possible alternatives.

## The scoring matrix (Needleman/Waterman)

AlnSeq also reduces the scoring matrix down to one row,
  which holds the last scores or the previously updated
  scores. This removes the scoring matrix, but also removes
  the ability to scan the scoring matrix for other high
  scores, which is sometimes done for an Waterman Smith
  alignment. The best score is found while building the
  matrix.

AlnSeq also supports alternative alignments with
  -query-ref-scan by storing the best score for each
  reference base and each query base (Starting positions,
  ending positions, score). The score, starting reference
  position, starting query position, ending reference
  position and ending query position are then printed
  out (currently in the same file before the alignment).
  There is no filter, so this will print out everything
  that is at or above -min-score. This includes duplicate
  scores.

## My memory efficient Waterman

The more memory efficient Waterman is slow, but it also
  uses memory in linear time. It is like the alternative
  alignment step (-query-ref-scan), except that it is
  designed to work on a single directional row. It returns
  the best score, start of the alignment, and the end of
  the alignment. I then use the Hirschberg to find the
  best alignment. This means for you have to pay the time
  cost of the alternative alignment and Hirschberg steps.

# Some light benchmarking

## Benchmarking setup

For each benchmark I am using four different lengths of
  genomes. The small genome is ~1700 bases, the Mid genome
  is ~10000 bases, the large genome is roughly ~27000
  bases, and the huge genome (or ramkiller) is around
  ~199980 bases. However, I am currently removing the huge
  genomes from my graphs (in datasets), because not all
  programs can align them in under 32 Gb.

Each test size has a reference and query, which are similar
  and aligned to each other 10 times. I am also aligning
  each size to all the other genome sizes. In these cases
  I do two runs, one with the smaller size as the reference
  and the other with the larger size as the reference. Each
  test/run is replicated 10 times.

For the python test I am using Align.PairwiseAlign from
  biopython 1.82, sequence_aligner
  [https://github.com/kensho-technologies/sequence_align](
   https://github.com/kensho-technologies/sequence_align),
  and the python library for alnSeq (this repository).
  I am using pairwiseAlign from biopython because it is was
  the faster aligner in the benchmark done in the
  sequence_aligner repository. I am using sequence_aligner,
  because it is one of the more popular python Needleman
  libraries.

For the standalone test I am using EMBOSS version 6.6.0.0,
  bioalignment, ssw_test, and alnSeq.  EMBOSS was used
  because it is a large toolkit. Bioalignment was used
  because it was a popular Hirschberg, ssw_test was used
  because it is one of the faster Smith Waterman.

Time, memory, and CPU usage were recorded with gnu time
  (/usr/bin/time -f "%e\t%M\t%P"). Programs were
  benchmarked on a computer with an Intel i7-6700K CPU
  (4.00GHz) with 32 Gb of memory.

## Python library results

![Time usage in seconds to align genomes with each
  python library.
](analysis/20231029-benchmark/20231029-alnSeqPy-time.svg)

Figure 1: Time usage in seconds to align genomes with
  each python library.

We compared the time usage it took biopython's pairwise
  aligner, sequence aligner, and alnSeq to align 1700
  nucleotide to 27000 nucleotide long genomes. We found
  that alnSeq was slightly faster than biopython, which
  was often faster than sequence_aligner (Figure 1).
  However, biopython was slower than sequence_aligner for
  the small/small genome alignment (Figure 1).

![Memory usage in seconds to align genomes with each
  python library.
](analysis/20231029-benchmark/20231029-alnSeqPy-memory.svg)

Figure 2: Memory usage in seconds to align genomes with
  each python library.

We compared the memory usage of pairwise aligner from
  biopython, sequence_aligner, and alnSeq. We found that
  alnSeq and biopython used less memory than
  sequence_aligner when the same algorithms were compared
  (Needleman to Needleman or Hirschberg to Hirschbergj)
  (Figure 2). Biopython used more memory for the smaller
  alignments than alnSeq, but the memory difference became
  very small as larger genomes were aligned (Figure 2). 

As expected the Needleman aligners used more memory than
  the Hirschberg aligners (Figure 2).

![Percent CPU used to align genomes with each python
  library.
](analysis/20231029-benchmark/20231029-alnSeqPy-CPU.svg)

Figure 3: Percentage of CPU used to align genomes with each
  python library.

The CPU usage was compared for all three python libraries.
  As expected both alnSeq and sequence_aligner used 100%
  (one thread) (Figure 3). However, biopython used at least
  two threads. The impact of the extra threads used by
  biopython seemed to decrease as larger genomes were
  aligned (Figure 3).

Overall, it looks like the pairwise aligner from biopython
  was better than sequence_aligner for the non-Hirschberg 
  alignment. Biopython also has settings for gap extension
  and gap opening scores, while alnSeq only applies the
  gap extension penalties, and sequence_aligner applies
  neither. Compared to alnSeq, biopython had similar
  speed, which means that biopython would and would likely
  be as fast or faster than alnSeq if it did not use gap
  ending penalties.

What it really comes down to is this. If memory is a
  concern, then use the Hirschberg from alnSeq. However,
  if memory is not a concern, then use biopython.

I will note there could be better, programs I did not test
  on github. Also, I did not benchmark local alignments
  (both biopython and alnSeq have this option). However,
  for local aligners, I expect ssw_test to be better than
  biopython (unless it uses ssw) or alnSeq.

## Standalone benchmark

**I need to make new graphs to reflect the 20231022
  update**

I picked out three programs to compare alnSeq to. The first
  is emboss, which is a more commonly used toolkit. The
  second is bio-alignment (bio or ssw), which had a
  Hirschberg. The last was the
  Complete-Striped-Smith-Waterman-Library (bio or ssw),
  which supports a vectorized, memory efficient, striped
  Smith Wateman alignment. This is not an exhaustive list,
  nor does it include the best programs. It is just a
  simple and quick list.

Also my figures show two different runs of alnSeq. One
  run is the default compile, with the Hirschberg using
  byte arrays and the Waterman and Needleman using two
  bit arrays.
The other run (alnSeqFast) is `make fast` and has the
  gap opening penalty disabled, uses byte arrays, and 
  is hardcoded to prefer insertions, then snps, then
  deletions when scores are equal.
This is the closest you can get to have a similar
  comparison to bio-alignment.

![Memory usage of alnSeq compared to the Waterman Smith,
  Needleman Wunsch, and Hirschberg aligners from
  bio-alignment (bio-align), emobss, and
  the Complete-Striped-Smith-Waterman-Library.
  alnSeq-scan is printing alternative bases.
](
  analysis/20230827-benchmarking/20230827-alnSeq-memory.svg
)

As expected the memory usage was much lower for the
  Hirschberg aligners and striped Smith Waterman
  (local facet; bio_or_ssw), while the non-Hirschberg
  aligners for alnSeq, emboss, and bio-alginment needed
  large amounts of memory.
For the non-Hirschberg's and striped smith waterman
  alignments, alnSeq needed less memory than emboss or
  bio-alignment.
With alnSeq using two bit arrays taking the least amount
  of memory.
When compared to bio-alignments Hirschberg, alnSeq's
  Hirschberg used less memory than bio-alignments
  Hirschberg, however, the memory usage for both
  Hirschbers is small and so the memory saving has no real
  impact.

Bio-alignment's Hirschberg used more memory when two
  different alignments (small-huge, mid-huge, large-huge,
  small-large) were aligned.
From a glance at the code I suspect this was due to 
  bio-alignment allocating memory for its returned scoring
  arrays at each recursion call.
For highly different alignments this would result in
  the midpoint being closer to 1, which would result in
  an nearly identical returned scoring row size in the next
  recursion call.
Since, bio-alignment is not freeing these arrays right
  away, it is possible that these arrays would continue to
  build up at each call, which results in increased memory
  usage.
However, a local alignment should be used instead of an
  global alignment (Needleman/Hirschberg) in these cases.

We found that ssw had lower memory usage when aligning
  similar genomes and a higher memory usage when the
  genomes were very different.
I am not sure why this is happening, but it may be due to
  how much of the matrix it has to construct.
Despite this, ssw still uses less memory than alnSeq's
  Waterman alignment and so, is a better option.

![Time usage of alnSeq compared to the Waterman Smith,
    Needleman Wunsch, and Hirschberg aligners from
    , emobss, and
  Complete-Striped-Smith-Waterman-Library.
](analysis/20230827-benchmarking/20230827-alnSeq-time.svg)

For time usage we found that emboss, alnSeq scan, or
  alnSeqs two bit Hirschberg (not shown) were the slowest
  algorithms.
The fastest programs were ssw (30x faster then alnSeq and
  10x faster than alnSeqFast Waterman), followed by the
  Needleman and Hirschberg from alnSeq-fast, which uses a
  byte matrix and like bio-alignment, ignores gap opening
  penalties.
Some of the extra time needed for Emboss could be due to
  its calculating both a gap opening and gap ending
  penalty.

## Final notes

Overall alnSeq is not the best alignment program and other
  then the Striped Smith Waterman, it has not been compared
  to the more efficient programs that utilize vectors, 
  GPUS, or multiple threads. It does ok on memory usage,
  but this often comes at a steep time cost.

Overall this was a good learning project that I will use
  in the future just because I coded it. That being said,
  I am glad that this is almost over.

AlnSeq does not use decimals, so if you want decimals for
  the gap extension penalty, so you will have to multiply
  all scores by 10 to 1000.

# Thanks

- To my dad for being willing to listen and give advice.
- To the people (never met) who coded baba
  [https://baba.sourceforge.net/](https://baba.sourceforge.net/).
  It gave me a great visual on how the Needleman Wunsch
  algorithm worked. Their were many other sources, but this
  was the one that was the most useful to me.
- Wikipedia's entry about the Hirschberg. It is helping me
  understand how the Hirschberg works
- Bio-alignment coded a Hirsberg, which I used to help
  understand how the Hirsberg alignment worked. 
- To all the sights providing guidance on how aligners
  worked. There are to many to list, but each one I used
  helped me understand these algorithms a little better.
