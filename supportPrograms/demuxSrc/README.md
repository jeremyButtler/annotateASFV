# Use:

Find amplicons in sequences, such as reference genomes or
  reads.

Note: demux will not work well for cases were two primers
  could be nested in each other, such as the spoligotyping
  primer tuberculosis. In these cases nothing will be
  printed out (except to the log).

The origin of the name came from my thought of using it
  as a barcode demuxer, but I found out that the first
  100 bases in a Minion sequenced read are to noisy.
  Adapters are always missed and barcodes need low score
  thresholds. Use `dorado demux`. The reason for the
  attempt was me demuxing with the wrong kit.

# License:

The primary is the Unlicense, but if that does not work
  treat this as under the MIT license.

# Install:

## Unix (Linux/Mac/BSD):

```
cd "${HOME}/Downloads;

if [ ! -d bioTools ]; then
   git clone https://github.com/jeremyButtler/bioTools bioTools;
fi;

cd bioTools/demuxSrc;
make -f mkfile.unix;
sudo make -f mkfile.unix install;
```

# Usage:

You can print the help message with `demux -h`.

## amplicon trimming

You will need to provide a tsv file with `-prim-tsv`. The
  first row is the header and is ignored. The next rows
  should have the primers.

- TSV format:
  - Column 1: name of primers (no spaces or tabs)
  - Column 2:
    - Yes, Y, True, T, or 1 if both primers need to be
      present
    - No, N, False, F, or 0 if primers are not paired
      (looking for single primer entries)
      - For single primers, demux will not print cases
        were a forward primer is followed by a reverse,
        because it is hard to tell if these primers were
        paired or some other problem
  - Column 3: forward primer sequence
  - Column 4: reverse primer sequence

Here is an example of a tsv for two primers. In this case
  I am using spaces instead of tabs for column separators
  (demux supports tabs, spaces, and a mix of tabs and
   spaces).

```
id    paired  for_seq               rev_seq
katG  T       ACCCACCCATTACAGAAACC  AGCAGGTTCACGAAGAAGTC
embA  T       CATCCTCACCGCCCTTA     CAACCTGTGGCTTCTTCTC
```

Other then the table you will need to provide one or more
  fastx (fasta or fastq) files at the end. The file(s) can
  be gz compressed (fasta.gz or fastq.gz).

```
demux -prim-tsv primers.tsv reads.fasq > trimmed-reads.fq;
```

The default output is a log and the trimmed reads to
  stdout (terminal). You can rename the log with
  the `-prefix <your_name>` command.

```
demux -prim-tsv primers.tsv reads.fasq -prefix reads-log \
  > trimmed-reads.fq;
```

## gene detection

You can search for gene coordinates using the `-gene`
  flag. The system searches for one gene at a time so that
  the window size is correct for each gene. The default
  settings are slow (setup for primers).

```
demux \
    -gene genes.fasta \
    reads.fastq \
  >  gene-coordinates.tsv;
```

Since genes are at least 100 base pairs, I would recomend
  setting the kmer length (`-kmer-len`) to 7 and the
  percentage of kmers needed to do an alignment
  (`-kmer-perc`) to 60%. This provides a good speed boost.

```
demux \
    -gene genes.fasta \
    -kmer-len 7 \
    -kmer-perc 0.6 \
    assembly.fasta \
  >  gene-coordinates.tsv;
```
