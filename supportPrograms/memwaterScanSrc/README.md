# Use:

Does a waterman alignment between two sequences and
  returns the coordiantes and scores of the best
  alignments.

# Install

## Unix:

```
if [ ! -d "${HOME}/Downloads/bioTools" ];
then
   git \
       clone https://github.com/jeremybuttler/bioTools \
       ~/Downloads/bioTools;
fi;

cd ~/Downloads/bioTools/memwaterScanSrc;
make -f mkfile.unix;
sudo make -f mkfile.unix install;
```

# Use:

You can get the help message with `memwaterScan -h`. For
  basic usage, inpute a reference with
  the `-ref reference.fasta` and the query with
  the `-qry query.fasta` command.

Example:

```
memwaterScan -ref reference.fasta -qry query.fasta > coordinates.tsv
```

# System

The scan part of memwaterScan means it keeps the
  coordiantes and score of the best mapping alignment for
  each base in the reference sequence and the query
  sequence. These kept alignment coordiantes are then
  filtered to remove low scoring alignments and alignments
  that are nested (last alignment is the last stage of the
  current alignment).

The water in memwaterScan means it is doing a Smith
  Waterman alignment.

The mem part of memwaterScan means that it is only keeping
  track of the score, direction, start, and end of one
  row. So, it has linear memory usage, but n*m time usage.
  It also means it can not output an alignment, but only
  coordiantes and a score.
