# Use:

Here to give an idea of how to use the seqST libraries in
  genBio.

Functions and structure to extract sequences from fastx
  files and reverse complement. Can not handle .gz files,
  use gzSeqST for that.

Files: seqST.c and seqST.h

- Depencencies:
  - ulCp from genLib
  - fileFun from genLib

# seqST

seqST is a library designed to read in sequences from
  fasta and fastq files. It also includes a reverse
  complement function.

For reading `.gz` compressed files, use the seqST
  structure in seqST with `gzSeqST`. The gz size of
  sequence reading is separate because adding gz support
  adds roughly 18kb to the musl static binary size.

## seqST structure

The seqST structure is designed to hold sequence data
  from fasta or fastq files. The variables you migth work
  with are:

- idStr holds the read id
  - the read id length is stored in idLenSL
- seqStr holds the sequence
  - the sequence length is stored in seqLenSL
- qStr (fastq only) holds the quality score entry
  - the q-score length is stored in qLenSL
  - the first character set to 0/null if no q-score entry
- offsetSL is the index of first base to align
  - for pairwise aligners (needle, water, memwater)
  - index 0
- endAlnSL is the index of last base to align
  - for pairwise aligners (needle, water, memwater)
  - index 0

## seqST functions

The workflow for a seqST structure is to initialize, then
  read sequences. You can then manipulate the seqST
  structer how you want.

### seqST setup

To initialize a seqST structure use init\_seqST (fun07).
  The input is a pointer to a seqST structure to
  initialize.

### seqST get sequences

You can read a fasta sequence using getFa\_seqST (fun03).
  The input is a pointer to a seqST structure to hold the
  sequence and a FILE poitner to the file to read from.
  The output is 0 for no errors, def\_EOF\_seqST if at end
  of file (no sequence read in), def\_fileErr\_seqST for
  file errors,
  and `def\_fileErr\_seqST | def\_badLine\_seqST` for
  non-fasta entries.

You can read a fastq sequence using getFq\_seqST (fun03).
  The input is a pointer to a seqST structure to hold the
  sequence and a FILE poitner to the file to read from.
  The output is 0 for no errors, def\_EOF\_seqST if at end
  of file (no sequence read in), def\_fileErr\_seqST for
  file errors,
  and `def\_fileErr\_seqST | def\_badLine\_seqST` for
  non-fasta entries.

### seqST reverse complement

You can reverse complement a sequence using
  revComp\_seqST.  The input is a pionter to a seqST
  structure to reverse complement.

### seqST read fasta file

You can get every sequence from a fasta file using
  readFaFile\_seqST. The output is a pointer to a seqST
  structure array (sorted by id) or 0 for any kind of
  error.

- Input
  - FILE pointer to fasta file to get sequences from
  - signed long pointer get number of sequences read
  - signed long pointer get number of seqST structs in
    the returned array
  - signed char pointer to get any errors
    - set to def\_memErr\_seqST for memory errors
    - set to def\_fileErr\_seqST for file errors
    - set to `def\_fileErr\_seqST | def\_badLine\_seqST`
      a invalid entry

You can then find sequences by name/id using
  search\_seqST.  The input is the sorted seqST structure
  array to search, the name/id to search for (c-string),
  and the number of sequences in the array.

When finished with the seqST structure array, you can
  free it using freeHeapAry\_seqST. You will need to
  provide the array and the number of seqST structers in
  the array.

## seqST example

Here is a quick example the reverse complements every
  read in a fastq file. Not very usefull, but gives an
  idea of the workflow.

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "seqST.h"

int
main(
){
   signed char errorSC = 0;
   struct seqST sequenceStruct;
   FILE *fastqFILE = fopen("file.fastq", "r");

   init_seqST(sequenceStruct);

   if(! fastqFILE)
      /*HANDEL FILE ERRORS*/

   while(! errorSC)
   { /*Loop: read in fastq sequences*/
      errorSC = getFq_seqST(fastqFILE, &sequenceStruct);
      if(errorSC)
         break;
      revComp_seqST(&sequenceStruct);
      fprintf(
          stdout,
          "@%s\n%s\n+\n%s",
          sequenceStruct.idStr,
          sequenceStruct.seqStr,
          sequenceStruct.qStr
      );
   } /*Loop: read in fastq sequences*/

   if(errorSC == def_EOF_seqST)
      ;
   else if(errorSC != def_memErr_seqST)
      /*handel memory errors*/
   else
      /*handel file errors*/

   freeStack_seqST(&sequenceStruct);
   fclose(fastqFILE);
   fastqFILE = 0;
   return 0;
}
```
