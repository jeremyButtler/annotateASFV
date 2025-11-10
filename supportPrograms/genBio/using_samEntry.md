# Use:

Here to give an idea of how to use the samEntry library in
  genBio. This file is/will be very long, so leverage the
  table of contents to find what you want.

Functions that deal with reading and printing sam files.

Files: samEntry.c and samEntry.h

- Dependencies:
  - base10str from genLib
  - numToStr from genLib
  - ulCp from genLib
  - fileFun from genLib
  - endLine.h from genLib
  - ntTo5Bit.h from genBio

# samEntry

Has the samEntry structure, which is for reading in
  read mappings from sam files. This structure can be
  used to print sam file entries as stats (filtsam
  -out-stats), as fastq's, as sam files, or as fasta
  files.

## Error types:

| error                  | meaning                     |
|:-----------------------|:----------------------------|
| 0                      | no error or memory error    |
| def\_EOF\_samEntry     | at end of file              |
| def\_memErr\_samEntry  | memory error                |
| def\_fileErr\_samEntry | error reading input file    |
| def\_expand\_samEntry  | resized input array         |

Table: error types returned by samEntry functions.
  means return message depends on fuction.

Not all functions return all error types. Also, 0 can
  mean no errors (if error type returned) or an error if
  a value, such as a pointer is returned.

## samEntry struct:

The samEntry struct is large, but should have enough
  in it to work with sam file entries.

- read and reference ID/name variables
   - qryIdStr: is a c-string with the read id for the
     mapped read/sequence
     - qryIdLenUC holds the length
   - refIdStr: is a c-string with the read id for the
     reference the read/sequence mapped to
     - refIdLenUC holds the length
- variables with the cigar entry
  - cigTypeStr: is the symbols (=, X, M, I, D, S, H) in
    the cigar entry
  - cigArySI: has the number of times each symbol in the
    cigar entry is repeated
  - cigLenUI is the length of cigTypeStr and cigArySI
- sequence variables
  - seqStr is a c-string with the sequence
  - qStr is a c-string with the q-score entry
    - for no q-score it can be null, `\0`, or `*\0`
  - readLenUI is the length of the sequence
- reference coordinates
  - refStartUI is the first reference base in alingment
  - refEndUI is the last reference base in alingment
  - alnReadLenUI is the number of reference bases that
    were aligned (matches + snps + deletions)
- stats:
  - mapqUC has the mapping quality
  - flagUS has the flags (2nd sam file entry)
  - medianQF has the median q-score
  - meanQF has the mean q-score
  - matchCntUI has the number of matches
  - snpCntUI has the number of SNPs
  - insCntUI has the number of insertions
  - delCntUI has the number of deletions
  - maskCntUI has the number of soft masked bases
- extra entries:
  - extraStr is a c-strig the extra entries (beyond the
    first 11) in the sam file
    - it also gets the header or comment when an `@line`
      is read
   - extraLenUI is the length (number of bytes) in
     extraStr
- other entries:
  - rNextStr has the RNEXT entry
    - length is in rnextLenUC
  - pNextSI has the PNEXT entry
  - tLenSI has the TLEN entry

### samEntry struct workflow:

To use a samEntry struct you must first initalize it,
  then setup it up (allocate initial memory). After this
  it is ready to use. Finally, you will need to free your
  samEntry struct.

### samEntry initalize and setup:

You can initializing a samEntry struct with init\_samEntry
  (fun02 samEntry.c/h). This takes in a point to a
  samEntry struct to initialize. There is no return.

To setup a samEntry structure use setup\_samEntry (fun03
  samEntry.c/h). This takes in a samEntry pionter and
  returns 0 for no errors and memErr\_samEntry for memory
  errors.

You can skip these two steps by using makeSamEntry (fun06
  samEntry.c/h). It returns 0 for memory error or a
  pointer to a initialized and setup heap allocated
  samEntry struct on success.

At times you may need to blank all values in a samEntry
  struct. To do this use blank\_samEntry (fun01
  samEntry.c/h). This will lazly blank all c-strings and
  the cigar values (first values set to 0), but will fully
  clear the quality score histogram (not mentioned in
  core section because not imporant).

### samEntry freeing:

you can free samEntry structs using freeStack\_samEntry
  (fun04 samEntry.c/h) to free all internal variables or
  samEntry structs on the stack (call setup\_samEntry if
  you plan to reuse).

For heap allocated samEntry structs use freeHeap\_samEntry
  (fun05 samEntry.c/h). You must set the samEntry
  structure pointer to null/0.

### samEntry manipulation:

#### samEntry: get sam file entry

You can use get\_samEntry (fun12 samEntry.c/h) to get a
  sam file entry from a file. This returns 0 for success,
  def\_EOF\_samEntry for end of file,
  def\_fileErr\_samEntry for an invalid line, and
   def\_memErr\_samEntry for memory errors.

- Input:
  - samEntry structure pionter to get new entry
    - this function calls blank\_samEntry, so do not worry
      about this
  - FILE pointer with sam file to get entry from
    - advanced to next line in sam file

Be aware, comments go in extra entry.

```
errorSC = get_samEntry(samEntryStructPointer, samFILE);
if(errorSC == def_memErr_samEntry)
   /*deal with memory error*/
else if(errorSC == def_fileErr_samEntry)
   /*deal with file error*/
else if(errorSC == def_EOF_samEntry)
{ /*Else If: at end of file*/
   freeHeap_samEntry(samEntryStructPointer);
   return 1;
} /*Else If: at end of file*/

if(samEntryStructPointer->extraStr[0] == '@')
   /*if is comment (header line)*/;
else
   /*else is read*/
```

You can also read in a line and convert it to a sam file
  entry using lineTo\_samEntry (fun11 samEntry.c/h). The
  input is a pointer to a samEntry struct and the c-string
  to convert. The return values are def\_memErr\_samEntry
  for memory errors or 0 for success.

```
lengthSL =
   getFullLine_fileFun(
      samFILE,
      &bufferStr,
      &bufferSizeSL,
      &bytesSL,  /*real bytes read from file*/
      0          /*copy at index 0 in buffStr*/
   );

if(lineTo_samEntry(samStructPionter, bufferStr))
   /*deal with memory error*/
```

#### samEntry: move to target position

You can find the sequence position of a reference position
  in a samEntry using findRefPos\_samEntry (fun13
  samEntry.c/h). This returns the index of the last
  cigar entry before the current one. It also modifies
  the input to be on the correct base.

- Input:
  - samEntry structure pointer with sequence to find
    refernce position
  - int pointer to cigar entry currently on
  - int pointer to number of bases left in the cigar entry
    (what have I not moved past)
  - reference position to find (target)
    - can be forwards or backwards
  - int pointer to the current reference position
    - this is updated to the target position or the last
      reference base in the sequence
  - int poiter to the current sequence position
    - this is updated to the target position or the last
      base

Here is an example to get base 100 in the query sequence.
  Just remember to check for deletions (missing in query).

```
signed int sequencePosition = 0;
signed int referencePosition = samEntryStruct->refStartUI;
signed int cigarPosition = 0;
signed int basesInCigar = samEntryStruct->cigArySI[0];
signed int target = 100;
signed int lastCigar = 0;

lastCigar =
   findRefPos_samEntry(
     samEntryStruct,
     &cigarPosition,
     &basesInCigar,
     target,
     &referencePosition,
     &sequencePosition,
);

if(samEntryStruct->cigTypeStr[cigarPosition] == 'D')
   /*reference base at target is a deltion in the query*/;
else
  /*do somenthing with position*/

/*now find position 1000 from current position*/
lastCigar =
   findRefPos_samEntry(
     samEntryStruct,
     &cigarPosition,
     &basesInCigar,
     1000,
     &referencePosition,
     &sequencePosition,
);
```

#### samEntry: swap entries

You can swap the entries in two samEntry structures using
  swap\_samEntry (fun14 samEntry.c/h). This takes pointers
  to two samEntry structures to sawp and returns nothing.

Swapping is done lazely, so pointers are swapped instead
  of copying the sequence, q-score entry, cigar arrays,
  and the extra entry, pointers are swapped. The query id,
  reference id, RNEXT, and qHistUI are not pointers, and
  so are copied.

#### samEntry: less usefull functions

**q-score update**

If you modify the q-score entry in a samEntry struct, then
  you can update thm with findQScores\_samEntry (fun08
  samEntry.c/h). This takes in a pointer to a samEntry
  struct to update and returns nothing. It will blank and
  then update all q-score values in the structure.

**q-score copy**

If you need to copy a q-score entry into you samEntry
  structure, use cpQEntry\_samEntry (fun09 samEntry.c/h).
  This function will both copy the q-score entry and
  also update all q-score stats. You must make sure there
  is enough memory before updating the q-scores.

- Input:
  - samEntry structure pointer to copy q-scores to
  - c-string with q-scores to copy
  - 1 (always provide this)

The return value is the length of the copied c-string.

**get sam file header**

You can read the entire header for a samEntry struct
  into a c-string using getHead\_samEntry
  (fun20 samEntry.c/h).

- Input:
  - samEntry structure pointer to read sam file with
    - will be set to frist non-header line
  - FILE poitner to sam file to get header from
- Returns:
  - pointer to c-string with header
  - 0 for memory errors

### samEntry printing:

#### samEntry printing; sam files

You can print a samEntry struct to a sam file using
  p\_samEntry (fun15 samEntry.c/h).

- Input:
  - samEntry structure pointer to print to sam file
  - 1 or 0 marking if you want a line break at the end
    - 0: put line break at end
    - 1: do not add a line break
    - Note: line breaks will vary by OS (`\n`, `\r`,
      `\n\r`, or `\r\n`).
  - file to write the samEntry to

#### samEntry printing; fastx files

You can also print to a fastq file using pfq\_samEntry (
  fun16 samEntry.c/h) or to a fasta file with
  pfa\_samEntry (fun17 samEntry.c/h). Both fuctions take
  a samEntry structure to print and a file to save the
  fastx entry to. There is no return value.

#### samEntry printing; stats

Finally, you can print out the stats in a spread sheet
  format (filtsam -out-stats) with pstats\_samEntry. There
  is not return value.

- Input:
  - samEntry structure pointer with entry to print
  - signed char pointer telling if to print the header
    - 1 print the header (set to 0 afterwards)
    - 0 do not print the header
  - 0 or 1:
    - 1: print counts for anonymous bases
    - 0: do not print counts for anonymous bases
  - FILE pointer with file to print to

#### samEntry example

Here is a quick example of using a samEntry struct. In
  this example I am reading each entry and printing out
  the stats for it (like filtsam does).

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #inlcude <stdlib.h>
#endif

#include <stdio.h>
#include "../genLib/endLine.h" /*gives me str_endLine*/
#include "../genBio/samEntry.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char headBl = 1; /*print header for stats file*/
   signed char errorSC = 0;
   struct samEntry samStruct;

   singed char *fileStr = (signed char *) "file.sam";
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ setup and open file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_samEntry(&samStruct);

   if(setup_samEntry(&samStruct))
      goto memErr_main;

   samFILE = fopen((char *) fileStr, "r");
   if(! samFILE)
      goto invalidFile_main;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ get past headers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = get_samEntry(&samStruct, samFILE);

   if(errorSC == def_memErr_samEntry)
      goto memErr_main;
   else if(errorSC == def_fileErr_samEntry)
      goto fileErr_main;
   else if(errorSC)
      goto empytFile_main;

   while(! errorSC)
   { /*Loop: get headers*/
      if(samStruct.extraStr[0] != '@')
         break; /*end of header*/

      errorSC = get_samEntry(&samStruct, samFILE);
   } /*Loop: get headers*/

   if(errorSC == def_memErr_samEntry)
      goto memErr_main;
   else if(errorSC == def_fileErr_samEntry)
      goto fileErr_main;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ print stats for each read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errorSC)
   { /*Loop: get reads in sam file*/
      if(samStruct.flagUS & 4)
         continue; /*umapped read*/
      else if(samStruct.seqStr[0] == '\0')
         continue; /*no sequence*/
      else if(samStruct.seqStr[0] == '*')
         continue; /*no sequence (version 2)*/

      pstats_samEntry(&samStruct, &headBl, 1, stdout);
           /*1 is to get anonymous base counts*/
           /*headBl is set to 0 after frist call*/

      errorSC = get_samEntry(&samStruct, samFILE);
   } /*Loop: get reads in sam file*/

   if(errorSC == def_memErr_samEntry)
      goto memErr_main;
   if(errorSC == def_fileErr_samEntry)
      goto fileErr_main;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ clean up, deal with errors, and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto cleanUp_main;

   memErr_main:;
      errorSC = 1;
      fprintf(stderr, "memory error%s", str_endLine);
      goto cleanUp_main;

   invalidFile_main:;
      errorSC = 2;
      fprintf(
         stderr,
         "could not open %s%s",
         fileStr,
         str_endLine
      );
      goto cleanUp_main;

   fileErr_main:;
      errorSC = 3;
      fprintf(
         stderr,
         "%s is not a sam file%s",
         fileStr,
         str_endLine
      );
      goto cleanUp_main;

   empytFile_main:;
      errorSC = 4;
      fprintf(stderr, "empty file%s", str_endLine);
      goto cleanUp_main;

   cleanUp_main:;
      if(samFILE)
         fclose(samFILE);
      samFILE = 0;

      freeStack_samEntry(&samStruct);

      return errorSC;
}
```
