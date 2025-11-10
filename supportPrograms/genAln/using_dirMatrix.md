# Goal:

Go over how to use the dirMatrix library to convert
  Needleman and Waterman matrix's into sam file
  alignments. The dirMatrix structure is also used in
  the Needleman and Waterman alignments.

- Files: dirMatrix.h and dirMatrix.c

- Dependencies:
  - alnSet
  - indexToCoord
  - alnDefs.h
  - samEntry from genBio
  - seqST from genBio
  - ulCp from genLib
  - charCp from genLib
  - numToStr from genLib
- Hidden Depdenencies
  - ntTo5Bit.h from genBio
  - base10str from genLib
  - fileFun from genLib
  - endLine.h from genLib

You should be familar with the alnSet (see
  using\_alnSet.md) library to hold alignment settings.
  Also, the seqST (see using\_seqST.md) to hold fasta
  and fastq sequences. And the samEntry (see
  using\_samEntry.md) library to hold the returned
  alignment.

This library is intended to be a part of the Needleman
  and Waterman libraries (see using\_needleman.md and
  using\_waterman.md).

# Overview

The dirMatrix library includes the dirMatrix structure and
  also has the functions to convert the dirMatrix
  structure into an alignment.

The dirMatrix structure is used will the Needleman and
  Waterman (not memWater) pairwise alignment libraries.

# dirMatrix structure

The dirMatrix structure holds the direction matrix and
  score row for a Needleman or Waterman alignment. It also
  has variables to record stats of the alignment.

- Variables:
  1. dirMatrixSC: directional matrix for the alignment
     - the different directions can be found in alnDef.h
     - the size of the matrix is stored in lenMatrixSL
  2. scoreSL: the score (or for Waterman; bestScore) of
     the alignment
  3. indexSL: is the index of the best scoring alignment
     (Waterman only)
  4. scoreArySL: has one row (the final row) of scores
     - the size is stored in lenScoreSL
  5. refLenSL: length of the reference genome
  6. refOffsetSL: first base aligned in reference
  7. refEndSL: last base aligned in reference
  8. qryLenSL: length of the query genome
  9. qryOffsetSL: first base aligned in query
  10. qryEndSL: last base aligned in query
  11. errSC: holds any errors from the Waterman or
      Needleman alignment
  
## setting up and freeing

The setup for a dirMatrix invovles making/declaring the
  dirMatrix structure and then initializing it. The
  Needleman and Waterman pairwise alingers will allocate
  the memory. When finished you shoudl free the dirMatrix
  with the free functions.

You can initialize the dirMatrix structure with
  the `init_dirMatrix()` function.

- Input:
  1. dirMatrix structure pointer to initialize

You can then free the dirMatrix structure using
  the `freeStack_dirMatrix()` function to free the
  varaibles in a dirMatrix structure. You can free a heap
  allocated dirMatrix structure
  with `freeHeap_dirMatrix()`, but you have to remember to
  set your dirMatrix structure pointer to 0/null
  afterwards.

- Input (freeStack and freeHeap):
  1. dirMatrix structure pointer

**Example 1; stack**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "dirMatrix.h"

int
main(
){
   struct dirMatrix matrixStackST;
   init_dirMatrix(&matrixStackST);

   /*code to use matrix goes here*/

   freeStack_dirMatrix(&matrixStackST);
} /*main*/
```

**Example 2; heap**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "dirMatrix.h"

int
main(
){
   struct dirMatrix *matrixHeapST = 0;
   matrixHeapST = malloc(sizeof(struct dirMatrix));

   if(! matrixHeapST)
      return 1; /*memory error*/

   init_dirMatrix(matrixHeapST);

   /*code to use matrix goes here*/

   freeHeap_dirMatrix(matrixHeapST);
   matrixHeapST = 0;
} /*main*/
```

# functions:

There are two functions that you can use. The first will
  convert an alignment into a sam file entry. The second
  will convert an alignment to a cigar.

## get sam file alignment

You can use the `getAln_dirMatrix()` function to convert
  an alignment from a Needleman or Waterman into a
  samEntry structure with the alignment.

- Input:
  1. dirMatrix structure pointer with the pairwise
     alignment
  2. index of the aligment or you can use 0 to use the
     index in the dirMatrix struct (input 1)
  3. 1 if the reference sequence is reverse complement to
       the read
     - use 0 if foward
  4. seqST structure pointer with the query sequence
  5. seqST structure pointer with the reference sequence
  6. samEntry structure pointer to store the alignment in
  7. unsigned int pionter to get the number of anonymous
     bases in the alignment
  8. alnSet structure pointer with the alignment settings
- Output:
  1. Modifies: input 7 to have the number of anonymous
     bases
  2. Modifies: input 6 to have the alignment
  3. Returns:
     - 0 for no errors
     - 1 for memory errors

**Example**

Would be command to compile the example if it was saved
  in a file name `tmp.c`.

```
gcc \
    -o tmp.o \
    tmp.c \
    needle.c \
    dirMatrix.c \
    alnSet.c \
    ../genBio/seqST.c \
    ../genBio/samEntry.c \
    ../genLib/ulCp.c \
    ../genLib/base10str.c \
    ../genLib/numToStr.c \
    ../genLib/charCp.c \
    ../genLib/fileFun.c \
    indexToCoord.c
```

Code for the example

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>

#include "needle.h"
#include "alnSet.h"
#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"

#include "dirMatrix.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int errSI = 0;
   unsigned int anonymousUI = 0;

   struct seqST queryStackST;
   struct seqST refStackST;

   struct alnSet settingsStackST;
   struct dirMatrix matrixStackST;
   struct samEntry samStackST;

   FILE *inFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_samEntry(&samStackST);
   init_seqST(&queryStackST);
   init_seqST(&refStackST);

   init_alnSet(&settingsStackST);
   init_dirMatrix(&matrixStackST);

   if(setup_samEntry(&samStackST))
      goto memErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: get sequences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*get the reference sequence*/
   inFILE = fopen("reference.fasta", "r");

   if(! inFILE)
      goto refFileErr_main_sec05;
   if( getFa_seqST(inFILE, &refStackST) )
      goto refFileErr_main_sec05;

   seqToIndex_alnSet(refStackST.seqStr);
      /*need to convert to index for alignment*/

   fclose(inFILE);
   inFILE = 0;


   /*get the query sequence*/
   inFILE = fopen("query.fasta", "r");

   if(! inFILE)
      goto queryFileErr_main_sec05;
   if( getFa_seqST(inFILE, &queryStackST) )
      goto queryFileErr_main_sec05;

   seqToIndex_alnSet(queryStackST.seqStr);
      /*need to convert to index for alignment*/

   fclose(inFILE);
   inFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: get alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    needle(
       &queryStackST,
       &refStackST,
       &matrixStackST,
       &settingsStackST
    ); /*fill the direction matrix*/

   if(matrixStackST.errSC)
      goto alignmentErr_main_sec05;

   /*convert the filled matrix into an alignment*/
   if(
      getAln_dirMatrix(
         &matrixStackST,
         0,             /*use index in matrixStackST*/
         0,             /*foward alignment*/
         &queryStackST, /*has query sequence*/
         &refStackST,   /*has reference sequence*/
         &samStackST,   /*gets the alignment*/
         &anonymousUI,  /*gets number anonymous bases*/
         &settingsStackST /*alignment settings*/
      )
   ) goto alignmentErr_main_sec05; /*memory error*/

   /*print out alignment*/
   p_samEntry(&samStackST, 0, stdout);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSI = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      fprintf(stderr, "memory error\n");
      errSI = 1;
      goto ret_main_sec05;

   refFileErr_main_sec05:;
      fprintf(
          stderr,
          "error getting reference sequence\n"
      );
      errSI = 2;
      goto ret_main_sec05;

   queryFileErr_main_sec05:;
      fprintf(stderr, "error getting query sequence\n");
      errSI = 2;
      goto ret_main_sec05;

   alignmentErr_main_sec05:;
      fprintf(stderr, "error aliging sequences\n");
      errSI = 3;
      goto ret_main_sec05;

   ret_main_sec05:;
      freeStack_samEntry(&samStackST);
      freeStack_seqST(&queryStackST);
      freeStack_seqST(&refStackST);

      freeStack_alnSet(&settingsStackST);
      freeStack_dirMatrix(&matrixStackST);

      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else
         fclose(inFILE);
      inFILE = 0;

      return errSI;
} /*main*/
```

## get cigar

You can use the `getCig_dirMatrix()` function to convert
  an alignment from a Needleman or Waterman into a
  cigar.

- Input:
  1. dirMatrix structure pointer with the pairwise
     alignment
  2. index of the aligment or you can use 0 to use the
     index in the dirMatrix struct (input 1)
  3. 1 if the reference sequence is reverse complement to
       the read
     - use 0 if foward
  4. seqST structure pointer with the query sequence
  5. seqST structure pointer with the reference sequence
  6. c-string pointer to get cigar letters (mutations)
     - gets "=XDIS" symbol
  7. unsigned int pionter to get the number of times each
     symbol repeats
  8. current position in the cigar
     - this is in case you are adding to a cigar
  9. size of cigar buffer
     - this is for resizing
  10. unsigned int pointer to get the first reference
      base in the cigar
  11. unsigned int pionter to get the number of anonymous
     bases in the alignment
  12. unsigned int pionter to get the number of masked
     bases (at ends; for Waterman)
  13. tells if to add deletions at end if reference is a
      bit sort
      - use 0 (disable)
  14. alnSet structure pointer with the alignment settings
- Output:
  1. Modifies: input 6 to have the cigar symbols
  2. Modifies: input 7 to have the cigar lengths
  3. Modifies: input 9 to have the new cigar size if is
     resized
  4. Modifies: input 10 to have first reference base in
     alignment
  5. Modifies: input 11 to have the number of anonymous
     bases
  6. Modifies: input 12 to have the number of masked bases
  7. Returns:
     - end of the cigar (length - 1; so is index 0)
     - 1 for memory errors

**Example**

Would be command to compile the example if it was saved
  in a file name `tmp.c`.

```
gcc \
    -o tmp.o \
    tmp.c \
    needle.c \
    dirMatrix.c \
    alnSet.c \
    ../genBio/seqST.c \
    ../genBio/samEntry.c \
    ../genLib/ulCp.c \
    ../genLib/base10str.c \
    ../genLib/numToStr.c \
    ../genLib/charCp.c \
    ../genLib/fileFun.c \
    indexToCoord.c
```

Code for the example

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "needle.h"
#include "alnSet.h"
#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"

#include "dirMatrix.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int errSI = 0;
   signed long lenSL = 0;

   signed char *cigHeapStr = 0;
   signed int *cigCountHeapSI = 0;
   unsigned int cigSizeUI = 0;
   unsigned int firstRefBaseUI = 0;

   unsigned int anonymousUI = 0;
   unsigned int maskUI = 0;

   struct seqST queryStackST;
   struct seqST refStackST;

   struct alnSet settingsStackST;
   struct dirMatrix matrixStackST;

   FILE *inFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_seqST(&queryStackST);
   init_seqST(&refStackST);

   init_alnSet(&settingsStackST);
   init_dirMatrix(&matrixStackST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: get sequences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*get the reference sequence*/
   inFILE = fopen("reference.fasta", "r");

   if(! inFILE)
      goto refFileErr_main_sec05;
   if( getFa_seqST(inFILE, &refStackST) )
      goto refFileErr_main_sec05;

   seqToIndex_alnSet(refStackST.seqStr);
      /*need to convert to index for alignment*/

   fclose(inFILE);
   inFILE = 0;


   /*get the query sequence*/
   inFILE = fopen("query.fasta", "r");

   if(! inFILE)
      goto queryFileErr_main_sec05;
   if( getFa_seqST(inFILE, &queryStackST) )
      goto queryFileErr_main_sec05;

   seqToIndex_alnSet(queryStackST.seqStr);
      /*need to convert to index for alignment*/

   fclose(inFILE);
   inFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: get alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    needle(
       &queryStackST,
       &refStackST,
       &matrixStackST,
       &settingsStackST
    ); /*fill the direction matrix*/

   if(matrixStackST.errSC)
      goto alignmentErr_main_sec05;

   /*convert the filled matrix into an alignment*/
   lenSL =
      getCig_dirMatrix(
         &matrixStackST,
         0,             /*use index in matrixStackST*/
         0,             /*foward alignment*/
         &queryStackST, /*has query sequence*/
         &refStackST,   /*has reference sequence*/
         &cigHeapStr,   /*gets cigar symbols*/
         &cigCountHeapSI, /*gets number times each symbol
                          ` is repeated
                          */
         0,              /*cigar is empty*/
         &cigSizeUI,
         &firstRefBaseUI, /*first ref base in alignment*/
         &anonymousUI,  /*gets number anonymous bases*/
         &maskUI,         /*number of masked bases*/
         0,               /*do not add deletions*/
         &settingsStackST /*alignment settings*/
      );

    if(lenSL < 0)
       goto memErr_main_sec05; /*memory error*/

    printf(
       "alignment starts at reference position %u\n",
       firstRefBaseUI
    );

    printf(
      "  cigar: %i%c",
      cigCountHeapSI[0],
      cigHeapStr[0]
   );

   for(errSI = 1; errSI <= lenSL; ++errSI)
   { /*Loop: print out the cigar*/
      printf(
         "%i%c",
         cigCountHeapSI[errSI],
         cigHeapStr[errSI]
      );
   } /*Loop: print out the cigar*/

   printf("\n");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSI = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      fprintf(stderr, "memory error\n");
      errSI = 1;
      goto ret_main_sec05;

   refFileErr_main_sec05:;
      fprintf(
          stderr,
          "error getting reference sequence\n"
      );
      errSI = 2;
      goto ret_main_sec05;

   queryFileErr_main_sec05:;
      fprintf(stderr, "error getting query sequence\n");
      errSI = 2;
      goto ret_main_sec05;

   alignmentErr_main_sec05:;
      fprintf(stderr, "error aliging sequences\n");
      errSI = 3;
      goto ret_main_sec05;

   ret_main_sec05:;
      if(cigHeapStr)
         free(cigHeapStr);
      cigHeapStr = 0;

      if(cigCountHeapSI)
         free(cigCountHeapSI);
      cigCountHeapSI = 0;

      freeStack_seqST(&queryStackST);
      freeStack_seqST(&refStackST);

      freeStack_alnSet(&settingsStackST);
      freeStack_dirMatrix(&matrixStackST);

      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else
         fclose(inFILE);
      inFILE = 0;

      return errSI;
} /*main*/
```
