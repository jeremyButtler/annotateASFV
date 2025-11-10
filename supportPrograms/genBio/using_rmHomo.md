# Goal

Describe how to use the rmHomo library to deal with
  homopolymers. This currently is homopolymer indel
  removal or homopolymer left align alignment.

Files: rmHomo.c rmHomo.h

- Dependenices:
  - samEntry
  - ulCp from genLib
- Hidden depencencies
  - fileFun from genLib
  - base10str from genLib
  - numToStr from genLib
  - ulCp from genLib
  - ntTo5Bit.h
  - endLine.h from genLib

You should be familar with the samEntry library (see
  using\_samEntry.md). It is not needed, but will help to
  be familar with the seqST library (see using\_seqST.md).

# Overview

The rmHomo library has two functions currently. The first
  is to remove index from homopolymers. The second is to
  left align indels in homopolymers.

## indel removal

Indel removal in homopolymers for reads in sam files is
  done using the `indel_rmHomo()` function.

- Input:
  1. samEntry structure pointer with read/sequence to
     remove homopolymers from
  2. c-string with sequence of reference the read/sequence
     (input 1) was mapped to
  3. minimum length of the homopolymer to remove indels
     from
  4. maximum indel size to remove
     - idea is large indels are very unlikely
  5. base to mask deletions with
     - if 0, then the reference base is used
  6. 1: to scan neighbor bases to see if indel is next to
     a homopolymer (0 do nothing)
  7. c-string buffer to use while removing indels
     - for the sequence
  8. (signed int \*) current size of sequence c-string
     buffer
  9. c-string buffer to use while removing indels
     - for the q-score entry
  10. (signed int \*) current size of quality score
      c-string buffer
- Output:
  - Modifies samEntry (input 1) to have indels in
    homopolymers removed
  - Modifies: input 7 and 9 to be resized as needed
  - Modifies: input 8 and 10 to have new buffer sizes if
    they were resized
  - Returns:
    - 0 for no errors
    - 1 for memory errors

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>
#include "rmHomo.h"
#include "samEntry.h"
#include "seqST.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: declare variables and initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;

   /*buffers for indel clean up step*/
   signed char *sequenceBufferStr = 0;
   unsigned int sequenceSizeUI = 0;
   signed char *qScoreBufferStr = 0;
   unsigned int qScoreSizeUI = 0;

   struct samEntry samStackST; /*gets reads*/
   struct seqST refStackST;    /*gets reference sequence*/

   FILE *inFILE = 0;

   /*initialize*/
   init_samST(&samStackST);
   init_seqST(&refStackST);

   if( setup_samST(&samStackST) )
      goto memoryErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: open files and get reference/first read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   inFILE = fopen("reference.fasta", "r");

   if(! inFILE)
      goto fileErr_main_sec04;

   errorSC = getFa_seqST(inFILE, &refStackST);
   if(! errorSC)
      ;
   else if(errorSC = def_memErr_seqST)
      goto memoryErr_main_sec04;
   else
      goto fileErr_main_sec04;

   fclose(inFILE);
   inFILE = 0;

   inFILE = fopen("reads.sam", "r");
   if(! inFILE)
      goto fileErr_main_sec04;
   errorSC = get_samEntry(&samStackST, inFILE);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: remove indels from reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errorSC)
   { /*Loop: read in sam file*/
      if(samStackST.extraStr[0] == '\0')
         p_samEntr(&samStackST, 0, stdout); /*header*/
      else if(samStackST.flagUS & 4)
         goto nextEntry_main_sec02; /*unmapped read*/
      else if(! samStackST.seqStr[0])
         goto nextEntry_main_sec02; /*no sequence*/
      else if(! samStackST.seqStr[0] == '*')
         goto nextEntry_main_sec02; /*no sequence*/

      errorSC =
         indel_rmHomo(
            &samStackST,
            refStackST.seqStr,
            3,    /*homopolymer is at least 3 bases long*/
            5,    /*indel is less then 5 bases long*/
            'N',  /*mask removed deletions with N*/
            0,    /*do not scan neighbor bases for
                    homopolymers*/
            &sequenceBufferStr, /*do not mess with*/
            &sequenceSizeUI,    /*do not mess with*/
            &qScoreBufferStr,   /*do not mess with*/
            &qScoreSizeUI       /*do not mess with*/
         );

      if(errorSC)
         goto memoryErr_main_sec04;

      nextEntry_main_sec02:;
         errorSC = get_samEntry(&samStackST, inFILE);
   } /*Loop: read in sam file*/

   if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else if(errorSC != def_EOF_samEntry)
      goto fileErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto return_main_sec04;

   memoryErr_main_sec04:;
      fprintf(stderr, "memory error\n");
      errorSC = 1;
      goto return_main_sec04;

   fileErr_main_sec04:;
      fprintf(
         stderr,
         "error in reference.fasta or reads.sam\n"
      );
      errorSC = 2;
      goto return_main_sec04;

   return_main_sec04:;
      freeStack_samST(&samStackST);
      freeStack_seqST(&refStackST);

      if(sequenceBuffStr)
         free(sequenceBufferStr);
      sequenceBufferStr = 0;

      if(qSocreBuffStr)
         free(qSocreBufferStr);
      qSocreBufferStr = 0;

      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else if(inFILE == stdout) ;
      else if(inFILE == stderr) ;
      else
         fclose(inFILE);
      inFILE = 0;

      return errorSC;
} /*main*/
```

## left alignment

The `leftAlnIndel_rmHomo` will left align deletions in
  homopolymers. It will also left align insertions that
  are the same as the homopolymer. Think of it as shifting
  indels around.

- Inputs:
  1. samEntry structure pointer with read/sequence to
     left align
  2. reference sequence the read was mapped to
- Output:
  - Modifies the samEntry structure (input 1) to have a
    left aligned sequence
  - Returns:
    - 0: for no errors
    - 1: for memory errors

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>
#include "rmHomo.h"
#include "samEntry.h"
#include "seqST.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: declare variables and initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct samEntry samStackST; /*gets reads*/
   struct seqST refStackST;    /*gets reference sequence*/
   FILE *inFILE = 0;

   /*initialize*/
   init_samST(&samStackST);
   init_seqST(&refStackST);

   if( setup_samST(&samStackST) )
      goto memoryErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: open files and get reference/first read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   inFILE = fopen("reference.fasta", "r");

   if(! inFILE)
      goto fileErr_main_sec04;

   errorSC = getFa_seqST(inFILE, &refStackST);
   if(! errorSC)
      ;
   else if(errorSC = def_memErr_seqST)
      goto memoryErr_main_sec04;
   else
      goto fileErr_main_sec04;

   fclose(inFILE);
   inFILE = 0;


   inFILE = fopen("reads.sam", "r");
   if(! inFILE)
      goto fileErr_main_sec04;
   errorSC = get_samEntry(&samStackST, inFILE);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: left align sequences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errorSC)
   { /*Loop: read in sam file*/
      if(samStackST.extraStr[0] == '\0')
         p_samEntr(&samStackST, 0, stdout); /*header*/
      else if(samStackST.flagUS & 4)
         goto nextEntry_main_sec02; /*unmapped read*/
      else if(! samStackST.seqStr[0])
         goto nextEntry_main_sec02; /*no sequence*/
      else if(! samStackST.seqStr[0] == '*')
         goto nextEntry_main_sec02; /*no sequence*/

      errorSC =
         leftAlnIndel_rmHomo(
            &samStackST,
            refStackST.seqStr
         ); /*is the left alignment step*/

      if(errorSC)
         goto memoryErr_main_sec04;

      nextEntry_main_sec02:;
         errorSC = get_samEntry(&samStackST, inFILE);
   } /*Loop: read in sam file*/

   if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else if(errorSC != def_EOF_samEntry)
      goto fileErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto return_main_sec04;

   memoryErr_main_sec04:;
      fprintf(stderr, "memory error\n");
      errorSC = 1;
      goto return_main_sec04;

   fileErr_main_sec04:;
      fprintf(
         stderr,
         "error in reference.fasta or reads.sam\n"
      );
      errorSC = 2;
      goto return_main_sec04;

   return_main_sec04:;
      freeStack_samST(&samStackST);
      freeStack_seqST(&refStackST);

      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else if(inFILE == stdout) ;
      else if(inFILE == stderr) ;
      else
         fclose(inFILE);
      inFILE = 0;

      return errorSC;
} /*main*/
```
