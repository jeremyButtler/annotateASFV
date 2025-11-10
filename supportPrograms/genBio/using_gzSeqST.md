# Goal

Describe how to use the gzSeqST library to read in gzipped
  compressed fastq and fasta files.

Files: gzSeqST.c gzSeqST.h

- Dependencies:
  - seqST
  - ulCp from genLib
  - fileFun from genLib
  - inflate from genLib
- Hidden Dependencies
  -   checkSum from genLib
  -   genMath from genLib
  -   endin from genLib
  -   endLine.h

You should be familar with the seqST (seqST.c/h and
  using\_seqST.md) library and the inflate
  (../genLib/inflate.c/h and ../genLib/using\_inflate.md)
  library.

# Overview

gzSeqST is separated from seqST because the gz inflate
  system provides a good boos to binary size. There is
  only one function in gzSeqST and that is `get_gzSeqST`.

- Input:
  1. a `file_inflate` structure pointer to hold the input
     file
  2. signed char pionter to be set to the file type input
     - get\_gzSeqST will figure out if it is a fasta,
       fastq, fasta.gz, or fastq.gz, so do not mess with
       this
  3. seqST structure pointer to get the next sequence from
     the file
  4. the file to work on or 0/null if working on a
     previous file
     - never use `fclose()` on FILES provided by input 4,
       because the file is stored in the file\_inflate
       structure (input 4). It will be closed when a new
       file is input
- Output:
  - 0 for no errors
  - def\_EOF\_seqST is at the end of the file
  - def\_memErr\_seqST for memory errors
  - def\_badLine\_seqST | def\_fileErr\_seqST for an
    invalid line
  - def\_fileErr\_seqST for any file errors
- Note: seqST gets the next sequence

**Example:**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>

#include "../genLib/inflate.h"
#include "seqST.h"
#include "gzSeqST.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed long readCountSL = 0;
   signed char fileTypeSC = 0;

   struct seqST readStackST;
   struct file_inflate inflateStackST;

   FILE inFILE = fopen("reads.fastq.gz", "rb");
      /*this file needs to be read in as a binary file
      `   otherwise widows will through errors
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initilize and get first read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_seqST(&readStackST);
   init_file_inflate(&inflateStackST);

   errSC =
      get_gzSeqST(
         &inflateStackST,
         &fileTypeSC, /*do not mess with this*/
         &readStackST,
         inFILE
      );
   inFILE = 0;
      /*the internal inflate functions will
      `   handle closing this file
      */

   if(! errSC)
      ; /*no errors; do nothing*/
   else if(errSC == def_EOF_seqST)
      goto cleanUp_main_sec04;
   else if(errSC == def_memErr_seqST)
      goto cleanUp_main_sec04;
   else if(errSC & def_badLine_seqST)
      goto cleanUp_main_sec04;
   else
      goto cleanUp_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec03: read rest of file and print as fasta file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errSC)
   { /*Loop: get reads from fastx file*/
       ++readCountSL;

       fprintf(
           stdout,
           ">%s\n%s\n",
           readStackST.idStr,
           readStackST.seqStr,
       ); /*convert fastq file to fasta file*/

      errSC =
         get_gzSeqST(
            &inflateStackST,
            &fileTypeSC, /*do not mess with this*/
            &readStackST,
            0            /*inflateStackST has file in it*/
         );
         /*to change the file, just replace the 0 with the
         `   file. get_gzSeqST will close the older file
         `   before adding the new file
         */
   } /*Loop: get reads from fastx file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec04: clean up the variables and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   cleanUp_main_sec04:;
      freeStack_seqST(&readStackST);
      freeStack_file_inflate(&inflateStackST);

      if(! readCountSL)
      { /*If: no reads in file*/
         fprintf(stderr, "no reads in reads.fastq.gz\n");
         return 4;
      } /*If: no reads in file*/
 
      else if(errSC == def_EOF_seqST)
         return 0; /*end of file; do nothing*/

      else if(errSC == def_memErr_seqST)
      { /*Else If: memory error*/
         fprintf(
            stderr,
            "ran out of memory on read %li in %s\n",
            readCountSL,
            "reads.fastq.gz"
         );
         return 1;
      } /*Else If: memory error*/

      else if(errSC & def_badLine_seqST)
      { /*Else If: bad line*/
         fprintf(
            stderr,
            "line %li in  is not a fastx entry\n",
            readCountSL,
            "reads.fastq.gz"
         );
         return 2;
      } /*Else If: bad line*/

      else
      { /*Else If: error reading file*/
         fprintf(
            stderr,
            "error reading reads.fastq.gz on line %li\n",
            readCountSL
         );
         return 3;
      } /*Else If: error reading file*/
} /*main*/
```
