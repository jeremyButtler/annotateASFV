# Goal

Describe how to use the cigToEqx libary to convert sam
  file cigars `MDISH` to eqx cigars `=XDISH`.

Files: cigToEqx.c and cigToEqx.h

- Dependencies:
  - samEnty.h
  - ntTo5Bit.h


- Hidden Dependencies:
  - base10str from genLib
  - numToStr from genLib
  - ulCp from genLib
  - fileFun from genLib
  - endLine.h from genLib

You should be familar with the samEntry library
  (samEntry.c/h and using\_samEntry.md). It might be
  also helpfull to have some familiarity with seqST
  (using\_seqST.md), but is not needed.

# cigToEqx

cigToEqx has the cigToEqx function to convert non-eqx
  cigars to eqx cigars. Beyond that, this does nothing
  else.

- Input:
  1. samEntry structure pointer with sam file entry to
     convert to an eqx cigar
  2. c-string with reference sequence to compare the cigar
     to
     - this does not need to be in "ATGC" format, but it
       does need to match the sequence in the samEntry
       structure
       - for example, I have used my indexToSeq\_alnSet()
         function (genAln/alnSet.c/h) to convert the
         sequences to numbers and then run them through
         this. You could also convert the sequences two
         two bit (0-3) numbers with ntTo2Bit.h and then
         compare. The only requirement is that they be
         in the same format.
  3. c-string pointer to be used in storing the temporary
     cigar. This is resized as needed
     - you will need to free this
  4. signed int array pointer to be used in storing the
     temporary cigar entry. This is resized as needed.
     - you will need to free this
  5. signed int pointer with current size of the cigar
     buffers (Input 3 and 4)
  6. if you want to treat anonymous bases as SNPs even
     when then can match
     - 1 all anonymous bases are mis-matches
     - 0 treads anonymous bases as matches if they match
- Output:
  - changes the cigar in the samEntry structure (input 1)
    to be an eqx cigar
  - Returns 0 for success
  - Returns 1 for memory errors
  - Returns 2 if the samEntry structure (input 1) had
    no sequence

**Example:** convert cigar entries to eqx cigar

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "seqST.h"

int
main(
){
   signed int errSI = 0;
   signed char *cigarString = 0;
   signed int *cigarValuesArySI = 0;
   signed int sizeBufferSI = 0;

   struct seqST referenceStackST;
   struct samEntry readStackST;

   FILE *refFILE = fopen("reference.fasta", "r");
   FILE *samFILE = fopen("reads.sam", "r");

   init_seqST(&refStackST);
   init_samEntry(&readStachST);

   if( setup_samEntry(&readStackST) )
      /*deal with memory errors*/
   
   if( getFa_seqST(refFILE, &referenceStackST))
      /*deal with no reference sequence
      `   either memory error (1) or no sequence (2)
      */
   fclose(refFILE);
   refFILE = 0;

   errSI = get_samEntry(&readStackST, samFILE);

   if(errSI)
      /*memory error or blank file*/

   while(! errSI)
   { /*Loop: read in sam file*/
      errSI =
         cigToEqx(
            &readStackST,
            referenceStackST.seqStr,
            &cigString,
            &cigarValuesArySI,
            &sizeBufferSI,
            0 /*treat anonymous bases as matches
         );

      if(errSI == 1)
         /*deal with memory errors*/
      else if(! errSI)
         ; /*no error, do nothing*/
      else if(readStackST.extraStr[0] 1= '@')
         goto nextEntry; /*entery had no sequence*/

      p_samEntry(&readStackST, 0, stdout);
         /*print header or converted sequence*/

      nextEntry:;
         errSI = get_samEntry(&readStackST, samFILE);
   } /*Loop: read in sam file*/

   fclose(samFILE);
   samFILE = 0;

   if(errSI != def_EOF_samEntry)
       /*deal with file error or memory error*/

   if(cigarString)
      free(cigarString);
   cigarString = 0;

   if(cigarValuesArySI)
      free(cigarValuesArySI);
   cigarValuesArySI = 0;

   freeStack_samEntry(&readStackST);
   freeStack_seqST(&referenceStackST);
   return 0;
} /*main*/
```
