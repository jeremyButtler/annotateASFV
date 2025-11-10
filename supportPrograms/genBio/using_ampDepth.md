# Goal

Describe how to use the ampDepth libary to get and
  print read depth histograms.

Files: ampDepth.c and ampDepth.h

- Dependencies:
  - samEnty.h
  - geneCoord.h


- Hidden Dependencies:
  - base10str from genLib
  - numToStr from genLib
  - ulCp from genLib
  - charCp from genLib
  - fileFun from genLib
  - ntTo5Bit.h
  - endLine.h from genLib
  - genMath.h from genLib

You should be familar with the samEntry library
  (samEntry.c/h and using\_samEntry.md) and the geneCoord
  (geneCoord.c/h and using\_geneCoord.md) library.

# ampDepth

The amp depth functions are here to build histograms
  with read depths and then print out the historgrams
  in a format for the user.

# Flow:

The flow is you first create an array the size of your
  reference sequence. You can realloc this array on the
  flye if needed, but you need to make sure all new
  entries are always set to 0.

Then you add the reads to the histogram.

Once all reads are added you can then print out the
  histogram with one of the functions.

## addRead\_ampDepth

The addRead\_ampDepth adds a read to the histogram.

- Input:
  1. samEntry structure pointer with read to add to the
     histogram
  2. geneCoord structure pointer with genes or genome
     targets for filtering
     - use 0/null to not filter by genome targets
  3. number of genes in the geneCoord structure (input 2)
  4. signed int array holding the depth for each base
     - this needs to be at least as long as the
       last reference coordinate (endRefUI) in the
       reads samEntry structure (input 1)
     - idealy, make as long as reference or if array is
       to small, realloc array to size and initialize
       the new entries with 0's
  5. signed int pointer to get the off target gene count
- Output:
  - adds read to the depth array
  - incurments off target input if filters removed the
    read

## printing the summary

The summary is a table that includes the fagrment
  size (areas connected with enough coverage) and the
  mean/median read depths for these fragments.

You can print the header for the file with
  `phead_ampDepth()`. The only input is the output file.

After that you can print the summary with
  the `phist_ampDepth()` function.

- Input:
  - signed int array with the histgram to print
    - built using addRead\_ampDepth
  - minimum depth to include the position in the summary
  - geneCoords structure pointer with coordinates used
    to filter reads (input 2 used with addRead\_ampDepth)
    - this can not be 0/null
  - number of genes in the geneCoords structure (input 3
    for addRead\_ampDepth)
  - number of off target reads (input 5 from
    addRead\_ampDepth)
  - number of unmapped reads
  - c-string with the run/output id
    - no spaces or tabs
  - file to save summary to

**Example:**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>

#include "geneCoord.h"
#include "samEntry.h"
#include "seqST.h"

int
main(
){
   signed int errSI = 0;
   signed int offTargetSI = 0;
   signed int unmappedSI = 0;

   signed int *depthHeapArySI = 0;

   struct geneCoord *geneCoordsHeapST = 0;
   signed int numerGenesSI = 0;
   signed char *geneTableFileStr =
      (signed char *) "genes.tsv";

   struct samEntry readStackST;

   FILE *samFILE = fopen("reads.sam", "r");

   init_samEntry(&readStachST);
   if( setup_samEntry(&readStackST) )
      /*deal with memory errors*/
   

   geneCoordsHeapST =
      getCoords_geneCoord(
         geneTableFileStr,
         &numberGenesSI, /*is index 0; so last gene*/
         &errorUL /*gets error type*/
      );
   if(! geneCoordsHeapST)
      /*use errorUL to find error and deal with it*/


   depthHeapArySI =
      calloc(
         genesHeapST->endAryUI[numberGenesSI] + 2,
         sizeof(signed int)
      );
      /*the last entry in the genesHeapST structure is
      `  always the last position (sorted by start
      `  position)
      */
   if(! depthHeapArySI)
      /*deal with memory error*/

   errSI = get_samEntry(&readStackST, samFILE);

   if(errSI)
      /*memory error or blank file*/

   while(! errSI)
   { /*Loop: read in sam file*/
      if(readStackST.extraStr[0] == '@')
         goto nextEntry; /*header*/

      else if(readStackST.flagUS & 4)
      { /*Else If: unmapped read*/
         ++unmappedSI;
         goto nextEntry; /*header*/
      } /*Else If: unmapped read*/

      else if(readStackST.flagUS & (256 | 2048))
         goto nextEntry;
         /*I do not want supplemental or secondary
         `   alignments. It means the same read mapped
         `   to multiple locations/references
         */

      addRead_ampDepth(
         &readStackST,
         geneCoordsHeapST,
         numberGenesSI,
         depthHeapArySI,
         &offTargetSI
      ); /*add read to the depth array*/

      nextEntry:;
         errSI = get_samEntry(&readStackST, samFILE);
   } /*Loop: read in sam file*/

   fclose(samFILE);
   samFILE = 0;

   if(errSI != def_EOF_samEntry)
       /*deal with file error or memory error*/

   phead_ampDepth(stdout);
   phist_ampDepth(
      depthHeapArySI,
      10,             /*minimum read depth*/
      geneCoordsHeapST,
      numberGenesSI,
      offTargetSI,
      unmappedSI,
      (signed char *) "out", /*id column of run*/
      stdout
   );

   if(depthHeapArySI)
      free(depthHeapArySI);
   depthHeapArySI = 0;

   freeStack_samEntry(&readStackST);

   if(geneCoordsHeapST)
      freeHeap_geneCoord(geneCoordsHeapST);
   geneCoordsHeapST = 0;

   return 0;
} /*main*/
```

## printing the depths

The read depth reports prints the number of bases the
  supported each reference position. It is filtered
  by the minimum read depth.

First add the reads to the read depth histogram, like in
  the print the summary step. Then
  call `pDeptHead_ampDepth()` to print the header. After
  that call `pdepth_ampDepth()` to print the read depth
  file.

The `pDepthHead_ampDepth()` function takes the output
  FILE (as pointer) as input and prints the header to the
  output file.

The `pdepth_ampDepth()` function takes a bit more input.
  The output is the read depths for each position that
  either has at least the minimum read depth.

- Input:
  1. signed int array with read depths for each position
    (input 4 for addRead\_ampDepth)
  2. size of the read depth array (input 1)
  3. minium read depth to keep a position
  4. c-string with reference name/id
     - do not use tabs or spaces
  5. c-string with identifier for run
     - do not use tabs or spaces
  6. FILE pointer with file to save printed depths to
- Output:
  - prints read depths above minimum read depth (input 3)
    to output file (input 6)

**Example:**

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
   signed char *referenceNameStr = "reference";
   signed int errSI = 0;
   signed int offTargetSI = 0;
   signed int unmappedSI = 0;

   signed int *swapArySI = 0;
   signed int *depthHeapArySI = 0;
   unsigned int *depthSizeUI = 0;
   struct samEntry readStackST;

   FILE *samFILE = fopen("reads.sam", "r");

   init_samEntry(&readStachST);
   if( setup_samEntry(&readStackST) )
      /*deal with memory errors*/
   
   /*get some initial memory*/
   depthSizeUI = 1024;
   depthHeapArySI =
      calloc((depthSizeUI + 8), sizeof(unsigned int));
   if(! depthHeapArySI)
      /*deal with memory error*/

   /*get first entry from sam file*/
   errSI = get_samEntry(&readStackST, samFILE);
   if(errSI)
      /*memory error or blank file*/

   while(! errSI)
   { /*Loop: read in sam file*/
      if(readStackST.extraStr[0] == '@')
         goto nextEntry; /*header*/

      else if(readStackST.flagUS & 4)
      { /*Else If: unmapped read*/
         ++unmappedSI;
         goto nextEntry; /*header*/
      } /*Else If: unmapped read*/

      else if(readStackST.flagUS & (256 | 2048))
         goto nextEntry;
         /*I do not want supplemental or secondary
         `   alignments. It means the same read mapped
         `   to multiple locations/references
         */

      /*make sure array covers full reference length*/
      if(readStackST.refEndUI > readMapSizeUI)
      { /*If: array needs to be larger*/
         swapArySI =
           realloc(
              depthHeapArySI,
              (readStackST.refEndUI + 1032)
                * sizeof(signed int)
           );
         if(! swapArySI)
         { /*If: had a memory error*/
            fprintf(
               stderr,
               "Memory error (main sec04 sub02)%s",
               str_endLine
            );

            goto err_main_sec06_sub02;
         } /*If: had a memory error*/

         depthHeapArySI = swapArySI;
         swapArySI = 0;
         errSI = depthSizeUI;
         depthSizeUI = readStackST.refEndUI + 1024;

         while(errSI < (signed int) depthSizeUI)
            depthHeapArySI[errSI++] = 0;
         errSI = 0;
      }  /*If: array needs to be larger*/

      addRead_ampDepth(
         &readStackST,
         geneCoordsHeapST,
         numberGenesSI,
         depthHeapArySI,
         &offTargetSI
      ); /*add read to the depth array*/

      nextEntry:;
         errSI = get_samEntry(&readStackST, samFILE);
   } /*Loop: read in sam file*/

   fclose(samFILE);
   samFILE = 0;

   if(errSI != def_EOF_samEntry)
       /*deal with file error or memory error*/

   pDepthHead_ampDepth(stdout);
   pdepth_ampDepth(
      depthHeapArySI,
      depthSizeUI,
      10,             /*minimum read depth*/
      referenceNameStr,
      (signed char *) "out", /*id column of run*/
      stdout
   );

   if(depthHeapArySI)
      free(depthHeapArySI);
   depthHeapArySI = 0;

   freeStack_samEntry(&readStackST);

   return 0;
} /*main*/
```
