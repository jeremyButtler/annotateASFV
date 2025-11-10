# Goal:

Describe how to use the adjCoords (adjust coordinates)
  library to adjust sam files mapped to genes to reference
  global coordinates.

Files: adjCoords.c adjCoords.h

- Dependencies:
  - geneCoord
  - samEntry
  - ulCp from genLib

- Hidden Dependencies:
  - base10str from genLib
  - numToStr from genLib
  - charCp from genLib
  - fileFun from genLib
  - ntTo5Bit.h
  - endLine.h from genLib
  - genMath.h from genLib

You should be familar with the samEntry library
  (samEntry.c/h and using\_samEntry.md) and the geneCoord
  (geneCoord.c/h and using\_geneCoord.md) library.

# Overview:

The ajdCoords library is designed to read in sam files
  that are mapped to genes from a reference and ajust to
  coordinates to the correct location on the reference.
  Not very usefull, but can save some time in some cases
  when only a few genes are needed.

There is only one function and that is `adjCoords`.

- Input
  1. samEntry structure pointer with the read to adjust
  2. geneCoord structure pointer with gene mapping
     coordinates on the reference genome
  3. number of genes (index 0) in the geneCoord structure
     (input 2)
- Output:
  - the reference position in the samEntry structure is
    adjusted to the reference position

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>

#include "adjCoords.h"
#include "geneCoord.h"
#include "samEntry.h"

int
main(
){
   unsigned long errUL = 0;
   signed int errSI = 0;

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
   sortName_geneCoord(coordsHeapST, 0, numGenesSI);
      /*need to sort by name for gene lookups*/


   errSI = get_samEntry(&readStackST, samFILE);

   if(errSI)
      /*memory error or blank file*/

   while(! errSI)
   { /*Loop: read in sam file*/
      if(readStackST.extraStr[0] == '@')
      { /*If: header*/
         p_samEntry(&readStackST, 0, stdout);
         goto nextEntry; /*header*/
      } /*If: header*/

      errSI =
         adjCoords(
            &readStackST,
            geneCoordsHeapST,
            numberGenesSI
         );

      if(! errSI)
         p_samEntry(&readStackST, 0, stdout);
         /*read has been adjusted*/
      /*else read is not in list, do not print*/

      nextEntry:;
         errSI = get_samEntry(&readStackST, samFILE);
   } /*Loop: read in sam file*/

   fclose(samFILE);
   samFILE = 0;

   if(errSI != def_EOF_samEntry)
       /*deal with file error or memory error*/

   freeStack_samEntry(&readStackST);

   if(geneCoordsHeapST)
      freeHeap_geneCoord(geneCoordsHeapST);
   geneCoordsHeapST = 0;

   return 0;
} /*main*/
```
