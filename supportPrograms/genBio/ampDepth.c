/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ampDepth SOF: Start Of File
'   - Had functions for updating and printing an read
'     depth histogram (array of integers)
'   o header:
'     - indcluded libraries
'   o fun01: addRead_ampDepth
'     - adds a read to an amp depth arrary (histogram)
'   o fun02: phead_ampDepth
'     - Prints the header for an ampDepth tsv file
'   o fun03: pDepthHead_ampDepth
'     - prints the header for the read depth output
'   o fun04: phist_ampDepth
'     - Prints an histogram of read depths as an tsv with
'       individual genes for ampDepth
'   o fun05: pDepthHead_ampDepth
'     - prints the header for the read depth output
'   o fun06: pdepth_ampDepth
'     - prints the read depth of each base
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - indcluded libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "ampDepth.h"

#include <stdio.h>

#include "geneCoord.h"
#include "samEntry.h"

/*No .c files*/
#include "../genLib/endLine.h"
#include "../genLib/genMath.h" /*using .h macros*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/ulCp.h"
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: addRead_ampDepth
|   - adds a read to an amp depth arrary (histogram)
| Input:
|   - samSTPtr:
|     o pointer to an samEntry structure with a read to
|       add to the histogram
|   - coordsSTPtr:
|     o geneCoord struct pointer with gene/target
|       coordinates want to extract
|     o input 0/null to not filter by genes
|   - numGenesSI:
|     o number of genes in coordsSTPtr (index 1)
|   - depthArySI:
|     o histogram (int array) to add each base to
|   - numOffTargSI:
|     o number of reads not mapping to a gene coordiante,
|       but are mapped to the genome
| Output:
|   - Modifies:
|     o ampDepthIAry to be incurmented for each reference
|       base in the read
\-------------------------------------------------------*/
void
addRead_ampDepth(
   struct samEntry *samSTPtr,
   struct geneCoord *coordsSTPtr, /*list of genes*/
   signed int numGenesSI,   /*number genes in list*/
   signed int *depthArySI,  /*depth array to update*/
   signed int *numOffTargSI /*number reads not in list*/
){
   signed int siBase = 0;
   signed int siIndex = 0;
   signed int endSI = 0;

   if(coordsSTPtr)
   { /*If: doing coordinate filtering*/
      siIndex =
         findRange_geneCoord(
            coordsSTPtr,
            samSTPtr->refStartUI,
            samSTPtr->refEndUI,
            numGenesSI
         );

      if(siIndex < 0)
         ++*numOffTargSI;
      else
      { /*Else: gene has some on target coordiantes*/
         addBases_fun01:;
         siBase =
            max_genMath(
                (signed int) samSTPtr->refStartUI,
                (signed int)
                   coordsSTPtr->startAryUI[siIndex]
            );

         endSI =
            min_genMath(
               (signed int) samSTPtr->refEndUI,
               (signed int) coordsSTPtr->endAryUI[siIndex]
            );

         while(siBase <= endSI)
            ++depthArySI[siBase++];


         /*see if read has muttiple genes*/
         ++siIndex;

         if(siIndex >= numGenesSI)
            ; /*end of genes list*/
         else if(
              samSTPtr->refEndUI
            > coordsSTPtr->startAryUI[siIndex]
         ) goto addBases_fun01; /*have another gene*/
      } /*Else: gene has some on target coordiantes*/
   } /*If: doing coordinate filtering*/

   else
   { /*Else: not doing coordinate filtering*/
      siIndex = samSTPtr->refStartUI;
      while(siIndex <= (signed int) samSTPtr->refEndUI)
        ++depthArySI[siIndex++];
   } /*Else: not doing coordinate filtering*/
} /*addRead_ampDepth*/

/*-------------------------------------------------------\
| Fun02: phead_ampDepth
|   - prints the header for an ampDepth tsv file
| Input:
|   - outFILE:
|     o file to print the header to
| Output:
|   - Prints:
|     o ampDepth header to outFILE
\-------------------------------------------------------*/
void
phead_ampDepth(
   void *outFILE
){
   fprintf((FILE *) outFILE, "flag\tampNumber\trefStart");
   fprintf((FILE *) outFILE, "\trefEnd\tampStart");
   fprintf((FILE *) outFILE, "\tampEnd\tavgAmpDepth");
   fprintf((FILE *) outFILE, "\tminAmpDepth\tmaxAmpDepth");
   fprintf((FILE *) outFILE, "\tgeneId\trefGeneStart");
   fprintf((FILE *) outFILE, "\trefGeneEnd");
   fprintf((FILE *) outFILE, "\tampGeneStart");
   fprintf((FILE *) outFILE, "\tampGeneEnd");
   fprintf((FILE *) outFILE, "\tfirstBaseDepth");
   fprintf((FILE *) outFILE, "\tlastBaseDepth");
   fprintf((FILE *) outFILE, "\tavgDepth\tminDepth");
   fprintf(
      (FILE *) outFILE,
      "\tmaxDepth%s",
      str_endLine
   );
} /*phead_ampDepth*/

/*-------------------------------------------------------\
| Fun04: phist_ampDepth
|   - prints an histogram of read depths as an tsv with
|     individual genes for ampDepth
| Input:
|   - histArySI:
|     o integer array with the histogram to print out
|   - minDepthSI:
|     o integer with the min depth to keep an histogram
|       entry
|   - geneCoordSTPtr:
|     o pointer to an geneCoord structure array with the
|       locations of each gene
|   - numGenesSI:
|     o integer with number of genes in geneCoordSTPtr
|   - offTargSI:
|     o number off of target reads to print out
|   - noMapSI:
|     o number off of unmapped reads to print out
|   - extraColStr:
|     o c-string with the first column to print out
|   - outFILE:
|     o file to print to
| Output:
|   - Prints:
|     o stats for each gene to outFILE
\-------------------------------------------------------*/
void
phist_ampDepth(
   signed int *histArySI,
   signed int minDepthSI,
   struct geneCoord *geneCoordSTPtr,
   signed int numGenesSI,
   signed int offTargSI,
   signed int noMapSI,
   signed char *extraColStr,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04: phist_ampDepth
   '    - prints out an ampDepth histogram to file
   '    o fun04 sec01:
   '      - variable declerations
   '    o fun04 sec02:
   '      - find the mapping positions
   '    o fun04 sec03:
   '      - print out the off targets an unampped reads
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^  - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    signed int mapStartSI = 0;
    signed int tmpStartSI = 0;

    signed int ampNumSI = -1;
    signed int geneIndexSI = 0;
    signed int tmpSI = 0;

    signed int readsAtStartSI = 0;
    signed int readsAtEndSI = 0;

    /*Stats for individual genes*/
    signed int maxReadsSI = 0;
    signed int minReadsSI = 0;
    signed long avgDepthSL = 0;

    /*Stats for the amplicon*/
    signed int ampStartSI = 0;
    signed int ampEndSI = 0;
    signed int ampGeneEndSI = 0;
    signed int ampGeneStartSI = 0;

    signed long ampAvgDepthSL = 0;
    signed int ampMaxReadsSI = 0;
    signed int ampMinReadSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - find the mapping positions
   ^   o fun04 sec02 sub01:
   ^     - check if positon has enough depth to keep
   ^   o fun04 sec02 sub02:
   ^     - find the amplicon stats
   ^  o fun04 sec02 sub03:
   ^    - find the stats for each gene in the amplicon
   ^  o fun04 sec02 sub04:
   ^    - print the stats for each gene in the amplicon
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec02 Sub01:
   *   - check if positon has enough depth to keep
   \*****************************************************/

   mapStartSI = 0;

   while(
         mapStartSI
      <= (signed int)
         geneCoordSTPtr->endAryUI[numGenesSI] + 1
   ){ /*Loop: Get the gene positions that mapped*/
      ++mapStartSI;

      if(
           mapStartSI
         < (signed int) geneCoordSTPtr->startAryUI[0]
      ) continue;

      geneIndexSI = 
         findStart_geneCoord(
            geneCoordSTPtr,
            mapStartSI,
            numGenesSI
         ); /*Find the index of the gene at the position*/

      if(geneIndexSI < 0)
         continue; /*Unmapped base*/

      if(histArySI[mapStartSI] < minDepthSI)
         continue;

      /**************************************************\
      * Fun04 Sec02 Sub02:
      *   - find amplicon stats
      \**************************************************/

      ++ampNumSI;

      ampGeneEndSI = geneIndexSI;

      ampStartSI = mapStartSI;

      ampGeneStartSI =
         geneCoordSTPtr->startAryUI[geneIndexSI];

      ampAvgDepthSL = 0;
      ampMaxReadsSI = histArySI[mapStartSI];
      ampMinReadSI = ampMaxReadsSI;

      for(
         ampEndSI = mapStartSI;
         histArySI[ampEndSI] >= minDepthSI;
         ++ampEndSI
      ){ /*Loop: Find end of region*/
          ampGeneEndSI +=
            (    ampEndSI
               > (signed int)
                 geneCoordSTPtr->endAryUI[ampGeneEndSI]
            );

          ampMaxReadsSI = 
             max_genMath(
                ampMaxReadsSI,
                histArySI[ampEndSI]
             ); /*Find the maximum*/

          ampMinReadSI = 
             min_genMath(
                ampMinReadSI,
                histArySI[ampEndSI]
             ); /*Find the maximum*/
          
          ampAvgDepthSL += histArySI[ampEndSI];
      } /*Loop: Find end of region*/

      --ampEndSI; /*Account for overcounting*/

      ampGeneEndSI =
         geneCoordSTPtr->endAryUI[ampGeneEndSI];

      /*Using interger mean, not decimal*/
      tmpSI = (ampEndSI - ampStartSI);

      if(tmpSI)
         ampAvgDepthSL /= tmpSI;

      /**************************************************\
      * Fun04 Sec02 Sub03:
      *   - find the stats for each gene in the amplicon
      \**************************************************/

      tmpSI = 0;

      nextGene_fun04_sec02_sub03:

      ++tmpSI;
      tmpStartSI = mapStartSI;
      avgDepthSL = 0;
      readsAtStartSI = histArySI[mapStartSI];
      maxReadsSI = readsAtStartSI;
      minReadsSI = maxReadsSI;

      while(
            mapStartSI
         <= (signed int)
            geneCoordSTPtr->endAryUI[geneIndexSI]
      ){ /*Loop: Check if gene is complete*/
          if(histArySI[mapStartSI] < minDepthSI) break;

          maxReadsSI = 
             max_genMath(
                maxReadsSI,
                histArySI[mapStartSI]
             ); /*Find the maximum*/

          minReadsSI = 
             min_genMath(
                minReadsSI,
                histArySI[mapStartSI]
             ); /*Find the maximum*/
          
          avgDepthSL += histArySI[mapStartSI];

         ++mapStartSI;
      } /*Loop: Check if gene is complete*/

      readsAtEndSI = histArySI[mapStartSI - 1];

      /*Using interger mean, not decimal*/
      tmpSI = (mapStartSI - 1 - tmpStartSI);

      if(tmpSI)
         avgDepthSL /= tmpSI;

      /**************************************************\
      * Fun04 Sec02 Sub04:
      *   - print the stats for each gene in the amplicon
      \**************************************************/


      fprintf(
         (FILE *) outFILE,
         "%s\t%i\t%i\t%i\t%i\t%i\t%li\t%i\t%i",
         extraColStr,
         ampNumSI,
         ampGeneStartSI + 1,
         ampGeneEndSI + 1,
         ampStartSI + 1,
         ampEndSI + 1,
         ampAvgDepthSL,
         ampMinReadSI,
         ampMaxReadsSI
      ); /*Print out the starting and ending position*/

      fprintf(
         (FILE *) outFILE,
         "\t%s\t%u\t%u\t%i\t%i\t%i\t%i\t%li\t%i\t%i%s",
         geneCoordSTPtr->idStrAry[geneIndexSI],
         geneCoordSTPtr->startAryUI[geneIndexSI] + 1,
         geneCoordSTPtr->endAryUI[geneIndexSI] + 1,
         (signed int)
             max_genMath(
                (signed int)
                  geneCoordSTPtr->startAryUI[geneIndexSI],
                ampStartSI
             ) + 1, /*first gene base in amplicon*/
         (signed int)
             min_genMath(
                (signed int)
                   geneCoordSTPtr->endAryUI[geneIndexSI],
                ampEndSI
             ) + 1, /*last gene base in amplicon*/
         readsAtStartSI,
         readsAtEndSI,
         avgDepthSL,
         minReadsSI,
         maxReadsSI,
         str_endLine
      );

      if(histArySI[mapStartSI] > minDepthSI)
      { /*If: I have another gene*/
         ++geneIndexSI;
         if(geneIndexSI > numGenesSI) continue;;

         /*Make sure I am on the next gene*/
         while(
              mapStartSI
            < (signed int)
              geneCoordSTPtr->startAryUI[geneIndexSI]
         ){ /*Loop: Move to next gene*/
            if(histArySI[mapStartSI] < minDepthSI)
               break;

            ++mapStartSI;
         } /*Loop: Move to next gene*/

         if(histArySI[mapStartSI] < minDepthSI)
            continue;
 
         goto nextGene_fun04_sec02_sub03;
      } /*If: I have another gene*/
   } /*Loop: Get the gene positions that mapped*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec03:
   ^  - print off targets and unampped reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ++ampNumSI;

   fprintf(
      (FILE *) outFILE,
      "%s\t%i\tNA\tNA\tNA\tNA\t%i\t%i\t%i",
      extraColStr,
      ampNumSI,
      offTargSI,
      offTargSI,
      offTargSI
   ); /*Print out the off target reads*/

   fprintf(
      (FILE *) outFILE,
      "\tx-off-target\tNA\tNA\tNA\tNA\t%i\t%i",
      offTargSI,
      offTargSI
   );

   fprintf(
      (FILE *) outFILE,
      "\t%i\t%i\t%i%s",
      offTargSI,
      offTargSI,
      offTargSI,
      str_endLine
   );

   ++ampNumSI;

   fprintf(
      (FILE *) outFILE,
      "%s\t%i\tNA\tNA\tNA\tNA\t%i\t%i\t%i",
      extraColStr,
      ampNumSI,
      noMapSI,
      noMapSI,
      noMapSI
   ); /*Print out the off target reads*/

   fprintf(
      (FILE *) outFILE,
      "\tz-unmapped\tNA\tNA\tNA\tNA\t%i\t%i",
      noMapSI,
      noMapSI
   );

   fprintf(
      (FILE *) outFILE,
      "\t%i\t%i\t%i%s",
      noMapSI,
      noMapSI,
      noMapSI,
      str_endLine
   );
} /*phist_ampDepth*/

/*-------------------------------------------------------\
| Fun05: pDepthHead_ampDepth
|   - prints the header for the read depth output
| Input:
|   - outFILE:
|     o file to print the header to
| Output:
|   - Prints:
|     o read depth header to outFILE
\-------------------------------------------------------*/
void
pDepthHead_ampDepth(
   void *outFILE
){
   fprintf((FILE *) outFILE, "flag\treference\tbase");
   fprintf((FILE *) outFILE, "\tdepth%s", str_endLine);
} /*pDepthHead_ampDepth*/

/*-------------------------------------------------------\
| Fun06: pdepth_ampDepth
|   - prints the read depth of each base
| Input:
|   - depthArySI:
|     o integer array with the depthogram to print out
|   - depthLenSI:
|     o number of bases (length) of depthArySI (index 1)
|   - minDepthSI:
|     o integer with the min depth to keep an depthogram
|       entry
|   - refStr:
|     o c-string with name of reference sequence
|     o if 0/null then "reference"
|   - flagStr:
|     o c-string with flag to go in left column
|     o if 0/null then "out"
|   - outFILE:
|     o file to print to
| Output:
|   - Prints:
|     o depth for each base to outFILE
\-------------------------------------------------------*/
void
pdepth_ampDepth(
   signed int *depthArySI,  /*has read depths*/
   signed int depthLenSI,   /*number bases in depthArySI*/
   signed int minDepthSI,   /*minimum read depth*/
   signed char *refStr,     /*name of reference*/
   signed char *flagStr,    /*flag to go in left column*/
   void *outFILE            /*output file*/
){

   signed int indexSI = 0;

   if(! flagStr)
      flagStr = (signed char *) "out";
   if(! refStr)
      flagStr = (signed char *) "reference";

   while(indexSI < depthLenSI)
   { /*Loop: print bases*/
      if(depthArySI[indexSI] < minDepthSI)
      { /*Else If: position is filtered by read depth*/
         ++indexSI;
         continue;
      } /*Else If: position is filtered by read depth*/

      fprintf(
         (FILE *) outFILE,
         "%s\t%s\t%i\t%i%s",
         flagStr,
         refStr,
         indexSI,
         depthArySI[indexSI],
         str_endLine
      );
      ++indexSI;
   } /*Loop: print bases*/
} /*pdepth_ampDepth*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy, modify, publish, use, compile,
:   sell, or distribute this software, either in source
:   code form or as a compiled binary, for any purpose,
:   commercial or non-commercial, and by any means.
: 
: In jurisdictions that recognize copyright laws, the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
:   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information, please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted, free of charge, to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"), to
:   deal in the Software without restriction, including
:   without limitation the rights to use, copy, modify,
:   merge, publish, distribute, sublicense, and/or sell
:   copies of the Software, and to permit persons to whom
:   the Software is furnished to do so, subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
