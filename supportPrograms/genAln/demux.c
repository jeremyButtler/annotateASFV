/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' demux SOF: Start Of File
'   - has functions to demultiplexes reads by barcode
'   o header:
'     - included libraries
'   o .c fun01: sortBarcodes_demux
'     - sorts the mapped barcodes by their barcode and
'       then by start coordiante
'   o .c fun02: sortStart_demux
'     - sorts the mapped barcodes by starting coordinates
'   o fun03: barcodeCoords_demux
'     - demux a read
'   o fun04: read_demux
'     - convert barcode coordinates to demuxed reads
'   o fun05: primer_demux
'     - get primer target regins from the input sequence
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genBio/seqST.h"
#include "alnSet.h"
#include "kmerFind.h"

/*.h files only*/
#include "../genLib/endLine.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .c  #include "../genLib/ulCp.h"
!   - .c  #include "../genLib/base10Str.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .c  #include "../genLib/shellSort.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .c  #include "memwater.h"
!   - .h  #include "../genLib/genMath.h"
!   - .h  #include "alnDefs.h"
!   - .h  #include "../genBio/kmerBit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: sortBarcodes_demux
|   - sorts the mapped barcodes by their barcode and then
|     by start coordiante
| Input:
|   - barAryUS:
|     o unsigned short array with barcodes
|   - startArySI:
|     o signed int array with starting positions
|   - endArySI:
|     o signed int array with ending coordinates
|   - scoreArySI:
|     o signed int array with mapping scores
|   - lenSI:
|     o number of items to sort (index 1)
| Output:
|   - Modifies:
|     o barAryUS, startArySI, endArySI, and scoreArySI
|       to be sorted by barAryUS and then startArySI
\-------------------------------------------------------*/
void
sortBarcodes_demux(
   unsigned short *barAryUS,
   signed int *startArySI,
   signed int *endArySI,
   signed int *scoreArySI,
   signed int lenSI
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - sorts the mapped barcodes by their barcode and
   '     then by start coordiante
   '   o fun01 sec01:
   '     - variable declerations
   '   o fun01 sec02:
   '     - find the number of rounds to sort for
   '   o fun01 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   signed int subSI = 0;
   signed int nextSI = 0;
   signed int lastSI = 0;
   signed int onSI = 0;

   /*Variables to incurment loops*/
   signed int siIndex = 0;
   signed int siElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   --lenSI;
   if(! lenSI)
      return;

   /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSI = 1; /*Initialzie first array*/

   while(subSI < lenSI - 1)
      subSI = (3 * subSI) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSI > 0)
   { /*Loop: all rounds*/
      for(
         siIndex = 0;
         siIndex <= subSI;
         ++siIndex
      ){ /*Loop: though sub array*/
         siElm = siIndex;

         for(
            siElm = siIndex;
            siElm + subSI <= lenSI;
            siElm += subSI
         ){ /*Loop: swap elements in subarray*/
            nextSI = siElm + subSI;

            if(barAryUS[siElm] >  barAryUS[nextSI])
               goto swap_fun01_sec03;
            else if(barAryUS[siElm] <  barAryUS[nextSI])
               ;
            else if(startArySI[siElm] >startArySI[nextSI])
            { /*Else If I need to swap an element*/
               swap_fun01_sec03:;

               barAryUS[siElm] ^= barAryUS[nextSI];
               barAryUS[nextSI] ^= barAryUS[siElm];
               barAryUS[siElm] ^= barAryUS[nextSI];

               startArySI[siElm] ^= startArySI[nextSI];
               startArySI[nextSI] ^= startArySI[siElm];
               startArySI[siElm] ^= startArySI[nextSI];

               endArySI[siElm] ^= endArySI[nextSI];
               endArySI[nextSI] ^= endArySI[siElm];
               endArySI[siElm] ^= endArySI[nextSI];

               scoreArySI[siElm] ^= scoreArySI[nextSI];
               scoreArySI[nextSI] ^= scoreArySI[siElm];
               scoreArySI[siElm] ^= scoreArySI[nextSI];

               lastSI = siElm;
               onSI = siElm;

               while(lastSI >= subSI)
               { /*loop: move swapped element back*/
                  lastSI -= subSI;

                  if(barAryUS[onSI] >  barAryUS[lastSI])
                      break; /*Positioned the element*/
                  else if(
                     barAryUS[onSI] <  barAryUS[lastSI]
                  ) break;
                  else if(
                    startArySI[onSI] > startArySI[lastSI]
                  ) break; /*Positioned the element*/
   
                  barAryUS[onSI] ^= barAryUS[lastSI];
                  barAryUS[lastSI] ^= barAryUS[onSI];
                  barAryUS[onSI] ^= barAryUS[lastSI];

                  startArySI[onSI] ^= startArySI[lastSI];
                  startArySI[lastSI] ^= startArySI[onSI];
                  startArySI[onSI] ^= startArySI[lastSI];

                  endArySI[onSI] ^= endArySI[lastSI];
                  endArySI[lastSI] ^= endArySI[onSI];
                  endArySI[onSI] ^= endArySI[lastSI];

                  scoreArySI[onSI] ^= scoreArySI[lastSI];
                  scoreArySI[lastSI] ^= scoreArySI[onSI];
                  scoreArySI[onSI] ^= scoreArySI[lastSI];

                  onSI = lastSI;
               } /*Loop: move swapped element back*/
            } /*Else If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSI = (subSI - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sortBarcodes_demux*/

/*-------------------------------------------------------\
| Fun02: sortStart_demux
|   - sorts the mapped barcodes by starting coordinates
| Input:
|   - barArySS:
|     o signed short array with barcodes
|   - startArySI:
|     o signed int array with starting positions
|   - endArySI:
|     o signed int array with ending coordinates
|   - scoreArySI:
|     o signed int array with mapping scores
|   - lenSI:
|     o number of items to sort (index 1)
| Output:
|   - Modifies:
|     o barArySS, startArySI, endArySI, and scoreArySI
|       to be sorted by startArySI
\-------------------------------------------------------*/
void
sortStart_demux(
   signed short *barArySS,
   signed int *startArySI,
   signed int *endArySI,
   signed int *scoreArySI,
   signed int lenSI
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - sorts the mapped barcodes by starting coordinates
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - find the number of rounds to sort for
   '   o fun02 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   signed int subSI = 0;
   signed int nextSI = 0;
   signed int lastSI = 0;
   signed int onSI = 0;

   /*Variables to incurment loops*/
   signed int siIndex = 0;
   signed int siElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   --lenSI;
   if(! lenSI)
      return;

   /*Recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSI = 1; /*Initialzie first array*/

   while(subSI < lenSI - 1)
      subSI = (3 * subSI) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSI > 0)
   { /*Loop: all rounds*/
      for(
         siIndex = 0;
         siIndex <= subSI;
         ++siIndex
      ){ /*Loop: though sub array*/
         siElm = siIndex;

         for(
            siElm = siIndex;
            siElm + subSI <= lenSI;
            siElm += subSI
         ){ /*Loop: swap elements in subarray*/
            nextSI = siElm + subSI;

            if(startArySI[siElm] >  startArySI[nextSI])
            { /*If: I need to swap an element*/
               barArySS[siElm] ^= barArySS[nextSI];
               barArySS[nextSI] ^= barArySS[siElm];
               barArySS[siElm] ^= barArySS[nextSI];

               startArySI[siElm] ^= startArySI[nextSI];
               startArySI[nextSI] ^= startArySI[siElm];
               startArySI[siElm] ^= startArySI[nextSI];

               endArySI[siElm] ^= endArySI[nextSI];
               endArySI[nextSI] ^= endArySI[siElm];
               endArySI[siElm] ^= endArySI[nextSI];

               scoreArySI[siElm] ^= scoreArySI[nextSI];
               scoreArySI[nextSI] ^= scoreArySI[siElm];
               scoreArySI[siElm] ^= scoreArySI[nextSI];

               lastSI = siElm;
               onSI = siElm;

               while(lastSI >= subSI)
               { /*loop: move swapped element back*/
                  lastSI -= subSI;

                  if(startArySI[onSI] >startArySI[lastSI])
                     break; /*Positioned the element*/
   
                  barArySS[onSI] ^= barArySS[lastSI];
                  barArySS[lastSI] ^= barArySS[onSI];
                  barArySS[onSI] ^= barArySS[lastSI];

                  startArySI[onSI] ^= startArySI[lastSI];
                  startArySI[lastSI] ^= startArySI[onSI];
                  startArySI[onSI] ^= startArySI[lastSI];

                  endArySI[onSI] ^= endArySI[lastSI];
                  endArySI[lastSI] ^= endArySI[onSI];
                  endArySI[onSI] ^= endArySI[lastSI];

                  scoreArySI[onSI] ^= scoreArySI[lastSI];
                  scoreArySI[lastSI] ^= scoreArySI[onSI];
                  scoreArySI[onSI] ^= scoreArySI[lastSI];

                  onSI = lastSI;
               } /*Loop: move swapped element back*/
            } /*If: I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSI = (subSI - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sortStart_demux*/

/*-------------------------------------------------------\
| Fun03: barcodeCoords_demux
|   - find the coordinates of barcodes in a read
| Input:
|   - lenSIPtr:
|     o signed int pointer to get length of returned array
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - barAryST:
|     o refST_kmerFind struct array with barcodes to demux
|   - barLenSI:
|     o number of barcodes in demux
|   - tblSTPtr:
|     o tblST_kmerFind struct pointer to use
|   - minPercScoreF:
|     o float with minimum percent score to keep a primer
|       mapping or discard
|     o recomend 0.9 (90%) to 0.95 (95%)
|     o do not go lower then 80% (otherwise to many
|       mappings kept)
|   - minScoreSL:
|     o signed long with minimum score to keep a primer
|       mapping or discard
|     o disable by using 0
|   - alnSetPtr:
|     o alnSet struct pointer with alignment settings
| Output:
|   - Modifies:
|     o lenSIPtr
|       - length of the returned array
|       - -1 if had two different barcodes overlapping
|       - -2 for memory errors
|   - Returns:
|     o array with the barcode index, start, end, and
|       score of each barcode
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + is negative if reverse, positive if forward
|       * ex: [index_1, start_1, end_1, score_1,
|              index_2, start_2 end_2, score_2,
|              index_3, start_3, end_3, score_3,
|              ...
|             ]
|     o for memory errors or no sequences (check lenSI)
\-------------------------------------------------------*/
signed int *
barcodeCoords_demux(
   signed int *lenSIPtr,    /*gets returned array length*/
   struct seqST *seqSTPtr,          /*has read to demux*/
   struct refST_kmerFind *barSTPtr, /*has barcodes*/
   signed int barLenSI,             /*number barcodes*/
   struct tblST_kmerFind *tblSTPtr, /*for searching*/
   float minPercScoreF,             /*min % score*/
   signed long minScoreSL,          /*min score to keep*/
   struct alnSet *alnSetPtr         /*settings*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - demux a read
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 Sec02:
   '     - detect primers
   '   o fun03 sec03:
   '     - filter barcodes and duplicate mappings
   '   o fun03 sec05:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *dirArySC = 0;
   signed short *primArySS = 0;
   signed int *scoreArySI = 0;
   signed int *seqStartArySI = 0;
   signed int *seqEndArySI = 0;
   signed short *primStartArySS = 0;
   signed short *primEndArySS = 0;
   signed int *retArySI = 0;
   signed int primSizeSI = 0;

   signed int lenSI = 0;
   signed int endSI = 0;
   signed int siMap = 0;
   signed int scoreSI = 0;
   signed int lastScoreSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - detect primers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenSIPtr = 0;
   lenSI =
      fxAllFindPrims_kmerFind(
         tblSTPtr,
         barSTPtr,
         barLenSI,
         seqSTPtr,
         minPercScoreF,
         &dirArySC,
         &primArySS,
         &scoreArySI,
         &seqStartArySI,
         &seqEndArySI,
         &primStartArySS,
         &primEndArySS,
         &primSizeSI,
         alnSetPtr
      );

   if(! lenSI)
      goto noPrim_fun03_sec05;
   else if(lenSI < 0)
      goto memErr_fun03_sec05;

   for(endSI = 0; endSI < lenSI; ++endSI)
      if( (dirArySC[endSI] & ~32) == 'R' )
         scoreArySI[endSI] *= -1; /*so have direction*/
   endSI = 0;

   free(dirArySC);
   dirArySC = 0;

   free(primStartArySS);
   primStartArySS = 0;

   free(primEndArySS);
   primEndArySS = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - filter barcodes and duplicate mappings
   ^   o fun03 sec03 sub01:
   ^     - filter by minimum scores
   ^   o fun03 sec03 sub02:
   ^     - flag duplicate mappings
   ^   o fun03 sec03 sub03:
   ^     - remove flagged duplicate mappings
   ^   o fun03 sec03 sub04:
   ^     - sort mappings by starting coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - filter by minimum scores
   \*****************************************************/

   endSI = lenSI;
   lenSI = 0;
   if(minScoreSL > 0)
   { /*If: have a min score to filter by*/
      for(siMap = 0; siMap < endSI; ++siMap)
      { /*Loop: remove low scores (non-percent check)*/
         if(scoreArySI[siMap] < 0)
            scoreSI = scoreArySI[siMap] * -1;
         else
            scoreSI = scoreArySI[siMap];

         if(scoreSI < minScoreSL)
            primArySS[siMap] = -1;
         else
            ++lenSI;
      } /*Loop: remove low scores (non-percent check)*/
   } /*If: have a min score to filter by*/

   else
      lenSI = endSI;

   sortBarcodes_demux(
      (unsigned short *) primArySS,
      seqStartArySI,
      seqEndArySI,
      scoreArySI,
      siMap
   ); /*remove filtered barcodes and group barcodes*/

   /*****************************************************\
   * Fun03 Sec03 Sub02:
   *   - flag duplicate mappings
   \*****************************************************/

   endSI = 0;
   lastScoreSI = 0;
   scoreSI = 0;
   for(siMap = 1; siMap < lenSI; ++siMap)
   { /*Loop: remove duplicate mappings*/
      if(seqStartArySI[siMap] < seqEndArySI[siMap-1])
      { /*If: barcode overlap or double map*/
         if(primArySS[siMap] != primArySS[siMap - 1])
            goto barcodeOverlap_fun03_sec05;

         if(scoreArySI[siMap] < 0)
            scoreSI = scoreArySI[siMap] * -1;
         else
            scoreSI = scoreArySI[siMap];

         if(scoreArySI[siMap - 1] < 0)
            lastScoreSI = scoreArySI[siMap - 1] * -1;
         else
            lastScoreSI = scoreArySI[siMap - 1];

         /*mapped the same position twice*/
         if(scoreSI >= lastScoreSI)
            primArySS[siMap - 1] = -1;
         else
            primArySS[siMap] = -1;

         if(siMap < 2)
            ++endSI; /*first comparision*/
      } /*If: barcode overlap or double map*/

      else
         ++endSI; /*make sure have at least one hit*/
   } /*Loop: remove duplicate mappings*/

   if(lenSI == 1)
      endSI = 1; /*no duplicate mappings*/
   else if(endSI <= 0)
      goto noPrim_fun03_sec05;

   /*****************************************************\
   * Fun03 Sec03 Sub03:
   *   - remove flagged duplicate mappings
   \*****************************************************/

   siMap = 0;
   while(siMap < lenSI)
   { /*Loop: remove duplicate mappings*/
      if(primArySS[siMap] < 0)
      { /*If: need to swap*/
         --lenSI;
         primArySS[siMap] ^= primArySS[lenSI];
         primArySS[lenSI] ^= primArySS[siMap];
         primArySS[siMap] ^= primArySS[lenSI];

         seqStartArySI[siMap] ^= seqStartArySI[lenSI];
         seqStartArySI[lenSI] ^= seqStartArySI[siMap];
         seqStartArySI[siMap] ^= seqStartArySI[lenSI];

         seqEndArySI[siMap] ^= seqEndArySI[lenSI];
         seqEndArySI[lenSI] ^= seqEndArySI[siMap];
         seqEndArySI[siMap] ^= seqEndArySI[lenSI];

         scoreArySI[siMap] ^= scoreArySI[lenSI];
         scoreArySI[lenSI] ^= scoreArySI[siMap];
         scoreArySI[siMap] ^= scoreArySI[lenSI];
      } /*If: need to swap*/

      else
         ++siMap;
   } /*Loop: remove duplicate mappings*/

   /*****************************************************\
   * Fun03 Sec03 Sub04:
   *   - sort mappings by starting coordinates
   \*****************************************************/

   sortStart_demux(
      primArySS,
      seqStartArySI,
      seqEndArySI,
      scoreArySI,
      siMap
   ); /*get positions*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - build the return array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retArySI = malloc((lenSI << 2) * sizeof(signed int));
   if(! retArySI)
      goto memErr_fun03_sec05;

   siMap = 0;
   endSI = lenSI << 2;
   lenSI = 0;

   while(siMap < endSI)
   { /*Loop: build the return array*/
      retArySI[siMap++] = primArySS[lenSI];
      retArySI[siMap++] = seqStartArySI[lenSI];
      retArySI[siMap++] = seqEndArySI[lenSI];
      retArySI[siMap++] = scoreArySI[lenSI];
      ++lenSI;
   } /*Loop: build the return array*/

   lenSI = siMap;

   free(primArySS);
   primArySS = 0;

   free(seqStartArySI);
   seqStartArySI = 0;

   free(seqEndArySI);
   seqEndArySI = 0;

   free(scoreArySI);
   scoreArySI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec05:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


   *lenSIPtr = lenSI;
   goto ret_fun03_sec05;
 
   noPrim_fun03_sec05:;
      goto errClean_fun03_sec05;

   barcodeOverlap_fun03_sec05:;
      lenSI = -1;
      goto errClean_fun03_sec05;

   memErr_fun03_sec05:;
      lenSI = -2;
      goto errClean_fun03_sec05;

   errClean_fun03_sec05:;
      if(retArySI)
         free(retArySI);
      retArySI = 0;

      goto ret_fun03_sec05;

   ret_fun03_sec05:;      
      if(dirArySC)
         free(dirArySC);
      dirArySC = 0;

      if(primArySS)
         free(primArySS);
      primArySS = 0;

      if(scoreArySI)
         free(scoreArySI);
      scoreArySI = 0;

      if(seqStartArySI)
         free(seqStartArySI);
      seqStartArySI = 0;

      if(seqEndArySI)
         free(seqEndArySI);
      seqEndArySI = 0;

      if(primStartArySS)
         free(primStartArySS);
      primStartArySS = 0;

      if(primEndArySS)
         free(primEndArySS);
      primEndArySS = 0;

      return retArySI;
} /*barcodeCoords_demux*/

/*-------------------------------------------------------\
| Fun04: read_demux
|   - convert barcode coordinates to demuxed reads
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - maxSplitSI:
|     o maximum number of splits before tossing read
|   - minDistSI:
|     o minimum distance between barcodes
|   - trimBl:
|     o 1: trim barcodes
|     o 0: do not trim barcodes
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + a negative score means a reverse mapping
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - barFileStrAry:
|     o c-string array with output file names
|     o you need to check and make sure these files are
|       writeable before hand (I assume here)
| Output:
|   - Prints:
|     o demuxed and if requested split sequence to the
|       correct file in barFileStrAry
|   - Returns:
|     o 0 for no errors
|     o 1 if read had no barcodes
|     o 2 if read had to many barcodes
|     o 3 if read had barcodes at both ends; no easy split
|     o 4 if read had to barcodes right next to each other
\-------------------------------------------------------*/
signed int
read_demux(
   struct seqST *seqSTPtr,   /*read to demux*/
   signed int maxSplitSI,    /*max splits before discard*/
   signed int minDistSI, /*min distance between barcodes*/
   signed char trimBl,       /*1: trim barcodes*/
   signed int coordArySI[],  /*has barcode mappings*/
   signed int coordLenSI,    /*length of coordLenSI*/
   signed char *barFileStrAry[]/*output file names*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - convert barcode coordinates to demuxed reads
   '   o fun04 sec01:
   '     - variable declarations
   '   o fun04 sec02:
   '     - check if to many splits and print no split case
   '   o fun04 sec03:
   '     - check if splits are to close and ending barcode
   '   o fun04 sec04:
   '     - print split sequences
   '   o fun04 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char tmpSC = 0;
   signed char tmpQSC = 0;
   signed char atEndBl = 0;/*first barcode at end of seq*/
   signed int posSI = 0;
   signed char tmpIdStr[4096];
   signed int headSI = 0;

   signed int scoreSI = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - check if to many splits and print no split case
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(coordLenSI < 4)
      goto noBarcodes_fun04_sec05;
   else if(coordLenSI / 4 > maxSplitSI)
      goto toManyBarcodes_fun04_sec05;

   if(coordLenSI < 5)
   { /*If: only have one barcode*/
      outFILE =
         fopen(
            (char *) barFileStrAry[coordArySI[0]],
            "a"
         );
         /*user should have checked this before hand*/

      if(trimBl)
      { /*If: trimming barcodes*/
         if(coordArySI[1] > seqSTPtr->seqLenSL / 2)
         { /*If: barcode at end*/
            tmpSC = seqSTPtr->seqStr[coordArySI[1]];
            seqSTPtr->seqStr[coordArySI[1]] = 0;
            posSI = 0;
         } /*If: barcode at end*/

         else
            posSI = coordArySI[2] + 1; /*at start*/
      } /*If: trimming barcodes*/

      if(coordArySI[3] < 0)
         scoreSI = coordArySI[3] * -1;
      else
         scoreSI = coordArySI[3];

      fprintf(
         outFILE,
         "@%s\tscore=%i\tstart=%i\tend=%i%s%s%s+%s%s%s",
         seqSTPtr->idStr,
         scoreSI,
         coordArySI[1] + 1,
         coordArySI[2] + 1,
         str_endLine,
         &seqSTPtr->seqStr[posSI],
         str_endLine,
         str_endLine,
         seqSTPtr->qStr,
         str_endLine
      );

      if(tmpSC)
         seqSTPtr->seqStr[coordArySI[1]] = tmpSC;
 
      fclose(outFILE);
      outFILE = 0;
      goto done_fun04_sec05;
   } /*If: only have one barcode*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec03:
   ^   - check if splits are to close and ending barcode
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(coordArySI[1] > minDistSI)
   { /*If: nearby barcode is not near the start*/
      if(
          seqSTPtr->seqLenSL - coordArySI[coordLenSI - 1]
        > minDistSI
      ) goto barcodesAtEnds_fun04_sec05;

      atEndBl = 1;
   } /*If: nearby barcode is not near the start*/

   else if(
        seqSTPtr->seqLenSL - coordArySI[coordLenSI - 1]
      < minDistSI
   ) goto barcodesAtEnds_fun04_sec05;
     /*ending barcodes are at the ends (no way to tell
     `  were to split
     */
   else
      atEndBl = 0;

   maxSplitSI = 4;
   while(maxSplitSI < coordLenSI)
   { /*Loop: check for odd errors*/
      if(
           coordArySI[maxSplitSI + 1] /*this start*/
         - coordArySI[maxSplitSI - 2] /*last end*/
         < minDistSI
      ) goto barcodesToClose_fun04_sec05;
        /*barcodes are to close together*/
      maxSplitSI += 4; /*move to next barcode map*/
   } /*Loop: check for odd errors*/

   headSI = 0;
   while(seqSTPtr->idStr[headSI] > 32)
   { /*Loop: copy id until no white space*/
      tmpIdStr[headSI] = seqSTPtr->idStr[headSI];
      ++headSI;
   } /*Loop: copy id until no white space*/

   tmpIdStr[headSI] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec04:
   ^   - print split sequences
   ^   o fun04 sec04 sub01:
   ^     - print cases were starting barcode is at end
   ^   o fun04 sec04 sub02:
   ^     - print cases were starting barcode is at start
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec04 Sub01:
   *   - print cases were starting barcode is at end
   \*****************************************************/

   if(atEndBl)
   { /*If: working from end to start*/
      for(
         maxSplitSI = coordLenSI - 1;
         maxSplitSI >= 0;
         maxSplitSI -= 4
      ){ /*Loop: print out the split sequences*/
         outFILE =
            fopen(
               (char *)
                  barFileStrAry[coordArySI[maxSplitSI-2]],
               "a"
            );
         posSI = coordArySI[maxSplitSI - 5]; /*last end*/

         if(trimBl)
         { /*If: trimming barcodes*/
            if(maxSplitSI < coordLenSI - 2)
            { /*If: not on last sequence (first loop)*/
               tmpSC =
                 seqSTPtr->seqStr[
                    coordArySI[maxSplitSI-2]/*this start*/
                 ];
               seqSTPtr->seqStr[
                  coordArySI[maxSplitSI - 2]/*this start*/
               ] = 0;

               tmpQSC =
                 seqSTPtr->qStr[
                    coordArySI[maxSplitSI-2]/*this start*/
                 ];
               seqSTPtr->qStr[
                  coordArySI[maxSplitSI - 2]/*this start*/
               ] = 0;
            } /*If: not on last sequence (first loop)*/
         } /*If: trimming barcodes*/

         else
         { /*Else: need to get barcode end*/
            tmpSC =
              seqSTPtr->seqStr[coordArySI[maxSplitSI-3]];
            seqSTPtr->seqStr[coordArySI[maxSplitSI-3]] =0;

            tmpQSC =
              seqSTPtr->qStr[coordArySI[maxSplitSI - 3]];
            seqSTPtr->qStr[coordArySI[maxSplitSI - 3]] =0;
         } /*Else: need to get barcode end*/

         fprintf(outFILE,"@%s-%i",tmpIdStr,maxSplitSI>>2);

         fprintf(
            outFILE,
            "%s\tscore=%i\tstart=%i\tend=%i\t%s%s%s",
            &seqSTPtr->idStr[headSI],
            coordArySI[maxSplitSI],
            coordArySI[maxSplitSI - 2] + 1,
            coordArySI[maxSplitSI - 3] + 1,
            str_endLine,
            &seqSTPtr->seqStr[posSI],
            str_endLine
         );

         fprintf(
            outFILE,
            "+%s%s%s",
            str_endLine,
            &seqSTPtr->qStr[posSI],
            str_endLine
         );

         if(tmpSC)
         { /*If: need to trim and did trim*/
            seqSTPtr->seqStr[coordArySI[maxSplitSI - 2]] =
               tmpSC;
            seqSTPtr->qStr[coordArySI[maxSplitSI - 2]] =
               tmpQSC;
            tmpSC = 0;
            tmpQSC = 0;
         } /*If: need to trim and did trim*/

         fclose(outFILE);
         outFILE = 0;
      } /*Loop: print out the split sequences*/
   } /*If: working from end to start*/

   /*****************************************************\
   * Fun04 Sec04 Sub02:
   *   - print cases were starting barcode is at start
   \*****************************************************/

   else
   { /*Else: working from start to end*/
      for(
         maxSplitSI = 0;
         maxSplitSI < coordLenSI;
         maxSplitSI += 4
      ){ /*Loop: print out the split sequences*/
         outFILE =
            fopen(
               (char *)
                  barFileStrAry[coordArySI[maxSplitSI]],
               "a"
            );

         if(trimBl)
         { /*If: trimming barcodes*/
            posSI = coordArySI[maxSplitSI + 2]; /*end*/

         } /*If: trimming barcodes*/

         else
         { /*Else: not trimming*/
            if(! maxSplitSI)
               posSI = 0;
            else
               posSI = coordArySI[maxSplitSI + 1];
         } /*Else: not trimming*/

         if(maxSplitSI < coordLenSI - 4)
         { /*If: not on the last barcode*/
            /*set start of next barcode to 0*/
            tmpSC =
               seqSTPtr->seqStr[
                  coordArySI[maxSplitSI + 5]
               ];
            seqSTPtr->seqStr[
               coordArySI[maxSplitSI + 5]
            ] = 0;

            tmpQSC =
               seqSTPtr->qStr[
                  coordArySI[maxSplitSI + 5]
               ];
            seqSTPtr->qStr[
               coordArySI[maxSplitSI + 5]
            ] = 0;
         } /*If: not on the last barcode*/

         fprintf(outFILE,"@%s-%i",tmpIdStr,maxSplitSI>>2);

         fprintf(
            outFILE,
            "%s\tscore=%i\tstart=%i\tend=%i%s%s%s",
            &seqSTPtr->idStr[headSI],
            coordArySI[maxSplitSI + 3],
            coordArySI[maxSplitSI + 1] + 1,
            coordArySI[maxSplitSI + 2] + 1,
            str_endLine,
            &seqSTPtr->seqStr[posSI],
            str_endLine
         );

         fprintf(
            outFILE,
            "+%s%s%s",
            str_endLine,
            &seqSTPtr->qStr[posSI],
            str_endLine
         );

         if(tmpSC)
         { /*If: need to trim and did trim*/
            seqSTPtr->seqStr[
               coordArySI[maxSplitSI + 5]
            ] = tmpSC;
            seqSTPtr->qStr[
               coordArySI[maxSplitSI + 5]
            ] = tmpQSC;

            tmpSC = 0;
            tmpQSC = 0;
         } /*If: need to trim and did trim*/

         fclose(outFILE);
         outFILE = 0;
      } /*Loop: print out the split sequences*/
   } /*Else: working from start to end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   done_fun04_sec05:;
      posSI = 0;
      goto ret_fun04_sec05;

   noBarcodes_fun04_sec05:;
      posSI = 1;
      goto ret_fun04_sec05;

   toManyBarcodes_fun04_sec05:;
      posSI = 2;
      goto ret_fun04_sec05;

   barcodesAtEnds_fun04_sec05:;
      posSI = 3;
      goto ret_fun04_sec05;

   barcodesToClose_fun04_sec05:;
      posSI = 4;
      goto ret_fun04_sec05;

   ret_fun04_sec05:;
      if(outFILE)
         fclose(outFILE);
      outFILE = 0;

      return posSI;
} /*read_demux*/

/*-------------------------------------------------------\
| Fun05: primer_demux
|   - get primer target regins from the input sequence
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - minDistSI:
|     o minimum distance between primers
|   - maxDistSI:
|     o maximum distance between primers
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + a negative score means a reverse mapping
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - barSTPtr:
|     o refST_kmerFind struct pionter with the primer
|       ids and the index of their mates (if paired)
|   - outFILE:
|     o FILE pionter to print reads to
| Output:
|   - Prints:
|     o amplicons to outFILE
|   - Returns:
|     o number of amplicons found
|     o 0 if no amplicons
\-------------------------------------------------------*/
signed int
primer_demux(
   struct seqST *seqSTPtr,/*read to split into amplicons*/
   signed int minDistSI,  /*min distance between primers*/
   signed int maxDistSI,  /*max distance between primers*/
   signed int coordArySI[],/*has barcode mappings*/
   signed int coordLenSI, /*length of coordLenSI*/
   struct refST_kmerFind *barSTPtr,/*primer ids & mates*/
   void *outFILE          /*print sequences to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - get primer target regins from the input sequence
   '   o fun05 sec01:
   '     - variable declarations
   '   o fun05 sec02:
   '     - check if to many splits and print no split case
   '   o fun05 sec03:
   '     - check if splits are to close and ending barcode
   '   o fun05 sec04:
   '     - print split sequences
   '   o fun05 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char tmpSeqSC = 0; /*seting sequence end to 0*/
   signed char tmpQSC = 0;/*for setting q-score end to 0*/

   signed int startSI = 0; /*first base in amplicon*/
   signed int endSI = 0;   /*last base in amplicon*/
   signed int lenSI = 0;   /*amplicon length*/
   signed int cntSI = 0;   /*number primers found*/
   signed int siCoord = 0; /*primer on*/

   signed int barSI = 0;  /*id of primer on*/
   signed int mateSI = 0; /*id of the mate primer*/
   signed int siNext = 0;/*finding next primer in a pair*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - check if to many splits and print no split case
   ^   o fun05 sec02 sub01:
   ^     - check if have primers and start primer loop
   ^   o fun05 sec02 sub02:
   ^     - deal with primer pairs
   ^   o fun05 sec02 sub03:
   ^     - non-paired primers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec02 Sub01:
   *   - check if have primers and start primer loop
   \*****************************************************/

   if(coordLenSI < 4)
      goto noBarcodes_fun05_sec05;

   for(siCoord = 0; siCoord < coordLenSI; siCoord += 4)
   { /*Loop: print primers*/
      barSI = coordArySI[siCoord];
      mateSI = barSTPtr[barSI].mateSI;

      /**************************************************\
      * Fun05 Sec02 Sub02:
      *   - deal with primer pairs
      *   o fun05 sec02 sub02 cat01:
      *     - find the mate for the current primer
      *   o fun05 sec02 sub02 cat02:
      *     - print the sequence from primer start to
      *       the mates primer end
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun05 Sec02 Sub02 Cat01:
      +   - find the mate for the current primer
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(mateSI >= 0)
      { /*If: have a reverse barcode*/
         if(coordArySI[siCoord + 3] < 0)
            continue;
            /*should have already found the mate, since
            `   this requires searching backwards
            */
         for(
            siNext = siCoord + 4;
            siNext < coordLenSI;
            siNext += 4
         ){ /*Loop: find the mate primer*/
            if(coordArySI[siNext] != mateSI)
               ; /*not the mate of the primer pair*/
            else if(coordArySI[siNext + 3] > 0)
               ; /*wrong direction to amplify*/
            else if(coordArySI[siNext + 3] < 0)
               break;
         }  /*Loop: find the mate primer*/

         if(siNext >= coordLenSI)
            continue;

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun05 Sec02 Sub02 Cat02:
         +   - print the sequence from primer start to
         +     the mates primer end
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         lenSI =
           coordArySI[siNext +2] - coordArySI[siCoord +1];
         ++lenSI;
         startSI = coordArySI[siCoord + 1];
         endSI = coordArySI[siNext + 2];

         if(lenSI < minDistSI)
            continue;
         else if(lenSI > maxDistSI && maxDistSI)
            continue;
         else
         { /*Else: have match with a good length*/
            tmpSeqSC =
               seqSTPtr->seqStr[coordArySI[siNext +2] +1];
            seqSTPtr->seqStr[coordArySI[siNext +2]+1] = 0;

            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
            { /*If: have quality score entry*/
               tmpQSC =
                  seqSTPtr->qStr[coordArySI[siNext +2]+1];
               seqSTPtr->qStr[coordArySI[siNext +2] +1] =
                  0;
               fprintf(
                  (FILE *) outFILE,
                  "@%s-%03i",
                  seqSTPtr->idStr,
                  cntSI++
               );
            } /*If: have quality score entry*/

            else
               fprintf(
                  (FILE *) outFILE,
                  ">%s-%03i",
                  seqSTPtr->idStr,
                  cntSI++
               );

            if(mateSI < barSI)
            { /*If: reverse complemet*/
               fprintf(
                  (FILE *) outFILE,
                  "\tdir=rev\tlen=%i\tstart=%i\tend=%i",
                  lenSI,
                  startSI + 1,
                  endSI + 1
               );
               fprintf(
                  (FILE *) outFILE,
                  "\tforPrim=%s\trevPrim=%s",
                  barSTPtr[mateSI].forSeqST->idStr,
                  barSTPtr[barSI].forSeqST->idStr
               );
               fprintf(
                  (FILE *) outFILE,
                  "\tforStart=%i\tforEnd=%i\tforScore=%i",
                  coordArySI[siNext + 1] - startSI + 1,
                  coordArySI[siNext + 2] - startSI + 1,
                  coordArySI[siNext + 3] * -1
               );
               fprintf(
                  (FILE *) outFILE,
                  "\trevStart=%i\trevEnd=%i\trevScore=%i",
                  coordArySI[siCoord + 1] - startSI + 1,
                  coordArySI[siCoord + 2] - startSI + 1,
                  coordArySI[siCoord + 3]
               );
            } /*If: reverse complemet*/
               
            else
            { /*Else: forward*/
               fprintf(
                  (FILE *) outFILE,
                  "\tdir=for\tlen=%i\tstart=%i\tend=%i",
                  lenSI,
                  startSI + 1,
                  endSI + 1
               );
               fprintf(
                  (FILE *) outFILE,
                  "\tforPrim=%s\trevPrim=%s",
                  barSTPtr[barSI].forSeqST->idStr,
                  barSTPtr[mateSI].forSeqST->idStr
               );
               fprintf(
                  (FILE *) outFILE,
                  "\tforStart=%i\tforEnd=%i\tforScore=%i",
                  coordArySI[siCoord + 1] - startSI + 1,
                  coordArySI[siCoord + 2] - startSI + 1,
                  coordArySI[siCoord + 3]
               );
               fprintf(
                  (FILE *) outFILE,
                  "\trevStart=%i\trevEnd=%i\trevScore=%i",
                  coordArySI[siNext + 1] - startSI + 1,
                  coordArySI[siNext + 2] - startSI + 1,
                  coordArySI[siNext + 3] * -1
               );
            } /*Else: forward*/

            if(! seqSTPtr->qStr || ! seqSTPtr->qStr[0])
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[coordArySI[siCoord+1]],
                 str_endLine
               );
            else
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s+%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[coordArySI[siCoord+1]],
                 str_endLine, /*end of sequence*/
                 str_endLine, /*for +*/
                 &seqSTPtr->qStr[coordArySI[siCoord + 1]],
                 str_endLine  /*for q-score entry*/
               );
               
            seqSTPtr->seqStr[coordArySI[siNext + 2] + 1] =
               tmpSeqSC;

            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
               seqSTPtr->qStr[coordArySI[siNext + 2] + 1]=
                  tmpQSC;
         } /*Else: have match with a good length*/
      } /*If: have a reverse barcode*/

      /**************************************************\
      * Fun05 Sec02 Sub03:
      *   - non-paired primers
      *   o fun05 sec02 sub03 cat01:
      *     - non-paired else + deal with reverse primers
      *   o fun05 sec02 sub03 cat02:
      *     - deal with single forward primers
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun05 Sec02 Sub03 Cat01:
      +   - non-paired else + deal with reverse primers
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: solo barcode*/
         if(coordArySI[siCoord + 3] < 0)
         { /*If: need to go backwards*/
            lenSI = coordArySI[siCoord + 2];
            endSI = lenSI;

            if(siCoord)
            { /*If: have a pervious primer*/
               lenSI -= coordArySI[siCoord - 2];
               ++lenSI;
               startSI = coordArySI[siCoord - 2];

               if(coordArySI[siCoord - 1] >= 0)
                  continue; /*the last and current
                            `  primers are in oposite
                            `  directions and so overlap,
                            `  so I can not make a call
                            */
            } /*If: have a pervious primer*/

            else
               startSI = 0;

            if(lenSI < minDistSI)
               continue; /*to short*/
            else if(lenSI > maxDistSI && maxDistSI)
               continue; /*to long*/

            tmpSeqSC =
               seqSTPtr->seqStr[coordArySI[siCoord +2]+1];
            seqSTPtr->seqStr[coordArySI[siCoord + 2]+1]=0;

            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
            { /*If: have quality score entry*/
               tmpQSC =
                  seqSTPtr->qStr[coordArySI[siCoord+2]+1];
               seqSTPtr->qStr[coordArySI[siCoord+2]+1]=0;

               fprintf(
                  (FILE *) outFILE,
                  "@%s-%03i\tdir=rev\tlen=%i\t",
                  seqSTPtr->idStr,
                  cntSI++,
                  lenSI
               );
            } /*If: have quality score entry*/

            else
               fprintf(
                  (FILE *) outFILE,
                  ">%s-%03i\tdir=rev\tlen=%i\t",
                  seqSTPtr->idStr,
                  cntSI++,
                  lenSI
               );

            fprintf(
               (FILE *) outFILE,
               "start=%i\tend=%i\tforPrim=NA\trevPrim=%s",
               startSI + 1,
               endSI + 1,
               barSTPtr[barSI].forSeqST->idStr
            );
            fprintf(
               (FILE *) outFILE,
               "\tforStart=NA\tforEnd=NA\tforScore=NA"
            );
            fprintf(
               (FILE *) outFILE,
               "\trevStart=%i\trevEnd=%i\trevScore=%i",
               coordArySI[siCoord + 1] - startSI + 1,
               coordArySI[siCoord + 2] - startSI + 1,
               coordArySI[siCoord + 3] * -1
            );

            if(! seqSTPtr->qStr || ! seqSTPtr->qStr[0])
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[
                    coordArySI[siCoord - 2] + 1
                 ],
                 str_endLine
               );
            else
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s+%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[
                    coordArySI[siCoord - 2] + 1
                 ],
                 str_endLine, /*end of sequence*/
                 str_endLine, /*for +*/
                 &seqSTPtr->qStr[coordArySI[barSI + 1]],
                 str_endLine  /*for q-score entry*/
               );

            seqSTPtr->seqStr[coordArySI[siCoord + 2]+1] =
               tmpSeqSC;
            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
               seqSTPtr->qStr[coordArySI[siCoord + 2]+1] =
                  tmpQSC;
         } /*If: need to go backwards*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun05 Sec02 Sub03 Cat02:
         +   - deal with single forward primers
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: forward*/
            lenSI = coordArySI[siCoord + 1];

            if(siCoord < coordLenSI - 4)
            { /*If: have a next primer*/
               lenSI = coordArySI[siCoord + 6] - lenSI +1;

               if(coordArySI[siCoord + 7] >= 0)
                  continue; /*the next primer and this
                            `  primers are in oposite
                            `  directions and so you can
                            `  not make a call
                            */
            } /*If: have a next primer*/

            else
               lenSI = seqSTPtr->seqLenSL - lenSI;

            if(lenSI < minDistSI)
               continue; /*to short*/
            else if(lenSI > maxDistSI && maxDistSI)
               continue; /*to long*/

            tmpSeqSC =
               seqSTPtr->seqStr[coordArySI[siCoord +6]+1];
            seqSTPtr->seqStr[coordArySI[siCoord + 6]+1]=0;
               /*6 = 4 + 2; 4 moves to next primer and
               `  2 is end of the primer
               */

            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
            { /*If: have a quality score entry*/
               tmpQSC =
                  seqSTPtr->qStr[coordArySI[siCoord+6]+1];
               seqSTPtr->qStr[coordArySI[siCoord + 6]+1] =
                  0;

               fprintf(
                  (FILE *) outFILE,
                  ">%s-%03i\tdir=for\tlen=%i\t",
                  seqSTPtr->idStr,
                  cntSI++,
                  lenSI
               );
            } /*If: have a quality score entry*/

            else
            { /*Else: is a fastq file*/
               fprintf(
                  (FILE *) outFILE,
                  "@%s-%03i\tdir=for\tlen=%i\t",
                  seqSTPtr->idStr,
                  cntSI++,
                  lenSI
               );
            } /*Else: is a fastq file*/

            fprintf(
               (FILE *) outFILE,
               "start=%i\tend=%i\tforPrim=%s\trevPrim=NA",
               startSI + 1,
               endSI + 1,
               barSTPtr[siNext].forSeqST->idStr
            );
            fprintf(
               (FILE *) outFILE,
               "\tforStart=%i\tforEnd=%i\tforScore=%i",
               coordArySI[siCoord + 1] - startSI + 1,
               coordArySI[siCoord + 2] - startSI + 1,
               coordArySI[siCoord + 3]
            );
            fprintf(
               (FILE *) outFILE,
               "\trevStart=NA\trevEnd=NA\trevScore=NA"
            );

            if(! seqSTPtr->qStr || ! seqSTPtr->qStr[0])
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[coordArySI[siCoord+1]],
                 str_endLine
               );
            else
               fprintf(
                 (FILE *) outFILE,
                 "%s%s%s+%s%s%s",
                 str_endLine, /*for end of header*/
                 &seqSTPtr->seqStr[coordArySI[siCoord+1]],
                 str_endLine, /*end of sequence*/
                 str_endLine, /*for +*/
                 &seqSTPtr->qStr[coordArySI[barSI + 1]],
                 str_endLine  /*for q-score entry*/
               );

            seqSTPtr->seqStr[coordArySI[siCoord + 6]+1] =
               tmpSeqSC;

            if(seqSTPtr->qStr && seqSTPtr->qStr[0])
               seqSTPtr->qStr[coordArySI[siCoord + 6]+1] =
                  tmpQSC;
         } /*Else: forward*/
      } /*Else: solo barcode*/
   } /*Loop: print primers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun05_sec05;

   noBarcodes_fun05_sec05:;
      cntSI = 0;
      goto ret_fun05_sec05;

   ret_fun05_sec05:;
      return cntSI;
} /*primer_demux*/

/*-------------------------------------------------------\
| Fun06: pGeneCoord_demux
|   - get coordinates of all genes found in target
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to find genes
|       in
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + a negative score means a reverse mapping
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - geneSTPtr:
|     o refST_kmerFind struct pionter with the genes
|   - headBlPtr:
|     o 1: print header and then set to 0
|     o 0: do not print the header
|   - geneSTPtr
|     o refST_kmerFind struct pionter with the genes
|       to find
|   - outFILE:
|     o FILE pionter to print gene coordinates
| Output:
|   - Prints:
|     o gene coordinates to outFILE
|   - Returns:
|     o number of amplicons found
|     o 0 if no genes
\-------------------------------------------------------*/
signed int
pGeneCoord_demux(
   struct seqST *seqSTPtr,/*read to split into amplicons*/
   signed int coordArySI[],/*has barcode mappings*/
   signed int coordLenSI, /*length of coordLenSI*/
   signed char *headBlPtr,/*1: print header + set to 0*/
   struct refST_kmerFind *geneSTPtr,/*genes searched*/
   void *outFILE          /*print sequences to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC:
   '   - get coordinates of all genes found in target
   '   o fun06 sec01:
   '     - variable declarations
   '   o fun06 sec02:
   '     - print header if needed
   '   o fun06 sec03:
   '     - print gene coordinates
   '   o fun06 sec04:
   '     - print split sequences
   '   o fun06 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int lenSI = 0;   /*amplicon length*/
   signed int cntSI = 0;   /*number primers found*/
   signed int siCoord = 0; /*primer on*/

   signed int geneSI = 0; /*index of gene on*/

   signed int *coordHeapArySI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - print header if needed
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(coordLenSI < 4)
      goto noGenes_fun06_sec04;

   if(*headBlPtr)
   { /*If: printing the header*/
      fprintf(
         (FILE *) outFILE,
         "id-count\tgene\tdir\tscore\tmaxScore\tlength"
      );
      fprintf(
         (FILE *) outFILE,
         "\tmaxLength\tstart\tend%s",
         str_endLine
      );
      *headBlPtr = 0;
   } /*If: printing the header*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^   - print gene coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(siCoord = 0; siCoord < coordLenSI; siCoord += 4)
   { /*Loop: print gene coordinates*/
      geneSI = coordArySI[siCoord];
      lenSI =
          coordArySI[siCoord+2] - coordArySI[siCoord+1]+1;
      if(coordArySI[siCoord + 3] < 0)
         fprintf(
           (FILE *) outFILE,
           "%s-%03i\t%s\tR\t%i\t%0.2f\t%i\t%li\t%i\t%i%s",
           seqSTPtr->idStr,
           cntSI++,                 /*number of genes*/
           geneSTPtr[geneSI].forSeqST->idStr,
           coordArySI[siCoord + 3], /*score*/
           geneSTPtr[geneSI].maxRevScoreF,
           lenSI,
           geneSTPtr[geneSI].forSeqST->seqLenSL,
           coordArySI[siCoord + 1] + 1,     /*start*/
           coordArySI[siCoord + 2] + 1,     /*end*/
           str_endLine
         );
      else
         fprintf(
           (FILE *) outFILE,
           "%s-%03i\t%s\tF\t%i\t%0.2f\t%i\t%li\t%i\t%i%s",
           seqSTPtr->idStr,
           cntSI++,                       /*number genes*/
           geneSTPtr[geneSI].forSeqST->idStr,
           coordArySI[siCoord + 3],   /*score*/
           geneSTPtr[geneSI].maxForScoreF,
           lenSI,
           geneSTPtr[geneSI].forSeqST->seqLenSL,
           coordArySI[siCoord +1] +1, /*start*/
           coordArySI[siCoord +2] +1, /*end*/
           str_endLine
         );
   } /*Loop: print gene coordinates*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun06_sec04;

   noGenes_fun06_sec04:;
      cntSI = 0;
      goto ret_fun06_sec04;

   ret_fun06_sec04:;
      if(coordHeapArySI)
         free(coordHeapArySI);
      coordHeapArySI = 0;

      return cntSI;
} /*pGeneCoord_demux*/

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
: Copyright (c) 2025 jeremyButtler
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
