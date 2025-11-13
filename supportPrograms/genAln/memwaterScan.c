/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memwaterScan SOF: Start Of Functions
'   - has memory efficent Waterman Smith alignment scan
'     (best alignment for each reference/query base)
'   o header:
'     - included libraries
'   o .h st01: aln_memwaterScan
'     - stores the alignment results from a waterman scan
'   o fun01: blank_aln_memwaterScan
'     - blanks (sets non-array values to 0) an
'   o fun02: init_aln_memwaterScan
'     - initializes (sets arrays to 0 + blanks) an
'       aln_memwaterScan structure
'   o fun03: freeStack_aln_memwaterScan
'     - frees the variables in an aln_memwaterScan struct
'   o fun04: freeHeap_aln_memwaterScan
'     - frees an aln_memwaterScan struct
'   o fun05: getCoord_aln_memwaterScan
'     - gets the coordinates for input base index for an
'       aln_memwaterScan structure
'   o fun06: swap_memwaterScan
'     - swaps position in reference array
'   o fun07: refCoordSort_aln_memwaterScan
'     - sorts best bases alignments by reference start,
'       query start, reference end (closest first), then
'       and query end (closest first)
'     - all negatives (non-alignment found) are pushed to
'       the end
'   o fun08: filter_memwaterScan
'     - removes low scoring alignments and alignments that
'       are nested alignments
'   o fun09: merge_aln_memwaterScan
'     - merges one aln_memwaterScan struct into another
'       aln_memwaterScan struct
'   o fun10 memwaterScan:
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference base)
'       alignment on a pair of sequences
'   o fun11: simple_memwaterScan
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference
'       base) alignment on a pair of sequences
'     - simple means no settings or seqST structure used
'       and no match matrix
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

#include "memwaterScan.h"

#include "../genBio/seqST.h"

#include "alnSet.h"
#include "indexToCoord.h"

/*.h files only*/
#include "../genLib/genMath.h" /*only using .h portion*/
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o std #include <stdio.h>
!   o .c  #include "../genLib/base10str.h"
!   o .c  #include "../genLib/ulCp.h"
!   o .c  #include "../genLib/fileFun.h"
!   o .h  #include "../genLib/endLine.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_aln_memwaterScan
|   - blanks (sets non-array values to 0) an
|     aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to blank
| Output:
|   - Modifies:
|     o sets refLenSI, qryLenSI, refOffsetSI, qryOffsetSI
|       in alnsTPtr to 0
\-------------------------------------------------------*/
void
blank_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){
   if(! alnSTPtr)
      return;

   alnSTPtr->refLenSI = 0;
   alnSTPtr->qryLenSI = 0;
   alnSTPtr->refOffsetSI = 0;
   alnSTPtr->qryOffsetSI = 0;

   alnSTPtr->outLenSL = 0;
} /*blank_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun02: init_aln_memwaterScan
|   - initializes (sets arrays to 0 + blanks) an
|     aln_memwaterScan structure
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to initialize
| Output:
|   - Modifies:
|     o sets all values to 0
\-------------------------------------------------------*/
void
init_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){
   if(! alnSTPtr)
      return;

   alnSTPtr->startArySL = 0;
   alnSTPtr->endArySL = 0;
   alnSTPtr->scoreArySL = 0;
   alnSTPtr->outSizeSL = 0;

   alnSTPtr->scoreRowSL = 0;
   alnSTPtr->indexRowSL = 0;
   alnSTPtr->dirRowSC = 0;
   alnSTPtr->rowSizeSI = 0;

   blank_aln_memwaterScan(alnSTPtr);
} /*init_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun03: freeStack_aln_memwaterScan
|   - frees the variables in an aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer with variable to
|       free
| Output:
|   - Modifies:
|     o frees all arrays and initializes
\-------------------------------------------------------*/
void
freeStack_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){
   if(! alnSTPtr)
      return;

   if(alnSTPtr->startArySL)
      free(alnSTPtr->startArySL);
   if(alnSTPtr->endArySL)
      free(alnSTPtr->endArySL);
   if(alnSTPtr->scoreArySL)
      free(alnSTPtr->scoreArySL);

   if(alnSTPtr->scoreRowSL)
      free(alnSTPtr->scoreRowSL);
   if(alnSTPtr->indexRowSL)
      free(alnSTPtr->indexRowSL);
   if(alnSTPtr->dirRowSC)
      free(alnSTPtr->dirRowSC);

   init_aln_memwaterScan(alnSTPtr);
} /*freeStack_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun04: freeHeap_aln_memwaterScan
|   - frees an aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to free
| Output:
|   - Modifies:
|     o frees alnSTPtr; you must set the pointer to 0/null
\-------------------------------------------------------*/
void
freeHeap_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){
   if(! alnSTPtr)
      return;
   freeStack_aln_memwaterScan(alnSTPtr);
   free(alnSTPtr);
} /*freeHeap_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun05: getCoord_aln_memwaterScan
|   - gets the coordinates for input base index for an
|     aln_memwaterScan structure
| Input:
|   - baseSI:
|     o base to get coordinates for
|   - refStartSIPtr:
|     o signed int pointer that gets the reference start
|   - refEndSIPtr:
|     o signed int pointer that gets the reference end
|   - qryStartSIPtr:
|     o signed int pointer that gets the query start
|   - qryEndSIPtr:
|     o signed int pointer that gets the query end
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer with coordinates
| Output:
|   - Modifies:
|     o refStartSIPtr to have reference starting position
|     o refEndSIPtr to have reference ending position
|     o qryStartSIPtr to have query starting position
|     o qryEndSIPtr to have query ending position
|   - Returns:
|     o score for success
|     o 0 if no position at index
|     o -1 for out of bounds error
\-------------------------------------------------------*/
signed char
getCoord_aln_memwaterScan(
   signed long baseSL,        /*index to grab*/
   signed int *refStartSIPtr, /*gets reference start*/
   signed int *refEndSIPtr,   /*gets reference end*/
   signed int *qryStartSIPtr, /*gets query start*/
   signed int *qryEndSIPtr,   /*gets query end*/
   struct aln_memwaterScan *alnSTPtr
){
   unsigned long tmpRefUL = 0;
   unsigned long tmpQryUL = 0;

   if(baseSL > alnSTPtr->outLenSL)
      goto boundsErr_fun05_sec04;
   else if(alnSTPtr->scoreArySL[baseSL] < 0)
      goto noCoord_fun05_sec04;

   indexToCoord(
      alnSTPtr->refLenSI,
      alnSTPtr->startArySL[baseSL],
      &tmpRefUL,
      &tmpQryUL
   );

   *refStartSIPtr = tmpRefUL + alnSTPtr->refOffsetSI;
   *qryStartSIPtr = tmpQryUL + alnSTPtr->qryOffsetSI;

   indexToCoord(
      alnSTPtr->refLenSI,
      alnSTPtr->endArySL[baseSL],
      &tmpRefUL,
      &tmpQryUL
   );
   *refEndSIPtr = tmpRefUL + alnSTPtr->refOffsetSI;
   *qryEndSIPtr = tmpQryUL + alnSTPtr->qryOffsetSI;

   return alnSTPtr->scoreArySL[baseSL];

   noCoord_fun05_sec04:;
      return 0;
   boundsErr_fun05_sec04:;
      return -1;
} /*getCoord_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun06: swap_memwaterScan
|   - swaps position in reference array
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to swap positions
|   - firstSI:
|     o first index to swap
|   - secSI:
|     o second index to swap
| Output:
|   - Modifies:
|     o scoreArySL, startArySL, endArySL in alnSTPtr to
|       have firstSI and secSI index's swapped
\-------------------------------------------------------*/
void
swap_memwaterScan(
   struct aln_memwaterScan *alnSTPtr,
   signed int firstSI,
   signed int secSI
){
   alnSTPtr->scoreArySL[firstSI] ^=
      alnSTPtr->scoreArySL[secSI];
   alnSTPtr->scoreArySL[secSI] ^=
      alnSTPtr->scoreArySL[firstSI];
   alnSTPtr->scoreArySL[firstSI] ^=
      alnSTPtr->scoreArySL[secSI];

   alnSTPtr->startArySL[firstSI] ^=
      alnSTPtr->startArySL[secSI];
   alnSTPtr->startArySL[secSI] ^=
      alnSTPtr->startArySL[firstSI];
   alnSTPtr->startArySL[firstSI] ^=
      alnSTPtr->startArySL[secSI];

   alnSTPtr->endArySL[firstSI] ^=
      alnSTPtr->endArySL[secSI];
   alnSTPtr->endArySL[secSI] ^=
      alnSTPtr->endArySL[firstSI];
   alnSTPtr->endArySL[firstSI] ^=
      alnSTPtr->endArySL[secSI];
} /*swap_memwaterScan*/

/*-------------------------------------------------------\
| Fun07: refCoordSort_aln_memwaterScan
|   - sorts best bases alignments by reference start,
|     query start, reference end (closest first), then
|     and query end (closest first)
|   - all negatives (non-alignment found) are pushed to
|     the end
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pionter to sort
| Output:
|   - Modifies:
|     o qryStartArySL in alnSTPtr to be sorted by
|       reference coodrinate
|     o qryEndArySL and qryScoreArySL in alnSTPtr to be
|       kept in sync with qryStartArySL
\-------------------------------------------------------*/
void
refCoordSort_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC:
   '   - sorts best bases alignments by reference start,
   '     query start, reference end (closest first), then
   '     query end (closest first), and finally score.
   '     all negatives (non-alignment found) are pushed to
   '     the end
   '   o fun07 sec01:
   '     - variable declerations
   '   o fun07 sec02:
   '     - find the number of rounds to sort for
   '   o fun07 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of elements to sort*/
   signed long numElmSL = alnSTPtr->outLenSL - 1;

   /*Number of sorting rounds*/
   signed int subSI = 0;
   signed int nextSI = 0;
   signed int lastSI = 0;
   signed int onSI = 0;

   unsigned long onRefStartUL;
   unsigned long onRefEndUL;
   unsigned long cmpRefStartUL;
   unsigned long cmpRefEndUL;

   unsigned long onQryStartUL;
   unsigned long onQryEndUL;
   unsigned long cmpQryStartUL;
   unsigned long cmpQryEndUL;

   /*Variables to incurment loops*/
   signed int siIndex = 0;
   signed int siElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numElmSL <= 0)
      return;

   /*Recursion formsla: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSI = 1; /*Initialzie first array*/

   while(subSI < numElmSL - 1)
      subSI = (3 * subSI) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - sort arrays
   ^   o fun07 sec03 sub01:
   ^     - outer loops for shell short
   ^   o fun07 sec03 sub02:
   ^     - compare and check if need to do first swap
   ^   o fun07 sec03 sub03:
   ^     - do inner (insertion) swaps
   ^   o fun07 sec03 sub04:
   ^     - move to the next element in the group
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec03 Sub01:
   *   - outer loops for shell short
   \*****************************************************/

   while(subSI > 0)
   { /*Loop: all rounds*/
      for(siIndex = 0; siIndex <= subSI; ++siIndex)
      { /*Loop: though sub array*/
         siElm = siIndex;

         for(
            siElm = siIndex;
            siElm + subSI <= numElmSL;
            siElm += subSI
         ){ /*Loop: swap elements in subarray*/
            nextSI = siElm + subSI;

            /********************************************\
            * Fun07 Sec03 Sub02:
            *   - compare & check if need to do first swap
            \********************************************/

            if(alnSTPtr->startArySL[siElm] < 0)
               goto swapOuter_fun07_sec03_sub02;
            else if(alnSTPtr->startArySL[nextSI] < 0)
               goto nextElm_fun07_sec03_sub04;

            indexToCoord(
               alnSTPtr->refLenSI,
               alnSTPtr->startArySL[siElm],
               &onRefStartUL,  
               &onQryStartUL
            );

            indexToCoord(
               alnSTPtr->refLenSI,
               alnSTPtr->startArySL[nextSI],
               &cmpRefStartUL,  
               &cmpQryStartUL
            );

            indexToCoord(
               alnSTPtr->refLenSI,
               alnSTPtr->endArySL[siElm],
               &onRefEndUL,  
               &onQryEndUL
            );

            indexToCoord(
               alnSTPtr->refLenSI,
               alnSTPtr->endArySL[nextSI],
               &cmpRefEndUL,  
               &cmpQryEndUL
            );

            if(onRefStartUL > cmpRefStartUL)
               goto swapOuter_fun07_sec03_sub02;
            else if(onRefStartUL < cmpRefStartUL)
               goto nextElm_fun07_sec03_sub04;

            else if(onQryStartUL > cmpQryStartUL)
               goto swapOuter_fun07_sec03_sub02;
            else if(onQryStartUL < cmpQryStartUL)
               goto nextElm_fun07_sec03_sub04;

            else if(onRefEndUL > cmpRefEndUL)
               goto swapOuter_fun07_sec03_sub02;
            else if(onRefEndUL < cmpRefEndUL)
               goto nextElm_fun07_sec03_sub04;

            else if(onQryEndUL > cmpQryEndUL)
            { /*Else If: need to swap (qry is later)*/
               swapOuter_fun07_sec03_sub02:;
               swap_memwaterScan(alnSTPtr, siElm, nextSI);

               /*****************************************\
               * Fun07 Sec03 Sub03:
               *   - do inner (insertion) swaps
               \*****************************************/

               lastSI = siElm;
               onSI = siElm;

               while(lastSI >= subSI)
               { /*loop: move swapped element back*/
                  lastSI -= subSI;

                  if(alnSTPtr->startArySL[onSI] < 0)
                     break;
                  else if(alnSTPtr->startArySL[lastSI] <0)
                     goto swapInner_fun07_sec03_sub04;

                  indexToCoord(
                     alnSTPtr->refLenSI,
                     alnSTPtr->startArySL[onSI],
                     &onRefStartUL,  
                     &onQryStartUL
                  );

                  indexToCoord(
                     alnSTPtr->refLenSI,
                     alnSTPtr->startArySL[lastSI],
                     &cmpRefStartUL,  
                     &cmpQryStartUL
                  );

                  indexToCoord(
                     alnSTPtr->refLenSI,
                     alnSTPtr->endArySL[onSI],
                     &onRefEndUL,  
                     &onQryEndUL
                  );

                  indexToCoord(
                     alnSTPtr->refLenSI,
                     alnSTPtr->endArySL[lastSI],
                     &cmpRefEndUL,  
                     &cmpQryEndUL
                  );

                  if(onRefStartUL > cmpRefStartUL)
                     break;
                  else if(onRefStartUL < cmpRefStartUL)
                     ;

                  else if(onQryStartUL > cmpQryStartUL)
                     break;
                  else if(onQryStartUL < cmpQryStartUL)
                     ;

                  else if(onRefEndUL > cmpRefEndUL)
                     break;
                  else if(onRefEndUL < cmpRefEndUL)
                     ;

                  else if(onQryEndUL > cmpQryEndUL)
                     break;

                  swapInner_fun07_sec03_sub04:;
                  swap_memwaterScan(alnSTPtr,onSI,lastSI);
                  onSI = lastSI;
               } /*Loop: move swapped element back*/
            }  /*Else If: need to swap (qry is later)*/

            /********************************************\
            * Fun07 Sec03 Sub04:
            *   - move to the next element in the group
            \********************************************/

            nextElm_fun07_sec03_sub04:;
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSI = (subSI - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*refCoordSort_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun08: filter_memwaterScan
|   - removes low scoring alignments and alignments that
|     are nested alignments
| Input:
|   - alnSTPtr:
|     o pointer to aln_memwaterScan struct with alignments
|       to filter
|   - minScoreSL:
|     o minimum score for an alignment
|   - minScoreTwoSL:
|     o second minimum score for an alignment
|       * use when you are putting a percent minimum score
|         for minScoreSL
|       * allows for a hard cutoff and a percent cuttoff
| Output:
|   - Modifies:
|     o startArySL, endArySL, and scoreArySL in alnSTPtr
|       to only have high scoring, non-nested alignments
|     o outLenSL to have the new number of alignments
|   - Returns:
|     o number of kept alignments
\-------------------------------------------------------*/
signed long
filter_memwaterScan(
   struct aln_memwaterScan *alnSTPtr,
   signed long minScoreSL,
   signed long minScoreTwoSL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun08 TOC:
   '   - removes low scoring alignments and alignments
   '     that are nested alignments
   '   o fun08 sec01:
   '     - variable declarations
   '   o fun08 sec02:
   '     - remove low scoring alignments & sort position
   '   o fun08 sec03:
   '     - mark nested alignments
   '   o fun08 sec04:
   '     - remove nested alignments
   '   o fun08 sec05:
   '     - resort by coordinate and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long refStartUL = 0;
   unsigned long refEndUL = 0;
   unsigned long refLastStartUL = 0;
   unsigned long refLastEndUL = 0;

   unsigned long qryStartUL = 0;
   unsigned long qryEndUL = 0;
   unsigned long qryLastStartUL = 0;
   unsigned long qryLastEndUL = 0;

   signed long keptPosSL = -1;
   signed long posSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec02:
   ^   - remove low scoring alignments & sort by position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSL = 0;
   while(posSL < alnSTPtr->outLenSL)
   { /*Loop: filter out low scoring alignments*/
      if(
            alnSTPtr->scoreArySL[posSL] < minScoreSL
         || alnSTPtr->scoreArySL[posSL] < minScoreTwoSL
      ){ /*If: under minimum score*/
         swap_memwaterScan(
            alnSTPtr,
            posSL,
            alnSTPtr->outLenSL - 1
         );

         --alnSTPtr->outLenSL;
      } /*If: under minimum score*/

      else
         ++posSL;
   } /*Loop: filter out low scoring alignments*/

   /*need to sort so nested alignments are next to each
   `  other (in most cases)
   */
   refCoordSort_aln_memwaterScan(alnSTPtr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec03:
   ^   - mark nested alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   keptPosSL = -1;
   posSL = 0;
   while(posSL < alnSTPtr->outLenSL)
   { /*Loop: filter out alignments*/
      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->startArySL[posSL],
         &refStartUL,  
         &qryStartUL
      );

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->endArySL[posSL],
         &refEndUL,  
         &qryEndUL
      );

      if(refStartUL != refLastStartUL)
         keptPosSL = posSL;

      else if(qryStartUL != qryLastStartUL)
         keptPosSL = posSL;

      else if(
            refEndUL - 1 != refLastEndUL
         && refEndUL != refLastEndUL
      ) keptPosSL = posSL;

      else if(
            qryEndUL - 1 != qryLastEndUL
         && qryEndUL != qryLastEndUL
      ) keptPosSL = posSL;

      else
      { /*Else: this is an overlapping alignment*/
         if(keptPosSL < 0)
            keptPosSL = posSL;

         else if(
               alnSTPtr->scoreArySL[keptPosSL]
            >= alnSTPtr->scoreArySL[posSL]
         ) alnSTPtr->startArySL[posSL] *= -1;

         else
         { /*Else: keeping new position*/
            alnSTPtr->startArySL[keptPosSL] *= -1;
            keptPosSL = posSL;
         } /*Else: keeping new position*/

         goto nextAln_fun08_sec0x_sub0x;
      } /*Else: this is an overlapping alignment*/

      nextAln_fun08_sec0x_sub0x:;
         refLastStartUL = refStartUL; 
         refLastEndUL = refEndUL; 

         qryLastStartUL = qryStartUL; 
         qryLastEndUL = qryEndUL; 
         ++posSL;
   } /*Loop: filter out alignments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec04:
   ^   - remove nested alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSL = 0;
   while(posSL < alnSTPtr->outLenSL)
   { /*Loop: remove flagged alignments*/
      if(alnSTPtr->startArySL[posSL] < 0)
      { /*If: flagged this for removal*/
         alnSTPtr->startArySL[posSL] *= -1;
         swap_memwaterScan(
            alnSTPtr,
            posSL,
            alnSTPtr->outLenSL - 1
         );

         --alnSTPtr->outLenSL;
      } /*If: flagged this for removal*/

      else
         ++posSL;
   } /*Loop: remove flagged alignments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec05:
   ^   - resort by coordinate and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refCoordSort_aln_memwaterScan(alnSTPtr);
   return alnSTPtr->outLenSL;
} /*filter_memwaterScan*/

/*-------------------------------------------------------\
| Fun09: merge_aln_memwaterScan
|   - merges one aln_memwaterScan struct into another
|     aln_memwaterScan struct
| Input:
|   - mergeSTPtr:
|     o aln_memwaterScan struct to merge into the final
|       structure
|   - outSTPtr:
|     o aln_memwaterScan struct to merge mergeSTPtr into
|   - minScoreSL:
|     o minumum score to keep score in mergeSTPtr
| Output:
|   - Modifies:
|     o variables in outSTPtr to include variables in
|       mergeSTPtr and outSTPtr
|   - Returns:
|     o 0 for success
|     o -1 if reference lengths differ (no way to merge)
|     o -2 for memory errors
\-------------------------------------------------------*/
signed char
merge_aln_memwaterScan(
   struct aln_memwaterScan *mergeSTPtr,
   struct aln_memwaterScan *outSTPtr,
   signed long minScoreSL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun09 TOC:
   '   - merges one aln_memwaterScan struct into another
   '     aln_memwaterScan struct
   '   o fun09 sec01:
   '     - variable declarations
   '   o fun09 sec02:
   '     - allocate memory for arrays
   '   o fun09 sec03:
   '     - merge alignments
   '   o fun09 sec04:
   '     - return the results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long *tmpSLPtr = 0;
   signed long newOutLenSL =
      mergeSTPtr->outLenSL + outSTPtr->outLenSL;
   signed long slAln = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec02:
   ^   - allocate memory for arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! outSTPtr->refLenSI)
      outSTPtr->refLenSI = mergeSTPtr->refLenSI;
   else if(outSTPtr->refLenSI != mergeSTPtr->refLenSI)
      goto differentRefsErr_fun09_sec0x;

   if(newOutLenSL >= outSTPtr->outSizeSL)
   { /*If: need to resize the output struct*/
      if(! outSTPtr->startArySL)
         tmpSLPtr =
            malloc(
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      else
         tmpSLPtr =
            realloc(
               outSTPtr->startArySL,
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      if(! tmpSLPtr)
         goto memErr_fun09_sec0x;
      outSTPtr->startArySL = tmpSLPtr;
      tmpSLPtr = 0;

      if(! outSTPtr->endArySL)
         tmpSLPtr =
            malloc(
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      else
         tmpSLPtr =
            realloc(
               outSTPtr->endArySL,
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      if(! tmpSLPtr)
         goto memErr_fun09_sec0x;
      outSTPtr->endArySL = tmpSLPtr;
      tmpSLPtr = 0;

      if(! outSTPtr->scoreArySL)
         tmpSLPtr =
            malloc(
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      else
         tmpSLPtr =
            realloc(
               outSTPtr->scoreArySL,
               (newOutLenSL + 16) * sizeof(unsigned long)
            );
      if(! tmpSLPtr)
         goto memErr_fun09_sec0x;
      outSTPtr->scoreArySL = tmpSLPtr;
      tmpSLPtr = 0;

      outSTPtr->outSizeSL = newOutLenSL + 16;
   } /*If: need to resize the output struct*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec03:
   ^   - merge alignments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(outSTPtr->qryLenSI < mergeSTPtr->qryLenSI)
      outSTPtr->qryLenSI = mergeSTPtr->qryLenSI;
      /*make sure always have the max query lengths*/

   for(slAln = 0; slAln < mergeSTPtr->outLenSL; ++slAln)
   { /*Loop: merge alignments*/
      if(mergeSTPtr->scoreArySL[slAln] < minScoreSL)
         continue;

      outSTPtr->startArySL[outSTPtr->outLenSL] =
         mergeSTPtr->startArySL[slAln];
      outSTPtr->endArySL[outSTPtr->outLenSL] =
         mergeSTPtr->endArySL[slAln];
      outSTPtr->scoreArySL[outSTPtr->outLenSL] =
         mergeSTPtr->scoreArySL[slAln];
      ++outSTPtr->outLenSL;
   } /*Loop: merge alignments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec04:
   ^   - return the results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   newOutLenSL = 0;
   goto ret_fun09_sec0x;

   differentRefsErr_fun09_sec0x:;
      newOutLenSL = -1;
      goto ret_fun09_sec0x;
  
   memErr_fun09_sec0x:;
      newOutLenSL = -1;
      goto ret_fun09_sec0x;

   ret_fun09_sec0x:;
      if(tmpSLPtr)
         free(tmpSLPtr);
      tmpSLPtr = 0;

      return newOutLenSL;
} /*merge_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun10: memwaterScan
|   - performs a memory efficent Smith Waterman scan
|     (keep best alignment for each query/reference base)
|     alignment on a pair of sequences
| Input;
|   - qrySeqSTPtr:
|     o pointer to an seqST with the query sequence and
|       index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - refSeqSTPtr:
|     o pointer to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - alnSTPtr:
|     o pointer to aln_memwaterScan structure to hold the
|       results of the alignment
|   - settings:
|     o pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o variables in alnSTPtr to have the new alignment
|  - Returns:
|    o score of the best alignment found
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
memwaterScan(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct aln_memwaterScan *alnSTPtr,/*gets alignment*/
   struct alnSet *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - performs a memory efficent Smith Waterman scan
   '     (keep best alignment for each query/reference
   '     base) alignment on a pair of sequences
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun10 sec01:
   '    - Variable declerations
   '  o fun10 sec02:
   '    - Allocate memory for alignment
   '  o fun10 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun10 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01: Variable declerations
   ^  o fun10 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun10 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun10 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun10 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   signed long scoreSL = 0;     /*score to return*/
   signed long bestScoreSL = 0;

   /*Get start & end of query and reference sequences*/
   signed char *refSeqStr = 0;
   signed char *qrySeqStr = 0;

   signed long refLenSL =
      refSTPtr->endAlnSL - refSTPtr->offsetSL + 1;

   signed long qryLenSL =
      qrySTPtr->endAlnSL - qrySTPtr->offsetSL + 1;

   /*Iterators for loops*/
   signed long slRef = 0;
   signed long slQry = 0;

   /*****************************************************\
   * Fun10 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   signed long snpScoreSL = 0;
   signed long insScoreSL = 0;
   signed long delScoreSL = 0;   /*score for deletion*/
   signed long nextSnpScoreSL = 0;/*next match/snp score*/

   /*****************************************************\
   * Fun10 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Keeping track of alignment starting positions*/
   signed long indexSL = 0;    /*index at in matrix*/
   signed long snpIndexSL = 0; /*last snp index*/
   signed long tmpIndexSL = 0; /*for getting snp index*/

   /*first query base in scoring array*/
   signed long startSL = refLenSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun10 sec02 sub01:
   ^    - get lengths and offsets
   ^  o fun10 sec02 sub02:
   ^    - output alignment coodinates memory allocation
   ^  o fun10 sec02 sub03:
   ^    - alignment rows memory allocate
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec02 Sub01:
   *   - get lengths and offsets
   \****************************************************/

   alnSTPtr->refLenSI = refLenSL;
   alnSTPtr->qryLenSI = qryLenSL;
   alnSTPtr->refOffsetSI = refSTPtr->offsetSL;
   alnSTPtr->qryOffsetSI = qrySTPtr->offsetSL;

   alnSTPtr->outLenSL = refLenSL + qryLenSL;

   /*****************************************************\
   * Fun10 Sec02 Sub02:
   *   - output alignment coodinates memory allocation
   \****************************************************/

   if(alnSTPtr->outSizeSL < alnSTPtr->outLenSL)
   { /*If: need to resize the alignment coordinates*/
      if(alnSTPtr->startArySL)
         free(alnSTPtr->startArySL);
      alnSTPtr->startArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->startArySL)
         goto memErr_fun10_sec05_sub03;

      if(alnSTPtr->endArySL)
         free(alnSTPtr->endArySL);
      alnSTPtr->endArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->endArySL)
         goto memErr_fun10_sec05_sub03;

      if(alnSTPtr->scoreArySL)
         free(alnSTPtr->scoreArySL);
      alnSTPtr->scoreArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->scoreArySL)
         goto memErr_fun10_sec05_sub03;

      alnSTPtr->outSizeSL = alnSTPtr->outLenSL;
   } /*If: need to resize the alignment coordinates*/

   /*****************************************************\
   * Fun10 Sec02 Sub03:
   *   - alignment rows memory allocate
   \****************************************************/

   if(alnSTPtr->rowSizeSI < refLenSL)
   { /*If: need more memory for alignment rows*/
      if(alnSTPtr->indexRowSL)
         free(alnSTPtr->indexRowSL);
      alnSTPtr->indexRowSL =
         malloc((refLenSL + 1) * sizeof(signed long));
      if(! alnSTPtr->indexRowSL)
         goto memErr_fun10_sec05_sub03;

      if(alnSTPtr->scoreRowSL)
         free(alnSTPtr->scoreRowSL);
      alnSTPtr->scoreRowSL =
         malloc((refLenSL + 1) * sizeof(signed long));
      if(! alnSTPtr->scoreRowSL)
         goto memErr_fun10_sec05_sub03;

      if(alnSTPtr->dirRowSC)
         free(alnSTPtr->dirRowSC);
      alnSTPtr->dirRowSC =
         malloc((refLenSL + 1) * sizeof(signed long));
      if(! alnSTPtr->dirRowSC)
         goto memErr_fun10_sec05_sub03;

      alnSTPtr->rowSizeSI = refLenSL;
   } /*If: need more memory for alignment rows*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^  - initialize all values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(indexSL=0; indexSL < alnSTPtr->outLenSL; ++indexSL)
   { /*loop; till have initalized the query scores*/
      alnSTPtr->startArySL[indexSL] = -1;
      alnSTPtr->endArySL[indexSL] = -1;
      alnSTPtr->scoreArySL[indexSL] = -1;
   } /*loop; till have initalized the query scores*/

   for(indexSL = 0; indexSL <= refLenSL; ++indexSL)
   { /*loop; till have initalized the first row*/
      alnSTPtr->dirRowSC[indexSL] = def_mvStop_alnDefs;
      alnSTPtr->indexRowSL[indexSL] = indexSL + 1;
      alnSTPtr->scoreRowSL[indexSL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun10 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun10 sec04 sub02:
   ^    - get snp and ins scores + start loop
   ^  o fun10 sec04 sub03:
   ^    - find high score
   ^  o fun10 sec04 sub04:
   ^    - check if keep score (score > 0)
   ^  o fun10 sec04 sub05:
   ^    - find next deletion score and move to next index
   ^  o fun10 sec04 sub07:
   ^    - prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   alnSTPtr->indexRowSL[0] = indexSL;

   /*These are always negative*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   snpIndexSL = alnSTPtr->indexRowSL[0];

   /*incurment to frist base*/
   ++indexSL;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetSL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetSL - 1;
      /*offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun10 Sec04 Sub02:
   *  - get snp and ins scores + start loop
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(slQry = 0; slQry < qryLenSL; ++slQry)
   { /*loop; compare query base against all ref bases*/
      for(slRef = 1; slRef <= refLenSL; ++slRef)
      { /*loop; compare one query to one reference base*/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[slQry],
               refSeqStr[slRef],
               settings
            ); /*find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = alnSTPtr->scoreRowSL[slRef];

         /*insertion score*/
         #ifdef NOEXTEND
            insScoreSL =
                 alnSTPtr->scoreRowSL[slRef]
               + settings->gapSS;
         #else
            insScoreSL = alnSTPtr->scoreRowSL[slRef];
            insScoreSL +=
               settings->insArySS[
                  alnSTPtr->dirRowSC[slRef]
               ];
         #endif

         /***********************************************\
         * Fun10 Sec04 Sub03:
         *   - find high score
         \***********************************************/

         alnSTPtr->scoreRowSL[slRef] =
            max_genMath(insScoreSL, snpScoreSL);
            /*find if ins/snp is best (5 Op)*/

         tmpIndexSL = alnSTPtr->indexRowSL[slRef];

         alnSTPtr->indexRowSL[slRef] =
            ifmax_genMath(
               insScoreSL,
               snpScoreSL,
               alnSTPtr->indexRowSL[slRef], /*ins index*/
               snpIndexSL                   /*snp index*/
            ); /*get index of high score*/

         snpIndexSL = tmpIndexSL;

         /*find direction (5 Op)*/
         alnSTPtr->dirRowSC[slRef] =
            alnSTPtr->scoreRowSL[slRef] > delScoreSL;
         alnSTPtr->dirRowSC[slRef] +=
            (
                 (snpScoreSL <= insScoreSL)
               & alnSTPtr->dirRowSC[slRef]
            );
         ++alnSTPtr->dirRowSC[slRef];

         /*Logic:
         `   - noDel: maxSC > delSc:
         `     o 1 if deletion not max score
         `     o 0 if deletion is max score
         `   - type: noDel + ((snpSc < insSc) & noDel):
         `     o 1 + (1 & 1) = 2 if insertion is maximum
         `     o 1 + (0 & 1) = 1 if snp is maximum
         `     o 0 + (0 & 0) = 0 if del is max; snp > ins
         `     o 0 + (1 & 0) = 0 if del is max, ins >= snp
         `   - dir: type + 1
         `     o adds 1 to change from stop to direction
         */

         /*finish finding max's*/
         alnSTPtr->indexRowSL[slRef] =
            ifmax_genMath(
               delScoreSL,
               alnSTPtr->scoreRowSL[slRef],
               alnSTPtr->indexRowSL[slRef-1],/*del index*/
               alnSTPtr->indexRowSL[slRef]/*current best*/
            ); /*get index of high score*/

         alnSTPtr->scoreRowSL[slRef] =
            max_genMath(
               delScoreSL,
               alnSTPtr->scoreRowSL[slRef]
         ); /*find if del is best (5 Op)*/
            
         /***********************************************\
         * Fun10 Sec04 Sub04:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(alnSTPtr->scoreRowSL[slRef] <= 0)
         {
            alnSTPtr->dirRowSC[slRef] = 0;
            alnSTPtr->scoreRowSL[slRef] = 0;
            alnSTPtr->indexRowSL[slRef] = indexSL + 1;
            /*always one off for new index's*/
            /* branchless method is slower here*/
         }

         else
         { /*Else: check if have new high score*/
            scoreSL = alnSTPtr->scoreRowSL[slRef];

            if(scoreSL > alnSTPtr->scoreArySL[slRef-1])
            { /*If: new best outerence position score*/
               alnSTPtr->startArySL[slRef - 1] =
                   alnSTPtr->indexRowSL[slRef];
               alnSTPtr->endArySL[slRef - 1] = indexSL;
               alnSTPtr->scoreArySL[slRef-1] = scoreSL;

               if(scoreSL > bestScoreSL)
                  bestScoreSL = scoreSL;
            } /*If: new best outerence position score*/

            if(
                 scoreSL
               > alnSTPtr->scoreArySL[slQry+startSL]
            ){ /*If: new best query position score*/
               alnSTPtr->startArySL[slQry + startSL] =
                  alnSTPtr->indexRowSL[slRef];
               alnSTPtr->endArySL[slQry + startSL] =
                  indexSL;
               alnSTPtr->scoreArySL[slQry + startSL] =
                  scoreSL;

               if(scoreSL > bestScoreSL)
                  bestScoreSL = scoreSL;
            } /*If: new best query position score*/
         } /*Else: check if have new high score*/

         /***********************************************\
         * Fun10 Sec04 Sub05:
         *   - find next deletion score and move to next
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
                 alnSTPtr->scoreRowSL[slRef]
               + settings->gapSS;
         #else
            delScoreSL = alnSTPtr->scoreRowSL[slRef];
            delScoreSL +=
               settings->delArySS[
                  alnSTPtr->dirRowSC[slRef]
               ];
         #endif

         ++indexSL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun10 Sec04 Sub07:
     *   - prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     alnSTPtr->indexRowSL[0] = indexSL; /*next index*/
     snpIndexSL = alnSTPtr->indexRowSL[0];

     ++indexSL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec05:
   ^  - set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto cleanUp_fun10_sec05;

   memErr_fun10_sec05_sub03:;
      bestScoreSL = -1;
      goto cleanUp_fun10_sec05;

   cleanUp_fun10_sec05:;
      return bestScoreSL;
} /*memwaterScan*/

/*-------------------------------------------------------\
| Fun11: simple_memwaterScan
|   - performs a memory efficent Smith Waterman scan
|     (keep best alignment for each query/reference base)
|     alignment on a pair of sequences
|   - simple means no settings or seqST structure used and
|     no match matrix
| Input;
|   - qrySeqStr:
|     o c-string with query sequence
|   - qryLenSI:
|     o length of qrySeqStr (index 1)
|   - refSeqStr:
|     o c-string with reference sequence
|   - refLenSI:
|     o length of refSeqStr (index 1)
|   - alnSTPtr:
|     o pointer to aln_memwaterScan structure to hold the
|       results of the alignment
| Output:
|  - Modifies:
|    o variables in alnSTPtr to have the new alignment
|  - Returns:
|    o score of the best alignment found
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
simple_memwaterScan(
   signed char *qrySeqStr,   /*query sequence*/
   signed int qryLenSI,      /*length of query sequence*/
   signed char *refSeqStr,   /*reference sequence*/
   signed int refLenSI,      /*reference sequence length*/
   struct aln_memwaterScan *alnSTPtr /*gets alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC: simple_memwaterScan
   '   - performs a memory efficent Smith Waterman scan
   '     (keep best alignment for each query/reference
   '     base) alignment on a pair of sequences
   '   - simple means no settings or seqST structure used
   '     and no match matrix
   '  o fun11 sec01:
   '    - Variable declerations
   '  o fun11 sec02:
   '    - Allocate memory for alignment
   '  o fun11 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun11 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01: Variable declerations
   ^  o fun11 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun11 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun11 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun11 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   signed long scoreSL = 0;     /*holds score*/
   signed long bestScoreSL = 0; /*score to return*/

   /*Iterators for loops*/
   signed int siRef = 0;
   signed int siQry = 0;

   /*****************************************************\
   * Fun11 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   signed long snpScoreSL = 0;
   signed long insScoreSL = 0;
   signed long delScoreSL = 0;   /*score for deletion*/
   signed long nextSnpScoreSL = 0;/*next match/snp score*/

   /*****************************************************\
   * Fun11 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Keeping track of alignment starting positions*/
   signed long indexSL = 0;    /*index at in matrix*/
   signed long snpIndexSL = 0; /*last snp index*/
   signed long tmpIndexSL = 0; /*for getting snp index*/

   /*first query base in scoring array*/
   signed long startSL = refLenSI;

   signed int insArySI[4];
   signed int delArySI[4];

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun11 sec02 sub01:
   ^    - get lengths and offsets
   ^  o fun11 sec02 sub02:
   ^    - output alignment coodinates memory allocation
   ^  o fun11 sec02 sub03:
   ^    - alignment rows memory allocate
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec02 Sub01:
   *   - get lengths and offsets
   \****************************************************/

   alnSTPtr->refLenSI = refLenSI;
   alnSTPtr->qryLenSI = qryLenSI;
   alnSTPtr->refOffsetSI = 0;
   alnSTPtr->qryOffsetSI = 0;

   alnSTPtr->outLenSL = refLenSI + qryLenSI;

   insArySI[def_mvStop_alnDefs] = def_gapOpen_alnDefs;
   insArySI[def_mvDel_alnDefs] = def_gapOpen_alnDefs;
   insArySI[def_mvSnp_alnDefs] = def_gapOpen_alnDefs;
   insArySI[def_mvIns_alnDefs] = def_gapExtend_alnDefs;

   delArySI[def_mvStop_alnDefs] = def_gapOpen_alnDefs;
   delArySI[def_mvDel_alnDefs] = def_gapExtend_alnDefs;
   delArySI[def_mvSnp_alnDefs] = def_gapOpen_alnDefs;
   delArySI[def_mvIns_alnDefs] = def_gapOpen_alnDefs;

   /*****************************************************\
   * Fun11 Sec02 Sub02:
   *   - output alignment coodinates memory allocation
   \****************************************************/

   if(alnSTPtr->outSizeSL < alnSTPtr->outLenSL)
   { /*If: need to resize the alignment coordinates*/
      if(alnSTPtr->startArySL)
         free(alnSTPtr->startArySL);
      alnSTPtr->startArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->startArySL)
         goto memErr_fun11_sec05_sub03;

      if(alnSTPtr->endArySL)
         free(alnSTPtr->endArySL);
      alnSTPtr->endArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->endArySL)
         goto memErr_fun11_sec05_sub03;

      if(alnSTPtr->scoreArySL)
         free(alnSTPtr->scoreArySL);
      alnSTPtr->scoreArySL =
         malloc(alnSTPtr->outLenSL * sizeof(signed long));
      if(! alnSTPtr->scoreArySL)
         goto memErr_fun11_sec05_sub03;

      alnSTPtr->outSizeSL = alnSTPtr->outLenSL;
   } /*If: need to resize the alignment coordinates*/

   /*****************************************************\
   * Fun11 Sec02 Sub03:
   *   - alignment rows memory allocate
   \****************************************************/

   if(alnSTPtr->rowSizeSI < refLenSI)
   { /*If: need more memory for alignment rows*/
      if(alnSTPtr->indexRowSL)
         free(alnSTPtr->indexRowSL);
      alnSTPtr->indexRowSL =
         malloc((refLenSI + 1) * sizeof(signed long));
      if(! alnSTPtr->indexRowSL)
         goto memErr_fun11_sec05_sub03;

      if(alnSTPtr->scoreRowSL)
         free(alnSTPtr->scoreRowSL);
      alnSTPtr->scoreRowSL =
         malloc((refLenSI + 1) * sizeof(signed long));
      if(! alnSTPtr->scoreRowSL)
         goto memErr_fun11_sec05_sub03;

      if(alnSTPtr->dirRowSC)
         free(alnSTPtr->dirRowSC);
      alnSTPtr->dirRowSC =
         malloc((refLenSI + 1) * sizeof(signed long));
      if(! alnSTPtr->dirRowSC)
         goto memErr_fun11_sec05_sub03;

      alnSTPtr->rowSizeSI = refLenSI;
   } /*If: need more memory for alignment rows*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^  - initialize all values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(indexSL=0; indexSL < alnSTPtr->outLenSL; ++indexSL)
   { /*loop; till have initalized the query scores*/
      alnSTPtr->startArySL[indexSL] = -1;
      alnSTPtr->endArySL[indexSL] = -1;
      alnSTPtr->scoreArySL[indexSL] = -1;
   } /*loop; till have initalized the query scores*/

   for(indexSL = 0; indexSL <= refLenSI; ++indexSL)
   { /*loop; till have initalized the first row*/
      alnSTPtr->dirRowSC[indexSL] = def_mvStop_alnDefs;
      alnSTPtr->indexRowSL[indexSL] = indexSL + 1;
      alnSTPtr->scoreRowSL[indexSL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun11 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun11 sec04 sub02:
   ^    - get snp and ins scores + start loop
   ^  o fun11 sec04 sub03:
   ^    - find high score
   ^  o fun11 sec04 sub04:
   ^    - check if keep score (score > 0)
   ^  o fun11 sec04 sub05:
   ^    - find next deletion score and move to next index
   ^  o fun11 sec04 sub07:
   ^    - prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   alnSTPtr->indexRowSL[0] = indexSL;

   /*These are always negative*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   snpIndexSL = alnSTPtr->indexRowSL[0];

   /*incurment to frist base*/
   ++indexSL;
   --refSeqStr;
      /*offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun11 Sec04 Sub02:
   *  - get snp and ins scores + start loop
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(siQry = 0; siQry < qryLenSI; ++siQry)
   { /*loop; compare query base against all ref bases*/
      for(siRef = 1; siRef <= refLenSI; ++siRef)
      { /*loop; compare one query to one reference base*/

         snpScoreSL =
              def_matchScore_alnDefs
            - (
               (
                    def_matchScore_alnDefs
                  - def_snpScore_alnDefs
               ) & -(
                          (qrySeqStr[siQry] & ~32)
                       != (refSeqStr[siRef] & ~32)
                    )
            );

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = alnSTPtr->scoreRowSL[siRef];

         /*insertion score*/
         #ifdef NOEXTEND
            insScoreSL =
              scoreHeapArySL[siRef] + def_gapOpen_alnDefs;
         #else
            insScoreSL =
                 alnSTPtr->scoreRowSL[siRef]
               + insArySI[alnSTPtr->dirRowSC[siRef]];
         #endif

         /***********************************************\
         * Fun11 Sec04 Sub03:
         *   - find high score
         \***********************************************/

         alnSTPtr->scoreRowSL[siRef] =
            max_genMath(insScoreSL, snpScoreSL);
            /*find if ins/snp is best (5 Op)*/

         tmpIndexSL = alnSTPtr->indexRowSL[siRef];

         alnSTPtr->indexRowSL[siRef] =
            ifmax_genMath(
               insScoreSL,
               snpScoreSL,
               alnSTPtr->indexRowSL[siRef], /*ins index*/
               snpIndexSL                   /*snp index*/
            ); /*get index of high score*/

         snpIndexSL = tmpIndexSL;

         /*find direction (5 Op)*/
         alnSTPtr->dirRowSC[siRef] =
            alnSTPtr->scoreRowSL[siRef] > delScoreSL;
         alnSTPtr->dirRowSC[siRef] +=
            (
                 (snpScoreSL <= insScoreSL)
               & alnSTPtr->dirRowSC[siRef]
            );
         ++alnSTPtr->dirRowSC[siRef];

         /*Logic:
         `   - noDel: maxSC > delSc:
         `     o 1 if deletion not max score
         `     o 0 if deletion is max score
         `   - type: noDel + ((snpSc < insSc) & noDel):
         `     o 1 + (1 & 1) = 2 if insertion is maximum
         `     o 1 + (0 & 1) = 1 if snp is maximum
         `     o 0 + (0 & 0) = 0 if del is max; snp > ins
         `     o 0 + (1 & 0) = 0 if del is max, ins >= snp
         `   - dir: type + 1
         `     o adds 1 to change from stop to direction
         */

         /*finish finding max's*/
         alnSTPtr->indexRowSL[siRef] =
            ifmax_genMath(
               delScoreSL,
               alnSTPtr->scoreRowSL[siRef],
               alnSTPtr->indexRowSL[siRef-1],/*del index*/
               alnSTPtr->indexRowSL[siRef]/*current best*/
            ); /*get index of high score*/

         alnSTPtr->scoreRowSL[siRef] =
            max_genMath(
               delScoreSL,
               alnSTPtr->scoreRowSL[siRef]
         ); /*find if del is best (5 Op)*/
            
         /***********************************************\
         * Fun11 Sec04 Sub04:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(alnSTPtr->scoreRowSL[siRef] <= 0)
         {
            alnSTPtr->dirRowSC[siRef] = 0;
            alnSTPtr->scoreRowSL[siRef] = 0;
            alnSTPtr->indexRowSL[siRef] = indexSL + 1;
            /*always one off for new index's*/
            /* branchless method is slower here*/
         }

         else
         { /*Else: check if have new high score*/
            scoreSL = alnSTPtr->scoreRowSL[siRef];

            if(scoreSL > alnSTPtr->scoreArySL[siRef - 1])
            { /*If: new best outerence position score*/
               alnSTPtr->startArySL[siRef - 1] =
                   alnSTPtr->indexRowSL[siRef];
               alnSTPtr->endArySL[siRef - 1] = indexSL;
               alnSTPtr->scoreArySL[siRef-1] = scoreSL;

               if(scoreSL > bestScoreSL)
                  bestScoreSL = scoreSL;
            } /*If: new best outerence position score*/

            if(
                 scoreSL
               > alnSTPtr->scoreArySL[siQry + startSL]
            ){ /*If: new best query position score*/
               alnSTPtr->startArySL[siQry + startSL] =
                  alnSTPtr->indexRowSL[siRef];
               alnSTPtr->endArySL[siQry + startSL] =
                  indexSL;
               alnSTPtr->scoreArySL[siQry + startSL] =
                  scoreSL;

               if(scoreSL > bestScoreSL)
                  bestScoreSL = scoreSL;
            } /*If: new best query position score*/
         } /*Else: check if have new high score*/

         /***********************************************\
         * Fun11 Sec04 Sub05:
         *   - find next deletion score and move to next
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
              scoreHeapArySL[siRef] + def_gapOpen_alnDefs;
         #else
            delScoreSL =
                 alnSTPtr->scoreRowSL[siRef]
               + delArySI[alnSTPtr->dirRowSC[siRef]];
         #endif

         ++indexSL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun11 Sec04 Sub07:
     *   - prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     alnSTPtr->indexRowSL[0] = indexSL; /*next index*/
     snpIndexSL = alnSTPtr->indexRowSL[0];

     ++indexSL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^  - set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto cleanUp_fun11_sec05;

   memErr_fun11_sec05_sub03:;
      bestScoreSL = -1;
      goto cleanUp_fun11_sec05;

   cleanUp_fun11_sec05:;
      return bestScoreSL;
} /*simple_memwaterScan*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible, this code is under the
:   MIT license.
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
