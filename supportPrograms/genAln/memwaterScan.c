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
'   o fun06: qryByRefSort_aln_memwaterScan
'     - sorts the query coordinates in a aln_memwaterScan
'       structure by the reference coordinate
'     - the sorting is first smallest query state, then
'       highest score, then longest (ending base) length
'   o fun07 memwaterScan:
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference base)
'       alignment on a pair of sequences
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

   alnSTPtr->refStartArySL = 0;
   alnSTPtr->refEndArySL = 0;
   alnSTPtr->refScoreArySL = 0;

   alnSTPtr->qryStartArySL = 0;
   alnSTPtr->qryEndArySL = 0;
   alnSTPtr->qryScoreArySL = 0;
   alnSTPtr->qrySizeSI = 0;

   alnSTPtr->scoreArySL = 0;
   alnSTPtr->indexArySL = 0;
   alnSTPtr->dirArySC = 0;
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

   if(alnSTPtr->refStartArySL)
      free(alnSTPtr->refStartArySL);
   if(alnSTPtr->refEndArySL)
      free(alnSTPtr->refEndArySL);
   if(alnSTPtr->refScoreArySL)
      free(alnSTPtr->refScoreArySL);

   if(alnSTPtr->qryStartArySL)
      free(alnSTPtr->qryStartArySL);
   if(alnSTPtr->qryEndArySL)
      free(alnSTPtr->qryEndArySL);
   if(alnSTPtr->qryScoreArySL)
      free(alnSTPtr->qryScoreArySL);

   if(alnSTPtr->scoreArySL)
      free(alnSTPtr->scoreArySL);
   if(alnSTPtr->indexArySL)
      free(alnSTPtr->indexArySL);
   if(alnSTPtr->dirArySC)
      free(alnSTPtr->dirArySC);

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
|   - qryBl:
|     o 1 if index is for the query array
|     o 0 if index is for the reference array
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
|     o frees alnSTPtr; you must set the pointer to 0/null
|   - Returns:
|     o 0 for success
|     o 1 if no position at index
|     o -1 for out of bounds error
\-------------------------------------------------------*/
signed char
getCoord_aln_memwaterScan(
   signed int baseSI,
   signed char qryBl,         /*1: get from query array*/
   signed int *refStartSIPtr, /*gets reference start*/
   signed int *refEndSIPtr,   /*gets reference end*/
   signed int *qryStartSIPtr, /*gets query start*/
   signed int *qryEndSIPtr,   /*gets query end*/
   struct aln_memwaterScan *alnSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - gets the coordinates for input base index for an
   '     aln_memwaterScan structure
   '   o fun05 sec01:
   '     - variable declarations
   '   o fun05 sec02:
   '     - get query array coordaintes (if requested)
   '   o fun05 sec03:
   '     - get reference array coordaintes (if requested)
   '   o fun05 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long tmpRefUL = 0;
   unsigned long tmpQryUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - get query array coordaintes (if requested)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   if(qryBl)
   { /*If: using the query array*/
      if(baseSI > alnSTPtr->qryLenSI)
         goto boundsErr_fun05_sec04;
      else if(alnSTPtr->qryScoreArySL[baseSI] < 0)
         goto noCoord_fun05_sec04;

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->qryStartArySL[baseSI],
         &tmpRefUL,
         &tmpQryUL
      );
      *qryStartSIPtr = tmpQryUL;
      *refStartSIPtr = tmpRefUL;

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->qryEndArySL[baseSI],
         &tmpRefUL,
         &tmpQryUL
      );
      *qryEndSIPtr = tmpQryUL;
      *refEndSIPtr = tmpRefUL;
   } /*If: using the query array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - get reference array coordaintes (if requested)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else
   { /*Else: using the reference array*/
      if(baseSI > alnSTPtr->refLenSI)
         goto boundsErr_fun05_sec04;
      else if(alnSTPtr->refScoreArySL[baseSI] < 0)
         goto noCoord_fun05_sec04;

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->refStartArySL[baseSI],
         &tmpRefUL,
         &tmpQryUL
      );
      *qryStartSIPtr = tmpQryUL;
      *refStartSIPtr = tmpRefUL;

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->refEndArySL[baseSI],
         &tmpRefUL,
         &tmpQryUL
      );
      *qryEndSIPtr = tmpQryUL;
      *refEndSIPtr = tmpRefUL;
   } /*Else: using the reference array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *refStartSIPtr += alnSTPtr->refOffsetSI;
   *refEndSIPtr += alnSTPtr->refOffsetSI;

   *qryStartSIPtr += alnSTPtr->qryOffsetSI;
   *qryEndSIPtr += alnSTPtr->qryOffsetSI;

   return 0;

   noCoord_fun05_sec04:;
      return 1;
   boundsErr_fun05_sec04:;
      return -1;
} /*getCoord_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun06: qryByRefSort_aln_memwaterScan
|   - sorts the query coordinates in a aln_memwaterScan
|     structure by the reference coordinate
|   - the sorting is first smallest query state, then
|     highest score, then longest (ending base) length
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
qryByRefSort_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC:
   '   - sorts the query coordinates in a aln_memwaterScan
   '     structure by the reference coordinate
   '   o fun06 sec01:
   '     - variable declerations
   '   o fun06 sec02:
   '     - find the number of rounds to sort for
   '   o fun06 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of elements to sort*/
   signed int numElmSI = alnSTPtr->qryLenSI - 1;

   /*Number of sorting rounds*/
   signed int subSI = 0;
   signed int nextSI = 0;
   signed int lastSI = 0;
   signed int onSI = 0;

   signed int onQrySI = 0;
   signed int compareQrySI = 0;

   /*Variables to incurment loops*/
   signed int siIndex = 0;
   signed int siElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numElmSI <= 0)
      return;

   /*Recursion formsla: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSI = 1; /*Initialzie first array*/

   while(subSI < numElmSI - 1)
      subSI = (3 * subSI) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^   - sort arrays
   ^   o fun06 sec03 sub01:
   ^     - outer loops for shell short
   ^   o fun06 sec03 sub02:
   ^     - compare and check if need to do first swap
   ^   o fun06 sec03 sub03:
   ^     - do inner (insertion) swaps
   ^   o fun06 sec03 sub04:
   ^     - move to the next element in the group
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec03 Sub01:
   *   - outer loops for shell short
   \*****************************************************/

   while(subSI > 0)
   { /*Loop: all rounds*/
      for(siIndex = 0; siIndex <= subSI; ++siIndex)
      { /*Loop: though sub array*/
         siElm = siIndex;

         for(
            siElm = siIndex;
            siElm + subSI <= numElmSI;
            siElm += subSI
         ){ /*Loop: swap elements in subarray*/
            nextSI = siElm + subSI;

            /********************************************\
            * Fun06 Sec03 Sub02:
            *   - compare & check if need to do first swap
            \********************************************/

            if(alnSTPtr->qryScoreArySL[siElm] < 0)
               goto swapOuter_fun06_sec03_sub02;
            else if(alnSTPtr->qryScoreArySL[nextSI] < 0)
               goto nextElm_fun06_sec03_sub04;
 
            onQrySI =
               qryCoord_indexToCoord(
                  alnSTPtr->qryLenSI,
                  alnSTPtr->qryStartArySL[siElm]
               );

            compareQrySI =
               qryCoord_indexToCoord(
                  alnSTPtr->qryLenSI,
                  alnSTPtr->qryStartArySL[nextSI]
               );

            if(onQrySI > compareQrySI)
               ;
            else if(onQrySI > compareQrySI)
               goto swapOuter_fun06_sec03_sub02;
            else if(
                 alnSTPtr->qryScoreArySL[siElm]
               < alnSTPtr->qryScoreArySL[nextSI]
            ) goto swapOuter_fun06_sec03_sub02;

            else if(
                  alnSTPtr->qryScoreArySL[siElm]
               == alnSTPtr->qryScoreArySL[nextSI]
            ){ /*Else If: had equal scores*/
                  onQrySI =
                     qryCoord_indexToCoord(
                        alnSTPtr->qryLenSI,
                        alnSTPtr->qryEndArySL[siElm]
                     );

                  compareQrySI =
                     qryCoord_indexToCoord(
                        alnSTPtr->qryLenSI,
                        alnSTPtr->qryEndArySL[nextSI]
                     );

               if(onQrySI > compareQrySI)
                  goto swapOuter_fun06_sec03_sub02;
            }  /*Else If: had equal scores*/

            else
            { /*Else: check if have length*/
               swapOuter_fun06_sec03_sub02:;

               alnSTPtr->qryStartArySL[siElm] ^=
                  alnSTPtr->qryStartArySL[nextSI];
               alnSTPtr->qryStartArySL[nextSI] ^=
                  alnSTPtr->qryStartArySL[siElm];
               alnSTPtr->qryStartArySL[siElm] ^=
                  alnSTPtr->qryStartArySL[nextSI];

               alnSTPtr->qryEndArySL[siElm] ^=
                  alnSTPtr->qryEndArySL[nextSI];
               alnSTPtr->qryEndArySL[nextSI] ^=
                  alnSTPtr->qryEndArySL[siElm];
               alnSTPtr->qryEndArySL[siElm] ^=
                  alnSTPtr->qryEndArySL[nextSI];

               alnSTPtr->qryScoreArySL[siElm] ^=
                  alnSTPtr->qryScoreArySL[nextSI];
               alnSTPtr->qryScoreArySL[nextSI] ^=
                  alnSTPtr->qryScoreArySL[siElm];
               alnSTPtr->qryScoreArySL[siElm] ^=
                  alnSTPtr->qryScoreArySL[nextSI];

               /*****************************************\
               * Fun06 Sec03 Sub03:
               *   - do inner (insertion) swaps
               \*****************************************/

               lastSI = siElm;
               onSI = siElm;

               while(lastSI >= subSI)
               { /*loop: move swapped element back*/
                  lastSI -= subSI;

                  if(alnSTPtr->qryScoreArySL[onSI] < 0)
                     break;
                  else if(
                     alnSTPtr->qryScoreArySL[lastSI] < 0
                  ) goto swapInner_fun06_sec03_sub03;
                   
                  onQrySI =
                     qryCoord_indexToCoord(
                        alnSTPtr->qryLenSI,
                        alnSTPtr->qryStartArySL[onSI]
                     );

                  compareQrySI =
                     qryCoord_indexToCoord(
                        alnSTPtr->qryLenSI,
                        alnSTPtr->qryStartArySL[lastSI]
                     );

                  if(onQrySI > compareQrySI)
                     break; /*Positioned the element*/
                  else if(onQrySI == compareQrySI)
                  { /*Else If: same reference position*/
                     if(
                          alnSTPtr->qryScoreArySL[onSI]
                        < alnSTPtr->qryScoreArySL[lastSI]
                     ) break;

                     else if(
                           alnSTPtr->qryScoreArySL[onSI]
                        == alnSTPtr->qryScoreArySL[lastSI]
                     ){ /*Else If: have same score*/

                        onQrySI =
                           qryCoord_indexToCoord(
                              alnSTPtr->qryLenSI,
                              alnSTPtr->qryEndArySL[onSI]
                           );

                        compareQrySI =
                           qryCoord_indexToCoord(
                             alnSTPtr->qryLenSI,
                             alnSTPtr->qryEndArySL[lastSI]
                           );
                        if(onQrySI <= compareQrySI)
                           break;
                     }  /*Else If: have same score*/
                  } /*Else If: same reference position*/

                  swapInner_fun06_sec03_sub03:;
                     alnSTPtr->qryStartArySL[onSI] ^=
                        alnSTPtr->qryStartArySL[lastSI];
                     alnSTPtr->qryStartArySL[lastSI] ^=
                        alnSTPtr->qryStartArySL[onSI];
                     alnSTPtr->qryStartArySL[onSI] ^=
                        alnSTPtr->qryStartArySL[lastSI];

                     alnSTPtr->qryEndArySL[onSI] ^=
                        alnSTPtr->qryEndArySL[lastSI];
                     alnSTPtr->qryEndArySL[lastSI] ^=
                        alnSTPtr->qryEndArySL[onSI];
                     alnSTPtr->qryEndArySL[onSI] ^=
                        alnSTPtr->qryEndArySL[lastSI];

                     alnSTPtr->qryScoreArySL[onSI] ^=
                        alnSTPtr->qryScoreArySL[lastSI];
                     alnSTPtr->qryScoreArySL[lastSI] ^=
                        alnSTPtr->qryScoreArySL[onSI];
                     alnSTPtr->qryScoreArySL[onSI] ^=
                        alnSTPtr->qryScoreArySL[lastSI];

                     onSI = lastSI;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/

            /********************************************\
            * Fun06 Sec03 Sub04:
            *   - move to the next element in the group
            \********************************************/

            nextElm_fun06_sec03_sub04:;
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSI = (subSI - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*qryByRefSort_aln_memwaterScan*/

/*-------------------------------------------------------\
| Fun07: memwaterScan
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
|    o score for aligment
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
memwaterScan(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct aln_memwaterScan *alnSTPtr,/*gets alignment*/
   struct alnSet *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC: memwaterAln
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun07 sec01:
   '    - Variable declerations
   '  o fun07 sec02:
   '    - Allocate memory for alignment
   '  o fun07 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun07 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01: Variable declerations
   ^  o fun07 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun07 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun07 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun07 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   signed long scoreSL = 0;     /*score to return*/

   /*Get start & end of query and reference sequences*/
   signed char *refSeqStr = 0;
   signed char *qrySeqStr = 0;

   signed long lenRefSL =
      refSTPtr->endAlnSL - refSTPtr->offsetSL + 1;

   signed long lenQrySL =
      qrySTPtr->endAlnSL - qrySTPtr->offsetSL + 1;

   /*Iterators for loops*/
   signed long slRef = 0;
   signed long slQry = 0;

   /*****************************************************\
   * Fun07 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   signed long snpScoreSL = 0;
   signed long insScoreSL = 0;
   signed long delScoreSL = 0;   /*score for deletion*/
   signed long nextSnpScoreSL = 0;/*next match/snp score*/

   /*****************************************************\
   * Fun07 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Keeping track of alignment starting positions*/
   signed long indexSL = 0;    /*index at in matrix*/
   signed long snpIndexSL = 0; /*last snp index*/
   signed long tmpIndexSL = 0; /*for getting snp index*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun07 sec02 sub01:
   ^    - allocate memory for the reference values
   ^  o fun07 sec02 sub02:
   ^    - allocate memory for the query values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec02 Sub01:
   *   - allocate memory for the reference values
   \****************************************************/

   if(alnSTPtr->rowSizeSI < refSTPtr->seqLenSL)
   { /*If: need more memory for the reference*/
      if(alnSTPtr->refStartArySL)
         free(alnSTPtr->refStartArySL);
      alnSTPtr->refStartArySL =
         malloc(refSTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->refStartArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->refEndArySL)
         free(alnSTPtr->refEndArySL);
      alnSTPtr->refEndArySL =
         malloc(refSTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->refEndArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->refScoreArySL)
         free(alnSTPtr->refScoreArySL);
      alnSTPtr->refScoreArySL =
         malloc(refSTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->refScoreArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->indexArySL)
         free(alnSTPtr->indexArySL);
      alnSTPtr->indexArySL =
         malloc(
            (refSTPtr->seqLenSL + 1) * sizeof(signed long)
         );
      if(! alnSTPtr->indexArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->scoreArySL)
         free(alnSTPtr->scoreArySL);
      alnSTPtr->scoreArySL =
         malloc(
            (refSTPtr->seqLenSL + 1) * sizeof(signed long)
         );
      if(! alnSTPtr->scoreArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->dirArySC)
         free(alnSTPtr->dirArySC);
      alnSTPtr->dirArySC =
         malloc(
            (refSTPtr->seqLenSL + 1) * sizeof(signed long)
         );
      if(! alnSTPtr->dirArySC)
         goto memErr_fun07_sec05_sub03;

      alnSTPtr->rowSizeSI = refSTPtr->seqLenSL;
   } /*If: need more memory for the reference*/

   /*****************************************************\
   * Fun07 Sec02 Sub02:
   *   - allocate memory for the query values
   \****************************************************/

   if(alnSTPtr->qrySizeSI < qrySTPtr->seqLenSL)
   { /*If: need more memory for the query*/
      if(alnSTPtr->qryStartArySL)
         free(alnSTPtr->qryStartArySL);
      alnSTPtr->qryStartArySL =
         malloc(qrySTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->qryStartArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->qryEndArySL)
         free(alnSTPtr->qryEndArySL);
      alnSTPtr->qryEndArySL =
         malloc(qrySTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->qryEndArySL)
         goto memErr_fun07_sec05_sub03;

      if(alnSTPtr->qryScoreArySL)
         free(alnSTPtr->qryScoreArySL);
      alnSTPtr->qryScoreArySL =
         malloc(qrySTPtr->seqLenSL * sizeof(signed long));
      if(! alnSTPtr->qryScoreArySL)
         goto memErr_fun07_sec05_sub03;

      alnSTPtr->qrySizeSI = qrySTPtr->seqLenSL;
   } /*If: need more memory for the query*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^  - initialize all values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnSTPtr->refOffsetSI = refSTPtr->offsetSL;
   alnSTPtr->qryOffsetSI = qrySTPtr->offsetSL;

   alnSTPtr->refLenSI = lenRefSL;
   alnSTPtr->qryLenSI = lenQrySL;

   for(indexSL = 0; indexSL < lenQrySL; ++indexSL)
   { /*loop; till have initalized the query scores*/
      alnSTPtr->qryStartArySL[indexSL] = -1;
      alnSTPtr->qryEndArySL[indexSL] = -1;
      alnSTPtr->qryScoreArySL[indexSL] = -1;
   } /*loop; till have initalized the query scores*/

   for(indexSL = 0; indexSL <= lenRefSL; ++indexSL)
   { /*loop; till have initalized the first row*/
      alnSTPtr->dirArySC[indexSL] = def_mvStop_alnDefs;
      alnSTPtr->indexArySL[indexSL] = indexSL + 1;
      alnSTPtr->scoreArySL[indexSL] = 0;

      if(indexSL < lenRefSL)
      { /*If: reference position*/
         alnSTPtr->refStartArySL[indexSL] = -1;
         alnSTPtr->refEndArySL[indexSL] = -1;
         alnSTPtr->refScoreArySL[indexSL] = -1;
      } /*If: reference position*/
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun07 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun07 sec04 sub02:
   ^    - get snp and ins scores + start loop
   ^  o fun07 sec04 sub03:
   ^    - find high score
   ^  o fun07 sec04 sub04:
   ^    - check if keep score (score > 0)
   ^  o fun07 sec04 sub05:
   ^    - find next deletion score and move to next index
   ^  o fun07 sec04 sub07:
   ^    - prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   alnSTPtr->indexArySL[0] = indexSL;

   /*These are always negative*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   snpIndexSL = alnSTPtr->indexArySL[0];

   /*incurment to frist base*/
   ++indexSL;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetSL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetSL - 1;
      /*offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun07 Sec04 Sub02:
   *  - get snp and ins scores + start loop
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(slQry = 0; slQry < lenQrySL; ++slQry)
   { /*loop; compare query base against all ref bases*/
      for(slRef = 1; slRef <= lenRefSL; ++slRef)
      { /*loop; compare one query to one reference base*/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[slQry],
               refSeqStr[slRef],
               settings
            ); /*find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = alnSTPtr->scoreArySL[slRef];

         /*insertion score*/
         #ifdef NOEXTEND
            insScoreSL =
                 alnSTPtr->scoreArySL[slRef]
               + settings->gapSS;
         #else
            insScoreSL = alnSTPtr->scoreArySL[slRef];
            insScoreSL +=
               settings->insArySS[
                  alnSTPtr->dirArySC[slRef]
               ];
         #endif

         /***********************************************\
         * Fun07 Sec04 Sub03:
         *   - find high score
         \***********************************************/

         alnSTPtr->scoreArySL[slRef] =
            max_genMath(insScoreSL, snpScoreSL);
            /*find if ins/snp is best (5 Op)*/

         tmpIndexSL = alnSTPtr->indexArySL[slRef];

         alnSTPtr->indexArySL[slRef] =
            ifmax_genMath(
               insScoreSL,
               snpScoreSL,
               alnSTPtr->indexArySL[slRef], /*ins index*/
               snpIndexSL                   /*snp index*/
            ); /*get index of high score*/

         snpIndexSL = tmpIndexSL;

         /*find direction (5 Op)*/
         alnSTPtr->dirArySC[slRef] =
            alnSTPtr->scoreArySL[slRef] > delScoreSL;
         alnSTPtr->dirArySC[slRef] +=
            (
                 (snpScoreSL <= insScoreSL)
               & alnSTPtr->dirArySC[slRef]
            );
         ++alnSTPtr->dirArySC[slRef];

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
         alnSTPtr->indexArySL[slRef] =
            ifmax_genMath(
               delScoreSL,
               alnSTPtr->scoreArySL[slRef],
               alnSTPtr->indexArySL[slRef-1],/*del index*/
               alnSTPtr->indexArySL[slRef]/*current best*/
            ); /*get index of high score*/

         alnSTPtr->scoreArySL[slRef] =
            max_genMath(
               delScoreSL,
               alnSTPtr->scoreArySL[slRef]
         ); /*find if del is best (5 Op)*/
            
         /***********************************************\
         * Fun07 Sec04 Sub04:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(alnSTPtr->scoreArySL[slRef] <= 0)
         {
            alnSTPtr->dirArySC[slRef] = 0;
            alnSTPtr->scoreArySL[slRef] = 0;
            alnSTPtr->indexArySL[slRef] = indexSL + 1;
            /*always one off for new index's*/
            /* branchless method is slower here*/
         }

         else
         { /*Else: check if have new high score*/
            scoreSL = alnSTPtr->scoreArySL[slRef];

            if(scoreSL > alnSTPtr->refScoreArySL[slRef-1])
            { /*If: new best reference position score*/
               alnSTPtr->refStartArySL[slRef - 1] =
                   alnSTPtr->indexArySL[slRef];
               alnSTPtr->refEndArySL[slRef - 1] = indexSL;
               alnSTPtr->refScoreArySL[slRef-1] = scoreSL;
            } /*If: new best reference position score*/

            if(scoreSL > alnSTPtr->qryScoreArySL[slQry])
            { /*If: new best query position score*/
               alnSTPtr->qryStartArySL[slQry] =
                   alnSTPtr->indexArySL[slRef];
               alnSTPtr->qryEndArySL[slQry] = indexSL;
               alnSTPtr->qryScoreArySL[slQry] = scoreSL;
            } /*If: new best query position score*/
         } /*Else: check if have new high score*/

         /***********************************************\
         * Fun07 Sec04 Sub05:
         *   - find next deletion score and move to next
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
                 alnSTPtr->scoreArySL[slRef]
               + settings->gapSS;
         #else
            delScoreSL = alnSTPtr->scoreArySL[slRef];
            delScoreSL +=
               settings->delArySS[
                  alnSTPtr->dirArySC[slRef]
               ];
         #endif

         ++indexSL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun07 Sec04 Sub07:
     *   - prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     alnSTPtr->indexArySL[0] = indexSL; /*next index*/
     snpIndexSL = alnSTPtr->indexArySL[0];

     ++indexSL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec05:
   ^  - set up for returing the matrix (clean up/wrap up)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto cleanUp_fun07_sec05;

   memErr_fun07_sec05_sub03:;
      scoreSL = -1;
      goto cleanUp_fun07_sec05;

   cleanUp_fun07_sec05:;
      return scoreSL;
} /*memwaterScan*/

#ifdef NEW_DEBUG
/*-------------------------------------------------------\
| Fun08: trsGetLongest_memwaterScan
|   - find longest trs (tandum repeat sequence) repeat
| Input;
|   - readSTPtr:
|     o pointer to an seqST with the read to scan for trs
|       index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - trsSeqSTPtr:
|     o pointer to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - startSIPtr:
|     o signed int pionter set to 1st reference base in
|       the trs
|   - endSIPtr:
|     o signed int pionter set to last reference base in
|       the trs
|   - alnSTPtr:
|     o pointer to aln_memwaterScan structure to hold the
|       results of the alignment (and used to align)
|   - settings:
|     o pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o variables in alnSTPtr to have the new alignment
|    o startSIPtr to have the 1st reference coorindate in
|      the longest found tandum repeat sequence
|    o endSIPtr to have the 1st reference coorindate in
|      the longest found tandum repeat sequence
|  - Returns:
|    o score of longest trs
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
trsGetLongest_memwaterScan(
   struct seqST *readSTPtr, /*read/sequence to scan*/
   struct seqST *trsSTPtr,  /*ref sequence and data*/
   signed int *startSIPtr,  /*gets start of longest map*/
   signed int *endSIPtr,    /*gets end of longest map*/
   struct aln_memwaterScan *alnSTPtr,/*gets alignment*/
   struct alnSet *settings
){
   signed long scoreSL = 0;

   unsigned long refCoordAryUL[2]; /*0 = start; 1 = end*/
   unsigned long qryCoordAryUL[2]; /*0 = start; 1 = end*/

   unsigned long lastRefCoordAryUL[2];/*0=start; 1=end*/
   unsigned long lastQryCoordAryUL[2];/*0=start; 1=end*/
 
   signed int trsLenSI = 0;
   signed int indexSI = 0;
   signed int nextSI = 0;
   signed int qryNextSI = 0;
   signed int endSI = 0;
   signed short gapSS = trsSTPtr->seqLenUL / 10;
      /*limits to a 10% gap*/

   signed short indelSS = trsSTPtr->seqLenUL / 10;
      /*limits for gaps between bases*/

   *startSIPtr = 0;
   *endSIPtr = 0;

   scoreSL =
      memwaterScan(readSTPtr,trsSTPtr,alnSTPtr,settings);
   if(scoreSL < 0)
      goto memErr_fun08_sec0x;

   for(indexSI=0; indexSI < alnSTPtr->refLenSI; ++indexSI)
      alnSTPtr->scoreArySL[indexSI] = 0;

   for(indexSI=0; indexSI < alnSTPtr->refLenSI; ++indexSI)
   { /*Loop: find longest length*/
      if(
           settings->minScoreSL
         < alnSTPtr->qryScoreSL[indexSI]
      ){ /*If: no valid score*/
         alnSTPtr->scoreArySL[indexSI] = -2;
         continue;
      }  /*If: no valid score*/

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->qryStartArySL[indexSI],
         &lastRefCoordAryUL[0],
         &lastQryCoordAryUL[0]
      );

      indexToCoord(
         alnSTPtr->refLenSI,
         alnSTPtr->qryEndArySL[indexSI],
         &lastRefCoordAryUL[1],
         &lastQryCoordAryUL[1]
      );

      trsLenSI =
         lastRefCoordAryUL[1] - lastRefCoordAryUL[0];
      nextSI = lastRefCoordAryUL[1] + 1;
      endSI = nextSI + gapSS;

      while(nextSI <= endSI)
      { /*Loop: see if can extend trs*/
         indexToCoord(
            alnSTPtr->refLenSI,
            alnSTPtr->qryStartArySL[indexSI],
            &refCoordAryUL[0],
            &qryCoordAryUL[0]
         );

         indexToCoord(
            alnSTPtr->refLenSI,
            alnSTPtr->qryEndArySL[indexSI],
            &refCoordAryUL[1],
            &qryCoordAryUL[1]
         );

         if(
                (signed long) refCoordAryUL[0]
              - lastRefCoordUL[0]
            > gapUS
         ){ /*If: end of the repeat*/
            alnSTPtr->scoreArySL[indexSI] = trsLenSI;

            nextSI = lastRefCoordAry[0] + 1;
            qryNextSI = lastQryCoordAry[0] + 1;
            while(nextSI < lastRefCoordAry[1])
            { /*Loop: fill in duplicate entries*/
               indexToCoord(
                  alnSTPtr->refLenSI,
                  alnSTPtr->qryStartArySL[indexSI],
                  &refCoordAryUL[0],
                  &qryCoordAryUL[0]
               );

               if(qryCoordAryUL[0] < qryNextSI - indelSS)
                  ;
               else if(
                  qryCoordAryUL[0] > qryNextSI + indelSS
               ) ;
               else
                  alnSTPtr->scoreArySL[nextSI] = -1;

               ++nextSI;
               ++qryNextSI;
            } /*Loop: fill in duplicate entries*/

            break;
         } /*If: end of the repeat*/

         else if(
            qryCoordAryUL[0] - lastEndAryUL[1] > gapUS
         ){ /*Else If: end of the repeat*/
         } /*Else If: end of the repeat*/

      } /*Loop: see if can extend trs*/

      curStartSI = 
   } /*Loop: find longest length*/

   return scoreSL;

   memErr_fun08_sec0x:;
      return scoreSL;
} /*trsGetLongest_memwaterScan*/
#endif

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
