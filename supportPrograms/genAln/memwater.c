/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memwater SOF: Start Of Functions
'   - has memory efficent Waterman Smith alignment
'   o header:
'     - Included libraries
'   o fun01 memwater:
'     - Run a memory efficent Waterman Smith alignment on
'       input sequences
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "memwater.h"

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
| Fun01: memwater
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qrySeqSTVoidPtr:
|     o Point to an seqST with the query sequence and
|       index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - refSeqSTVoidPtr:
|     o Point to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - refStartSL:
|     o Pointer to unsigned long to hold the frist
|       reference base in the alignment
|   - refEndSL:
|     o Pointer to unsigned long to hold the last
|       reference base in the alignment
|   - qryStartSL:
|     o Pointer to unsigned long to hold the frist query 
|       base in the alignment
|   - qryEndSL:
|     o Pointer to unsigned long to hold the last query
|       base in the alignment
|   - alnSetVoidPtr:
|     o Pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o refStartSL to have 1st reference base in alignment
|    o refEndSL to have last reference base in alignment
|    o qryStartSL to have first query base in alignment
|    o qryEndSL to have last query base in alignment
|  - Returns:
|    o score for aligment
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
memwater(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   signed long *refStartSL,
   signed long *refEndSL,
   signed long *qryStartSL,
   signed long *qryEndSL,
   struct alnSet *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: memwaterAln
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun01 sec01:
   '    - Variable declerations
   '  o fun01 sec02:
   '    - Allocate memory for alignment
   '  o fun01 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun01 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01: Variable declerations
   ^  o fun01 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun01 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun01 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun01 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   signed long scoreSL = 0;     /*score to return*/
   signed long bestStartSL = 0; /*best starting index*/
   signed long bestEndSL = 0;   /*best ending index*/

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
   * Fun01 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   signed long snpScoreSL = 0;
   signed long insScoreSL = 0;
   signed long delScoreSL = 0;   /*score for deletion*/
   signed long nextSnpScoreSL = 0;/*next match/snp score*/
   signed long *scoreHeapArySL = 0; /*alignment matrix*/

   /*****************************************************\
   * Fun01 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Direction matrix (1 cell holds a single direction)*/
   signed char *dirRowHeapSC = 0;  /*Holds directions*/

   /*Keeping track of alignment starting positions*/
   signed long indexSL = 0;    /*index at in matrix*/
   signed long snpIndexSL = 0; /*last snp index*/
   signed long tmpIndexSL = 0; /*for getting snp index*/

   signed long *indexHeapSL=0; /*row; starting indexes*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun01 sec02 sub01:
   ^    - Allocate memory for the alignment
   ^  o fun01 sec02 sub02:
   ^    - Allocate memory for keeping track of indexes
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *  - Allocate memory for the alignment
   \****************************************************/

   dirRowHeapSC =
      calloc((lenRefSL + 1), sizeof(signed char));
   if(dirRowHeapSC == 0)
      goto memErr_fun01_sec05_sub03;

   scoreHeapArySL =
      calloc((lenRefSL + 1), sizeof(signed long));
      /*+ 1 is for the indel column*/
   if(scoreHeapArySL == 0)
      goto memErr_fun01_sec05_sub03;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *  - Get memory for keeping track of starting indexes
   \*****************************************************/

   /*Set up the first row of starting indexes*/
   indexHeapSL =
      calloc((lenRefSL + 1), sizeof(unsigned long));
   if(indexHeapSL == 0)
      goto memErr_fun01_sec05_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      indexSL = 0;
      indexSL <= lenRefSL;
      ++indexSL
   ){ /*loop; till have initalized the first row*/
      dirRowHeapSC[indexSL] = def_mvStop_alnDefs;
      indexHeapSL[indexSL] = indexSL + 1;
      scoreHeapArySL[indexSL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun01 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun01 sec04 sub02:
   ^    - get snp and ins scores + start loop
   ^  o fun01 sec04 sub03:
   ^    - find high score
   ^  o fun01 sec04 sub04:
   ^    - check if keep score (score > 0)
   ^  o fun01 sec04 sub05:
   ^    - find deletion score and move to next index
   ^  o fun01 sec04 sub07:
   ^    - prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   indexHeapSL[0] = indexSL;

   /*These are always negative*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   snpIndexSL = indexHeapSL[0];

   /*incurment to frist base*/
   ++indexSL;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetSL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetSL - 1;
      /*offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *  - get snp and ins scores + start loop
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      slQry = 0;
      slQry < lenQrySL;
      ++slQry
   ){ /*loop; compare query base against all ref bases*/

      for(
         slRef = 1;
         slRef <= lenRefSL;
         ++slRef
      ){ /*loop; compare one query to one reference base*/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[slQry],
               refSeqStr[slRef],
               settings
            ); /*find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = scoreHeapArySL[slRef];

         /*insertion score*/
         #ifdef NOEXTEND
            insScoreSL =
               scoreHeapArySL[slRef] + settings->gapSS;
         #else
            insScoreSL = scoreHeapArySL[slRef];
            insScoreSL +=
               settings->insArySS[dirRowHeapSC[slRef]];
         #endif

         /***********************************************\
         * Fun02 Sec04 Sub03:
         *   - find high score
         \***********************************************/

         scoreHeapArySL[slRef] =
            max_genMath(
               insScoreSL,
               snpScoreSL
         ); /*find if ins/snp is best (5 Op)*/

         tmpIndexSL = indexHeapSL[slRef];

         indexHeapSL[slRef] =
            ifmax_genMath(
               insScoreSL,
               snpScoreSL,
               indexHeapSL[slRef], /*insertion index*/
               snpIndexSL          /*snp index*/
            ); /*get index of high score*/

         snpIndexSL = tmpIndexSL;

         /*find direction (5 Op)*/
         dirRowHeapSC[slRef] =
            scoreHeapArySL[slRef] > delScoreSL;
         dirRowHeapSC[slRef] +=
            (
                 (snpScoreSL <= insScoreSL)
               & dirRowHeapSC[slRef]
            );
         ++dirRowHeapSC[slRef];

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
         indexHeapSL[slRef] =
            ifmax_genMath(
               delScoreSL,
               scoreHeapArySL[slRef],
               indexHeapSL[slRef - 1],    /*del index*/
               indexHeapSL[slRef]         /*current best*/
            ); /*get index of high score*/

         scoreHeapArySL[slRef] =
            max_genMath(
               delScoreSL,
               scoreHeapArySL[slRef]
         ); /*find if del is best (5 Op)*/
            
         /***********************************************\
         * Fun02 Sec04 Sub04:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(scoreHeapArySL[slRef] <= 0)
         {
            dirRowHeapSC[slRef] = 0;
            scoreHeapArySL[slRef] = 0;
            indexHeapSL[slRef] = indexSL + 1;
            /*always one off for new index's*/
         }

         else if(scoreSL < scoreHeapArySL[slRef])
         { /*Else If: this is the best score*/
            scoreSL = scoreHeapArySL[slRef];
            bestStartSL = indexHeapSL[slRef];
            bestEndSL = indexSL;
         } /*Else If: this was an snp or match*/


         /* branchless method is slower here
         keepSL =
             (signed long) -(scoreHeapArySL[slRef] > 0);
         dirRowHeapSC[slRef] &= keepSL;
         scoreHeapArySL[slRef] &= keepSL;

         indexHeapSL[slRef] &= keepSL;
         indexHeapSL[slRef] |= ( indexSL & (~keepSL) );
         */

         /***********************************************\
         * Fun02 Sec04 Sub05:
         *   - find deletion score and move to next index
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
               scoreHeapArySL[slRef] + settings->gapSS;
         #else
            delScoreSL = scoreHeapArySL[slRef];
            delScoreSL +=
               settings->delArySS[dirRowHeapSC[slRef]];
         #endif

         ++indexSL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun01 Sec04 Sub07:
     *   - prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     indexHeapSL[0] = indexSL; /*next index*/
     snpIndexSL = indexHeapSL[0];

     ++indexSL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun01 sec05 sub01:
   ^    - get coordinates
   ^  o fun01 sec05 sub02:
   ^    - no error clean up
   ^  o fun01 sec05 sub03:
   ^    - memory error clean up
   ^  o fun01 sec05 sub04:
   ^    - general clean up (no error or error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec05 Sub01:
   *  - get coordinates
   \*****************************************************/

   *refStartSL =
      refCoord_indexToCoord(
         lenRefSL,
         bestStartSL
      ); /*Convert the starting index to coordinates*/

   *qryStartSL =
      qryCoord_indexToCoord(
         lenRefSL,
         bestStartSL
      ); /*Convert the starting index to coordinates*/

   *refStartSL += refSTPtr->offsetSL;
   *qryStartSL += qrySTPtr->offsetSL;

   *refEndSL =
      refCoord_indexToCoord(
         lenRefSL,
         bestEndSL
      ); /*Convert ending index to coordinates*/

   *qryEndSL =
      qryCoord_indexToCoord(
         lenRefSL,
         bestEndSL
      ); /*Convert ending index to coordinates*/

   *refEndSL += refSTPtr->offsetSL;
   *qryEndSL += qrySTPtr->offsetSL;

   /*****************************************************\
   * Fun01 Sec05 Sub02:
   *  - no error clean up
   \*****************************************************/

   goto cleanUp_fun01_sec05_sub04;

   /*****************************************************\
   * Fun01 Sec05 Sub03:
   *  - memory error clean up
   \*****************************************************/

   memErr_fun01_sec05_sub03:;
      scoreSL = -1;
      goto cleanUp_fun01_sec05_sub04;

   /*****************************************************\
   * Fun01 Sec05 Sub04:
   *  - general clean up (no error or error)
   \*****************************************************/

   cleanUp_fun01_sec05_sub04:;
      free(dirRowHeapSC);
      dirRowHeapSC = 0;

      free(scoreHeapArySL);
      scoreHeapArySL = 0;

      free(indexHeapSL);
      indexHeapSL = 0;

      return scoreSL;
} /*memwater*/

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
