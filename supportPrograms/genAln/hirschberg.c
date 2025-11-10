/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' hirschberg SOF: Start Of File
'   o header:
'     - included libraries
'   o .c fun01: alnOneNt_hirschberg
'     - align a single base to a sequence
'   o .c fun02: forwardScore_hirschberg
'     - scores the input query and reference sequence in a
'        forward direction
'   o .c fun03: reverseScore_hirschberg
'     - scores input query and reference sequence
'       backwards for a hirschberg
'   o .c fun04: recursive_hirschberg
'     - does the recurisve part of a hirschberg alignment
'   o .c fun05: alnToSam_hirschberg
'     - converts an alignment from recursive_hirschberg to
'       a sam file entry
'   o fun06: aln_hirschberg
'     - does hischberg alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "../genLib/ulCp.h"
#include "../genLib/numToStr.h"

#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"
#include "hirschberg.h"
#include "alnSet.h"

/*.h files only*/
#include "../genLib/genMath.h"
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - std #include <stdio.h>
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "../genLib/endLine.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*these are for the alignment arrays, not the final
`   output
*/
#define def_end_hirschberg 0
#define def_gap_hirschberg 1
#define def_nt_hirschberg 2

/*-------------------------------------------------------\
| Fun01: alnOneNt_hirschberg
|  - align a single base to a sequence
| Input:
|  - baseSC:
|    o Base to align to a sequence
|  - baseIndexSL:
|    o Position of base in its seqequence (index 0)
|  - seqStr:
|    o Sequence to align the base to
|  - startOfSeqSL:
|    o Fist base to align in seqStr (index 0)
|  - lenSeqSL:
|    o Number of bases to align in seqStr (index 1)
|  - baseAlnStr:
|    o Holds the alignment of baseC's sequence
|  - seqAlnStr:
|    o Holds the alignment of seqStr
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o baseAlnStr to hold the alignment for the single
|      base (baseC)
|    o seqAlnStr to hold the alignment for the sequence
|      (seqStr)
\-------------------------------------------------------*/
void
alnOneNt_hirschberg(
   signed char baseSC,      /*single base to align*/
   signed long baseIndexSL, /*index baseSC is at*/
   signed char *seqStr,     /*seq to align baseSC to*/
   signed long startOfSeqSL,/*1st base to align; seqStr*/
   signed long lenSeqSL,    /*number bases in seqStr*/
   signed char *baseAlnStr, /*holds baseC alignment*/
   signed char *seqAlnStr,  /*holds seqStr alignment*/
   struct alnSet *settings  /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '  - Align a single base to a sequence
   '  o fun01 sec01:
   '    - variable declerations
   '  o fun01 sec02:
   '    - find the reference bases position on the query
   '  o fun01 sec03:
   '    - fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^  - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long insScoreSL = 0;   /*score for a insertion*/
   signed long matchScoreSL = 0; /*score for a match*/
   signed long lastInsSL = 0;    /*score for next snp*/

   signed long slBase = 0;       /*for the for loop*/
   signed long snpIndexSL = 0;   /*index of last match*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^  - find the reference bases position on the query
   ^  o fun01 sec02 sub01:
   ^    - find the first scores for the loop
   ^  o fun01 sec02 sub02:
   ^    - find the remaing scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *  - find the first scores for the loop
   \*****************************************************/

   lastInsSL = 0;
   insScoreSL = settings->gapSS;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *  - Find the remaing scores
   \*****************************************************/
 
   for(
      slBase = startOfSeqSL;
      slBase < startOfSeqSL + lenSeqSL;
      ++slBase
   ){ /*Loop: align baseC to seqStr*/
     matchScoreSL =
         lastInsSL
       + getScore_alnSet(seqStr[slBase],baseSC,settings);

     #ifdef NOEXTEND
        insScoreSL += settings->gapSS;
     #else
        insScoreSL += settings->extendSS;
     #endif

     lastInsSL = insScoreSL;

     if(matchScoreSL > insScoreSL)
        snpIndexSL = slBase;
   } /*Loop: align baseC to seqStr*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /* A series of deletions and insertions are prefered
   `  over matches and smps. In this case put the base
   `  at the first base. There is no good position
   */
   if(snpIndexSL == 0)
      snpIndexSL = startOfSeqSL;

   /*add in the insertions at the start*/
   for(
      slBase = startOfSeqSL;
      slBase < snpIndexSL;
      ++slBase
   ) seqAlnStr[slBase] = def_gap_hirschberg;
   
   /*Add in the position of the base*/
   baseAlnStr[baseIndexSL] = def_nt_hirschberg;
   seqAlnStr[slBase] = def_nt_hirschberg;

   /*Finish adding in the insertions at the end*/
   for(
      slBase = slBase + 1;
      slBase < startOfSeqSL + lenSeqSL;
      ++slBase
   ) seqAlnStr[slBase] = def_gap_hirschberg;

   return;
} /*alnOneNt_hirschberg*/

/*-------------------------------------------------------\
| Fun02: forwardScore_hirschberg
|  - scores the input query and reference sequence in a
|    forward direction
| Input:
|  - refStr:
|    o C-string with the reference sequence to align
|  - refStartSL:
|    o Position of the first base to align in the
|      reference sequence (is index 0)
|  - refLenSL:
|    o Number of bases to align in the reference sequence
|      (is index 1)
|  - qryStr:
|    o C-string with the query sequence to align
|  - qryStartSL:
|    o Position of the first base to align in the query
|      sequence (is index 0)
|  - qryLenSL:
|    o Number of bases to align in the query sequence
|      (is index 1)
|  - scoreArySC:
|    o Array to fill with scores
|  - dirArySC:
|    o array to use in finding directions
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o scoreArySC to hold the last row of scores in a
|      Needleman Wunsch
|    o dirArySC to hold the last row of directions in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
signed long
forwardScore_hirschberg(
  signed char *refStr,    /*reference sequence*/
  signed long refStartSL, /*index of first ref base*/
  signed long refLenSL,   /*length of reference*/

  signed char *qryStr,    /*query sequence*/
  signed long qryStartSL, /*index of first query base*/
  signed long qryLenSL,   /*length of query*/

  signed long *scoreArySL,/*array of scores to fill*/
  signed char *dirArySC,  /*direction row (gap extend)*/
  struct alnSet *alnSetPtr /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC: forwardScore_hirschberg
   '  - does the scoring step for a hirschberg alignment
   '    (forward direction)
   '  o fun02 sec01:
   '    - variable declerations
   '  o fun02 sec02:
   '    - set up the first row (indel row) of scores
   '  o fun02 sec03:
   '    - find the scores
   '  o fun02 sec04:
   '    - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^  - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long refNtSL = 0;
   signed long qryNtSL = 0;
   signed long gapColSL = 0;
   signed long nextNtScoreSL = 0;

   signed long delScoreSL = 0;
   signed long insScoreSL = 0;
   signed long ntScoreSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^  - set up the first row (indel row) of scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   gapColSL = 0;
   scoreArySL[refStartSL] = alnSetPtr->gapSS;
   dirArySC[refStartSL] = -1;

   for(
      refNtSL = refStartSL + 1;
      refNtSL < refStartSL + refLenSL;
      ++refNtSL
   ){ /*Loop:Set the initial blank scores*/
      dirArySC[refNtSL] = -1;

      #ifdef NOEXTEND
         scoreArySL[refNtSL] =
            scoreArySL[refNtSL - 1] + alnSetPtr->gapSS;
      #else
         scoreArySL[refNtSL] =
            scoreArySL[refNtSL - 1] + alnSetPtr->extendSS;
      #endif
   } /*Loop:Set the initial blank scores*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - find the scores
   ^   o fun02 sec03 sub01:
   ^     - setup initial scores and start scoring loops
   ^   o fun02 sec03 sub02:
   ^     - find snp and insertion scores
   ^   o fun02 sec03 sub03:
   ^     - find the best score
   ^   o fun02 sec03 sub04:
   ^     - find next deletion score
   ^   o fun02 sec03 sub05:
   ^     - find score at rows end
   ^   o fun02 sec03 sub06:
   ^     - prepare for the next round
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec03 Sub01:
   *   - setup initial scores and start scoring loops
   \*****************************************************/

   /*Get the value I need to find the next snp's score*/
   nextNtScoreSL = 0;
   gapColSL = alnSetPtr->gapSS;

   #ifdef NOEXTEND
      delScoreSL = gapColSL + alnSetPtr->gapSS;
   #else
      delScoreSL = gapColSL + alnSetPtr->extendSS;
   #endif

   for(
      qryNtSL = qryStartSL;
      qryNtSL < qryStartSL + qryLenSL;
      ++qryNtSL
   ){ /*Loop: score all query bases (rows)*/

      for(
         refNtSL = refStartSL;
         refNtSL < refStartSL + refLenSL - 1;
         ++refNtSL
      ){ /*Loop: Find the max scores for a single row*/

        /************************************************\
        * Fun02 Sec03 Sub02:
        *   - find snp and insertion scores
        \************************************************/

        ntScoreSL =
           getScore_alnSet(
              refStr[refNtSL],
              qryStr[qryNtSL],
              alnSetPtr
           ); /*find score of aligned bases*/
        ntScoreSL += nextNtScoreSL;

        nextNtScoreSL = scoreArySL[refNtSL];
           /*snp score for next round*/

        #ifdef NOEXTEND
           insScoreSL = scoreArySL[refNtSL];
           insScoreSL += alnSetPtr->gapSS;
        #else
           insScoreSL = scoreArySL[refNtSL];
           insScoreSL +=
              alnSetPtr->insArySS[ dirArySC[refNtSL] ];
              /*find gap extension or open & add penalty*/
        #endif

        /************************************************\
        * Fun02 Sec03 Sub03:
        *   - find the best score
        \************************************************/

         scoreArySL[refNtSL] =
            max_genMath(insScoreSL, ntScoreSL);
            /*find if ins/snp is best (5 Op)*/
 
         /*find direction (5 Op)*/
         dirArySC[refNtSL] =
            scoreArySL[refNtSL] > delScoreSL;
         dirArySC[refNtSL] +=
            (
                 (ntScoreSL <= insScoreSL)
               & dirArySC[refNtSL]
            );
         ++dirArySC[refNtSL];

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

         /*finish max*/
         scoreArySL[refNtSL] =
            max_genMath(delScoreSL, scoreArySL[refNtSL]);
            /*find if del is best (5 Op)*/

        /************************************************\
        * Fun02 Sec03 Sub04:
        *   - find next deletion score
        \************************************************/

        #ifdef NOEXTEND
           delScoreSL = scoreArySL[refNtSL];
           delScoreSL += alnSetPtr->gapSS;
        #else
           delScoreSL = scoreArySL[refNtSL];
           delScoreSL +=
              alnSetPtr->delArySS[ dirArySC[refNtSL] ];
              /*find gap extension or open & add penalty*/
        #endif
      } /*Loop: Find the max scores for a single row*/

      /**************************************************\
      * Fun02 Sec03 Sub05:
      *   - find score at rows end
      *   o fun02 sec03 sub05 cat01:
      *     - find the sNP and insertion score
      *   o fun02 sec03 sub05 cat02:
      *     - find the best score
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec03 Sub05 Cat01:
      +   - find the SNP and insertion score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      ntScoreSL =
         getScore_alnSet(
            refStr[refNtSL],
            qryStr[qryNtSL],
            alnSetPtr
         ); /*find score of aligned bases*/
      ntScoreSL += nextNtScoreSL;

      #ifdef NOEXTEND
         insScoreSL = scoreArySL[refNtSL];
         insScoreSL += alnSetPtr->gapSS;
      #else
         insScoreSL = scoreArySL[refNtSL];
         insScoreSL +=
            alnSetPtr->insArySS[ dirArySC[refNtSL] ];
            /*find gap extension or open & add penalty*/
      #endif

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec03 Sub05 Cat02:
      +   - find the best score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      scoreArySL[refNtSL] =
         max_genMath(insScoreSL, ntScoreSL);
         /*find if ins/snp is best (5 Op)*/
 
      /*find direction (5 Op)*/
      dirArySC[refNtSL] =
         scoreArySL[refNtSL] > delScoreSL;
      dirArySC[refNtSL] +=
         ((ntScoreSL <= insScoreSL) & dirArySC[refNtSL]);
      ++dirArySC[refNtSL];

      scoreArySL[refNtSL] =
         max_genMath(delScoreSL, scoreArySL[refNtSL]);
         /*find if del is best (5 Op)*/

      /**************************************************\
      * Fun02 Sec03 Sub06:
      *   - prepare for the next round
      \**************************************************/

      nextNtScoreSL = gapColSL;

      gapColSL += alnSetPtr->gapSS;
         /*not applying extensions because I do not want
         `  to apply an extension penalty when an
         `  insertion changes to a deletion
         */

      #ifdef NOEXTEND
         delScoreSL = gapColSL + alnSetPtr->gapSS;
      #else
         delScoreSL = gapColSL + alnSetPtr->extendSS;
      #endif
   } /*Loop: score all query bases (rows)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^  - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   gapColSL -= alnSetPtr->gapSS;
      /*correct for being on the next row*/

   return gapColSL;
} /*forwardScore_hirschberg*/

/*-------------------------------------------------------\
| Fun03: reverseScore_hirschberg
|  - scores input query and reference sequence backwards
|    for a hirschberg
| Input:
|  - refStr:
|    o C-string with the reference sequence to align
|  - refStartSL:
|    o Position of the first base to align in the
|      reference sequence (is index 0)
|  - refLenSL:
|    o Number of bases to align in the reference sequence
|      (is index 1)
|  - qryStr:
|    o C-string with the query sequence to align
|  - qryStartSL:
|    o Position of the first base to align in the query
|      sequence (is index 0)
|  - qryLenSL:
|    o Number of bases to align in the query sequence
|      (is index 1)
|  - scoreArySC:
|    o Array to fill with scores
|  - dirArySC:
|    o array to use in finding directions
|  - alnSetPtr:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o scoreArySC to hold the last row of scores in a
|      Needleman Wunsch
|    o dirArySC to hold the last row of directions in a
|      Needleman Wunsch
|  - Returns:
|    o The gap column score
\-------------------------------------------------------*/
signed long
reverseScore_hirschberg(
  signed char *refStr,    /*reference sequence*/
  signed long refStartSL, /*index of first ref base*/
  signed long refLenSL,   /*length of reference*/

  signed char *qryStr,    /*query sequence*/
  signed long qryStartSL, /*index of first query base*/
  signed long qryLenSL,   /*length of query*/

  signed long *scoreArySL, /*array of scores to fill*/
  signed char *dirArySC,  /*direction row (gap extend)*/
  struct alnSet *alnSetPtr /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC: scoreHirschRev
   '  - does a single round of scoring for a hirschberg
   '    alignment (backwards direction)
   '  o fun03 sec01:
   '    - variable declerations
   '  o fun03 sec02:
   '    - set up the first row (indel row) of scores
   '  o fun03 sec03:
   '    - find the scores
   '  o fun03 sec04:
   '    - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^  - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long refNtSL = 0;
   signed long qryNtSL = 0;
   signed long gapColSL = 0;
   signed long nextNtScoreSL = 0;

   signed long delScoreSL = 0;
   signed long insScoreSL = 0;
   signed long ntScoreSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^  - set up the first row (indel row) of scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   dirArySC[refStartSL + refLenSL - 1] = -1;
   gapColSL = 0;
   scoreArySL[refStartSL + refLenSL - 1] =
      alnSetPtr->gapSS;

   for(
      refNtSL = refStartSL + refLenSL - 2;
      refNtSL > refStartSL;
      --refNtSL
   ){ /*Loop:Set the initial blank scores*/
      dirArySC[refNtSL] = -1;

      #ifdef NOEXTEND
         scoreArySL[refNtSL] =
            scoreArySL[refNtSL + 1] + alnSetPtr->gapSS;
      #else
         scoreArySL[refNtSL] =
            scoreArySL[refNtSL + 1] + alnSetPtr->extendSS;
      #endif
   } /*Loop:Set the initial blank scores*/


   dirArySC[refNtSL] = -1;

   #ifdef NOEXTEND
      scoreArySL[refNtSL] =
         scoreArySL[refStartSL + 1] + alnSetPtr->gapSS;
   #else
      scoreArySL[refNtSL] =
         scoreArySL[refStartSL + 1] + alnSetPtr->extendSS;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^  - find the scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Get the value I need to find the next snp's score*/
   nextNtScoreSL = 0;
   gapColSL = alnSetPtr->gapSS;

   #ifdef NOEXTEND
      delScoreSL = gapColSL + alnSetPtr->gapSS;
   #else
      delScoreSL = gapColSL + alnSetPtr->extendSS;
   #endif

   for(
      qryNtSL = qryStartSL + qryLenSL - 1;
      qryNtSL >= qryStartSL;
      --qryNtSL
   ){ /*Loop: score all query bases (rows)*/

      for(
         refNtSL = refStartSL + refLenSL - 1;
         refNtSL > refStartSL;
         --refNtSL
      ){ /*Loop: Find the max scores for a single row*/

        /************************************************\
        * Fun03 Sec03 Sub02:
        *   - find snp and insertion scores
        \************************************************/

        ntScoreSL =
           getScore_alnSet(
              refStr[refNtSL],
              qryStr[qryNtSL],
              alnSetPtr
           ); /*find score of aligned bases*/
        ntScoreSL += nextNtScoreSL;

        nextNtScoreSL = scoreArySL[refNtSL];
           /*snp score for next round*/

        #ifdef NOEXTEND
           insScoreSL = scoreArySL[refNtSL];
           insScoreSL += alnSetPtr->gapSS;
        #else
           insScoreSL = scoreArySL[refNtSL];
           insScoreSL +=
              alnSetPtr->insArySS[ dirArySC[refNtSL] ];
              /*find gap extension or open & add penalty*/
        #endif

        /************************************************\
        * Fun03 Sec03 Sub03:
        *   - find the best score
        \************************************************/

         scoreArySL[refNtSL] =
            max_genMath(insScoreSL, ntScoreSL);
            /*find if ins/snp is best (5 Op)*/
 
         /*find direction (5 Op)*/
         dirArySC[refNtSL] =
            scoreArySL[refNtSL] > delScoreSL;
         dirArySC[refNtSL] +=
            (
                 (ntScoreSL <= insScoreSL)
               & dirArySC[refNtSL]
            );
         ++dirArySC[refNtSL];

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

         /*finish max*/
         scoreArySL[refNtSL] =
            max_genMath(delScoreSL, scoreArySL[refNtSL]);
            /*find if del is best (5 Op)*/

        /************************************************\
        * Fun03 Sec03 Sub04:
        *   - find next deletion score
        \************************************************/

        #ifdef NOEXTEND
           delScoreSL = scoreArySL[refNtSL];
           delScoreSL += alnSetPtr->gapSS;
        #else
           delScoreSL = scoreArySL[refNtSL];
           delScoreSL +=
              alnSetPtr->delArySS[ dirArySC[refNtSL] ];
              /*find gap extension or open & add penalty*/
        #endif
      } /*Loop: Find the max scores for a single row*/

      /**************************************************\
      * Fun03 Sec03 Sub05:
      *   - find score at rows end
      *   o fun03 sec03 sub05 cat01:
      *     - find the sNP and insertion score
      *   o fun03 sec03 sub05 cat02:
      *     - find the best score
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub05 Cat01:
      +   - find the SNP and insertion score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      ntScoreSL =
         getScore_alnSet(
            refStr[refNtSL],
            qryStr[qryNtSL],
            alnSetPtr
         ); /*find score of aligned bases*/
      ntScoreSL += nextNtScoreSL;

      #ifdef NOEXTEND
         insScoreSL = scoreArySL[refNtSL];
         insScoreSL += alnSetPtr->gapSS;
      #else
         insScoreSL = scoreArySL[refNtSL];
         insScoreSL +=
            alnSetPtr->insArySS[ dirArySC[refNtSL] ];
            /*find gap extension or open & add penalty*/
      #endif

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub05 Cat02:
      +   - find the best score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      scoreArySL[refNtSL] =
         max_genMath(insScoreSL, ntScoreSL);
         /*find if ins/snp is best (5 Op)*/
 
      /*find direction (5 Op)*/
      dirArySC[refNtSL] =
         scoreArySL[refNtSL] > delScoreSL;
      dirArySC[refNtSL] +=
         ((ntScoreSL <= insScoreSL) & dirArySC[refNtSL]);
      ++dirArySC[refNtSL];

      scoreArySL[refNtSL] =
         max_genMath(delScoreSL, scoreArySL[refNtSL]);
         /*find if del is best (5 Op)*/

      /**************************************************\
      * Fun03 Sec03 Sub06:
      *   - prepare for the next round
      \**************************************************/

      nextNtScoreSL = gapColSL;

      gapColSL += alnSetPtr->gapSS;
         /*not applying extensions because I do not want
         `  to apply an extension penalty when an
         `  insertion changes to a deletion
         */

      #ifdef NOEXTEND
         delScoreSL = gapColSL + alnSetPtr->gapSS;
      #else
         delScoreSL = gapColSL + alnSetPtr->extendSS;
      #endif
   } /*Loop: score all query bases (rows)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^  - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   gapColSL -= alnSetPtr->extendSS;
      /*correct for being on the next row*/

   return gapColSL;
} /*reverseScore_hirschberg*/

/*-------------------------------------------------------\
| Fun04: recursive_hirschberg
|  - does the recurisve part of a hirschberg alignment
| Input
|  - refSeqStr:
|    o reference sequence to align
|  - refStartUL:
|    o first base to align in the reference (index 0)
|  - refLenUL:
|    o number of bases to align in the reference (index 1)
|  - qrySeqStr:
|    o query sequence to align
|  - qryStartUL:
|    o first base to align in the query (index 0)
|  - qryLenUL:
|    o number of bases to align in the query (index 1)
|  - forScoreRowL:
|    o row holding the forward (1st half of query) scores
|    o must be the size of the full length reference
|  - revScoreRowL:
|    o row holding the reverse (last half of query) scores
|    o must be the size of the full length reference
|  - refAlnST:
|    o holds the reference alignment and is a temporary
|      row for finding directions
|  - qryAlnST:
|    o holds the query alignment
|  - settings:
|    o pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Modifies:
|    o refAlnST to hold the reference alignment
|    o qryAlnST to hold the query alignment
\-------------------------------------------------------*/
void
recursive_hirschberg(
  signed char *refSeqStr,/*reference sequence*/
  signed long refStartSL,/*1st reference base to align*/
  signed long refLenSL,  /*number of ref bases to align*/

  signed char *qrySeqStr,/*query sequence*/
  signed long qryStartSL,/*1st query base to align*/
  signed long qryLenSL,  /*number query bases to align*/

  signed long *forScoreRowSL, /*holds final forward row*/
  signed long *revScoreRowSL, /*finding reverse scores*/
  
  signed char *refAlnStr, /*output reference alignment*/
  signed char *qryAlnStr, /*output query alignment*/

  struct alnSet *settings /*Settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '  - does the recursive part of a Hirschberg alignment
   '  o fun04 sec01:
   '    - Variable declerations
   '  o fun04 sec02:
   '    - Check if on a leaf (final part of alignment
   '  o fun04 sec03:
   '    - Get scores
   '  o fun04 sec04:
   '    - Find the midpoint
   '  o fun04 sec05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun04 Sec01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long forGapColSL = 0;
   signed long revGapColSL = 0;
   signed long midPointSL = 0;
   signed long slFor = 0;/*Loop iterator (forward score)*/

   signed long slGapOn = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^  - Check if on a leaf (final part of alignment
   ^  o fun04 sec02 sub01:
   ^    - Handle cases were I have just insertions
   ^  o fun04 sec02 sub02:
   ^    - Handle cases were I have just deletions
   ^  o fun04 sec02 sub03:
   ^    - Handle cases were I have to align last ref base
   ^  o fun04 sec02 sub04:
   ^  - Handle cases were I have to align last query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec02 Sub01:
   *  - Handle cases were I have just insertions
   \*****************************************************/

   if(refLenSL == 0)
   { /*If all remaing bases are query insertions*/
    for(
       slGapOn = qryStartSL;
       slGapOn < qryLenSL + qryStartSL;
       ++slGapOn
     ) qryAlnStr[slGapOn] = def_gap_hirschberg;

     return; /*Nothing else to do*/
   } /*If all remaing bases are query insertions*/

   /*****************************************************\
   * Fun04 Sec02 Sub02:
   *  - Handle cases were I have just deletions
   \*****************************************************/

   if(qryLenSL == 0)
   { /*If all remaing bases are query deletions*/
    for(
       slGapOn = refStartSL;
       slGapOn < refLenSL + refStartSL;
       ++slGapOn
    ) refAlnStr[slGapOn] = def_gap_hirschberg;

     return; /*Nothing else to do*/
   } /*If all remaing bases are query deletions*/

   /*****************************************************\
   * Fun04 Sec02 Sub03:
   *  - Handle cases were I have to align last ref base
   \*****************************************************/

   if(refLenSL == 1)
   { /*If I have to align the last reference base*/
     if(qryLenSL == 0)
     { /*If bases are aligned (one ref & one query)*/
        refAlnStr[refStartSL] = def_gap_hirschberg;
        return; 
     } /*If bases are aligned (one ref & one query)*/

     if(qryLenSL == 1)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnStr[qryStartSL] = def_nt_hirschberg;
        refAlnStr[refStartSL] = def_nt_hirschberg;
        return;
     } /*If bases are aligned (one ref & one query)*/

     alnOneNt_hirschberg(
       *(refSeqStr + refStartSL),/*ref base*/
       refStartSL,             /*Position of ref base*/
       qrySeqStr,              /*first base of query*/
       qryStartSL,             /*positoin of query*/
       qryLenSL,               /*Length of the query*/
       refAlnStr,              /*array to hold alignment*/
       qryAlnStr,              /*array to hold alignment*/
       settings                /*Has Scoring variables*/
     );

     return; /*This base is now aligned*/
   } /*If I have to align the last reference base*/

   /*****************************************************\
   * Fun04 Sec02 Sub04:
   *  - Handle cases were I have to align last query base
   \*****************************************************/

   if(qryLenSL == 1)
   { /*If I have to align the last query base*/
     if(refLenSL == 0)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnStr[qryStartSL] = def_gap_hirschberg;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     if(refLenSL == 1)
     { /*If bases are aligned (one ref & one query)*/
        qryAlnStr[qryStartSL] = def_nt_hirschberg;
        refAlnStr[refStartSL] = def_nt_hirschberg;
        return; /*Finished*/
     } /*If bases are aligned (one ref & one query)*/

     alnOneNt_hirschberg(
       *(qrySeqStr + qryStartSL),/*ref base*/
       qryStartSL,             /*Position of ref base*/
       refSeqStr,              /*first base of reference*/
       refStartSL,             /*positoin of query*/
       refLenSL,               /*Length of the query*/
       qryAlnStr,              /*array to hold alignment*/
       refAlnStr,              /*array to hold alignment*/
       settings                /*Has Scoring variables*/
     );

     return; /*Finshed aligning this query base*/
   } /*If I have to align the last query base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec03:
   ^  - Get scores
   ^  - for -DNOEXTEND, refAlnStr is ignored
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    forGapColSL = 
      forwardScore_hirschberg(
        refSeqStr,       /*Entire reference sequence*/
        refStartSL,      /*Starting base of ref target*/
        refLenSL,        /*length of ref target region*/
        qrySeqStr,       /*Query seq with coordinates*/
        qryStartSL,      /*Starting base of query target*/
        qryLenSL / 2,    /*Length of query target region*/
        forScoreRowSL,   /*Array of scores to fill*/
        refAlnStr,       /*direction row for gap extend*/
        settings         /*setttings to use*/
    ); /*Get the scores for the forward direction*/

    revGapColSL = 
      reverseScore_hirschberg(
        refSeqStr,     /*Entire reference sequence*/
        refStartSL,    /*Starting base of ref target*/
        refLenSL,      /*length of ref target region*/
        qrySeqStr,     /*Query seq with coordinates*/
        qryStartSL + (qryLenSL / 2),/*new query start*/
        qryLenSL - (qryLenSL / 2),  /*New query length*/
        revScoreRowSL,  /*Array of scores to fill*/
        refAlnStr,       /*direction row for gap extend*/
        settings       /* setttings to use*/
      ); /*Get the scores for the reverse direction*/
      /* I can get away with queryLen/2 here, because 
      `  queryLen is index 1 and the function takes in
      `  an lenth 1 argument
      `  This part could be made thread safe by using
      `  refAlnST for the forward score drection row and
      `  another array for the reverse direction row.
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun04 Sec04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   forScoreRowSL[refStartSL + refLenSL -1] += revGapColSL;
   midPointSL = refStartSL + refLenSL - 1;

   for(
      slFor = refStartSL;
      slFor < refStartSL + refLenSL - 1;
      ++slFor
   ){ /*Loop; add up all scores*/
     forScoreRowSL[slFor] += revScoreRowSL[slFor + 1];
       /*The reverse row is already reversed*/

     if(forScoreRowSL[slFor] > forScoreRowSL[midPointSL])
        midPointSL = slFor;
   } /*Loop; add up all scores*/

   forGapColSL += revScoreRowSL[refStartSL];

   if(forGapColSL > forScoreRowSL[midPointSL])
      midPointSL = 0;
   else
      midPointSL = midPointSL + 1 - refStartSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun04 Sec05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   recursive_hirschberg(
     refSeqStr,      /*Full reference sequence*/
     refStartSL,     /*Full reference sequence*/
     midPointSL,     /*Length of new reference sequence*/
     qrySeqStr,      /*Full query sequence*/
     qryStartSL,     /*Start of queyr target region*/
     qryLenSL / 2,   /*Length of query target region*/
     forScoreRowSL,  /*For scoring*/
     revScoreRowSL,  /*Has last line of scores*/
     refAlnStr,      /*direction row for gap extend*/
     qryAlnStr,      /*Holds the alignment codes*/
     settings        /*Settings for the alignment*/
   );

   recursive_hirschberg(
     refSeqStr,             /*Full reference sequence*/
     refStartSL + midPointSL,/*New reference start*/
     refLenSL - midPointSL, /*New reference end*/
     qrySeqStr,             /*Full query sequence*/
     qryStartSL + (qryLenSL / 2),/*New query start*/
     qryLenSL - (qryLenSL / 2),  /*New query length*/
     forScoreRowSL,         /*For scoring*/
     revScoreRowSL,         /*Has last line of scores*/
     refAlnStr,             /*Holds reference alingment*/
     qryAlnStr,             /*Holds query alingment*/
     settings               /*Settings for alignment*/
   );

   return;
} /*recursive_hirschberg*/

/*-------------------------------------------------------\
| Fun05: alnToSam_hirschberg
|  - converts an alignment from recursive_hirschberg to a
|    sam file entry
| Input:
|  - samSTPtr:
|    o samEntry struct pointer to add alignment to
|  - refST:
|    o pointer to seqST struct with the reference sequence
|      and the first base to align in the reference
|      sequence (refST->offsetUL; as index 0)
|  - qryST:
|    o pointer to seqST struct with the query sequence and
|      the first base to align in the query sequence
|      (qryST->offsetUL; as index 0)
|  - refAlnArySC:
|    o Has the reference alignment output by HirschbergFun
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - qryAlnArySC:
|    o Has the query alignment output by HirschbergFun
|    o if -DHIRSCHTWOBIT is two bit array, else char array
|  - settings:
|    o pointer to alnSet structure with the match matrix
|      used to determine if bases were a match or snp
| Output:
|  - Modifies:
|    o samSTPtr (everything) to have the alignment
|  - Returns:
|    o score of alignment for success
|    o def_memErr_hirschberg for memory errors
\-------------------------------------------------------*/
signed long
alnToSam_hirschberg(
  struct samEntry *samSTPtr, /*gets alignment*/
  struct seqST *refST,     /*reference sequence & offset*/
  struct seqST *qryST,     /*query sequence & offset*/
  signed char *refAlnArySC,/*has reference alignment*/
  signed char *qryAlnArySC,/*has query alignment*/
  struct alnSet *settings  /*has score/match matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC: hirschToAlnST
   '   - converts an alignment from recursive_hirschberg
   '     to a sam file entry
   '   o fun05 sec01:
   '     - variable declerations
   '   o fun05 sec02:
   '     - copy q-score and sequence entries
   '   o fun05 sec03:
   '     - build the cigar
   '   o fun05 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *refSeqStr = 0;
   signed char *qrySeqStr = 0;
   unsigned char matchBl = 0;
   unsigned int uiCig = 0;

   signed long scoreSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - copy q-score and sequence entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);


   samSTPtr->qryIdLenUC =
      cpWhite_ulCp(samSTPtr->qryIdStr, qryST->idStr);

   samSTPtr->refIdLenUC =
      cpWhite_ulCp(samSTPtr->refIdStr, refST->idStr);

   if(
      samSTPtr->seqSizeUI < (unsigned int) qryST->seqLenSL
   ){ /*If: need more memory for sequence*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr =
         malloc(
           (qryST->seqLenSL + 8) * sizeof(signed char)
         );
      if(! samSTPtr->seqStr)
         goto memErr_fun05_sec0x;
      samSTPtr->seqSizeUI = qryST->seqLenSL;
   } /*If: need more memory for sequence*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      qryST->seqStr,
      qryST->seqLenSL
   );
   samSTPtr->readLenUI = qryST->seqLenSL;
   indexToSeq_alnSet(samSTPtr->seqStr);


   if(! qryST->qStr)
      ;
   else if(! qryST->qStr[0])
      ;
   else if(qryST->qStr[0] && qryST->qStr[1] < 33)
      ;
   else
   { /*Else: have q-score entry to copy*/
      if(samSTPtr->qSizeUI < (unsigned int) qryST->qLenSL)
      { /*If: need more memory for q-score entry*/
         free(samSTPtr->qStr);
         samSTPtr->qStr =
            malloc(
              (qryST->qLenSL + 8) * sizeof(signed char)
            );
         if(! samSTPtr->qStr)
            goto memErr_fun05_sec0x;
         samSTPtr->qSizeUI = qryST->qLenSL;
      } /*If: need more memory for q-score entry*/

      cpQEntry_samEntry(samSTPtr, qryST->qStr, 1);
   } /*Else: have q-score entry to copy*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - build the cigar
   ^   o fun05 sec03 sub01:
   ^     - move to 1st base & get the first aligned base
   ^   o fun05 sec03 sub02:
   ^     - start loop and check if have insertions
   ^   o fun05 sec03 sub03:
   ^     - check for stop and snps/matches
   ^   o fun05 sec03 sub04:
   ^     - check if have deletions
   ^   o fun05 sec03 sub05:
   ^     - record the end and start of the alignmnent
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec03 Sub01:
   *   - move to first bases and get first aligned base
   \*****************************************************/

   refSeqStr = refST->seqStr + refST->offsetSL;
   qrySeqStr = qryST->seqStr + qryST->offsetSL;
   samSTPtr->cigTypeStr[0] = 0;

   if(qryST->offsetSL)
   { /*If: user wanted an offset*/
      samSTPtr->cigTypeStr[uiCig] = 'S';
      samSTPtr->cigArySI[uiCig] = qryST->offsetSL;

      samSTPtr->maskCntUI = qryST->offsetSL;
      samSTPtr->refStartUI = qryST->offsetSL;

      ++uiCig;
   } /*If: user wanted an offset*/

   else
      samSTPtr->refStartUI = 0;

   /*****************************************************\
   * Fun05 Sec03 Sub02:
   *   - start loop and check if have insertions
   \*****************************************************/

   while(
         *refAlnArySC != def_end_hirschberg
      || *qryAlnArySC != def_end_hirschberg
   ){ /*Loop: Add reference aligned bases to alnStruct*/
      if(*qryAlnArySC == def_gap_hirschberg)
      { /*If: I had an deltion*/
         ++samSTPtr->insCntUI;

         if(samSTPtr->cigTypeStr[uiCig] == 'I')
         { /*If: extending an insertion*/
            ++samSTPtr->cigArySI[uiCig];

            #ifdef NOEXTEND
               scoreSL += settings->gapSS;
            #else
               scoreSL += settings->extendSS;
            #endif
         } /*If: extending an insertion*/

         else
         { /*Else: need to add insertion entry*/
            if(samSTPtr->cigTypeStr[0])
               ++uiCig;
            if( checkCigSize_samEntry(samSTPtr, uiCig) )
               goto memErr_fun05_sec0x;

            samSTPtr->cigTypeStr[uiCig] = 'I';
            samSTPtr->cigArySI[uiCig] = 1;

            scoreSL += settings->gapSS;
         } /*Else: need to add insertion entry*/

         ++qryAlnArySC;
         ++qrySeqStr;
         continue;
      } /*If: I had an deltion*/

      /***************************************************\
      * Fun05 Sec03 Sub03:
      *   - check for stop and snps/matches
      \***************************************************/

      switch(*refAlnArySC)
      { /*Switch: Check the error type*/
         case 0: break;

         case def_nt_hirschberg:
         /*Case: snp or match*/
            matchBl =
               getMatch_alnSet(
                  *qrySeqStr,
                  *refSeqStr,
                  settings
               );
            scoreSL +=
               getScore_alnSet(
                  *qrySeqStr,
                  *refSeqStr,
                  settings
               );

            if(matchBl & def_ntEql_alnDefs)
            { /*If: I had a match*/
               ++samSTPtr->matchCntUI;

               if(samSTPtr->cigTypeStr[uiCig] == '=')
                  ++samSTPtr->cigArySI[uiCig];
               else
               { /*Else: need to add insertion entry*/
                  if(samSTPtr->cigTypeStr[0])
                     ++uiCig;
                  if(
                    checkCigSize_samEntry(samSTPtr, uiCig)
                  ) goto memErr_fun05_sec0x;

                  samSTPtr->cigTypeStr[uiCig] = '=';
                  samSTPtr->cigArySI[uiCig] = 1;
               } /*Else: need to add insertion entry*/
            } /*If: I had a match*/

            else
            { /*Else: I had an snp (mismatch)*/
               ++samSTPtr->snpCntUI;

               if(samSTPtr->cigTypeStr[uiCig] == 'X')
                  ++samSTPtr->cigArySI[uiCig];
               else
               { /*Else: need to add insertion entry*/
                  if(samSTPtr->cigTypeStr[0])
                     ++uiCig;
                  if(
                    checkCigSize_samEntry(samSTPtr, uiCig)
                  ) goto memErr_fun05_sec0x;

                  samSTPtr->cigTypeStr[uiCig] = 'X';
                  samSTPtr->cigArySI[uiCig] = 1;
               } /*Else: need to add insertion entry*/
            } /*Else: I had an snp (mismatch)*/

            ++refAlnArySC;
            ++qryAlnArySC;

            ++refSeqStr;
            ++qrySeqStr;

            break;
         /*Case: snp or match*/

         /************************************************\
         * Fun05 Sec03 Sub04:
         *   - check for deletions
         \************************************************/

         case def_gap_hirschberg:
         /*Case: Deletions*/
            ++samSTPtr->delCntUI;

            if(samSTPtr->cigTypeStr[uiCig] == 'D')
            { /*If: extending a deletion*/
               ++samSTPtr->cigArySI[uiCig];

               #ifdef NOEXTEND
                  scoreSL += settings->gapSS;
               #else
                  scoreSL += settings->extendSS;
               #endif
            } /*If: extending a deletion*/

            else
            { /*Else: need to add deletion entry*/
               if(samSTPtr->cigTypeStr[0])
                  ++uiCig;

               if(
                 checkCigSize_samEntry(samSTPtr, uiCig)
               ) goto memErr_fun05_sec0x;

               samSTPtr->cigTypeStr[uiCig] = 'D';
               samSTPtr->cigArySI[uiCig] = 1;
               scoreSL += settings->gapSS;
            } /*Else: need to add deletion entry*/
   
            ++refAlnArySC;
            ++refSeqStr;
            break;
         /*Case: Deletions*/
      } /*Switch: Check the error type*/
   } /*Loop: Add reference aligned bases to alnStruct*/

   /*****************************************************\
   * Fun05 Sec03 Sub05:
   *   - record the end and start of the alignmnent
   \*****************************************************/

   ++uiCig;
   samSTPtr->cigLenUI = uiCig;
   samSTPtr->cigTypeStr[uiCig] = 0;
   samSTPtr->cigArySI[uiCig] = 0;

   samSTPtr->alnReadLenUI = samSTPtr->matchCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->snpCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->delCntUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   /*****************************************************\
   * Fun05 Sec03 Sub06:
   *   - add edit distance and scores
   \*****************************************************/

   refSeqStr = samSTPtr->extraStr;

   /*edit distance*/
   *refSeqStr++ = 'N';
   *refSeqStr++ = 'M';
   *refSeqStr++ = ':';
   *refSeqStr++ = 'i';
   *refSeqStr++ = ':';

   refSeqStr +=
      numToStr(
         refSeqStr,
         samSTPtr->snpCntUI
            + samSTPtr->insCntUI
            + samSTPtr->delCntUI
      );
 

   /*assembler score*/
   *refSeqStr++ = '\t';
   *refSeqStr++ = 'A';
   *refSeqStr++ = 'S';
   *refSeqStr++ = ':';
   *refSeqStr++ = 'i';
   *refSeqStr++ = ':';

   refSeqStr +=
      signed_numToStr(refSeqStr, scoreSL);

   samSTPtr->extraLenUI = refSeqStr - samSTPtr->extraStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun05_sec0x;

   memErr_fun05_sec0x:;
      scoreSL = def_memErr_hirschberg;
      goto ret_fun05_sec0x;

   ret_fun05_sec0x:;
      return scoreSL;
} /*alnToSam_hirschberg*/

/*-------------------------------------------------------\
| Fun06: aln_hirschberg
|  - does hischberg alignment
| Input:
|  - refST:
|    o Pointer to seqST structure with the reference
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - qryST:
|    o Pointer to seqST structure with the query
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - samSTPtr:
|    o pointer to samEntry struct to store alignment in
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Returns:
|    o score of alignment for success
|    o def_memErr_hirschberg for memory errors
\-------------------------------------------------------*/
signed long
aln_hirschberg(
  struct seqST *refST, /*reference sequence to align*/
  struct seqST *qryST, /*qeury sequence to align*/
  struct samEntry *samSTPtr,/*gets alignment*/
  struct alnSet *settings   /*settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC: sam_hirschberg
   '   - hischberg alignment and returns samEntry struct
   '   o fun06 sec01:
   '     - variable declerations
   '   o fun06 sec02:
   '     - memory allocation (set up for Hirschberg)
   '   o fun06 sec03:
   '     - run hirschberg alignment
   '   o fun06 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun06 Sec01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   signed long lenRefSL =
      refST->endAlnSL - refST->offsetSL + 1;
   signed long lenQrySL =
      qryST->endAlnSL - qryST->offsetSL + 1;
     /*+ 1 to convert to index 1 (values are index 0)*/

   signed long *forScoreArySL = 0;
   signed long *revScoreArySL = 0;
   signed long scoreSL = 0;

   signed char *refAlnStr = 0;
   signed char *qryAlnStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun06 Sec02:
   ^    - memory allocation (set up for Hirschberg)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refAlnStr = calloc(lenRefSL + 1, sizeof(signed char));
   if(! refAlnStr)
      goto memErr_fun06_sec04;
   refAlnStr[lenRefSL] = def_end_hirschberg;

   qryAlnStr = calloc(lenQrySL + 1, sizeof(signed char));
   if(qryAlnStr == 0)
      goto memErr_fun06_sec04;
   qryAlnStr[lenQrySL] = def_end_hirschberg;

   forScoreArySL = malloc(sizeof(signed long) * lenRefSL);
   if(forScoreArySL == 0)
      goto memErr_fun06_sec04;

   revScoreArySL = malloc(sizeof(signed long) * lenRefSL);
   if(revScoreArySL == 0)
      goto memErr_fun06_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun06 Sec03:
   ^    - run hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*sending in offset values, because alignment array is
   ` sized to the alignment region
   */
   recursive_hirschberg(
     refST->seqStr + refST->offsetSL,
     0,                /*1st reference base to align*/
     lenRefSL,         /*length of ref region to align*/
     qryST->seqStr + qryST->offsetSL,
     0,                /*1st query base to align*/
     lenQrySL,         /*length of query target region*/
     forScoreArySL, /*for scoring*/
     revScoreArySL, /*for scoring*/
     refAlnStr,     /*holds the reference alignment*/
     qryAlnStr,     /*holds the query alignment*/
     settings       /*settings for the alignment*/
   );

   free(forScoreArySL);
   forScoreArySL = 0;

   free(revScoreArySL);
   revScoreArySL = 0;

   scoreSL =
     alnToSam_hirschberg(
        samSTPtr,
        refST,
        qryST,
        refAlnStr,
        qryAlnStr,
        settings
     ); /*convert alignment to a sam file entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun06 Sec04:
   ^    - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun06_sec04;

   memErr_fun06_sec04:;
      scoreSL = def_memErr_hirschberg;
      scoreSL = 1;
      goto ret_fun06_sec04;

   ret_fun06_sec04:;
      if(forScoreArySL)
         free(forScoreArySL);
      forScoreArySL = 0;

      if(revScoreArySL)
         free(revScoreArySL);
      revScoreArySL = 0;
 
      if(refAlnStr)
         free(refAlnStr);
      refAlnStr = 0;

      if(qryAlnStr)
         free(qryAlnStr);
      qryAlnStr = 0;

      return scoreSL;
} /*aln_hirschberg*/
