/*#########################################################
# Name: hirschberg
# Use:
#  - Holds functions for doing a hirschberg global
#    alignment
# Libraries:
#  - "alnStruct.h"
#  o "generalAlnFun.h"
#  o "alnMatrixStruct.h"
#  o "twoBitArrays.h"
#  o "scoresST.h"
#  o "seqStruct.h"
#  o "alnSetStruct.h"
#  o "alnSeqDefaults.h"
# C Standard Libraries:
#  o <stdlib.h>
#  o <stdint.h>
#  o <stdio.h>  // by alnSetStructure.h
#  o <string.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o fun-01 Hirschberg:
'    - Sets up for and calls the recursvie function to
'      run a Hirschberg alignment
'  o fun-02 HirschbergFun:
'    - Does the recursive part of a Hirschberg alignment
'  o fun-03 scoreForwardHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (forward direction)
'  o fun-04 scoreReverseHirsch:
'    - Does a single round of scoring for a hirschberg
'      alignment (reverse direction)
'  o fun-05 positionSingleRefBase:
'    - Align a single reference base to a query sequence
'  o fun-06 twoBitAlnToAlnST:
'    - Converts a two bit array with an alignment to an
'      alnStruct structure
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

#include "alnStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o A twoBitAry structure with the alignment.
|    o 0 For memory errors
\--------------------------------------------------------*/
struct alnStruct * Hirschberg(
  struct seqStruct *refST, /*Reference sequence to align*/
  struct seqStruct *qryST, /*Qeury sequence to align*/
    /* For refST and qryST, use seqStruct->offsetUL to set
    `  the starting point for the alignmnet and
    `  seqStruct->endAlnUL to set the ending point
    */
  struct alnSet *settings /*Settings to use for alignment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Hirschberg
   '  - Sets up for and calls the recursvie function to
   '    run a Hirschberg alignment
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-01 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-01 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the output alignment
\--------------------------------------------------------*/
void HirschbergFun(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,
    /*index 0: 1st bast to algin in reference*/
  unsigned long refLenUL,
    /*index 1 length of region to Align*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL,/*index 0 Starting query base*/
  unsigned long qryLenUL,/*index 1 Length of query region*/

  long *forwardScoreRowL,   /*Holds final forward row*/
  long *reverseScoreRowL,   /*For finding reverse scores*/
    /* both the forward and reverse scoring rows must be
    `  the size of the full length reference.
    */
  
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAlnST,/*Holds ref alignment*/
     struct twoBitAry *qryAlnST,/*Holds query alignment*/
     struct twoBitAry *dirRow,
       /*spare direction row (makes thread safe)*/
  #else
     char *refAlnST,   /*Holds output reference alignment*/
     char *qryAlnST,   /*Holds the output query alignment*/
     char *dirRow,/*direction row for thread safe scoring*/
  #endif
  struct alnSet *settings /*Settings to use for alignment*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: HirschbergFun
   '  - Does the recursive part of a Hirschberg alignment
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Check if on a leaf (final part of alignment
   '  o fun-02 sec-03:
   '    - Get scores
   '  o fun-02 sec-04:
   '    - Find the midpoint
   '  o fun-02 sec-05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch / Smith Waterman alignment
|    o dirRowSt to hold the last row of directions in a
|      Needleman Wunsch / Smith Waterman alignment
\--------------------------------------------------------*/
long scoreForwardHirsch(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,  /*index 0 starting ref base*/
  unsigned long refLenUL,    /*index 1 Length of target*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL, /*Index 0 Starting query base*/
  unsigned long qryLenUL,    /*index 1 length of target*/

  long *scoreRowPtrL,        /*Array of scores to fill*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *dirRowST,/*direction row*/
  #else
     char *dirRowST,      /*direction row, for gap extend*/
  #endif
  struct alnSet *settings    /*setttings to use*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreForwardHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (forward direction)
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-03 sec-03:
   '    - Score till on the last row
   '  o fun-03 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      backwards Needleman Wunsch /Smith Waterman alignment
|    o dirRowST to hold the last row of directions in a
|      backwards Needleman Wunsch /Smith Waterman alignment
\--------------------------------------------------------*/
long scoreReverseHirsch(
  char *refSeqCStr,          /*Reference sequence*/
  unsigned long refStartUL,  /*index 0 starting ref base*/
  unsigned long refLenUL,    /*index 1 Length of target*/

  char *qrySeqCStr,          /*Query sequence*/
  unsigned long qryStartUL, /*Index 0 Starting query base*/
  unsigned long qryLenUL,    /*index 1 length of target*/

  long *scoreRowPtrL,        /*Array of scores to fill*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *dirRowST,/*direction row*/
  #else
     char *dirRowST,      /*direction row, for gap extend*/
  #endif
  struct alnSet *settings    /*setttings to use*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreReverseHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (reverse direction)
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-04 sec-03:
   '    - Score till on the last row
   '  o fun-04 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the alignment for the single
|      base aligned to the sequence
\--------------------------------------------------------*/
void positionSingleBase(
  char baseC,             /*Single base to align to a seq*/
  unsigned long baseIndexUL,/*Index base is at*/
  char *seqCStr,            /*Sequence position on*/
  unsigned long startOfSeqUL,
    /*Index 0 of first base to align bascC to in seqCStr*/
  unsigned long lenSeqUL,
    /*Index 1; Length of the aligned region in seqCStr*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *baseCAlnST,
       /* Two bit alingment array for the sequence having
       `  baseC
       */
     struct twoBitAry *seqAlnST, 
       /* Two bit alignment array for the sequence alinging
       `  baseC to
       */
  #else
     char *baseCAlnST,
     char *seqAlnST, 
  #endif
  struct alnSet *settings        /*setttings to use*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single base to a sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: for error
|    o pointer to alnStruct with alignment
\--------------------------------------------------------*/
struct alnStruct * twoBitAlnToAlnST(
  struct seqStruct *refST,
   /*Has reference alignment start/end & reference length*/
  struct seqStruct *qryST,
   /*Has query alignment start/end and query length*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAlignment,
       /*Two bit array with the reference alignment*/
     struct twoBitAry *qryAlignment
       /*Two bit array with the query alignment*/
  #else
     char *refAlignment, /*has reference alignment*/
     char *qryAlignment /*has query alignment*/
  #endif
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: twoBitAlnToAlnST
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Set-01: maxGap (non-vector)
|  - Finds the best score and determines if the score was
|    an gap (del/ins) or snp/match.
|  - Each macro in the set prioritizes insertions,
|    deletions, and snps differently.
|  o set-01 macro-09: maxGapInsDelSnp
|  o set-01 macro-10: maxGapDelInsSnp
|  o set-01 macro-11: maxGapInsSnpDel
|  o set-01 macro-12: maxGapDelSnpIns
|  o set-01 macro-13: maxGapSnpInsDel
|  o set-01 macro-14: maxGapSnpDelIns
| Input:
|  - retMax;
|    o Holds the maximum value
|  - isGap:
|    o Holds if kept score was an snp (-1) or an gap (0)
|  - snp:
|    o Score for an snp
|  - ins:
|    o Score for an insertion
|  - del:
|    o Score for a deletion
| Output:
|  - Modifies:
|    o isGap to hold if best score was an snp or gap
|    o retMax to hold the best scores
\--------------------------------------------------------*/

/*--------------------------------------------------------\
| Set-01 Macro-09: maxGapDelInsSnp
\--------------------------------------------------------*/
#define maxGapInsDelSnp(retMax,isGap,snp,ins,del)\
{ /*maxGapInsDelSnp*/ \
   max((isGap), (ins), (del));         /*5 op*/ \
   max((retMax), (isGap), (snp));      /*5 op*/\
   (isGap) = -( (isGap) == (retMax) ); /*2 op*/\
      /*Logic:
      ` isGap == del
      `    1 if selected an insertion/deletion, else 0
      ` -(isGap == del)
      `    Changes 1 to -1 (111...1) (ins/del), else 0
      */\
} /*maxGapInsDelSnp*/

/*--------------------------------------------------------\
| Set-01 Macro-10: maxGapDelInsSnp
\--------------------------------------------------------*/
#define maxGapDelInsSnp(retMax,isGap,snp,ins,del)\
{ /*maxGapDelInsSnp*/\
   max((isGap), (del), (ins)); /*5 op*/ \
   max((del), (isGap), (snp)); /*5 op*/\
   (isGap) = -( (isGap) == (del) ); /*2 op*/\
      /*Logic:
      ` isGap == del
      `    1 if selected an insertion/deletion, else 0
      ` -(isGap == del)
      `    Changes 1 to -1 (111...1) (ins/del), else 0
      */\
} /*maxGapDelInsSnp*/

/*--------------------------------------------------------\
| Set-01 Macro-11: maxGapInsSnpDel
\--------------------------------------------------------*/
#define maxGapInsSnpDel(retMax,isGap,snp,ins,del)\
{ /*maxGapInsSnpDel*/\
   max((retMax), (snp), (del))               /*5 op*/\
   (isGap) = ((snp) < (retMax));             /*1 op*/\
   \
   max((retMax), (ins), (retMax));           /*5 op*/ \
   (isGap) = -((ins == (retMax)) | (isGap)); /*3op*/\
      /*Logic
      ` isGap == 1
      `    is if a deletion is prefered (del > snp), else 0
      ` ins == del
      `    is 1 if a insertion is prefered, else 0
      ` (ins == del) | isGap
      `    is 1 if an insertion or deletion was prefered
      ` -((in == del) | isGap)
      `    is -1 (111...1) if an insertion or deletion
      `    was preffered, else 0
      */\
} /*maxGapInsSnpDel*/

/*--------------------------------------------------------\
| Set-01 Macro-12: maxGapDelSnpIns
\--------------------------------------------------------*/
#define maxGapDelSnpIns(retMax,isGap,snp,ins,del)\
{ /*maxGapDelSnpIns*/\
   max((retMax), (snp), (ins));              /*5 op*/\
   (isGap) = ((snp) < (retMax));             /*1 op*/\
   max((retMax), (del), (retMax));           /*5 op*/\
   (isGap) = -((del == (retMax)) | (isGap)); /*3 op*/\
      /*Logic
      ` isGap
      `    1 if a deletion is prefered (snp <= del), else 0
      ` ins == del
      `    1 if a insertion is prefered, else 0
      ` (ins == del) | isGap
      `    is 1 if an insertion or deletion was prefered
      ` -((in == del) | isGap)
      `    is -1 (111...1) if an insertion or deletion
      `    was preffered, else 0
      */\
} /*maxGapDelSnpIns*/

/*--------------------------------------------------------\
| Set-01 Macro-13: maxGapSnpInsDel
\--------------------------------------------------------*/
#define maxGapSnpInsDel(retMax,isGap,snp,ins,del)\
{\
   max((retMax), (snp), (ins));   /*5 op*/ \
   max((retMax), (retMax), (del));   /*5 op*/\
   (isGap) = -((snp) < (retMax)); /*2 op*/\
      /*Logic:
      ` snp < del
      `    is 1 if deletion selected, else 0
      ` -(snp < del)
      `    is -1 (111...1) if deletion/insertion, else 0
      */\
} /*maxGapSnpInsDel*/

/*--------------------------------------------------------\
| Set-01 Macro-14: maxGapSnpDelIns
\--------------------------------------------------------*/
#define maxGapSnpDelIns(retMax,isGap,snp,ins,del)\
{\
   max((retMax), (snp), (del)); /*5 op*/ \
   max((retMax), (retMax), (ins)); /*5 op*/\
   (isGap) = -((snp) < (retMax));    /*2 op*/\
      /*Logic:
      ` snp < del
      `    is 1 if deletion selected, else 0
      ` -(snp < del)
      `    is -1 (111...1) if deletion/insertion, else 0
      */\
} /*maxGapSnpDelIns*/

/*--------------------------------------------------------\
| Name: maxGapInsDel (Macro-15:)
| Use:
|  - Chooses which maxGapX function to call to maximize and
|    get the gaps.
| Input:
|  o snp:
|    - Is the score for the snp
|  o ins:
|    -  for the insertion
|  o del:
|    -  for an deltion. Will hold the max score
|  o isGap:
|    - Is set to -1 if gap chosen, 0 for snp/match
|  o preflag:
|    - Flag having preferance for direction
|    - options: defInsDelSnp, defDelInsSnp, defInsSnpDel,
|               defDelSnpIns, defSnpInsDel, defSnpDelIns
| Output:
|  - Modifies::
|    o del to hold the maximum value
|    o isGap to if was a gap or not
\--------------------------------------------------------*/
/*Check if harcoding a preferared direction*/
#ifdef INSDELSNP
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapInsDelSnp((retMax),(isGap),(snp),(ins),(del))}

#elif defined DELINSSNP
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapDelInsSnp((retMax),(isGap),(snp),(ins),(del))}

#elif defined INSSNPDEL
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapInsSnpDel((retMax),(isGap),(snp),(ins),(del))}

#elif defined DELSNPINS
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapDelSnpIns((retMax),(isGap),(snp),(ins),(del))}

#elif defined SNPINSDEL
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapSnpInsDel((retMax),(isGap),(snp),(ins),(del))}

#elif defined SNPDELINS
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
      {maxGapSnpdelIns((retMax),(isGap),(snp),(ins),(del))}

#else
   #define maxGapScore(retMax,isGap,snp,ins,del,prefFlag)\
   {\
    switch((prefFlag))\
    { /*Switch; get an snp/match priority*/\
     case defSnpInsDel:\
      maxGapSnpInsDel((retMax),(isGap),(snp),(ins),(del));\
      break;\
     case defSnpDelIns:\
      maxGapSnpDelIns((retMax),(isGap),(snp),(ins),(del));\
      break;\
     case defInsSnpDel:\
      maxGapInsSnpDel((retMax),(isGap),(snp),(ins),(del));\
      break;\
     case defDelSnpIns:\
      maxGapDelSnpIns((retMax),(isGap),(snp),(ins),(del));\
      break;\
     case defInsDelSnp:\
      maxGapInsDelSnp((retMax),(isGap),(snp),(ins),(del));\
      break;\
     case defDelInsSnp:\
      maxGapDelInsSnp((retMax),(isGap),(snp),(ins),(del));\
      break;\
    } /*Switch; get an snp/match priority*/\
   }
#endif /*Check if harcoding a preferared direction*/

#endif
