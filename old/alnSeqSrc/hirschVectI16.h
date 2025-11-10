/*#########################################################
# Name: hirschVectI16
# Use:
#  - Holds functions and macros for running a hirschberg
#    with 16 bit vector support
# Libraries:
#  - "vectorWrap.h"
#  o Compiler flag determined
#    i) or "vectWrapAVX512.h" 
#    ii) or "vectWrapAVX2.h"
#    iii) or "vectWrapSSE.h"
#    iv) or "vectWrapNEON.h"
#    v) or "vectWrapNEON64.h"
# C Standard Libraries:
#   o <limits.h>
#   o <stdint.h>
#   o <stdlib.h>
#   o Compiler flaged determined
#     i) or <immintrin.h>
#     ii) or <arm_neon.h>
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' hirschVectFun SOH: Start Of Header
'  - header and macros for the 16 bit hirschberg
'  o fun-07 vectAddGaps:
'    - Adds gaps into a array of characters. Uses vectors
'      when adding a large number of gaps(>= vector length)
'  o set-01 vI16MaxGap: [macros 01 to 06] (16 bit vector)
'    - Set of macros for finding the maximum value and if
'      the maximum was a gap (-1) or snp/match (0).
'  o macro-07 vI16FindMaxGap:
'    - Determines macro to call for finding the best score
'      and then stores the best scores in delVectI16 and
'      the gap/no-gap in gapVectI16. This functions
'      determines which vI16MaxGapXXX function should be
'      called.
'    - 4 vectors + 4 operations or 5 vectors + 5 operations
'  o set-02 maxGap: [macros-09 to 10] (non-vector)
'    - Set of macros for finding the maximum value and if
'      the maximum was a gap (-1) or snp/match (0).
'  o macro-15 maxGapInsDel:
'    - Uses macros from set-02 maxGap to maximize the value
'      and decide if a gap was preffered.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef HIRSCHVECTI16_H /*Guard for header*/
#define HIRSCHVECTI16_H /*Guard for header*/

#include "vectorWrap.h"

/*--------------------------------------------------------\
| Name: vectAddGaps
| Use:
|  - Adds gaps into a array of characters. Uses vectors
|    when large enough.
| Input:
|  - alnAryC:
|    o Array to add gaps to
|  - lenGapL:
|    o number of gaps to add
|  - gapC:
|    o Symbol of gap to add
| Output:
|  - Modifies:
|    o alnAryC to have lenGapL gaps (gapC)
\--------------------------------------------------------*/
void vectAddGaps(
   char *alnAryC, /*Has alignment array to add gaps to*/
   unsigned long lenGapL, /*Size of gap*/
   char gapC,     /*Symbol for a gap*/
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: vectAddGaps
   '  - Adds gaps to an array (by vectors if large enough)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*--------------------------------------------------------\
| Set-01: vI16MaxGap (16 bit vector)
|  - Finds the best score and determines if the score was
|    an gap (del/ins) or snp/match.
|  - Each macro in the set prioritizes insertions,
|    deletions, and snps differently.
|  o set-01 macro-01: vI16MaxGapInsDelSnp
|  o set-01 macro-02: vI16MaxGapDelInsSnp
|  o set-01 macro-03: vI16MaxGapInsSnpDel
|  o set-01 macro-04: vI16MaxGapDelSnpIns
|  o set-01 macro-05: vI16MaxGapSnpInsDel
|  o set-01 macro-06: vI16MaxGapSnpDelIns
| Input:
|  - insVectI16:
|    o Vector with scores for insertions
|  - delVectI16:
|    o Vector with scores for deletions. This will hold
|      the best scores.
|  - snpVectI16:
|    o Vector with scores for snps/matches
|  - gapVectI16:
|    o Vector told hold if kept score was an snp (0)
|      or a gap (-1)
| Output:
|  - Modifies:
|    o gapVectI16 to hold if best score was an snp or gap
|    o delVectI16 to hold the best scores
| Note:
|  - I would recomend finding the next rounds snp scores
|    before calling this function. That way the last rounds
|    gap scores can be used.
\--------------------------------------------------------*/

/*--------------------------------------------------------\
| Set-01 Macro-01: vI16MaxGapInsDelSnp
\--------------------------------------------------------*/
#define vI16MaxGapInsDelSnp(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   mask16 snpGapMask16; /*For comparisons*/\
   \
   (gapVectI16) = \
      max_vectI16_retVectI16((insVectI16), (delVectI16));\
      /*See if an insertion is better than a deletion*/\
   \
   (delVectI16) = \
      max_vectI16_retVectI16(\
         (gapVectI16),\
         (snpVectI16)\
   ); /*See if a snp is preffered*/\
   \
   /*See if this move was a gap (del is gap extend/open)*/\
   snpGapMask16 = \
      cmpeq_vectI16_retMask16((delVectI16), (gapVectI16));\
      /*See if snp was favored over an gap*/\
   \
   /*Convert the masks to vectors*/\
   gapVectI16 = mask16_retVectI16(snpGapMask16);\
} /*vI16MaxGapInsDelSnp*/

/*--------------------------------------------------------\
| Set-01 Macro-02: vI16MaxGapDelInsSnp
\--------------------------------------------------------*/
#define vI16MaxGapDelInsSnp(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   mask16 snpGapMask16; /*For comparisons*/\
   \
   (gapVectI16) = \
      max_vectI16_retVectI16((delVectI16), (insVectI16));\
      /*See if an insertion is better than a deletion*/\
   \
   (delVectI16) = \
      max_vectI16_retVectI16(\
         (gapVectI16),\
         (snpVectI16)\
   ); /*See if a snp is preffered*/\
   \
   /*See if this move was a gap (del is gap extend/open)*/\
   snpGapMask16 = \
      cmpeq_vectI16_retMask16((delVectI16), (gapVectI16));\
      /*See if snp was favored over an gap*/\
   \
   /*Convert the masks to vectors*/\
   gapVectI16 = mask16_retVectI16(snpGapMask16);\
} /*vI16MaxGapDelInsSnp*/

/*--------------------------------------------------------\
| Set-01 Macro-03: vI16MaxGapInsSnpDel
\--------------------------------------------------------*/
#define vI16MaxGapInsSnpDel(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   vectI16 tmpVI16;\
   \
   (delVectI16) = \
      max_vectI16_retVectI16((snpVectI16), (delVectI16));\
   \
   /*Checking to see if a deletion was selected*/\
   gapVectI16 = \
      mask16_retVectI16(\
         cmpgt_vectI16_retMask16(\
            (delVectI16),\
            (snpVectI16)\
      )\
   ); /*This cast comand is only compiled for AVX512*/\
   \
   /*Finish finding the maxium*/\
   (delVectI16) = \
      max_vectI16_retVectI16((insVectI16), (delVectI16));\
   \
   tmpVI16 =\
      mask16_retVectI16(\
         cmpeq_vectI16_retMask16(\
            (insVectI16), /*Score for insertions*/\
            (delVectI16)  /*Has maximum value*/\
         ) /*is -1 if insertion selected, else 0*/\
   );  /*cast comand is only compiled for AVX512*/\
   \
   /*Finish determining if a gap was selected*/\
   gapVectI16=or_vectI16_retVectI16((gapVectI16),tmpVI16);\
} /*vI16MaxGapInsSnpDel*/

/*--------------------------------------------------------\
| Set-01 Macro-04: vI16MaxGapDelSnpIns
\--------------------------------------------------------*/
#define vI16MaxGapDelSnpIns(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   vectI16 tmpVI16;\
   \
   tmpVI16 = \
      max_vectI16_retVectI16((delVectI16), (snpVectI16));\
   \
   /*Check to see if a deletion was selected*/\
   gapVectI16 = \
      mask16_retVectI16(\
         cmpeq_vectI16_retMask16(\
            (delVectI16),\
            tmpVI16\
      )\
   ); /*This cast comand is only compiled for AVX512*/\
   \
   /*Finish finding the maximum*/\
   (delVectI16) = \
      max_vectI16_retVectI16(tmpVI16, (insVectI16));\
   \
   /*Finsh finding if a gap was selected*/\
   tmpVI16 =\
      mask16_retVectI16(\
         cmpgt_vectI16_retMask16(\
            (insVectI16), /*Score for insertions*/\
            (delVectI16)  /*Has maximum value*/\
         ) /*is -1 if insertion selected, else 0*/\
   );  /*cast comand is only compiled for AVX512*/\
   \
   /*Finish determining if a gap was selected*/\
   gapVectI16=or_vectI16_retVectI16((gapVectI16),tmpVI16);\
} /*vI16MaxGapDelSnpIns*/

/*--------------------------------------------------------\
| Set-01 Macro-05: vI16MaxGapSnpInsDel
\--------------------------------------------------------*/
#define vI16MaxGapSnpInsDel(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   mask16 snpGapMask16; /*For comparisons*/\
   \
   (gapVectI16) = \
      max_vectI16_retVectI16((insVectI16), (delVectI16));\
      /*See if an insertion is better than a deletion*/\
   \
   (delVectI16) = \
      max_vectI16_retVectI16(\
         (snpVectI16),\
         (gapVectI16)\
   ); /*See if a snp is preffered*/\
   \
   /*See if this move was a gap (del is gap extend/open)*/\
   snpGapMask16 = \
      cmplt_vectI16_retMask16((snpVectI16), (gapVectI16));\
      /*See if snp was favored over an gap*/\
   \
   /*Convert the masks to vectors*/\
   gapVectI16 = mask16_retVectI16(snpGapMask16);\
} /*vI16MaxGapSnpInsDel*/

/*--------------------------------------------------------\
| Set-01 Macro-06: vI16MaxGapSnpDelIns
\--------------------------------------------------------*/
#define vI16MaxGapSnpDelIns(\
   insVectI16, /*insertion move scores (used as temp)*/\
   delVectI16, /*has deletions, ends with best scores*/\
   snpVectI16, /*snp move scores (set to insVectI16)*/\
   gapVectI16  /*Will hold if best score was snp or gap*/\
) {\
   mask16 snpGapMask16; /*For comparisons*/\
   \
   (gapVectI16) = \
      max_vectI16_retVectI16((delVectI16), (insVectI16));\
      /*See if an insertion is better than a deletion*/\
   \
   (delVectI16) = \
      max_vectI16_retVectI16(\
         (snpVectI16),\
         (gapVectI16)\
   ); /*See if a snp is preffered*/\
   \
   /*See if this move was a gap (del is gap extend/open)*/\
   snpGapMask16 = \
      cmplt_vectI16_retMask16((snpVectI16), (gapVectI16));\
      /*See if snp was favored over an gap*/\
   \
   /*Convert the masks to vectors*/\
   gapVectI16 = mask16_retVectI16(snpGapMask16);\
} /*vI16MaxGapSnpDelIns*/

/*--------------------------------------------------------\
| Name: vI16FindMaxGap (Macro-07:)
| Use:
|  - Determines macro to call for finding the best score
|    and then stores the best scores in delVectI16 and
|    the gap/no-gap in gapVectI16. This functions
|    determines which vI16MaxGapXXX function should be
|    called.
| Input:
|  - insVectI16:
|    o Vector with scores for insertions
|  - delVectI16:
|    o Vector with scores for deletions. This will hold
|      the best scores.
|  - snpVectI16:
|    o Vector with scores for snps/matches
|  - gapVectI16:
|    o Vector told hold if if kept score was an snp (-1)
|      or a gap (0)
|  - priorityFlagC:
|    o Tells which direction (ins/del/snp) to prioritize
|      when scores are equal.
| Output:
|  - Modifies:
|    o gapVectI16 to hold if best score was an snp or gap
|    o delVectI16 to hold the best scores
| Note:
|  - insSnpDel and delSnpIns require one extra vector and
|    one extra operation.
|  - I would recomend finding the next rounds snp scores
|    before calling this function. That way the last rounds
|    gap scores can be used.
\--------------------------------------------------------*/
#if defined INSDELSNP
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapInsDelSnp(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#elif defined DELINSSNP
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapDelInsSnp(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#elif defined INSSNPDEL
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapInsSnpDel(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#elif defined DELSNPINS
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapDelSnpIns(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#elif defined SNPINSDEL
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapSnpInsDel(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#elif defined SNPDELINS
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   )(\
      vI16MaxGapSnpDelIns(\
         (insVectI16),\
         (delVectI16),\
         (snpVectI16),\
         (gapVectI16)\
      );\
   )

#else
   #define vI16FindMaxGap(\
      insVectI16,\
      delVectI16,\
      snpVectI16,\
      gapVectI16,\
      priorityFlagC\
   ){\
     switch(alnSetST->bestDirC)\
     { /*Switch; get an snp/match priority*/\
       case defInsDelSnp:\
           vI16MaxGapInsDelSnp(\
              (insVectI16),\
              (delVectI16),\
              (snpVectI16),\
              (gapVectI16)\
           );\
           break;\
       case defDelInsSnp:\
          vI16MaxGapDelInsSnp(\
             (insVectI16),\
             (delVectI16),\
             (snpVectI16),\
             (gapVectI16)\
          );\
          break;\
       case defInsSnpDel:\
          vI16MaxGapInsSnpDel(\
             (insVectI16),\
             (delVectI16),\
             (snpVectI16),\
             (gapVectI16)\
          );\
          break;\
       case defDelSnpIns:\
          vI16MaxGapDelSnpIns(\
             (insVectI16),\
             (delVectI16),\
             (snpVectI16),\
             (gapVectI16)\
          );\
          break;\
       case defSnpInsDel:\
          vI16MaxGapSnpInsDel(\
             (insVectI16),\
             (delVectI16),\
             (snpVectI16),\
             (gapVectI16)\
          );\
          break;\
       case defSnpDelIns:\
          vI16MaxGapSnpDelIns(\
             (insVectI16),\
             (delVectI16),\
             (snpVectI16),\
             (gapVectI16)\
          );\
          break;\
     } /*Switch; get an snp/match priority*/\
   } /*vI16FindMaxGap*/
#endif /*Macro-07 (max best scores)*/


#endif /*guard for header*/


