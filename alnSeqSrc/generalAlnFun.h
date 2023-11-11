/*########################################################
# Name generalAlnFun
# Use:
#  o Holds general functions used in my Needleman Wunsch
#    or Waterman Smith alignment.
# Libraries:
#   - "alnSetStruct.h"
#   o "alignmentSettings.h"
# C Standard libraries:
#   o <stdlib.h>
#   o <stdint.h>
#   o <stdio.h>  // Used by alnSetStructures.h
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
' o fun-01 twoBitMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for two bit arrays.
'   - Function is inlied and 2nd to last function in header
' o fun-02 charMaxScore:
'   - Picks the best score and direction for the current
'     base pairs being compared in an alignment. This
'     function is set up for charters
'   - Function is inlined and is last function in header
'  o fun-03 alnMaxScore:
'    - Picks the best score for the current base pairs
'      being compared in an alignment.
'  o fun-04 checkIfBasesMatch:
'    - Are two bases are same? (includes anonymous bases)
'  o macro-01 indelScore:
'    - Calculates the score for an indel
'  o macro-02 max:
'    - Find the maximum value (branchless)
'  o macro-03 ifMax:
'    - Set a value (ret) to a value based on which value
'      is greater.
'  o macro-04 snpInsDel:
'    - Selects the max score and direction for max score.
'    - snps and then insertions when max scores are equal.
'  o macro-05 snpDelIns:
'    - Selects the max score and direction for max score.
'    - perfers snps and then deletions when max scores are
'      equal.
'  o macro-06 insSnpDel:
'    - Selects the max score and direction for max score.
'    - Perfers insertions and then snps when max scores
'      are equal.
'  o macro-07 insDelSnp:
'    - Selects the max score and direction for max score.
'    - Perfers indels and then deletions when max scores
'      are equal.
'  o macro-08: delSnpIns
'    - Selects the max score and direction for max score.
'    - perfers deletions and then snps when max scores are
'      equal.
'  o macro-09: delInsSnp
'    - Selects the max score and direction for max score.
'    - Perfers insertions and then deletions when max
'      scores are equal.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef GENERALLALNFUN_H
#define GENERALLALNFUN_H

#include "alnSetStruct.h"

/*--------------------------------------------------------\
| Output:
|  - Returns
|    o 1: if bases were a match
|    o 0 if bases do not mach
\--------------------------------------------------------*/
char checkIfBasesMatch(
    char *queryBaseC,// Query base to compare to reference
    char *refBaseC   // Reference base to compare to query
); /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: Sec-1 Sub-1: checkIfBasesMatch
   '  - Are two bases are same? (includes anonymous bases)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Macro-01: indelScore
| Use:
|  - Calculates the score for an indel
| Input:
|   - retScore:
|     o Holds the score for the indel
|   - dirC:
|     o Direction of the last cell
|   - oldScore:
|     o The previous score to add to the gap penalty
|   - settings:
|     o Has the gap opeing and extension penalties.
| Output:
|  - Modifies:
|    o retScore to hold the score for an indel
\--------------------------------------------------------*/
#define indelScore(retScore, dirC, oldScore, settings){\
   (retScore) = 0 - ((dirC) == defMvSnp); \
   (retScore) = \
        (oldScore) \
      + ((settings)->gapOpenS & (retScore)) \
      + ((settings)->gapExtendS & (~(retScore))); \
    /* Branchless, apply a gap opening and extension
    `   penalty. This is faster then my old switch
    `   statment method.
    ` tmpDirUL = 0 - (direction != diagnol)
    `   direction != diagnol
    `     Becomes 1 if the move results in a gap
    `   0 - direction != diagnol
    `     Becomes 0 - 1 if have a gap. This results in
    `       all bits being set to 1.
    `     Becomes 0 - 0 if have no gaps. This is 0.
    ` score + (gapOpen & (!dir)) + (gapExtend & dir)
    `   gapOpen & (~dir)
    `     ~dir gives zero if all bits are set to one
    `       (gap), but all bits set to one if 0.
    `     If ~dir is zero, gapOpen = 0, else applies the
    `       gap opening penalty.
    `   gapExtend & (dir)
    `     Becomes zero when there are no gaps (dir = 0).
    */ \
}

/*--------------------------------------------------------\
| Macro-02: max
| Use:
|  - Find the maximum value (branchless)
| Input:
|  o ret:
|    - Value to hold the maximum value
|  o x:
|    - First value to compare, ret is set to x if x >= y
|  o y:
|    - second value to compare, ret is set to y if y > x
| Output:
|  - Sets:
|    - Sets ret to x if x is greater than or equal to y
|    - Sets ret to y if y is greather than x
| From:
|  - https://graphics.stanford.edu/~seander/bithacks.html#IntegerMinOrMax
\--------------------------------------------------------*/
#define max(ret, x, y){\
   (ret) = (x) ^ (((x) ^ (y)) & (-((x) < (y)))); \
   /*Logic:
   `  x < y:
   `    if x < y I get 0 (x > y); else I get 1
   `  -(x < y):
   `    If x < y I get -1 (111...)
   `    If x >= y I get 0
   `  x ^ ((x ^ y) & 111...): (only when x < y)
   `    This gives me x
   `  x ^ (x ^ y) & 0: (only when x >= y)
   `    This gives me y
   */ \
}

/*--------------------------------------------------------\
| Macro-03: ifMax
| Use:
|  - Set a value (ret) to a value based on which value
|    is greater.
| Input:
|  o ret:
|    - This will hold the return value
|  o x:
|    - First value to compare, (if x >= y)
|  o y:
|    - second value to compare, (if y > x)
|  o xRet:
|    - Value to set ret of x is >= y
|  o yRet:
|    - Value to set ret of y is > x
| Output:
|  - Sets:
|    - ret to xRet if x is greater than or equal to y
|    - ret to yRet if y is greater than x
\--------------------------------------------------------*/
#define ifMax(ret, x, y, xRet, yRet){\
   (ret) = (xRet) ^ (((xRet) ^ (yRet)) & (-((x) < (y)))); \
   /*This follows the same logic as max(ret, x, y), except
   ` instead of setting ret to the highest value, I set
   ` ret to xRet if x is >= to y or yRet if y > x.
   */ \
}

/*--------------------------------------------------------\
| Set-01: Maximize score/direction (non-vector)
|  - Selects the max score and direction selected for the
|    max score.
|  o set-01 macro-04 a: insDelSnp (branched)
|  o set-01 macro-05 a: delInsSnp (branched)
|  o set-01 macro-06 a: insSnpDel (branched)
|  o set-01 macro-07 a: delSnpIns (branched)
|  o set-01 macro-08 a: snpInsDel (branched)
|  o set-01 macro-09 a: snpDelIns (branched)
|
|  o Set-01 Macro-04 b: insDelSnp (no braching)
|  o set-01 macro-05 b: delInsSnp (no branching)
|  o set-01 macro-06 b: insSnpDel (no branching)
|  o set-01 macro-07 b: delSnpIns (no branching)
|  o set-01 macro-08 b: snpInsDel (no branching)
|  o set-01 macro-09 b: snpDelIns (no branching)
| Variations:
|  - The default is branchless (no ifs, set b), but a
|    branched variation (if statements, set a) can be used
|    with -DBRANCHED 
| Input:
|  o maxSc:
|    - This will hold the max score
|  o maxDir
|    - This will hold the direction of the max score
|  o insSc:
|    - Score for having an insertion at this position
|  o snpSc
|    - Score for having an SNP/match at this position
|  o delSc:
|    - Score for having an deletion at this position
| Output:
|  - Sets:
|    o Sets maxDir to the direction of the max score
|    - Sets maxSc to the max score
\--------------------------------------------------------*/

/*branched or banchless score/direction maximize*/
#ifdef BRANCHED
   /*-----------------------------------------------------\
   | Set-01 Macro-04 a: insDelSnp
   \-----------------------------------------------------*/
   #define insDelSnp(maxSc, maxDir, snpSc, insSc, delSc){\
      if((delSc) > (insSc))\
         {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
      else {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
      if((snpSc) > (maxSc))\
         {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-05 a: delInsSnp
   \-----------------------------------------------------*/
   #define delInsSnp(maxSc, maxDir, snpSc, insSc, delSc){\
      if((insSc) > (delSc))\
         {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
      else {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
      if((snpSc) > (maxSc))\
         {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-06 a: insSnpDel
   \-----------------------------------------------------*/
   #define insSnpDel(maxSc, maxDir, snpSc, insSc, delSc){\
      if((snpSc) > (insSc))\
         {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
      else {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
      if((delSc) > (maxSc))\
         {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-07 a: delSnpIns
   \-----------------------------------------------------*/
   #define delSnpIns(maxSc, maxDir, snpSc, insSc, delSc){\
      if((snpSc) > (delSc))\
         {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
      else {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
      if((insSc) > (maxSc))\
         {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-08 a: snpInsDel
   \-----------------------------------------------------*/
   #define snpInsDel(maxSc, maxDir, snpSc, insSc, delSc){\
      if((insSc) > (snpSc))\
         {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
      else {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
      if((delSc) > (maxSc))\
         {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-09 a: snpDelIns
   \-----------------------------------------------------*/
   #define snpDelIns(maxSc, maxDir, snpSc, insSc, delSc){\
      if((delSc) > (snpSc))\
         {(maxSc) = (delSc); (maxDir) = (defMvDel);}\
      else {(maxSc) = (snpSc); (maxDir) = (defMvSnp);}\
      if((insSc) > (maxSc))\
         {(maxSc) = (insSc); (maxDir) = (defMvIns);}\
   }

#else
   /*-----------------------------------------------------\
   | Set-01 Macro-04 b: insDelSnp
   \-----------------------------------------------------*/
   #define insDelSnp(maxSc, maxDir, snpSc, insSc, delSc){\
      max((maxSc), (insSc), (delSc));            /*5 Op*/\
      (maxDir) = (snpSc) > (maxSc);              /*1 Op*/\
      max((maxSc), (maxSc), (snpSc));            /*5 Op*/\
      (maxDir) +=                                /*4 Op*/\
         ( ( (insSc) == (maxSc) ) | (maxDir) ) + defMvDel;\
      /*Logic:
      ` maxDir = snp > max
      `   1 if snp was selected, else 0 (ins/del)
      ` defMvDel = 1
      ` (ins == max) | maxDir
      `   1 if insertion or snp was selected, else 0 (del)
      ` For an snp 1 += (0 | 1) + 1 = 3
      ` For an ins 0 += (1 | 0) + 1 = 2
      ` For an del 0 += (0 | 0) + 1 = 1
      */\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-05 b: delInsSnp
   \-----------------------------------------------------*/
   #define delInsSnp(maxSc, maxDir, snpSc, insSc, delSc){\
      max((maxSc), (delSc), (insSc));             /*5 Op*/\
      maxDir = (snpSc) > (maxSc);                 /*1 Op*/\
      max((maxSc), (maxSc), (snpSc));             /*5 Op*/\
      (maxDir) +=                                 /*3 Op*/\
         ( (delSc) != (maxSc) ) + defMvDel;\
      /*Logic:
      ` maxDir = snp > max
      `   maxDir is 1 if snp was selected, else 0
      ` defMvDel = 1
      ` del != maxSc
      `   1 if a deletion (ins/snp) was selected
      ` For an snp 1 += (1) + 1 = 3
      ` For an ins 0 += (1) + 1 = 2
      ` For an del 0 += 0 + 1 = 1
      */\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-06 a: insSnpDel
   \-----------------------------------------------------*/
   #define insSnpDel(maxSc, maxDir, snpSc, insSc, delSc){\
      max(maxSc, (insSc), (snpSc));               /*5 Op*/\
      (maxDir) = (delSc) <= (maxSc);              /*1 Op*/\
        /*1 if kept score is not a deletion, else 0*/\
      max((maxSc), (maxSc), (delSc));             /*5 Op*/\
      (maxDir) +=                                 /*4 Op*/\
         ( ( (snpSc) > (insSc) ) & maxDir ) + defMvDel;\
      /*Logic:
      ` maxDir del <= max
      `   1 if an insertion or snp, else 0 (deletion)
      ` defMvDel = 1
      ` (snp > ins) & maxDir
      `   1 if an snp is selected, else 0 (ins/del)
      ` For an snp 1 += (1 & 1) + 1 = 3
      ` For an ins 1 += (0 & 1) + 1 = 2
      ` For an del 0 += (0/1 & 0) + 1 = 1
      */\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-07 b: delSnpIns
   \-----------------------------------------------------*/
   #define delSnpIns(maxSc, maxDir, snpSc, insSc, delSc){\
      max((maxSc), (snpSc), (insSc));             /*5 Op*/\
      (maxDir) = (maxSc) <= (delSc);              /*1 Op*/\
      max((maxSc), (delSc), (maxSc));             /*5 Op*/\
      (maxDir) =                                  /*4 Op*/\
          defMvSnp\
        - ( ( (insSc) > (snpSc) ) | (maxDir) )\
        - (maxDir);\
      /*Logic:
      `  maxDir is 1 if a deletion is selected (max <= del)
      `  defMvSnp = 3
      `  (ins > snp) | (maxDir)
      `     1 if insertion or deletion selected, else 0
      `  -maxDir
      `     Is 1 if deletion was selected
      ` For an del I get 3 - (0/1 | 1) - 1 = 1
      ` For an ins I get 3 - (1 | 0) - 0 = 2
      ` For an snp I get 3 - (0 | 0) - 0 = 3
      */\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-08 b: snpInsDel
   \-----------------------------------------------------*/
   #define snpInsDel(maxSc, maxDir, snpSc, insSc, delSc){\
      max((maxSc), (snpSc), (insSc));             /*5 Op*/\
      (maxDir) = ((delSc) <= (maxSc));            /*1 Op*/\
        /*1 if deletion not choosen, else 0*/\
      max((maxSc), (maxSc), (delSc));             /*5 Op*/\
      (maxDir) +=                                 /*3 Op*/\
         ( (snpSc) == (maxSc) ) + defMvDel;\
      /*Logic
      ` maxDir = delSc <= maxSc
      `   1 if an insertion or snp was selected, else 0 (del)
      ` snp == max
      `   1 if the snp was selected.
      ` For an snp 1 += (1) + 1 = 3
      ` For an ins 1 += (0) + 1 = 2
      ` For an del 0 += (0) + 1 = 1
      */\
   }

   /*-----------------------------------------------------\
   | Set-01 Macro-09 b: snpDelIns
   \-----------------------------------------------------*/
   #define snpDelIns(maxSc, maxDir, snpSc, insSc, delSc){\
      max((maxSc), (delSc), (insSc));             /*5 Op*/\
      maxDir = (snpSc) < (maxSc);                 /*1 Op*/\
      max((maxSc), (snpSc), (maxSc));             /*5 Op*/\
      (maxDir) =                                  /*4 Op*/\
           defMvSnp /*Always 3*/\
         - (maxDir) \
         - ( ( (delSc) == (maxSc) ) & (maxDir) );\
        /*Logic
        ` maxDir = snp < max
        `   1 if deletion/insertion selected, else 0 (snp)
        ` defMvSnp is always 3
        ` (del == max) & maxDir
        `   1 if a deletin was selected, else 0 (snp/ins)
        ` For an snp (over ins/snp)   3 - 0 - (0/1 & 0) = 3
        ` For an ins (ignore if eqal) 3 - 1 - (0 & 0) = 2
        ` For an del (over ins)       3 - 1 - (1 & 1) = 1
        */\
   }
#endif /*branched or banchless score/diretion maximize*/

/*--------------------------------------------------------\
| Name: charMaxScore (Macro-10:)
| Use:
|  - This finds the max score based on the users
|    preferences for direction
| Input:
|  - maxScore:
|    o This will hold the maximum score
|  - maxDir:
|    o This will hold the direction of the max score
|  - snp:
|    o Score for making an snp/match move
|  - ins:
|    o Score for making an insertion move
|  - del:
|    o Score for making an deletion move
|   - pref:
|     o Direction combination that is prefered when scores
|       are equal. This can be hardcoded in.
| Output:
|  - Modifies
|    o maxScore to hold the maxium score
|    o maxDir to hold the direction of the max score
\--------------------------------------------------------*/
#ifdef INSDELSNP
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {insDelSnp((maxScore),(maxDir),(snp),(ins),(del))}

#elif defined DELINSSNP
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {delInsSnp((maxScore),(maxDir),(snp),(ins),(del))}

#elif defined INSSNPDEL 
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {insSnpDel((maxScore),(maxDir),(snp),(ins),(del))}

#elif defined DELSNPINS
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {delSnpIns((maxScore),(maxDir),(snp),(ins),(del))}

#elif defined SNPINSDEL
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {snpInsDel((maxScore),(maxDir),(snp),(ins),(del))}

#elif defined SNPDELINS
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
     {snpDelIns((maxScore),(maxDir),(snp),(ins),(del))}

#else
   #define charMaxScore(maxScore,maxDir,snp,ins,del,pref)\
{\
   switch(pref)\
   { /*Switch; get an snp/match priority*/\
      case defInsDelSnp:\
        insDelSnp((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
      case defDelInsSnp:\
        delInsSnp((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
      case defInsSnpDel:\
        insSnpDel((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
      case defDelSnpIns:\
        delSnpIns((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
      case defSnpInsDel:\
        snpInsDel((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
      case defSnpDelIns:\
        snpDelIns((maxScore),(maxDir),(snp),(ins),(del));\
        break;\
   } /*Switch; get an snp/match priority*/\
} /*charMaxScore*/

#endif /*For charMaxScore variations*/

/*--------------------------------------------------------\
| Name: alnMaxScore (Macro-11:)
| Use: 
|  - Picks the best score for the current base pairs
|     being compared in an alignment.
| Variations:
|  - This can be complied as a switch statement or a single
|    direction can be hardcoded in
| Input:
|   - maxScore:
|     o Holds the maximum score
|   - snp:
|     o Score for an snp/match
|   - ins:
|     o Score for an insertion
|   - del:
|     o Score for an deletion
|   - prefFlag:
|     o Direction combination that is prefered when scores
|       are equal. This can be hardcoded in.
| Output:
|   - Modifies
|     o maxScore to hold the maximum score
\--------------------------------------------------------*/
 /*Check how maximizes scores*/
#if defined INSDELSNP
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (ins), (del));\
      max((maxScore), (maxScore), (snp));\
   }

#elif defined DELINSSNP
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (del), (ins));\
      max((maxScore), (maxScore), (snp));\
   }

#elif defined INSSNPDEL
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (ins), (snp));\
      max((maxScore), (maxScore), (del));\
   }

#elif defined DELSNPINS
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (del), (snp));\
      max((maxScore), (maxScore), (ins));\
   }

#elif defined SNPINSDEL
   #define alnMaxScoreore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (snpScore), (insScore));\
      max((maxScore), (maxScore), (delScore));\
   }

#elif defined SNPDELINS
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      max((maxScore), (snp), (del));\
      max((maxScore), (maxScore), (ins));\
   }

#else
   #define alnMaxScore(maxScore,snp,ins,del,prefFlag){\
      switch(prefFlag)\
      { /*Switch; get an snp/match priority*/\
         case defSnpInsDel:\
           max((maxScore), (snp), (ins));\
           max((maxScore), (maxScore), (del));\
           break;\
         case defSnpDelIns:\
           max((maxScore), (snp), (del));\
           max((maxScore), (maxScore), (ins));\
           break;\
         case defInsSnpDel:\
           max((maxScore), (ins), (snp));\
           max((maxScore), (maxScore), (del));\
           break;\
         case defInsDelSnp:\
           max((maxScore), (ins), (del));\
           max((maxScore), (maxScore), (snp));\
           break;\
         case defDelSnpIns:\
           max((maxScore), (del), (snp));\
           max((maxScore), (maxScore), (ins));\
           break;\
         case defDelInsSnp:\
           max((maxScore), (del), (ins));\
           max((maxScore), (maxScore), (snp));\
           break;\
      } /*Switch; get an snp/match priority*/\
   } /*alnMaxScoreore*/
#endif /*Check how maximizes scores*/

#endif

