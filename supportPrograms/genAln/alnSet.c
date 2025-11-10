/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' alnSet SOF: Start Of File
'  - holds settings structures and supporting functions
'    for my pairwise aligners
'  o header:
'    - included libraries
'  o .h st01 alnSet:
'     o Holds settings for my alignment program
'  o .h fun01 setScore_alnSet:
'    - Sets the score for a base pair (reference/query)
'  o .h fun02 setMatch_alnSet:
'    - Sets if two bases are a match or mismtach
'  o fun03 freeStack_alnSet:
'    o Frees variables inside alnSet
'  o fun04 freeHeap_alnSet:
'    o Frees an alnSet structure (and sets to 0)
'  o .h fun05 getScore_alnSet:
'    - Get the score for a pair of bases from an alignment
'  o .h fun06 getMatch_alnSet:
'    - Check if two bases were a match or mismatch
'  o fun07 readScoreFile_alnSet
'     - Reads in a file of scores for a scoring matrix
'  o fun08 readMatchFile_alnSet:
'    - Reads in a file of matches
'  o fun09 seqToIndex_alnSet:
'    - Converts a sequence to a look up table index
'  o fun10 indexToSeq_alnSet:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
'  o fun11: revCmpIndex_alnSet
'    - reverse complement sequence in lookup index format
'      (from seqToIndex_alnSet)
'  o fun12: changeGap_alnSet
'    - changes the gap penalties in an alnSet structure
'  o fun13: maxScore_alnSet
'    - finds maximum score possible for a sequence
'  o fun14 init_alnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun15: pDefMatchMatrix_alnSet
'    - print out the default match matrix
'  o license:
'    - Licensing for this code (public domain / mit)
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

#include "alnSet.h"

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/ulCp.h"
#include "../genLib/fileFun.h"

/*no .c files*/
#include "../genLib/endLine.h"
#include "alnDefs.h"

/*-------------------------------------------------------\
| Fun03: freeStack_alnSet
|  - Does a stack free of an alnSet structer
| Input:
|  - alnSetSTPtr:
|    o alnSetSTPtr to free internal variables in
| Output:
|  - Free:
|    o Nothing, there are no heap allocated variables.
|      Is here in case I decide to have heap allocated
|      variables on day.
\-------------------------------------------------------*/
void
freeStack_alnSet(
   struct alnSet *alnSetSTPtr
){
   init_alnSet(alnSetSTPtr); /*avoids error messages*/
} /*freeStack_alnSet*/

/*-------------------------------------------------------\
| Fun04: freeHeap_alnSet
|  - Frees and alnSet (alignment settings) structure
| Input:
|  - alnSetSTPtr:
|    o Pionter to alnSetST to free
| Output:
|  - Free:
|    o alnSetSTPtr
|  - Set:
|    o alnSetSTPtr to 0 (NULL)
\-------------------------------------------------------*/
void
freeHeap_alnSet(
   struct alnSet *alnSetSTPtr
){
   freeStack_alnSet(alnSetSTPtr);
   free(alnSetSTPtr);
} /*freeHeap_alnSet*/

/*-------------------------------------------------------\
| Fun07: readScoreFile_alnSet
|  - Reads in a file of scores for a scoring matrix
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the score matrix to modify
|  - scoreFILE:
|    o File to get scores from
| Output:
|  - Modifies:
|    o Score matrix in alngSetST to hold the scores from
|      the file (scoreFILE)
|    o scoreFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readScoreFile_alnSet(
    struct alnSet *alnSetSTPtr, /*score matrix to change*/
    void *inFILE         /*File scoring matrix scores*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC: readScoreFile_alnSet
   '  o fun07 sec01:
   '    - Variable declerations & set up
   '  o fun07 sec02:
   '    - Blank the scoring matrix
   '  o fun07 sec03:
   '    - Read in line and check if comment
   '  o fun07 sec04:
   '    - Convert score & add to matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun07 1024
   signed char buffStr[def_lenBuff_fun07];
   signed char *tmpStr = 0;
   signed short scoreSS = 0;

   unsigned char colUC = 0;
   unsigned char rowUC = 0;
   signed long lenSL = 0;

   buffStr[def_lenBuff_fun07 - 1] = '\0';
   buffStr[def_lenBuff_fun07 - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetSTPtr->scoreMatrixSS[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(1)
   { /*Loop: read in scores*/
      lenSL =
         getLine_fileFun(inFILE, buffStr, 1024, &lenSL);
      if(! lenSL)
         break;
       
      if(buffStr[0] == '/' && buffStr[1] == '/')
      { /*On a comment, move onto the next line*/
          --lenSL;

          while(buffStr[lenSL] == '\0')
          { /*Loop: read in more buffer*/
             lenSL =
                getLine_fileFun(
                   inFILE,
                   buffStr,
                   1024,
                   &lenSL
                );
             --lenSL;
          } /*Loop: read in more buffer*/

          continue;
      } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun07 Sec04:
       ^  - Convert score & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffStr[0] == '\r')
           continue;                        /*Blank line*/
       else if(buffStr[0] == '\n')
           continue;                        /*Blank line*/
       else if(buffStr[0] == '\0')
           continue;                        /*Blank line*/

       if(buffStr[0] < 64 && buffStr[2] < 64)
           return ftell(inFILE);  /*Invalid character*/
       
       tmpStr +=
           strToSS_base10str(
              &buffStr[4],
              &scoreSS
           );

       setScore_alnSet(
         buffStr[0],
         buffStr[2],
         scoreSS,
         alnSetSTPtr
       ); /*Add the score to the matrix*/

       if(tmpStr == &buffStr[3])
           return ftell(inFILE);         /*No score*/

       scoreSS = endLine_ulCp(buffStr);

       while(
              buffStr[lenSL] != '\n'
           && buffStr[lenSL] != '\r'
       ){ /*Loop: get next line*/
          lenSL =
             getLine_fileFun(inFILE,buffStr,1024,&lenSL);
           if(! lenSL)
              break;
           --lenSL;
       }  /*Loop: get next line*/
   } /*Loop: read in scores*/

   return 0;
} /*readScoreFile_alnSet*/

/*-------------------------------------------------------\
| Fun08: readMatchFile_alnSet
|  - Reads in a file of matches
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the match matrix to modify
|  - matchFILE:
|    o File to get matchs from
| Output:
|  - Modifies:
|    o Match matrix in alngSetST to hold the matchs from
|      the file (matchFILE)
|    o matchFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readMatchFile_alnSet(
    struct alnSet *alnSetSTPtr,
    void *matchFILE  /*File of matchs for scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun08 TOC: readMatchFile_alnSet
   '  o fun08 sec01:
   '    - Variable declerations & set up
   '  o fun08 sec02:
   '    - Blank the match matrix
   '  o fun08 sec03:
   '    - Read in line and check if comment
   '  o fun08 sec04:
   '    - Add match/snp (mismatch) to match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun08 1024
   signed char buffStr[def_lenBuff_fun08];
   signed long lenSL = 0;

   unsigned char colUC = 0;
   unsigned char rowUC = 0;

   buffStr[def_lenBuff_fun08 - 1] = '\0';
   buffStr[def_lenBuff_fun08 - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetSTPtr->matchMatrixSC[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(1)
   { /*While I have matchs to read in*/
      lenSL =
         getLine_fileFun(
            matchFILE,
            buffStr,
            1024,
            &lenSL
          );
       if(! lenSL)
          goto ret_fun08_sec04;
       --lenSL;
       
       if(buffStr[0] == '/' && buffStr[1] == '/')
       { /*On a comment, move onto the next line*/

           while(! buffStr[lenSL])
           { /*While have more buffer to read in*/
              lenSL =
                 getLine_fileFun(
                    matchFILE,
                    buffStr,
                    1024,
                    &lenSL
                  );
               if(! lenSL)
                  goto ret_fun08_sec04;
              --lenSL;
           } /*While have more buffer to read in*/

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun08 Sec04:
       ^  - Convert match & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffStr[0] == '\r')
           continue;                        /*Blank line*/
       else if(buffStr[0] == '\n')
           continue;                        /*Blank line*/
       else if(buffStr[0] == '\0')
           continue;                        /*Blank line*/

       if(buffStr[4] == '1')
          ; /*match*/
       else if( (buffStr[4] & ~32) == 'Y' )
          buffStr[4] = '1'; /*match*/
       else if( (buffStr[4] & ~32) == 'T' )
          buffStr[4] = '1'; /*match*/
       else if(buffStr[4] == '0')
          ; /*not match*/
       else if( (buffStr[4] & ~32) == 'N')
          buffStr[4] = '0'; /*mismatch*/
       else if( (buffStr[4] & ~32) == 'F')
          buffStr[4] = '0'; /*mismatch*/
       else if( (buffStr[4] & ~32) == 'N')
          buffStr[4] = def_anonymous_alnDefs;
          /*anonymosus base*/
       else
           return ftell((FILE *) matchFILE);
           /*This error means I do not know if match/snp*/

       setMatch_alnSet(
         buffStr[0],
         buffStr[2],
         buffStr[4],
         alnSetSTPtr
       ); /*Add the match to the matrix*/

       while(
              buffStr[lenSL] != '\n'
           && buffStr[lenSL] != '\r'
       ){ /*Loop: get next line*/
          lenSL =
             getLine_fileFun(
                matchFILE,
                buffStr,
                1024,
                &lenSL
              );
           if(! lenSL)
              goto ret_fun08_sec04;
           --lenSL;
       }  /*Loop: get next line*/
   } /*While I have matchs to read in*/

   ret_fun08_sec04:;
      return 0;
} /*readMatchFile_alnSet*/

/*-------------------------------------------------------\
| Fun09: seqToIndex_alnSet
|  - Converts a sequence to a look up table index
|    (table is in alnSetStruct.c/h)
| Input:
|  - seqStr:
|    o pointer c-string with the sequence to convert
| Output:
|  - Modifies:
|    o seqST->seqStr to have look up table indexs (1-27, 
|      with null as 0) instead of bases
\-------------------------------------------------------*/
void
seqToIndex_alnSet(
   signed char *seqStr
){
   while(*seqStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *seqStr &= def_ntToCode_alnSet;
      ++seqStr;
   } /*Loop: convert bases to lookup table values*/
} /*seqToIndex_alnSet*/

/*-------------------------------------------------------\
| Fun10: indexToSeq_alnSet
|  - Converts a sequence of lookup indexs back into
|    uppercase characters (a-z)
| Input:
|  - seqStructPtr:
|    o Pointer to sequence structer with converte sequence
|      to deconvert (lookup index to base)
| Output:
|  - Modifies:
|    o seqST->seqStr to have bases instead of look up
|      table indexs
\-------------------------------------------------------*/
void
indexToSeq_alnSet(
   signed char *seqStr
){
   while(*seqStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *seqStr |= 64;
      ++seqStr;
   } /*Loop: convert bases to lookup table values*/
} /*indexToSeq_alnSet*/

/*-------------------------------------------------------\
| Fun11: revCmpIndex_alnSet
|   - reverse complement sequence in lookup index format
|     (from seqToIndex_alnSet)
| Input:
|   - seqStr:
|     o c-string with look up index sequence to
|       reverse complement
|   - syncStr:
|     o c-string to keep in sync with seqStr
|     o use 0/null for no sequence
|   - lenSeqSL:
|     o length of sequence to convert (index 1)
| Output:
|   - Modifies:
|     o seqStr to be reverse complemented
\-------------------------------------------------------*/
void
revCmpIndex_alnSet(
   signed char *seqStr,
   signed char *syncStr,
   unsigned long lenSeqUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - reverse complement sequence in lookup index
   '     format (from seqToIndex_alnSet)
   '   o fun11 sec01:
   '     - variable declarations and reverse complement
   '   o fun11 sec02:
   '     - reverse complement final base
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declarations and reverse complement
   ^   o fun11 sec01 sub01:
   ^     - variable declaration
   ^   o fun11 sec01 sub02:
   ^     - complement start base + start loop
   ^   o fun11 sec01 sub03:
   ^     - complement end base
   ^   o fun11 sec01 sub04:
   ^     - swap start and end bases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /**************************************************\
   * Fun11 Sec01 Sub01:
   *   - variable declaration
   \**************************************************/

   signed char *endStr = seqStr + lenSeqUL - 1;
   signed char *syncEndStr = 0;

   if(! lenSeqUL)
      return; /*nothing to do*/

   if(syncStr)
      syncEndStr = syncStr + lenSeqUL - 1;
  

   /**************************************************\
   * Fun11 Sec01 Sub02:
   *   - complement start base + start loop
   \**************************************************/

    while(endStr > seqStr)
    { /*Loop: reverse complement*/

       switch(*seqStr)
       { /*switch, reverse complement*/

           case ('A' & def_ntToCode_alnSet):
              *seqStr = 'T' & def_ntToCode_alnSet;
               break;
           case ('C' & def_ntToCode_alnSet):
              *seqStr = 'G' & def_ntToCode_alnSet;
               break;
           case ('G' & def_ntToCode_alnSet):
              *seqStr = 'C' & def_ntToCode_alnSet;
               break;
           case ('T' & def_ntToCode_alnSet):
              *seqStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('U' & def_ntToCode_alnSet):
              *seqStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('W' & def_ntToCode_alnSet):
              *seqStr = 'W' & def_ntToCode_alnSet;
               break;
           case ('S' & def_ntToCode_alnSet):
              *seqStr = 'S' & def_ntToCode_alnSet;
               break;
           case ('M' & def_ntToCode_alnSet):
              *seqStr = 'K' & def_ntToCode_alnSet;
               break;
           case ('K' & def_ntToCode_alnSet):
              *seqStr = 'M' & def_ntToCode_alnSet;
               break;
           case ('R' & def_ntToCode_alnSet):
              *seqStr = 'Y' & def_ntToCode_alnSet;
               break;
           case ('Y' & def_ntToCode_alnSet):
              *seqStr = 'R' & def_ntToCode_alnSet;
               break;
           case ('B' & def_ntToCode_alnSet):
              *seqStr = 'V' & def_ntToCode_alnSet;
               break;
           case ('D' & def_ntToCode_alnSet):
              *seqStr = 'H' & def_ntToCode_alnSet;
               break;
           case ('H' & def_ntToCode_alnSet):
              *seqStr = 'D' & def_ntToCode_alnSet;
               break;
           case ('V' & def_ntToCode_alnSet):
              *seqStr = 'B' & def_ntToCode_alnSet;
               break;
           case ('N' & def_ntToCode_alnSet):
              *seqStr = 'N' & def_ntToCode_alnSet;
               break;
           default:
              *seqStr = 0;

       } /*switch, reverse complement*/

       /*************************************************\
       * Fun11 Sec01 Sub03:
       *   - complement end base
       \*************************************************/

       switch(*endStr)
       { /*switch, reverse complement*/

           case ('A' & def_ntToCode_alnSet):
              *endStr = 'T' & def_ntToCode_alnSet;
              break;
           case ('C' & def_ntToCode_alnSet):
              *endStr = 'G' & def_ntToCode_alnSet;
               break;
           case ('G' & def_ntToCode_alnSet):
              *endStr = 'C' & def_ntToCode_alnSet;
               break;
           case ('T' & def_ntToCode_alnSet):
              *endStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('U' & def_ntToCode_alnSet):
              *endStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('W' & def_ntToCode_alnSet):
              *endStr = 'W' & def_ntToCode_alnSet;
               break;
           case ('S' & def_ntToCode_alnSet):
              *endStr = 'S' & def_ntToCode_alnSet;
               break;
           case ('M' & def_ntToCode_alnSet):
              *endStr = 'K' & def_ntToCode_alnSet;
               break;
           case ('K' & def_ntToCode_alnSet):
              *endStr = 'M' & def_ntToCode_alnSet;
               break;
           case ('R' & def_ntToCode_alnSet):
              *endStr = 'Y' & def_ntToCode_alnSet;
               break;
           case ('Y' & def_ntToCode_alnSet):
              *endStr = 'R' & def_ntToCode_alnSet;
               break;
           case ('B' & def_ntToCode_alnSet):
              *endStr = 'V' & def_ntToCode_alnSet;
               break;
           case ('D' & def_ntToCode_alnSet):
              *endStr = 'H' & def_ntToCode_alnSet;
               break;
           case ('H' & def_ntToCode_alnSet):
              *endStr = 'D' & def_ntToCode_alnSet;
               break;
           case ('V' & def_ntToCode_alnSet):
              *endStr = 'B' & def_ntToCode_alnSet;
               break;
           case ('N' & def_ntToCode_alnSet):
              *endStr = 'N' & def_ntToCode_alnSet;
               break;
           default:
              *endStr = 0;

       } /*switch, reverse complement*/

       /*************************************************\
       * Fun11 Sec01 Sub04:
       *   - swap start and end bases
       \*************************************************/

       *seqStr ^= *endStr;
       *endStr ^= *seqStr;
       *seqStr ^= *endStr;
       
       ++seqStr;
       --endStr;

       if(syncStr)
       { /*If: keeping c-string in sync*/
          *syncStr ^= *syncEndStr;
          *syncEndStr ^= *syncStr;
          *syncStr ^= *syncEndStr;
          
          ++syncStr;
          --syncEndStr;
       } /*If: keeping c-string in sync*/
    } /*Loop: reverse complement*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun11 Sec02:
    ^   - reverse complement final base
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if ended on same base; if so complement base*/
    if(endStr == seqStr)
    { /*If: ending on same base (not swaped)*/

       switch(*seqStr & (~32)) /*Make sure upper case*/
       { /*switch: last base*/

           case ('A' & def_ntToCode_alnSet):
               *seqStr = 'T' & def_ntToCode_alnSet;
               break;
           case ('C' & def_ntToCode_alnSet):
               *seqStr = 'G' & def_ntToCode_alnSet;
               break;
           case ('G' & def_ntToCode_alnSet):
               *seqStr = 'C' & def_ntToCode_alnSet;
               break;
           case ('T' & def_ntToCode_alnSet):
               *seqStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('U' & def_ntToCode_alnSet):
               *seqStr = 'A' & def_ntToCode_alnSet;
               break;
           case ('W' & def_ntToCode_alnSet):
               *seqStr = 'W' & def_ntToCode_alnSet;
               break;
           case ('S' & def_ntToCode_alnSet):
               *seqStr = 'S' & def_ntToCode_alnSet;
               break;
           case ('M' & def_ntToCode_alnSet):
               *seqStr = 'K' & def_ntToCode_alnSet;
               break;
           case ('K' & def_ntToCode_alnSet):
               *seqStr = 'M' & def_ntToCode_alnSet;
               break;
           case ('R' & def_ntToCode_alnSet):
               *seqStr = 'Y' & def_ntToCode_alnSet;
               break;
           case ('Y' & def_ntToCode_alnSet):
               *seqStr = 'R' & def_ntToCode_alnSet;
               break;
           case ('B' & def_ntToCode_alnSet):
               *seqStr = 'V' & def_ntToCode_alnSet;
               break;
           case ('D' & def_ntToCode_alnSet):
               *seqStr = 'H' & def_ntToCode_alnSet;
               break;
           case ('H' & def_ntToCode_alnSet):
               *seqStr = 'D' & def_ntToCode_alnSet;
               break;
           case ('V' & def_ntToCode_alnSet):
               *seqStr = 'B' & def_ntToCode_alnSet;
               break;
           case ('N' & def_ntToCode_alnSet):
               *seqStr = 'N' & def_ntToCode_alnSet;
               break;
           default:
              *seqStr = 0;

       } /*switch: last base*/

    } /*If: ending on same base (not swaped)*/
} /*revCmpIndex_alnSet*/

/*-------------------------------------------------------\
| Fun12: changeGap_alnSet
|   - changes the gap penalties in an alnSet structure
| Input:
|   - alnSetSTPtr:
|     o pointer to alnSet structure to update gap scores
|   - gapSS:
|     o new gap opening penalty
|   - extendSS:
|     o new gap extension penalty
| Output:
|   - Modifies:
|     o delArySS to have gap extension/opening penalties
|     o insArySS to have gap extension/opening penalties
|     o gapSS in alnSetSTPtr to have gap opening penalty
|     o extendSS to have gap extension penalty
\-------------------------------------------------------*/
void
changeGap_alnSet(
   struct alnSet *alnSetSTPtr,
   signed short gapSS,
   signed short extendSS
){

   alnSetSTPtr->gapSS = gapSS;
   alnSetSTPtr->extendSS = extendSS;
   
   alnSetSTPtr->delArySS[def_mvStop_alnDefs] = gapSS;
   alnSetSTPtr->delArySS[def_mvDel_alnDefs] = extendSS;
   alnSetSTPtr->delArySS[def_mvIns_alnDefs] = gapSS;
   alnSetSTPtr->delArySS[def_mvSnp_alnDefs] = gapSS;

   /*insetion gap penalty matrix*/
   alnSetSTPtr->insArySS[def_mvStop_alnDefs] = gapSS;
   alnSetSTPtr->insArySS[def_mvDel_alnDefs] = gapSS;
   alnSetSTPtr->insArySS[def_mvIns_alnDefs] = extendSS;
   alnSetSTPtr->insArySS[def_mvSnp_alnDefs] = gapSS;
} /*changeGap_alnSet*/

/*-------------------------------------------------------\
| Fun13: maxScore_alnSet
|  - finds maximum score possible for a sequence
| Input:
|  - seqStr:
|    o c-string with sequence to find max score for
|    o needs to be converted to index with
|      indexToSeq_alnSet (prepared for alignment)
|  - startUL:
|    o index of first base to check (index 0)
|  - lenSeqUL:
|    o length of the sequence to check
|  - alnSetPtr:
|    o pointer to an alnSet (settings) structure to
|      initialize
| Output:
|  o Returns:
|    - maximum possible score for alignment
\-------------------------------------------------------*/
signed long
maxScore_alnSet(
   signed char *seqStr,
   unsigned long startUL,
   unsigned long lenSeqUL,
   struct alnSet *alnSetPtr
){
   signed long scoreSL = 0;

   while(startUL < lenSeqUL)
   { /*Loop: align sequences*/
      scoreSL +=
         getScore_alnSet(
            seqStr[startUL],
            seqStr[startUL],
            alnSetPtr
         );

      ++startUL;
   } /*Loop: align sequences*/

   return scoreSL;
} /*maxScore_alnSet*/

/*-------------------------------------------------------\
| Fun14: init_alnSet
|  - Set values in altSet (alingment settings) structure
|    to default values
| Input:
|  - alnSetSTPtr:
|    o poineter to an alnSet (settings) structure to
|      initialize
| Output:
|  o Modifies:
|    - alnSetST to have default alignment settings values
\-------------------------------------------------------*/
void
init_alnSet(
    struct alnSet *alnSetST /*Has settings to initialize*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC: init_alnSet
   '  - Set values in altSet (alingment settings)
   '    structure to defaults
   '  o fun14 sec01:
   '    - Set non-matrix variables
   '  o fun14 sec02:
   '    - Initialize scoring matrix
   '  o fun14 sec03:
   '    - Initialize match matrix
   '  o fun14 sec04:
   '    - set up scoring matrix for nucleotides
   '  o fun14 sec05:
   '    - set up matching matrix for nucleotides
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^  - Set non-matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Variables for my for loop*/
   unsigned char colUC = 0;
   unsigned char rowUC = 0;


   alnSetST->lineWrapUI = 58; /*wrap every 58 characters*/
   alnSetST->pFullAlnBl = 0;  /*remove masked ends*/
   alnSetST->pBasePosBl = 1;  /*print base coordinates*/

   changeGap_alnSet(
      alnSetST,
      def_gapOpen_alnDefs,
      def_gapExtend_alnDefs
   ); /*add in gap scoring*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^  - Initialize scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*loop for all columns in the comparison matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetST->scoreMatrixSS[colUC][rowUC] = 0;
           /*Most of these cells will never be used*/
           /*But are needed to build the table*/
   } /*loop for all columns in the comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^  - Initialize match matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Both words, DNA, and AA all are the same when both
    ` characters/bases/amino acids are the same.
    */
    for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
    { /*Loop: fill in comparison matrix*/
        for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
        { /*Loop: Fill in the entire matrix*/ 
           if(colUC == rowUC)
             alnSetST->matchMatrixSC[colUC][rowUC] =
                def_ntEql_alnDefs;
           else
             alnSetST->matchMatrixSC[colUC][rowUC] =
                def_ntNotEql_alnDefs;
        } /*Loop: Fill in the entire matrix*/ 
    } /*Loop: fill in comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec04:
   ^   - set up scoring matrix for nucleotides
   ^   o fun14 sec04 sub01:
   ^     - a as first base
   ^   o fun14 sec04 sub02:
   ^     - t as first base
   ^   o fun14 sec04 sub03:
   ^     - u (t) as first base
   ^   o fun14 sec04 sub04:
   ^     - g as first base
   ^   o fun14 sec04 sub05:
   ^     - c as first base
   ^   o fun14 sec04 sub06:
   ^     - w (anonymous) as first base
   ^   o fun14 sec04 sub07:
   ^     - s (anonymous) as first base
   ^   o fun14 sec04 sub08:
   ^     - m (anonymous) as first base
   ^   o fun14 sec04 sub09:
   ^     - k (anonymous) as first base
   ^   o fun14 sec04 sub10:
   ^     - r (anonymous) as first base
   ^   o fun14 sec04 sub11:
   ^     - y (anonymous) as first base
   ^   o fun14 sec04 sub12:
   ^     - b (anonymous) as first base
   ^   o fun14 sec04 sub13:
   ^     - d (anonymous) as first base
   ^   o fun14 sec04 sub14:
   ^     - h (anonymous) as first base
   ^   o fun14 sec04 sub15:
   ^     - v (anonymous) as first base
   ^   o fun14 sec04 sub16:
   ^     - n (anonymous) as first base
   ^   o fun14 sec04 sub17:
   ^     - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec04 Sub01:
   *   - a as first base
   \*****************************************************/

   setScore_alnSet('a','a',def_AToA_alnDefs,alnSetST);
   setScore_alnSet('a','t',def_AToT_alnDefs,alnSetST);
   setScore_alnSet('a','u',def_AToU_alnDefs,alnSetST);
   setScore_alnSet('a','g',def_AToG_alnDefs,alnSetST);
   setScore_alnSet('a','c',def_AToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('a','w',def_AToW_alnDefs,alnSetST);
   setScore_alnSet('a','s',def_AToS_alnDefs,alnSetST);
   setScore_alnSet('a','m',def_AToM_alnDefs,alnSetST);
   setScore_alnSet('a','k',def_AToK_alnDefs,alnSetST);
   setScore_alnSet('a','r',def_AToR_alnDefs,alnSetST);
   setScore_alnSet('a','y',def_AToY_alnDefs,alnSetST);
   setScore_alnSet('a','b',def_AToB_alnDefs,alnSetST);
   setScore_alnSet('a','d',def_AToD_alnDefs,alnSetST);
   setScore_alnSet('a','h',def_AToH_alnDefs,alnSetST);
   setScore_alnSet('a','v',def_AToV_alnDefs,alnSetST);
   setScore_alnSet('a','n',def_AToN_alnDefs,alnSetST);
   setScore_alnSet('a','x',def_AToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub02:
   *   - t as first base
   \*****************************************************/

   setScore_alnSet('t','a',def_TToA_alnDefs,alnSetST);
   setScore_alnSet('t','t',def_TToT_alnDefs,alnSetST);
   setScore_alnSet('t','u',def_TToU_alnDefs,alnSetST);
   setScore_alnSet('t','g',def_TToG_alnDefs,alnSetST);
   setScore_alnSet('t','c',def_TToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('t','w',def_TToW_alnDefs,alnSetST);
   setScore_alnSet('t','s',def_TToS_alnDefs,alnSetST);
   setScore_alnSet('t','m',def_TToM_alnDefs,alnSetST);
   setScore_alnSet('t','k',def_TToK_alnDefs,alnSetST);
   setScore_alnSet('t','r',def_TToR_alnDefs,alnSetST);
   setScore_alnSet('t','y',def_TToY_alnDefs,alnSetST);
   setScore_alnSet('t','b',def_TToB_alnDefs,alnSetST);
   setScore_alnSet('t','d',def_TToD_alnDefs,alnSetST);
   setScore_alnSet('t','h',def_TToH_alnDefs,alnSetST);
   setScore_alnSet('t','v',def_TToV_alnDefs,alnSetST);
   setScore_alnSet('t','n',def_TToN_alnDefs,alnSetST);
   setScore_alnSet('t','x',def_TToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setScore_alnSet('u','a',def_UToA_alnDefs,alnSetST);
   setScore_alnSet('u','t',def_UToT_alnDefs,alnSetST);
   setScore_alnSet('u','u',def_UToU_alnDefs,alnSetST);
   setScore_alnSet('u','g',def_UToG_alnDefs,alnSetST);
   setScore_alnSet('u','c',def_UToC_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setScore_alnSet('u','w',def_UToW_alnDefs,alnSetST);
   setScore_alnSet('u','s',def_UToS_alnDefs,alnSetST);
   setScore_alnSet('u','m',def_UToM_alnDefs,alnSetST);
   setScore_alnSet('u','k',def_UToK_alnDefs,alnSetST);
   setScore_alnSet('u','r',def_UToR_alnDefs,alnSetST);
   setScore_alnSet('u','y',def_UToY_alnDefs,alnSetST);
   setScore_alnSet('u','b',def_UToB_alnDefs,alnSetST);
   setScore_alnSet('u','d',def_UToD_alnDefs,alnSetST);
   setScore_alnSet('u','h',def_UToH_alnDefs,alnSetST);
   setScore_alnSet('u','v',def_UToV_alnDefs,alnSetST);
   setScore_alnSet('u','n',def_UToN_alnDefs,alnSetST);
   setScore_alnSet('u','x',def_UToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub04:
   *   - g as first base
   \*****************************************************/

   setScore_alnSet('g','a',def_GToA_alnDefs,alnSetST);
   setScore_alnSet('g','t',def_GToT_alnDefs,alnSetST);
   setScore_alnSet('g','u',def_GToU_alnDefs,alnSetST);
   setScore_alnSet('g','g',def_GToG_alnDefs,alnSetST);
   setScore_alnSet('g','c',def_GToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('g','w',def_GToW_alnDefs,alnSetST);
   setScore_alnSet('g','s',def_GToS_alnDefs,alnSetST);
   setScore_alnSet('g','m',def_GToM_alnDefs,alnSetST);
   setScore_alnSet('g','k',def_GToK_alnDefs,alnSetST);
   setScore_alnSet('g','r',def_GToR_alnDefs,alnSetST);
   setScore_alnSet('g','y',def_GToY_alnDefs,alnSetST);
   setScore_alnSet('g','b',def_GToB_alnDefs,alnSetST);
   setScore_alnSet('g','d',def_GToD_alnDefs,alnSetST);
   setScore_alnSet('g','h',def_GToH_alnDefs,alnSetST);
   setScore_alnSet('g','v',def_GToV_alnDefs,alnSetST);
   setScore_alnSet('g','n',def_GToN_alnDefs,alnSetST);
   setScore_alnSet('g','x',def_GToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub05:
   *   - c as first base
   \*****************************************************/

   setScore_alnSet('c','a',def_CToA_alnDefs,alnSetST);
   setScore_alnSet('c','t',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','u',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','g',def_CToG_alnDefs,alnSetST);
   setScore_alnSet('c','c',def_CToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('c','w',def_CToW_alnDefs,alnSetST);
   setScore_alnSet('c','s',def_CToS_alnDefs,alnSetST);
   setScore_alnSet('c','m',def_CToM_alnDefs,alnSetST);
   setScore_alnSet('c','k',def_CToK_alnDefs,alnSetST);
   setScore_alnSet('c','r',def_CToR_alnDefs,alnSetST);
   setScore_alnSet('c','y',def_CToY_alnDefs,alnSetST);
   setScore_alnSet('c','b',def_CToB_alnDefs,alnSetST);
   setScore_alnSet('c','d',def_CToD_alnDefs,alnSetST);
   setScore_alnSet('c','h',def_CToH_alnDefs,alnSetST);
   setScore_alnSet('c','v',def_CToV_alnDefs,alnSetST);
   setScore_alnSet('c','n',def_CToN_alnDefs,alnSetST);
   setScore_alnSet('c','x',def_CToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('w','a',def_WToA_alnDefs,alnSetST);
   setScore_alnSet('w','c',def_WToC_alnDefs,alnSetST);
   setScore_alnSet('w','g',def_WToG_alnDefs,alnSetST);
   setScore_alnSet('w','t',def_WToT_alnDefs,alnSetST);
   setScore_alnSet('w','u',def_WToT_alnDefs,alnSetST);

   setScore_alnSet('w','w',def_WToW_alnDefs,alnSetST);
   setScore_alnSet('w','s',def_WToS_alnDefs,alnSetST);
   setScore_alnSet('w','m',def_WToM_alnDefs,alnSetST);
   setScore_alnSet('w','k',def_WToK_alnDefs,alnSetST);
   setScore_alnSet('w','r',def_WToR_alnDefs,alnSetST);
   setScore_alnSet('w','y',def_WToY_alnDefs,alnSetST);
   setScore_alnSet('w','b',def_WToB_alnDefs,alnSetST);
   setScore_alnSet('w','d',def_WToD_alnDefs,alnSetST);
   setScore_alnSet('w','h',def_WToH_alnDefs,alnSetST);
   setScore_alnSet('w','v',def_WToV_alnDefs,alnSetST);
   setScore_alnSet('w','n',def_WToN_alnDefs,alnSetST);
   setScore_alnSet('w','x',def_WToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('s','a',def_SToA_alnDefs,alnSetST);
   setScore_alnSet('s','c',def_SToC_alnDefs,alnSetST);
   setScore_alnSet('s','g',def_SToG_alnDefs,alnSetST);
   setScore_alnSet('s','t',def_SToT_alnDefs,alnSetST);
   setScore_alnSet('s','u',def_SToT_alnDefs,alnSetST);

   setScore_alnSet('s','w',def_SToW_alnDefs,alnSetST);
   setScore_alnSet('s','s',def_SToS_alnDefs,alnSetST);
   setScore_alnSet('s','m',def_SToM_alnDefs,alnSetST);
   setScore_alnSet('s','k',def_SToK_alnDefs,alnSetST);
   setScore_alnSet('s','r',def_SToR_alnDefs,alnSetST);
   setScore_alnSet('s','y',def_SToY_alnDefs,alnSetST);
   setScore_alnSet('s','b',def_SToB_alnDefs,alnSetST);
   setScore_alnSet('s','d',def_SToD_alnDefs,alnSetST);
   setScore_alnSet('s','h',def_SToH_alnDefs,alnSetST);
   setScore_alnSet('s','v',def_SToV_alnDefs,alnSetST);
   setScore_alnSet('s','n',def_SToN_alnDefs,alnSetST);
   setScore_alnSet('s','x',def_SToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('m','a',def_MToA_alnDefs,alnSetST);
   setScore_alnSet('m','c',def_MToC_alnDefs,alnSetST);
   setScore_alnSet('m','g',def_MToG_alnDefs,alnSetST);
   setScore_alnSet('m','t',def_MToT_alnDefs,alnSetST);
   setScore_alnSet('m','u',def_MToT_alnDefs,alnSetST);

   setScore_alnSet('m','w',def_MToW_alnDefs,alnSetST);
   setScore_alnSet('m','s',def_MToS_alnDefs,alnSetST);
   setScore_alnSet('m','m',def_MToM_alnDefs,alnSetST);
   setScore_alnSet('m','k',def_MToK_alnDefs,alnSetST);
   setScore_alnSet('m','r',def_MToR_alnDefs,alnSetST);
   setScore_alnSet('m','y',def_MToY_alnDefs,alnSetST);
   setScore_alnSet('m','b',def_MToB_alnDefs,alnSetST);
   setScore_alnSet('m','d',def_MToD_alnDefs,alnSetST);
   setScore_alnSet('m','h',def_MToH_alnDefs,alnSetST);
   setScore_alnSet('m','v',def_MToV_alnDefs,alnSetST);
   setScore_alnSet('m','n',def_MToN_alnDefs,alnSetST);
   setScore_alnSet('m','x',def_MToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('k','a',def_KToA_alnDefs,alnSetST);
   setScore_alnSet('k','c',def_KToC_alnDefs,alnSetST);
   setScore_alnSet('k','g',def_KToG_alnDefs,alnSetST);
   setScore_alnSet('k','t',def_KToT_alnDefs,alnSetST);
   setScore_alnSet('k','u',def_KToT_alnDefs,alnSetST);

   setScore_alnSet('k','w',def_KToW_alnDefs,alnSetST);
   setScore_alnSet('k','s',def_KToS_alnDefs,alnSetST);
   setScore_alnSet('k','m',def_KToM_alnDefs,alnSetST);
   setScore_alnSet('k','k',def_KToK_alnDefs,alnSetST);
   setScore_alnSet('k','r',def_KToR_alnDefs,alnSetST);
   setScore_alnSet('k','y',def_KToY_alnDefs,alnSetST);
   setScore_alnSet('k','b',def_KToB_alnDefs,alnSetST);
   setScore_alnSet('k','d',def_KToD_alnDefs,alnSetST);
   setScore_alnSet('k','h',def_KToH_alnDefs,alnSetST);
   setScore_alnSet('k','v',def_KToV_alnDefs,alnSetST);
   setScore_alnSet('k','n',def_KToN_alnDefs,alnSetST);
   setScore_alnSet('k','x',def_KToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('r','a',def_RToA_alnDefs,alnSetST);
   setScore_alnSet('r','c',def_RToC_alnDefs,alnSetST);
   setScore_alnSet('r','g',def_RToG_alnDefs,alnSetST);
   setScore_alnSet('r','t',def_RToT_alnDefs,alnSetST);
   setScore_alnSet('r','u',def_RToT_alnDefs,alnSetST);

   setScore_alnSet('r','w',def_RToW_alnDefs,alnSetST);
   setScore_alnSet('r','s',def_RToS_alnDefs,alnSetST);
   setScore_alnSet('r','m',def_RToM_alnDefs,alnSetST);
   setScore_alnSet('r','k',def_RToK_alnDefs,alnSetST);
   setScore_alnSet('r','r',def_RToR_alnDefs,alnSetST);
   setScore_alnSet('r','y',def_RToY_alnDefs,alnSetST);
   setScore_alnSet('r','b',def_RToB_alnDefs,alnSetST);
   setScore_alnSet('r','d',def_RToD_alnDefs,alnSetST);
   setScore_alnSet('r','h',def_RToH_alnDefs,alnSetST);
   setScore_alnSet('r','v',def_RToV_alnDefs,alnSetST);
   setScore_alnSet('r','n',def_RToN_alnDefs,alnSetST);
   setScore_alnSet('r','x',def_RToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('y','a',def_YToA_alnDefs,alnSetST);
   setScore_alnSet('y','c',def_YToC_alnDefs,alnSetST);
   setScore_alnSet('y','g',def_YToG_alnDefs,alnSetST);
   setScore_alnSet('y','t',def_YToT_alnDefs,alnSetST);
   setScore_alnSet('y','u',def_YToT_alnDefs,alnSetST);

   setScore_alnSet('y','w',def_YToW_alnDefs,alnSetST);
   setScore_alnSet('y','s',def_YToS_alnDefs,alnSetST);
   setScore_alnSet('y','m',def_YToM_alnDefs,alnSetST);
   setScore_alnSet('y','k',def_YToK_alnDefs,alnSetST);
   setScore_alnSet('y','r',def_YToR_alnDefs,alnSetST);
   setScore_alnSet('y','y',def_YToY_alnDefs,alnSetST);
   setScore_alnSet('y','b',def_YToB_alnDefs,alnSetST);
   setScore_alnSet('y','d',def_YToD_alnDefs,alnSetST);
   setScore_alnSet('y','h',def_YToH_alnDefs,alnSetST);
   setScore_alnSet('y','v',def_YToV_alnDefs,alnSetST);
   setScore_alnSet('y','n',def_YToN_alnDefs,alnSetST);
   setScore_alnSet('y','x',def_YToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('b','a',def_BToA_alnDefs,alnSetST);
   setScore_alnSet('b','c',def_BToC_alnDefs,alnSetST);
   setScore_alnSet('b','g',def_BToG_alnDefs,alnSetST);
   setScore_alnSet('b','t',def_BToT_alnDefs,alnSetST);
   setScore_alnSet('b','u',def_BToT_alnDefs,alnSetST);

   setScore_alnSet('b','w',def_BToW_alnDefs,alnSetST);
   setScore_alnSet('b','s',def_BToS_alnDefs,alnSetST);
   setScore_alnSet('b','m',def_BToM_alnDefs,alnSetST);
   setScore_alnSet('b','k',def_BToK_alnDefs,alnSetST);
   setScore_alnSet('b','r',def_BToR_alnDefs,alnSetST);
   setScore_alnSet('b','y',def_BToY_alnDefs,alnSetST);
   setScore_alnSet('b','b',def_BToB_alnDefs,alnSetST);
   setScore_alnSet('b','d',def_BToD_alnDefs,alnSetST);
   setScore_alnSet('b','h',def_BToH_alnDefs,alnSetST);
   setScore_alnSet('b','v',def_BToV_alnDefs,alnSetST);
   setScore_alnSet('b','n',def_BToN_alnDefs,alnSetST);
   setScore_alnSet('b','x',def_BToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('d','a',def_DToA_alnDefs,alnSetST);
   setScore_alnSet('d','c',def_DToC_alnDefs,alnSetST);
   setScore_alnSet('d','g',def_DToG_alnDefs,alnSetST);
   setScore_alnSet('d','t',def_DToT_alnDefs,alnSetST);
   setScore_alnSet('d','u',def_DToT_alnDefs,alnSetST);

   setScore_alnSet('d','w',def_DToW_alnDefs,alnSetST);
   setScore_alnSet('d','s',def_DToS_alnDefs,alnSetST);
   setScore_alnSet('d','m',def_DToM_alnDefs,alnSetST);
   setScore_alnSet('d','k',def_DToK_alnDefs,alnSetST);
   setScore_alnSet('d','r',def_DToR_alnDefs,alnSetST);
   setScore_alnSet('d','y',def_DToY_alnDefs,alnSetST);
   setScore_alnSet('d','b',def_DToB_alnDefs,alnSetST);
   setScore_alnSet('d','d',def_DToD_alnDefs,alnSetST);
   setScore_alnSet('d','h',def_DToH_alnDefs,alnSetST);
   setScore_alnSet('d','v',def_DToV_alnDefs,alnSetST);
   setScore_alnSet('d','n',def_DToN_alnDefs,alnSetST);
   setScore_alnSet('d','x',def_DToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('h','a',def_HToA_alnDefs,alnSetST);
   setScore_alnSet('h','c',def_HToC_alnDefs,alnSetST);
   setScore_alnSet('h','g',def_HToG_alnDefs,alnSetST);
   setScore_alnSet('h','t',def_HToT_alnDefs,alnSetST);
   setScore_alnSet('h','u',def_HToT_alnDefs,alnSetST);

   setScore_alnSet('h','w',def_HToW_alnDefs,alnSetST);
   setScore_alnSet('h','s',def_HToS_alnDefs,alnSetST);
   setScore_alnSet('h','m',def_HToM_alnDefs,alnSetST);
   setScore_alnSet('h','k',def_HToK_alnDefs,alnSetST);
   setScore_alnSet('h','r',def_HToR_alnDefs,alnSetST);
   setScore_alnSet('h','y',def_HToY_alnDefs,alnSetST);
   setScore_alnSet('h','b',def_HToB_alnDefs,alnSetST);
   setScore_alnSet('h','d',def_HToD_alnDefs,alnSetST);
   setScore_alnSet('h','h',def_HToH_alnDefs,alnSetST);
   setScore_alnSet('h','v',def_HToV_alnDefs,alnSetST);
   setScore_alnSet('h','n',def_HToN_alnDefs,alnSetST);
   setScore_alnSet('h','x',def_HToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('v','a',def_VToA_alnDefs,alnSetST);
   setScore_alnSet('v','c',def_VToC_alnDefs,alnSetST);
   setScore_alnSet('v','g',def_VToG_alnDefs,alnSetST);
   setScore_alnSet('v','t',def_VToT_alnDefs,alnSetST);
   setScore_alnSet('v','u',def_VToT_alnDefs,alnSetST);

   setScore_alnSet('v','w',def_VToW_alnDefs,alnSetST);
   setScore_alnSet('v','s',def_VToS_alnDefs,alnSetST);
   setScore_alnSet('v','m',def_VToM_alnDefs,alnSetST);
   setScore_alnSet('v','k',def_VToK_alnDefs,alnSetST);
   setScore_alnSet('v','r',def_VToR_alnDefs,alnSetST);
   setScore_alnSet('v','y',def_VToY_alnDefs,alnSetST);
   setScore_alnSet('v','b',def_VToB_alnDefs,alnSetST);
   setScore_alnSet('v','d',def_VToD_alnDefs,alnSetST);
   setScore_alnSet('v','h',def_VToH_alnDefs,alnSetST);
   setScore_alnSet('v','v',def_VToV_alnDefs,alnSetST);
   setScore_alnSet('v','n',def_VToN_alnDefs,alnSetST);
   setScore_alnSet('v','x',def_VToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('n','a',def_NToA_alnDefs,alnSetST);
   setScore_alnSet('n','c',def_NToC_alnDefs,alnSetST);
   setScore_alnSet('n','g',def_NToG_alnDefs,alnSetST);
   setScore_alnSet('n','t',def_NToT_alnDefs,alnSetST);
   setScore_alnSet('n','u',def_NToT_alnDefs,alnSetST);

   setScore_alnSet('n','w',def_NToW_alnDefs,alnSetST);
   setScore_alnSet('n','s',def_NToS_alnDefs,alnSetST);
   setScore_alnSet('n','m',def_NToM_alnDefs,alnSetST);
   setScore_alnSet('n','k',def_NToK_alnDefs,alnSetST);
   setScore_alnSet('n','r',def_NToR_alnDefs,alnSetST);
   setScore_alnSet('n','y',def_NToY_alnDefs,alnSetST);
   setScore_alnSet('n','b',def_NToB_alnDefs,alnSetST);
   setScore_alnSet('n','d',def_NToD_alnDefs,alnSetST);
   setScore_alnSet('n','h',def_NToH_alnDefs,alnSetST);
   setScore_alnSet('n','v',def_NToV_alnDefs,alnSetST);
   setScore_alnSet('n','n',def_NToN_alnDefs,alnSetST);
   setScore_alnSet('n','x',def_NToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec04 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setScore_alnSet('x','a',def_XToA_alnDefs,alnSetST);
   setScore_alnSet('x','c',def_XToC_alnDefs,alnSetST);
   setScore_alnSet('x','g',def_XToG_alnDefs,alnSetST);
   setScore_alnSet('x','t',def_XToT_alnDefs,alnSetST);
   setScore_alnSet('x','u',def_XToT_alnDefs,alnSetST);

   setScore_alnSet('x','w',def_XToW_alnDefs,alnSetST);
   setScore_alnSet('x','s',def_XToS_alnDefs,alnSetST);
   setScore_alnSet('x','m',def_XToM_alnDefs,alnSetST);
   setScore_alnSet('x','k',def_XToK_alnDefs,alnSetST);
   setScore_alnSet('x','r',def_XToR_alnDefs,alnSetST);
   setScore_alnSet('x','y',def_XToY_alnDefs,alnSetST);
   setScore_alnSet('x','b',def_XToB_alnDefs,alnSetST);
   setScore_alnSet('x','d',def_XToD_alnDefs,alnSetST);
   setScore_alnSet('x','h',def_XToH_alnDefs,alnSetST);
   setScore_alnSet('x','v',def_XToV_alnDefs,alnSetST);
   setScore_alnSet('x','n',def_XToN_alnDefs,alnSetST);
   setScore_alnSet('x','x',def_XToX_alnDefs,alnSetST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^  - set up matching matrix for nucleotides
   ^  o fun14 sec05 sub01:
   ^    - a as first base
   ^  o fun14 sec05 sub02:
   ^    - t as first base
   ^  o fun14 sec05 sub03:
   ^    - u (t) as first base
   ^  o fun14 sec05 sub04:
   ^    - g as first base
   ^  o fun14 sec05 sub05:
   ^    - c as first base
   ^  o fun14 sec05 sub06:
   ^    - w (anonymous) as first base
   ^  o fun14 sec05 sub07:
   ^    - s (anonymous) as first base
   ^  o fun14 sec05 sub08:
   ^    - m (anonymous) as first base
   ^  o fun14 sec05 sub09:
   ^    - k (anonymous) as first base
   ^  o fun14 sec05 sub10:
   ^    - r (anonymous) as first base
   ^  o fun14 sec05 sub11:
   ^    - y (anonymous) as first base
   ^  o fun14 sec05 sub12:
   ^    - b (anonymous) as first base
   ^  o fun14 sec05 sub13:
   ^    - d (anonymous) as first base
   ^  o fun14 sec05 sub14:
   ^    - h (anonymous) as first base
   ^  o fun14 sec05 sub15:
   ^    - v (anonymous) as first base
   ^  o fun14 sec05 sub16:
   ^    - n (anonymous) as first base
   ^  o fun14 sec05 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec05 Sub01:
   *   - a as first base
   \*****************************************************/

   setMatch_alnSet('a','a',def_AEqlA_alnDefs,alnSetST);
   setMatch_alnSet('a','t',def_AEqlT_alnDefs,alnSetST);
   setMatch_alnSet('a','u',def_AEqlU_alnDefs,alnSetST);
   setMatch_alnSet('a','g',def_AEqlG_alnDefs,alnSetST);
   setMatch_alnSet('a','c',def_AEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('a','w',def_AEqlW_alnDefs,alnSetST);
   setMatch_alnSet('a','s',def_AEqlS_alnDefs,alnSetST);
   setMatch_alnSet('a','m',def_AEqlM_alnDefs,alnSetST);
   setMatch_alnSet('a','k',def_AEqlK_alnDefs,alnSetST);
   setMatch_alnSet('a','r',def_AEqlR_alnDefs,alnSetST);
   setMatch_alnSet('a','y',def_AEqlY_alnDefs,alnSetST);
   setMatch_alnSet('a','b',def_AEqlB_alnDefs,alnSetST);
   setMatch_alnSet('a','d',def_AEqlD_alnDefs,alnSetST);
   setMatch_alnSet('a','h',def_AEqlH_alnDefs,alnSetST);
   setMatch_alnSet('a','v',def_AEqlV_alnDefs,alnSetST);
   setMatch_alnSet('a','n',def_AEqlN_alnDefs,alnSetST);
   setMatch_alnSet('a','x',def_AEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub02:
   *   - t as first base
   \*****************************************************/

   setMatch_alnSet('t','a',def_TEqlA_alnDefs,alnSetST);
   setMatch_alnSet('t','t',def_TEqlT_alnDefs,alnSetST);
   setMatch_alnSet('t','u',def_TEqlU_alnDefs,alnSetST);
   setMatch_alnSet('t','g',def_TEqlG_alnDefs,alnSetST);
   setMatch_alnSet('t','c',def_TEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('t','w',def_TEqlW_alnDefs,alnSetST);
   setMatch_alnSet('t','s',def_TEqlS_alnDefs,alnSetST);
   setMatch_alnSet('t','m',def_TEqlM_alnDefs,alnSetST);
   setMatch_alnSet('t','k',def_TEqlK_alnDefs,alnSetST);
   setMatch_alnSet('t','r',def_TEqlR_alnDefs,alnSetST);
   setMatch_alnSet('t','y',def_TEqlY_alnDefs,alnSetST);
   setMatch_alnSet('t','b',def_TEqlB_alnDefs,alnSetST);
   setMatch_alnSet('t','d',def_TEqlD_alnDefs,alnSetST);
   setMatch_alnSet('t','h',def_TEqlH_alnDefs,alnSetST);
   setMatch_alnSet('t','v',def_TEqlV_alnDefs,alnSetST);
   setMatch_alnSet('t','n',def_TEqlN_alnDefs,alnSetST);
   setMatch_alnSet('t','x',def_TEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setMatch_alnSet('u','a',def_UEqlA_alnDefs,alnSetST);
   setMatch_alnSet('u','g',def_UEqlG_alnDefs,alnSetST);
   setMatch_alnSet('u','c',def_UEqlC_alnDefs,alnSetST);
   setMatch_alnSet('u','t',def_UEqlT_alnDefs,alnSetST);
   setMatch_alnSet('u','u',def_UEqlU_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setMatch_alnSet('u','w',def_UEqlW_alnDefs,alnSetST);
   setMatch_alnSet('u','s',def_UEqlS_alnDefs,alnSetST);
   setMatch_alnSet('u','m',def_UEqlM_alnDefs,alnSetST);
   setMatch_alnSet('u','k',def_UEqlK_alnDefs,alnSetST);
   setMatch_alnSet('u','r',def_UEqlR_alnDefs,alnSetST);
   setMatch_alnSet('u','y',def_UEqlY_alnDefs,alnSetST);
   setMatch_alnSet('u','b',def_UEqlB_alnDefs,alnSetST);
   setMatch_alnSet('u','d',def_UEqlD_alnDefs,alnSetST);
   setMatch_alnSet('u','h',def_UEqlH_alnDefs,alnSetST);
   setMatch_alnSet('u','v',def_UEqlV_alnDefs,alnSetST);
   setMatch_alnSet('u','n',def_UEqlN_alnDefs,alnSetST);
   setMatch_alnSet('u','x',def_UEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub04:
   *   - g as first base
   \*****************************************************/

   setMatch_alnSet('g','a',def_GEqlA_alnDefs,alnSetST);
   setMatch_alnSet('g','t',def_GEqlT_alnDefs,alnSetST);
   setMatch_alnSet('g','u',def_GEqlU_alnDefs,alnSetST);
   setMatch_alnSet('g','g',def_GEqlG_alnDefs,alnSetST);
   setMatch_alnSet('g','c',def_GEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('g','w',def_GEqlW_alnDefs,alnSetST);
   setMatch_alnSet('g','s',def_GEqlS_alnDefs,alnSetST);
   setMatch_alnSet('g','m',def_GEqlM_alnDefs,alnSetST);
   setMatch_alnSet('g','k',def_GEqlK_alnDefs,alnSetST);
   setMatch_alnSet('g','r',def_GEqlR_alnDefs,alnSetST);
   setMatch_alnSet('g','y',def_GEqlY_alnDefs,alnSetST);
   setMatch_alnSet('g','b',def_GEqlB_alnDefs,alnSetST);
   setMatch_alnSet('g','d',def_GEqlD_alnDefs,alnSetST);
   setMatch_alnSet('g','h',def_GEqlH_alnDefs,alnSetST);
   setMatch_alnSet('g','v',def_GEqlV_alnDefs,alnSetST);
   setMatch_alnSet('g','n',def_GEqlN_alnDefs,alnSetST);
   setMatch_alnSet('g','x',def_GEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub05:
   *   - c as first base
   \*****************************************************/

   setMatch_alnSet('c','a',def_CEqlA_alnDefs,alnSetST);
   setMatch_alnSet('c','t',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','u',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','g',def_CEqlG_alnDefs,alnSetST);
   setMatch_alnSet('c','c',def_CEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('c','w',def_CEqlW_alnDefs,alnSetST);
   setMatch_alnSet('c','s',def_CEqlS_alnDefs,alnSetST);
   setMatch_alnSet('c','m',def_CEqlM_alnDefs,alnSetST);
   setMatch_alnSet('c','k',def_CEqlK_alnDefs,alnSetST);
   setMatch_alnSet('c','r',def_CEqlR_alnDefs,alnSetST);
   setMatch_alnSet('c','y',def_CEqlY_alnDefs,alnSetST);
   setMatch_alnSet('c','b',def_CEqlB_alnDefs,alnSetST);
   setMatch_alnSet('c','d',def_CEqlD_alnDefs,alnSetST);
   setMatch_alnSet('c','h',def_CEqlH_alnDefs,alnSetST);
   setMatch_alnSet('c','v',def_CEqlV_alnDefs,alnSetST);
   setMatch_alnSet('c','n',def_CEqlN_alnDefs,alnSetST);
   setMatch_alnSet('c','x',def_CEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('w','a',def_WEqlA_alnDefs,alnSetST);
   setMatch_alnSet('w','c',def_WEqlC_alnDefs,alnSetST);
   setMatch_alnSet('w','g',def_WEqlG_alnDefs,alnSetST);
   setMatch_alnSet('w','t',def_WEqlT_alnDefs,alnSetST);
   setMatch_alnSet('w','u',def_WEqlT_alnDefs,alnSetST);

   setMatch_alnSet('w','w',def_WEqlW_alnDefs,alnSetST);
   setMatch_alnSet('w','s',def_WEqlS_alnDefs,alnSetST);
   setMatch_alnSet('w','m',def_WEqlM_alnDefs,alnSetST);
   setMatch_alnSet('w','k',def_WEqlK_alnDefs,alnSetST);
   setMatch_alnSet('w','r',def_WEqlR_alnDefs,alnSetST);
   setMatch_alnSet('w','y',def_WEqlY_alnDefs,alnSetST);
   setMatch_alnSet('w','b',def_WEqlB_alnDefs,alnSetST);
   setMatch_alnSet('w','d',def_WEqlD_alnDefs,alnSetST);
   setMatch_alnSet('w','h',def_WEqlH_alnDefs,alnSetST);
   setMatch_alnSet('w','v',def_WEqlV_alnDefs,alnSetST);
   setMatch_alnSet('w','n',def_WEqlN_alnDefs,alnSetST);
   setMatch_alnSet('w','x',def_WEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('s','a',def_SEqlA_alnDefs,alnSetST);
   setMatch_alnSet('s','c',def_SEqlC_alnDefs,alnSetST);
   setMatch_alnSet('s','g',def_SEqlG_alnDefs,alnSetST);
   setMatch_alnSet('s','t',def_SEqlT_alnDefs,alnSetST);
   setMatch_alnSet('s','u',def_SEqlT_alnDefs,alnSetST);

   setMatch_alnSet('s','w',def_SEqlW_alnDefs,alnSetST);
   setMatch_alnSet('s','s',def_SEqlS_alnDefs,alnSetST);
   setMatch_alnSet('s','m',def_SEqlM_alnDefs,alnSetST);
   setMatch_alnSet('s','k',def_SEqlK_alnDefs,alnSetST);
   setMatch_alnSet('s','r',def_SEqlR_alnDefs,alnSetST);
   setMatch_alnSet('s','y',def_SEqlY_alnDefs,alnSetST);
   setMatch_alnSet('s','b',def_SEqlB_alnDefs,alnSetST);
   setMatch_alnSet('s','d',def_SEqlD_alnDefs,alnSetST);
   setMatch_alnSet('s','h',def_SEqlH_alnDefs,alnSetST);
   setMatch_alnSet('s','v',def_SEqlV_alnDefs,alnSetST);
   setMatch_alnSet('s','n',def_SEqlN_alnDefs,alnSetST);
   setMatch_alnSet('s','x',def_SEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('m','a',def_MEqlA_alnDefs,alnSetST);
   setMatch_alnSet('m','c',def_MEqlC_alnDefs,alnSetST);
   setMatch_alnSet('m','g',def_MEqlG_alnDefs,alnSetST);
   setMatch_alnSet('m','t',def_MEqlT_alnDefs,alnSetST);
   setMatch_alnSet('m','u',def_MEqlT_alnDefs,alnSetST);

   setMatch_alnSet('m','w',def_MEqlW_alnDefs,alnSetST);
   setMatch_alnSet('m','s',def_MEqlS_alnDefs,alnSetST);
   setMatch_alnSet('m','m',def_MEqlM_alnDefs,alnSetST);
   setMatch_alnSet('m','k',def_MEqlK_alnDefs,alnSetST);
   setMatch_alnSet('m','r',def_MEqlR_alnDefs,alnSetST);
   setMatch_alnSet('m','y',def_MEqlY_alnDefs,alnSetST);
   setMatch_alnSet('m','b',def_MEqlB_alnDefs,alnSetST);
   setMatch_alnSet('m','d',def_MEqlD_alnDefs,alnSetST);
   setMatch_alnSet('m','h',def_MEqlH_alnDefs,alnSetST);
   setMatch_alnSet('m','v',def_MEqlV_alnDefs,alnSetST);
   setMatch_alnSet('m','n',def_MEqlN_alnDefs,alnSetST);
   setMatch_alnSet('m','x',def_MEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('k','a',def_KEqlA_alnDefs,alnSetST);
   setMatch_alnSet('k','c',def_KEqlC_alnDefs,alnSetST);
   setMatch_alnSet('k','g',def_KEqlG_alnDefs,alnSetST);
   setMatch_alnSet('k','t',def_KEqlT_alnDefs,alnSetST);
   setMatch_alnSet('k','u',def_KEqlT_alnDefs,alnSetST);

   setMatch_alnSet('k','w',def_KEqlW_alnDefs,alnSetST);
   setMatch_alnSet('k','s',def_KEqlS_alnDefs,alnSetST);
   setMatch_alnSet('k','m',def_KEqlM_alnDefs,alnSetST);
   setMatch_alnSet('k','k',def_KEqlK_alnDefs,alnSetST);
   setMatch_alnSet('k','r',def_KEqlR_alnDefs,alnSetST);
   setMatch_alnSet('k','y',def_KEqlY_alnDefs,alnSetST);
   setMatch_alnSet('k','b',def_KEqlB_alnDefs,alnSetST);
   setMatch_alnSet('k','d',def_KEqlD_alnDefs,alnSetST);
   setMatch_alnSet('k','h',def_KEqlH_alnDefs,alnSetST);
   setMatch_alnSet('k','v',def_KEqlV_alnDefs,alnSetST);
   setMatch_alnSet('k','n',def_KEqlN_alnDefs,alnSetST);
   setMatch_alnSet('k','x',def_KEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('r','a',def_REqlA_alnDefs,alnSetST);
   setMatch_alnSet('r','c',def_REqlC_alnDefs,alnSetST);
   setMatch_alnSet('r','g',def_REqlG_alnDefs,alnSetST);
   setMatch_alnSet('r','t',def_REqlT_alnDefs,alnSetST);
   setMatch_alnSet('r','u',def_REqlT_alnDefs,alnSetST);

   setMatch_alnSet('r','w',def_REqlW_alnDefs,alnSetST);
   setMatch_alnSet('r','s',def_REqlS_alnDefs,alnSetST);
   setMatch_alnSet('r','m',def_REqlM_alnDefs,alnSetST);
   setMatch_alnSet('r','k',def_REqlK_alnDefs,alnSetST);
   setMatch_alnSet('r','r',def_REqlR_alnDefs,alnSetST);
   setMatch_alnSet('r','y',def_REqlY_alnDefs,alnSetST);
   setMatch_alnSet('r','b',def_REqlB_alnDefs,alnSetST);
   setMatch_alnSet('r','d',def_REqlD_alnDefs,alnSetST);
   setMatch_alnSet('r','h',def_REqlH_alnDefs,alnSetST);
   setMatch_alnSet('r','v',def_REqlV_alnDefs,alnSetST);
   setMatch_alnSet('r','n',def_REqlN_alnDefs,alnSetST);
   setMatch_alnSet('r','x',def_REqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('y','a',def_YEqlA_alnDefs,alnSetST);
   setMatch_alnSet('y','c',def_YEqlC_alnDefs,alnSetST);
   setMatch_alnSet('y','g',def_YEqlG_alnDefs,alnSetST);
   setMatch_alnSet('y','t',def_YEqlT_alnDefs,alnSetST);
   setMatch_alnSet('y','u',def_YEqlT_alnDefs,alnSetST);

   setMatch_alnSet('y','w',def_YEqlW_alnDefs,alnSetST);
   setMatch_alnSet('y','s',def_YEqlS_alnDefs,alnSetST);
   setMatch_alnSet('y','m',def_YEqlM_alnDefs,alnSetST);
   setMatch_alnSet('y','k',def_YEqlK_alnDefs,alnSetST);
   setMatch_alnSet('y','r',def_YEqlR_alnDefs,alnSetST);
   setMatch_alnSet('y','y',def_YEqlY_alnDefs,alnSetST);
   setMatch_alnSet('y','b',def_YEqlB_alnDefs,alnSetST);
   setMatch_alnSet('y','d',def_YEqlD_alnDefs,alnSetST);
   setMatch_alnSet('y','h',def_YEqlH_alnDefs,alnSetST);
   setMatch_alnSet('y','v',def_YEqlV_alnDefs,alnSetST);
   setMatch_alnSet('y','n',def_YEqlN_alnDefs,alnSetST);
   setMatch_alnSet('y','x',def_YEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('b','a',def_BEqlA_alnDefs,alnSetST);
   setMatch_alnSet('b','c',def_BEqlC_alnDefs,alnSetST);
   setMatch_alnSet('b','g',def_BEqlG_alnDefs,alnSetST);
   setMatch_alnSet('b','t',def_BEqlT_alnDefs,alnSetST);
   setMatch_alnSet('b','u',def_BEqlT_alnDefs,alnSetST);

   setMatch_alnSet('b','w',def_BEqlW_alnDefs,alnSetST);
   setMatch_alnSet('b','s',def_BEqlS_alnDefs,alnSetST);
   setMatch_alnSet('b','m',def_BEqlM_alnDefs,alnSetST);
   setMatch_alnSet('b','k',def_BEqlK_alnDefs,alnSetST);
   setMatch_alnSet('b','r',def_BEqlR_alnDefs,alnSetST);
   setMatch_alnSet('b','y',def_BEqlY_alnDefs,alnSetST);
   setMatch_alnSet('b','b',def_BEqlB_alnDefs,alnSetST);
   setMatch_alnSet('b','d',def_BEqlD_alnDefs,alnSetST);
   setMatch_alnSet('b','h',def_BEqlH_alnDefs,alnSetST);
   setMatch_alnSet('b','v',def_BEqlV_alnDefs,alnSetST);
   setMatch_alnSet('b','n',def_BEqlN_alnDefs,alnSetST);
   setMatch_alnSet('b','x',def_BEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('d','a',def_DEqlA_alnDefs,alnSetST);
   setMatch_alnSet('d','c',def_DEqlC_alnDefs,alnSetST);
   setMatch_alnSet('d','g',def_DEqlG_alnDefs,alnSetST);
   setMatch_alnSet('d','t',def_DEqlT_alnDefs,alnSetST);
   setMatch_alnSet('d','u',def_DEqlT_alnDefs,alnSetST);

   setMatch_alnSet('d','w',def_DEqlW_alnDefs,alnSetST);
   setMatch_alnSet('d','s',def_DEqlS_alnDefs,alnSetST);
   setMatch_alnSet('d','m',def_DEqlM_alnDefs,alnSetST);
   setMatch_alnSet('d','k',def_DEqlK_alnDefs,alnSetST);
   setMatch_alnSet('d','r',def_DEqlR_alnDefs,alnSetST);
   setMatch_alnSet('d','y',def_DEqlY_alnDefs,alnSetST);
   setMatch_alnSet('d','b',def_DEqlB_alnDefs,alnSetST);
   setMatch_alnSet('d','d',def_DEqlD_alnDefs,alnSetST);
   setMatch_alnSet('d','h',def_DEqlH_alnDefs,alnSetST);
   setMatch_alnSet('d','v',def_DEqlV_alnDefs,alnSetST);
   setMatch_alnSet('d','n',def_DEqlN_alnDefs,alnSetST);
   setMatch_alnSet('d','x',def_DEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('h','a',def_HEqlA_alnDefs,alnSetST);
   setMatch_alnSet('h','c',def_HEqlC_alnDefs,alnSetST);
   setMatch_alnSet('h','g',def_HEqlG_alnDefs,alnSetST);
   setMatch_alnSet('h','t',def_HEqlT_alnDefs,alnSetST);
   setMatch_alnSet('h','u',def_HEqlT_alnDefs,alnSetST);

   setMatch_alnSet('h','w',def_HEqlW_alnDefs,alnSetST);
   setMatch_alnSet('h','s',def_HEqlS_alnDefs,alnSetST);
   setMatch_alnSet('h','m',def_HEqlM_alnDefs,alnSetST);
   setMatch_alnSet('h','k',def_HEqlK_alnDefs,alnSetST);
   setMatch_alnSet('h','r',def_HEqlR_alnDefs,alnSetST);
   setMatch_alnSet('h','y',def_HEqlY_alnDefs,alnSetST);
   setMatch_alnSet('h','b',def_HEqlB_alnDefs,alnSetST);
   setMatch_alnSet('h','d',def_HEqlD_alnDefs,alnSetST);
   setMatch_alnSet('h','h',def_HEqlH_alnDefs,alnSetST);
   setMatch_alnSet('h','v',def_HEqlV_alnDefs,alnSetST);
   setMatch_alnSet('h','n',def_HEqlN_alnDefs,alnSetST);
   setMatch_alnSet('h','x',def_HEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('v','a',def_VEqlA_alnDefs,alnSetST);
   setMatch_alnSet('v','c',def_VEqlC_alnDefs,alnSetST);
   setMatch_alnSet('v','g',def_VEqlG_alnDefs,alnSetST);
   setMatch_alnSet('v','t',def_VEqlT_alnDefs,alnSetST);
   setMatch_alnSet('v','u',def_VEqlT_alnDefs,alnSetST);

   setMatch_alnSet('v','w',def_VEqlW_alnDefs,alnSetST);
   setMatch_alnSet('v','s',def_VEqlS_alnDefs,alnSetST);
   setMatch_alnSet('v','m',def_VEqlM_alnDefs,alnSetST);
   setMatch_alnSet('v','k',def_VEqlK_alnDefs,alnSetST);
   setMatch_alnSet('v','r',def_VEqlR_alnDefs,alnSetST);
   setMatch_alnSet('v','y',def_VEqlY_alnDefs,alnSetST);
   setMatch_alnSet('v','b',def_VEqlB_alnDefs,alnSetST);
   setMatch_alnSet('v','d',def_VEqlD_alnDefs,alnSetST);
   setMatch_alnSet('v','h',def_VEqlH_alnDefs,alnSetST);
   setMatch_alnSet('v','v',def_VEqlV_alnDefs,alnSetST);
   setMatch_alnSet('v','n',def_VEqlN_alnDefs,alnSetST);
   setMatch_alnSet('v','x',def_VEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('n','a',def_NEqlA_alnDefs,alnSetST);
   setMatch_alnSet('n','c',def_NEqlC_alnDefs,alnSetST);
   setMatch_alnSet('n','g',def_NEqlG_alnDefs,alnSetST);
   setMatch_alnSet('n','t',def_NEqlT_alnDefs,alnSetST);
   setMatch_alnSet('n','u',def_NEqlT_alnDefs,alnSetST);

   setMatch_alnSet('n','w',def_NEqlW_alnDefs,alnSetST);
   setMatch_alnSet('n','s',def_NEqlS_alnDefs,alnSetST);
   setMatch_alnSet('n','m',def_NEqlM_alnDefs,alnSetST);
   setMatch_alnSet('n','k',def_NEqlK_alnDefs,alnSetST);
   setMatch_alnSet('n','r',def_NEqlR_alnDefs,alnSetST);
   setMatch_alnSet('n','y',def_NEqlY_alnDefs,alnSetST);
   setMatch_alnSet('n','b',def_NEqlB_alnDefs,alnSetST);
   setMatch_alnSet('n','d',def_NEqlD_alnDefs,alnSetST);
   setMatch_alnSet('n','h',def_NEqlH_alnDefs,alnSetST);
   setMatch_alnSet('n','v',def_NEqlV_alnDefs,alnSetST);
   setMatch_alnSet('n','n',def_NEqlN_alnDefs,alnSetST);
   setMatch_alnSet('n','x',def_NEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun14 Sec05 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setMatch_alnSet('x','a',def_XEqlA_alnDefs,alnSetST);
   setMatch_alnSet('x','c',def_XEqlC_alnDefs,alnSetST);
   setMatch_alnSet('x','g',def_XEqlG_alnDefs,alnSetST);
   setMatch_alnSet('x','t',def_XEqlT_alnDefs,alnSetST);
   setMatch_alnSet('x','u',def_XEqlT_alnDefs,alnSetST);

   setMatch_alnSet('x','w',def_XEqlW_alnDefs,alnSetST);
   setMatch_alnSet('x','s',def_XEqlS_alnDefs,alnSetST);
   setMatch_alnSet('x','m',def_XEqlM_alnDefs,alnSetST);
   setMatch_alnSet('x','k',def_XEqlK_alnDefs,alnSetST);
   setMatch_alnSet('x','r',def_XEqlR_alnDefs,alnSetST);
   setMatch_alnSet('x','y',def_XEqlY_alnDefs,alnSetST);
   setMatch_alnSet('x','b',def_XEqlB_alnDefs,alnSetST);
   setMatch_alnSet('x','d',def_XEqlD_alnDefs,alnSetST);
   setMatch_alnSet('x','h',def_XEqlH_alnDefs,alnSetST);
   setMatch_alnSet('x','v',def_XEqlV_alnDefs,alnSetST);
   setMatch_alnSet('x','n',def_XEqlN_alnDefs,alnSetST);
   setMatch_alnSet('x','x',def_XEqlX_alnDefs,alnSetST);

   return;
} /*init_alnSet*/

/*-------------------------------------------------------\
| Fun15: pDefMatchMatrix_alnSet
|   - print out the default match matrix
| Input:
|   - outFILE:
|     o file to print match matrix to
| Output:
|   - Prints:
|     o default match matrix in alnDefs.h to outFILE
\-------------------------------------------------------*/
void
pDefMatchMatrix_alnSet(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC:
   '   - print out the default match matrix
   '   o fun15 sec01:
   '     - print out header
   '   o fun15 sec02:
   '     - print out match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^  - print out header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   fprintf(
      (FILE *) outFILE,
      "// default match matrix%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "// legend:%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is match%s",
      def_ntEql_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is mismatch%s",
      def_ntNotEql_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is anonymous match%s",
      def_anonMatch_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i flag for anonymous; alone is mismatch%s",
      def_anonymous_alnDefs,
      str_endLine
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^  - print out match matrix
   ^  o fun15 sec02 sub01:
   ^    - a as first base
   ^  o fun15 sec02 sub02:
   ^    - t as first base
   ^  o fun15 sec02 sub03:
   ^    - u (t) as first base
   ^  o fun15 sec02 sub04:
   ^    - g as first base
   ^  o fun15 sec02 sub05:
   ^    - c as first base
   ^  o fun15 sec02 sub06:
   ^    - w (anonymous) as first base
   ^  o fun15 sec02 sub07:
   ^    - s (anonymous) as first base
   ^  o fun15 sec02 sub08:
   ^    - m (anonymous) as first base
   ^  o fun15 sec02 sub09:
   ^    - k (anonymous) as first base
   ^  o fun15 sec02 sub10:
   ^    - r (anonymous) as first base
   ^  o fun15 sec02 sub11:
   ^    - y (anonymous) as first base
   ^  o fun15 sec02 sub12:
   ^    - b (anonymous) as first base
   ^  o fun15 sec02 sub13:
   ^    - d (anonymous) as first base
   ^  o fun15 sec02 sub14:
   ^    - h (anonymous) as first base
   ^  o fun15 sec02 sub15:
   ^    - v (anonymous) as first base
   ^  o fun15 sec02 sub16:
   ^    - n (anonymous) as first base
   ^  o fun15 sec02 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec02 Sub01:
   *   - a as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "a a %i%s",
      def_AEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a t %i%s",
      def_AEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a u %i%s",
      def_AEqlU_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a g %i%s",
      def_AEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a c %i%s",
      def_AEqlC_alnDefs,
      str_endLine
   );

   /*anonymous matches*/
   fprintf(
      (FILE *) outFILE,
      "a w %i%s",
      def_AEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a s %i%s",
      def_AEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a m %i%s",
      def_AEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a k %i%s",
      def_AEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a r %i%s",
      def_AEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a y %i%s",
      def_AEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a b %i%s",
      def_AEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a d %i%s",
      def_AEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a h %i%s",
      def_AEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a v %i%s",
      def_AEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a n %i%s",
      def_AEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "a x %i%s",
      def_AEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub02:
   *   - t as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "t a %i%s",
      def_TEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t t %i%s",
      def_TEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t u %i%s",
      def_TEqlU_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t g %i%s",
      def_TEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t c %i%s",
      def_TEqlC_alnDefs,
      str_endLine
   );

   /*anonymous matches*/
   fprintf(
      (FILE *) outFILE,
      "t w %i%s",
      def_TEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t s %i%s",
      def_TEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t m %i%s",
      def_TEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t k %i%s",
      def_TEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t r %i%s",
      def_TEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t y %i%s",
      def_TEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t b %i%s",
      def_TEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t d %i%s",
      def_TEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t h %i%s",
      def_TEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t v %i%s",
      def_TEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t n %i%s",
      def_TEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "t x %i%s",
      def_TEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "u a %i%s",
      def_UEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u g %i%s",
      def_UEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u c %i%s",
      def_UEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u t %i%s",
      def_UEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u u %i%s",
      def_UEqlU_alnDefs,
      str_endLine
   );

   /*Set u & t to same scores (U is RNA version of T)*/
   fprintf(
      (FILE *) outFILE,
      "u w %i%s",
      def_UEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u s %i%s",
      def_UEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u m %i%s",
      def_UEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u k %i%s",
      def_UEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u r %i%s",
      def_UEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u y %i%s",
      def_UEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u b %i%s",
      def_UEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u d %i%s",
      def_UEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u h %i%s",
      def_UEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u v %i%s",
      def_UEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u n %i%s",
      def_UEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "u x %i%s",
      def_UEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub04:
   *   - g as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "g a %i%s",
      def_GEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g t %i%s",
      def_GEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g u %i%s",
      def_GEqlU_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g g %i%s",
      def_GEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g c %i%s",
      def_GEqlC_alnDefs,
      str_endLine
   );

   /*anonymous matches*/
   fprintf(
      (FILE *) outFILE,
      "g w %i%s",
      def_GEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g s %i%s",
      def_GEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g m %i%s",
      def_GEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g k %i%s",
      def_GEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g r %i%s",
      def_GEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g y %i%s",
      def_GEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g b %i%s",
      def_GEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g d %i%s",
      def_GEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g h %i%s",
      def_GEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g v %i%s",
      def_GEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g n %i%s",
      def_GEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "g x %i%s",
      def_GEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub05:
   *   - c as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "c a %i%s",
      def_CEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c t %i%s",
      def_CEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c u %i%s",
      def_CEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c g %i%s",
      def_CEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c c %i%s",
      def_CEqlC_alnDefs,
      str_endLine
   );

   /*anonymous matches*/
   fprintf(
      (FILE *) outFILE,
      "c w %i%s",
      def_CEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c s %i%s",
      def_CEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c m %i%s",
      def_CEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c k %i%s",
      def_CEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c r %i%s",
      def_CEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c y %i%s",
      def_CEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c b %i%s",
      def_CEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c d %i%s",
      def_CEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c h %i%s",
      def_CEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c v %i%s",
      def_CEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c n %i%s",
      def_CEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "c x %i%s",
      def_CEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "w a %i%s",
      def_WEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w c %i%s",
      def_WEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w g %i%s",
      def_WEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w t %i%s",
      def_WEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w u %i%s",
      def_WEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "w w %i%s",
      def_WEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w s %i%s",
      def_WEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w m %i%s",
      def_WEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w k %i%s",
      def_WEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w r %i%s",
      def_WEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w y %i%s",
      def_WEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w b %i%s",
      def_WEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w d %i%s",
      def_WEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w h %i%s",
      def_WEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w v %i%s",
      def_WEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w n %i%s",
      def_WEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "w x %i%s",
      def_WEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "s a %i%s",
      def_SEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s c %i%s",
      def_SEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s g %i%s",
      def_SEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s t %i%s",
      def_SEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s u %i%s",
      def_SEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "s w %i%s",
      def_SEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s s %i%s",
      def_SEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s m %i%s",
      def_SEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s k %i%s",
      def_SEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s r %i%s",
      def_SEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s y %i%s",
      def_SEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s b %i%s",
      def_SEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s d %i%s",
      def_SEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s h %i%s",
      def_SEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s v %i%s",
      def_SEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s n %i%s",
      def_SEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "s x %i%s",
      def_SEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "m a %i%s",
      def_MEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m c %i%s",
      def_MEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m g %i%s",
      def_MEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m t %i%s",
      def_MEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m u %i%s",
      def_MEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "m w %i%s",
      def_MEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m s %i%s",
      def_MEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m m %i%s",
      def_MEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m k %i%s",
      def_MEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m r %i%s",
      def_MEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m y %i%s",
      def_MEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m b %i%s",
      def_MEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m d %i%s",
      def_MEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m h %i%s",
      def_MEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m v %i%s",
      def_MEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m n %i%s",
      def_MEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "m x %i%s",
      def_MEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "k a %i%s",
      def_KEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k c %i%s",
      def_KEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k g %i%s",
      def_KEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k t %i%s",
      def_KEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k u %i%s",
      def_KEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "k w %i%s",
      def_KEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k s %i%s",
      def_KEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k m %i%s",
      def_KEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k k %i%s",
      def_KEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k r %i%s",
      def_KEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k y %i%s",
      def_KEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k b %i%s",
      def_KEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k d %i%s",
      def_KEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k h %i%s",
      def_KEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k v %i%s",
      def_KEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k n %i%s",
      def_KEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "k x %i%s",
      def_KEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "r a %i%s",
      def_REqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r c %i%s",
      def_REqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r g %i%s",
      def_REqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r t %i%s",
      def_REqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r u %i%s",
      def_REqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "r w %i%s",
      def_REqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r s %i%s",
      def_REqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r m %i%s",
      def_REqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r k %i%s",
      def_REqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r r %i%s",
      def_REqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r y %i%s",
      def_REqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r b %i%s",
      def_REqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r d %i%s",
      def_REqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r h %i%s",
      def_REqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r v %i%s",
      def_REqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r n %i%s",
      def_REqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "r x %i%s",
      def_REqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "y a %i%s",
      def_YEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y c %i%s",
      def_YEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y g %i%s",
      def_YEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y t %i%s",
      def_YEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y u %i%s",
      def_YEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "y w %i%s",
      def_YEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y s %i%s",
      def_YEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y m %i%s",
      def_YEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y k %i%s",
      def_YEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y r %i%s",
      def_YEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y y %i%s",
      def_YEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y b %i%s",
      def_YEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y d %i%s",
      def_YEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y h %i%s",
      def_YEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y v %i%s",
      def_YEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y n %i%s",
      def_YEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "y x %i%s",
      def_YEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "b a %i%s",
      def_BEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b c %i%s",
      def_BEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b g %i%s",
      def_BEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b t %i%s",
      def_BEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b u %i%s",
      def_BEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "b w %i%s",
      def_BEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b s %i%s",
      def_BEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b m %i%s",
      def_BEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b k %i%s",
      def_BEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b r %i%s",
      def_BEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b y %i%s",
      def_BEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b b %i%s",
      def_BEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b d %i%s",
      def_BEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b h %i%s",
      def_BEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b v %i%s",
      def_BEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b n %i%s",
      def_BEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "b x %i%s",
      def_BEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "d a %i%s",
      def_DEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d c %i%s",
      def_DEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d g %i%s",
      def_DEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d t %i%s",
      def_DEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d u %i%s",
      def_DEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "d w %i%s",
      def_DEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d s %i%s",
      def_DEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d m %i%s",
      def_DEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d k %i%s",
      def_DEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d r %i%s",
      def_DEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d y %i%s",
      def_DEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d b %i%s",
      def_DEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d d %i%s",
      def_DEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d h %i%s",
      def_DEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d v %i%s",
      def_DEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d n %i%s",
      def_DEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "d x %i%s",
      def_DEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "h a %i%s",
      def_HEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h c %i%s",
      def_HEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h g %i%s",
      def_HEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h t %i%s",
      def_HEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h u %i%s",
      def_HEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "h w %i%s",
      def_HEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h s %i%s",
      def_HEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h m %i%s",
      def_HEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h k %i%s",
      def_HEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h r %i%s",
      def_HEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h y %i%s",
      def_HEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h b %i%s",
      def_HEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h d %i%s",
      def_HEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h h %i%s",
      def_HEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h v %i%s",
      def_HEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h n %i%s",
      def_HEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "h x %i%s",
      def_HEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "v a %i%s",
      def_VEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v c %i%s",
      def_VEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v g %i%s",
      def_VEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v t %i%s",
      def_VEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v u %i%s",
      def_VEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "v w %i%s",
      def_VEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v s %i%s",
      def_VEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v m %i%s",
      def_VEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v k %i%s",
      def_VEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v r %i%s",
      def_VEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v y %i%s",
      def_VEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v b %i%s",
      def_VEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v d %i%s",
      def_VEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v h %i%s",
      def_VEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v v %i%s",
      def_VEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v n %i%s",
      def_VEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "v x %i%s",
      def_VEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "n a %i%s",
      def_NEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n c %i%s",
      def_NEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n g %i%s",
      def_NEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n t %i%s",
      def_NEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n u %i%s",
      def_NEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "n w %i%s",
      def_NEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n s %i%s",
      def_NEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n m %i%s",
      def_NEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n k %i%s",
      def_NEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n r %i%s",
      def_NEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n y %i%s",
      def_NEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n b %i%s",
      def_NEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n d %i%s",
      def_NEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n h %i%s",
      def_NEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n v %i%s",
      def_NEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n n %i%s",
      def_NEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "n x %i%s",
      def_NEqlX_alnDefs,
      str_endLine
   );

   /*****************************************************\
   * Fun15 Sec02 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "x a %i%s",
      def_XEqlA_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x c %i%s",
      def_XEqlC_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x g %i%s",
      def_XEqlG_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x t %i%s",
      def_XEqlT_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x u %i%s",
      def_XEqlT_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "x w %i%s",
      def_XEqlW_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x s %i%s",
      def_XEqlS_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x m %i%s",
      def_XEqlM_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x k %i%s",
      def_XEqlK_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x r %i%s",
      def_XEqlR_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x y %i%s",
      def_XEqlY_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x b %i%s",
      def_XEqlB_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x d %i%s",
      def_XEqlD_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x h %i%s",
      def_XEqlH_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x v %i%s",
      def_XEqlV_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x n %i%s",
      def_XEqlN_alnDefs,
      str_endLine
   );
   fprintf(
      (FILE *) outFILE,
      "x x %i%s",
      def_XEqlX_alnDefs,
      str_endLine
   );
} /*pDefMatchMatrix_alnSet*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However,for cases were the public domain is not
:   suitable,such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible,this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy,modify,publish,use,compile,
:   sell,or distribute this software,either in source
:   code form or as a compiled binary,for any purpose,
:   commercial or non-commercial,and by any means.
: 
: In jurisdictions that recognize copyright laws,the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY,WHETHER IN AN ACTION OF
:   CONTRACT,TORT OR OTHERWISE,ARISING FROM,OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information,please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted,free of charge,to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"),to
:   deal in the Software without restriction,including
:   without limitation the rights to use,copy,modify,
:   merge,publish,distribute,sublicense,and/or sell
:   copies of the Software,and to permit persons to whom
:   the Software is furnished to do so,subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY,WHETHER IN
:   AN ACTION OF CONTRACT,TORT OR OTHERWISE,ARISING
:   FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
