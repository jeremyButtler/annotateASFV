/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' dirMatrix SOF: Start Of File
'   - holds functions for dealing with the dirMatrix
'     returned by water and needle
'   o header:
'     - included libraries
'   o .h st01: alnMatrixStruct
'     - Holds the direction matrix and best score(s) for a
'       single aligment
'   o fun01: blank_dirMatrix
'     - blanks all score, index, and length variables in a
'       dirMatrix structure
'   o fun02: init_dirMatrix
'     - initializes a dirMatrix structure
'   o fun03: freeStack_dirMatrix
'     - frees heap allocated variables in a dirMatrix
'       structure
'   o fun04: freeHeap_dirMatrix
'     - frees a dirMatrix structure
'   o fun05: getAln_dirMatrix
'     - gets a sam file entry (alignment) from a direction
'       matrix (inside the dirMatrix structure)
'   o fun06: getCig_dirMatrix
'     - gets a cigar for an alignment
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

#include "dirMatrix.h"

#include "../genLib/ulCp.h"
#include "../genLib/charCp.h"
#include "../genLib/numToStr.h"

#include "../genBio/samEntry.h"
#include "../genBio/seqST.h"

#include "alnSet.h"
#include "indexToCoord.h"

/*.h files only*/
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .c  #include "../genLib/base10str.h"
!   o .c  #include "../genLib/fileFun.h"
!   o .h  #include "../genLib/endLine.h"
!   o .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_dirMatrix
|   - blanks all score, index, and length variables in a
|     dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix to blank
| Output:
|   - Sets:
|     o all score, length, and index variables are set to
|       0 (direction matrix is not touched)
\-------------------------------------------------------*/
void
blank_dirMatrix(
  struct dirMatrix *dirMatrixSTPtr
){
   dirMatrixSTPtr->scoreSL = 0;
   dirMatrixSTPtr->indexSL = 0;

   dirMatrixSTPtr->refLenSL = 0;
   dirMatrixSTPtr->refOffsetSL = 0;
   dirMatrixSTPtr->refEndSL = 0;

   dirMatrixSTPtr->qryLenSL = 0;
   dirMatrixSTPtr->qryOffsetSL = 0;
   dirMatrixSTPtr->qryEndSL = 0;

   dirMatrixSTPtr->errSC = 0;
} /*blank_dirMatrixSC*/

/*-------------------------------------------------------\
| Fun02: init_dirMatrix
|   - initializes a dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix structure initialize
| Output:
|   - Sets:
|     o all variables in matrixSTPtr to 0
\-------------------------------------------------------*/
void
init_dirMatrix(
  struct dirMatrix *dirMatrixSTPtr
){
   blank_dirMatrix(dirMatrixSTPtr);
  
   dirMatrixSTPtr->dirMatrixSC = 0;
   dirMatrixSTPtr->lenMatrixSL = 0;
  
   dirMatrixSTPtr->scoreArySL = 0;
   dirMatrixSTPtr->lenScoreSL = 0;
} /*init_dirMatrixSC*/

/*-------------------------------------------------------\
| Fun03: freeStack_dirMatrix
|   - frees heap allocated variables in a dirMatrix
|     structure
| Input:
|   - matrixSTPtr
|     o pointer to dirMatrix structure with variables to
|       free
| Output:
|   - Frees:
|     o dirMatrix->dirMatrixSC
|   - Sets:
|     o all non-freeded variables to 0
\-------------------------------------------------------*/
void
freeStack_dirMatrix(
   struct dirMatrix *matrixSTPtr
){
   if(matrixSTPtr)
   { /*If: not null*/
      free(matrixSTPtr->dirMatrixSC);
      free(matrixSTPtr->scoreArySL);

      init_dirMatrix(matrixSTPtr); /*sets to null*/
   } /*If: not null*/

   return;
} /*freeAlnMatrixStack*/

/*-------------------------------------------------------\
| Fun04: freeHeap_dirMatrix
|   - frees a dirMatrix structure
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to free
| Output:
|   - Frees:
|     o matrixSTPtr
\-------------------------------------------------------*/
void
freeHeap_dirMatrix(
   struct dirMatrix *matrixSTPtr
){
   freeStack_dirMatrix(matrixSTPtr);
   free(matrixSTPtr);
} /*freeAlnMatrix*/

/*-------------------------------------------------------\
| Fun05: getAln_dirMatrix
|   - gets a sam file entry (alignment) from a direction
|     matrix (inside the dirMatrix structure)
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to get alignment
|       from
|   - indexSL:
|     o index of last base in the alignment
|     o 0 to use index from matirxSTPtr
|   - revBl:
|     o 1: reverse alignment (sam flag is 16)
|       - this means I had to reverse complement the
|         reference sequence
|     o 0: foward alignment (sam flag is 0)
|   - qrySTPtr:
|     o pointer to a seqST with the query sequence
|   - refSTPtr:
|     o pointer to a seqST with the reference sequence
|   - samSTPtr:
|     o pointer to a samEntry struct to hold the alignment
|   - numAnonUI:
|     o pointer to unsigned in to hold the number of
|       anonymous bases (matches only)
|   - alnSetSTPtr:
|     o pointer to alnSet structure with the match matrix
| Output:
|   - Modifies:
|     o samSTPtr to have the aligned sequence
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error (only error possible)
\-------------------------------------------------------*/
signed char
getAln_dirMatrix(
   struct dirMatrix *matrixSTPtr,
   signed long indexSL,
   signed char revBl,
   struct seqST *qrySTPtr,
   struct seqST *refSTPtr,
   struct samEntry *samSTPtr,
   unsigned int *numAnonUI,
   struct alnSet *alnSetSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - gets an alignment from a dirMatrix structure
   '   o fun05 sec01:
   '     - variable declerations
   '   o fun05 sec02:
   '     - find start and ending positions
   '   o fun05 sec03:
   '     - allocate memroy and copy query
   '   o fun05 sec04:
   '     - get alignment form matrix
   '   o fun05 sec05:
   '     - add starting softmasked bases and invert cigar
   '   o fun05 sec06:
   '     - add tags (NM, AS, nn)
   '   o fun05 sec07:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define expandCig_fun05 256

   unsigned int qryPosUI = 0;
   unsigned int refPosUI = 0;
   unsigned int lenRefUI = matrixSTPtr->refLenSL;

   signed char *qrySeqStr =
               (signed char *) qrySTPtr->seqStr;
   signed char *refSeqStr =
               (signed char *) refSTPtr->seqStr;
 
   signed char *tmpStr = 0;
   unsigned char tmpUC = 0;
   signed char matchBl = 0; /*check if had match or snp*/

   signed char *dirMatrixSC = matrixSTPtr->dirMatrixSC;

   struct seqST seqDoNotFreeST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - find sequence start and ending positions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numAnonUI = 0;

   if(! indexSL)
      indexSL = matrixSTPtr->indexSL; /*no index input*/

   qryPosUI =
      (unsigned int)
      qryCoord_indexToCoord(
         (unsigned long) lenRefUI,
         indexSL
      );

   refPosUI =
      (unsigned int)
      refCoord_indexToCoord(
         (unsigned long) lenRefUI,
         indexSL
      );

   /*bases not aligned by user*/
   refSeqStr += matrixSTPtr->refOffsetSL;
   qrySeqStr += matrixSTPtr->qryOffsetSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - allocate memroy and copy query
   ^   o fun05 sec03 sub01:
   ^     - add ending soft masked bases to cigar
   ^   o fun05 sec03 sub02:
   ^     - copy read id
   ^   o fun05 sec03 sub03:
   ^     - copy the reference id
   ^   o fun05 sec03 sub04:
   ^     - copy query sequence
   ^   o fun05 sec03 sub05:
   ^     - copy query q-score entry
   ^   o fun05 sec03 sub06:
   ^     - set flag and reference end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec03 Sub01:
   *   - add ending soft masked bases to cigar
   \*****************************************************/

   blank_samEntry(samSTPtr);

   if(qryPosUI + 1 < qrySTPtr->seqLenSL)
   { /*If: need to add ending softmasked bases*/
      samSTPtr->cigTypeStr[0] = 'S';

      samSTPtr->cigArySI[0] =
         qrySTPtr->seqLenSL - qryPosUI - 1;
         /*-1 to account for seqLenSL being index 1 and
         `   qryPosUI being index 0
         */
      samSTPtr->maskCntUI =
         qrySTPtr->seqLenSL - qryPosUI - 1;
   } /*If: need to add ending softmasked bases*/

   /*else is set to null*/

   /*****************************************************\
   * Fun05 Sec03 Sub02:
   *   - copy read id
   \*****************************************************/

   /*check for header marks*/
   tmpUC = (qrySTPtr->idStr[0] == '>');
   tmpUC += (qrySTPtr->idStr[0] == '@');

   tmpStr = (signed char *) qrySTPtr->idStr + tmpUC;

   while(*tmpStr++ > 32) ;

   --tmpStr;
   matchBl = *tmpStr;
   *tmpStr = '\0';

   samSTPtr->qryIdLenUC = 
      cpStr_ulCp(
         (signed char *) samSTPtr->qryIdStr, 
         (signed char *) qrySTPtr->idStr + tmpUC
      ); /*copy read id*/

   *tmpStr = matchBl;

   /*****************************************************\
   * Fun05 Sec03 Sub03:
   *   - copy the reference id
   \*****************************************************/

   /*check for header marks*/
   tmpUC = (refSTPtr->idStr[0] == '>');
   tmpUC += (refSTPtr->idStr[0] == '@');

   tmpStr = (signed char *) refSTPtr->idStr;

   while(*tmpStr++ > 32) ;

   --tmpStr;
   matchBl = *tmpStr;
   *tmpStr = '\0';

   samSTPtr->refIdLenUC = 
      cpStr_ulCp(
         (signed char *) samSTPtr->refIdStr, 
         (signed char *) refSTPtr->idStr + tmpUC
      ); /*copy read id*/

   *tmpStr = matchBl;

   /*****************************************************\
   * Fun05 Sec03 Sub04:
   *   - copy query sequence
   \*****************************************************/

   if(samSTPtr->seqSizeUI < qrySTPtr->seqLenSL)
   { /*If: I need more memory for the sequence*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         (signed char *)
         malloc(
              (qrySTPtr->seqLenSL+1024)
            * sizeof(signed char)
         );

      if(! samSTPtr->seqStr)
         goto memErr_fun05_sec07;

      samSTPtr->seqSizeUI = qrySTPtr->seqLenSL + 1023;
      tmpStr = 0;
   } /*If: I need more memory for the sequence*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      qrySTPtr->seqStr,
      qrySTPtr->seqLenSL
   );

   samSTPtr->readLenUI = qrySTPtr->seqLenSL;

   /*convert the index to a real sequence*/
   indexToSeq_alnSet(samSTPtr->seqStr);

   /*****************************************************\
   * Fun05 Sec03 Sub05:
   *   - copy query q-score entry
   \*****************************************************/

   if(qrySTPtr->qLenSL)
   { /*If: have q-score entry*/
      if(samSTPtr->qSizeUI < qrySTPtr->qLenSL)
      { /*If: I need more memory for the sequence*/
         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;

         samSTPtr->qStr =
            (signed char *)
            malloc(
                 (qrySTPtr->qLenSL + 1024)
               * sizeof(signed char)
            );

         if(! samSTPtr->qStr)
            goto memErr_fun05_sec07;

         samSTPtr->qSizeUI = qrySTPtr->qLenSL + 1023;
      } /*If: I need more memory for the sequence*/

      cpQEntry_samEntry(
         samSTPtr,
         (signed char *) qrySTPtr->qStr,
         1 /*need to make sure histogram is blanked*/
      );

      samSTPtr->readLenUI = qrySTPtr->seqLenSL;
   } /*If: have q-score entry*/

   else
   { /*Else: set the q-score entry to nothing*/
      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
   } /*Else: set the q-score entry to nothing*/

   /*****************************************************\
   * Fun05 Sec03 Sub06:
   *   - set flag and reference end
   \*****************************************************/

   if(revBl)
   { /*If: read is reverse complment forward*/
      samSTPtr->flagUS = 16;
      samSTPtr->refStartUI = lenRefUI - refPosUI - 1;
   } /*If: read is reverse complment forward*/

   else
   { /*Else: read is forward*/
      samSTPtr->flagUS = 0;
      samSTPtr->refEndUI = refPosUI;
   } /*Else: read is forward*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - get alignment form matrix
   ^   o fun05 sec04 sub01:
   ^     - check cig memory + start loop cigar fill loop
   ^   o fun05 sec04 sub02:
   ^     - insertion cases + cigar switch statment
   ^   o fun05 sec04 sub03:
   ^     - snp or match cases
   ^   o fun05 sec04 sub04:
   ^     - deletion cases
   ^   o fun05 sec04 sub05:
   ^     - note ending referencde position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec04 Sub01:
   *   - check cig memory + start loop cigar fill loop
   \*****************************************************/

   dirMatrixSC += indexSL;

   while(*dirMatrixSC != def_mvStop_alnDefs)
   { /*Loop: trace alignment path*/

      /*I do not expect sequence long cigars, so resizing
      `   in small chunks is prefered, expecially since
      `   I often keep the same samEntry struct, so
      `   one resize effects the entire program
      */
      if(
         samSTPtr->cigLenUI + 3 >= samSTPtr->cigSizeUI
      ){ /*If: want more cigar memory*/
         tmpStr =
            realloc(
               samSTPtr->cigTypeStr,
               (samSTPtr->cigLenUI + expandCig_fun05)
                  * sizeof(signed char)
            );
         if(! tmpStr)
            goto memErr_fun05_sec07;
         samSTPtr->cigTypeStr = tmpStr;

         tmpStr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (samSTPtr->cigLenUI + expandCig_fun05)
                  * sizeof(signed int)
            );
         if(! tmpStr)
            goto memErr_fun05_sec07;
         samSTPtr->cigArySI = (signed int *) tmpStr;

         samSTPtr->cigSizeUI = samSTPtr->cigLenUI;
         samSTPtr->cigSizeUI += (expandCig_fun05 - 8);
      } /*If: want more cigar memory*/

      /**************************************************\
      * Fun05 Sec04 Sub02:
      *   - insertion cases + cigar switch statment
      \**************************************************/

      switch(*dirMatrixSC)
      { /*Switch: find direction (snp/del/ins)*/
         case def_mvStop_alnDefs:
            break;

         case def_mvIns_alnDefs:
         /*Case: is a insertion*/
            if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               != 'I'
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               ) ++samSTPtr->cigLenUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->cigLenUI
               ] = 'I';

               samSTPtr->cigArySI[
                  samSTPtr->cigLenUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigArySI[samSTPtr->cigLenUI];

            ++samSTPtr->insCntUI;

            --qryPosUI;
            dirMatrixSC -= (lenRefUI + 1);
            break;
         /*Case: is a insertion*/

         /***********************************************\
         * Fun05 Sec04 Sub03:
         *   - snp or match cases
         \***********************************************/

         case def_mvSnp_alnDefs:
         /*Case: is a snp or match*/
            matchBl =
               getMatch_alnSet(
                  qrySeqStr[qryPosUI],
                  refSeqStr[refPosUI],
                  alnSetSTPtr
               ); /*find if had a match*/

            if(matchBl & def_ntEql_alnDefs)
            { /*If: had a match*/
               /*count number of anonymous matches*/
               *numAnonUI +=
                  !!(matchBl & def_anonymous_alnDefs);

               matchBl = '='; /*match*/
               ++samSTPtr->matchCntUI;

            } /*If: had a match*/

            else
            { /*Else: had a mismatch*/
               matchBl = 'X'; /*mismatch*/
               ++samSTPtr->snpCntUI;
            } /*Else: had a mismatch*/

            if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               != matchBl
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               ) ++samSTPtr->cigLenUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->cigLenUI
               ] = matchBl;

               samSTPtr->cigArySI[
                  samSTPtr->cigLenUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigArySI[samSTPtr->cigLenUI];

            --refPosUI;
            --qryPosUI;

            ++samSTPtr->alnReadLenUI;

            dirMatrixSC -= (lenRefUI + 2);
            break;
         /*Case: is a snp or match*/

         /***********************************************\
         * Fun05 Sec04 Sub04:
         *   - deletion cases
         \***********************************************/

         case def_mvDel_alnDefs:
         /*Case: is a deletion*/
            if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               != 'D'
            ){ /*If: need to make new cigar entry*/
               if(
                  samSTPtr->cigTypeStr[samSTPtr->cigLenUI]
               ) ++samSTPtr->cigLenUI; /*if not null*/

               samSTPtr->cigTypeStr[
                  samSTPtr->cigLenUI
               ] = 'D';

               samSTPtr->cigArySI[
                  samSTPtr->cigLenUI
               ] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++samSTPtr->cigArySI[samSTPtr->cigLenUI];

            ++samSTPtr->delCntUI;
            ++samSTPtr->alnReadLenUI;

            --refPosUI;
            --dirMatrixSC;
            break;
         /*Case: is a deletion*/
      } /*Switch: find direction (snp/del/ins)*/
   } /*Loop: trace alignment path*/

   /*****************************************************\
   * Fun05 Sec04 Sub06:
   *   - note ending reference position
   \*****************************************************/

   /*account for ending on one index behind*/
   ++refPosUI;
   ++qryPosUI;

   if(revBl)
      samSTPtr->refEndUI = lenRefUI - refPosUI - 1;
   else
      samSTPtr->refStartUI = refPosUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec05:
   ^   - add starting softmasked bases and invert cigar
   ^   o fun05 sec05 sub01:
   ^     - add starting softmasked bases to cigar
   ^   o fun05 sec05 sub02:
   ^     - if forward sequence invert cigar; is backwards
   ^   o fun05 sec05 sub03:
   ^     - if reverse sequence; reverse complement
   ^   o fun05 sec05 sub04:
   ^     - make sure cigar length is index 1 (not 0)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec05 Sub01:
   *   - add starting softmasked bases to cigar
   \*****************************************************/

   if(qryPosUI > 0)
   { /*If: need to add starting softmasked bases*/
      ++samSTPtr->cigLenUI;
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = 'S';

      samSTPtr->cigArySI[samSTPtr->cigLenUI] = qryPosUI;
      samSTPtr->maskCntUI += qryPosUI;
   } /*If: need to add starting softmasked bases*/

   /*****************************************************\
   * Fun05 Sec05 Sub02:
   *   - if forward sequence invert cigar; is backwards
   \*****************************************************/

   if(! revBl)
   { /*If: need to reverse (ivert) cigar*/
      qryPosUI = 0;
      refPosUI = samSTPtr->cigLenUI;

      while(qryPosUI < refPosUI)
      { /*Loop: invert the cigar to forward direction*/
         samSTPtr->cigTypeStr[qryPosUI] ^=
            samSTPtr->cigTypeStr[refPosUI];
         samSTPtr->cigTypeStr[refPosUI] ^=
            samSTPtr->cigTypeStr[qryPosUI];
         samSTPtr->cigTypeStr[qryPosUI] ^=
            samSTPtr->cigTypeStr[refPosUI];

         samSTPtr->cigArySI[qryPosUI] ^=
            samSTPtr->cigArySI[refPosUI];
         samSTPtr->cigArySI[refPosUI] ^=
            samSTPtr->cigArySI[qryPosUI];
         samSTPtr->cigArySI[qryPosUI] ^=
            samSTPtr->cigArySI[refPosUI];

         ++qryPosUI;
         --refPosUI;
      } /*Loop: invert the cigar to forward direction*/

      /*null/0 or softmask at start flipped to end*/
   } /*If: need to reverse (ivert) cigar*/

   /*****************************************************\
   * Fun05 Sec05 Sub03:
   *   - if reverse sequence; reverse complement
   \*****************************************************/

   else
   { /*Else: reverse complement sequence*/
      seqDoNotFreeST.seqStr = samSTPtr->seqStr;
      seqDoNotFreeST.seqLenSL = samSTPtr->readLenUI;

      if(
            samSTPtr->qStr[0] != '*'
         && samSTPtr->qStr[1] != '\0'
         && samSTPtr->qStr[0] != '\0'
      ){ /*If: have q-score entry*/
         seqDoNotFreeST.qStr = samSTPtr->qStr;
         seqDoNotFreeST.qLenSL = samSTPtr->readLenUI;
      } /*If: have q-score entry*/

      else
      { /*Else: do not have q-score entry*/
         seqDoNotFreeST.qStr = 0;
         seqDoNotFreeST.qLenSL = 0;
      } /*Else: do not have q-score entry*/

      revComp_seqST(&seqDoNotFreeST);

      seqDoNotFreeST.seqStr = 0;
      seqDoNotFreeST.seqLenSL = 0;

      seqDoNotFreeST.qStr = 0;
      seqDoNotFreeST.qLenSL = 0;
   } /*Else: reverse complement sequence*/


   /*****************************************************\
   * Fun05 Sec05 Sub04:
   *   - make sure cigar length is index 1 (not 0)
   \*****************************************************/

   if(samSTPtr->cigTypeStr[samSTPtr->cigLenUI] != '\0')
   { /*If: is index 0*/
      ++samSTPtr->cigLenUI; /*convert to index 1*/
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
   } /*If: is index 0*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec06:
   ^   - add tags (NM, AS, nn)
   ^   o fun05 sec06 sub01:
   ^     - add NM (edit distance/number differences) flag
   ^   o fun05 sec06 sub02:
   ^     - add AS (score) flag
   ^   o fun05 sec06 sub03:
   ^     - add nn (anonymoys bases) flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec06 Sub01:
   *   - add NM (edit distance/number differences) flag
   \*****************************************************/

   tmpStr = (signed char *) samSTPtr->extraStr;

   samSTPtr->extraLenUI +=
      cpDelim_charCp(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         (signed char *) "NM:i:",
         '\0'
      ); /*copy edit distance*/

   samSTPtr->extraLenUI +=
      numToStr(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         samSTPtr->snpCntUI
            + samSTPtr->delCntUI
            + samSTPtr->insCntUI
      ); /*find edit distance*/

   /*****************************************************\
   * Fun05 Sec06 Sub02:
   *   - add AS (score) flag
   \*****************************************************/

   tmpStr[samSTPtr->extraLenUI++] = '\t';

   samSTPtr->extraLenUI +=
      cpDelim_charCp(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         (signed char *) "AS:i:",
          '\0'
      ); /*copy assembler score*/

   samSTPtr->extraLenUI +=
      numToStr(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         matrixSTPtr->scoreSL /  def_scoreAdj_alnDefs
      ); /*get alignment score*/

   /*****************************************************\
   * Fun05 Sec06 Sub03:
   *   - add nn (anonymoys bases) flag
   \*****************************************************/

   tmpStr[samSTPtr->extraLenUI++] = '\t';

   samSTPtr->extraLenUI +=
      cpDelim_charCp(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         (signed char *) "nn:i:",
         '\0'
      ); /*copy number anonymous bases*/

   samSTPtr->extraLenUI +=
      numToStr(
         (signed char *) &tmpStr[samSTPtr->extraLenUI],
         *numAnonUI
      ); /*get number anonymous bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec07:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun05_sec07:;

   return 1;
} /*getAln_dirMatrix*/

/*-------------------------------------------------------\
| Fun06: getCig_dirMatrix
|   - gets a cigar for an alignment
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to get alignment
|       from
|   - indexSL:
|     o index of last base in the alignment
|     o 0 to use index from matirxSTPtr
|   - revBl:
|     o 1: reverse alignment (sam flag is 16)
|       - this means I had to reverse complement the
|         reference sequence
|     o 0: foward alignment (sam flag is 0)
|   - qrySTPtr:
|     o pointer to a seqST with the query sequence
|   - refSTPtr:
|     o pointer to a seqST with the reference sequence
|   - cigTypeStr:
|     o pointer to c-string with cigar entry type array
|   - cigArySI:
|     o pointer to signed int array with the length of
|       each cigar entry
|   - cigPosUI:
|     o position at in cigar
|   - cigLenUI:
|     o pointer to unsigned int with length of the cigar
|       buffer
|   - refStartUI:
|     o unsigned int pointer to point to first reference
|       base in cigar
|   - numAnonUI:
|     o pointer to unsigned in to hold the number of
|       anonymous bases (matches only)
|   - maskCntUI:
|     o pointer to unsigned long to hold number of
|       maksed bases
|   - delAtEndBl:
|     o 1: add deltions if reference is short at end
|     o 0: ignore
|   - alnSetSTPtr:
|     o pointer to alnSet structure with the match matrix
| Output:
|   - Modifies:
|     o cigTypeStr to have the cigar entry types
|     o cigArySI to have the length of each cigar entry
|     o cigLenUI if cigTypeStr and cigArySI needed to be
|       resized
|     o refStartUI to have first reference base in cigar
|     o numAnonUI to have number of matching anonymous
|       bases
|     o maskCntUI to have number of masked bases
|   - Returns:
|     o new position in cigar
|     o -1 for memory error (only error possible)
\-------------------------------------------------------*/
signed long
getCig_dirMatrix(
   struct dirMatrix *matrixSTPtr,
   signed long indexSL,
   signed char revBl,
   struct seqST *qrySTPtr,
   struct seqST *refSTPtr,
   signed char **cigTypeStr,
   signed int **cigArySI,
   unsigned int cigPosUI,
   unsigned int *cigLenUI,
   unsigned int *refStartUI,
   unsigned int *numAnonUI,
   unsigned int *maskCntUI,
   signed char delAtEndBl,
   struct alnSet *alnSetSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC:
   '   - gets an alignment from a dirMatrix structure
   '   o fun06 sec01:
   '     - variable declerations
   '   o fun06 sec02:
   '     - find start and ending positions
   '   o fun06 sec03:
   '     - add ending soft masked bases to cigar
   '   o fun06 sec04:
   '     - get alignment form matrix
   '   o fun06 sec05:
   '     - add starting softmasked bases and invert cigar
   '   o fun06 sec06:
   '     - add tags (NM, AS, nn)
   '   o fun06 sec07:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define expandCig_fun06 256

   unsigned int qryPosUI = 0;
   unsigned int refPosUI = 0;
   unsigned int lenRefUI = matrixSTPtr->refLenSL;
   unsigned int startUI = cigPosUI;

   signed char *qrySeqStr =
      (signed char *) qrySTPtr->seqStr;
   signed char *refSeqStr =
      (signed char *) refSTPtr->seqStr;
 
   signed char *tmpStr = 0;
   signed char matchBl = 0; /*check if had match or snp*/

   signed char *dirMatrixSC = matrixSTPtr->dirMatrixSC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - find sequence start and ending positions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(*refStartUI < matrixSTPtr->refOffsetSL)
      *refStartUI += matrixSTPtr->refOffsetSL;
      /*likely user did not account for offset*/

   *numAnonUI = 0;

   if(! indexSL)
      indexSL = matrixSTPtr->indexSL; /*no index input*/

   qryPosUI =
      (unsigned int)
      qryCoord_indexToCoord(
         (unsigned long) lenRefUI,
         indexSL
      );

   refPosUI =
      (unsigned int)
      refCoord_indexToCoord(
         (unsigned long) lenRefUI,
         indexSL
      );

   /*bases not aligned by user*/
   refSeqStr += matrixSTPtr->refOffsetSL;
   qrySeqStr += matrixSTPtr->qryOffsetSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^   - add ending soft masked bases to cigar
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! *cigTypeStr || ! *cigArySI)
   { /*If: need memory*/
      if(*cigTypeStr)
         free(*cigTypeStr);
      *cigTypeStr = 0;

      *cigTypeStr =malloc((256 +8) * sizeof(signed char));

      if(! *cigTypeStr)
         goto memErr_fun06_sec07;


      if(*cigArySI)
         free(*cigArySI);
      *cigArySI = 0;

      *cigArySI = malloc((256 + 8) * sizeof(signed int));
      if(! *cigArySI)
         goto memErr_fun06_sec07;

      *cigLenUI = 256;

      *cigTypeStr[0] = 0;
      *cigArySI[0] = 0;
   } /*If: need memory*/

   if(
         ! revBl
      && (*cigTypeStr)[cigPosUI] > 0
   ){ /*If: cigar entry already here*/
       ++cigPosUI;
       ++startUI;
       /*forward sequence mean backwards cigars*/
   } /*If: cigar entry already here*/

   if( refPosUI < (matrixSTPtr->refLenSL - 1) )
   { /*If: missing bases at end*/
      if(revBl)
      { /*If: reverse complement sequence*/
         /*find start offset*/
         *refStartUI = matrixSTPtr->refEndSL - refPosUI;
         *refStartUI += matrixSTPtr->refOffsetSL;
      } /*If: reverse complement sequence*/

      else if(delAtEndBl)
      { /*Else: foward sequence (deletions at end)*/
         if((*cigTypeStr)[cigPosUI] > 32)
            ++cigPosUI;

         (*cigTypeStr)[cigPosUI] = 'D';

         (*cigArySI)[cigPosUI] =
            matrixSTPtr->refEndSL
               - matrixSTPtr->refOffsetSL
               - refPosUI;
      } /*Else: foward sequence (deletions at end)*/
   } /*If: missing bases at end*/

   if(qryPosUI + 1 < matrixSTPtr->qryLenSL)
   { /*If: need to add ending softmasked bases*/
      if((*cigTypeStr)[cigPosUI] > 32)
         ++cigPosUI;

      (*cigTypeStr)[cigPosUI] = 'S';
      (*cigArySI)[cigPosUI] =
         matrixSTPtr->qryLenSL - qryPosUI;

      *maskCntUI += matrixSTPtr->qryLenSL - qryPosUI;
   } /*If: need to add ending softmasked bases*/

   /*else is set to null*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec04:
   ^   - get alignment form matrix
   ^   o fun06 sec04 sub01:
   ^     - check cig memory + start loop cigar fill loop
   ^   o fun06 sec04 sub02:
   ^     - insertion cases + cigar switch statment
   ^   o fun06 sec04 sub03:
   ^     - snp or match cases
   ^   o fun06 sec04 sub04:
   ^     - deletion cases
   ^   o fun06 sec04 sub05:
   ^     - move to first base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec04 Sub01:
   *   - check cig memory + start loop cigar fill loop
   \*****************************************************/

   dirMatrixSC += indexSL;

   while(*dirMatrixSC != def_mvStop_alnDefs)
   { /*Loop: trace alignment path*/

      /*I do not expect sequence long cigars, so resizing
      `   in small chunks is prefered, expecially since
      `   I often keep the same samEntry struct, so
      `   one resize effects the entire program
      */
      if((cigPosUI + 3) >= *cigLenUI )
      { /*If: want more cigar memory*/
         tmpStr =
            realloc(
               *cigTypeStr,
               (cigPosUI + expandCig_fun06)
                  * sizeof(signed char)
            );
         if(! tmpStr)
            goto memErr_fun06_sec07;
         *cigTypeStr = tmpStr;

         tmpStr =
            (signed char *)
            realloc(
               *cigArySI,
               (cigPosUI + expandCig_fun06)
                  * sizeof(signed int)
            );
         if(! tmpStr)
            goto memErr_fun06_sec07;
         *cigArySI = (signed int *) tmpStr;

         *cigLenUI += (expandCig_fun06 - 8);
      } /*If: want more cigar memory*/

      /**************************************************\
      * Fun06 Sec04 Sub02:
      *   - insertion cases + cigar switch statment
      \**************************************************/

      switch(*dirMatrixSC)
      { /*Switch: find direction (snp/del/ins)*/
         case def_mvStop_alnDefs:
            break;

         case def_mvIns_alnDefs:
         /*Case: is a insertion*/
            if( (*cigTypeStr)[cigPosUI] != 'I' )
            { /*If: need to make new cigar entry*/

               if( (*cigTypeStr)[cigPosUI] )
                  ++cigPosUI; /*if not null*/

               (*cigTypeStr)[cigPosUI] = 'I';
               (*cigArySI)[cigPosUI] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++(*cigArySI)[cigPosUI];

            --qryPosUI;
            dirMatrixSC -= (lenRefUI + 1);
            break;
         /*Case: is a insertion*/

         /***********************************************\
         * Fun06 Sec04 Sub03:
         *   - snp or match cases
         \***********************************************/

         case def_mvSnp_alnDefs:
         /*Case: is a snp or match*/
            matchBl =
               getMatch_alnSet(
                  qrySeqStr[qryPosUI],
                  refSeqStr[refPosUI],
                  alnSetSTPtr
               ); /*find if had a match*/

            if(matchBl & def_ntEql_alnDefs)
            { /*If: had a match*/
               /*count number of anonymous matches*/
               *numAnonUI +=
                  !!(matchBl & def_anonymous_alnDefs);

               matchBl = '='; /*match*/
            } /*If: had a match*/

            else
               matchBl = 'X'; /*mismatch*/

            if( (*cigTypeStr)[cigPosUI] != matchBl )
            { /*If: need to make new cigar entry*/
               if( (*cigTypeStr)[cigPosUI] )
                  ++cigPosUI; /*if not null*/

               (*cigTypeStr)[cigPosUI] = matchBl;
               (*cigArySI)[cigPosUI] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++(*cigArySI)[cigPosUI];

            --refPosUI;
            --qryPosUI;

            dirMatrixSC -= (lenRefUI + 2);
            break;
         /*Case: is a snp or match*/

         /***********************************************\
         * Fun06 Sec04 Sub04:
         *   - deletion cases
         \***********************************************/

         case def_mvDel_alnDefs:
         /*Case: is a deletion*/
            if( (*cigTypeStr)[cigPosUI] != 'D' )
            { /*If: need to make new cigar entry*/
               if( (*cigTypeStr)[cigPosUI])
                  ++cigPosUI; /*if not null*/

               (*cigTypeStr)[cigPosUI] = 'D';
               (*cigArySI)[cigPosUI] = 1;
            } /*If: need to make new cigar entry*/

            else
               ++(*cigArySI)[cigPosUI];

            --refPosUI;
            --dirMatrixSC;
            break;
         /*Case: is a deletion*/
      } /*Switch: find direction (snp/del/ins)*/
   } /*Loop: trace alignment path*/

   /*****************************************************\
   * Fun06 Sec04 Sub05:
   *   - move to first base
   \*****************************************************/

   /*account for ending on one index behind*/
   ++refPosUI;
   ++qryPosUI;

   if(! revBl)
      *refStartUI = refPosUI; /*index 0*/

   if(
         revBl
      &&
           (refPosUI + matrixSTPtr->refOffsetSL)
         < matrixSTPtr->refEndSL
      && delAtEndBl
   ){ /*If: deletions at end*/
      
      if((*cigTypeStr)[cigPosUI] > 32)
         ++cigPosUI;

      (*cigTypeStr)[cigPosUI] = 'D';
      (*cigArySI)[cigPosUI] = refPosUI;
   } /*If: deletions at end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec05:
   ^   - add starting softmasked bases and invert cigar
   ^   o fun06 sec05 sub01:
   ^     - add starting softmasked bases to cigar
   ^   o fun06 sec05 sub02:
   ^     - if forward sequence invert cigar; is backwards
   ^   o fun06 sec05 sub04:
   ^     - set return to index 0 and add null at end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec05 Sub01:
   *   - add starting softmasked bases to cigar
   \*****************************************************/

   if(qryPosUI > 0)
   { /*If: need to add starting softmasked bases*/
      ++cigPosUI;
      (*cigTypeStr)[cigPosUI] = 'S';

      (*cigArySI)[cigPosUI] = qryPosUI;
      *maskCntUI += qryPosUI;
   } /*If: need to add starting softmasked bases*/

   /*****************************************************\
   * Fun06 Sec05 Sub02:
   *   - if forward sequence invert cigar; is backwards
   \*****************************************************/

   if(! revBl)
   { /*If: need to reverse (invert) cigar*/
      qryPosUI = startUI;
      refPosUI = cigPosUI;

      while(qryPosUI < refPosUI)
      { /*Loop: invert the cigar to forward direction*/
         /*swap was from standfords bithacking guide*/

         (*cigTypeStr)[qryPosUI]^=(*cigTypeStr)[refPosUI];
         (*cigTypeStr)[refPosUI]^=(*cigTypeStr)[qryPosUI];
         (*cigTypeStr)[qryPosUI]^=(*cigTypeStr)[refPosUI];

         (*cigArySI)[qryPosUI] ^= (*cigArySI)[refPosUI];
         (*cigArySI)[refPosUI] ^= (*cigArySI)[qryPosUI];
         (*cigArySI)[qryPosUI] ^= (*cigArySI)[refPosUI];

         ++qryPosUI;
         --refPosUI;
      } /*Loop: invert the cigar to forward direction*/

      /*null/0 or softmask at start flipped to end*/
   } /*If: need to reverse (invert) cigar*/

   if(startUI > 0)
   { /*If: have pervious entries*/
      qryPosUI = startUI - 1;
      refPosUI = startUI;

      while(refPosUI <= cigPosUI)
      { /*Loop: move entries back*/
         if(
            (*cigTypeStr)[qryPosUI]
               == (*cigTypeStr)[refPosUI]
         ) (*cigArySI)[qryPosUI] += (*cigArySI)[refPosUI];
           /*same entry*/

         else if(! delAtEndBl)
         { /*Else If: masking wanted*/
            if(
                  (*cigTypeStr)[qryPosUI] == 'S'
               && (*cigTypeStr)[refPosUI] == 'D'
            ) (*cigArySI)[qryPosUI] +=
                  (*cigArySI)[refPosUI];
   
            else if(
                  (*cigTypeStr)[qryPosUI] == 'D'
               && (*cigTypeStr)[refPosUI] == 'S'
            ){ /*Else If: deletion next to softmasked*/
               (*cigTypeStr)[qryPosUI] = 'S';
               (*cigArySI)[qryPosUI] +=
                  (*cigArySI)[refPosUI];
            } /*Else If: deletion next to softmasked*/
   
            else if(
                  (*cigTypeStr)[qryPosUI] == 'S'
               && (*cigTypeStr)[refPosUI] == 'I'
            ) (*cigArySI)[qryPosUI] +=
                 (*cigArySI)[refPosUI];
   
            else if(
                  (*cigTypeStr)[qryPosUI] == 'I'
               && (*cigTypeStr)[refPosUI] == 'S'
            ){ /*If: softmask next to insertion*/
               (*cigTypeStr)[qryPosUI] = 'S';
               (*cigArySI)[qryPosUI] +=
                  (*cigArySI)[refPosUI];
            } /*If: softmask next to insertion*/
   
            else if(
                  (*cigTypeStr)[qryPosUI] == 'S'
               && (*cigTypeStr)[refPosUI] == 'X'
            ) (*cigArySI)[qryPosUI] +=
                 (*cigArySI)[refPosUI];
   
            else if(
                  (*cigTypeStr)[qryPosUI] == 'X'
               && (*cigTypeStr)[refPosUI] == 'S'
            ){ /*If: softmask next to snp*/
               (*cigTypeStr)[qryPosUI] = 'S';
               (*cigArySI)[qryPosUI] +=
                  (*cigArySI)[refPosUI];
            } /*If: softmask next to snp*/

            else
            { /*Else: no issues*/
               ++qryPosUI;

               (*cigTypeStr)[qryPosUI] =
                  (*cigTypeStr)[refPosUI];

               (*cigArySI)[qryPosUI] =
                  (*cigArySI)[refPosUI];
            } /*Else: no issues*/
         } /*Else If: masking wanted*/

         else
         { /*Else: not softmasking ends*/
            ++qryPosUI;

            (*cigTypeStr)[qryPosUI] =
               (*cigTypeStr)[refPosUI];

            (*cigArySI)[qryPosUI] =
               (*cigArySI)[refPosUI];
         } /*Else: not softmasking ends*/

         ++refPosUI;
      } /*Loop: move entries back*/

      (*cigTypeStr)[qryPosUI + 1] = '\0';
      cigPosUI = qryPosUI;
   } /*If: have pervious entries*/

   /*****************************************************\
   * Fun06 Sec05 Sub03:
   *   - set return to index 0 and add null at end
   \*****************************************************/

   if( (*cigTypeStr)[cigPosUI] != '\0' )
      (*cigTypeStr)[cigPosUI + 1] = '\0';

   else
      --cigPosUI; /*is index 1*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec07:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return cigPosUI;

   memErr_fun06_sec07:;
      return -1;
} /*getCig_dirMatrix*/

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
