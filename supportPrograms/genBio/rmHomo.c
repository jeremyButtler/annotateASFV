/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' rmHomo SOF: Start Of File
'   - removes and or adjusts indels in homopolymers
'   o header:
'     - included libraries
'   o fun01: indel_rmHomo
'     - remove indels from homopolymers
'   o fun02: leftAlnIndel_rmHomo
'     - forces all homopolymer indels to be left or right
'       aligned
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

#include "rmHomo.h"

#include "../genLib/ulCp.h"
#include "samEntry.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - std #include <stdio.h>
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "../genLib/endLine.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: indel_rmHomo
|   - remove indels from homopolymers
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to sequence to remove
|       indels from homopolymers
|   - refStr:
|     o c-string with reference sequence in samSTPtr is
|       mapped to
|   - minHomoSI:
|     o minimum homopolymer length to remove indel from
|   - maxIndelSI:
|     o maximum indel size to remove
|   - maskSC:
|     o base to mask deletions with
|     o if 0 if(! maskSC), then uses reference base
|   - scanBl:
|     o scan both neighboring bases next to indel to
|       see if the indel is next to or in middle of large
|       homopolymer
|   - seqBuffStrPtr:
|     o c-string pointer to buffer to use in copy sequence
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeSeqUIPtr:
|     o unsigned long pointer to with size of
|       seqBuffStrPtr
|   - qBuffStrPtr:
|     o c-string pointer to buffer to use in copy q-score
|       entry
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeQUIPtr:
|     o unsigned long pointer to with size of qBuffStrPtr
| Output:
|   - Modifies:
|     o seqStr in samSTPtr to have no homopolymer indels
|     o seqBuffStrPtr and qBuffStrPtr to be resized (or
|       allocated) to sequence length + deletions
|     o sizeQUIPtr and sizeSeqUIPtr to have their buffer
|       size if seqBuffStrPtr or qBuffStrPtr are resized
|   - Returns:
|     o 0 for no errors
|     o 1 for memory errors
\-------------------------------------------------------*/
signed char
indel_rmHomo(
   struct samEntry *samSTPtr,/*sequence with indels*/
   signed char *refStr,      /*reference sequence*/
   signed int minHomoSI,     /*min homoplymer length*/
   signed int maxIndelSI,    /*indel to large to remove*/
   signed char maskSC,       /*base to mask dels with*/
   signed char scanBl,       /*scan neighbor bases*/
   signed char **seqBuffStrPtr, /*sequence buffer*/
   unsigned int *sizeSeqUIPtr,  /*size of seqBuffStrPtr*/
   signed char **qBuffStrPtr,   /*q-score buffer*/
   unsigned int *sizeQUIPtr     /*size of qBuffStrPtr*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - remove indels from homopolymers
   '   o fun01 sec01:
   '     - variable declarations 
   '   o fun01 sec02:
   '     - memory allocation
   '   o fun01 sec03:
   '     - indel removal
   '   o fun01 sec04:
   '     - add null to end of c-strings (sequence/cigar)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char qUC = 0;

   unsigned int uiCig = 0;
   unsigned int cpCigUI = 0;

   unsigned int seqPosUI = 0;
   unsigned int dupPosUI = 0;

   unsigned int refPosUI = 0;
   unsigned int tmpPosUI = 0;
   unsigned int lenHomoUI = 0;

   /*for finding homopolymer around base*/
   unsigned int lenForHomoUI = 0;
   unsigned int lenBackHomoUI = 0;
   signed char forSC = 0;
   signed char backSC = 0;

   unsigned int lenBuffUI = 0;
   signed char qBl = 0; /*1: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - memory allocation
   ^   o fun01 sec02 sub01:
   ^     - get maximum possible length
   ^   o fun01 sec02 sub02:
   ^     - resize sequence buffers
   ^   o fun01 sec02 sub03:
   ^     - resize q-score buffer (if have)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - get maximum possible length
   \*****************************************************/

   lenBuffUI = samSTPtr->readLenUI;
   lenBuffUI += samSTPtr->delCntUI;
      /*deletions mean adding in a base, so need to have
      `  some extra room
      */

   lenBuffUI += 8;
      /*prevents read overflows from ulCp functions*/

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - resize sequence buffers
   \*****************************************************/

   if(*sizeSeqUIPtr < lenBuffUI)
   { /*If: need to resize buffer*/
      if(*seqBuffStrPtr)
         free(*seqBuffStrPtr);
      *seqBuffStrPtr = 0;
   } /*If: need to resize buffer*/

   if(! *seqBuffStrPtr)
   { /*If: need new sequence buffer*/
      *sizeSeqUIPtr = 0;

      *seqBuffStrPtr =
         malloc(lenBuffUI * sizeof(signed char));
      if(! *seqBuffStrPtr)
         goto memErr_fun01_sec05;

      *sizeSeqUIPtr = lenBuffUI;
   } /*If: need new sequence buffer*/

   /*****************************************************\
   * Fun01 Sec02 Sub03:
   *   - resize q-score buffer (if have)
   \*****************************************************/

   if(! samSTPtr->qStr)
      qBl = 1;
   else if(samSTPtr->qStr[0] == '\0')
      qBl = 1;
   else if(
         samSTPtr->qStr[0] == '*'
      && samSTPtr->qStr[1] == '\0'
   ) ; /*no q-score entry*/

   else
      qBl = 1;

   if(qBl)
   { /*If: have q-score entry*/
      if(*sizeQUIPtr < lenBuffUI)
      { /*If: need to resize buffer*/
         if(*qBuffStrPtr)
            free(*qBuffStrPtr);
         *qBuffStrPtr = 0;
      } /*If: need to resize buffer*/

      if(! *qBuffStrPtr)
      { /*If: need new quence buffer*/
         *sizeQUIPtr = 0;

         *qBuffStrPtr =
            malloc(lenBuffUI * sizeof(signed char));
         if(! *qBuffStrPtr)
            goto memErr_fun01_sec05;

         *sizeQUIPtr = lenBuffUI;
      } /*If: need new quence buffer*/
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - indel removal
   ^   o fun01 sec03 sub01:
   ^     - copy matches/snps/soft masked bases
   ^   o fun01 sec03 sub02:
   ^     - deal with hard masking cases
   ^   o fun01 sec03 sub03:
   ^     - find homopolymer length
   ^   o fun01 sec03 sub04:
   ^     - check if discarding indel
   ^   o fun01 sec03 sub05:
   ^     - remove indels
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec03 Sub01:
   *  - copy matches/snps/soft masked bases
   \*****************************************************/

   refPosUI = samSTPtr->refStartUI;

   for(
      cpCigUI = 0;
      cpCigUI < samSTPtr->cigLenUI;
      ++cpCigUI
   ){ /*Loop: remove indels*/

      samSTPtr->cigTypeStr[uiCig] =
         samSTPtr->cigTypeStr[cpCigUI];
      samSTPtr->cigArySI[uiCig] =
         samSTPtr->cigArySI[cpCigUI];

      if(
            samSTPtr->cigTypeStr[uiCig] == 'M'
         || samSTPtr->cigTypeStr[uiCig] == '='
         || samSTPtr->cigTypeStr[uiCig] == 'X'
         || samSTPtr->cigTypeStr[uiCig] == 'S'
      ){ /*If: non-indel mutation*/
         cpLen_ulCp(
            &(*seqBuffStrPtr)[dupPosUI],
            &samSTPtr->seqStr[seqPosUI],
            samSTPtr->cigArySI[uiCig]
         );

         if(qBl)
         { /*If: have q-score entry*/
            cpLen_ulCp(
               &(*qBuffStrPtr)[dupPosUI],
               &samSTPtr->qStr[seqPosUI],
               samSTPtr->cigArySI[uiCig]
            );
         } /*If: have q-score entry*/

         dupPosUI += samSTPtr->cigArySI[uiCig];
         seqPosUI += samSTPtr->cigArySI[uiCig];

         if(samSTPtr->cigTypeStr[uiCig] != 'S')
            refPosUI += samSTPtr->cigArySI[uiCig];

         ++uiCig;

         continue;
      }  /*If: non-indel mutation*/

      /**************************************************\
      * Fun01 Sec03 Sub02:
      *  - deal with hard masking cases
      \**************************************************/

      if(
            samSTPtr->cigTypeStr[uiCig] != 'D'
         && samSTPtr->cigTypeStr[uiCig] != 'I'
      ){ /*If: not an indel (hard masking)*/
         ++uiCig;
         continue;
      }  /*If: not an indel (hard masking)*/

      /**************************************************\
      * Fun01 Sec03 Sub03:
      *   - find homopolymer length
      *   o fun01 sec03 sub03 cat01:
      *     - find strict homopolyer length
      *   o fun01 sec03 sub03 cat02:
      *     - find length of neighbor homopolymers
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun01 Sec03 Sub03 Cat01:
      +   - find strict homopolyer length
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(samSTPtr->cigArySI[uiCig] > maxIndelSI)
         goto keepIndel_fun01_sec03_sub04;
         /*indel is to large to remove*/


      lenHomoUI = 1;

      if(refStr[refPosUI] != '\0')
      { /*If: have next bases (potentail homopolymer)*/
         tmpPosUI = refPosUI + 1;

         while(refStr[tmpPosUI] == refStr[refPosUI])
         { /*Loop: find end of homopolymer*/
            ++lenHomoUI;
            ++tmpPosUI;
         } /*Loop: find end of homopolymer*/
      } /*If: have next bases (potentail homopolymer)*/



      if(refPosUI > 0)
      { /*If: have previous bases (homopolymer?)*/
         tmpPosUI = refPosUI - 1;

         while(refStr[tmpPosUI] == refStr[refPosUI])
         { /*Loop: find end of homopolymer*/
            ++lenHomoUI;

            if(! tmpPosUI)
               break;

            --tmpPosUI;
         } /*Loop: find end of homopolymer*/
      } /*If: have previous bases (homopolymer?)*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun01 Sec03 Sub03 Cat02:
      +   - find length of neighbor homopolymers
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(! scanBl)
         ;
      else if(lenHomoUI < (unsigned int) minHomoSI)
      { /*Else If: homopolymer will be to small*/
         lenForHomoUI = 0;
         lenBackHomoUI = 0;
         forSC = 0; 
         backSC = -1; /*never want these to match*/

         if(refStr[refPosUI] != '\0')
         { /*If: have next bases (potentail homopolymer)*/
            tmpPosUI = refPosUI + 1;
            forSC = refStr[tmpPosUI];
            ++tmpPosUI;

            while(refStr[tmpPosUI] == forSC)
            { /*Loop: find end of homopolymer*/
               ++lenForHomoUI;
               ++tmpPosUI;
            } /*Loop: find end of homopolymer*/
         } /*If: have next bases (potentail homopolymer)*/

         if(refPosUI > 0)
         { /*If: have previous bases (homopolymer?)*/
            tmpPosUI = refPosUI - 1;
            backSC = refStr[tmpPosUI];
            --tmpPosUI;

            while(refStr[tmpPosUI] == backSC)
            { /*Loop: find end of homopolymer*/
               ++lenBackHomoUI;

               if(! tmpPosUI)
                  break;

               --tmpPosUI;
            } /*Loop: find end of homopolymer*/
         } /*If: have previous bases (homopolymer?)*/

         if(forSC == backSC)
            lenForHomoUI += lenBackHomoUI;
         else if(lenForHomoUI < lenBackHomoUI)
            lenForHomoUI = lenBackHomoUI;

         if(lenHomoUI < lenForHomoUI)
            lenHomoUI = lenForHomoUI;
      } /*Else If: homopolymer will be to small*/

      /**************************************************\
      * Fun01 Sec03 Sub04:
      *  - check if discarding indel
      \**************************************************/

      if(lenHomoUI < (unsigned int) minHomoSI)
      { /*If: homoplymer to small*/
         keepIndel_fun01_sec03_sub04:;

         if(samSTPtr->cigTypeStr[uiCig] == 'I')
         { /*If: need to copy insertion*/

            tmpPosUI = seqPosUI;
            tmpPosUI += samSTPtr->cigArySI[uiCig];

            while(seqPosUI < tmpPosUI)
            { /*Loop: copy insertion*/
               if(qBl)
                  (*qBuffStrPtr)[dupPosUI] =
                     samSTPtr->qStr[seqPosUI];

               (*seqBuffStrPtr)[dupPosUI++] =
                  samSTPtr->seqStr[seqPosUI++];
            } /*Loop: copy insertion*/
         } /*If: need to copy insertion*/

         else
            refPosUI += samSTPtr->cigArySI[uiCig];
            /*need to update reference position*/

         ++uiCig; /*keep cigar entry*/
         continue;
      } /*If: homoplymer to small*/
 
      /**************************************************\
      * Fun01 Sec03 Sub05:
      *  - remove indels
      \**************************************************/

      if(samSTPtr->cigTypeStr[uiCig] == 'D')
      { /*If: removing deletions*/
         tmpPosUI = refPosUI;
         tmpPosUI += samSTPtr->cigArySI[uiCig];

         while(refPosUI < tmpPosUI)
         { /*Loop: add reference bases in (remove del)*/
            if(qBl)
            { /*If: have q-score entries to add to*/
               /*I am using 0, becuase replacing a
               `  deletion with a reference base means no
               `  confidence in the call
               */
               ++samSTPtr->qHistUI[0];
               (*qBuffStrPtr)[dupPosUI] =
                  def_adjQ_samEntry; /*adjQ is q-score 0*/
            } /*If: have q-score entries to add to*/

            if(maskSC)
            { /*If: user wanted base masked*/
               (*seqBuffStrPtr)[dupPosUI++] = maskSC;
               ++refPosUI;
            } /*If: user wanted base masked*/

            else
               (*seqBuffStrPtr)[dupPosUI++] =
                  refStr[refPosUI++];
         } /*Loop: add reference bases in (remove del)*/

         samSTPtr->delCntUI -= samSTPtr->cigArySI[uiCig];
         samSTPtr->matchCntUI +=
            samSTPtr->cigArySI[uiCig];
         samSTPtr->readLenUI += samSTPtr->cigArySI[uiCig];


         tmpPosUI += samSTPtr->cigTypeStr[uiCig] = '=';

         if(samSTPtr->cigTypeStr[uiCig - 1] == '=')
         { /*If: this  cigar entry can be merged*/
            samSTPtr->cigArySI[uiCig - 1] +=
               samSTPtr->cigArySI[uiCig];
         }  /*If: this cigar entry can be merged*/

         else
            ++uiCig; /*can not merge*/
      } /*If: removing deletions*/


      else
      { /*Else: insertion*/
         samSTPtr->insCntUI -= samSTPtr->cigArySI[uiCig];
         samSTPtr->readLenUI -= samSTPtr->cigArySI[uiCig];

         tmpPosUI = seqPosUI;
         tmpPosUI += samSTPtr->cigArySI[uiCig];


         if(qBl)
         { /*If: have quality score entry*/
            while(seqPosUI < tmpPosUI)
            { /*Loop: remove insertion quality scores*/
               qUC =
                  (unsigned char)
                  samSTPtr->qStr[seqPosUI++];
               qUC -= def_adjQ_samEntry;

               --samSTPtr->qHistUI[qUC];
               samSTPtr->sumQUL -= qUC;
            } /*Loop: remove insertion quality scores*/
         } /*If: have quality score entry*/

         else
         { /*Else: no quality score entry*/
            seqPosUI = tmpPosUI;
            continue;
         } /*Else: no quality score entry*/
      } /*Else: insertion*/
   } /*Loop: remove indels*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - add null to end of c-strings (sequence/cigar)
   ^   o fun01 sec04 sub01:
   ^     - merge duplicate cigar entries
   ^   o fun01 sec04 sub02:
   ^     - copy sequence and q-score entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *   - merge duplicate cigar entries
   \*****************************************************/

   samSTPtr->cigLenUI = uiCig;
   uiCig = 0;

   for(
     cpCigUI = 1;
     cpCigUI < samSTPtr->cigLenUI;
     ++cpCigUI
   ){ /*Loop: merge duplicate cigar entries*/
      if(
            samSTPtr->cigTypeStr[uiCig]
         == samSTPtr->cigTypeStr[cpCigUI]
      ){ /*If: can merge cigar entries*/
         samSTPtr->cigArySI[uiCig] +=
            samSTPtr->cigArySI[cpCigUI];
      }  /*If: can merge cigar entries*/

      else
      { /*Else: moving to next cigar entry*/
         ++uiCig;
         samSTPtr->cigTypeStr[uiCig] =
            samSTPtr->cigTypeStr[cpCigUI];
         samSTPtr->cigArySI[uiCig] =
            samSTPtr->cigArySI[cpCigUI];
      } /*Else: moving to next cigar entry*/
   }  /*Loop: merge duplicate cigar entries*/

   ++uiCig;
   samSTPtr->cigLenUI = uiCig;
   samSTPtr->cigArySI[uiCig] = '\0';
   samSTPtr->cigTypeStr[uiCig] = '\0';

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *   - copy sequence and q-score entries
   \*****************************************************/

   (*seqBuffStrPtr)[dupPosUI] = '\0';
   (*qBuffStrPtr)[dupPosUI] = '\0';

   if(samSTPtr->seqSizeUI < lenBuffUI)
   { /*If: new sequence buffer is larger*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc(lenBuffUI * sizeof(signed char));
      if(! samSTPtr->seqStr)
         goto memErr_fun01_sec05;

      samSTPtr->seqSizeUI = lenBuffUI - 8;
   } /*If: new sequence buffer is larger*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      (*seqBuffStrPtr),
      samSTPtr->readLenUI
   );

   if(qBl)
   { /*If: have q-score entry*/
      if(samSTPtr->qSizeUI < lenBuffUI)
      { /*If: current q-score buffer is larger*/
         free(samSTPtr->qStr);

         samSTPtr->qStr =
            malloc(lenBuffUI * sizeof(signed char));
         if(! samSTPtr->qStr)
            goto memErr_fun01_sec05;

         samSTPtr->qSizeUI = lenBuffUI - 8;
      } /*If: current q-score buffer is larger*/

      cpLen_ulCp(
         samSTPtr->qStr,
         (*qBuffStrPtr),
         samSTPtr->readLenUI
      );

      findQScores_samEntry(samSTPtr);
         /*need to re-adjust q-scores for removed indels*/
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^   - return and clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   qUC = 0;
   goto ret_fun01_sec05;

   memErr_fun01_sec05:;
      qUC = 1;
      goto ret_fun01_sec05;

   ret_fun01_sec05:;
      return (signed char) qUC;
} /*indel_rmHomo*/

/*-------------------------------------------------------\
| Fun02: leftAlnIndel_rmHomo
|   - forces all homopolymer indels to be left or right
|     aligned
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to sequence to remove
|       indels from homopolymers
|   - refStr:
|     o c-string with reference sequence in samSTPtr is
|       mapped to
| Output:
|   - Modifies:
|     o seqStr and qStr in samSTPtr to have indel
|       homopolymers left or right aligned
|   - Returns:
|     o 0 for no errors
|     o 1 for memory errors
\-------------------------------------------------------*/
signed char
leftAlnIndel_rmHomo(
   struct samEntry *samSTPtr,/*sequence with indels*/
   signed char *refStr       /*reference sequence*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - forces all homopolymer indels to be left or right
   '     aligned
   '   o fun02 sec01:
   '     - variable declarations 
   '   o fun02 sec02:
   '     - memory allocation and q-score entry detection
   '   o fun02 sec03:
   '     - indel removal
   '   o fun02 sec04:
   '     - add null to end of c-strings (sequence/cigar)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *cigHeapStr = 0;
   signed int *cigHeapArySI = 0;

   unsigned char qUC = 0;

   unsigned int uiCig = 0;
   signed int tmpCigSI = 0;
   signed int tmpNtSI = 0;
   unsigned int cpCigUI = 0;
   unsigned int newCigUI = 0; /*position in new cigar*/
   signed int extraNtSI = 0;

   unsigned int seqPosUI = 0;
   unsigned int dupPosUI = 0;

   unsigned int refPosUI = 0;
   unsigned int tmpPosUI = 0;

   /*homopolyer lengths*/
   signed char ntSC = 0;          /*nucleotide of indel*/
   unsigned int tmpUI = 0;
   unsigned int backHomoLenUI = 0;
   unsigned int forHomoLenUI = 0;

   unsigned int lenBuffUI = 0;
   signed char qBl = 0; /*1: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - memory allocation and q-score entry detection
   ^   o fun02 sec02 sub01:
   ^     - get maximum possible length
   ^   o fun02 sec02 sub02:
   ^     - detect if have q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - get maximum possible length
   \*****************************************************/

   if(! samSTPtr->seqStr[0])
      goto noSeq_fun02_sec05;
   else if(samSTPtr->seqStr[0] == '*')
      goto noSeq_fun02_sec05;

   lenBuffUI = samSTPtr->readLenUI;
   lenBuffUI += samSTPtr->delCntUI;
      /*deletions mean adding in a base, so need to have
      `  some extra room
      */

   lenBuffUI += 8;
      /*prevents read overflows from ulCp functions*/

   cigHeapStr = malloc(lenBuffUI * sizeof(signed char));
   if(! cigHeapStr)
      goto memErr_fun02_sec05;
   cigHeapStr[0] = 0;

   cigHeapArySI = malloc(lenBuffUI * sizeof(signed int));
   if(! cigHeapArySI)
      goto memErr_fun02_sec05;
   cigHeapArySI[0] = 0;

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - detect if have q-score entry
   \*****************************************************/

   if(! samSTPtr->qStr)
      qBl = 1;
   else if(samSTPtr->qStr[0] == '\0')
      qBl = 1;
   else if(
         samSTPtr->qStr[0] == '*'
      && samSTPtr->qStr[1] == '\0'
   ) ; /*no q-score entry*/

   else
      qBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - indel removal
   ^   o fun02 sec03 sub01:
   ^     - copy matches/snps/soft masked bases
   ^   o fun02 sec03 sub02:
   ^     - deal with hard masking cases
   ^   o fun02 sec03 sub03:
   ^     - find homopolymer length
   ^   o fun02 sec03 sub04:
   ^     - shift indels
   ^   o fun02 sec03 sub04:
   ^     - check if can shift indels
   ^   o fun02 sec03 sub05:
   ^     - shift deletions
   ^   o fun02 sec03 sub06:
   ^     - shift insertions
   ^   o fun02 sec03 sub07:
   ^     - move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec03 Sub01:
   *  - copy matches/snps/soft masked bases
   \*****************************************************/

   newCigUI = 0;
   refPosUI = samSTPtr->refStartUI;

   for(
      cpCigUI = 0;
      cpCigUI < samSTPtr->cigLenUI;
      ++cpCigUI
   ){ /*Loop: remove indels*/

      if(
            samSTPtr->cigTypeStr[uiCig] == 'M'
         || samSTPtr->cigTypeStr[uiCig] == '='
         || samSTPtr->cigTypeStr[uiCig] == 'X'
         || samSTPtr->cigTypeStr[uiCig] == 'S'
      ){ /*If: non-indel mutation*/
         cigHeapStr[newCigUI] =
            samSTPtr->cigTypeStr[uiCig];
         cigHeapArySI[newCigUI] =
            samSTPtr->cigArySI[uiCig];

         dupPosUI += samSTPtr->cigArySI[uiCig];
         seqPosUI += samSTPtr->cigArySI[uiCig];

         if(samSTPtr->cigTypeStr[uiCig] != 'S')
            refPosUI += samSTPtr->cigArySI[uiCig];

         ++newCigUI;
         goto nextCigarEntry_fun02_sec03_sub07;
      }  /*If: non-indel mutation*/

      /**************************************************\
      * Fun02 Sec03 Sub02:
      *  - deal with hard masking cases
      \**************************************************/

      if(
            samSTPtr->cigTypeStr[uiCig] != 'D'
         && samSTPtr->cigTypeStr[uiCig] != 'I'
      ){ /*If: not an indel (hard masking)*/
         cigHeapStr[newCigUI] =
            samSTPtr->cigTypeStr[uiCig];
         cigHeapArySI[newCigUI] =
            samSTPtr->cigArySI[uiCig];

         refPosUI += samSTPtr->cigArySI[uiCig];

         ++newCigUI;
         goto nextCigarEntry_fun02_sec03_sub07;
      }  /*If: not an indel (hard masking)*/

      /**************************************************\
      * Fun02 Sec03 Sub03:
      *   - detect if in homopolymer
      *   o fun02 sec03 sub03 cat01:
      *     - check if insertion matches bases around it
      *   o fun02 sec03 sub03 cat02:
      *     - find forward and reverse homopolmer lengths
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec03 Sub03 Cat01:
      +   - check if insertion matches bases around it
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(samSTPtr->cigTypeStr[uiCig] == 'D')
         ntSC = refStr[refPosUI]; /*if deletion*/

      else if(
         samSTPtr->seqStr[seqPosUI] == refStr[refPosUI -1]
      ) ntSC = samSTPtr->seqStr[seqPosUI];

      else if(
         samSTPtr->seqStr[seqPosUI] == refStr[refPosUI]
      ) ntSC = samSTPtr->seqStr[seqPosUI];

      else
      { /*Else: inserton is different*/
         cigHeapStr[newCigUI] =
            samSTPtr->cigTypeStr[uiCig];
         cigHeapArySI[newCigUI] =
            samSTPtr->cigArySI[uiCig];

         seqPosUI += samSTPtr->cigArySI[uiCig];

         ++newCigUI;
         goto nextCigarEntry_fun02_sec03_sub07;
      } /*Else: insertion is different*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec03 Sub03 Cat02:
      +   - find homopolmer length
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      backHomoLenUI = 0;
      forHomoLenUI = 0;

      if(refPosUI > 0)
      { /*If: may have previous bases (homopolymer?)*/
         tmpPosUI = refPosUI - 1;

         while(refStr[tmpPosUI] == ntSC)
         { /*Loop: find end of homopolymer*/
            ++backHomoLenUI;

            if(! tmpPosUI)
               break;

            --tmpPosUI;
         } /*Loop: find end of homopolymer*/

      } /*If: may have previous bases (homopolymer?)*/

      if(refStr[refPosUI])
      { /*If: not end of sequence*/
         tmpPosUI = refPosUI + 1;

         while(refStr[tmpPosUI++] == ntSC)
            ++forHomoLenUI;
      } /*If: not end of sequence*/

      /**************************************************\
      * Fun02 Sec03 Sub04:
      *   - check if can shift indels
      \**************************************************/

      if(! backHomoLenUI)
      { /*Else If: no homo or indel already left aligned*/
         cigHeapStr[newCigUI] =
            samSTPtr->cigTypeStr[uiCig];
         cigHeapArySI[newCigUI] =
            samSTPtr->cigArySI[uiCig];

         if(samSTPtr->cigTypeStr[uiCig] == 'D')
            refPosUI += samSTPtr->cigArySI[uiCig];
         else
            seqPosUI += samSTPtr->cigArySI[uiCig];

         ++newCigUI;
         goto nextCigarEntry_fun02_sec03_sub07;
      } /*Else If: no homo or indel already left aligned*/

      /**************************************************\
      * Fun02 Sec03 Sub05:
      *   - shift deletions
      \**************************************************/

      else if(samSTPtr->cigTypeStr[uiCig] == 'D')
      { /*Else: need to align a deletion*/
         tmpCigSI = cigHeapArySI[newCigUI - 1];
         tmpPosUI = refPosUI;
         tmpUI = 0;

         if(
               cigHeapStr[newCigUI - 1] == '='
            || cigHeapStr[newCigUI - 1] == 'M'
         ){ /*If: may have a match*/
            while(tmpCigSI > 0)
            { /*Loop: find insert position*/
               if(refStr[tmpPosUI] != ntSC)
                  break;
               ++tmpUI;
               --tmpPosUI;
               --tmpCigSI;

               if(tmpPosUI <= 0)
                  break;
               if(tmpUI == backHomoLenUI)
                  break; /*at start of indel*/
            }  /*Loop: find insert position*/
         }  /*If: may have a match*/

         cigHeapArySI[newCigUI - 1] -= tmpUI;
         cigHeapStr[newCigUI] = 'D';

         if(
              forHomoLenUI
            < (unsigned int) samSTPtr->cigArySI[uiCig]
         ) cigHeapArySI[newCigUI] = forHomoLenUI + 1;
               /*+1 to account for not counting deletion*/

         else
            cigHeapArySI[newCigUI] =
               samSTPtr->cigArySI[uiCig];

         /*check if did not shift all deletions*/
         if(
              cigHeapArySI[newCigUI]
            < samSTPtr->cigArySI[uiCig]
         ) extraNtSI = 
                samSTPtr->cigArySI[uiCig]
              - cigHeapArySI[newCigUI];
         else
            extraNtSI = 0;

         ++newCigUI;

         /*finish the new cigar entry for the old del*/
         if(tmpUI > 0)
         { /*If: need to shift bases around*/
            if(cigHeapStr[newCigUI - 2] == '=')
               cigHeapStr[newCigUI] = '=';
            else
               cigHeapStr[newCigUI] = 'M';

            cigHeapArySI[newCigUI] = tmpUI;
            ++newCigUI;
         } /*If: need to shift bases around*/

         if(extraNtSI)
         { /*If: had extra deletions*/
            cigHeapStr[newCigUI] = 'D';
            cigHeapArySI[newCigUI] = extraNtSI;
            ++newCigUI;
         } /*If: had extra deletions*/

         refPosUI += samSTPtr->cigArySI[uiCig];
         goto nextCigarEntry_fun02_sec03_sub07;
      } /*Else: need to align a deletion*/

      /**************************************************\
      * Fun02 Sec03 Sub06:
      *   - shift insertions
      *   o fun02 sec03 sub06 cat01:
      *     - find new insertion position and shift length
      *   o fun02 sec03 sub06 cat02:
      *     - swap quality scores around
      *   o fun02 sec03 sub06 cat03:
      *     - insert insertion into new cigar position
      *     - shift insertion cigar entry
      *     * ins matches homopolymer so safe to do after
      *       swaping q-scores
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec03 Sub06 Cat01:
      +   - find new insertion position and shift length
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: need to shift insertion*/
         tmpCigSI = cigHeapArySI[newCigUI - 1];
         tmpPosUI = refPosUI - 1;
            /*ins is always one ahead*/
         tmpUI = 0;

         if(
               cigHeapStr[newCigUI - 1] == '='
            || cigHeapStr[newCigUI - 1] == 'M'
         ){ /*If: may have a match*/
            while(tmpCigSI && refStr[tmpPosUI] == ntSC)
            { /*Loop: find insert position*/
               ++tmpUI;
               --tmpPosUI;
               --tmpCigSI;

               if(tmpPosUI <= 0)
                  break;
               if(tmpUI == backHomoLenUI)
                  break; /*at start of indel*/
            }  /*Loop: find insert position*/
         }  /*If: may have a match*/


         /*find insertion length to shift*/
         extraNtSI = 1;

         if(samSTPtr->cigArySI[uiCig] > 1)
         { /*If: need to find number of insertions*/
            tmpCigSI = samSTPtr->cigArySI[uiCig] - 1;
               /*-1 to account for first insertion*/
            tmpNtSI = seqPosUI + 1;
               /*+1 to get of first insertion*/

            while(
                  tmpCigSI
               && samSTPtr->seqStr[tmpNtSI] == ntSC
            ){ /*Loop: find insertion length*/
               ++extraNtSI;
               ++tmpNtSI;
               --tmpNtSI;
               --tmpCigSI;
            } /*Loop: find insertion length*/
         } /*If: need to find number of insertions*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun02 Sec03 Sub06 Cat02:
         +   - swap quality scores around
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(qBl)
         { /*If: have q-score entry to swap*/
            tmpPosUI = seqPosUI - backHomoLenUI;

            for(
               tmpNtSI = 0;
               tmpNtSI < extraNtSI;
               ++tmpNtSI
            ){ /*Loop: swap quality score entries*/
               samSTPtr->qStr[tmpPosUI + tmpNtSI] ^=
                  samSTPtr->qStr[seqPosUI + tmpNtSI];
               samSTPtr->qStr[seqPosUI + tmpNtSI] ^=
                  samSTPtr->qStr[tmpPosUI + tmpNtSI];
               samSTPtr->qStr[tmpPosUI + tmpNtSI] ^=
                  samSTPtr->qStr[seqPosUI + tmpNtSI];
            } /*Loop: swap quality score entries*/
         } /*If: have q-score entry to swap*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun02 Sec03 Sub06 Cat03:
         +   - shift insertion cigar entry
         +   * ins matches homopolymer so safe to do after
         +     swaping q-scores
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         cigHeapArySI[newCigUI - 1] -= tmpUI;
         cigHeapStr[newCigUI] = 'I';
         cigHeapArySI[newCigUI] = extraNtSI;

         /*check if did not shift all insertions*/
         if(extraNtSI < samSTPtr->cigArySI[uiCig])
            extraNtSI = 
               samSTPtr->cigArySI[uiCig] - extraNtSI;
         else
            extraNtSI = 0;

         ++newCigUI;

         /*finish the new cigar entry for the old ins*/
         if(tmpUI > 0)
         { /*If: need to shift bases around*/
            if(cigHeapStr[newCigUI - 2] == '=')
               cigHeapStr[newCigUI] = '=';
            else
               cigHeapStr[newCigUI] = 'M';

            cigHeapArySI[newCigUI] = tmpUI;
            ++newCigUI;
         } /*If: need to shift bases around*/

         if(extraNtSI)
         { /*If: had extra insertions*/
            cigHeapStr[newCigUI] = 'I';
            cigHeapArySI[newCigUI] = extraNtSI;
            ++newCigUI;
         } /*If: had extra insertions*/

         seqPosUI += samSTPtr->cigArySI[uiCig];
         goto nextCigarEntry_fun02_sec03_sub07;
      } /*Else: need to shift insertion*/
 
      /**************************************************\
      * Fun02 Sec03 Sub07:
      *   - move to next cigar entry
      \**************************************************/

      nextCigarEntry_fun02_sec03_sub07:;
         if(newCigUI < 2)
            ;
         else if(
               cigHeapStr[newCigUI - 1]
            == cigHeapStr[newCigUI - 2]
         ){ /*Else If: can merge cigar entries*/
            cigHeapArySI[newCigUI - 2] +=
               cigHeapArySI[newCigUI - 1];
            --newCigUI;

            cigHeapArySI[newCigUI] = 0;
            cigHeapStr[newCigUI] = 0;
         }  /*Else If: can merge cigar entries*/

         ++uiCig;
   } /*Loop: remove indels*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - add null to end of c-strings (sequence/cigar)
   ^   o fun02 sec04 sub01:
   ^     - merge duplicate cigar entries
   ^   o fun02 sec04 sub02:
   ^     - copy new cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec04 Sub01:
   *   - merge duplicate cigar entries & remove 0 entries
   \*****************************************************/

   uiCig = 0;

   for(
     cpCigUI = 1;
     cpCigUI < newCigUI;
     ++cpCigUI
   ){ /*Loop: merge duplicate cigar entries*/
      if(cigHeapStr[uiCig] == cigHeapStr[cpCigUI])
         cigHeapArySI[uiCig] += cigHeapArySI[cpCigUI];
         /*can merge entries*/

      else if(cigHeapArySI[cpCigUI] <= 0)
         continue; /*cigar entry is 0 (no bases)*/

      else
      { /*Else: moving to next cigar entry*/
         ++uiCig;
         cigHeapStr[uiCig] = cigHeapStr[cpCigUI];
         cigHeapArySI[uiCig] = cigHeapArySI[cpCigUI];
      } /*Else: moving to next cigar entry*/
   }  /*Loop: merge duplicate cigar entries*/

   ++uiCig;
   cigHeapStr[uiCig] = 0;
   cigHeapArySI[uiCig] = 0;
   newCigUI = uiCig;

   /*****************************************************\
   * Fun02 Sec04 Sub02:
   *   - copy new cigar entry
   \*****************************************************/

   if(newCigUI > samSTPtr->cigSizeUI)
   { /*If: new cigar is longer then old cigar memory*/
      free(samSTPtr->cigTypeStr);
      free(samSTPtr->cigArySI);

      samSTPtr->cigTypeStr = cigHeapStr;
      samSTPtr->cigArySI = cigHeapArySI;

      cigHeapStr = 0;
      cigHeapArySI = 0;

      samSTPtr->cigSizeUI = newCigUI;
   } /*If: new cigar is longer then old cigar memory*/

   else
   { /*Else: need to copy in new cigar entries*/
      cpLen_ulCp(
         samSTPtr->cigTypeStr,
         cigHeapStr,
         newCigUI
      );

      cpLen_ulCp(
         (signed char *) samSTPtr->cigArySI,
         (signed char *) cigHeapArySI,
         newCigUI * (sizeof(signed int))
      );
   } /*Else: need to copy in new cigar entries*/

   samSTPtr->cigLenUI = newCigUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec05:
   ^   - return and clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   qUC = 0;
   goto ret_fun02_sec05;

   noSeq_fun02_sec05:;
      qUC = 0;
      goto ret_fun02_sec05;
 
   memErr_fun02_sec05:;
      qUC = 1;
      goto ret_fun02_sec05;

   ret_fun02_sec05:;
      if(cigHeapStr)
         free(cigHeapStr);
      cigHeapStr = 0;

      if(cigHeapArySI)
         free(cigHeapArySI);
      cigHeapArySI = 0;

      return (signed char) qUC;
} /*leftAlnIndel_rmHomo*/

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
