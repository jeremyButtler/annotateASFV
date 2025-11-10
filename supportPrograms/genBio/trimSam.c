/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' trimSam SOF: Start Of File
'   - holds functions for trimming softmasked regoins at
'     the start and end of sequences in an sam file
'   o header:
'     - Included libraries and defined variables
'   o fun01 seq_trimSam:
'     o Trim soft mask regions off end of sam entry
'   o fun02: coords_trimSam
'     - Trim an sam file entry by coordinates
'   o fun03: ends_trimSam
'     - trims x (user input) bases off ends of read
'   o fun04 trimReads_trimSam:
'     o Trims soft mask regions for all reads with a
'       sequence in a sam file
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries and defined variables
|   o header Sec01:
|     - Included libraries
|   o header Sec02:
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec01:
^   - Included libraries
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "trimSam.h"

#include <stdio.h>

#include "samEntry.h"
#include "../genLib/ulCp.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden depenencies
!   o .c  #include "../genLib/base10StrToNum.h"
!   o .c  #include "../genLib/numToStr.h"
!   o .c  #include "../genLib/fileFun.h"
!   o .h  #include "../genLib/endLine.h"
!   o .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec02:
^   - Defined varialbes
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*This .c version is slightly slower, but it works*/
/*This will work so long as the integer is beneath
`   128 bits, which is non-realistic value.
*/
#define def_uiToUC_trimSam \
 ((sizeof(unsigned int) >>1) - (sizeof(unsigned int) >>3))

/*Logic:
`   - sizeof(unsigned int) >> 1:
`     o This gets the number of bytes to shift to get
`       an 32 bit or 16 bit number to one bytes. Basically
`       it is divide by 2
`   - (sizeof(unsigned int) >> 3)
`     o This adds in 1 if I have eight bytes. It fails
`       for ints and shorts
`   - (sizeof(int) >> 1) - (sizeof(int) >> 3)
`     o Gives me the number of bits I need to shift to
`       get the number of characters in an int
`   - (sizeof(int) >> 3) << (sizeof(int) >> 4)
`     o This is not used, but would correct for a 128 bit
`       number. Each time I increase the bit size I need
`       to add a correction in
*/

/*-------------------------------------------------------\
| Fun01: seq_trimSam
| Use:
|  - Trims off the soft masked regions of a sam entry
| Input:
|  - samSTPtr:
|    o samEntry struct with sam entry to trim.
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_header_trimSam if header (invalid and ignored)
|    o def_noMap_trimSam if unmapped read (no reference)
|    o def_noSeq_trimSam if no sequence line
|  - Modifies:
|    o Trims cigar, sequence, & q-score entries in samST.
\-------------------------------------------------------*/
signed char
seq_trimSam(
    struct samEntry *samSTPtr   /*has sam line to trim*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: seq_trimSam
   '  - Trims soft masked regions at start & end of sam
   '    entry
   '  o fun01 sec01:
   '    - Variable declerations
   '  o fun01 sec02:
   '    - Check if I have something to trim
   '  o fun01 sec03:
   '    - Finding trimming postions and trim the cigar
   '  o fun01 sec04:
   '    - trim reads
   '  o fun01 sec05:
   '    - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    signed char errSC = 0;
    signed int startPosSI = 0;
    signed int endTrimSI = 0;
        
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec02:
    ^   - Check if I have something to trim
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if is an header is an header entry*/
    if(samSTPtr->extraStr[0] == '@')
       goto header_fun01_sec05;

    /*Check if is an unmapped read*/
    if(samSTPtr->flagUS & 4)
       goto noMap_fun01_sec05;

    /*Check if is an unmapped read*/
    if(samSTPtr->cigTypeStr[0] == '*')
       goto noMap_fun01_sec05;

    /*Check if has a sequence to trim*/
    if(samSTPtr->seqStr[0] == '*')
       goto noSeq_fun01_sec05;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec03:
    ^   - Finding trimming postions and trim the cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Position my self at the start of the sequence*/
    if(samSTPtr->cigTypeStr[0] == 'S')
    { /*If: need to remove starting softmask bases*/
       startPosSI = samSTPtr->cigArySI[0];
       --samSTPtr->cigLenUI;
       samSTPtr->readLenUI -= startPosSI;

       cpLen_ulCp(
          samSTPtr->cigTypeStr,
          &samSTPtr->cigTypeStr[1],
          samSTPtr->cigLenUI
       ); /*Trim the cigar types array*/

       cpLen_ulCp(
          (signed char *) samSTPtr->cigArySI,
          (signed char *) &samSTPtr->cigArySI[1],
          samSTPtr->cigLenUI * sizeof(unsigned int)
       ); /*Trim the cigar types array*/
    } /*If: need to remove starting softmask bases*/

    if(samSTPtr->cigTypeStr[samSTPtr->cigLenUI -1] == 'S')
    { /*If: need to remove ending softmasked bases*/
       --samSTPtr->cigLenUI;
       endTrimSI = samSTPtr->cigArySI[samSTPtr->cigLenUI];
       samSTPtr->readLenUI -= endTrimSI;
       samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = 0;
       samSTPtr->cigArySI[samSTPtr->cigLenUI] = 0;
    } /*If: need to remove ending softmasked bases*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec04:
    ^   - trim reads
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(startPosSI)
    { /*If: trimming of starting bases*/
       cpLen_ulCp(
          samSTPtr->seqStr,
          &samSTPtr->seqStr[startPosSI],
          samSTPtr->readLenUI
       );

       if(! samSTPtr->qStr[0]) ;
       else if(! samSTPtr->qStr[1]) ;
       else
          cpQEntry_samEntry(
             samSTPtr,
             &samSTPtr->qStr[startPosSI],
             1 /*blank q-score histograms*/
          ); /*also finds mean/median q-score*/
    } /*If: trimming of starting bases*/

    else if(endTrimSI)
    { /*Else If: only trimming ending bases*/
       samSTPtr->seqStr[samSTPtr->readLenUI] = 0;

       if(! samSTPtr->qStr[0]) ;
       else if(! samSTPtr->qStr[1]) ;
       else
       { /*Else: have q-score entry*/
          samSTPtr->qStr[samSTPtr->readLenUI] = 0;
          findQScores_samEntry(samSTPtr);
       } /*Else: have q-score entry*/
    } /*Else If: only trimming ending bases*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec05:
    ^   - return 
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errSC = 0;
    goto cleanUp_fun01_sec05;

    header_fun01_sec05:;
       errSC = def_header_trimSam;
       goto cleanUp_fun01_sec05;

    noMap_fun01_sec05:;
       errSC = def_noMap_trimSam;
       goto cleanUp_fun01_sec05;

    noSeq_fun01_sec05:;
       errSC = def_noSeq_trimSam;
       goto cleanUp_fun01_sec05;

    cleanUp_fun01_sec05:;
       return errSC;
} /*seq_trimSam*/

/*-------------------------------------------------------\
| Fun02: coords_trimSam
|   - Trim an sam file entry by coordinates
| Input:
|   - samSTPtr:
|     o Pointer to an sam entry structure with an read to
|       trim
|   - startSI:
|     o Singed integer with the starting coordinate
|   - endSI:
|     o Singed integer with the ending coordinate
|     o < 1 is treated as do not trim
| Output:
|   - Modifies:
|     o samSTPtr to be trimmed
|   - Returns:
|     o 0 for no errors
|     o def_rangeErr_trimSam for coordinates out of range
|     o def_header_trimSam if header (invalid and ignored)
|     o def_noMap_trimSam if unmapped read (no reference)
|     o def_noSeq_trimSam if no sequence line
\-------------------------------------------------------*/
signed char
coords_trimSam(
   struct samEntry *samSTPtr,
   signed int startSI,
   signed int endSI
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - Trim an sam file entry by coordinates
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - check ranges and set up variables
   '   o fun02 sec03:
   '     - trim read
   '   o fun02 sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;

   signed int siCig = 0;
   signed int cigBaseOnSI = 0;
   signed int refPosSI = 0;
   signed int seqPosSI = 0;

   signed int seqEndPosSI = 0;
   signed int siEndCig = 0;
   signed int cigEndBaseSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - check ranges and set up variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if is an header is an header entry*/
    if(samSTPtr->extraStr[0] == '@')
       goto header_fun02_sec05;

    /*Check if is an unmapped read*/
    if(samSTPtr->flagUS & 4)
       goto noMap_fun02_sec05;

    if(samSTPtr->cigTypeStr[0] == '*')
       goto noMap_fun02_sec05;

    /*Check if has a sequence to trim*/
    if(samSTPtr->seqStr[0] == '*')
       goto noSeq_fun02_sec05;

   if(samSTPtr->refStartUI > (unsigned int) endSI)
      goto outOfRange_fun02_sec05;

   if(samSTPtr->refEndUI < (unsigned int) startSI)
      goto outOfRange_fun02_sec05;

   cigBaseOnSI = samSTPtr->cigArySI[0];
   refPosSI = (signed int) samSTPtr->refStartUI;

   /*Setting to zero becuase it is easier to recount*/
   samSTPtr->maskCntUI = 0;
   samSTPtr->insCntUI = 0;
   samSTPtr->delCntUI = 0;
   samSTPtr->snpCntUI = 0;
   samSTPtr->matchCntUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - trim read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   findRefPos_samEntry(
      samSTPtr,
      &siCig,
      &cigBaseOnSI,
      startSI,
      &refPosSI,
      &seqPosSI
   );

   samSTPtr->cigArySI[siCig] = cigBaseOnSI;
   samSTPtr->refStartUI = (unsigned int) refPosSI;

   seqEndPosSI = seqPosSI;
   siEndCig = siCig;
   cigEndBaseSI = cigBaseOnSI;

   if(samSTPtr->refEndUI < (unsigned int) endSI)
      endSI = samSTPtr->refEndUI;
   else
      samSTPtr->refEndUI = (unsigned int) endSI;

   findRefPos_samEntry(
      samSTPtr,
      &siEndCig,
      &cigEndBaseSI,
      endSI,
      &refPosSI,
      &seqEndPosSI
   );

   samSTPtr->cigArySI[siEndCig] -= cigEndBaseSI;
   samSTPtr->readLenUI = seqEndPosSI - seqPosSI;

   if(! siEndCig) ; /*at start (nothing kept)*/
   else if(! samSTPtr->cigArySI[siEndCig])
      --siEndCig;
      

   samSTPtr->cigLenUI = siEndCig - siCig + 1;

   cpLen_ulCp(
      samSTPtr->cigTypeStr,
      &samSTPtr->cigTypeStr[siCig],
      samSTPtr->cigLenUI
   );

   cpLen_ulCp(
      (signed char *) samSTPtr->cigArySI,
      (signed char *) &samSTPtr->cigArySI[siCig],
      samSTPtr->cigLenUI * sizeof(signed int)
   );

   cpLen_ulCp(
      samSTPtr->seqStr,
      &samSTPtr->seqStr[seqPosSI],
      samSTPtr->readLenUI
   );

   if(! samSTPtr->qStr[0]) ;
   else if(! samSTPtr->qStr[1]) ;
   else
      cpQEntry_samEntry(
         samSTPtr,
         &samSTPtr->qStr[seqPosSI],
         1 /*blank q-score histograms*/
      ); /*also finds mean/median q-score*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - find alignment stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      siCig = 0;
      siCig < (signed int) samSTPtr->cigLenUI;
      ++siCig
   ){ /*Loop: refind counts*/
      switch(samSTPtr->cigTypeStr[siCig])
      { /*Switch: count entries*/
         case 'M':
         case '=':
            samSTPtr->matchCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         case 'X':
            samSTPtr->snpCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         case 'I':
            samSTPtr->insCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         case 'D':
            samSTPtr->delCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         case 'S':
            samSTPtr->maskCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
      } /*Switch: count entries*/
   }  /*Loop: refind counts*/

   samSTPtr->alnReadLenUI = samSTPtr->matchCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->snpCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->delCntUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec05:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun02_sec05;

   outOfRange_fun02_sec05:;
      errSC = def_rangeErr_trimSam;
      goto cleanUp_fun02_sec05;

   header_fun02_sec05:;
      errSC = def_header_trimSam;
      goto cleanUp_fun02_sec05;

   noMap_fun02_sec05:;
      errSC = def_noMap_trimSam;
      goto cleanUp_fun02_sec05;

   noSeq_fun02_sec05:;
      errSC = def_noSeq_trimSam;
      goto cleanUp_fun02_sec05;

   cleanUp_fun02_sec05:;
      return errSC;
} /*coords_trimSam*/

/*-------------------------------------------------------\
| Fun03: ends_trimSam
|   - trims x (user input) bases off ends of read
| Input:
|   - samSTPtr:
|     o Pointer to an sam entry structure with an read to
|       trim
|   - startSI:
|     o number bases to trim off the start of the read 
|   - endSI:
|     o number bases to trim off the end of the read
|   - strictBl:
|     o 3: trim extactly only both ends
|     o 2: trim extactly endSI bases (end only)
|     o 1: trim extactly startSI bases (start only)
|     o 0: trim at until match is found after startSI and
|          endSI
| Output:
|   - Modifies:
|     o seqSTPtr, qStr, cigTypeStr, and cigArySI in
|       samSTPtr to be trimmed
|   - Returns:
|     o 0 for no errors
|     o def_rangeErr_trimSam for coordinates out of range
|     o def_noSeq_trimSam if no sequence line
\-------------------------------------------------------*/
signed char
ends_trimSam(
   struct samEntry *samSTPtr,
   signed int startSI,
   signed int endSI,
   signed char strictBl
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - trims x (user input) bases off ends of read
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - trim of sequence and quality score bases
   '   o fun03 sec03:
   '     - trim sequence, q-score, and cigar entries
   '   o fun03 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int ntCntSI = 0; /*number nucleotides trimming*/
   signed int rmSI = 0;

   signed int cigStartSI = 0;
   signed int startLeftSI = 0;
      /*bases left in start entry at end of cigar*/
   signed int cigEndSI = 0;
   signed int endLeftSI = 0;
      /*bases left in end entry at end of cigar*/
  
   unsigned int matchCntUI = samSTPtr->matchCntUI;
   unsigned int snpCntUI = samSTPtr->snpCntUI;
   unsigned int insCntUI = samSTPtr->insCntUI;
   unsigned int delCntUI = samSTPtr->delCntUI;
   unsigned int maskCntUI = samSTPtr->maskCntUI;
   unsigned int readLenUI = samSTPtr->readLenUI;
   unsigned int refStartUI = samSTPtr->refStartUI;
   unsigned int refEndUI = samSTPtr->refEndUI;
   unsigned int alnReadLenUI = samSTPtr->alnReadLenUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if can trim and find new length
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ntCntSI = startSI + endSI;

   if(samSTPtr->readLenUI <= (unsigned int) ntCntSI)
      goto rangeErr_fun03_sec04;
   else if(*samSTPtr->seqStr == '*')
      goto noSeq_fun03_sec04;
   else if(! *samSTPtr->seqStr)
      goto noSeq_fun03_sec04;
   else if(startSI)
      ;
   else if(! endSI)
      goto done_fun03_sec04; /*nothing to trim*/

   samSTPtr->readLenUI -= startSI;
   samSTPtr->readLenUI -= endSI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - trim cigar entries
   ^   o fun03 sec03 sub01:
   ^     - trim starting cigar entries
   ^   o fun03 sec03 sub02:
   ^     - trim ending cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - trim starting cigar entries
   *   o fun03 sec03 sub01 cat01:
   *     - check if at end of cigar + setup for loop
   *   o fun03 sec03 sub01 cat02:
   *     - get trim count for non-deletion cases
   *   o fun03 sec03 sub01 cat03:
   *     - deal with deletions and hard masks
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun03 Sec03 Sub01 Cat01:
   +   - check if at end of cigar + setup for loop
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(startSI <= 0)
      goto trimEnd_fun03_sec04_sub02;

   ntCntSI = 0;
   cigStartSI = 0;

   while(ntCntSI < startSI || ! (strictBl & 1) )
   { /*Loop: find number entries to trim*/
      if(cigStartSI >= (signed int) samSTPtr->cigLenUI)
         goto rangeErr_fun03_sec04;
         /*trimmed entire cigar/read*/

      if(
            samSTPtr->cigTypeStr[cigStartSI] == '='
         || samSTPtr->cigTypeStr[cigStartSI] == 'M'
      ){ /*If: on match or snp/match*/
         if(ntCntSI >= startSI)
            goto trimEnd_fun03_sec04_sub02;
         else if(
              ntCntSI + samSTPtr->cigArySI[cigStartSI]
            > startSI
         ){  /*Else If: only part of match trimmed*/
            ntCntSI += samSTPtr->cigArySI[cigStartSI];
            rmSI = ntCntSI;
            rmSI -= startSI;
            rmSI= samSTPtr->cigArySI[cigStartSI] -rmSI;

            startLeftSI = samSTPtr->cigArySI[cigStartSI];
            startLeftSI -= rmSI;

            refStartUI += rmSI;
            alnReadLenUI -= rmSI;
            readLenUI -= rmSI;
            matchCntUI -= rmSI;

            goto trimEnd_fun03_sec04_sub02;
         }  /*Else If: only part of match trimmed*/
      }  /*If: on match or snp/match*/
 
      /*~~~++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub01 Cat02:
      +   - get trim count for non-deletion cases
      \+++~~~+++++++++++++++++++++++++++++++++++++++++++*/

      switch(samSTPtr->cigTypeStr[cigStartSI])
      { /*Switch: find number bases in cigar entry*/
         case '=':
         case 'M':
         case 'X':
         case 'I':
         case 'S':
         /*Case: match, snp/match, snp, insertion, mask*/
            ntCntSI += samSTPtr->cigArySI[cigStartSI];

            if(! (strictBl & 1) )
               rmSI = samSTPtr->cigArySI[cigStartSI];

            else if(ntCntSI >= startSI)
            { /*Else If: overshot on match or mask entry*/
               rmSI = ntCntSI;
               rmSI -= startSI;
               rmSI= samSTPtr->cigArySI[cigStartSI] -rmSI;

               startLeftSI=samSTPtr->cigArySI[cigStartSI];
               startLeftSI -= rmSI;
            } /*Else If: overshot on match or mask entry*/

            else
               rmSI = samSTPtr->cigArySI[cigStartSI];

            if(samSTPtr->cigTypeStr[cigStartSI] == '=')
            { /*If: match*/
               refStartUI += rmSI;
               alnReadLenUI -= rmSI;
               matchCntUI -= rmSI;
            } /*If: match*/

            else if(
               samSTPtr->cigTypeStr[cigStartSI] == 'X'
            ){ /*Else If: snp*/
               refStartUI += rmSI;
               alnReadLenUI -= rmSI;
               snpCntUI -= rmSI;
            } /*Else If: snp*/

            else if(
               samSTPtr->cigTypeStr[cigStartSI] == 'M'
            ){ /*Else If: match or snp*/
               refStartUI += rmSI;
               alnReadLenUI -= rmSI;
               matchCntUI -= rmSI;
            } /*Else If: match or snp*/

            else if(
               samSTPtr->cigTypeStr[cigStartSI] == 'I'
            ) insCntUI -= rmSI;

            else
               maskCntUI -= rmSI;

            readLenUI -= rmSI;

            if(ntCntSI < startSI)
               ++cigStartSI;
            else if(! (strictBl & 1) )
               ++cigStartSI;
            break;
         /*Case: match, snp/match, snp, insertion, mask*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun03 Sec03 Sub01 Cat03:
         +   - deal with deletions and hard masks
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'D':
         /*Case: deletion*/
            rmSI = samSTPtr->cigArySI[cigStartSI];
            refStartUI += rmSI;
            alnReadLenUI -= rmSI;
            delCntUI -= rmSI;

            ++cigStartSI;
            break; /*deletion or hard masked*/
         /*Case: deletion*/

         default:
         /*Case: some weird hard mask*/
            ++cigStartSI;
            break; /*deletion or hard masked*/
         /*Case: some weird hard mask*/
      } /*Switch: find number bases in cigar entry*/
   } /*Loop: find number entries to trim*/

   /*****************************************************\
   * Fun03 Sec03 Sub02:
   *   - trim ending cigar entries
   *   o fun03 sec03 sub02 cat01:
   *     - check if at start of cigar + setup for loop
   *   o fun03 sec03 sub02 cat02:
   *     - check if match case is ending case
   *   o fun03 sec03 sub02 cat03:
   *     - handle non-deleltion cases (including match)
   *   o fun03 sec03 sub02 cat04:
   *     - handle deletion and hardmask cases
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun03 Sec03 Sub02 Cat01:
   +   - check if at start of cigar + setup for loop
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   trimEnd_fun03_sec04_sub02:;
   cigEndSI = samSTPtr->cigLenUI - 1;
   ntCntSI = 0;

   if(endSI <= 0)
      goto trim_fun03_sec04;

   while(ntCntSI < endSI || ! (strictBl & 2) )
   { /*Loop: trim end of sam entry*/

      if(cigEndSI == cigStartSI)
      { /*If: on last entry*/
         if(! startLeftSI)
            goto match_fun03_sec03_sub02;

         else if(ntCntSI + startLeftSI >= endSI)
         { /*Else If: found end*/
            rmSI = ntCntSI + startLeftSI;
            rmSI -= endSI;
            rmSI = startLeftSI - rmSI;

            endLeftSI = startLeftSI;
            endLeftSI -= rmSI;

            alnReadLenUI -= rmSI;
            readLenUI -= rmSI;
            matchCntUI -= rmSI;

            ntCntSI = endSI;
            goto trim_fun03_sec04;
         } /*Else If: found end*/

         else
            goto rangeErr_fun03_sec04;
            /*trimmed entire cigar/read*/
      } /*If: on last entry*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat02:
      +   - check if match case is ending case
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(
            samSTPtr->cigTypeStr[cigEndSI] == '='
         || samSTPtr->cigTypeStr[cigEndSI] == 'M'
      ){ /*If: could be a match*/
         match_fun03_sec03_sub02:;

         if(ntCntSI >= endSI)
            goto trim_fun03_sec04;
         else if(
              ntCntSI + samSTPtr->cigArySI[cigEndSI]
            > endSI
         ){  /*Else If: only part of match trimmed*/
            ntCntSI += samSTPtr->cigArySI[cigEndSI];
            rmSI = ntCntSI;
            rmSI -= endSI;
            rmSI = samSTPtr->cigArySI[cigEndSI] - rmSI;

            endLeftSI = samSTPtr->cigArySI[cigEndSI];
            endLeftSI -= rmSI;

            alnReadLenUI -= rmSI;
            readLenUI -= rmSI;
            matchCntUI -= rmSI;

            goto trim_fun03_sec04;
         }  /*Else If: only part of match trimmed*/
      }  /*If: could be a match*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat03:
      +   - handle non-deleltion cases; includes match
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      switch(samSTPtr->cigTypeStr[cigEndSI])
      { /*Switch: find number bases in cigar entry*/
         case '=':
         case 'M':
         case 'X':
         case 'I':
         case 'S':
         /*Case: match or snp/match entry*/
            ntCntSI += samSTPtr->cigArySI[cigEndSI];

            if(! (strictBl & 2) )
               rmSI = samSTPtr->cigArySI[cigEndSI];

            else if(ntCntSI > endSI)
            { /*Else If: overshot on match or mask entry*/
               rmSI = ntCntSI;
               rmSI -= endSI;
               rmSI = samSTPtr->cigArySI[cigEndSI] - rmSI;

               endLeftSI = samSTPtr->cigArySI[cigEndSI];
               endLeftSI -= rmSI;
            } /*Else If: overshot on match or mask entry*/

            else
               rmSI = samSTPtr->cigArySI[cigEndSI];

            if(samSTPtr->cigTypeStr[cigEndSI] == '=')
            { /*If: match*/
               matchCntUI -= rmSI;
               refEndUI -= rmSI;
               alnReadLenUI -= rmSI;
            } /*If: match*/

            else if(samSTPtr->cigTypeStr[cigEndSI] == 'X')
            { /*Else If: snp*/
               snpCntUI -= rmSI;
               refEndUI -= rmSI;
               alnReadLenUI -= rmSI;
            } /*Else If: snp*/

            else if(samSTPtr->cigTypeStr[cigEndSI] == 'M')
            { /*Else If: match or snp*/
               matchCntUI -= rmSI;
               refEndUI -= rmSI;
               alnReadLenUI -= rmSI;
            } /*Else If: match or snp*/

            else if(samSTPtr->cigTypeStr[cigEndSI] == 'I')
               insCntUI -= rmSI;
            else
               maskCntUI -= rmSI;

            readLenUI -= rmSI;

            if(ntCntSI < endSI)
               --cigEndSI;
            else if(! (strictBl & 2) )
               --cigEndSI;
            break;
         /*Case: match or snp/match entry*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun03 Sec03 Sub02 Cat04:
         +   - handle deletion and hardmask cases
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'D':
         /*Case: deletion*/
            rmSI = samSTPtr->cigArySI[cigEndSI];
            refEndUI -= rmSI;
            alnReadLenUI -= rmSI;
            delCntUI -= rmSI;

            --cigEndSI;
            break; /*deletion or hard masked*/
         /*Case: deletion*/

         default:
            --cigEndSI;
            break; /*deletion or hard masked*/
      } /*Switch: find number bases in cigar entry*/
   } /*Loop: trim end of sam entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - trim sequence, q-score, and cigar entries
   ^   o fun03 sec04 sub01:
   ^     - trim cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec04 Sub01:
   *   - trim cigar entries
   \*****************************************************/

   trim_fun03_sec04:;

   if(! readLenUI)
      goto rangeErr_fun03_sec04;

   if(startLeftSI)
      samSTPtr->cigArySI[cigStartSI] = startLeftSI;

   if(endLeftSI)
      samSTPtr->cigArySI[cigEndSI] = endLeftSI;


   if(cigStartSI)
   { /*If: need to trim cigar entries from start*/
      samSTPtr->cigLenUI = cigEndSI - cigStartSI;
      ++samSTPtr->cigLenUI; /*convert index 0 to index 1*/
      cpLen_ulCp(
         samSTPtr->cigTypeStr,
         &samSTPtr->cigTypeStr[cigStartSI],
         samSTPtr->cigLenUI
      );

      cpLen_ulCp(
         (signed char *) samSTPtr->cigArySI,
         (signed char *) &samSTPtr->cigArySI[cigStartSI],
         samSTPtr->cigLenUI * sizeof(signed int)
      );
   } /*If: need to trim cigar entries from start*/

   else if(cigEndSI < (signed int) samSTPtr->cigLenUI)
   { /*Else If: only trimming end of cigar*/
      samSTPtr->cigLenUI = cigEndSI;
      ++samSTPtr->cigLenUI;
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = 0;
      samSTPtr->cigArySI[samSTPtr->cigLenUI] = 0;
   } /*Else If: only trimming end of cigar*/

   /*****************************************************\
   * Fun03 Sec04 Sub02:
   *   - trim sequence and q-score entries
   \*****************************************************/

   samSTPtr->refStartUI = refStartUI;
   samSTPtr->refEndUI = refEndUI;
   samSTPtr->readLenUI = readLenUI;
   samSTPtr->alnReadLenUI = alnReadLenUI;

   samSTPtr->matchCntUI = matchCntUI;
   samSTPtr->snpCntUI = snpCntUI;
   samSTPtr->insCntUI = insCntUI;
   samSTPtr->delCntUI = delCntUI;
   samSTPtr->maskCntUI = maskCntUI;

   if(! startSI)
   { /*If: not trimming the start*/
      samSTPtr->seqStr[samSTPtr->readLenUI] = 0;

      if(! samSTPtr->qStr)
         ;
      else if(! samSTPtr->qStr[1])
         ;
      else
      { /*Else: have q-score entry*/
         samSTPtr->qStr[samSTPtr->readLenUI] = 0;
         findQScores_samEntry(samSTPtr);
      } /*Else: have q-score entry*/
   } /*If: not trimming the start*/

   else
   { /*Else: trimming bases from the start*/
      cpLen_ulCp(
         samSTPtr->seqStr,
         &samSTPtr->seqStr[startSI],
         samSTPtr->readLenUI
      );

      if(samSTPtr->qStr[0])
         ;
      else if(! samSTPtr->qStr[1])
         ;
      else
         cpQEntry_samEntry(
            samSTPtr,
            &samSTPtr->qStr[startSI],
            1
         );
   } /*Else: trimming bases from the start*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   done_fun03_sec04:;
      return 0;

   rangeErr_fun03_sec04:;
      return def_rangeErr_trimSam;

   noSeq_fun03_sec04:;
      return def_noSeq_trimSam;
} /*ends_trimSam*/

/*-------------------------------------------------------\
| Fun04: trimReads_trimSam
| Use:
|  - Goes though sam file and calls seq_trimSam for each
|    entry
| Input:
|  - samFILE:
|    o Pointer to FILE. Sam file with entries to trim
|  - outFILE:
|    o File to write the trimmed sam file entries to
|  - keepUmapBl:
|    o Also print out entries with unmapped reads
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_memErr_trimSam for memory errors
|  - Prints:
|    o trimmed sam entries with sequences to outFILE, but
|      ignores sam entries without sequences
\-------------------------------------------------------*/
signed char
trimReads_trimSam(
    void *samFILE,           /*sam file to trim*/
    void *outFILE,           /*file to store output*/
    signed char keepUmapBl   /*1: keep unmapped reads*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ' Fun04 TOC: trimReads_trimSam
   '  - Goes though sam file and calls seq_trimSam for
   '    each entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    signed char errSC = 0;  /*Tells if memory error*/
    struct samEntry samStackST;
    
    init_samEntry(&samStackST);
    errSC = setup_samEntry(&samStackST);
    
    if(errSC)
       goto memErr_fun04;
   
    errSC = get_samEntry(&samStackST, (FILE *) samFILE);
    
    while(! errSC)
    { /*While there are lines in sam file to convert*/
        if(*(samStackST.extraStr) == '@')
        { /*If was a header*/
           p_samEntry(&samStackST, 0, (FILE *) outFILE);

           if(errSC)
              goto memErr_fun04;
            
           errSC =
              get_samEntry(&samStackST, (FILE *) samFILE);
           
           if(errSC)
              goto memErr_fun04;

           continue; /*header line, move to next line*/
        } /*If was a header*/
        
        /*Convert & print out sam file entry*/
        errSC = seq_trimSam(&samStackST);

        if(errSC)
        { /*If: had problematic input*/
           if(errSC & def_noMap_trimSam && keepUmapBl)
             p_samEntry(&samStackST, 0, (FILE *) outFILE);
             /*printing umapped reads*/

           /*remainder is no sequence or header*/
        } /*If: had problematic input*/

        else
           p_samEntry(&samStackST, 0, (FILE *) outFILE);
           /*print the trimed entry*/
        
        errSC =
           get_samEntry(&samStackST, (FILE *) samFILE);
    } /*While there are lines in sam file to convert*/
    
    if(errSC == def_memErr_samEntry)
       goto memErr_fun04;
     
    errSC = 0;
    goto cleanUp_fun04;

    memErr_fun04:;
    errSC = def_memErr_trimSam;
    goto cleanUp_fun04;

    cleanUp_fun04:;

    freeStack_samEntry(&samStackST);

    return errSC;
} /*trimReads_trimSam*/

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
