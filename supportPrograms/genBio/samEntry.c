/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' samEntry SOF: Start Of File
'   - struct to hold sam file entry and support functions
'   o header:
'     - Included libraries
'   o defined variables:
'     - system line break variables
'   o .h st01 samEntry:
'     - holds a single samfile entry
'   o .h st02: refs_samEntry
'     - holds reference length & name from sam file header
'   o fun01 blank_samEntry:
'     - sets all non-alloacted variables in samEntry to 0
'   o fun02 init_samEntry:
'     - Initalize a samEntry struct to 0's
'   o fun03: setup_samEntry
'     - allocates memory for a samEntry structure (call
'       after init_samEntry (fun02))
'   o fun04 freeStack_samEntry:
'     - Frees heap allocations in a stack allocated
'       samEntry struct
'   o fun05 freeHeap_samEntry:
'     - Frees a samEntry structer (and sets to null)
'   o fun06: mk_samEntry
'     - Makes an heap allocated samEntry structure
'   o fun07: qhistToMed_samEntry
'     - Gets the median q-score for an histogram of
'       q-scores in a samStruct
'   o fun08: findQScores_samEntry
'      - Gets the median and mean q-scores from a samEntry
'        Structure.
'   o fun09: cpQEntry_samEntry
'     - Copies q-scores from a string into a samEntry
'       structure
'   o fun11: lineTo_samEntry
'     - adds a sam file line to a samEntry struct
'   o fun12: get_samEntry
'     - reads in a single line from a sam file
'   o fun13: findRefPos_samEntry
'     - Find an reference coordinate in an sequence in
'       an sam entry structure
'   o fun14: swap_samEntry
'     - swaps two samEntry structs
'   o fun15: p_samEntry
'     - Prints the sam file entry to a file. This does not
'       print any extra stats that were found.
'   o fun16: pfq_samEntry
'     - Prints the sam entry as a fastq entry to a fastq
'       file
'   o fun17: pfa_samEntry
'     - Prints the sam entry as a fasta entry to a fasta
'       file
'   o fun18: pstats_samEntry
'     - Prints out the stats in a samEntry struct to a file
'   o fun19: revCmp_samEntry
'     - reverse complements a sam file sequence entry
'   o fun20: getHead_samEntry
'     - get header for a sam file
'   o fun21: checkCigSize_samEntry
'     - checks to see if I need to resize cigar and if
'       needed resizes the cigar
'   o .h note01:
'      - Notes about the sam file format from the sam file
'        pdf
'   o license:
'     - licensing for this code (public domain / mit)
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

#include <stdio.h>

#include "samEntry.h"

#include "../genLib/base10str.h"
#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"
#include "../genLib/fileFun.h"

/*These have no .c files*/
#include "../genLib/endLine.h" /*for OS line breaks*/
#include "ntTo5Bit.h"
   /*look up table to see if have anonymous bases*/

/*-------------------------------------------------------\
| Fun01: blank_samEntry
|   - Sets all values to 0, or for c-strings to '\0'
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure to blank
| Output:
|   - Modifies:
|     o Sets every variable to 0
|     o the first character in every c-string to be '\0'
\-------------------------------------------------------*/
void
blank_samEntry(
   struct samEntry *samSTPtr
){
    unsigned int uiIter = 0;

    samSTPtr->qryIdStr[0] = '\0'; /*query id/name*/
    samSTPtr->qryIdLenUC = 0;     /*Length of query id*/

    samSTPtr->refIdStr[0] = '\0'; /*reference id/name*/
    samSTPtr->refIdLenUC = 0;     /*Length; reference id*/

    if(samSTPtr->cigTypeStr)
       samSTPtr->cigTypeStr[0] = '\0';
    if(samSTPtr->cigArySI)
       samSTPtr->cigArySI[0] = 0;
    samSTPtr->cigLenUI = 0;       /*Length; cigar entry*/

    samSTPtr->rNextStr[0] = '*';  /*rNext entry*/
    samSTPtr->rNextStr[1] = '\0'; /*rNext entry*/
    samSTPtr->rnextLenUC = 1;     /*Length of rNext*/

    if(samSTPtr->seqStr)
       samSTPtr->seqStr[0] = '\0';
    if(samSTPtr->qStr)
       samSTPtr->qStr[0] = '\0';

    if(samSTPtr->extraStr)
       samSTPtr->extraStr[0] = '\0';
    samSTPtr->extraLenUI = 0;

    /*Flags/single numeric values in the sam entry*/
    samSTPtr->mapqUC = 0;         /*mapping quality*/
    samSTPtr->flagUS = 0;         /*flag in samEntryCStr*/
    samSTPtr->pNextSI = 0;         /*PNEXT postion*/
    samSTPtr->tLenSI = 0;          /*TLEN*/

    /*Stats*/

    /*The q-scores excluded the softmasked regions*/
    samSTPtr->medianQF = 0;       /*median read q-score*/
    samSTPtr->meanQF = 0;         /*mean read q-score*/

    samSTPtr->refStartUI = 0;     /*First reference base*/
    samSTPtr->refEndUI = 0;       /*Last reference base*/
    samSTPtr->readLenUI = 0;      /*Holds read length*/
    samSTPtr->alnReadLenUI = 0;   /*Aligned length*/

    samSTPtr->matchCntUI = 0;     /*number of matches*/
    samSTPtr->snpCntUI = 0;       /*number of mismatches*/
    samSTPtr->insCntUI = 0;       /*number of insertions*/
    samSTPtr->delCntUI = 0;       /*number of deletions*/
    samSTPtr->maskCntUI = 0;      /*number soft masked*/

    /*These variables are used in finding the q-scores*/
    for(
       uiIter = 0;
       uiIter < def_maxQ_samEntry;
       ++uiIter
    ) samSTPtr->qHistUI[uiIter] = 0;

    samSTPtr->sumQUL = 0;
} /*blank_samEntry*/

/*-------------------------------------------------------\
| Fun02: init_samEntry
|   - Initializes a samEntry structure for use. This 
|     function should only ever be called once per
|     structure or after freeStack_samEntry has been used.
|     Use blank_samEntry to reset the values in a samEntry
|     struct after initialization.
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to initialize
| Output:
|   - Modifies:
|     o all variables in samEntry to be 0 or null
|   - Returns:
|     o 0 for success
|     o 1 for memory error
\-------------------------------------------------------*/
void
init_samEntry(
   struct samEntry *samSTPtr
){
    (samSTPtr)->seqStr = 0;
    (samSTPtr)->seqSizeUI = 0;

    (samSTPtr)->qStr = 0;
    (samSTPtr)->qSizeUI = 0;

    (samSTPtr)->cigTypeStr = 0;
    (samSTPtr)->cigSizeUI = 0;
    (samSTPtr)->cigArySI = 0;

    (samSTPtr)->extraStr = 0;
    (samSTPtr)->extraSizeUI = 0;

    blank_samEntry(samSTPtr);
} /*init_samEntry*/

/*-------------------------------------------------------\
| Fun03: setup_samEntry
|   - allocates memory for a samEntry structure (call
|     after init_samEntry (fun02))
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to allocate memory to
| Output:
|   - Allocates:
|     o memory for seqStr, qStr, cigTypStr, cigArySI,
|       and extraStr
|   - Returns:
|     o 0 for success
|     o def_memErr_samEntry for memory error
\-------------------------------------------------------*/
signed char
setup_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03: setup_samEntry
   '   o fun03 sec01:
   '     - allocate sequence memory
   '   o fun03 sec02:
   '     - allocate q-score memory
   '   o fun03 sec03:
   '     - allocate cigar type (snp/match/ins/del) memory
   '   o fun03 sec04:
   '     - allocate cigar count array
   '   o fun03 sec05:
   '     - allocate extra array
   '   o fun03 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec01:
    ^   - allocate sequence memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->seqStr)
       free(samSTPtr->seqStr);
  
    samSTPtr->seqStr = 0;

    samSTPtr->seqStr =
       calloc(
          1024,
          sizeof(signed char)
       ); /*allocate some memory for sequence*/

    if(! samSTPtr->seqStr)
       goto memErr_fun03_sec06;

    samSTPtr->seqSizeUI = 1024 - 8;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec02:
    ^   - allocate q-score memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->qStr)
       free(samSTPtr->qStr);
  
    samSTPtr->qStr = 0;

    samSTPtr->qStr =
       calloc(
          1024,
          sizeof(signed char)
       ); /*allocate some memory for q-score*/

    if(! samSTPtr->qStr)
       goto memErr_fun03_sec06;

    samSTPtr->qSizeUI = 1024 - 8;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec03:
    ^   - allocate cigar type (snp/match/ins/del) memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigTypeStr)
       free(samSTPtr->cigTypeStr);
  
    samSTPtr->cigTypeStr = 0;

    samSTPtr->cigTypeStr =
       calloc(
          256,
          sizeof(signed char)
       ); /*allocate memory for cigary type*/

    if(! samSTPtr->cigTypeStr)
       goto memErr_fun03_sec06;

    samSTPtr->cigSizeUI = 256 - 8;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec04:
    ^   - allocate cigar count array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigArySI)
       free(samSTPtr->cigArySI);
  
    samSTPtr->cigArySI = 0;

    samSTPtr->cigArySI =
       calloc(
          256,
          sizeof(signed int)
       ); /*allocate memory for cigar counts*/

    if(! samSTPtr->cigArySI)
       goto memErr_fun03_sec06;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec05:
    ^   - allocate extra array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->extraStr)
       free(samSTPtr->extraStr);
  
    samSTPtr->extraStr = 0;

    samSTPtr->extraStr =
       calloc(
          1024,
          sizeof(signed char)
       ); /*allocate memory for extra entry*/

    if(! samSTPtr->extraStr)
       goto memErr_fun03_sec06;

    samSTPtr->extraSizeUI = 1024 - 8;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec06:
    ^   - return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    blank_samEntry(samSTPtr);
    return 0;

    memErr_fun03_sec06:;
       return def_memErr_samEntry;
} /*setup_samEntry*/

/*-------------------------------------------------------\
| Fun04: freeStack_samEntry
| Use:
|   - Frees all variables in samEntry, but not samEntry
| Input:
|   - samSTPtr
|     o Pointer to samEntry struct to free the interal
|       memory of
| Output:
|   - Frees:
|     o allocated memory in samSTPtr
\-------------------------------------------------------*/
void
freeStack_samEntry(
   struct samEntry *samSTPtr
){
    if(! samSTPtr)
       return; /*nothing to do*/

    blank_samEntry((samSTPtr));
    
    free((samSTPtr)->seqStr);
    (samSTPtr)->seqStr = 0;
    (samSTPtr)->seqSizeUI = 0;
    
    free((samSTPtr)->qStr);
    (samSTPtr)->qStr = 0;
    (samSTPtr)->qSizeUI = 0;
    
    free((samSTPtr)->cigTypeStr);
    (samSTPtr)->cigTypeStr = 0;
    (samSTPtr)->cigSizeUI = 0;
    
    free((samSTPtr)->cigArySI);
    (samSTPtr)->cigArySI = 0;
    
    free((samSTPtr)->extraStr);
    (samSTPtr)->extraStr = 0;
    (samSTPtr)->extraSizeUI = 0;
} /*freeStack_samEntry*/

/*-------------------------------------------------------\
| Fun05: freeHeap_samEntry
|   - Frees a samEntry struct
| Input:
|   - samSTPtr
|     o Pointer to Sam entry struct to free
| Output:
|   - Frees:
|     o samSTPtr and its variables with allocated memory
\-------------------------------------------------------*/
void freeHeap_samEntry(
   struct samEntry *samSTPtr
){
    if(! samSTPtr)
       return; /*nothing to do*/

    freeStack_samEntry(samSTPtr);
    free(samSTPtr);
} /*freeHeap_samEntry*/

/*-------------------------------------------------------\
| Fun06: mk_samEntry
|   - Makes an heap allocated samEntry structure
| Input:
| Output:
|   - Returns:
|     o Pointer to a samEntry structure
|     o 0 if had an memory error
\-------------------------------------------------------*/
struct samEntry *
mk_samEntry(
  void
){
  struct samEntry *samST=malloc(sizeof(struct samEntry));
  unsigned char errUC = 0;
  
  if(samST)
  { /*If: I did not have a memory error*/
     init_samEntry(samST);
     errUC = setup_samEntry(samST);

     if(errUC)
     { /*If: memory error*/
        freeHeap_samEntry(samST);
        samST = 0;
     } /*If: memory error*/
  } /*If: I did not have a memory error*/
  
  return samST;
} /*makeSamEntry*/

/*-------------------------------------------------------\
| Fun07: qhistToMed_samEntry
|   - Gets the median q-score for an histogram of q-scores
|     in a samStruct
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct with a q-score
|       histogram to find the median q-score for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
\-------------------------------------------------------*/
void
qhistToMed_samEntry(
   struct samEntry *samSTPtr 
){
    unsigned int numBasesUI = 0;
    unsigned int midPointUL = samSTPtr->readLenUI >> 1;
    unsigned int uiQ = 0;

    for(uiQ = 0; uiQ < def_maxQ_samEntry; uiQ++)
    { /*Loop: through q-score histogram for midpoint*/
        numBasesUI += samSTPtr->qHistUI[uiQ];

        if(numBasesUI >= midPointUL)
        { /*if found the midpoint, then find the median*/
            if(numBasesUI > midPointUL || numBasesUI & 1)
                samSTPtr->medianQF = (float) uiQ;

            else
            { /*Else: even & 2 differnt Q's at midpoint*/
                numBasesUI = uiQ;
                ++uiQ;

                while(samSTPtr->qHistUI[uiQ++] == 0) ;

                --uiQ; /*account for being one off*/

                samSTPtr->medianQF =
                    (numBasesUI + uiQ) / ((float) 2);
            } /*Else: even & 2 differnt Q's at midpoint*/

            break;
        } /*if found the midpoint, then find the median*/
    } /*Loop: through q-score histogram for midpoint*/
} /*qhistToMed_samEntry*/

/*-------------------------------------------------------\
| Fun08: findQScores_samEntry
|   - Gets the median and mean q-scores from a samEntry
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct to find the median and
|       mean q-scores for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
|     o samSTPtr->meanQF to have the mean q-score
|     o samSTPtr->qHistUI to have a histogram of q-scores
|     o samSTPtr->sumQUL to have the sum of all q-scores
\-------------------------------------------------------*/
void
findQScores_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' Fun08 Sec01:
    '   - variable declarations
    '   o fun08 sec02:
    '     - blank q-score variables
    '   o fun08 sec03:
    '     - unsigned long; find most of q-scores
    '   o fun08 sec04:
    '     - finsh of q-scores (byte method, find qscore)
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun08 Sec01:
    ^   - variable declarations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    unsigned int uiQScore = 0;
    unsigned int uiChar = 0;
    unsigned char *scoreAryUC = 0;
    
    #ifdef NOUL
    #else
       unsigned long qScoresUL = 0;
       ulong_ulCp qAdjustUL =
          mkDelim_ulCp((signed char) def_adjQ_samEntry);
       unsigned long *qPtrUL =
          (unsigned long *) samSTPtr->qStr;
    #endif

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun08 Sec02:
    ^   - blank q-score variables
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    uiQScore = 0;
    while(uiQScore < def_maxQ_samEntry)
       samSTPtr->qHistUI[uiQScore++] = 0;
    samSTPtr->sumQUL = 0;
    samSTPtr->medianQF = 0;
    samSTPtr->meanQF = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun08 Sec03:
    ^   - unsigned long; find most of q-scores
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    #ifdef NOUL
    #else
       /*Find the number of q-score characters in buffer*/
       for(
          uiQScore = 0;
          uiQScore <
             samSTPtr->readLenUI >> def_shiftULBy_ulCp;
          ++uiQScore
       ) { /*Loop: Update the q-score historgram and sum*/
          qScoresUL = qPtrUL[uiQScore] - qAdjustUL;
          scoreAryUC = (unsigned char *) &qScoresUL;

          for(
             uiChar = 0;
             uiChar < def_charInUL_ulCp;
             ++uiChar
          ){ /*Loop: Get the q-score entries*/

            ++samSTPtr->qHistUI[scoreAryUC[uiChar]];
            samSTPtr->sumQUL +=
                (unsigned char)
                scoreAryUC[uiChar];
          } /*Loop: Get the q-score entries*/
       } /*Loop: Update the q-score historgram and sum*/
    
          uiQScore = (samSTPtr)->readLenUI;
          scoreAryUC = (unsigned char *) samSTPtr->qStr;
    #endif

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun08 Sec04:
    ^   - finsh of q-scores (for byte method, find qscore)
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    #ifdef NOUL
       scoreAryUC = (unsigned char *) samSTPtr->qStr;
    #endif
       
    for(
       #ifdef NOUL
          uiQScore = 0;
       #else
          uiQScore -=
             (samSTPtr->readLenUI & def_modUL_ulCp);
       #endif
       uiQScore < samSTPtr->readLenUI;
       ++uiQScore
    ) { /*Loop: Copy the q-score entries*/
        uiChar = scoreAryUC[uiQScore];
        uiChar -= def_adjQ_samEntry;
        ++samSTPtr->qHistUI[uiChar];
        samSTPtr->sumQUL += uiChar;
    } /*Loop: Copy the q-score entries*/
    
    samSTPtr->meanQF =
       (float) samSTPtr->sumQUL /(samSTPtr)->readLenUI;
    qhistToMed_samEntry((samSTPtr));
} /*findQScores_samEntry*/

/*-------------------------------------------------------\
| Fun09: cpQEntry_samEntry
|   - Copies q-scores from a string into a samEntry
|     structure
| Input:
|   - samSTPtr:
|     o Pionter to sam entry struct to copys q-scores to
|   _ cpQStr:
|     o C-string with q-scores to copy to samSTPtr
|   - blankQHistBl:
|     o 1: Blank q-score vars (histogram/sum/mean/median)
|     o 0: do not blank the q-score variables
| Output:
|   - Mofidies:
|     o qStr in samSTPtry to have the q-scores
|     o medianQF in samSTPtr to have the median q-score
|     o meanQF in samSTPtr to have the mean q-score
|     o qHistUI in samSTPtr to have histogram of q-scores
|     o samQUL in samSTPtr to have sum off all q-scores
|   - Returns
|     o The value in samSTPtr->readLenUI
\-------------------------------------------------------*/
unsigned int
cpQEntry_samEntry(
   struct samEntry *samSTPtr, /*Copy q-scores to*/
   signed char *cpQStr,       /*q-scores to copy*/
   signed char blankQHistBl   /*1: to blank q-score hist*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun09 TOC: cpQEntry_samEntry
  '   o fun09 sec01:
  '     - Variable declerations
  '   o fun09 sec02:
  '     - Check and if asked blank the q-score values
  '   o fun09 sec03:
  '     - Copy q-scores using unsigned longs
  '   o fun09 sec04:
  '     - Copy last q-scores I could not copy with
  '       unsigned longs
  '   o fun09 sec05:
  '     - Find the median and mean q-scores
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec01:
  ^   - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  unsigned int uiQ = 0;
  unsigned int uiChar = 0;
  unsigned long qScoreUL = 0;

  #ifdef NOUL
  #else
     unsigned char *tmpStr = 0;
     unsigned long *cpPtrUL = (unsigned long *) cpQStr;
     ulong_ulCp
        qAdjustUL =
           mkDelim_ulCp((signed char) def_adjQ_samEntry);
     unsigned long
        *dupPtrUL = (unsigned long *) samSTPtr->qStr;
  #endif
     
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec02:
  ^   - Check and if asked blank the q-score values
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  if(blankQHistBl)
  { /*If: I need to blank the q-score histograms*/
    for(uiChar = 0; uiChar < def_maxQ_samEntry; ++uiChar)
       samSTPtr->qHistUI[uiChar] = 0;

    samSTPtr->sumQUL = 0;
    samSTPtr->medianQF = 0;
    samSTPtr->meanQF = 0;
  } /*If: I need to blank the q-score histograms*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec03:
  ^   - Copy q-scores using unsigned longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  #ifndef NOUL
     for(
        uiQ = 0;
        uiQ < (samSTPtr->readLenUI >> def_shiftULBy_ulCp);
        ++uiQ
     ) { /*Loop: Copy the q-score entries*/
        dupPtrUL[uiQ] = cpPtrUL[uiQ];
        qScoreUL = dupPtrUL[uiQ] - qAdjustUL;
        tmpStr = (unsigned char *) &qScoreUL;
        
        for(
           uiChar = 0;
           uiChar < def_charInUL_ulCp;
           ++uiChar
        ) { /*Loop: Get the q-score entries*/
           ++samSTPtr->qHistUI[tmpStr[uiChar]];
           samSTPtr->sumQUL += tmpStr[uiChar];
        } /*Loop: Get the q-score entries*/
     } /*Loop: Copy the q-score entries*/
  #endif
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec04:
  ^   - Copy last q-scores I could not copy with unsgined
  ^     longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  uiQ = (samSTPtr)->readLenUI;
  
  for(
     #ifdef NOUL
        uiQ = 0;
     #else
        uiQ -= (samSTPtr->readLenUI & def_modUL_ulCp);
     #endif
     uiQ < samSTPtr->readLenUI;
     ++uiQ
  ) { /*Loop: Copy the q-score entries*/
     samSTPtr->qStr[uiQ] = cpQStr[uiQ];
     qScoreUL = (unsigned char) cpQStr[uiQ];
     qScoreUL -= def_adjQ_samEntry;

     ++samSTPtr->qHistUI[qScoreUL];
     samSTPtr->sumQUL += qScoreUL;
  } /*Loop: Copy the q-score entries*/
  
  samSTPtr->qStr[uiQ] = '\0';
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec05:
  ^   - Find the median and mean q-scores
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  samSTPtr->meanQF = (float) samSTPtr->sumQUL;
  samSTPtr->meanQF /= (float) samSTPtr->readLenUI;
  
  qhistToMed_samEntry(samSTPtr);
  return uiQ;
} /*cpQEntry_samEntry*/

/*-------------------------------------------------------\
| Fun11: lineTo_samEntry
|  - adds a sam file line to a samEntry struct
| Input:
|  - samSTPtr:
|    o pionter to samEntry structure to store the line
|    o needs to be initialized and setup (first use only)
|  - buffStr:
|    o c-string with sam file entry for samSTPtr
| Output:
|  - Modifies:
|    o samSTPtr to have the next line
|      - comments are in extraStr
|  - Returns:
|    o 0 for success
|    o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
lineTo_samEntry(
   struct samEntry *samSTPtr,
   signed char *buffStr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - adds a sam file line to a samEntry struct
   '   o fun11 sec01:
   '     - variable declerations
   '   o fun11 sec02:
   '     - get query id from the buffer
   '   o fun11 sec03:
   '     - get flag
   '   o fun11 sec04:
   '     - read in reference name/id
   '   o fun11 sec05:
   '     - get reference position
   '   o fun11 sec06:
   '     - get mapping quality
   '   o fun11 sec07:
   '     - get cigar entry
   '   o fun11 sec08:
   '     - get RNEXT entry
   '   o fun11 sec09:
   '     - get PNEXT entry
   '   o fun11 sec10:
   '     - get TLEN entry
   '   o fun11 sec11:
   '     - get sequence entry
   '   o fun11 sec12:
   '     - get q-score entry
   '   o fun11 sec13:
   '     - copy extra entry; after strict sam entries
   '   o fun11 sec14:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *tmpStr = 0;
   signed int tmpSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - get the query id from the buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);

   /*This is a comment*/
   if(buffStr[0] == '@')
      goto extraEntry;
 
   /*Query id, single byte copy may be better here,
    ` hard to tell
   */

   tmpSI = lenStrNull_ulCp(buffStr, def_tab_ulCp, '\t');
   if(tmpSI > 120)
      goto memErr_fun11_sec14;

   cpLen_ulCp(samSTPtr->qryIdStr, buffStr, tmpSI);
   samSTPtr->qryIdLenUC = (unsigned char) tmpSI;
   buffStr += samSTPtr->qryIdLenUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - get the flag for the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr +=
      strToUS_base10str(
         buffStr,
         &samSTPtr->flagUS
      ); /*Get the flag for the alignment*/

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^   - read in the reference name/id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpSI = lenStrNull_ulCp(buffStr, def_tab_ulCp, '\t');
   if(tmpSI > 120)
      goto memErr_fun11_sec14;

   cpLen_ulCp(samSTPtr->refIdStr, buffStr, tmpSI);
   samSTPtr->refIdLenUC = (unsigned char) tmpSI;
   buffStr += samSTPtr->refIdLenUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^   - get reference position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr +=
      strToUI_base10str(
         buffStr,
         &samSTPtr->refStartUI
      ); /*Get the starting base in the reference*/

    /*Convert the starting positionto index 0*/
   samSTPtr->refStartUI -= (samSTPtr->refStartUI > 0);

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec06:
   ^   - get mapping quality
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr +=
      strToUC_base10str(
         buffStr,
         &samSTPtr->mapqUC
      ); /*Get the mapping quality of the alignment*/

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec07:
   ^   - get cigar entry
   ^   o fun11 sec07 sub01:
   ^     - check if there is a cigar entry
   ^   o fun11 sec07 sub02:
   ^     - read in the cigar entry
   ^   o fun11 sec07 sub03:
   ^     - count number of matchs/snps/dels/inss/masks in
   ^       the cigar entry
   ^   o fun11 sec07 sub04:
   ^     - get read lengths from the cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec07 Sub01:
   *   - check if there is a cigar entry
   \*****************************************************/

   if(buffStr[0] == '*')
   { /*If: the cigar entry was not present*/
      samSTPtr->cigTypeStr[0] = '*';
      samSTPtr->cigArySI[0] = 0;
      buffStr += 2;
      (samSTPtr)->cigLenUI = 1;

      goto rNextEntry_fun11_sec08;
   } /*If: the cigar entry was not present*/

   /*****************************************************\
   * Fun11 Sec07 Sub02:
   *   - read in cigar entry
   \*****************************************************/
   
   while(buffStr[0] > 32)
   { /*Loop: Read in the cigar entry*/

      /*Using -1 to account for the null I will add at the
      `end
      */
      if((samSTPtr)->cigLenUI >= samSTPtr->cigSizeUI-1)
      { /*If: I need to increase the cigar buff size*/
         samSTPtr->cigSizeUI <<= 1;

         tmpStr =
            realloc(
               samSTPtr->cigTypeStr,
               samSTPtr->cigSizeUI * sizeof(char)
            ); /*Resize the type cigar buffer*/

          if(tmpStr == 0)
             goto memErr_fun11_sec14;

          samSTPtr->cigTypeStr = tmpStr;

         tmpStr = (signed char *)
            realloc(
               samSTPtr->cigArySI,
               samSTPtr->cigSizeUI * sizeof(int)
            ); /*Resize the value cigar buffer*/

          if(tmpStr == 0)
             goto memErr_fun11_sec14;

          samSTPtr->cigArySI = (int *) tmpStr;
      } /*If: I need to increase the cigar buff size*/

      /*Record the cigar entry*/
      buffStr +=
         strToSI_base10str(
            buffStr,
            &samSTPtr->cigArySI[samSTPtr->cigLenUI]
          ); /*Get the number of bases for this type*/

      samSTPtr->cigTypeStr[samSTPtr->cigLenUI]=buffStr[0];

      /**************************************************\
      * Fun11 Sec07 Sub03:
      *   - count number of matchs/snps/dels/inss/masks in
      *     cigar entry
      \**************************************************/

      switch(buffStr[0])
      { /*Switch: Check the cigar entry type*/
         case '=':
         /*Case: This was a match*/
             samSTPtr->matchCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This was a match*/

         case 'M':
         /*Case: This is an snp or match*/
             samSTPtr->matchCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This is an snp or match*/

         case 'X':
         /*Case: This is an snp*/
             samSTPtr->snpCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This is an snp*/

         case 'I':
         /*Case: This is an insertion*/
             samSTPtr->insCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This is an insertion*/

         case 'D':
         /*Case: This is an deletion*/
             samSTPtr->delCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This is an deletion*/

         case 'S':
         /*Case: This is an softmasked region*/
             samSTPtr->maskCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];

             break;
         /*Case: This is an softmasked region*/
      } /*Switch: Check the cigar entry type*/

      ++buffStr; /*Get off the tab*/
      ++samSTPtr->cigLenUI;
   } /*Loop: Read in the cigar entry*/

   /*****************************************************\
   * Fun11 Sec07 Sub04:
   *   - get read lengths from the cigar entries
   \*****************************************************/

    samSTPtr->cigArySI[samSTPtr->cigLenUI] = 0;
    samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';

   samSTPtr->readLenUI =
        samSTPtr->matchCntUI
      + samSTPtr->snpCntUI
      + samSTPtr->insCntUI
      + samSTPtr->maskCntUI;

   samSTPtr->alnReadLenUI =
        samSTPtr->matchCntUI
      + samSTPtr->snpCntUI
      + samSTPtr->delCntUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   samSTPtr->refEndUI -= (samSTPtr->alnReadLenUI > 0);
      /*-1 from (alnReadLen > 0) converts to index 0*/

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec08:
   ^   - get RNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   rNextEntry_fun11_sec08:

   /*Not sure which is better here*/
   samSTPtr->rnextLenUC =
      (unsigned char)
      cpDelim_ulCp(
         samSTPtr->rNextStr,
         buffStr,
         def_tab_ulCp,
         '\t'
      ); /*Copy the query id/name*/

   buffStr += samSTPtr->rnextLenUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec09:
   ^   - get PNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr +=
      strToSI_base10str(
         buffStr,
         &samSTPtr->pNextSI
       ); /*Get the number of bases for this type*/

   samSTPtr->pNextSI -= (samSTPtr->pNextSI > 0);

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec10:
   ^   - get TLEN entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr +=
      strToSI_base10str(
         buffStr,
         &samSTPtr->tLenSI
       ); /*Get the number of bases for this type*/

   ++buffStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec11:
   ^   - get sequence entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->readLenUI == 0 && buffStr[0] != '*')
      samSTPtr->readLenUI =
         (unsigned int)
         lenStr_ulCp(buffStr, def_tab_ulCp, '\t');

   else if(buffStr[0] == '*')
   { /*Else If: There  is no sequence entry*/
      samSTPtr->seqStr[0] = '*';
      samSTPtr->seqStr[1] = '\0';
      buffStr += 2;

      goto noQEntry;
   } /*Else If: There  is no sequence entry*/

   if(samSTPtr->readLenUI + 8 > samSTPtr->seqSizeUI)
   { /*If: I need to resize sequence & q-score buffers*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc((samSTPtr->readLenUI + 9) * sizeof(char));

      if(!samSTPtr->seqStr)
         goto memErr_fun11_sec14;

      free(samSTPtr->qStr);
      samSTPtr->qStr = 0;

      samSTPtr->qStr =
         malloc((samSTPtr->readLenUI + 9) * sizeof(char));

      if(!samSTPtr->qStr)
         goto memErr_fun11_sec14;

      samSTPtr->seqSizeUI = samSTPtr->readLenUI;
      samSTPtr->qSizeUI = samSTPtr->readLenUI;
   } /*If: I need to resize sequence & q-score buffers*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      buffStr,
      samSTPtr->readLenUI
   );

   buffStr += samSTPtr->readLenUI + 1; /*+1 gets off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec12:
   ^   - Get q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(buffStr[0] == '*' && buffStr[1] < 32)
   { /*If: there is no q-score entry*/
      noQEntry:

      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
      buffStr += 2;
   } /*If: there is no q-score entry*/

   else
   { /*Else: is a q-score entry*/
      buffStr +=
         cpQEntry_samEntry(
            samSTPtr,
            buffStr,
            0
         ) + 1;
   } /*Else: is a q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec13:
   ^   - copy extra entry; after strict sam entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   extraEntry:

   if(! buffStr || *buffStr == '\0')
   { /*If: no extra entry*/
      if(samSTPtr->extraLenUI > 1)
      { /*If: I need to resize the buffer*/
         if(samSTPtr->extraStr)
            free(samSTPtr->extraStr);

         samSTPtr->extraStr = 0;

         samSTPtr->extraStr =
            malloc(10 * sizeof(signed char));
         
         if(samSTPtr->extraStr == 0)
            goto memErr_fun11_sec14;

         samSTPtr->extraStr[0] = '\0';

         samSTPtr->extraLenUI = 0;
      } /*If: I need to resize the buffer*/
   } /*If: no extra entry*/

   else
   { /*Else: have extra entry*/
      /*not sure if char or ul copy better here*/
      samSTPtr->extraLenUI = endStr_ulCp(buffStr);

      while(buffStr[samSTPtr->extraLenUI] < 33)
      { /*Loop: move off white space*/
         if(! samSTPtr->extraLenUI)
            break; /*extra entry is only white space*/

         --samSTPtr->extraLenUI;
      } /*Loop: move off white space*/

      if(! samSTPtr->extraLenUI)
         samSTPtr->extraStr[0] = '\0';

      else
      { /*Else: have extra entries*/
         ++samSTPtr->extraLenUI;
           /*convert index 0 position to index 1 length*/

         if(
              samSTPtr->extraLenUI
            > samSTPtr->extraSizeUI
         ){ /*If: I need to resize the buffer*/
            free(samSTPtr->extraStr);
            samSTPtr->extraStr = 0;

            samSTPtr->extraStr =
               calloc(
                  samSTPtr->extraLenUI + 9,
                  sizeof(signed char)
               );
            
            if(samSTPtr->extraStr == 0)
               goto memErr_fun11_sec14;
         } /*If: I need to resize the buffer*/

         cpLen_ulCp(
            samSTPtr->extraStr,
            buffStr,
            samSTPtr->extraLenUI
         ); /*Copy the extra entry*/
      } /*Else: have extra entries*/
   } /*Else: have extra entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec14:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun11_sec14:;
      return def_memErr_samEntry;
} /*lineTo_samEntry*/

/*-------------------------------------------------------\
| Fun12: get_samEntry
|  - reads in a single line from a sam file
| Input:
|  - samSTPtr:
|    o pionter to samEntry structure to store the sam file
|      line in. This structure should be initialized
|  - buffStr:
|    o buffer to read the sam file line temporarly into.
|      this is resized if needed. You can input NULL to
|      create a new buffer.
|  - lenBuffUL:
|    o length of buffStr (updated if buffStr is resized)
|  - samVoidFILE:
|    o sam file to read a line from. The void is so that
|      I can use samFILE in the function.
| Output:
|  - Modifies:
|    o samSTPtr to have the next line
|      - Comments are in extraStr
|      - not blanked if at end of file
|    o samFILE to be on the next line
|    o buffStr to hold a sam file line (resized if needed)
|    o lenBuffUL to hold the resized length of buffStr
|  - Returns:
|    o 0 for success
|    o def_EOF_samEntry for EOF (End Of File)
|    o def_memErr_samEntry for memory errors
|    o def_fileErr_samEntry for file errors
\-------------------------------------------------------*/
signed char
get_samEntry(
   struct samEntry *samSTPtr,
   void *samFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC:
   '   - adds a sam file line to a samEntry struct
   '   o fun12 sec01:
   '     - variable declerations
   '   o fun12 sec02:
   '     - get query id from the buffer
   '   o fun12 sec03:
   '     - get flag
   '   o fun12 sec04:
   '     - read in reference name/id
   '   o fun12 sec05:
   '     - get reference position
   '   o fun12 sec06:
   '     - get mapping quality
   '   o fun12 sec07:
   '     - get cigar entry
   '   o fun12 sec08:
   '     - get RNEXT entry
   '   o fun12 sec09:
   '     - get PNEXT entry
   '   o fun12 sec10:
   '     - get TLEN entry
   '   o fun12 sec11:
   '     - get sequence entry
   '   o fun12 sec12:
   '     - get q-score entry
   '   o fun12 sec13:
   '     - copy extra entry; after strict sam entries
   '   o fun12 sec14:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define size_fun12 512
   signed char buffStr[size_fun12];
   signed int posSI = 0;
   signed long lenSL = 0;
   signed long tmpSL = 0;
   signed long ignoreSL = 0;
   signed char *swapStr = 0;
   signed char commentBl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - get the query id from the buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);

   do{ /*Loop: get past white space*/
      lenSL =
         getLine_fileFun(
            samFILE,
            buffStr,
            size_fun12,
            &ignoreSL
         );
      if(! lenSL)
         goto eof_fun12_sec14;
   } while(buffStr[0] == '\n' || buffStr[0] == '\r');
     /*Loop: get past blank lines*/


   /*This is a comment*/
   if(buffStr[0] == '@')
   { /*If: comment*/
      commentBl = 1;
      goto extraEntry_fun12_sec13_sub01;
   } /*If: comment*/
 

   tmpSL = lenStrNull_ulCp(buffStr, def_tab_ulCp, '\t');
   if(tmpSL > 120)
      goto memErr_fun12_sec14;
   cpLen_ulCp( samSTPtr->qryIdStr, buffStr, tmpSL);

   samSTPtr->qryIdLenUC = (unsigned char) tmpSL;
   posSI = tmpSL;

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI; /*get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - get the flag for the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSI +=
     strToUS_base10str(&buffStr[posSI],&samSTPtr->flagUS);

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec04:
   ^   - read in the reference name/id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpSL = lenStrNull_ulCp(buffStr, def_tab_ulCp, '\t');
   if(tmpSL > 120)
      goto memErr_fun12_sec14;

   tmpSL =
     lenStrNull_ulCp(&buffStr[posSI], def_tab_ulCp, '\t');
   if(tmpSL > 120)
      goto memErr_fun12_sec14;
   cpLen_ulCp(samSTPtr->refIdStr, &buffStr[posSI], tmpSL);

   samSTPtr->refIdLenUC = (unsigned char) tmpSL;
   posSI += tmpSL;

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI; /*get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec05:
   ^   - get reference position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSI +=
      strToUI_base10str(
         &buffStr[posSI],
         &samSTPtr->refStartUI
      ); /*Get the starting base in the reference*/

   /*Convert the starting positionto index 0*/
   samSTPtr->refStartUI -= (samSTPtr->refStartUI > 0);

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec06:
   ^   - get mapping quality
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSI +=
      strToUC_base10str(
         &buffStr[posSI],
         &samSTPtr->mapqUC
      ); /*get mapping quality of alignment*/

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec07:
   ^   - get cigar entry
   ^   o fun12 sec07 sub01:
   ^     - check if there is a cigar entry
   ^   o fun12 sec07 sub02:
   ^     - read in the cigar entry
   ^   o fun12 sec07 sub03:
   ^     - count number of matchs/snps/dels/inss/masks in
   ^       the cigar entry
   ^   o fun12 sec07 sub04:
   ^     - get read lengths from the cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec07 Sub01:
   *   - check if there is a cigar entry
   \*****************************************************/

   if(buffStr[posSI] == '*')
   { /*If: the cigar entry was not present*/
      samSTPtr->cigTypeStr[0] = '*';
      samSTPtr->cigArySI[0] = 0;
      samSTPtr->cigTypeStr[1] = 0;
      samSTPtr->cigLenUI = 1;

      ++posSI;
      goto rNextEntry_fun12_sec08;
   } /*If: the cigar entry was not present*/

   /*****************************************************\
   * Fun12 Sec07 Sub02:
   *   - read in cigar entry
   \*****************************************************/
   
   while(buffStr[posSI] > 32)
   { /*Loop: get cigar entry*/
      if(posSI + 16 >= size_fun12)
      { /*If: need more buffer*/
         lenSL -= posSI;
         cpLen_ulCp(
            buffStr,
            &buffStr[posSI],
            lenSL
         );

         lenSL +=
            getLine_fileFun(
               samFILE,
               &buffStr[lenSL],
               size_fun12 - lenSL,
               &ignoreSL
            );

         if(! lenSL)
            goto fileErr_fun12_sec14;
         posSI = 0;
      } /*If: need more buffer*/


      /*Using -1 to account for the null I will add at the
      `end
      */
      if((samSTPtr)->cigLenUI >= samSTPtr->cigSizeUI-1)
      { /*If: I need to increase the cigar buff size*/
         samSTPtr->cigSizeUI <<= 1;

         swapStr =
            realloc(
               samSTPtr->cigTypeStr,
               samSTPtr->cigSizeUI * sizeof(char)
            ); /*Resize the type cigar buffer*/
          if(swapStr == 0)
             goto memErr_fun12_sec14;
          samSTPtr->cigTypeStr = swapStr;

         swapStr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               samSTPtr->cigSizeUI * sizeof(int)
            ); /*Resize the value cigar buffer*/
          if(swapStr == 0)
             goto memErr_fun12_sec14;
          samSTPtr->cigArySI = (signed int *) swapStr;
      } /*If: I need to increase the cigar buff size*/

      /*Record the cigar entry*/
      posSI +=
         strToSI_base10str(
            &buffStr[posSI],
            &samSTPtr->cigArySI[samSTPtr->cigLenUI]
          ); /*Get the number of bases for this type*/
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] =
         buffStr[posSI];

      /**************************************************\
      * Fun12 Sec07 Sub03:
      *   - count number of matchs/snps/dels/inss/masks in
      *     cigar entry
      \**************************************************/

      switch(buffStr[posSI])
      { /*Switch: Check the cigar entry type*/
         case '=':
         /*Case: This was a match*/
             samSTPtr->matchCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This was a match*/

         case 'M':
         /*Case: This is an snp or match*/
             samSTPtr->matchCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This is an snp or match*/

         case 'X':
         /*Case: This is an snp*/
             samSTPtr->snpCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This is an snp*/

         case 'I':
         /*Case: This is an insertion*/
             samSTPtr->insCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This is an insertion*/

         case 'D':
         /*Case: This is an deletion*/
             samSTPtr->delCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This is an deletion*/

         case 'S':
         /*Case: This is an softmasked region*/
             samSTPtr->maskCntUI +=
                samSTPtr->cigArySI[samSTPtr->cigLenUI];
             break;
         /*Case: This is an softmasked region*/
      } /*Switch: Check the cigar entry type*/

      ++posSI; /*move to next cigar entry*/
      ++samSTPtr->cigLenUI;
   } /*Loop: get cigar entry*/

   /*****************************************************\
   * Fun12 Sec07 Sub04:
   *   - get read lengths from the cigar entries
   \*****************************************************/

   samSTPtr->cigArySI[samSTPtr->cigLenUI] = 0;
   samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';

   samSTPtr->readLenUI =
        samSTPtr->matchCntUI
      + samSTPtr->snpCntUI
      + samSTPtr->insCntUI
      + samSTPtr->maskCntUI;

   samSTPtr->alnReadLenUI =
        samSTPtr->matchCntUI
      + samSTPtr->snpCntUI
      + samSTPtr->delCntUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   samSTPtr->refEndUI -= (samSTPtr->alnReadLenUI > 0);
      /*-1 from (alnReadLen > 0) converts to index 0*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec08:
   ^   - get RNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   rNextEntry_fun12_sec08:

   /*making sure have enough data for next three entries*/
   if(posSI + 256 < lenSL) ;

   else if(
         buffStr[lenSL - 1] == '\n'
      || buffStr[lenSL - 1] == '\r'
   ) ; /*have full line*/

   else
   { /*Else: need more data*/
      lenSL -= posSI;
      cpLen_ulCp(buffStr, &buffStr[posSI], lenSL);

      lenSL +=
         getLine_fileFun(
            samFILE,
            &buffStr[lenSL],
            size_fun12 - lenSL,
            &ignoreSL
         );

      if(! lenSL)
         goto fileErr_fun12_sec14;
      posSI = 0;
   } /*Else: need more data*/

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI;


   tmpSL =
     lenStrNull_ulCp(&buffStr[posSI], def_tab_ulCp, '\t');
   if(tmpSL > 120)
      goto memErr_fun12_sec14;

   cpLen_ulCp(samSTPtr->rNextStr, &buffStr[posSI], tmpSL);
   samSTPtr->rnextLenUC = (unsigned char) tmpSL;
   posSI += tmpSL;

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec09:
   ^   - get PNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSI +=
      strToSI_base10str(
         &buffStr[posSI],
         &samSTPtr->pNextSI
      );
   samSTPtr->pNextSI -= (samSTPtr->pNextSI > 0);

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec10:
   ^   - get TLEN entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   posSI +=
      strToSI_base10str(
         &buffStr[posSI],
         &samSTPtr->tLenSI
       ); /*Get the number of bases for this type*/

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec11:
   ^   - get sequence entry
   ^   o fun12 sec11 sub01:
   ^     - no sequence case
   ^   o fun12 sec11 sub02:
   ^     - get sequence in buffer
   ^   o fun12 sec11 sub03:
   ^     - get unkown length (umapped) sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec11 Sub01:
   *   - no sequence case
   \*****************************************************/

   if(buffStr[posSI] == '*')
   { /*If: no sequence in file*/
      samSTPtr->seqStr[0] = '*';
      samSTPtr->seqStr[1] = '\0';
      posSI += 2;
      samSTPtr->readLenUI = 0;
      goto noQEntry_fun12_sec12;
   } /*If: no sequence in file*/

   /*****************************************************\
   * Fun12 Sec11 Sub02:
   *   - get sequence in buffer
   \*****************************************************/

   tmpSL =
     lenStrNull_ulCp(&buffStr[posSI], def_tab_ulCp, '\t');

   if(tmpSL + posSI < lenSL)
   { /*If: went beyond sequence*/
      cpLen_ulCp(samSTPtr->seqStr,&buffStr[posSI],tmpSL);
      samSTPtr->readLenUI = tmpSL;
      tmpSL += posSI;
      lenSL -= tmpSL;
      cpLen_ulCp(buffStr, &buffStr[tmpSL], lenSL);

      posSI = 0;
      tmpSL = 0;
   } /*If: went beyond sequence*/

   /*****************************************************\
   * Fun12 Sec11 Sub03:
   *   - get unkown length (umapped) sequence
   \*****************************************************/

   else if(! samSTPtr->readLenUI)
   { /*Else If: umapped read (need to find length)*/
      lenSL -= posSI;
      cpLen_ulCp(samSTPtr->seqStr,&buffStr[posSI],lenSL);

      posSI = 0;
      samSTPtr->readLenUI = lenSL;

      while(samSTPtr->seqStr[posSI] != '\t')
      { /*Loop: get sequence*/
         if(
               samSTPtr->readLenUI + size_fun12
            >= samSTPtr->seqSizeUI
         ){ /*If: need more memory*/
            samSTPtr->seqSizeUI +=
               (samSTPtr->seqSizeUI >> 1);
            swapStr =
               realloc(
                  samSTPtr->seqStr,
                    (samSTPtr->seqSizeUI + 8)
                  * sizeof(signed char)
               );
            if(! swapStr)
               goto memErr_fun12_sec14;
            samSTPtr->seqStr = swapStr;
         } /*If: need more memory*/
  
         lenSL =
            getLine_fileFun(
               samFILE,
               &samSTPtr->seqStr[samSTPtr->readLenUI],
               size_fun12,
               &ignoreSL
            ); /*this may be a bit slower, but copying
               `  in buffer size lengths should allow me
               `  to never overrun my buffer
               */

         if(! lenSL)
            goto fileErr_fun12_sec14;
         tmpSL =
            lenStrNull_ulCp(
               &samSTPtr->seqStr[samSTPtr->readLenUI],
               def_tab_ulCp,
               '\t'
            );

         samSTPtr->readLenUI += tmpSL;
         if(samSTPtr->seqStr[samSTPtr->readLenUI] == '\t')
            break;
      } /*Loop: get sequence*/

      lenSL -= tmpSL;
      cpLen_ulCp(
         buffStr,
         &samSTPtr->seqStr[samSTPtr->readLenUI],
         lenSL
      );
      posSI = 0;
      samSTPtr->seqStr[samSTPtr->readLenUI] = 0;
   } /*Else If: umapped read (need to find length)*/

   /*****************************************************\
   * Fun12 Sec11 Sub04:
   *   - get known length sequence
   \*****************************************************/

   else
   { /*Else: know length of read*/
      if(samSTPtr->seqSizeUI < samSTPtr->readLenUI)
      { /*If: need more buffer*/
         free(samSTPtr->seqStr);
         samSTPtr->seqStr =
            malloc(
               (samSTPtr->readLenUI + 8)
                  * sizeof(signed char)
            );
         if(! samSTPtr->seqStr)
            goto memErr_fun12_sec14;
         samSTPtr->seqSizeUI = samSTPtr->readLenUI;
      } /*If: need more buffer*/

      cpLen_ulCp(samSTPtr->seqStr,&buffStr[posSI],tmpSL);
         /*tmpSL has current length of buffer*/
      lenSL = samSTPtr->readLenUI - tmpSL;

      tmpSL +=
         fread(
            &samSTPtr->seqStr[tmpSL],
            sizeof(signed char),
            lenSL,
            (FILE *) samFILE
         );

      if(tmpSL < samSTPtr->readLenUI)
         goto fileErr_fun12_sec14;
         /*missing a good chunk of the sequence*/
      lenSL = 0;
      tmpSL = 0;
      posSI = 0;
      samSTPtr->seqStr[samSTPtr->readLenUI] = 0;
   } /*Else: know length of read*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec12:
   ^   - Get q-score entry
   ^   o fun12 sec12 sub01:
   ^     - make sure have part of q-score in buffer
   ^   o fun12 sec12 sub02:
   ^     - no q-score entry case
   ^   o fun12 sec12 sub03:
   ^     - have q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec12 Sub01:
   *   - make sure have part of q-score in buffer
   \*****************************************************/

   if(posSI + 256 < lenSL) ;

   else if(
         buffStr[lenSL - 1] == '\n'
      || buffStr[lenSL - 1] == '\r'
   ) ; /*have full line*/

   else
   { /*Else: need to get more data*/
      lenSL -= posSI;
      cpLen_ulCp(buffStr, &buffStr[posSI], lenSL);

      lenSL +=
         getLine_fileFun(
            samFILE,
            &buffStr[lenSL],
            size_fun12 - lenSL,
            &ignoreSL
         );

      posSI = 0;
   } /*Else: need to get more data*/

   if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;
   ++posSI;

   /*****************************************************\
   * Fun12 Sec12 Sub02:
   *   - no q-score entry case
   \*****************************************************/

   if(buffStr[posSI] == '*' && buffStr[posSI + 1] < 32)
   { /*If: there is no q-score entry*/
      noQEntry_fun12_sec12:

      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
      posSI += 1; /*tab is checked for later*/
   } /*If: there is no q-score entry*/

   /*****************************************************\
   * Fun12 Sec12 Sub03:
   *   - have q-score entry
   \*****************************************************/

   else
   { /*Else: is a q-score entry*/
      if(samSTPtr->qSizeUI < samSTPtr->readLenUI)
      { /*If: need more buffer*/
         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;
         samSTPtr->qStr =
            malloc(
               (samSTPtr->readLenUI + 8)
                 * sizeof(signed char)
            );
         if(! samSTPtr->qStr)
            goto memErr_fun12_sec14;
         samSTPtr->qSizeUI = samSTPtr->readLenUI;
      } /*If: need more buffer*/


      tmpSL = lenSL - posSI;
      if(tmpSL >= (signed long) samSTPtr->readLenUI)
      { /*If: allof q-score entry is in buffer*/
         tmpSL = samSTPtr->readLenUI;
         cpLen_ulCp(
            samSTPtr->qStr,
            &buffStr[posSI],
            tmpSL
         );
         posSI += tmpSL;
      } /*If: allof q-score entry is in buffer*/

      else
      { /*Else: only part of q-score entry is in buffer*/
         cpLen_ulCp(
            samSTPtr->qStr,
            &buffStr[posSI],
            tmpSL
         );

         tmpSL +=
            fread(
               (char *) &samSTPtr->qStr[tmpSL],
               sizeof(signed char),
               samSTPtr->readLenUI - tmpSL,
               (FILE *) samFILE
            );

         if(tmpSL < samSTPtr->readLenUI)
            goto fileErr_fun12_sec14;
         samSTPtr->qStr[samSTPtr->readLenUI] = 0;
         posSI = 0;
         tmpSL = 0;
         lenSL = 0;
      } /*Else: only part of q-score entry is in buffer*/

      findQScores_samEntry(samSTPtr);
   } /*Else: is a q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec13:
   ^   - copy extra entry; after strict sam entries
   ^   o fun12 sec13 sub01:
   ^     - check if have extra entry in buffer
   ^   o fun12 sec13 sub02:
   ^     - copy extra entry in buffer
   ^   o fun12 sec13 sub03:
   ^     - check for line breaks
   ^   o fun12 sec13 sub04:
   ^     - finish getting extra entry line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /*****************************************************\
   * Fun12 Sec13 Sub01:
   *   - check if have extra entry in buffer
   \*****************************************************/

   extraEntry_fun12_sec13_sub01:

   samSTPtr->extraLenUI = 0;

   if(buffStr[posSI] == '\n')
      goto done_fun12_sec14;

   else if(buffStr[posSI] == '\r')
      goto done_fun12_sec14;

   else if(posSI >= lenSL)
   { /*Else If: need more buffer*/
      buffStr[0] = getc((FILE *) samFILE);
      if(buffStr[0] == '\n')
         goto checkNewline_fun12_sec03_sub03;
      else if(buffStr[0] == '\r')
         goto checkCarriage_fun12_sec03_sub03;
      else if(buffStr[0] == '\t')
         goto extraFinish_fun12_sec13_sub04;
      else
         goto fileErr_fun12_sec14;
   } /*Else If: need more buffer*/

   else if(commentBl)
      ;

   else if(buffStr[posSI] != '\t')
      goto fileErr_fun12_sec14;

   else
      ++posSI; /*get off tab*/

   /*****************************************************\
   * Fun12 Sec13 Sub02:
   *   - copy extra entry in buffer
   \*****************************************************/

   if(posSI >= lenSL)
      goto extraFinish_fun12_sec13_sub04;
      /*buffers last character was a tab*/
   else if(posSI < lenSL)
   { /*Else If: need to copy buffer*/
      tmpSL = lenSL - posSI;

      if(tmpSL >= samSTPtr->extraSizeUI)
      { /*If: need more buffer*/ 
         free(samSTPtr->extraStr);
         samSTPtr->extraStr = 0;
         samSTPtr->extraSizeUI += tmpSL;
         samSTPtr->extraSizeUI += (tmpSL >> 1);
         samSTPtr->extraStr =
            malloc(
               (samSTPtr->extraSizeUI + 8)
                 * sizeof(signed char)
            );
         if(! samSTPtr->extraStr)
            goto memErr_fun12_sec14;
      } /*If: need more buffer*/ 

      cpLen_ulCp(
         samSTPtr->extraStr,
         &buffStr[posSI],
         tmpSL
      );
      samSTPtr->extraLenUI = tmpSL;

      if(
            samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
         == '\n'
      ){ /*If: at end of line*/
         --samSTPtr->extraLenUI;
         goto checkNewline_fun12_sec03_sub03;
      }  /*If: at end of line*/

      else if(
            samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
         == '\r'
      ){ /*Else If: at end of line (carraige return)*/
         --samSTPtr->extraLenUI;
         goto checkCarriage_fun12_sec03_sub03;
      }  /*Else If: at end of line (carraige return)*/
   } /*Else If: need to copy buffer*/

   /*****************************************************\
   * Fun12 Sec13 Sub03:
   *   - check for line breaks
   \*****************************************************/


   else
   { /*Else If: get off tab*/
      buffStr[0] = fgetc((FILE *) samFILE);
      if(buffStr[0] == EOF)
         goto done_fun12_sec14;

      else if(buffStr[0] == '\n')
      { /*Else If: end of line*/
         checkNewline_fun12_sec03_sub03:;

         buffStr[0] = fgetc((FILE *) samFILE);
         if(buffStr[0] != '\r')
            ungetc(buffStr[0], (FILE *) samFILE);
         goto done_fun12_sec14;
      } /*Else If: end of line*/

      else if(buffStr[0] == '\r')
      { /*Else If: end of line*/
         checkCarriage_fun12_sec03_sub03:;

         buffStr[0] = fgetc((FILE *) samFILE);
         if(buffStr[0] != '\n')
            ungetc(buffStr[0], (FILE *) samFILE);
         goto done_fun12_sec14;
      } /*Else If: end of line*/

      else if(buffStr[0] != '\t')
         goto fileErr_fun12_sec14;
   } /*Else If: get off tab*/

   /*****************************************************\
   * Fun12 Sec13 Sub04:
   *   - finish getting extra entry line
   \*****************************************************/

   extraFinish_fun12_sec13_sub04:;
      tmpSL = samSTPtr->extraSizeUI;
      samSTPtr->extraLenUI +=
         getFullLine_fileFun(
            samFILE,
            &samSTPtr->extraStr,
            &tmpSL,    /*buffer size*/
            &ignoreSL, /*actual number bytes from file*/
            samSTPtr->extraLenUI /*were to start copy*/
         );
      if(lenSL < 0)
         goto memErr_fun12_sec14;
      samSTPtr->extraSizeUI = tmpSL;

      if(
            samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
         == '\n'
      ){ /*If: newline ending*/
        --samSTPtr->extraLenUI;
        if(
              samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
           == '\r'
        ) --samSTPtr->extraLenUI;
      }  /*If: newline ending*/

      else if(
            samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
         == '\r'
      ){ /*If: carraige return ending*/
        --samSTPtr->extraLenUI;
        if(
              samSTPtr->extraStr[samSTPtr->extraLenUI - 1]
           == '\n'
        ) --samSTPtr->extraLenUI;
      }  /*If: carraige return ending*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec14:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   done_fun12_sec14:;
      samSTPtr->extraStr[samSTPtr->extraLenUI] = '\0';
      return 0;
   eof_fun12_sec14:;
      return def_EOF_samEntry;
   memErr_fun12_sec14:;
      return def_memErr_samEntry;
   fileErr_fun12_sec14:;
      return def_fileErr_samEntry;
} /*get_samEntry*/

/*-------------------------------------------------------\
| Fun13: findRefPos_samEntry
|   - Find an reference coordinate in an sequence in
|     an sam entry structure
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure with sequence
|       (query) to find target reference position in
|   - siCig:
|     o pionter to signed int with index of the cigar
|       entry I am on
|   - cigNtSI:
|     o pionter to a singed int with the number of bases
|       left in the cigar entry
|   - targPosSI:
|     o Reference position I wish to find
|   - refPosSI:
|     o pointer to signed int with the current position
|       in reference sequence
|   - seqPosSI:
|     o pointer to signed int with the current position
|       in query sequence
| Output:
|   - Modifies:
|     o siCig to point to the next open cigar entry
|       - will be > samSTPtr->cigLenUI when the sequence
|         does not end at at targPosSI
|     o cigNtSI to have the number of bases remianing
|       in the current siCig entry
|     o refPosSI to be the index to the position in targSI
|     o seqPosSI to be the index of the sequence base at
|       refPosSI (check cigar to see if deletion)
|   - Returns:
|     o Index of last cigar entry I incurmented
\-------------------------------------------------------*/
signed int
findRefPos_samEntry(
   struct samEntry *samSTPtr,/*has cigar and sequence*/
   signed int *siCig,        /*current cigar entry on*/
   signed int *cigNtSI,      /*# of bases left in cigar*/
   signed int targPosSI,     /*position looking for*/
   signed int *refPosSI,     /*current ref position*/
   signed int *seqPosSI      /*current sequence position*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ` Fun13 TOC: findRefPos_samEntry
  `   - Find an reference coordinate in an sequence in
  `     an sam entry structure
  `   o fun13 sec01:
  `     - Start loop and check insertions/soft masking
  `   o fun13 sec02:
  `     - Move position in deletion cases
  `   o fun13 sec03:
  `     - Move position for snp/match cases
  `   o fun13 sec04:
  `     - Move to the next cigar entry
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - Start loop and check insertions/soft masking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int lastCigMacSI = *siCig;

   if(*refPosSI <= targPosSI)
   { /*If: need to move forwards*/
      if(samSTPtr->cigTypeStr[*siCig] == 'S')
         goto softMask_fun13_sec01;

      if(samSTPtr->cigTypeStr[*siCig] == 'H')
         goto hardMask_fun13_sec03;

      /*Check if I did a full cigar entry move*/
      while(*refPosSI < targPosSI)
      { /*Loop: till I am on the target base*/
         lastCigMacSI = *siCig;

         switch(samSTPtr->cigTypeStr[*siCig])
         { /*Switch: check what the next entry is*/
            case 'S':
            case 'I':
            /*Case: Softmasking or insertions*/
               softMask_fun13_sec01:;
               *seqPosSI += *cigNtSI;
               ++(*siCig);
               *cigNtSI = 0;
               break;
            /*Case: Softmasking or insertions*/

            /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
            ^ Fun13 Sec02:
            ^   - Move position in deletion cases
            \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

            case 'D':
            /*Case: Deletion*/
               *refPosSI += *cigNtSI;

               if(*refPosSI <= (targPosSI))
               { /*If: I have not found target position*/
                  ++(*siCig);
                  *cigNtSI = 0;
               } /*If: I have not found target position*/

               else
               { /*Else: I overshot the target*/
                  *cigNtSI = *refPosSI - targPosSI;
                  *refPosSI -= *cigNtSI;
                     /*correct overshot*/
               } /*Else: I overshot the target*/

               break;
            /*Case: Deletion*/

            /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
            ^ Fun13 Sec03:
            ^   - Move position for snp/match cases\
            \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

            case 'M':
            case 'X':
            case '=':
            /*Case: match (M or =) or snp (M or X)*/
               *refPosSI += *cigNtSI;
               *seqPosSI += *cigNtSI;

               if(*refPosSI <= targPosSI)
               { /*If: I have not found target position*/
                  *cigNtSI = 0;
                  ++(*siCig);
               } /*If: I have not found target position*/

               else 
               { /*Else: I overshot the target*/
                  *cigNtSI = *refPosSI - targPosSI;

                  *refPosSI -= *cigNtSI; /*overshot*/
                  *seqPosSI -= *cigNtSI; /*overshot*/
               } /*Else: I overshot the target*/

               break;
            /*Case: match (M or =) or snp (M or X)*/

            default:
            /*Case: hard mask of some kind*/
               hardMask_fun13_sec03:;
               ++(*siCig);
               *cigNtSI = 0;
            /*Case: hard mask of some kind*/
         } /*Switch: check what the next entry is*/

         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun13 Sec04:
         ^   - Move to the next cigar entry
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

         if(*siCig >= (signed int) samSTPtr->cigLenUI)
         { /*If: end of sequence*/
            if(samSTPtr->cigLenUI > 0)
               *siCig = samSTPtr->cigLenUI - 1;
            else
               *siCig = 0;
            *cigNtSI = 0;
            break;
         } /*If: end of sequence*/

         if(*cigNtSI == 0)
            *cigNtSI = samSTPtr->cigArySI[*siCig];
      } /*Loop: till I am on the target base*/
   } /*If: need to move forwards*/

   else
   { /*Else: need to move backwards*/
      if(samSTPtr->cigTypeStr[*siCig] == 'S')
         goto softMask_fun13_sec0x;

      if(samSTPtr->cigTypeStr[*siCig] == 'H')
         goto hardMask_fun13_sec0x;

      /*Check if I did a full cigar entry move*/
      while(*refPosSI > targPosSI)
      { /*Loop: till I am on the target base*/
         if(! *seqPosSI)
            break; /*end of seqeuence*/

         lastCigMacSI = *siCig;

         switch(samSTPtr->cigTypeStr[*siCig])
         { /*Switch: check what the next entry is*/
            case 'S':
            case 'I':
            /*Case: Softmasking or insertions*/
               softMask_fun13_sec0x:;
               *seqPosSI -=
                  samSTPtr->cigArySI[*siCig] - *cigNtSI;
               --*siCig;
               *cigNtSI = 0;
               break;
            /*Case: Softmasking or insertions*/

            /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
            ^ Fun13 Sec02:
            ^   - Move position in deletion cases
            \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

            case 'D':
            /*Case: Deletion*/
               *refPosSI -=
                  samSTPtr->cigArySI[*siCig] - *cigNtSI;

               if(*refPosSI >= targPosSI)
               { /*If: I have not found target position*/
                  --*siCig;
                  *cigNtSI = 0;
               } /*If: I have not found target position*/

               else
               { /*Else: I overshot the target*/
                  *cigNtSI = targPosSI - *refPosSI;
                  *refPosSI += *cigNtSI;
                     /*correct overshot*/
                  *cigNtSI =
                     samSTPtr->cigArySI[*siCig] -*cigNtSI;
               } /*Else: I overshot the target*/

               break;
            /*Case: Deletion*/

            /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
            ^ Fun13 Sec03:
            ^   - Move position for snp/match cases\
            \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

            case 'M':
            case 'X':
            case '=':
            /*Case: match (M or =) or snp (M or X)*/
               *refPosSI -=
                  samSTPtr->cigArySI[*siCig] - *cigNtSI;
               *seqPosSI -=
                  samSTPtr->cigArySI[*siCig] - *cigNtSI;

               if(*refPosSI >= targPosSI)
               { /*If: I have not found target position*/
                  *cigNtSI = 0;
                  --*siCig;
               } /*If: I have not found target position*/

               else 
               { /*Else: I overshot the target*/
                  *cigNtSI = targPosSI - *refPosSI;
                  *refPosSI += *cigNtSI;
                  *seqPosSI += *cigNtSI; /*overshot*/
                  *cigNtSI =
                     samSTPtr->cigArySI[*siCig] -*cigNtSI;
               } /*Else: I overshot the target*/

               break;
            /*Case: match (M or =) or snp (M or X)*/

            default:
            /*Case: hard mask of some kind*/
               hardMask_fun13_sec0x:;
               --*siCig;
               *cigNtSI = 0;
            /*Case: hard mask of some kind*/
         } /*Switch: check what the next entry is*/

         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun13 Sec04:
         ^   - Move to the next cigar entry
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

         if(*siCig < 0)
         { /*If: at start of seqence*/
            *cigNtSI = samSTPtr->cigArySI[0];
            *siCig = 0;
            break; /*start of the sequence*/
         } /*If: at start of seqence*/
      } /*Loop: till I am on the target base*/
   } /*Else: need to move backwards*/

   return lastCigMacSI;
} /*findRefPos_samEntry*/

/*-------------------------------------------------------\
| Fun14: swap_samEntry
|   - swaps two samEntry structs
| Input:
|   - firstSTPtr:
|     o pointer to first samEntry struct to swap
|   - secSTPtr:
|     o pointer to second samEntry struct to swap
| Output:
|   - Modifies:
|     o all variables in firstSTPtr to be same as secSTPtr
|     o all variables in sedSTPtr to be same as firstSTPtr
\-------------------------------------------------------*/
void
swap_samEntry(
   struct samEntry *firstSTPtr,
   struct samEntry *secSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - swap two samEntry structs
   '   o fun14 sec01:
   '     - swap query id + variable declaration
   '   o fun14 sec02:
   '     - swap reference id
   '   o fun14 sec03:
   '     - swap cigar entries
   '   o fun14 sec04:
   '     - swap rNext entries
   '   o fun14 sec05:
   '     - swap sequence entires
   '   o fun14 sec06:
   '     - swap q-score line
   '   o fun14 sec07:
   '     - swap extra line entires
   '   o fun14 sec08:
   '     - swap mapq, flag, PNEXT, and TLEN
   '   o fun14 sec09:
   '     - swap q-score stat entries
   '   o fun14 sec10:
   '     - swap lengths (read/aligned) and coordinates
   '   o fun14 sec11:
   '     - swap matches, SNPs, indels, and masked stats
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - swap query id + variable declaration
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *swapSCPtr = 0;
   float swapF = 0;
   unsigned short usIndex = 0;


   swapNull_ulCp(
      firstSTPtr->qryIdStr,
      secSTPtr->qryIdStr
   );

   firstSTPtr->qryIdLenUC ^= secSTPtr->qryIdLenUC;
   secSTPtr->qryIdLenUC ^= firstSTPtr->qryIdLenUC;
   firstSTPtr->qryIdLenUC ^= secSTPtr->qryIdLenUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - swap reference id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapNull_ulCp(
      firstSTPtr->refIdStr,
      secSTPtr->refIdStr
   );

   firstSTPtr->refIdLenUC ^= secSTPtr->refIdLenUC;
   secSTPtr->refIdLenUC ^= firstSTPtr->refIdLenUC;
   firstSTPtr->refIdLenUC ^= secSTPtr->refIdLenUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^   - swap cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapSCPtr = firstSTPtr->cigTypeStr;
   firstSTPtr->cigTypeStr = secSTPtr->cigTypeStr;
   secSTPtr->cigTypeStr = swapSCPtr;

   swapSCPtr = (signed char *) firstSTPtr->cigArySI;
   firstSTPtr->cigArySI = secSTPtr->cigArySI;
   secSTPtr->cigArySI = (signed int *) swapSCPtr;
  
   firstSTPtr->cigLenUI ^= secSTPtr->cigLenUI;
   secSTPtr->cigLenUI ^= firstSTPtr->cigLenUI;
   firstSTPtr->cigLenUI ^= secSTPtr->cigLenUI;

   firstSTPtr->cigSizeUI ^= secSTPtr->cigSizeUI;
   secSTPtr->cigSizeUI ^= firstSTPtr->cigSizeUI;
   firstSTPtr->cigSizeUI ^= secSTPtr->cigSizeUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec04:
   ^   - swap rNext entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapNull_ulCp(
      firstSTPtr->rNextStr,
      secSTPtr->rNextStr
   );

   firstSTPtr->rnextLenUC ^= secSTPtr->rnextLenUC;
   secSTPtr->rnextLenUC ^= firstSTPtr->rnextLenUC;
   firstSTPtr->rnextLenUC ^= secSTPtr->rnextLenUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^   - swap sequence entires
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapSCPtr = firstSTPtr->seqStr;
   firstSTPtr->seqStr = secSTPtr->seqStr;
   secSTPtr->seqStr = swapSCPtr;

   firstSTPtr->seqSizeUI ^= secSTPtr->seqSizeUI;
   secSTPtr->seqSizeUI ^= firstSTPtr->seqSizeUI;
   firstSTPtr->seqSizeUI ^= secSTPtr->seqSizeUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec06:
   ^   - swap q-score line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapSCPtr = firstSTPtr->qStr;
   firstSTPtr->qStr = secSTPtr->qStr;
   secSTPtr->qStr = swapSCPtr;

   firstSTPtr->qSizeUI ^= secSTPtr->qSizeUI;
   secSTPtr->qSizeUI ^= firstSTPtr->qSizeUI;
   firstSTPtr->qSizeUI ^= secSTPtr->qSizeUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec07:
   ^   - swap extra line entires
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapSCPtr = firstSTPtr->extraStr;
   firstSTPtr->extraStr = secSTPtr->extraStr;
   secSTPtr->extraStr = swapSCPtr;

   firstSTPtr->extraLenUI ^= secSTPtr->extraLenUI;
   secSTPtr->extraLenUI ^= firstSTPtr->extraLenUI;
   firstSTPtr->extraLenUI ^= secSTPtr->extraLenUI;

   firstSTPtr->extraSizeUI ^= secSTPtr->extraSizeUI;
   secSTPtr->extraSizeUI ^= firstSTPtr->extraSizeUI;
   firstSTPtr->extraSizeUI ^= secSTPtr->extraSizeUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec08:
   ^   - swap mapq, flag, PNEXT, and TLEN
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   firstSTPtr->mapqUC ^= secSTPtr->mapqUC;
   secSTPtr->mapqUC ^= firstSTPtr->mapqUC;
   firstSTPtr->mapqUC ^= secSTPtr->mapqUC;

   firstSTPtr->flagUS ^= secSTPtr->flagUS;
   secSTPtr->flagUS ^= firstSTPtr->flagUS;
   firstSTPtr->flagUS ^= secSTPtr->flagUS;

   firstSTPtr->pNextSI ^= secSTPtr->pNextSI;
   secSTPtr->pNextSI ^= firstSTPtr->pNextSI;
   firstSTPtr->pNextSI ^= secSTPtr->pNextSI;

   firstSTPtr->tLenSI ^= secSTPtr->tLenSI;
   secSTPtr->tLenSI ^= firstSTPtr->tLenSI;
   firstSTPtr->tLenSI ^= secSTPtr->tLenSI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec09:
   ^   - swap q-score stat entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   swapF = firstSTPtr->medianQF;
   firstSTPtr->medianQF = secSTPtr->medianQF;
   secSTPtr->medianQF = swapF;

   swapF = firstSTPtr->meanQF;
   firstSTPtr->meanQF = secSTPtr->meanQF;
   secSTPtr->meanQF = swapF;

   firstSTPtr->sumQUL ^= secSTPtr->sumQUL;
   secSTPtr->sumQUL ^= firstSTPtr->sumQUL;
   firstSTPtr->sumQUL ^= secSTPtr->sumQUL;

   for(
      usIndex = 0;
      usIndex < def_maxQ_samEntry;
      ++usIndex
   ){ /*Loop: swap q-score histogram entries*/
      firstSTPtr->qHistUI[usIndex] ^=
         secSTPtr->qHistUI[usIndex];

      secSTPtr->qHistUI[usIndex] ^=
         firstSTPtr->qHistUI[usIndex];

      firstSTPtr->qHistUI[usIndex] ^=
         secSTPtr->qHistUI[usIndex];
   }  /*Loop: swap q-score histogram entries*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec10:
   ^   - swap lengths (read/aligned) and coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   firstSTPtr->refStartUI ^= secSTPtr->refStartUI;
   secSTPtr->refStartUI ^= firstSTPtr->refStartUI;
   firstSTPtr->refStartUI ^= secSTPtr->refStartUI;

   firstSTPtr->refEndUI ^= secSTPtr->refEndUI;
   secSTPtr->refEndUI ^= firstSTPtr->refEndUI;
   firstSTPtr->refEndUI ^= secSTPtr->refEndUI;

   firstSTPtr->readLenUI ^= secSTPtr->readLenUI;
   secSTPtr->readLenUI ^= firstSTPtr->readLenUI;
   firstSTPtr->readLenUI ^= secSTPtr->readLenUI;

   firstSTPtr->alnReadLenUI ^= secSTPtr->alnReadLenUI;
   secSTPtr->alnReadLenUI ^= firstSTPtr->alnReadLenUI;
   firstSTPtr->alnReadLenUI ^= secSTPtr->alnReadLenUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec11:
   ^   - swap matches, SNPs, indels, and masked stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   firstSTPtr->matchCntUI ^= secSTPtr->matchCntUI;
   secSTPtr->matchCntUI ^= firstSTPtr->matchCntUI;
   firstSTPtr->matchCntUI ^= secSTPtr->matchCntUI;

   firstSTPtr->snpCntUI ^= secSTPtr->snpCntUI;
   secSTPtr->snpCntUI ^= firstSTPtr->snpCntUI;
   firstSTPtr->snpCntUI ^= secSTPtr->snpCntUI;

   firstSTPtr->insCntUI ^= secSTPtr->insCntUI;
   secSTPtr->insCntUI ^= firstSTPtr->insCntUI;
   firstSTPtr->insCntUI ^= secSTPtr->insCntUI;

   firstSTPtr->delCntUI ^= secSTPtr->delCntUI;
   secSTPtr->delCntUI ^= firstSTPtr->delCntUI;
   firstSTPtr->delCntUI ^= secSTPtr->delCntUI;

   firstSTPtr->maskCntUI ^= secSTPtr->maskCntUI;
   secSTPtr->maskCntUI ^= firstSTPtr->maskCntUI;
   firstSTPtr->maskCntUI ^= secSTPtr->maskCntUI;
} /*swap_samEntry*/

/*-------------------------------------------------------\
| Fun15: p_samEntry
|   - prints the sam file entry to a file. This does not
|     print any extra stats that were found.
| Input:
|   - samST
|     o Pointer to samEntry struct with sam entry to print
|  - pNoNewLineBl:
|     o 1: do not print a new line after; you will do this
|     o 0: end of sam entry, print a new line
|   - outFILE:
|     o File to print the sam entry to
| Output:
|  - Prints:
|    o sam file entry in samST to outFILE.
\-------------------------------------------------------*/
void
p_samEntry(
   struct samEntry *samSTPtr,
   signed char pNoNewLineBl,
   void *outFILE
){
   #define size_fun15 512
   #define char_fun15 sizeof(signed char)

   signed char buffStr[size_fun15];
   signed int lenSI = 0;

   unsigned int uiCig = 0;
   FILE *samFILE = (FILE *) outFILE;

   if(
         samSTPtr->extraStr[0] == '@'
      && samSTPtr->qryIdStr[0] == '\0'
   ){ /*If: have header*/
      fprintf(
         samFILE,
         "%s%s",
         samSTPtr->extraStr,
         str_endLine
      ); /*Print out the header (in extra since commnet)*/

      return;
   } /*If: have header*/
   

   /*print query id*/
   samSTPtr->qryIdStr[samSTPtr->qryIdLenUC] = '\t';
   fwrite(
      samSTPtr->qryIdStr,
      char_fun15,
      samSTPtr->qryIdLenUC + 1,
      samFILE
   );
   samSTPtr->qryIdStr[samSTPtr->qryIdLenUC] = '\0';

   /*Copy the flag over*/
   lenSI = numToStr(buffStr, samSTPtr->flagUS);
   buffStr[lenSI++] = '\t';
   fwrite(buffStr, char_fun15, lenSI, samFILE);

   /*Copy the referenced id to the buffer*/
   samSTPtr->refIdStr[samSTPtr->refIdLenUC] = '\t';
   fwrite(
      samSTPtr->refIdStr,
      char_fun15,
      samSTPtr->refIdLenUC + 1,
      samFILE
   );
   samSTPtr->refIdStr[samSTPtr->refIdLenUC] = '\0';

   /*Reference position*/
   if(samSTPtr->flagUS & 4)
   { /*If: unmapped read*/
      fputc('0', (FILE *) samFILE);
      fputc('\t', (FILE *) samFILE);
   } /*If: unmapped read*/

   else
   { /*Else: have starting position on reference*/
      lenSI = numToStr(buffStr, samSTPtr->refStartUI + 1);
      buffStr[lenSI++] = '\t';
      fwrite(buffStr, char_fun15, lenSI, samFILE);
   } /*Else: have starting position on reference*/

   /*mapq*/
   lenSI = numToStr(buffStr, samSTPtr->mapqUC);
   buffStr[lenSI++] = '\t';
   fwrite(buffStr, char_fun15, lenSI, samFILE);

   /*Check if there is a cigar entry*/
   if(
         ! samSTPtr->cigTypeStr[0]
      || samSTPtr->cigTypeStr[0]=='*'
   ){ /*If: no cigar entry*/
      fputc('*', samFILE);
      fputc('\t', samFILE);
   } /*If: no cigar entry*/

   else
   { /*Else: convert the cigar entry*/
      lenSI = 0;

      for(uiCig=0; uiCig < samSTPtr->cigLenUI; ++uiCig)
      { /*Loop: Convert cigar to string*/
         lenSI +=
            numToStr(
               &buffStr[lenSI],
               samSTPtr->cigArySI[uiCig]
            );
         buffStr[lenSI++] = samSTPtr->cigTypeStr[uiCig];

         if(lenSI + 16 >= size_fun15)
         { /*If: need to print buffer*/
            fwrite(buffStr, char_fun15, lenSI, samFILE);
            lenSI = 0;
         } /*If: need to print buffer*/
      } /*Loop: Convert cigar to string*/

      buffStr[lenSI++] = '\t';
      fwrite(buffStr, char_fun15, lenSI, samFILE);
   } /*Else: convert the cigar entry*/

   if(samSTPtr->rNextStr[0] != '\0')
   { /*If: need to write the RNEXT entry*/
      fwrite(
         samSTPtr->rNextStr,
         char_fun15,
         samSTPtr->rnextLenUC,
         samFILE
      );
   } /*If: need to write the RNEXT entry*/

   else
      fputc('*', (FILE *) samFILE);

   fputc('\t', (FILE *) samFILE);

   /*PNEXT*/
   lenSI = signed_numToStr(buffStr, samSTPtr->pNextSI);
   buffStr[lenSI++] = '\t';
   fwrite(buffStr, char_fun15, lenSI, samFILE);

   /*TLEN*/
   lenSI = signed_numToStr(buffStr, samSTPtr->tLenSI);
   buffStr[lenSI++] = '\t';
   fwrite(buffStr, char_fun15, lenSI, samFILE);

   /*sequence print*/
   if(samSTPtr->seqStr[0])
      fwrite(
         samSTPtr->seqStr,
         char_fun15,
         samSTPtr->readLenUI,
         samFILE
      );
   else
      fputc('*', (FILE *) samFILE);
   fputc('\t', (FILE *) samFILE);

   /*q-score entry print*/
   if(samSTPtr->qStr[0] && samSTPtr->qStr[1])
      fwrite(
         samSTPtr->qStr,
         char_fun15,
         samSTPtr->readLenUI,
         samFILE
      );
   else
      fputc('*', samFILE);

   if(! *samSTPtr->extraStr)
      ;
   else if(samSTPtr->extraLenUI)
   { /*Else If: have extra entry to print out*/
      fputc('\t', samFILE);
      fwrite(
         samSTPtr->extraStr,
         char_fun15,
         samSTPtr->extraLenUI,
         samFILE
      );
   } /*Else If: have extra entry to print out*/

   if(! pNoNewLineBl)
   { /*If: user wanted line ending*/
      fputc(str_endLine[0], samFILE);
      if(str_endLine[1])
         fputc(str_endLine[1], samFILE);
   } /*If: user wanted line ending*/
} /*p_samEntry*/

/*-------------------------------------------------------\
| Fun16: pfq_samEntry
|   - Prints the sam entry as fastq entry to a fastq file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fastq file to print the new fastq entry to
| Output:
|   - Prints:
|     o fastq entry from samST to outFILE
\-------------------------------------------------------*/
void
pfq_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if(   (samSTPtr)->seqStr[1] != '\0'
      && (samSTPtr)->qStr[1] != '\0'
   ){ /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        "@%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->matchCntUI,
        (samSTPtr)->snpCntUI,
        (samSTPtr)->insCntUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f%s",
        (samSTPtr)->delCntUI,
        (samSTPtr)->maskCntUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF,
        str_endLine
      );
      fprintf(
        (FILE *) (outFILE),
        "%s%s+%s%s%s",
        (samSTPtr)->seqStr,
        str_endLine,
        str_endLine,
        (samSTPtr)->qStr,
        str_endLine
      );
   } /*If: This entry can be printed out*/
} /*p_samEntryAsFq*/

/*-------------------------------------------------------\
| Fun17: pfa_samEntry
|   - Prints the sam entry as fasta entry to a fasta file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fasta file to print the new fasta entry to
| Output:
|   - Prints:
|     o fasta entry from samST to outFILE
\-------------------------------------------------------*/
void
pfa_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if((samSTPtr)->seqStr[1] != '\0')
   { /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        ">%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->matchCntUI,
        (samSTPtr)->snpCntUI,
        (samSTPtr)->insCntUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f%s",
        (samSTPtr)->delCntUI,
        (samSTPtr)->maskCntUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF,
        str_endLine
      );

      fprintf(
         (FILE *) (outFILE),
         "%s%s",
         (samSTPtr)->seqStr,
         str_endLine
      );
   } /*If: This entry can be printed out*/
} /*pfa_samEntry*/

/*-------------------------------------------------------\
| Fun18: pstats_samEntry
|   - Prints out the stats in a samEntry struct to a file
| Input:
|   - samEntryST:
|     o Pointer to samEntry struct to print stats for
|   - pHeadBl:
|     o 1: Print the header for the stats tsv file
|     o 0: Do not print the header
|   - pNsBl:
|     o 1: find and print out the anonymous base counts
|     o 0: do not print out anonymous base counts
|   - outFILE:
|     o TSV (tab deliminated) file to print stats to
| Output:
|   - Prints:
|     o line with stats from samEntryStruct to file (tsv)
|   - Modifies:
|     o printHeaderChar to be 0 if set to 1
\-------------------------------------------------------*/
void
pstats_samEntry(
   struct samEntry *samSTPtr,
   signed char *pHeadBl,
   signed char pNsBl,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun18 TOC:
   '   - Prints stats in a samEntry struct to file
   '   o fun18 sec01:
   '     - variable declerations
   '   o fun18 sec02:
   '     - check if comment, if not check if print header
   '   o fun18 sec03:
   '     - print out general stats
   '   o fun18 sec04:
   '     - print matches, snps, ins, dels, and masking.
   '     - if asked also include anonymous bases
   '   o fun18 sec05:
   '     - print out the accuracy stats
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int matchNCntUI = 0;
   unsigned int snpNCntUI = 0;
   unsigned int insNCntUI = 0;
   unsigned int maskNCntUI = 0;
   unsigned int *cntUI = 0;

   signed int siCig = 0;
   signed int numNtSI = 0;
   signed char *tmpStr = 0;
   unsigned char ntUC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec02:
   ^   - check if comment, if not check if print header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(   (samSTPtr)->extraStr[0] != '@'
      && (samSTPtr)->qryIdStr[0] != '\0'
   ){ /*If: This is not a comment*/

      if(*(pHeadBl))
      { /*If: I need to print the header*/
        fprintf((FILE *) (outFILE), "Read\tRef\tFlag");
        fprintf((FILE *) (outFILE), "\tMapQ\tRefPos"); 
        fprintf((FILE *) (outFILE), "\tReadLength");
        fprintf((FILE *) (outFILE), "\tRefAlnLength");

        if((pNsBl))
        { /*If: I am printing the anonymous base counts*/
           fprintf((FILE *) (outFILE), "\tmatch_total");
           fprintf((FILE *) (outFILE), "\tmatch");
           fprintf((FILE *)(outFILE),"\tmatch_anonymous");

           fprintf((FILE *) (outFILE), "\tsnp_total");
           fprintf((FILE *) (outFILE), "\tsnp");
           fprintf((FILE *) (outFILE), "\tsnp_anonymous");

           fprintf((FILE *) (outFILE), "\tins_total");
           fprintf((FILE *) (outFILE), "\tins");
           fprintf((FILE *) (outFILE), "\tins_anonymous");

           fprintf((FILE *) (outFILE), "\tdel");

           fprintf((FILE *) (outFILE), "\tmask_total");
           fprintf((FILE *) (outFILE), "\tmask");
           fprintf((FILE *) (outFILE),"\tmask_anonymous");
        } /*If: I am printing the anonymous base counts*/

        else
        { /*Else: Print out the normal header*/
           fprintf((FILE *) (outFILE), "\tMatches\tSnps");
           fprintf((FILE *) (outFILE), "\tInss\tDels");
           fprintf((FILE *) (outFILE), "\tmask");
        } /*Else: Print out the normal header*/

        fprintf(
           (FILE *) (outFILE),
           "\tMedianQ\tMeanQ%s",
           str_endLine
        );

        *(pHeadBl) = 0;
      } /*If: I need to print the header*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun18 Sec03:
      ^   - print out general stats
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      fprintf(
        (FILE *) (outFILE),
        "%s\t%s\t%u\t%u\t%u\t%u\t%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun18 Sec04:
      ^   - print matches, snps, ins, dels, and masking.
      ^   - if asked also include anonymous bases
      ^   o fun18 sec04 sub01:
      ^     - print anonymous bases; have sequence?
      ^   o fun18 sec04 sub02:
      ^     - find number anonymous bases
      ^   o fun18 sec04 sub03:
      ^     - print out anonymous base and other counts
      ^   o fun18 sec04 sub04:
      ^     - not printing anoynous bases do regular print
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun18 Sec04 Sub01:
      *   - (print anonymous bases) check if have sequence
      \**************************************************/

      if( pNsBl )
      { /*If: I am finding the anonymous base counts*/
         tmpStr = (signed char *) (samSTPtr)->seqStr;

         if(*tmpStr == '*')
         { /*If: no sequence to check for anonymous*/
            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->matchCntUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->snpCntUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA\t%u",
              (samSTPtr)->insCntUI,
              (samSTPtr)->delCntUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->maskCntUI
            );

            goto pqscore_fun18_sec05;
         } /*If: no sequence to check for anonymous*/

         /***********************************************\
         * Fun18 Sec04 Sub02:
         *   - find number anonymous bases
         \***********************************************/

         siCig = 0;

         while(siCig < (signed int) (samSTPtr)->cigLenUI)
         { /*Loop: count number of anonymous bases*/
            numNtSI = (samSTPtr)->cigArySI[siCig];

            switch((samSTPtr)->cigTypeStr[siCig])
            { /*Switch: check what type of cigar entry*/
               case 'I':
               /*Case: Insertions*/
                  cntUI = &insNCntUI;
                  break;
               /*Case: Insertions*/

               case 'X':
               /*Case: Mismatches*/
                  cntUI = &snpNCntUI;
                  break;
               /*Case: Mismatches*/

               case '=':
               case 'M':
               /*Case: Matches or non-eqx entry*/
                  cntUI = &matchNCntUI;
                  break;
               /*Case: Matches or non-eqx entry*/

               case 'S':
               /*Case: softmasking*/
                  cntUI = &maskNCntUI;
                  break;
               /*Case: softmasking*/

               /*Handle missing information cases*/
               default:
               /*Case: deletions or hard masks*/
                  cntUI = 0;
                  numNtSI = 0;
               /*Case: deletions or hard masks*/
            } /*Switch: check what type of cigar entry*/

            while(numNtSI > 0)
            { /*Loop: Check each nucleotide*/
               ntUC = ntTo5Bit[(unsigned char) *tmpStr];
               (*cntUI) +=
                  (!!(ntUC & def_n_fithBit_ntTo5Bit));

               ++tmpStr;
               --numNtSI;
            } /*Loop: Check each nucleotide*/

            ++siCig;
         } /*Loop: count number of anonymous bases*/

         /***********************************************\
         * Fun18 Sec04 Sub03:
         *   - print anonymous base and other counts
         \***********************************************/

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->matchCntUI,
           (samSTPtr)->matchCntUI - matchNCntUI,
           matchNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->snpCntUI,
           (samSTPtr)->snpCntUI - snpNCntUI,
           snpNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u",
           (samSTPtr)->insCntUI,
           (samSTPtr)->insCntUI - insNCntUI,
           insNCntUI,
           (samSTPtr)->delCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->maskCntUI,
           (samSTPtr)->maskCntUI - maskNCntUI,
           maskNCntUI
         );
      } /*If: I am finding the anonymous base counts*/

      /**************************************************\
      * Fun18 Sec04 Sub04:
      *   - not printing anoynous bases do regular print
      \**************************************************/

      else
      { /*Else: I am not counting anonymous bases*/
         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u\t%u",
           samSTPtr->matchCntUI,
           samSTPtr->snpCntUI,
           samSTPtr->insCntUI,
           samSTPtr->delCntUI,
           samSTPtr->maskCntUI
         );
      } /*Else: I am not counting anonymous bases*/
      
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun18 Sec05:
      ^   - print out the q-scores
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      pqscore_fun18_sec05:;

      fprintf(
        (FILE *) (outFILE),
        "\t%f\t%f%s",
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF,
        str_endLine
      );

   } /*If: This is not a comment*/
} /*pstats_samEntry*/

/*-------------------------------------------------------\
| Fun19: revCmp_samEntry
|   - reverse complements a sam file sequence entry
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to reverse complement
| Output:
|   - Modifies:
|     o samSTPtr (includes flag) to be reverse complment
\-------------------------------------------------------*/
void
revCmp_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun19 TOC:
   '   - reverse complements a sam file sequence entry
   '   o fun19 sec01:
   '     - variable declarations
   '   o fun19 sec02:
   '     - check if have sequence and qscore entries
   '   o fun19 sec03:
   '     - reverse complement sequence and q-scores
   '   o fun19 sec04:
   '     - reverse the cigar and flag
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *startStr = 0;
   signed char *endStr = 0;

   signed char *qStartStr = 0;
   signed char *qEndStr = 0;

   signed int *startSIPtr;
   signed int *endSIPtr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec02:
   ^   - check if have sequence and qscore entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! samSTPtr->seqStr
      || (
            samSTPtr->seqStr[0] == '*'
            && samSTPtr->seqStr[1] < 32
         )
   ) ; /*no sequence entry*/

   else
   { /*Else: have sequence entry*/
      startStr = samSTPtr->seqStr;
      endStr = startStr + samSTPtr->readLenUI - 1;
   } /*Else: have sequence entry*/

   if(
         ! samSTPtr->qStr
      || (
            samSTPtr->qStr[0] == '*'
            && samSTPtr->qStr[1] < 32
         )
   ) ;

   else
   { /*Else: have q-score entry*/
      qStartStr = samSTPtr->qStr;
      qEndStr = qStartStr + samSTPtr->readLenUI - 1;
   } /*Else: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec03:
   ^   - reverse complement sequence and q-scores
   ^   o fun19 sec03 sub01:
   ^     - reverse q-score entries + start loop
   ^   o fun19 sec03 sub02:
   ^     - reverse complement sequence entry
   ^   o fun19 sec03 sub03:
   ^     - make sure last sequence base is complement
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun19 Sec03 Sub01:
   *   - reverse q-score entries + start loop
   \*****************************************************/

   while(startStr < endStr)
   { /*Loop: reverse complement*/
      if(qStartStr)
      { /*If: have qscore entry*/
         *qStartStr ^= *qEndStr;
         *qEndStr ^= *qStartStr;
         *qStartStr ^= *qEndStr;

         ++qStartStr;
         --qEndStr;
      } /*If: have qscore entry*/

      /**************************************************\
      * Fun19 Sec03 Sub02:
      *   - reverse complement sequence entry
      \**************************************************/

      switch(*startStr & (~32))
      { /*Switch: complement*/
          case 'A': *startStr = 'T'; break;
          case 'C': *startStr = 'G'; break;
          case 'G': *startStr = 'C'; break;
          case 'T': *startStr = 'A'; break;
          case 'U': *startStr = 'A'; break;
          case 'W': *startStr = 'W'; break;
          case 'S': *startStr = 'S'; break;
          case 'M': *startStr = 'K'; break; /*A/C*/
          case 'K': *startStr = 'M'; break; /*T/G*/
          case 'R': *startStr = 'Y'; break; /*A/G*/
          case 'Y': *startStr = 'R'; break; /*T/C*/
          case 'B': *startStr = 'V'; break; /*C/G/T*/
          case 'D': *startStr = 'H'; break; /*G/T/A*/
          case 'H': *startStr = 'D'; break; /*C/A/T*/
          case 'V': *startStr = 'B'; break; /*A/C/G*/
          case 'N': *startStr = 'N'; break; /*A/C/G/T*/
          default: *startStr = 0;
      } /*Switch: complement*/

      switch(*endStr & (~32))
      { /*Switch: complement*/
          case 'A': *endStr = 'T'; break;
          case 'C': *endStr = 'G'; break;
          case 'G': *endStr = 'C'; break;
          case 'T': *endStr = 'A'; break;
          case 'U': *endStr = 'A'; break;
          case 'W': *endStr = 'W'; break;
          case 'S': *endStr = 'S'; break;
          case 'M': *endStr = 'K'; break; /*A/C*/
          case 'K': *endStr = 'M'; break; /*T/G*/
          case 'R': *endStr = 'Y'; break; /*A/G*/
          case 'Y': *endStr = 'R'; break; /*T/C*/
          case 'B': *endStr = 'V'; break; /*C/G/T*/
          case 'D': *endStr = 'H'; break; /*G/T/A*/
          case 'H': *endStr = 'D'; break; /*C/A/T*/
          case 'V': *endStr = 'B'; break; /*A/C/G*/
          case 'N': *endStr = 'N'; break; /*A/C/G/T*/
          default: *endStr = 0;
      } /*Switch: complement*/

      *startStr ^= *endStr;
      *endStr ^= *startStr;
      *startStr ^= *endStr;

      ++startStr;
      --endStr;
   } /*Loop: reverse complement*/

   /*****************************************************\
   * Fun19 Sec03 Sub03:
   *   - make sure last sequence base is complement
   \*****************************************************/

   if(
         startStr
      && startStr == endStr
   ){ /*If: have middle base*/
      switch(*startStr & (~32))
      { /*Switch: complement*/
          case 'A': *startStr = 'T'; break;
          case 'C': *startStr = 'G'; break;
          case 'G': *startStr = 'C'; break;
          case 'T': *startStr = 'A'; break;
          case 'U': *startStr = 'A'; break;
          case 'W': *startStr = 'W'; break;
          case 'S': *startStr = 'S'; break;
          case 'M': *startStr = 'K'; break; /*A/C*/
          case 'K': *startStr = 'M'; break; /*T/G*/
          case 'R': *startStr = 'Y'; break; /*A/G*/
          case 'Y': *startStr = 'R'; break; /*T/C*/
          case 'B': *startStr = 'V'; break; /*C/G/T*/
          case 'D': *startStr = 'H'; break; /*G/T/A*/
          case 'H': *startStr = 'D'; break; /*C/A/T*/
          case 'V': *startStr = 'B'; break; /*A/C/G*/
          case 'N': *startStr = 'N'; break; /*A/C/G/T*/
          default: *startStr = 0;
      } /*Switch: complement*/
   } /*If: have middle base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec04:
   ^   - reverse the cigar and flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->flagUS & 4) ;
      /*unmapped read; do nothing*/
   else if(samSTPtr->flagUS & 16)
      samSTPtr->flagUS &= ~16;
   else
      samSTPtr->flagUS |= 16;

   if(
         samSTPtr->cigTypeStr
      && samSTPtr->cigTypeStr[0] != '*'
   ){ /*If: have cigar to flip*/
      startStr = samSTPtr->cigTypeStr;
      endStr = startStr + samSTPtr->cigLenUI - 1;

      startSIPtr = samSTPtr->cigArySI;
      endSIPtr = startSIPtr + samSTPtr->cigLenUI - 1;

      while(startStr < endStr)
      { /*Loop: reverse cigar*/
         *startStr ^= *endStr;
         *endStr ^= *startStr;
         *startStr ^= *endStr;

         ++startStr;
         --endStr;

         *startSIPtr ^= *endSIPtr;
         *endSIPtr ^= *startSIPtr;
         *startSIPtr ^= *endSIPtr;

         ++startSIPtr;
         --endSIPtr;
      } /*Loop: reverse cigar*/
   } /*If: have cigar to flip*/
} /*revCmp_samEntry*/

/*-------------------------------------------------------\
| Fun20: getHead_samEntry
|   - get header for a sam file
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to get header with
|   - samFILE:
|     o sam file to get header from
| Output:
|   - Modifies:
|     o samSTPtr to have entry after header
|     o samFILE to point to second read entry
|   - Returns:
|     o pointer to c-string with header
|     o 0 for memory error
\-------------------------------------------------------*/
signed char *
getHead_samEntry(
   struct samEntry *samSTPtr,
   void *samFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun20 TOC:
   '   - get header for a sam file
   '   o fun20 sec01:
   '     - variable declarations
   '   o fun20 sec02:
   '     - memory allocation and get first line
   '   o fun20 sec03:
   '     - get header
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec01:
   ^   - varaible declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define len_fun20_samEntry 4096

   signed char errSC = 0;
   signed char *retHeapStr = 0;
   signed char *tmpStr = 0;
   unsigned int lenHeadUI = 0;
   unsigned int sizeRetUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec02:
   ^   - memory allocation and get first line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retHeapStr =
      malloc(len_fun20_samEntry * sizeof(signed char));
   sizeRetUI = len_fun20_samEntry;

   if(! retHeapStr)
      goto memErr_fun20_sec04;

   /*get first line of sam file*/
   errSC = get_samEntry(samSTPtr, samFILE);

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_samEntry)
         goto memErr_fun20_sec04;
      else
         goto memErr_fun20_sec04; /*EOF or file error*/
   } /*If: had error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec03:
   ^   - get header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errSC)
   { /*Loop: get header*/
      if(samSTPtr->extraStr[0] != '@')
         break; /*end of comments*/

      if(lenHeadUI + samSTPtr->extraLenUI > sizeRetUI - 1)
      { /*If: need to resize buffer*/
         sizeRetUI = lenHeadUI + len_fun20_samEntry;
         sizeRetUI += samSTPtr->extraLenUI;

         tmpStr =
            realloc(
               retHeapStr,
               sizeRetUI
            );

         if(! tmpStr)
            goto memErr_fun20_sec04;

         retHeapStr = tmpStr;
         tmpStr = 0;
      } /*If: need to resize buffer*/

      cpLen_ulCp(
         &retHeapStr[lenHeadUI],
         samSTPtr->extraStr,
         samSTPtr->extraLenUI
      );

      lenHeadUI += samSTPtr->extraLenUI;
      retHeapStr[lenHeadUI++] = *str_endLine;

      if(str_endLine[1])
         retHeapStr[lenHeadUI++] = str_endLine[1];
         /*windows has two line characters for line end*/

      retHeapStr[lenHeadUI] = '\0';

      /*get next line of sam file*/
      errSC = get_samEntry(samSTPtr, samFILE);
   } /*Loop: get header*/

   if(errSC == def_memErr_samEntry)
      goto memErr_fun20_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun20_sec04;

   memErr_fun20_sec04:;
   if(retHeapStr)
      free(retHeapStr);
   retHeapStr = 0;
   goto ret_fun20_sec04;

   ret_fun20_sec04:;
   return retHeapStr;
} /*getHead_samEntry*/

/*-------------------------------------------------------\
| Fun21: checkCigSize_samEntry
|   - checks to see if I need to resize cigar and if
|     needed resizes the cigar
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with cigar to check
|   - cigPosUI:
|     o position of next cigar entry to add
| Output:
|   - Modifies:
|     o cigTypeStr in samSTPtr to be resized if needed
|     o cigArySI in samSTPtr to be resized if needed
|     o cigSizeUI in samSTPtr to have new cigar size
|   - Returns:
|     o 0 for success
|     o 1 for memory error
\-------------------------------------------------------*/
signed char
checkCigSize_samEntry(
   struct samEntry *samSTPtr, /*has cigar*/
   unsigned int cigPosUI      /*positon at in cigar*/
){ /*checkCigSize_samEntry*/
   signed char *swapStr = 0;
   unsigned int newSizeUI = 0;

   if(cigPosUI >= samSTPtr->cigSizeUI)
   { /*If: need more memory*/
      newSizeUI = samSTPtr->cigSizeUI;
      newSizeUI += (newSizeUI >> 1);

      swapStr =
         realloc(
            samSTPtr->cigTypeStr,
            newSizeUI * sizeof(signed char)
         );
      if(! swapStr)
         goto memErr_fun21;
      samSTPtr->cigTypeStr = swapStr;

      swapStr =
         (signed char *)
         realloc(
            samSTPtr->cigArySI,
            newSizeUI * sizeof(signed int)
         );
      if(! swapStr)
         goto memErr_fun21;
      samSTPtr->cigArySI = (signed int *) swapStr;

      samSTPtr->cigSizeUI = newSizeUI;
   } /*If: need more memory*/

   return 0;

   memErr_fun21:;
      return 1;
}  /*checkCigSize_samEntry*/

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
